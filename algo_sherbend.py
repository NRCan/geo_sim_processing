#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
    This algorithm CIT-S implements the Wang Generalization algotithm with constraint checking

    This algorithm simplifies lines.  It detects for each line the bends.  It analyze the bend and
    remove the bends that are below a certain diameter. The point and lines that do not need 
    to be simplified are still used to enforce topology integrity between  those feature that need to be simplified
    
    Usage:
        import spike

    Limits and constraints
        Always works better when the line to process meet the OGC simple line.
          

"""

#####################################################################################################################################

import math, sys

from shapely.geometry import Point, LineString, LinearRing, Polygon
from shapely import affinity
from lib_geosim import GenUtil, SpatialContainer
                               
# Public key word contants
#MULTI_BENDS = "MULTI_BEND"
#SINGLE_BEND = "SINGLE_BEND"
#NO_VERTICE_ADD = "NO_VERTICE_ADD"

# Properties name
_DIAMETER = "diameter"
#_SIMPLIFY_FIRST_LAST = "simplify_first_last"

# Internal constant

_ALGO = "Sherbend"

# Internal constant ===> Should be modify with care...
_AREA_CMP_INDEX = .75                   # Compactness index factor applied to the adjusted area
_SIMILAR_BEND_ADJ_AREA_RATIO = .5       # Adjusted area ratio to detect similar bend; From [0..1] higher value mean 
                                        # more similar bends are accepted
_SIMILAR_BEND_BASE_RATIO = .5           # Base ratio to detect similar bend; From [0..1] Higher value mean more
                                        # similar bends are accepted
_SIMILAR_BEND_CMP_INDEX_RATIO = .6      # Compactness index ratio to detect similar bend; From [0..1] higher 
                                        # value mean more similar bends are accepted
_DEPTH_OFFSET_RATIO = 0.5               # When _BIG_BEND=True, it is the ratio by which the point in the
                                        # middle is offset in direction of the bend peak
_BIG_BEND_CMP_INDEX_RATIO = 0.6         # Compactness index ratio for a bend to be considered as a big bend
                                        # candidate; [0..1] A higher value will allow more bend to be BIG_BEND  candidate
_BIG_BEND_MAX_ADJ_AREA_RATIO = 0.6      # Minimal adjusted area ratio used to detect if a bend meets the minimum
                                        # area value; [0..1] A higher value will allow more bend to be BIG_BEND
                                        # candidate

#Internal key word constants
_BIG = "Big"
_BURNED = "Burned"
_SMALL = "Small"
_SIMPLIFIED = 'Simplified'
_NOT_SIMPLIFIED = 'NotSimplified'
_UNSIMPLIFIABLE = 'Unsimplifiable'


class LineStringSb(LineString):

    """LineString specialization for the SherBend algorithm"""

    def __init__(self, coords, fast_access=True):
        super().__init__(coords)
        self.fast_access = fast_access
        if self.fast_access:
            self.__lst_coords = list(super().coords)

        # Declaration of the instance variable
        self.sb_geom_type = self.geom_type # variable defined to avoid slower C calls with geom_type
        self.sb_is_simplest = False # The line is not at its simplest form
        self.sb_bends = [] # Holder for the bend of the line


    @property
    # Is the line string closed
    def sb_is_closed(self):
        try:
            return self._sb_is_closed
        except AttributeError:
            # A closed line need at least 4 vertex to be valid
            if len(self.coords) >=4 and  GenUtil.distance(self.coords[0], self.coords[-1]) <= GenUtil.ZERO:
                self._sb_is_closed = True
            else:
                self._sb_is_closed = False
            return self._sb_is_closed

    @property
    def coords(self):
        if self.fast_access:
            return self.__lst_coords
        else:
            return super().coords


    @coords.setter
    def coords(self, coords):
        LineString.coords.__set__(self, coords)
        if self.fast_access:
            self.__lst_coords = list(super().coords)
        # Delete variable that are now outdated. so they will be computed next time it will be accessed
        try:
            del self._vertex_orientation
        except AttributeError:
            pass


    @property
    def vertex_orientation(self):
        """List containing the orientation at each vertex of the line.
        -1: anti clockwise, +1 Clockwise; 0 Straight line
        For closed line the first and last vertice bear the same value
        For open line the first and last value are None"""
        try:
            return self._vertex_orientation
        except AttributeError:
            self._vertex_orientation = []
            for i in range(1, len(self.coords) - 1):  # '1' and 'cnt-1' to 'forget' first and last vertice
                orient = GenUtil.orientation(self.coords[i-1], self.coords[i], self.coords[i+1])
                self._vertex_orientation.append(orient)
            if self.is_closed:
                # Case of a closed line or polygon; we do not copy the first and lat even if they are the same
                orient = GenUtil.orientation(self.coords[-2], self.coords[0], self.coords[1])
                self._vertex_orientation = [orient] + self._vertex_orientation
            else:
                # Case of an open line; the first and last are None
                orient = None
                self._vertex_orientation = [orient] + self._vertex_orientation + [orient]
            return self._vertex_orientation


    def _remove_colinear_vertex(self):
        """This method remove the colinear verxtex in the line string. Also handles closed line"""
        if len(self.coords) <= 2:
            # Nothing to do with a line with 2 points
            pass
        else:
            # Detect the position of the colinear vertex
            vertex_to_del = [i for i, orient in (enumerate(self.vertex_orientation)) if orient == 0]
            if len(vertex_to_del) >= 1:
                # Delete the colinear vertex
                lst_coords = list(self.coords)
                for i in reversed(vertex_to_del):
                    del(lst_coords[i])
                if vertex_to_del[0] == 0:
                    # When delete the first vertex than we need to recopy the "new first" to the last vertice
                    lst_coords = lst_coords + [lst_coords[0]]
                self.coords = lst_coords


    def _rotate_start_bend(self):
        """Rotate a closed line string so the start of the line is also the start of a clockwise bend

        To be done on closed line only"""

        i = 0
        rotate = None
        max = len(self.vertex_orientation)
        for i in range(max):
            j = (i+1)%max
            if self.vertex_orientation[i] == GenUtil.CLOCKWISE and \
               self.vertex_orientation[j] == GenUtil.ANTI_CLOCKWISE:
                rotate = i
                break

        # Rotate the frist last vertex to the position of the biggest bend
        if rotate is None:
            # All the bend are clockwise.  Nothing to do
            pass
        elif rotate == 0:
            # The line string does not to be rotated
            pass
        else:
            lst_coord = self.coords[rotate:] + self.coords[1:rotate+1]
            self.coords = lst_coord # Update the LineString coordinate

    def _extract_coords(self, i,j):
        """Extract the coordinate between index [i,j]

        Return
            List of x,y coordinates

        If j is lower than i act like a circular array and avoid duplication of first/last vertice"""

        if i <= j:
            lst_coords = self.coords[i:j+1]
        else:
            lst_coords = self.coords[i:] + self.coords[0:j+1]

        return lst_coords

    def extract_sub_line_string(self, i, j):
        """Extract a sub line string

        Return
            Sub Linestring"""

        lst_coords = self._extract_coords(i, j)
        if len(lst_coords) >=2:
            line = LineString(lst_coords)
        else:
            # create an empty line
            line = LineString(None)

        return line


    def _change_inflexion(self, i,j):
        """Flag if there is an inflexion between the two vertice specidfied.

        There is inflexion when a change of orientation occurs from clock wise to anti clocwise or vice cersa"""

        max = len(self.vertex_orientation)
        if (self.vertex_orientation[i] == GenUtil.ANTI_CLOCKWISE and \
            self.vertex_orientation[(i+1)%max] == GenUtil.CLOCKWISE) or \
           (self.vertex_orientation[i] == GenUtil.CLOCKWISE and \
            self.vertex_orientation[(i+1)%max] == GenUtil.ANTI_CLOCKWISE):
            inflexion = True
        else:
            inflexion = False

        return inflexion

    def _add_bends(self, inflexions):
        """Add then to the line from the inflexion list"""

        for k in range(len(inflexions) - 1):
            i = inflexions[k][0]
            j = inflexions[k + 1][1]
            self.sb_bends.append(Bend(i, j, self._extract_coords(i, j)))


    def _create_bends(self):
        """Create the bends in the line"""

        #Delete any actual bend information
        self.sb_bends = []

        # Remove the colinear vertice in order to facilitate bend detection (moreover colinaer vertice are useless)
        self._remove_colinear_vertex()

        inflexions = []
        max = len(self.vertex_orientation)
        if self.is_closed:
            # Rotate the line to position at the start of a bend
            self._rotate_start_bend()
            # The vertex_oriention list is considered a circular list
            for i in range(max):
                j = (i + 1) % max
                if self._change_inflexion(i,j):
                    inflexions.append((i,j))
            # Create the bend from the inflexion point
            if (inflexions):
                if len(inflexions) >= 3:
                    # If there is more than 23 inflexions we add another circular inflexion
                    i = inflexions[-1][0]
                    j = inflexions[0][1]
                    inflexions.append((i, j))
                # Tranform the inflexion into bends
                self._add_bends(inflexions)

        else:
            # The vertex_oriention list is not considered a circular list
            if max == 3:
                # Special case there is only one bend to simplify
                j = len(self.coords)-1
                self.sb_bends.append(Bend(0, j, self._extract_coords(0, j)))
            elif max >=4:
                for i in range(1, max-2):
                    if self._change_inflexion(i, i+1):
                        inflexions.append((i, i+1))
                # Add inflexion to add the first and last bend
                inflexions = [(0, None)] + inflexions + [(None, max-1)]
                # Transform inflexion into bends
                self._add_bends(inflexions)

        return


    def _sort_bends(self, diameter):
        """Sort the bends from by order of min_adj_are"""

        lst_bends = []
        for i,bend in enumerate(self.sb_bends):
            if bend.adj_area <= self.sb_min_adj_area:
                # Only select the bend below the minimum adjusted area
                lst_bends.append((i, bend.adj_area))

        # Sort based of the adj_area from smalles to biggest
        lst_bends.sort(key=lambda tup: tup[1])  # sorts in place

        return lst_bends

    def _offset_bend_ij(self,i, j):
        """"Offset the value of the different bend i,j because some of the vertice of the line were removed"""

        if i < j:
            offset = j-i-1
        else:
            offset = j
        for bend in self.sb_bends:
            if bend.status == _NOT_SIMPLIFIED:
                if (bend.i < bend.j):
                    if bend.i >= j:
                        bend.i -= offset
                        bend.j -= offset
                else:
                    if bend.i >= j:
                        bend.i -= offset


    def _make_line_ccw(self):
        """Make sure the line is counter clockwise.

        Note: Only apply to closed line"""

        if self.sb_is_closed:
            tmp_ring = LinearRing(self.coords)
            if not tmp_ring.is_ccw:
                #The linear ring is clockwise. Reverse the coordinates to make it ccw
                self.coords = list(reversed(self.coords))


    def simplify(self, diameter, s_constraints=None):
        """Simplify the line by reducing each bend"""

        ray = diameter / 2.0
        self.sb_min_adj_area = _AREA_CMP_INDEX * math.pi * ray ** 2.0
        nbr_bend_simplified = 0

        # Make sure the line is counter clockwise
        #
        self._make_line_ccw()

        #Create the bend in the line
        self._create_bends()

        max_bends = len(self.sb_bends)
        sorted_bends = self._sort_bends(diameter)
        if len(sorted_bends) == 0:
            # No more bend to simplify.  Line is at its simplest form
            self.sb_is_simplest = True
        elif len(sorted_bends)>= 2:
            # Make the biggest bend (last one) unsimplifiable
            ind_last = sorted_bends[-1][0]
            self.sb_bends[ind_last].status = _UNSIMPLIFIABLE

        # Loop over each bend to simplify them
        for sorted_bend in sorted_bends:
            ind = sorted_bend[0]
            if self.sb_bends[ind].status == _NOT_SIMPLIFIED:
                ind_before = None
                ind_after = None
                if self.sb_is_closed:
                    if (max_bends >=2):
                        ind_before = (ind-1)%max_bends
                        ind_after = (ind+1)% max_bends
                else:
                   if ind > 0:
                        ind_before = ind-1
                   if ind < max_bends-1:
                        ind_after = ind+1

                # Validate the spatial constraints
                i = self.sb_bends[ind].i
                j = self.sb_bends[ind].j

                if i < j:
                    lst_coords = self.coords[0:i+1] + self.coords[j:]
                else:
                    # Manage circular list
                    lst_coords = self.coords[j:i+1] + self.coords[j:j+1]

                if s_constraints is not None:
                    in_conflict = s_constraints.check_constraints(self, self.sb_bends[ind])
                else:
                    in_conflict = False

                if not in_conflict:
                    # Update the coordinates
                    self.coords = lst_coords

                    # Bend before and after must no be simplified in this pass maybe a next pass
                    if ind_before is not None: self.sb_bends[ind_before].status = _UNSIMPLIFIABLE
                    if ind_after is not None: self.sb_bends[ind_after].status = _UNSIMPLIFIABLE

                    self.sb_bends[ind].status = _SIMPLIFIED
                    nbr_bend_simplified += 1

                    self._offset_bend_ij(i, j)

        return nbr_bend_simplified


class PointSb(Point):

    def __init__(self, coords, fast_access=True):
        super().__init__(coords)
        self.fast_access = fast_access
        if self.fast_access:
            self.__lst_coords = list(super().coords)

    @property
    def coords(self):
        if self.fast_access:
            return self.__lst_coords
        else:
            return super().coords

    @coords.setter
    def coords(self, coords):
        print ("Need to update the spatial container...")
        Point.coords.__set__(self, coords)
        if self.fast_access:
            self.__lst_coords = list(super().coords)


class SpatialConstraints(object):
    """
    """

    def __init__(self, simplicity=True, crossing=True, sidedness=True, s_container=None):
        """Constructor for the SpatialConstraint class"""

        self.simplicity = simplicity
        self.crossing = crossing
        self.sidedness = sidedness
        self.s_container = s_container
        self.nbr_err_simplicity = 0
        self.nbr_err_crossing = 0
        self.nbr_err_sidedness = 0

    def check_constraints(self, line, bend):
        """Validate the spatial constraint"""

        in_conflict = False

        if not in_conflict:
            in_conflict = self._check_simplicity(line, bend.i, bend.j, bend.replacement_line)

        if not in_conflict:
            in_conflict = self._check_crossing(line, bend.replacement_line)

        if not in_conflict:
            in_conflict = self._check_sidedness(line, bend.polygon)

        return in_conflict


    def _check_simplicity(self, line, i, j, new_sub_line):
        """Check if the new sub line creates a self intersection in the line

        Parameter:
            line -- LineString to verify for self intersection
            i -- Start index of the new sub line in the line
            j -- End index of the new sub line in the line
            new_sub_line -- Replacement line string

        Return
            Boolean True if the line is simple or False otherwise

        """

        # Create a very short line so that the line does not -touch the start and end line (increase performance)
        smaller_sub_line = affinity.scale(new_sub_line, xfact=1. - GenUtil.ZERO, yfact=1. - GenUtil.ZERO)

        # Creates two lines before and after the new_sub_line
        start_line = line.extract_sub_line_string(0, i)
        end_line = line.extract_sub_line_string(j, len(line.coords)-1)

        in_conflict = smaller_sub_line.crosses(start_line) or smaller_sub_line.crosses(end_line)
        if in_conflict:
            self.nbr_err_simplicity += 1

        return in_conflict


    def _check_crossing(self, line, new_sub_line):

        features = self.s_container.get_features(line.bounds, remove_features=[line._sb_sc_id])
        # Check that the new middle line does not cross any interior holes of the polygon
        gen_crosses = filter(new_sub_line.intersects, features)  # Creates a generator
        in_conflict = False
        for element in gen_crosses:
            in_conflict = True
            self.nbr_err_crossing += 1
            break

        return in_conflict

    def _check_sidedness(self, line, pol):

        features = self.s_container.get_features(pol.bounds, remove_features=[line._sb_sc_id])
        # Check that the new middle line does not cross any interior holes of the polygon
        gen_contains = filter(pol.contains, features)  # Creates a generator
        in_conflict = False
        for element in gen_contains:
            in_conflict = True
            self.nbr_err_sidedness += 1
            break

        return in_conflict


    def validateIntersection(self, sp_container, line, bend):
        """Check that the replacement line does not intersect against another line

        Keyword definition
           sp_container -- SpatialContainer containing all the features
           replacement_line -- LineString used to replace the bend
           line -- Linestring line to validate

        Return
            Boolean indicating if the line pass(True) or failed(False) the validation
        """

        lst_lines = sp_container.get_features(bounds=bend.replacement_line.bounds,
                                              remove_features=[line._sb_sc_id],
                                              filter=lambda feature: feature._gbt_geom_type == 'LineString')

        lst_line = list(lst_lines)
        if len(lst_line) !=0:
            print ("Le nombre de lignes est:", len(lst_line))

        gen_crosses = filter(bend.replacement_line.intersects, lst_lines)  # Creates a generator
        intersection = False
        for element in gen_crosses:
            intersection = True
            break # One crossing is enough... we can stop

        return intersection


class Bend(object):
    """
    This class defines attributes and operations for bends

    Attributes: None
    """

    def __init__(self, i, j, bend_coords):
        """
        Initialize Bend object


        """

        self.i = i  # Index of the start of the bend coordinate
        self.j = j  #  Index of the end of the bend coordinate
        self.status = _NOT_SIMPLIFIED  # Type of bend by default: UNTOUCHED
        self.bend_coords = bend_coords  # List of the coordinate forming the bend
        if len(bend_coords) <= 1:
            k = 0


    @property
    def polygon(self):  # Polygon formed by the bendbbbbbbb
        try:
            return self._polygon
        except AttributeError:
            self._polygon = Polygon(self.bend_coords)
            return self._polygon


    @property
    def area(self):  # Area formed by the bend
        try:
            return self._area
        except AttributeError:
            self._area = self.polygon.area
            if self._area <= GenUtil.ZERO: self._area = GenUtil.ZERO  # In case of area=0 we assume almost 0 area instead
            return self._area


    @property
    def base(self):  # The length of the base of the bend
        try:
            return self._base
        except AttributeError:
            self._base = GenUtil.distance(self.bend_coords[0], self.bend_coords[-1])
            if self._base <= GenUtil.ZERO: self._base = GenUtil.ZERO # Avois a case of division by zero
            return self._base


    @property
    def perimeter(self):  # The length of the base of the bend
        try:
            return self._perimeter
        except AttributeError:
            self._perimeter = self.polygon.length
            return self._perimeter


    @property
    def cmp_index(self):  # The compactness index of the bend
        try:
            return self._cmp_index
        except AttributeError:
            self._cmp_index = GenUtil.calculate_compactness_index(self.area, self.perimeter)
            return self._cmp_index


    @property
    def adj_area(self):  # The adjusted area of the bend
        try:
            return self._adj_area
        except AttributeError:
            self._adj_area = GenUtil.calculate_adjusted_area(self.area, self.cmp_index)
            return self._adj_area

    @property
    def replacement_line(self):  # The adjusted area of the bend
        try:
            return self._replacement_line
        except AttributeError:
            self._replacement_line = LineStringSb((self.bend_coords[0], self.bend_coords[-1]))
            return self._replacement_line

    def create_replacement_line (lst_coords, bend, diameter):
        """Calculate the replacement line for a bend"""

        # Extract the sub line containing the bend with one extra vertice on each side
        sub_line = LineStringSb(lst_coords[bend.i-1:bend.j+1])
        bend_i = 1
        bend_j = len(bend.j)-1

        # Translate to sub line so that the bend starts at 0,0
        xoff, yoff = lst_coords[bend.i][0], lst_coords[bend.i][1]
        line_translate = affinity.affine_transform(sub_line, [1, 0, 0, 1, -xoff, -yoff])

        # Extract the angle between the base of the bend (bendi, bendj) and the x axis
        lst_coord = list(line_translate.coords)
        p0 = (lst_coord[bend_j][0], lst_coord[bend_j][1])
        p1 = (lst_coord[bend_i][0], lst_coord[bend_i][1])
        p2 = (abs(p0[0])+1., 0)
        angle = GenUtil.angle_vecor(p0, p1, p2)
#        p0_x = line1_coord[bend_j][0]
#        p0_y = line1_coord[bend_j][1]
#        p1_x = abs(p0_x) + 1.  # In case x == 0
#        p1_y = 0.

#        dot = p0_x * p1_x + p0_y * p1_y
#        len_a = (p0_x ** 2 + p0_y ** 2) ** .5
#        len_b = (p1_x ** 2 + p1_y ** 2) ** .5

        angle = math.acos(dot / (len_a * len_b))
        angle = (angle * 180 / math.pi)

        if p0[1] >= 0.:
            angle = -angle  # Clockwise rotation
#        if p0_y >= 0.:
#            angle = -angle

        # Rotate the bend so it's on the x axis
        a = math.cos(angle)
        b = -math.sin(angle)
        d = math.sin(angle)
        e = math.cos(angle)
        line_rotate = affinity.rotate(line_translate, angle, origin=(0, 0))
        lst_coords = list(line_rotate.coords)


#        line_i = LineString(lst_coords[0:3])
#        line_j = LineString(lst_coords[-2:])
        # Calculate the angle between the base of the bend of segment before and after the bend
        theta_i = lib_geobato.GenUtil.compute_angle(lst_coords[0], lst_coords[1], lst_coords[bend_j])
        theta_j = lib_geobato.GenUtil.compute_angle(lst_coords[bend_j], lst_coords[-2], lst_coords[-1])

        # Determine if the
        bend_line = LineString(lst_coord[bend_i:bend_j+1])
        (minx, miny, maxx, maxy) = bend_line.bounds
        y_dynamic = (abs(miny) + abs(maxy)) * 10.
        x_middle = (lst_coords[bend_i][0] + lst_coords[bend_j][0]) / 2.
        line_y_positive = LineString(((x_middle, 0), (x_middle, y_dynamic)))
        line_y_negative = LineString(((x_middle, 0), (x_middle, -y_dynamic)))
        if line4.crosses(line_y_positive):
            bend_side = +1
        else:
            if line4.crosses(line_y_negative):
                bend_side = -1

        if lst_coords[0][1] >= 0.:
            start_line_side = 1
        else:
            start_line_side = -1

        if lst_coords[-1][1] >= 0.:
            end_line_side = 1
        else:
            end_line_side = -1

        if (start_line_side * end_line_side == -1):
            print("Nothing to do....")
            line5 = LineString(lst_coords[0:bend_i + 1] + lst_coords[bend_j:])
        else:
            # Both line are on the same side
            if start_line_side == 1 and end_line_side == 1:
                if bend_side == -1:
                    angle_bias = 2.
                    y_offset = -1
                else:
                    angle_bias = 3.
                    y_offset = 1
            if start_line_side == -1 and end_line_side == -1:
                if bend_side == 1:
                    angle_bias = 2.
                    y_offset = 1
                else:
                    angle_bias = 3.
                    y_offset = 1

            theta_i = (180. - theta_i) / angle_bias
            if theta_i >= 5.:
                hypothenus = x_middle / math.cos(theta_i * math.pi / 180.)
                y_height = math.sqrt(hypothenus ** 2 - x_middle ** 2)
                if bend_side == -1:
                    y_height *= y_offset
                new_coord = (x_middle, y_height)
                line5 = LineString(lst_coords[0:bend_i + 1] + [new_coord] + lst_coords[bend_j:])
            else:
                print("Nothing to do....")
                line5 = LineString(lst_coords[0:bend_i + 1] + lst_coords[bend_j:])


class AlgoSherbend(object):
    """Main class for the Sherbend algorithm
    
    Attributes:
        - None
        
    """

    def __init__(self, command, geo_content):
        """Initialize the attributes of an object of the class DPAlgorithm

        Keyword:
            command: dataclass containing all the commands for the sherbend line reduction algorithm
            geo_content: dataclass containing the geo information needed for the the sherbend line reduction algorithm

        Return value:
            None
        """

#        Algorithm.__init__(self)
        
        self.command = command
        self.geo_content = geo_content
        ray = command.diameter / 2.0
        self.min_adj_area = _AREA_CMP_INDEX * math.pi * ray ** 2.0
        self.nbr_bend_simplified = 0

    def _calculate_adj_area(self, exclude, coords):

        if (exclude):
            pol = Polygon(coords)
            cmp_index = GenUtil.calculate_compactness_index(pol.area, pol.length)
            adj_area = GenUtil.calculate_adjusted_area(pol.area, cmp_index)
        else:
            adj_area = sys.float_info.max

        return adj_area


    def load_features(self, geo_content, command):
        """Load the points, line strings and polygons in the spatial container.

        The Polygons are deconstructued into a list LineString with clockwise orientation and extra added information
        needed for the reconstruction of the original Polygon

        Args:
            geo_content (dataClass): Contains all the input#output geo spatial information
            command (object): Contains the parameters of the command line interface

        Return
            None
        """

        # Create the spatial container that will receive all the spatial features
        self.s_container = SpatialContainer()

        # Load all the features in the spatial container
        for feature in geo_content.in_features:
            if feature.geom_type == GenUtil.POLYGON:
                adj_area = self._calculate_adj_area(command.exclude_polygon, feature.exterior.coords)
                if (adj_area > self.min_adj_area):
                    # Only keep the polygon over the minimum adjusted area
                    # Deconstruct the Polygon into a list of LineString with supplementary information
                    # needed to reconstruct the original Polygon
                    ext_feature = LineStringSb(feature.exterior.coords)
                    interiors = feature.interiors
                    int_features = []
                    # Extract the interiors as LineString
                    for interior in interiors:
                        adj_area = self._calculate_adj_area(command.exclude_hole, interior.coords)
                        if (adj_area > self.min_adj_area):
                            # Keep the interior over the minimal adjusted area
                            interior = LineStringSb(interior.coords)  # Transform to LineString
                            interior.sb_original_type = GenUtil.POLYGON_INTERIOR
                            int_features.append(interior)
                        else:
                            geo_content.nbr_del_holes += len(feature.interiors)

                    # Add attributes needed for reconstruction
                    ext_feature.sb_interiors = int_features
                    ext_feature.sb_layer_name = feature.sb_layer_name
                    ext_feature.sb_properties = feature.sb_properties
                    ext_feature.sb_original_type = GenUtil.POLYGON_EXTERIOR

                    # Add the exterior and the interior independently
                    self.s_container.add_feature(ext_feature)  # Add the exterior
                    self.s_container.add_features(int_features)  # Add the interiorS
                else:
                    # Do not add the feature (exterior and interiors ) in the spatial container
                    # Update some stats
                    geo_content.nbr_del_polygons += 1
                    geo_content.nbr_del_holes += len(feature.interiors)
            else:  # Geometry is Point or LinseString
                feature = LineStringSb(feature.coords)
                feature.sb_geom_type = feature.geom_type  # For performance to avoid the C caller overhead
                feature.sb_original_type = feature.sb_geom_type

                self.s_container.add_feature(feature)  # Add the feature

        return

    def _manage_lines_simplification (self, s_constraints):
        """Main routine to simplify the lines
        
        For each line to simplify 
            For each valid bend to simplify
                check the consraints if the constraint are violated check alternative bends (only if the 
                number of bend to simplify is one.
                
        One of the costly operation specially for very long line string (like contour) is to rewrite the 
        coordinates into the Shapely structure.  This is why we updtade the shapely structure at the end 
        when the last bend of the line is processed
        
        Parameters: None
            
        Return value
                int: Total number of bend simplified
        """

        iter_nbr = 0
        total_nbr_bend_simplified = 0
        # Iterate until all the line are simplified or there are no more line have to be simplified
        while (True):
            iter_nbr_bend_simplified = 0
            print('Iteration # {}'.format(iter_nbr))
            for line in self.s_container.get_features(filter=lambda feature: feature.sb_geo_type == 'LineString' and not feature.sb_is_simplest):
                nbr_bend_simplified = line.simplify(self.command.diameter, s_constraints)
                iter_nbr_bend_simplified += nbr_bend_simplified
                total_nbr_bend_simplified += nbr_bend_simplified
            print('Number of bend simplified {}'.format(iter_nbr_bend_simplified))
            print('----------')
            iter_nbr += 1
            if iter_nbr_bend_simplified == 0:
                break

        print('Total number of bend simplified: {}'.format(total_nbr_bend_simplified))
        print ('Total number of simplicity error: {}'.format(s_constraints.nbr_err_simplicity))
        print('Total number of crossing error: {}'.format(s_constraints.nbr_err_crossing))
        print('Total number of sidedness error: {}'.format(s_constraints.nbr_err_sidedness))


        return total_nbr_bend_simplified


    def process(self):
        """Main routine for the Sherbend algorithm
        
        The algorithm will simplify the lines using the Sherbend algorithm. 
        It will iterate over the lines until there are no more bends to simplify.

        Keyword arguments:
            none

        Return value:
            geo_content: dataclass containing the output information

        """

        # Load the features into the spatial container
        self.load_features(self.geo_content, self.command)

        s_constraints = SpatialConstraints(s_container=self.s_container)

        self._manage_lines_simplification(s_constraints)

        for feature in self.s_container.get_features():
            if feature.sb_geom_type == GenUtil.POINT:
                self.geo_content.out_features.append(feature)
            elif feature.sb_geom_type == GenUtil.LINE_STRING:
                if feature.sb_original_type == GenUtil.LINE_STRING:
                    self.out_features.append(feature)
                else:
                    if feature.sb_original_type == GenUtil.POLYGON_EXTERIOR:
                        # The LineString was an exterior Polygon so reconstruct the originalPolygon
                        interiors = [list(interior.coords) for interior in feature.sb_interiors]
                        polygon = Polygon(feature.coords, interiors)
                        polygon.sb_layer_name = feature.sb_layer_name
                        polygon.sb_properties = feature.sb_properties
                        self.geo_content.out_features.append(polygon)
                    else:
                        pass  # Nothing to do with the holes here

        return
