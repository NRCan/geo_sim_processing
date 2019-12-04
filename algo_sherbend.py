#!/usr/local/bin/python
# -=- encoding: utf-8 -=-



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


import math, sys

from shapely.geometry import Point, LineString, LinearRing, Polygon
from shapely.prepared import prep
from shapely import affinity
from lib_geosim import GenUtil, PointSc, LineStringSc, SpatialContainer, GeoSimException

# Internal constant ===> Should be modify with care...
_AREA_CMP_INDEX = .75  # Compactness index factor applied to the adjusted area

#Internal key word constants
_BURNED = "Burned"
_DIAMETER = "diameter"
_SIMPLIFIED = 'Simplified'
_NOT_SIMPLIFIED = 'NotSimplified'
_UNSIMPLIFIABLE = 'Unsimplifiable'


class LineStringSb(LineStringSc):
    """
    A class to represent a LineString used by the SherBend algorithm

    Attributes
    ----------
    coords : List
        A list of coordinates (x,y)
    original_type: str
        The original type of the feature
    min_adj_are : float
        The minimal adjusted area below which the vends are deleted
    properties : dict
        The dictionary of the properties (attributes of the features)
    fast_access : Boolean
        A flag to indicate if we keep a copy od the coordinate in order to accelrate the access becase
        the access to the C function is slow

    """

    def __init__(self, coords, original_type, min_adj_area, layer_name, properties, fast_access=True):
        super().__init__(coords)
        self.sb_original_type = original_type
        self.sb_layer_name = layer_name
        self.sb_properties = properties
        self.sb_min_adj_area = min_adj_area
        self._sb_fast_access = fast_access
        if self._sb_fast_access:
            self.__lst_coords = list(super().coords)

        # Declaration of the instance variable
        self.sb_geom_type = self.geom_type  # variable defined to avoid slower C calls with geom_type
        self.sb_is_simplest = False  # The line is not at its simplest form
        self.sb_bends = []  # Holder for the bend of the line

    # Is the line string closed
    @property
    def sb_is_closed(self):
        try:
            return self._sb_is_closed
        except AttributeError:
            # A closed line need at least 4 vertex to be valid
            if len(self.coords) >= 4 and GenUtil.distance(self.coords[0], self.coords[-1]) <= GenUtil.ZERO:
                self._sb_is_closed = True
            else:
                self._sb_is_closed = False
            return self._sb_is_closed

    @property
    def coords(self):
        if self._sb_fast_access:
            return self.__lst_coords
        else:
            return super().coords

    @coords.setter
    def coords(self, coords):
        # Access the coord attribute in the parent class
        super(LineStringSb, self.__class__).coords.fset(self, coords)  # Odd writing but it's needed...
        if self._sb_fast_access:
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
        """This method remove the co linear vertex in the line string. Also handles closed line

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        if len(self.coords) <= 2:
            # Nothing to do with a line with 2 points
            pass
        else:
            # Detect the position of the colinear vertex
            vertex_to_del = [i for i, orient in (enumerate(self.vertex_orientation)) if orient == 0]
            if len(vertex_to_del) >= 1:
                # Delete the co linear vertex
                lst_coords = list(self.coords)
                for i in reversed(vertex_to_del):
                    del(lst_coords[i])
                if vertex_to_del[0] == 0:
                    # When delete the first vertex than we need to recopy the "new first" to the last vertice
                    lst_coords = lst_coords + [lst_coords[0]]
                self.coords = lst_coords

    def _rotate_start_bend(self):
        """Rotate a closed line string so the start of the line is also the start of a clockwise bend

        To be done on closed line only

        Parameters
        ----------
        None

        Return
        ------
        None

        """

        rotate = None
        max_v = len(self.vertex_orientation)
        for i in range(max_v):
            j = (i+1) % max_v
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

        If j is lower than i act like a circular array and avoid duplication of first/last vertice

        Parameters
        ----------
        i,j : int
            Index used to extract a sub list

        Return
        ------
        List
            list of (x,y) coordinates

        """

        if i <= j:
            lst_coords = self.coords[i:j+1]
        else:
            lst_coords = self.coords[i:] + self.coords[0:j+1]

        return lst_coords

    def _change_inflexion(self, i):
        """Flag if there is an inflexion between at the specified vertices.

        There is inflexion when a change of orientation occurs from clock wise to anti clocwise or vice cersa

        Parameters
        ----------
        i : int
            Index of for vertex orientation

        Return
        ------
        bool
            Flag indicating if an inflexion occurs or not
        """

        max_v = len(self.vertex_orientation)
        if (self.vertex_orientation[i] == GenUtil.ANTI_CLOCKWISE and
            self.vertex_orientation[(i+1) % max_v] == GenUtil.CLOCKWISE) or \
           (self.vertex_orientation[i] == GenUtil.CLOCKWISE and
            self.vertex_orientation[(i+1) % max_v] == GenUtil.ANTI_CLOCKWISE):
            inflexion = True
        else:
            inflexion = False

        return inflexion

    def _add_bends(self, inflexions):
        """Add Bend to the line from the inflexion list

        Parameters
        ----------
        inflexions : List
            List of the inflexions in the list

        Return
        ------
        None
        """

        for k in range(len(inflexions) - 1):
            i = inflexions[k][0]
            j = inflexions[k + 1][1]
            self.sb_bends.append(Bend(i, j, self._extract_coords(i, j)))

    def _create_bends(self):
        """Create the bends in the line

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        # Delete any actual bend information
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
                if self._change_inflexion(i):
                    inflexions.append((i, j))
            # Create the bend from the inflexion point
            if inflexions:
                if len(inflexions) >= 3:
                    # If there is more than 23 inflexions we add another circular inflexion
                    i = inflexions[-1][0]
                    j = inflexions[0][1]
                    inflexions.append((i, j))
                # Transform the inflexion into bends
                self._add_bends(inflexions)

        else:
            # The vertex_oriention list is not considered a circular list
            if max == 3:
                # Special case there is only one bend to simplify
                j = len(self.coords)-1
                self.sb_bends.append(Bend(0, j, self._extract_coords(0, j)))
            elif max >= 4:
                for i in range(1, max-2):
                    if self._change_inflexion(i):
                        inflexions.append((i, i+1))
                # Add inflexion to add the first and last bend
                inflexions = [(0, None)] + inflexions + [(None, max-1)]
                # Transform inflexion into bends
                self._add_bends(inflexions)

        return

    def _sort_bends(self):
        """Sort the bends by order of ascending min_adj_are

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        lst_bends = []
        for i, bend in enumerate(self.sb_bends):
            if bend.adj_area <= self.sb_min_adj_area:
                # Only select the bend below the minimum adjusted area
                lst_bends.append((i, bend.adj_area))

        # Sort based of the adj_area from smallest to biggest
        lst_bends.sort(key=lambda tup: tup[1])  # sorts in place

        return lst_bends

    def _offset_bend_ij(self, i, j):
        """"Offset the value of the different bend i,j because one or more vertice of the line were removed

        Handle circular list when j < i

        Parameters
        ----------
        i,j : int
            Index in the line where the vertice were removed

        Return
        ------
        None
        """

        if i < j:
            offset = j-i-1
        else:
            offset = j
        for bend in self.sb_bends:
            if bend.status == _NOT_SIMPLIFIED:
                if bend.i < bend.j:
                    if bend.i >= j:
                        bend.i -= offset
                        bend.j -= offset
                else:
                    if bend.i >= j:
                        bend.i -= offset

    def _make_line_ccw(self):
        """Make sure the line is counter clockwise.

        Only apply to closed line

        Parameters
        ----------
        None

        Return
        ------
        None"""

        if self.sb_is_closed:
            tmp_ring = LinearRing(self.coords)
            if not tmp_ring.is_ccw:
                # The linear ring is clockwise. Reverse the coordinates to make it ccw
                self.coords = list(reversed(self.coords))

    def simplify(self, diameter, s_constraints=None):
        """Simplify the line by reducing each bend

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        nbr_bend_simplified = 0

        # Make sure the line is counter clockwise
        #
        self._make_line_ccw()

        # Create the bend in the line
        self._create_bends()

        max_bends = len(self.sb_bends)
        sorted_bends = self._sort_bends()
        if len(sorted_bends) == 0:
            # No more bend to simplify.  Line is at its simplest form
            self.sb_is_simplest = True
        elif len(sorted_bends) >= 2:
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
                    if max_bends >= 2:
                        ind_before = (ind-1) % max_bends
                        ind_after = (ind+1) % max_bends
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

                if self.is_closed:
                    if len(lst_coords) >= 4:
                        if s_constraints is not None:
                            in_conflict = s_constraints.check_constraints(self, self.sb_bends[ind])
                        else:
                            in_conflict = False
                    else:
                        # A closed line cannot have less than 4 vertices
                        in_conflict = True
                else:
                    if len(lst_coords) >= 2:
                        if s_constraints is not None:
                            in_conflict = s_constraints.check_constraints(self, self.sb_bends[ind])
                        else:
                            in_conflict = False
                    else:
                        # An open line cannot have less than  3 vertices
                        in_conflict = True

                if not in_conflict:
                    # Update the coordinates
                    self.coords = lst_coords

                    # Bend before and after must no be simplified in this pass maybe a next pass
                    if ind_before is not None:
                        self.sb_bends[ind_before].status = _UNSIMPLIFIABLE
                    if ind_after is not None:
                        self.sb_bends[ind_after].status = _UNSIMPLIFIABLE

                    self.sb_bends[ind].status = _SIMPLIFIED
                    nbr_bend_simplified += 1

                    self._offset_bend_ij(i, j)

        return nbr_bend_simplified


class PointSb(PointSc):
    """
    A class to represent a Point used by the SherBend algorithm

    Attributes
    ----------
    coords : tuple
        A tuple (x,y) representing one coordinate
    properties : dict
        The dictionary of the properties (attributes of the features)
    fast_access : Boolean
        A flag to indicate if we keep a copy od the coordinate in order to accelrate the access becase
        the access to the C function is slow
    """

    def __init__(self, coords, layer_name, properties, fast_access=True):
        super().__init__(coords)
        self.sb_is_simplest = True
        self.sb_layer_name = layer_name
        self.sb_properties = properties
        self.sb_original_type = GenUtil.POINT
        self.sb_geom_type = GenUtil.POINT  # For faster access than calling C (geom_type)
        self._sb_fast_access = fast_access
        if self._sb_fast_access:
            self.__lst_coords = list(super().coords)

    @property
    def coords(self):
        if self._sb_fast_access:
            return self.__lst_coords
        else:
            return super().coords

    @coords.setter
    def coords(self, coords):
        Point.coords.__set__(self, coords)
        if self._sb_fast_access:
            self.__lst_coords = list(super().coords)


class SpatialConstraints(object):
    """
    A class to represent validation of spatial constraints

    Attributes
    ----------
    simplicity : bool
        Flag indicating if simplicity constraint (self crossing) is validated
    crossing : bool
        Flag indicating if crossing constraint (intersection between feature) is validated
    sidedness : bool
        Flag indicating if sidedness constraint (relative adjacency) is validated
    s_container : SpatialContainer
        Object containing all the feature
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
        """Validate the different spatial constraint

        Parameters
        ----------
        line : LineStringSb
            LineString to validate for spatial constraints
        bend : Bend
            Bend to validate for spatial constraints

        Returns
        -------
            bool
            Flag indicating if the spatial constrainst are valid or not"""

        in_conflict = False

        if not in_conflict:
            in_conflict = self._check_simplicity(line, bend.replacement_line)

        if not in_conflict:
            in_conflict = self._check_crossing(line, bend.replacement_line)

        if not in_conflict:
            in_conflict = self._check_sidedness(line, bend.polygon)

        return in_conflict

    def _check_simplicity(self, line, new_sub_line):
        """Check if the new sub line creates a self intersection in the line

        Parameter
        ---------
        line : LineStringSb
            LineString to validate for self intersection
        new_sub_line : LineString
            New LineString to validate for self intersection

        Returns
        -------
            Boolean
                Flag indicating if the line is simple or not
        """

        # Create a very short line so that the line does not -touch the start and end line (increase performance)
        smaller_sub_line = affinity.scale(new_sub_line, xfact=1. - GenUtil.ZERO, yfact=1. - GenUtil.ZERO)

        in_conflict = False
        prepared_smaller_sub_line = prep(smaller_sub_line)
        if prepared_smaller_sub_line.intersects(line):
            in_conflict = True
            self.nbr_err_simplicity += 1

        return in_conflict

    def _check_crossing(self, line, new_sub_line):
        """Check if the new sub line intersects other line

            Parameter
            ---------
            line : LineStringSb
                LineString to validate for intersection with other line
            new_sub_line : LineString
                New LineString to validate for intersection with other line

            Returns
            -------
                Boolean
                    Flag indicating if the line intersect with other line or not
        """

        features = self.s_container.get_features(new_sub_line.bounds, remove_features=(line,))

        # Check that the new middle line does not cross any interior holes of the polygon
        prepared_new_sub_line = prep(new_sub_line)
        in_conflict = False
        gen_crosses = filter(prepared_new_sub_line.intersects, features)
        for feature in gen_crosses:
            in_conflict = True
            self.nbr_err_crossing += 1
            break

        return in_conflict

    def _check_sidedness(self, line, pol):
        """Validate the line for adjacency constraints

        Parameter
        ---------
        line : LineStringSb
            LineString to validate for adjacency
        new_sub_line : LineString
            New Polygon to check for adjacency

        Returns
        -------
            Boolean
                Flag indicating if the line creates or not adjacency problem
        """

        features = self.s_container.get_features(pol.bounds, remove_features=(line,))
        # Check that the new middle line does not cross any interior holes of the polygon
        prepared_pol = prep(pol)
        gen_contains = filter(prepared_pol.contains, features)
        in_conflict = False
        for feature in gen_contains:
            in_conflict = True
            self.nbr_err_sidedness += 1
            break

        return in_conflict


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
            if self._area <= GenUtil.ZERO:
                self._area = GenUtil.ZERO  # In case of area=0 we assume almost 0 area instead
            return self._area

    @property
    def base(self):  # The length of the base of the bend
        try:
            return self._base
        except AttributeError:
            self._base = GenUtil.distance(self.bend_coords[0], self.bend_coords[-1])
            if self._base <= GenUtil.ZERO:
                self._base = GenUtil.ZERO # Avois a case of division by zero
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
            self._replacement_line = LineString((self.bend_coords[0], self.bend_coords[-1]))
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
        
        self.command = command
        self.geo_content = geo_content
        self.nbr_bend_simplified = 0

    def calculate_min_adj_area(self, diameter):

        return  (_AREA_CMP_INDEX * math.pi * (diameter/2.0)**2.0)

    def _calculate_adj_area(self, coords):

        pol = Polygon(coords)
        cmp_index = GenUtil.calculate_compactness_index(pol.area, pol.length)
        adj_area = GenUtil.calculate_adjusted_area(pol.area, cmp_index)

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
            diameter = command.dlayer_dict[feature.sb_layer_name]
            min_adj_area = self.calculate_min_adj_area(diameter)
            if feature.geom_type == GenUtil.POINT:
                out_feature = PointSb(feature.coords, feature.sb_layer_name, feature.sb_properties)
                # Add the feature
                self.s_container.add_feature(out_feature)
            elif feature.geom_type == GenUtil.LINE_STRING:
                out_feature = out_feature = LineStringSb(feature.coords, GenUtil.LINE_STRING, min_adj_area, feature.sb_layer_name,
                                                         feature.sb_properties)
                # Add the feature
                self.s_container.add_feature(out_feature)
            elif feature.geom_type == GenUtil.POLYGON:
                adj_area = self._calculate_adj_area(feature.exterior.coords)
                # Only keep the polygon over the minimum adjusted area
                if not command.exclude_polygon or adj_area > min_adj_area:
                    # Deconstruct the Polygon into a list of LineString with supplementary information
                    # needed to reconstruct the original Polygon
                    ext_feature = LineStringSb(feature.exterior.coords, GenUtil.POLYGON_EXTERIOR, min_adj_area,
                                               feature.sb_layer_name, feature.sb_properties)
                    interiors = feature.interiors
                    int_features = []
                    # Extract the interiors as LineString
                    for interior in interiors:
                        adj_area = self._calculate_adj_area(interior.coords)
                        # Only keep the interior (hole) over the minimal adjusted area
                        if not command.exclude_hole or adj_area > min_adj_area:
                            interior = LineStringSb(interior.coords, GenUtil.POLYGON_INTERIOR, min_adj_area, None, None)
                            int_features.append(interior)
                        else:
                            geo_content.nbr_del_holes += len(feature.interiors)

                    # Add interior features needed for Polygon reconstruction
                    ext_feature.sb_interiors = int_features

                    # Add the exterior and the interior independently
                    self.s_container.add_feature(ext_feature)  # Add the exterior
                    self.s_container.add_features(int_features)  # Add the interiorS
                else:
                    # Do not add the feature (exterior and interiors ) in the spatial container
                    # Update some stats
                    geo_content.nbr_del_polygons += 1
                    geo_content.nbr_del_holes += len(feature.interiors)
            else:
                raise GeoSimException ("Invalid geometry type: {}".format(feature.geometry))

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

            # Build line iterator
            lines = (feature for feature in self.s_container.get_features()
                                 if(not feature.sb_is_simplest and feature.sb_geom_type==GenUtil.LINE_STRING ))
            for line in lines:
                nbr_bend_simplified = line.simplify(self.command.diameter, s_constraints)
                iter_nbr_bend_simplified += nbr_bend_simplified
                total_nbr_bend_simplified += nbr_bend_simplified
            print('Number of bend simplified {}'.format(iter_nbr_bend_simplified))
            print('----------')
            iter_nbr += 1
            if iter_nbr_bend_simplified == 0:
                break

        print('Total number of bend simplified: {}'.format(total_nbr_bend_simplified))
        print('Total number of simplicity error: {}'.format(s_constraints.nbr_err_simplicity))
        print('Total number of crossing error: {}'.format(s_constraints.nbr_err_crossing))
        print('Total number of sidedness error: {}'.format(s_constraints.nbr_err_sidedness))

        return total_nbr_bend_simplified

    def process(self):
        """Main routine for the Sherbend algorithm
        
        The algorithm will simplify the lines using the Sherbend algorithm. 
        It will iterate over the lines until there are no more bends to simplify.

        Keyword arguments:
            None

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
                    self.geo_content.out_features.append(feature)
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
