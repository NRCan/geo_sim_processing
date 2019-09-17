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

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.geometry.polygon import orient
from shapely import affinity

from lib_geobato import GenUtil, SpatialContainer, PointSc, LineStringSc, Polygon
                               
# Public key word contants
MULTI_BENDS = "MULTI_BEND"
SINGLE_BEND = "SINGLE_BEND"
NO_VERTICE_ADD = "NO_VERTICE_ADD"

# Properties name
_DIAMETER = "diameter"
_SIMPLIFY_FIRST_LAST = "simplify_first_last"

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
_SMALL = "Small"
_ONE_BEND = 'OneBend'
_TWO_BENDS = 'TwoBends'
_THREE_BENDS = 'ThreeBends'
_FOUR_BENDS = 'FourBends'
_FIVE_BENDS = 'FiveBends'
_SIMPLIFIED = 'Simplified'
_UNKNOWN = 'Unknown'
_IN_CONFLICT = 'InConflict'

# class SherbendStatistics(GenStatistics):
#     """Class that contains the statistics for the Sherbend algorithm
#
#     Attributes
#         stat_names: Name of the statistics for the SherbendStatistics class. These name are
#                     used by the Statistics class
#
#     """
#
#     def __init__(self):
#         """Initialize the attributes of an object of the class Sherbend statistics
#
#         Parameter:
#             None
#
#         Return value
#             None
#
#         """
#
#         GenStatistics.__init__(self)
#         self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS,  \
#                              _BIG, _SMALL, _ONE_BEND, _TWO_BENDS, _THREE_BENDS, _FOUR_BENDS, _FIVE_BENDS ))
#
#     def get_stats (self, type=GenStatistics.SUMMARY):
#         """Extract the current statistics and build  a list of string that forms the statistical message"
#
#         Parameters:
#             type: Give the form of statistics to extract. Can take 2 values.
#                 SUMMARY: Summary information
#                 DETAILED: Detailed information
#
#         """
#
#         str_out = []
#         str_out.append( "Sherbend algorithm Statistics" )
#
#         str_out.append( "--------------------------" )
#         if (type == GenStatistics.DETAILED):
#             for i in xrange((self.get_nbr_iteration())):
#                 str_out.append("Detailed statistics")
#                 str_out.append("Iteration # " + str(i))
#                 str_out.append("Bend simplified: " + str(self.get_stats_name_count_iter( _ALGO, i)))
#                 str_out.append( "--------" )
#                 str_out.append( "Conflicts:" )
#                 str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_iter( GenUtil.SIMPLE_LINE, i)))
#                 str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
#                 str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
#                 str_out.append( "--------" )
#         str_out.append( "Summary statistics" )
#         str_out.append("Total bend simplified: " + str(self.get_stats_name_count_total(_ALGO)))
#         str_out.append( "--------" )
#         str_out.append("Statistics by bends")
#         str_out.append("     One bend      :  " +  str(self.get_stats_name_count_total(_ONE_BEND)))
#         str_out.append("     Two bends     :  " +  str(self.get_stats_name_count_total(_TWO_BENDS)))
#         str_out.append("     Three bends   :  " +  str(self.get_stats_name_count_total(_THREE_BENDS)))
#         str_out.append("     Four bends    :  " +  str(self.get_stats_name_count_total(_FOUR_BENDS)))
#         str_out.append("     Five bends    :  " +  str(self.get_stats_name_count_total(_FIVE_BENDS)))
#         str_out.append("     Small bends   :  " +  str(self.get_stats_name_count_total(_SMALL)))
#         str_out.append("     Big bends    :  " +  str(self.get_stats_name_count_total(_BIG)))
#         str_out.append( "--------" )
#         str_out.append( "Conflicts:" )
#         str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_total( GenUtil.SIMPLE_LINE)))
#         str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
#         str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
#         str_out.append( "--------" )
#         str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )
#
#         return str_out


class SpatialConstraints(object):
    """
    """

    _simplicity = []
    _intersection = []

    def validateSimplicity(self, bend, line, start_line, middle_line, end_line):
        """Check if the middle line is crossing the start_line or the end_line

        Keyword definition
           bend -- Bend to simplify
           line -- Line
           start_line -- LineString to check
           middle_line -- LineString to check
           end_line -- LineString to check

        Return
            Boolean indicating if the line pass(True) or failed(False) the validation
        """

        # Create a very short line so that the line does not touch the start and end line (increase performance)
        smaller_middle_line = affinity.scale(middle_line, xfact=1. - GenUtil.ZERO, yfact=1. - GenUtil.ZERO)
        in_conflict = smaller_middle_line.crosses(start_line) or smaller_middle_line.crosses(end_line)

        if not in_conflict:
            if line._gbt_original_type == 'Polygon-Exterior':
                # Manage the interiors of a polygon
                lst_holes = line._gbt_interiors

                # Check that the new middle line does not cross any interior holes of the polygon
                gen_crosses = list(filter(middle_line.intersects, lst_holes))  # Creates a generator
                in_conflict = False
                for element in gen_crosses:
                    in_conflict = True
                    break

                # Check that the bend to simplify does not contain any interior holes of the polygon
                if not in_conflict:
                    gen_inside = list(filter(bend.polygon.contains, lst_holes))
                    in_conflict = False
                    for element in gen_inside:
                        in_conflict = True
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
                                              remove_features=[line._gbt_sc_id],
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
        self.type = _UNKNOWN  # Type of bend by default:UNKNOWN
        self.bend_coords = bend_coords  # List of the coordinate forming the bend


    @property
    def polygon(self):  # Polygon formed by the bend
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
            self._cmp_index = 4*self._area*math.pi / (self.perimeter**2.0)
            return self._cmp_index


    @property
    def adj_area(self):  # The adjusted area of the bend
        try:
            return self._adj_area
        except AttributeError:
            self._adj_area = self.area * (0.75 / self.cmp_index)
            return self._adj_area

    @property
    def replacement_line(self):  # The adjusted area of the bend
        try:
            return self._replacement_line
        except AttributeError:
            self._replacement_line = LineStringSc((self.bend_coords[0], self.bend_coords[-1]))
            return self._replacement_line

    def create_replacement_line (lst_coords, bend, diameter):
        """Calculate the replacement line for a bend"""

        # Extract the sub line containing the bend with one extra vertice on each side
        sub_line = LineStringSc(lst_coords[bend.i-1:bend.j+1])
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

#        # Set the parameters according to the params.bend_mode parameter
#        if self.params.bend_mode == MULTI_BENDS:
#            self.params.multi_bend = True
#            self.params.big_bend = True
#        elif self.params.bend_mode == SINGLE_BEND:
#            self.params.multi_bend = False
#            self.params.big_bend = True
#        else:
#            self.params.multi_bend = False
#            self.params.big_bend = False

    def load_features(self, features):
        """Load the points, line strings and polygons in the spatial container.

        The Polygons are deconstructued into a list LineString with clockwise orientation and extra added information
        needed for the reconstruction of the original Polygon

        Keyword definition
            features: List of shapely features

        Return
            None
        """

        # Create the spatial container that will receive all the spatial features
        self.s_container = SpatialContainer()

        # Load all the features in the spatial container
        for feature in features:
            if feature.geom_type == GenUtil.POLYGON:
####                feature._gbt_geom_type = 'LineString'  # For performance to avoid the C caller overhead
                # Deconstruct the Polygon into a list of LineString with supplementary information
                # needed to reconstruct the original Polygon
                tmp_pol = orient(Polygon(feature.exterior.coords), GenUtil.CLOCKWISE) # Orient vertices clockwiswe
                ext_feature = LineString(tmp_pol.exterior.coords)
                interiors = feature.interiors
                int_features = []
                # Extract the interiors as LineString
                for interior in interiors:
                    tmp_pol = orient(Polygon(interior.coords), GenUtil.CLOCKWISE)  # Orient vertices clockwiswe
                    interior = LineString(tmp_pol.exterior.coords)  # Transform to LineString
                    interior._gbt_geom_type = GenUtil.LINE_STRING  # For performance to avoid the C caller overhead
                    interior._gbt_original_type = GenUtil.POLYGON_INTERIOR
                    int_features.append(interior)

                # Add attributes needed for reconstruction
                ext_feature._gbt_interiors = int_features
                ext_feature._gbt_layer_name = feature._gbt_layer_name
                ext_feature._gbt_properties = feature._gbt_properties
                ext_feature._gbt_geom_type = GenUtil.LINE_STRING  # For performance to avoid the C caller overhead
                ext_feature._gbt_original_type = GenUtil.POLYGON_EXTERIOR

                # Add the exterior and the interior independently
                self.s_container.add_feature(ext_feature)  # Add the exterior
                self.s_container.add_features(int_features)  # Add the interior
            else:  # Geometry is Point or LinseString
                feature._gbt_geom_type = feature.geom_type  # For performance to avoid the C caller overhead
                feature._gbt_original_type = feature._gbt_geom_type

                self.s_container.add_feature(feature)  # Add the feature

        # Empty the feature list
        features.clear()

        return


    def add_line_attributes (self, diameter):
        """This routine sets different attributes of the lines

        Keyword definition
           diameter -- diameter of the bend to simplify
            
        Return value: None

        """
        for line in self.s_container.get_features(filter= lambda feature : feature._gbt_geom_type==GenUtil.LINE_STRING):
            
            line._gbt_is_simplest = False
            line._gbt_bends = []
            # Set the minimal adjusted area
            ray = diameter/2.0
            line._gbt_min_adj_area = _AREA_CMP_INDEX * math.pi * ray**2.0
            
            # Check if the line is a closed line by comparing  first and vertice
            if GenUtil.distance(line.coords[0], line.coords[-1]) <= GenUtil.ZERO:
                line._gbt_is_closed = True
            else:
                line._gbt_is_closed = False
                
        return

#    def classify_bends (self, line):
#        """Classify the beds
#
#        High level routine to classify the bends from _ONE_BEND to _FIVE_BENDS and to calculate
#        the replacement line for the bends that are not of type "_UNKNOWN"
#
#        Parameter:
#            s_container: Object of type SpatialContainer containing the features nd the spatial index
#
#        Return value:
#            None
#
#        """
#
#        if (self.command.multi_bend and line.multi_bend):
#            self._detect_window_bend (line, _FIVE_BENDS)
#            self._detect_window_bend (line, _FOUR_BENDS)
#            self._detect_window_bend (line, _THREE_BENDS)
#            self._detect_window_bend (line, _TWO_BENDS)
#            self._detect_window_bend (line, _ONE_BEND)
#        else:
#            self._detect_window_bend (line, _ONE_BEND)
#        self._detect_window_bend(line, _ONE_BEND)
#
#        # Detect and flag lines at their simplest form
#        self._set_simplest_line(line)
#
#        # Compaction of the bends in the list of bends.  This compactness will speedup the search process
#        self._compact_bends(line)
#
#        # Calculates the replacement line for each bend found
#        self._reduce_bends(line)


    def _set_simplest_line(self, line):
        """ Determine if a line is at its simplest form
        
        Set the attribute of the simplest at true when a line as 0 bends (the bend list is empty) or
        all the bends of the line are of type "UNKNOWN"  this means that all the bends are over 
        the minimum adjusted size
        
        Parameter:
            line: LineString to check if it is at its simplest form
            
        Return value
            None
            
        """
        
        if line._gbt_bends:
            list_unknown = [bend.type for bend in line._gbt_bends if bend.type != _UNKNOWN]
            # If list is empty ==> All the bends are of type _UNKNOWN and the line is at its simplest form
            if list_unknown:
                pass
            else:
                line._gbt_simplest= True
        else:
            line._gbt_simplest = True


    def classify_bends(self, line):
        """This routine is sliding a window over the bends of a line in order to classify bend for simplification.
        
        It tries to find consecutives bends
        that meet the following criterias:
          - The bends must have the bend_type _UNKNOWN
          - The bends must surrounded by bend of type _UNKNOWN
          - The bends must have a adjusted area below min_adj_area
          - The first and last bend of the window must be smaller than the 
               previous and next bend of the window
               
        Parameter:
            line: line object to check for consecutives bend
            
        Return value:
            None
    
        """
        
#        window_length = self._number_of_bends_to_reduce(bend_type)
        window_length = 1
        bend_type = _ONE_BEND
#        if bend_type == _ONE_BEND:
#            window_length = 1
#        else:
#            window_length = 0

        nbr_bends = len(line._gbt_bends)
        last_bend = nbr_bends - window_length 
        
        # The last bend to process depends on the window length
        for i in range(0, last_bend + 1):
            # The sliding window must be surrounded by bend type _UNKNOWN 
            # otherwise we cannot check if they are similar bends
            if (i == 0):
                # Special case for the first bend
                surround_before = _UNKNOWN
            else:
                surround_before = line._gbt_bends[i - 1].type
            if (i == last_bend):
                #How the last bend we assume the next bend to _UNKNOWN
                surround_after = _UNKNOWN
            else:
                surround_after = line._gbt_bends[i + window_length].type
            
            # Check that all the bends over the window are of type UNKNOWN
            bend_unknown = True
            for j in range(window_length):
                if (line._gbt_bends[i + j].type != _UNKNOWN):
                    bend_unknown = False
            
            # Check previous and next bend are of type _UNKNOWN
            if (bend_unknown and
                 surround_before == _UNKNOWN  and
                 surround_after == _UNKNOWN):
                
                # Check that all bends are below min_adj_area
                window_adj_area = True
                for j in range(window_length):
                    if (line._gbt_bends[i + j].adj_area > line._gbt_min_adj_area):
                        window_adj_area = False
                    
                if (window_adj_area):
                    
                    if ( self._are_neighbours_bigger(line, i, i+window_length-1) ):
                        
                        lst_bends = []
                        for i_lst in range(window_length):
                            lst_bends.append(line._gbt_bends[i_lst+i])
                            
                        if ( self._are_bends_similar(lst_bends) ):
                            for i_lst in range(window_length):
                                line._gbt_bends[i_lst+i].type = bend_type
                        
        return
    
    def _are_neighbours_bigger(self, line, first_bend, last_bend):
        """Checks if the bends is surrounded by bigger area bends.
        
        Parameters:
            line: Line object to check
            first_bend: number of the first bend of the line to check
            last_bend: Number of the last bend of the line to check
            
        Return value:
            Boolean flag for the state of the bend in comparison to its neighbours
                True: The neighbours are bigger
                False:  The neighbours are smaller
        """
        
        # If it is the first bend of the line we assume previous bend as infinite
        if (first_bend == 0):
            previous_adj_area = 1.0e+99
        else:
            previous_adj_area = line._gbt_bends[first_bend - 1].adj_area
        
        # If it is the last bend of the line we assume the next bend as infinite    
        if (last_bend == len(line._gbt_bends) - 1):
            next_adj_area = 1.0e+99
        else:
            next_adj_area = line._gbt_bends[last_bend + 1].adj_area
            
        if (previous_adj_area >= line._gbt_bends[first_bend].adj_area and
             next_adj_area >= line._gbt_bends[last_bend].adj_area):
            
            bigger = True
        else:
            bigger = False
            
        return bigger

    
    # def _compact_bends (self, line):
    #     """This routine compacts the bends.
    #
    #     When there are multiple bends we compact them as one bend and with the those bends
    #     in the attribute multi_bends
    #
    #     Parameters:
    #         line: Line object to compact the bends
    #
    #     Return value:
    #         None
    #
    #     Note: We start the compaction by the end to avoid index problem when it remove
    #           number
    #     """
    #
    #     i = len(line._gbt_bends)-1 # i is on the last bend
    #     while (i >= 0):
    #         nbr_bends = self._number_of_bends_to_reduce (line.bends[i].type)
    #         if (nbr_bends >=1):
    #             lst_bends = []
    #             for j in range(nbr_bends):
    #                 lst_bends.append(line._gbt_bends[i-j])
    #             lst_bends.reverse()
    #
    #             i = i - (nbr_bends -1) # Go to the first bend of the multi bend
    #             line.bends[i].multi_bends = lst_bends
    #
    #             # Keep the first bend and delete the remaining bends
    #             for dummy in range(nbr_bends-1):
    #                 del line._gbt_bends[i+1]
    #
    #         i -= 1
    #
    # def _reduce_bends(self, line):
    #     """Detect the bends to simplify in the line
    #
    #     Parameters:
    #         line: Line object to detect bend
    #     """
    #
    #     for bend in line._gbt_bends:
    #         if bend.type == _ONE_BEND:
    #             if (self._are_bends_big(bend, line._gbt_min_adj_area)):
    #                 bend_size = _BIG
    #             else:
    #                 bend_size = _SMALL
    #
    #             bend.replacement_line = (line.coords[bend.i],line.coords[bend.j] )


    def create_bends(self, line):
        """Create the bends in a line"""

        # Determine the inflexion in the line
        lst_ij = GenUtil.locate_bends(line.coords)

        # Create bend attribute list or reset in the case  there are already some bends
        line._gbt_bends = []

        # Create the bends
        for (i, j) in lst_ij:
            line._gbt_bends.append(Bend(i, j, line.coords[i:j + 1]))

        return


    def rotate_coordinates(self, line):
        """Rotate the first and last vertice of a closed line on the place where the biggest bend because on a closed
        line bend located on the first and last bend are not simplified

        It is easier to move (rotate) the first/last vertice than to try to simplifiy the bend located
        of the first last vertice.  So we try to move the first/last vertice on a bend that do not need
        simplification
        
        Keyword definition
          line -- LineString object to calculate bends
        
        Return value: None  
        """

        if len(line._gbt_bends) >= 2:
            # Only process a line if there are two or more bends
            min_adj_area = -1.
            for bend in line._gbt_bends:
                if bend.adj_area > line._gbt_min_adj_area and bend.j-bend.i >= 4:
                    # A bend formed by 4 vertices (or more) is the ideal cases because if we split this bend
                    # in the middle it forms 2 smaller bends that are removed as we do not process the first and
                    # last bend of a closed line
                    i_j = (bend.i,bend.j)
                    break
                if bend.adj_area > min_adj_area:
                    # Acceptable bend found and favor the biggest bend
                    min_adj_area = bend.adj_area
                    i_j = (bend.i, bend.j)

            # Rotate the line at the center of the i_j
            i = i_j[0]
            j = i_j[1]
            if j-i == 2:
                mid_i_j = i+1
                del_first = True
                del_last  = True
            elif j-i == 3:
                mid_i_j = i + 1
                del_first = True
                del_last = True
            else:
                mid_i_j = int((i+j)/2)
                del_first = True
                del_last = True
            line.coords = line.coords[mid_i_j:] + line.coords[1:mid_i_j+1]
            print (line.coords[mid_i_j:])
            print (line.coords[1:mid_i_j+1])

            # Recreate the bends as the line as beed rotated with the first/last point on the biggest bend
            self.create_bends(line)

            # Delete the first bend if required
            if len(line._gbt_bends) >= 1 and del_first:
                del line._gbt_bends[0]

            # Delete the last bend if required
            if len(line._gbt_bends) >= 1 and del_last:
                del line._gbt_bends[-1]

        return


    def manage_lines_simplification (self):
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
                True: At least one line to process was simplified
                False: No line was simplified
        """
        
        line_simplified = False
        
        for line in self.s_container.get_features(filter= lambda feature: feature.geo_type == 'LineString' and not feature._gbt_is_simplest):

            # Create the bends in the line
            self.create_bends(line)

            if self.command.rotate_coord and line._gbt_is_closed:
                self.rotate_coordinates(line)

            # Classify each bend in the line
            self.classify_bends(line)
            
            last_bend = (len(line._gbt_bends))-1
            
            # Scan backwards all the bends of the lines from the last one to the first one
            i_bend = last_bend
            while (i_bend >= 0):

#                nbr_bends = self._number_of_bends_to_reduce(line._gbt_bends[i_bend].type)
                bend = line._gbt_bends[i_bend]
                if bend.type == _ONE_BEND:
                    bend_status = self._manage_bend_constraints(line, bend)
                    if (bend_status == _SIMPLIFIED):
                        line.coords = line.coords[0:bend.i] + bend.replacement_line.coords + line.coords[bend.j:]
                        line_simplified = True

                i_bend -= 1
                
            # Reset the bends to save some space...
            line._gbt_bends = []
            
#            # Update the Shapely structure using the coordinates in the temporary structure
#            if (line_simplified):
#                line.coords = new_coords
#            line.update_coords(in_hand_line_coords, self.s_container)
        
        return line_simplified

    def _manage_bend_constraints(self, line, bend):
        """Check if the bend to simplfy will violate the SIMPLE_LINE, LINE_CROSSING or SIDEDNESS constraint
        
        Parameters
            bend: Bend object to simplify
            line: line object to simplify

            
        Return value
            Flag indicating the status of the constraint verification
                _SIMPLIFIED: The bend simplification did not violate any constraint and was simplified
                _IN_CONFLICT:  The bend simplification did violate a constraint and was not simplified
        
        """
        in_conflict = False

        if self.command.simplicity and not in_conflict:
            start_line = GenUtil.create_LineString(line.coords[:bend.i+1])
            end_line = GenUtil.create_LineString(line.coords[bend.j:])

            in_conflict = self.spatial_constraints.validateSimplicity(bend, line, start_line, bend.replacement_line, end_line)

        if self.command.intersection and not in_conflict:
            in_conflict = self.spatial_constraints.validateIntersection(self.s_container, line, bend )

        if not in_conflict:
            # No conflict replace the bend by its replacement line
            status = _SIMPLIFIED
        else:
            status = _IN_CONFLICT
            bend.type = _IN_CONFLICT
                    
        return status

# The next lines are disable it was managing the alternate bends
# The alternate bends were adding far more complexity in the code for few better results
#    def _find_alternate_bends (self, line, i):
#        """Find an alternate bend when the bend to simplify is violating a constraint
#        
#        When a _ONE_BEND is in conflict we check if an alternate bend can be 
#        choosen  to be simplified instead. This strategy is implemented to avoid
#        that a dead lock occurs: the same bend is always selected and it is in 
#        conflict and nothing around get simplified.
#        
#        A bend is considered an alternate bend if it meets the following criterias
#        The bend[i] is a _ONE_BEND and 
#        The bend[i] is in conflict and
#        The bend[i-1] is < minimum adjusted area and the bend[i-2] is type UNKNOWN
#           ===> bend[i-1] is an alternate bend to be reducced 
#           
#        The bend[i] is a _ONE_BEND and 
#        The bend[i] is in conflict and
#        The bend[i+1] is < minimum adjusted area and the bend[i+2] is type UNKNOWN
#           ===> bend[i-1] is an alternate bend to be reducced
#           
#        Parameter:
#            line: line object to find an alternative bend
#            i: bend number for which to find alternate bend
#            
#        Return value
#            tuple containing 2 values; if different from None it corresponds to the bend number 
#            of the alternate bend
#           
#        """
#        
#        last_bend = len(line.bends)-1
#        
#        previous_bend = None
#        next_bend = None
#        
#        # Check if previous is an alternative bends
#        if (i !=0): # No alternate bends if we are on the first bend
#            if (i==1):
#                i_minus_two = _UNKNOWN
#            else:
#                i_minus_two = line.bends[i-2].type
#            
#            if (i_minus_two == _UNKNOWN):
#                previous_bend = self._validate_alternate_bends (line, i-1)
#        
#        # Check if the next bend is an alternative bend
#        if (i != last_bend):
#            if ( i == last_bend-1):
#                i_plus_two = _UNKNOWN
#            else:
#                i_plus_two = line.bends[i+2].type
#                
#            if (i_plus_two == _UNKNOWN):
#                next_bend = self._validate_alternate_bends (line, i+1)
#                
#        return (next_bend, previous_bend)
#            
#    def _validate_alternate_bends (self, line, i):
#        """Validate if the bend is a candidate for an alternate bends
#        
#        To be a candidate bend must conform to the rules of a bend to simplify which are the
#        area is smaller than the minimum adjusted area 
#        
#        Parameters:
#            line: line object to simplify
#            i: Bend number to check
#        
#        Return value
#             The position of the valid alternate bend or None if no alternate bend is possible
#             
#        """
#         
#        min_adj_area = line.min_adj_area
#        
#        # Check if the bend is below the minimum adjusted area
#        if (line.bends[i].adj_area <= min_adj_area):
#            line.bends[i].type = _ONE_BEND
#            lst_bends = [line.bends[i]]
#            
#            if (self._are_bends_big(lst_bends, line.min_adj_area)):
#                bend_size = _BIG
#            else:
#                bend_size = _SMALL
#                
#            self._reduce_one_bend (lst_bends, bend_size, line)
#            line.bends[i].multi_bends = lst_bends
#            line.bends[i].replacement_line = lst_bends[0].replacement_line
#        else:
#            i = None
#    
#        return i  
                
    # def _adjust_bend_special_cases (self, line):
    #     """Deal with bend sepcial cases.
    #
    #     If we don't simplify the first last bend than remove the first/last bend of the list
    #     If the line is closed and has only one bend than it is at its simplest form
    #
    #     Parameters
    #         line: line object to process
    #
    #     Return value
    #         None
    #
    #     """
    #
    #     # Option to keep or delete the first and last bend on the line
    #     if self.command.simplify_first_last:
    #         pass
    #     else:
    #         # Delete the first and last bend
    #         if (len(line._gbt_bends) >= 2):
    #             del line._gbt_bends[0]
    #             del line._gbt_bends[-1]
    #         elif (len(line._gbt_bends) == 1):
    #             del line._gbt_bends[0]
    #         else:
    #             pass
    #
    #     # Special case of a closed line with only one bend this line is at its simpliest form
    #     if (len(line._gbt_bends) == 1 and line._gbt_is_closed):
    #         line._gbt_bends = []
    #         line._gbt_simplest = True
    #
    #     return


#    def _create_replacement_line(self, bend, lst_coords):
#        """Create the replacement line for a bend, from a list of point"""
#
#        bend.replacement_line = MA_LineString(lst_coords)
#
#    def _reduce_one_bend(self, bend, line):
#        """Compute the replacement line for the case of a one bend
#
#        Parameter:
#            bend: bend to process (in the case of a one bend the list contains only one element)
#            line: line object to process
#
#        Return value:
#            None
#
#        """
#
#        bend._gbt_replacement_line = (line.coords[bend.i], line.coords[bend.j])
#
#
#        return
#
#     def _reduce_two_bends(self, lst_bends, bend_size, line):
#         """Compute the replacement line for the case of a two bend
#
#         Parameter:
#             lst_bends: list of bend to process
#             bend_size: Type of bend BIG or SMALL bend
#             line: line object to process
#
#         Return value:
#             None
#
#         """
#
#         bend1 = lst_bends[0]
#         bend2 = lst_bends[1]
#
#         # Calculate the mid position between the beginning of the first and second bend
#         base_mid_point = MA_Point(line.coords_dual[bend1.j]).mid_point(MA_Point(line.coords_dual[bend2.i]))
#
#         if (bend_size == _BIG):
#             # Calculate the middle coordinates of each bend
#             lst_coords_bend1 = self._reduce_bend(bend1, line)
#             lst_coords_bend2 = self._reduce_bend(bend2, line)
#             # Add the mid position and the middle coordinates of each bend
#             lst_coords = lst_coords_bend1 + list(base_mid_point.coords_dual) +  lst_coords_bend2
#
#         else:
#             # Bend is small
#             lst_coords = list(base_mid_point.coords_dual)
#
#         # Add the first/last vertice of the bends
#         lst_coords.insert(0, line.coords_dual[bend1.i])
#         lst_coords.append(line.coords_dual[bend2.j])
#
#         self._create_replacement_line(bend1, lst_coords)
#
#
#     def _reduce_three_bends(self, lst_bends, bend_size, line):
#         """Compute the replacement line for the case of a three bend
#
#         Parameter:
#             lst_bends: list of bend to process
#             bend_size: Type of bend BIG or SMALL bend
#             line: line object to process
#
#         Return value:
#             None
#
#         """
#
#         bend1 = lst_bends[0]
#         bend2 = lst_bends[1]
#         bend3 = lst_bends[2]
#
#         if (bend_size == _BIG):
#             # Calculate the middle coordinates of each bend
#             lst_coords_bend1 = self._reduce_bend (bend1, line)
#             lst_coords_bend2 = self._reduce_bend (bend2, line)
#             lst_coords_bend3 = self._reduce_bend (bend3, line)
#             lst_coords = lst_coords_bend1 + lst_coords_bend2 + lst_coords_bend3
#         else:
#             # Bend is small so only find the middle position of the second bend
#             lst_coords = self._reduce_bend (bend2, line)
#
#         # Add the first/last vertice of the bends
#         lst_coords.insert(0, line.coords_dual[bend1.i])
#         lst_coords.append(line.coords_dual[bend3.j])
#
#         self._create_replacement_line(bend1, lst_coords)
#
#
    # def _reduce_four_bends(self, lst_bends, bend_size, line):
    #     """Compute the replacement line for the case of a four bend
    #
    #     Parameter:
    #         lst_bends: list of bend to process
    #         bend_size: Type of bend BIG or SMALL bend
    #         line: line object to process
    #
    #     Return value:
    #         None
    #
    #     """
    #
    #     bend1 = lst_bends[0]
    #     bend2 = lst_bends[1]
    #     bend3 = lst_bends[2]
    #     bend4 = lst_bends[3]
    #
    #     # Calculate the mid position between the bends
    #     base_mid_point12 = MA_Point(line.coords_dual[bend1.j]).mid_point(MA_Point(line.coords_dual[bend2.i]))
    #     base_mid_point23 = MA_Point(line.coords_dual[bend2.j]).mid_point(MA_Point(line.coords_dual[bend3.i]))
    #     base_mid_point34 = MA_Point(line.coords_dual[bend3.j]).mid_point(MA_Point(line.coords_dual[bend4.i]))
    #
    #
    #     if (bend_size == _BIG):
    #         # Calculate the mid position of each bend
    #         lst_coords_bend1 = self._reduce_bend(bend1, line)
    #         lst_coords_bend2 = self._reduce_bend(bend2, line)
    #         lst_coords_bend3 = self._reduce_bend(bend3, line)
    #         lst_coords_bend4 = self._reduce_bend(bend4, line)
    #         lst_coords = lst_coords_bend1 + list(base_mid_point12.coords_dual) + \
    #                      lst_coords_bend2 + list(base_mid_point23.coords_dual) + \
    #                      lst_coords_bend3 + list(base_mid_point34.coords_dual) + \
    #                      lst_coords_bend4
    #
    #     else:
    #         # Bend is small so just use mid position between bends
    #         lst_coords = list(base_mid_point12.coords_dual) +\
    #                      list(base_mid_point23.coords_dual) +\
    #                      list(base_mid_point34.coords_dual)
    #
    #
    #
    #     # Add the first/last vertice of the bends
    #     lst_coords.insert(0, line.coords_dual[bend1.i])
    #     lst_coords.append(line.coords_dual[bend4.j])
    #
    #     self._create_replacement_line(bend1, lst_coords)
    #
    # def _reduce_five_bends(self, lst_bends, bend_size, line):
    #     """Compute the replacement line for the case of a five bend
    #
    #     Parameter:
    #         lst_bends: list of bend to process
    #         bend_size: Type of bend BIG or SMALL bend
    #         line: line object to process
    #
    #     Return value:
    #         None
    #
    #     """
    #
    #     bend1 = lst_bends[0]
    #     bend2 = lst_bends[1]
    # #    bend3 = lst_bends[2]
    #     bend4 = lst_bends[3]
    #     bend5 = lst_bends[4]
    #
    #     lst_coords = self._reduce_bend(bend2, line) + self._reduce_bend(bend4, line)
    #
    #     lst_coords.insert(0, line.coords_dual[bend1.i]) #
    #     lst_coords.append(line.coords_dual[bend5.j])   #
    #
    #     self._create_replacement_line(bend1, lst_coords)
     

    def _are_bends_big(self, bend, min_adj_area):
        """This routine if the bend are big bend
        
        To be considered as big bend the list of bends must have a adjusted area and the compactness index
        over a certain threshold     
        
        Parameters:
            bends: Bend to process
            min_adj_area: minimum adjusted area for the line
             
        Return value
            Boolean flag indicating if the bend are candidate for BIG bend
                True: The bends are big bend
                False: Otherwise
        """

        adj_area_ratio = 1.0 - (bend.adj_area/ min_adj_area)
        cmp_index_ratio = 1.0 - bend.cmp_index
        if (adj_area_ratio  < _BIG_BEND_MAX_ADJ_AREA_RATIO and
            cmp_index_ratio < _BIG_BEND_CMP_INDEX_RATIO ):
            is_big_bend = True
        else:
            is_big_bend = False

        return is_big_bend


    # def _number_of_bends_to_reduce (self, bend_type):
    #     """
    #     Output the number of bends depending on the type of bend
    #     _ONE_BEND ==> 1
    #     ...
    #     _FIVE_BENDS === 5
    #     If bend type is _IN_CONFLICT, UNKNOWN the value 0 is returned
    #     """
    #
    #     if (bend_type == _ONE_BEND):
    #         nbr_bends = 1
    #     elif (bend_type == _TWO_BENDS):
    #         nbr_bends = 2
    #     elif (bend_type == _THREE_BENDS):
    #         nbr_bends = 3
    #     elif (bend_type == _FOUR_BENDS):
    #         nbr_bends = 4
    #     elif (bend_type == _FIVE_BENDS):
    #         nbr_bends = 5
    #     else:
    #         nbr_bends = 0
    #
    #     return nbr_bends

    def _finalize_bend_information (self, bend, line, in_hand_line_coords):
        """This routine calculates extra attributes and informaion for the bends and line  that will be simplified.
        
        It calculates the following values:
        replacement_line bounding_box: Extent values for the replacement line
        bend_polygon: The list of points forming the polygon of all the bends (more than one bend with a multi bends)
        bend_polygon bounding box: Extent of the bend polygon
        
        Parameters:
            bend: Bend object to simplify
            line: LineString being processed
            in_hand_line_coords: In order to increase performance we are not rewriting the coordinates in
                                 the Sahpely structure we keep them in that list
        
        """   
        
        first_pt = bend.i  # First bend first point
        last_pt = bend.j  # Last bend last point
            
       
        # Calculate the complete new line with the bend removed
        new_line_coords = list(in_hand_line_coords)
        new_line_coords[first_pt:last_pt+1] = bend.replacement_line.coords_dual[0:]
        bend.new_line = LineString(new_line_coords)
        
        # Calculate the sidedness region to use for the constraints
#        old_sub_line = list(line.coords[first_pt:last_pt+1])
#        new_sub_line = list(bend.replacement_line.coords[0:])
#        bend.sidedness_polygon = GenUtil.calculate_sidedness_polygon(LineString(old_sub_line),
#                                                                     LineString(new_sub_line) )
        

    def _calculate_number_of_bends (self, lines):
        """
        This routine is calculating the total number of bends in the lines
        """
        
        nbr_bends = 0
        for line in lines:
            nbr_bends += len(line.bends)
            
        return nbr_bends
    
    # def _reduce_bend(self, bend, line, bend_depth=_DEPTH_OFFSET_RATIO):
    #
    #     """This routine reduce one bend"""
    #
    #     first_coord = line.coords[bend.i]
    #     last_coord = line.coords[bend.j]
    #
    #     if (bend.base / bend.depth < 1.0):
    #         offset_factor = bend.base / bend.depth
    #     else:
    #         offset_factor = bend.depth / bend.base
    #
    #     bend_depth = bend_depth * offset_factor
    #
    #     # Calculate the position of the middle of the bend base line
    #     base_mid_point = Point(first_coord).mid_point(Point(last_coord))
    #
    #     # Offset the point in direction of the bend peak
    #     scaled_point = self.rescale_vector (base_mid_point, MA_Point(bend.peak_coords),  bend_depth)
    #
    #     return [scaled_point.coords_dual[0]]

    def _are_bends_similar(self, lst_bends):
        
        """This routine determines determine if the list of bends are similar"""
        
        lst_bend_base      = []
        lst_bend_adj_area  = []
        lst_bend_cmp_index = []
        
        # Create the list for the base, adj_area and cmp_index
        for bend in lst_bends:
            lst_bend_base.append(bend.base)
            lst_bend_adj_area.append(bend.adj_area)
            lst_bend_cmp_index.append(bend.cmp_index)
        
        # Extract max list values
        max_base      = max(lst_bend_base)
        max_adj_area  = max(lst_bend_adj_area)
        max_cmp_index = max(lst_bend_cmp_index)
        
        # Extract min values
        min_base      = min(lst_bend_base)
        min_adj_area  = min(lst_bend_adj_area)
        min_cmp_index = min(lst_bend_cmp_index)
        
        delta_base      = 1.0 - (min_base / max_base)
        delta_cmp_index = 1.0 - (min_cmp_index / max_cmp_index)
        delta_adj_area  = 1.0 - (min_adj_area / max_adj_area)
        
        if (delta_base < _SIMILAR_BEND_BASE_RATIO and
             delta_cmp_index < _SIMILAR_BEND_CMP_INDEX_RATIO and
             delta_adj_area < _SIMILAR_BEND_ADJ_AREA_RATIO):
            similar_bends = True
        else:
            similar_bends = False
        
        return similar_bends
    
    # def _replace_bend (self, bend, line, in_hand_line_coords):
    #     """
    #     This function replace  the vertices of the bend (excluding the first and last point) by  the replacement
    #     line (excluding the first and last point which are the same as the first/last vertice of the bend)
    #
    #     Parameters:
    #        bend: Bend to process
    #        line: LineString object being processed
    #        in_hand_line_coords: In order to increase performance we are not rewriting the coordinates in
    #                             the Sahpely structure we keep them in that list
    #     """
    #     # Determine the position of the bend limits
    #     first = bend_to_check.i
    #     last = bend_to_check.j
    #
    #     # Replace the coordinate with the new line in the temporary line coordinates (but we do not create a new object)
    #     in_hand_line_coords[0:] = list(bend_to_check.new_line.coords)
    #
    #     if (self._are_bends_big([bend_to_check], line.min_adj_area)):
    #         self.stats.add_stats(_BIG)
    #     else:
    #         self.stats.add_stats(_SMALL)
    #
    #     self.stats.add_stats(_ALGO)
    #     self.stats.add_stats(bend_to_check.type)
    #
    #     return
        
#    def rescale_vector(self, p1, p2, scale_factor):
#        """This routine rescale the vector defined by the points P1 and P2 by a factor
#        of SCALE_FACTOR
#
#        P1: Point vector origin (first point)
#        P2: Point vector to rescale (second point)
#        scale_factor: factor to scale the vector (same for x and y)
#
#        """
#
#        x1 = p1.coords_dual[0][0]
#        y1 = p1.coords_dual[0][1]
#        x2 = p2.coords_dual[0][0]
#        y2 = p2.coords_dual[0][1]
#
#        vec_x = x2 - x1
#        vec_y = y2 - y1
#
#        vec_x = vec_x * scale_factor
#        vec_y = vec_y * scale_factor
#
#        x_out = vec_x + x1
#        y_out = vec_y + y1
#
#        return MA_Point([x_out, y_out])
    
    def _calculate_bend_depth(self, line, bend):
        """
        This routine extract the depth of a bend. 
        
        The depth is defined as the furthest point 
        from the base for the points that delineate the bend
        
        Returned values:
            a tuple containing
            (the depth as a distance, index of the point on the line where the peak is)
            
        """
        
        dist_depth = - 1.
        
#        bend_base_line = LineString((line.coords_dual[bend.i], line.coords_dual[bend.j]))
        base_line_p1 = line.coords[bend.i]
        base_line_p2 = line.coords[bend.j]
        
        for n in range(bend.i + 1, bend.j):
            
#            distance = bend_base_line.distance(Point(line.coords_dual[n]))
            distance = GenUtil.distance_line_point (base_line_p1, base_line_p2, line.coords[n])
            if (distance > dist_depth):
                dist_depth = distance
                pt_depth_index = n
    
        return (dist_depth, pt_depth_index)


#    def check_features(self):
#        """Check if the features passed in parameters are of the good class type and have the good attributes
#
#        Parameters: None
#
#        Return value: None
#
#        """
#
#        # Check the line string
#        properties_to_check_to_check = [_DIAMETER, _SIMPLIFY_FIRST_LAST]
#        for feature in self.features:
#            if isinstance(feature, MA_LineString):
#                GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check_to_check)
    

    def process(self):
        """Main routine for the Sherbend algorithm
        
        The algorithm will simplify the lines using the Sherbend algorithm. 
        It will iterate over the lines until there are no more bends to simplify.

        Keyword arguments:
            none

        Return value:
            geo_content: dataclass containing the output information

        """



        self.spatial_constraints = SpatialConstraints()

        # Load the features into the spatial container
        self.load_features(self.geo_content.features)
     
        self.add_line_attributes(self.command.diameter)
        
#        if (self.command.multi_bend):
#            nbr_step = 2
#        else:
        nbr_step = 1
        iter_nbr = 0
        for step in range(nbr_step):
            # The last step is always done with multi_bend to False.  We always process the last iterations of multi_bend to False
            # in order to process alternative bends correctly. An alternative bend is involved when a bend is in conflict and it tries
            # to simplify the bend just before or after the bend in conflict
            if (step == 1):
                self.command.multi_bend = False
            
            # Loop until no more lines are simplifiable
            line_simplified = True
            # Iterate until all the line are simplified or there are no more line have to be simplified
            while (line_simplified):
                print ('Start of iteration # {}'.format(iter_nbr ))
                line_simplified=False
                line_simplified = self.manage_lines_simplification()
                print('Number of bend simplified {}'.format(10))
                print('End of iteration # {}'.format(10))

                iter_nbr += 1
                    
        out_features = []
        for feature in self.s_container.get_features():
            if feature._gbt_geom_type == 'Point':
                out_features.append(feature)
            elif feature._gbt_geom_type == 'LineString':
                if feature._gbt_original_type == 'LineString':
                    out_features.append(feature)
                else:
                    if feature._gbt_original_type == 'Polygon-Exterior':
                        # The LineString was an exterior Polygon so reconstruct the originalPolygon
                        interiors = [list(interior.coords) for interior in feature._gbt_interiors]
                        polygon = Polygon(feature.coords, interiors)
                        polygon._gbt_layer_name = feature._gbt_layer_name
                        polygon._gbt_properties = feature._gbt_properties
                        out_features.append(polygon)
                    else:
                        pass  # Nothing to do with the holes here

        return out_features
