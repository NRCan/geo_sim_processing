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
from shapely.prepared import prep

from algo_bends import AlgoBends
from lib_geobato import MA_LineString, GenStatistics, Algorithm, GenUtil, \
                         SpatialContainer, MA_Point, Polygon, Parameters
                               
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

class SherbendStatistics(GenStatistics):
    """Class that contains the statistics for the Sherbend algorithm
    
    Attributes
        stat_names: Name of the statistics for the SherbendStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class Sherbend statistics
        
        Parameter: 
            None
            
        Return value
            None
            
        """
        
        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS,  \
                             _BIG, _SMALL, _ONE_BEND, _TWO_BENDS, _THREE_BENDS, _FOUR_BENDS, _FIVE_BENDS ))
        
    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information
        
        """
        
        str_out = []
        str_out.append( "Sherbend algorithm Statistics" )
        
        str_out.append( "--------------------------" )
        if (type == GenStatistics.DETAILED):
            for i in xrange((self.get_nbr_iteration())):
                str_out.append("Detailed statistics")
                str_out.append("Iteration # " + str(i))
                str_out.append("Bend simplified: " + str(self.get_stats_name_count_iter( _ALGO, i)))
                str_out.append( "--------" )
                str_out.append( "Conflicts:" )
                str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_iter( GenUtil.SIMPLE_LINE, i)))
                str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
                str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
                str_out.append( "--------" )
        str_out.append( "Summary statistics" )
        str_out.append("Total bend simplified: " + str(self.get_stats_name_count_total(_ALGO)))
        str_out.append( "--------" )
        str_out.append("Statistics by bends")
        str_out.append("     One bend      :  " +  str(self.get_stats_name_count_total(_ONE_BEND)))
        str_out.append("     Two bends     :  " +  str(self.get_stats_name_count_total(_TWO_BENDS)))
        str_out.append("     Three bends   :  " +  str(self.get_stats_name_count_total(_THREE_BENDS)))
        str_out.append("     Four bends    :  " +  str(self.get_stats_name_count_total(_FOUR_BENDS)))
        str_out.append("     Five bends    :  " +  str(self.get_stats_name_count_total(_FIVE_BENDS)))
        str_out.append("     Small bends   :  " +  str(self.get_stats_name_count_total(_SMALL)))
        str_out.append("     Big bends    :  " +  str(self.get_stats_name_count_total(_BIG)))
        str_out.append( "--------" )
        str_out.append( "Conflicts:" )
        str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_total( GenUtil.SIMPLE_LINE)))
        str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
        str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
        str_out.append( "--------" )
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )
        
        return str_out

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
                
    def add_line_attributes (self):
        """This routine sets different attributes of the lines

        Parameters: None
            
        Return value: None

        """
    
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            
            line.is_simplest = False
            line.bends = []
            line.min_adj_area = None
            line.multi_bend = False
            
            # Set the minimal adjusted area
            ray = line.ma_properties[_DIAMETER] / 2.0
            line.min_adj_area = _AREA_CMP_INDEX * math.pi * ray**2.0
            
            # Check if the line is a closed line by comparing  if the first 
            # and the last vertice are almost at the same place
            if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                line.is_closed = True
            else:
                line.is_closed = False
                
            # Set the multi bend attribute
            line.multi_bend = self.params.multi_bend
            if (line.is_closed):
                polygon = Polygon(line.coords_dual)
                # The next lines of codes are within "try .. except" as in some cases the closed line
                # may form a bad polygon.  The command line.polygon_line.exterior will create an exception
                # in that case we disable the multi bends
                try:
                    if (polygon.area <= 5. * line.min_adj_area):
                        # If the line is closed and the area of the polygon is near by the minimum adjusted area
                        # we disable the multi bend option because in that case the output is better with single bend
                        line.multi_bend = False
                except ValueError:
                    line.multi_bend = False
                            
        return

    def classify_bends (self, line):   
        """Classify the beds
        
        High level routine to classify the bends from _ONE_BEND to _FIVE_BENDS and to calculate
        the replacement line for the bends that are not of type "_UNKNOWN"
        
        Parameter: 
            s_container: Object of type SpatialContainer containing the features nd the spatial index
            
        Return value:
            None
        
        """ 
    
        if (self.params.multi_bend and line.multi_bend):
            self._detect_window_bend (line, _FIVE_BENDS)
            self._detect_window_bend (line, _FOUR_BENDS)
            self._detect_window_bend (line, _THREE_BENDS)
            self._detect_window_bend (line, _TWO_BENDS)
            self._detect_window_bend (line, _ONE_BEND)
        else:                
            self._detect_window_bend (line, _ONE_BEND)
            
        # Detect and flag lines at their simplest form
        self._set_simplest_line(line)
        
        # Compaction of the bends in the list of bends.  This compactness will speedup the search process
        self._compact_bends(line)
        
        # Calculates the replacement line for each bend found
        self._reduce_bends(line)
            
    def _set_simplest_line(self, line):
        """ Determine if a line is at its simplest form
        
        Set the attribute of the simplest at true when a line as 0 bends (the bend list is empty) or
        all the bends of the line are of type "UNKNOWN"  this means that all the bends are over 
        the minimum adjusted size
        
        Parameter:
            line: MA_LineString to check if it is at its simplest form
            
        Return value
            None
            
        """
        
        if (line.bends):
            list_unknown = [bend.type for bend in line.bends if bend.type != _UNKNOWN]
            # If list is empty ==> All the bends are of type _UNKNOWN and the line is at its simplest form
            if list_unknown:
                pass
            else:
                line.simplest= True
        else:
            line.simplest = True

    def _detect_window_bend(self, line, bend_type):
        """This routine is sliding a window over the bends of a line. 
        
        It tries to find consecutives bends
        that meet the following criterias:
          - The bends must have the bend_type _UNKNOWN
          - The bends must surrounded by bend of type _UNKNOWN
          - The bends must have a adjusted area below min_adj_area
          - The first and last bend of the window must be smaller than the 
               previous and next bend of the window
               
        Parameter:
            line: line object to check for consecutives bend
            bend_type: type of multi bend to check (from 1 to 5)
            
        Return value:
            None
    
        """
        
        window_length = self._number_of_bends_to_reduce(bend_type)
        min_adj_area = line.min_adj_area
        
        nbr_bends = len(line.bends)
        last_bend = nbr_bends - window_length 
        
        # The last bend to process depends on the window length
        for i in range(0, last_bend + 1):
            # The sliding window must be surrounded by bend type _UNKNOWN 
            # otherwise we cannot check if they are similar bends
            if (i == 0):
                # Special case for the first bend
                surround_before = _UNKNOWN
            else:
                surround_before = line.bends[i - 1].type
            if (i == last_bend):
                #How the last bend we assume the next bend to _UNKNOWN
                surround_after = _UNKNOWN
            else:
                surround_after = line.bends[i + window_length].type
            
            # Check that all the bends over the window are of type UNKNOWN
            bend_unknown = True
            for j in range(window_length):
                if (line.bends[i + j].type != _UNKNOWN):
                    bend_unknown = False
            
            # Check previous and next bend are of type _UNKNOWN
            if (bend_unknown and
                 surround_before == _UNKNOWN  and
                 surround_after == _UNKNOWN):
                
                # Check that all bends are below min_adj_area
                window_adj_area = True
                for j in range(window_length):
                    if (line.bends[i + j].adj_area > min_adj_area):
                        window_adj_area = False
                    
                if (window_adj_area):
                    
                    if ( self._are_neighbours_bigger(line, i, i+window_length-1) ):
                        
                        lst_bends = []
                        for i_lst in range(window_length):
                            lst_bends.append(line.bends[i_lst+i])
                            
                        if ( self._are_bends_similar(lst_bends) ):
                            for i_lst in range(window_length):
                                line.bends[i_lst+i].type = bend_type
                        
        return
    
    def _are_neighbours_bigger(self, line, first_bend, last_bend):
        """Checks if the bends is surrounded by bigger bends. This means bigger adjusted area.
        
        Parameters:
            line: Line object to check
            first_bend: number of the first bend of the line to check
            last_bend: Number of the last bend of the line to check
            
        Return value:
            Boolean flag for the state of the bend in comparison if its neighbours 
                True: The neighbours are bigger
                False:  The neighbours are smaller
        
        """
        
        # If it is the first bend of the line we assume previous bend as infinite
        if (first_bend == 0):
            previous_adj_area = 1.0e+99
        else:
            previous_adj_area = line.bends[first_bend - 1].adj_area
        
        # If it is the last bend of the line we assume the next bend as infinite    
        if (last_bend == len(line.bends) - 1):
            next_adj_area = 1.0e+99
        else:
            next_adj_area = line.bends[last_bend + 1].adj_area
            
        if (previous_adj_area >= line.bends[first_bend].adj_area and
             next_adj_area >= line.bends[last_bend].adj_area):
            
            bigger = True
        else:
            bigger = False
            
        return bigger

    
    def _compact_bends (self, line):
        """This routine compacts the bends.
        
        When there are multiple bends we compact them as one bend and with the those bends
        in the attribute multi_bends
        
        Parameters:
            line: Line object to compact the bends
            
        Return value:
            None
            
        Note: We start the compaction by the end to avoid index problem when it remove
              number
        """
        
        i = len(line.bends)-1 # i is on the last bend
        while (i >= 0):
            nbr_bends = self._number_of_bends_to_reduce (line.bends[i].type)
            if (nbr_bends >=1):
                lst_bends = []
                for j in range(nbr_bends):
                    lst_bends.append(line.bends[i-j])
                lst_bends.reverse()
                
                i = i - (nbr_bends -1) # Go to the first bend of the multi bend
                line.bends[i].multi_bends = lst_bends  
                 
                # Keep the first bend and delete the remaining bends
                for dummy in range(nbr_bends-1):
                    del line.bends[i+1]
    
            i -= 1
            
    def _reduce_bends(self, line):
        """Detect the bends to simplify in the line 
        
        Parameters: 
            line: Line object to detect bend
        """
        
        for bend in line.bends:
        
            # Choose the good method depending on the number of bend to simplify 
            if (bend.type == _ONE_BEND):
                sim_function = self._reduce_one_bend
            elif (bend.type == _TWO_BENDS):
                sim_function = self._reduce_two_bends
            elif (bend.type == _THREE_BENDS):
                sim_function = self._reduce_three_bends
            elif (bend.type == _FOUR_BENDS):
                sim_function = self._reduce_four_bends
            elif (bend.type == _FIVE_BENDS):
                sim_function = self._reduce_five_bends
                
            nbr_bends = self._number_of_bends_to_reduce (bend.type)
            
            if (nbr_bends >= 1):
                lst_bends = bend.multi_bends
                
                if (self._are_bends_big(lst_bends, line.min_adj_area)):
                    bend_size = _BIG
                else:
                    bend_size = _SMALL
                
                sim_function (lst_bends, bend_size, line)
                
                bend.replacement_line = lst_bends[0].replacement_line
                bend.i = lst_bends[0].i
                bend.j = lst_bends[-1].j
                
                bend.multi_bends = None

    def detect_bend_location(self, line): 
        """Calculates the position of each individual bends in the line
        
        The position of the bends are calculated according to the definition of the bencds 
        in the orginal paper Wang 1998.
        
        Parameter: LineString object to calculate bends
        
        Return value: None  
        """
        
        # Calculate the bends on the line
        algo_bends = AlgoBends()
        algo_bends.features.append( line )
        algo_bends.process()
        
        for bend in line.bends:
            # Add properties to the bends needed by the SherBend algorithm
            bend.area = None             # The area of the bend
            bend.perimeter = None        # The perimeter of the bend
            bend.cmp_index = None        # Compactness index (base ** 2) / (perimeter - base)
            bend.adj_area = None         # Adjusted area of the bend area * 0.75 /cmp_index)
            bend.type = _UNKNOWN         # Type of bend by default:UNKNOWN
            bend.base = None             # The length of the base of the bend
            bend.depth = None            # The distance from the bend peak to the base_line
            bend.peak_index = None       # Vertice number of the peak's bend 
            bend.polygon = None          # Polygon of the bend
            bend.replacement_line = None # LineString that will replace the bend   
            bend.multi_bends = None      # List of the bends contains in the multi bend 
            
        # Manage the bend for the first/last bend and for closed line
        self._adjust_bend_special_cases (line)
        
        # Calcultes extra attributes on the line
        self._add_bend_attributes(line)

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
        
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
                        
            # Copy the coordinates of the line in hand.  During the processing if the line in hand all the modifications
            # area written in the variable in_hand_line_coords.  at the end of the processing the coordinates are rewritten in the 
            # Shapely structure.  This work is done only for performance issue   
            in_hand_line_coords = list(line.coords_dual)
            
            # For the line in hand detect where are the bends and classify the type of bends
            self.detect_bend_location(line)
            self.classify_bends(line)
            
            last_bend = (len(line.bends))-1
            
            # Scan backwards all the bends of the lines from the last one to the first one
            i_bend = last_bend
            while (i_bend >= 0):
                            
                nbr_bends = self._number_of_bends_to_reduce(line.bends[i_bend].type)        
                if (nbr_bends >=1):
                    bend_status = self._manage_bend_constraints(i_bend, line, in_hand_line_coords)
                    if (bend_status == _SIMPLIFIED):
                        line_simplified = line_simplified or True
#                    # The next lines are now disable as well as the routine to find alternate bends
#                    # i.e. _find_alternate_bends  The routine was used determine the best alternate
#                    # bend to simplify when a bend to simplify was not simplify because of a contraint violation
#                    # This concept was not used very often and needed more complex code.
#                    else:
#                        # Alternate bend simplification for _ONE_BEND
#                        if (nbr_bends == 1):
#                            # Detect if the sibbling bends to the ONE_BEND in conflict could be a candidate 
#                            # bend for bend simplification
#                            (next_bend, preceding_bend) = self._find_alternate_bends (line, i_bend)
#                            for alt_bend in (next_bend, preceding_bend):
#                                if ( alt_bend is not None):
#                                    bend_status = self._manage_bend_constraints(alt_bend, line, in_hand_line_coords)
#                                    if (bend_status == _SIMPLIFIED):
#                                        line_simplified = line_simplified or True
#                            if (preceding_bend is not None):
#                                i_bend -= 1
                # Skip to the preceding bend
                i_bend -= 1
                
            # Reset the bends to save some space...
            line.bends = []
            
            # Update the Shapely structure using the coordinates in the temporary structure
            line.update_coords(in_hand_line_coords, self.s_container)
        
        return line_simplified

    def _manage_bend_constraints(self, i_bend, line, in_hand_line_coords):
        """Check if the bend to simplfy will violate the SIMPLE_LINE, LINE_CROSSING or SIDEDNESS constraint
        
        Parameters
            i_bend: Bend object to simplify
            line: line object to simplify
            in_hand_line_coords: In order to increase performance we are not rewriting the coordinates in
                                 the Sahpely structure we keep them in that list
            
        Return value
            Flag indicating the status of the constraint verification
                _SIMPLIFIED: The bend simplification did not violate any constraint and was simplified
                _IN_CONFLICT:  The bend simplification did violate a constraint and was not simplified
        
        """
        
        
        conflict_type = None
        bend_to_check = line.bends[i_bend]
                
        # Build the features used for constraint checking 
        self._finalize_bend_information(bend_to_check, line, in_hand_line_coords)
        line_simple_line = bend_to_check.new_line
        line_crossing_line = bend_to_check.replacement_line
        sideness_polygon = bend_to_check.sidedness_polygon
                
        conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                  sideness_polygon, self.s_container, line._sci_id)

        if (conflict_type is None):
            # No conflict replace the bend by its replacement line
            self._replace_bend (i_bend, line, in_hand_line_coords)
            status = _SIMPLIFIED
        else:
            status = _IN_CONFLICT
            line.bends[i_bend].type = _IN_CONFLICT
                    
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
                
    def _adjust_bend_special_cases (self, line):
        """Deal with bend sepcial cases.
        
        If we don't simplify the first last bend than remove the first/last bend of the list
        If the line is closed and has only one bend than it is at its simplest form
        
        Parameters
            line: line object to process
            
        Return value
            None

        """
    
        # Option to keep or delete the first and last bend on the line    
        if line.ma_properties[_SIMPLIFY_FIRST_LAST]:
            pass
        else:
            # Delete the first and last bend
            if (len(line.bends) >= 2):
                del line.bends[-1]
                del line.bends[0]
            elif (len(line.bends) == 1):
                del line.bends[0]
            else:
                pass
            
        # Special case of a closed line with only one bend this line is at its final state
        if (len(line.bends) == 1 and line.is_closed):
            line.bends = []
            line.simplest = True
            
        return
        
    def _add_bend_attributes(self, line):       
        """This routine computes the different attributes of the bends of a line.
        
        Parameter:
            line: line object to process
            
        Return value:
            None
        """
        
        bend_to_simplify = False
        for bend in line.bends:
            i = bend.i
            j = bend.j
            bend_polygon = Polygon(line.coords_dual[i:j+1])
            area = bend_polygon.area
            if area <= GenUtil.ZERO:
                area = GenUtil.ZERO # In case of area=0 we assume very small area
            else:
                area = area
            base = GenUtil.distance (line.coords_dual[i], line.coords_dual[j])
            if base <= GenUtil.ZERO: base = GenUtil.ZERO  # Avoid a case of division by zero
            perimeter = bend_polygon.length
            bend.base = base
            bend.perimeter = perimeter
            bend.cmp_index = 4 * area * math.pi / (perimeter ** 2.0)
            bend.adj_area = area * (0.75 / bend.cmp_index)
            
            # Because the evaluation on the next attributes are costly in time machine
            # It has to call many Shapely routines those attributes are only evaluated
            # for bend that are below the min_adj_are and are real candidates for 
            # simplification
            if (bend.adj_area <= line.min_adj_area):
                bend_to_simplify = True
                (depth, peak_index) = self._calculate_bend_depth(line, bend)
                bend.depth = depth
                bend.peak_index = peak_index
                bend.peak_coords = line.coords_dual[peak_index]
                
        # If there are no potential bend to simplify the line is at its simplest form
        if not bend_to_simplify:
            line.is_simplest = True
                
    def _create_replacement_line(self, bend, lst_coords):
        """Create the replacement line for a bend, from a list of point"""
        
        bend.replacement_line = MA_LineString(lst_coords)
    
    def _reduce_one_bend(self, lst_bends, bend_size, line):
        """Compute the replacement line for the case of a one bend
        
        Parameter:
            lst_bends: list of bend to process (in the case of a one bend the list contains only one element)
            bend_size: Type of bend BIG or SMALL bend
            line: line object to process
        
        Return value:
            None
        
        """
       
        bend = lst_bends[0]
        
        if (bend_size == _BIG):
            # Find the coordinate of the replacement line
            lst_coord = self._reduce_bend (bend, line)
        else:
            # Bend is small no middle coordinates
            lst_coord = []
        
        # Add first last coordinates
        lst_coord.insert(0, line.coords_dual[bend.i]) 
        lst_coord.append(line.coords_dual[bend.j])
        
        self._create_replacement_line(bend, lst_coord)
    
        return
    
    def _reduce_two_bends(self, lst_bends, bend_size, line):
        """Compute the replacement line for the case of a two bend
        
        Parameter:
            lst_bends: list of bend to process
            bend_size: Type of bend BIG or SMALL bend
            line: line object to process
        
        Return value:
            None
        
        """  
        
        bend1 = lst_bends[0]
        bend2 = lst_bends[1]
        
        # Calculate the mid position between the beginning of the first and second bend
        base_mid_point = MA_Point(line.coords_dual[bend1.j]).mid_point(MA_Point(line.coords_dual[bend2.i]))    
        
        if (bend_size == _BIG):
            # Calculate the middle coordinates of each bend
            lst_coords_bend1 = self._reduce_bend(bend1, line)
            lst_coords_bend2 = self._reduce_bend(bend2, line)
            # Add the mid position and the middle coordinates of each bend
            lst_coords = lst_coords_bend1 + list(base_mid_point.coords_dual) +  lst_coords_bend2
        
        else:
            # Bend is small
            lst_coords = list(base_mid_point.coords_dual)
    
        # Add the first/last vertice of the bends
        lst_coords.insert(0, line.coords_dual[bend1.i]) 
        lst_coords.append(line.coords_dual[bend2.j])   
    
        self._create_replacement_line(bend1, lst_coords)
      
      
    def _reduce_three_bends(self, lst_bends, bend_size, line):
        """Compute the replacement line for the case of a three bend
        
        Parameter:
            lst_bends: list of bend to process
            bend_size: Type of bend BIG or SMALL bend
            line: line object to process
        
        Return value:
            None
        
        """
    
        bend1 = lst_bends[0]
        bend2 = lst_bends[1]
        bend3 = lst_bends[2]
        
        if (bend_size == _BIG):
            # Calculate the middle coordinates of each bend
            lst_coords_bend1 = self._reduce_bend (bend1, line)
            lst_coords_bend2 = self._reduce_bend (bend2, line)
            lst_coords_bend3 = self._reduce_bend (bend3, line)
            lst_coords = lst_coords_bend1 + lst_coords_bend2 + lst_coords_bend3
        else:
            # Bend is small so only find the middle position of the second bend
            lst_coords = self._reduce_bend (bend2, line)
    
        # Add the first/last vertice of the bends
        lst_coords.insert(0, line.coords_dual[bend1.i]) 
        lst_coords.append(line.coords_dual[bend3.j]) 
            
        self._create_replacement_line(bend1, lst_coords)
            
    
    def _reduce_four_bends(self, lst_bends, bend_size, line):
        """Compute the replacement line for the case of a four bend
        
        Parameter:
            lst_bends: list of bend to process
            bend_size: Type of bend BIG or SMALL bend
            line: line object to process
        
        Return value:
            None
        
        """
        
        bend1 = lst_bends[0] 
        bend2 = lst_bends[1]
        bend3 = lst_bends[2]
        bend4 = lst_bends[3]
        
        # Calculate the mid position between the bends
        base_mid_point12 = MA_Point(line.coords_dual[bend1.j]).mid_point(MA_Point(line.coords_dual[bend2.i]))
        base_mid_point23 = MA_Point(line.coords_dual[bend2.j]).mid_point(MA_Point(line.coords_dual[bend3.i]))
        base_mid_point34 = MA_Point(line.coords_dual[bend3.j]).mid_point(MA_Point(line.coords_dual[bend4.i]))    
        
        
        if (bend_size == _BIG):
            # Calculate the mid position of each bend
            lst_coords_bend1 = self._reduce_bend(bend1, line)
            lst_coords_bend2 = self._reduce_bend(bend2, line)
            lst_coords_bend3 = self._reduce_bend(bend3, line)
            lst_coords_bend4 = self._reduce_bend(bend4, line)
            lst_coords = lst_coords_bend1 + list(base_mid_point12.coords_dual) + \
                         lst_coords_bend2 + list(base_mid_point23.coords_dual) + \
                         lst_coords_bend3 + list(base_mid_point34.coords_dual) + \
                         lst_coords_bend4
    
        else:
            # Bend is small so just use mid position between bends
            lst_coords = list(base_mid_point12.coords_dual) +\
                         list(base_mid_point23.coords_dual) +\
                         list(base_mid_point34.coords_dual) 
    
        
        
        # Add the first/last vertice of the bends
        lst_coords.insert(0, line.coords_dual[bend1.i]) 
        lst_coords.append(line.coords_dual[bend4.j]) 
    
        self._create_replacement_line(bend1, lst_coords)
           
    def _reduce_five_bends(self, lst_bends, bend_size, line):
        """Compute the replacement line for the case of a five bend
        
        Parameter:
            lst_bends: list of bend to process
            bend_size: Type of bend BIG or SMALL bend
            line: line object to process
        
        Return value:
            None
        
        """
    
        bend1 = lst_bends[0]
        bend2 = lst_bends[1]
    #    bend3 = lst_bends[2]
        bend4 = lst_bends[3]
        bend5 = lst_bends[4]
            
        lst_coords = self._reduce_bend(bend2, line) + self._reduce_bend(bend4, line)
        
        lst_coords.insert(0, line.coords_dual[bend1.i]) #
        lst_coords.append(line.coords_dual[bend5.j])   #
        
        self._create_replacement_line(bend1, lst_coords)
     

    def _are_bends_big(self, bends, min_adj_area):
        """This routine if the bend are big bend
        
        To be considered as big bend the list of bends must have a adjusted area and the compactness index
        over a certain threshold     
        
        Parameters:
            bends: List of bend to process
            min_adj_area: minimum adjusted area for the line
             
        Return value
            Boolean flag indicating if the bend are candidate for BIG bend
                True: The bends are big bend
                False: Otherwise
        """
        
        if (self.params.big_bend):
            nbr_bends = len(bends)
            list_adj_area = []
            list_cmp_index = []
            for i in range(nbr_bends):
                list_adj_area.append(bends[i].adj_area)
                list_cmp_index.append(bends[i].cmp_index)
            
            delta_adj_area  = 1.0 - (min(list_adj_area)/ min_adj_area)
            delta_cmp_index = 1.0 - (min(list_cmp_index))
            if (delta_adj_area  < _BIG_BEND_MAX_ADJ_AREA_RATIO and
                delta_cmp_index < _BIG_BEND_CMP_INDEX_RATIO ):  
                big = True
            else:
                big = False
        else:
            big = False
            
        return big

    def _number_of_bends_to_reduce (self, bend_type):
        """
        Output the number of bends depending on the type of bend
        _ONE_BEND ==> 1
        ...
        _FIVE_BENDS === 5
        If bend type is _IN_CONFLICT, UNKNOWN the value 0 is returned
        """
        
        if (bend_type == _ONE_BEND):
            nbr_bends = 1    
        elif (bend_type == _TWO_BENDS):
            nbr_bends = 2
        elif (bend_type == _THREE_BENDS):
            nbr_bends = 3
        elif (bend_type == _FOUR_BENDS):
            nbr_bends = 4
        elif (bend_type == _FIVE_BENDS):
            nbr_bends = 5
        else:
            nbr_bends = 0
                   
        return nbr_bends

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
        old_sub_line = list(line.coords_dual[first_pt:last_pt+1])
        new_sub_line = list(bend.replacement_line.coords_dual[0:])
        bend.sidedness_polygon = GenUtil.calculate_sidedness_polygon(LineString(old_sub_line),
                                                                     LineString(new_sub_line) )
        

    def _calculate_number_of_bends (self, lines):
        """
        This routine is calculating the total number of bends in the lines
        """
        
        nbr_bends = 0
        for line in lines:
            nbr_bends += len(line.bends)
            
        return nbr_bends
    
    def _reduce_bend(self, bend, line, bend_depth=_DEPTH_OFFSET_RATIO): 
        
        """This routine reduce one bend"""
    
        first_coord = line.coords_dual[bend.i]
        last_coord = line.coords_dual[bend.j]
            
        if (bend.base / bend.depth < 1.0):
            offset_factor = bend.base / bend.depth
        else:
            offset_factor = bend.depth / bend.base
        
        bend_depth = bend_depth * offset_factor 
               
            
        # Calculate the position of the middle of the bend base line
        base_mid_point = MA_Point(first_coord).mid_point(MA_Point(last_coord))
            
        # Offset the point in direction of the bend peak
        scaled_point = self.rescale_vector (base_mid_point, MA_Point(bend.peak_coords),  bend_depth)
            
        return [scaled_point.coords_dual[0]]

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
    
    def _replace_bend (self, i_bend, line, in_hand_line_coords):
        """
        This function replace  the vertices of the bend (excluding the first and last point) by  the replacement 
        line (excluding the first and last point which are the same as the first/last vertice of the bend)
        
        Parameters:
           i_bend: Bend number to process
           line: LineString object being processed
           in_hand_line_coords: In order to increase performance we are not rewriting the coordinates in
                                the Sahpely structure we keep them in that list
        """
    
        bend_to_check = line.bends[i_bend]
    
        # Determine the position of the bend limits
        first = bend_to_check.i
        last = bend_to_check.j    
    
        # Replace the coordinate with the new line in the temporary line coordinates (but we do not create a new object)
        in_hand_line_coords[0:] = list(bend_to_check.new_line.coords)
    
        if (self._are_bends_big([bend_to_check], line.min_adj_area)):
            self.stats.add_stats(_BIG)
        else:
            self.stats.add_stats(_SMALL)
            
        self.stats.add_stats(_ALGO)
        self.stats.add_stats(bend_to_check.type)
    
        return
        
    def rescale_vector(self, p1, p2, scale_factor):
        """This routine rescale the vector defined by the points P1 and P2 by a factor
        of SCALE_FACTOR
        
        P1: Point vector origin (first point)
        P2: Point vector to rescale (second point)
        scale_factor: factor to scale the vector (same for x and y)
        
        """ 
        
        x1 = p1.coords_dual[0][0]
        y1 = p1.coords_dual[0][1]
        x2 = p2.coords_dual[0][0]
        y2 = p2.coords_dual[0][1]
        
        vec_x = x2 - x1
        vec_y = y2 - y1
        
        vec_x = vec_x * scale_factor
        vec_y = vec_y * scale_factor
        
        x_out = vec_x + x1
        y_out = vec_y + y1
        
        return MA_Point([x_out, y_out])     
    
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
        base_line_p1 = line.coords_dual[bend.i]
        base_line_p2 = line.coords_dual[bend.j]
        
        for n in xrange(bend.i + 1, bend.j):
            
#            distance = bend_base_line.distance(Point(line.coords_dual[n]))
            distance = GenUtil.distance_line_point (base_line_p1, base_line_p2, line.coords_dual[n])
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
    

    def process(self, params):
        """Main routine for the Sherbend algorithm
        
        The algorithm will simplify the lines using the Sherbend algorithm. 
        It will iterate over the lines until there are no more bends to simplify.

        Keyword arguments:
            none

        Return value:
            geo_content: dataclass containing the output information

        """
        
#        print ("Start of sherbend  algorithm)")
#        print ("Parameter description:")
#        print ("  - Bend mode: {}".format(params.command.bend_mode))
#        print ("  - Simplicity constraint {}".format(params.command.simplicity))
#        print ("  - Crossing constraint {}".format(params.command.crossing))
#        print ("  - Adjacency constraint {}".format(params.command.adjacency))

#        # Check the feature's class and attributes
#        self.check_features()

        # Load the features into the spatial container
        self.s_container = self.load_features ()
     
        if (self.params.debug):
            #Test if print is needed before scanning the s_container for nothing... waste of time...
            nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING")) 
            nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT")) 
            GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
            GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
            
        self.add_line_attributes()
        
        if (self.params.multi_bend):
            nbr_step = 2
        else:
            nbr_step = 1
    
        for step in range(nbr_step):
        
            # The last step is always done with multi_bend to False.  We always process the last iterations to multi_bend to False
            # in order to process alternative bends correctly. An alternative bend is involved when a bend is in conflict and it tries
            # to simplify the bend just before or after the bend in conflict
            if (step == 1):
                self.params.multi_bend = False
            
            # Loop until no more lines are simplifiable
            line_simplified = True
            # Iterate until all the line are simplified or there are no more line have to be simplified
            while (line_simplified):
                # At each new iteration we create a new itearation
                self.stats.add_iteration()
                #Reset stats and error position
                self.stats.reset_stats_names([GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS]) # Reset error counter
                self.error_positions = []  # Reset error position we don't need the old one
                GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )  
                
                line_simplified = False
        
                line_simplified = self.manage_lines_simplification()
                                        
                GenUtil.print_debug (self.params, 'Number of bend simplified %s: ' %(self.stats.get_stats_name_count_iter(_ALGO)))
                GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )
                
                if (self.params.keep_iter_results):
                    self.iter_results.add_iteration()
                    self.iter_results.add_features(self.s_container.get_features(filter="feature_type==GenUtil.LINE_STRING"))
                    
        self.features = [feature for feature in self.s_container.get_features()] 
        
        GenUtil.print_debug (self.params, "End of %s" % (_ALGO))
