#!/usr/local/bin/python
# -=- encoding: utf-8 -=-

#####################################################################################################################################

"""
    This algorithm implement the classic Douglas Peucker line simplification with constraint.

    This algorithm simplify line using the classic Douglass Peucker algorithm.  It will disable line
    simplification if the line simplified cause one of the three constraint (simple line, line crossing, 
    sidedness) to be violated. If point feature are inputted they are used for constraint validation


    Usage:
        import algo_douglas_peucker

    Limits and constraints


"""
__revision__ = "--REVISION-- : $Id: algo_douglas_peucker.py 398 2011-07-26 15:02:51Z dpilon $"

#####################################################################################################################################

import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_genmetal import MA_Point, MA_LineString, MA_Polygon, InvalidParameterError, GenUtil,\
                         GenStatistics, PointErrorPosition, LineStringErrorPosition, IterationResults, Parameters,\
                         SpatialContainer, Algorithm

########################################################

#Internal key word constants
_ALGO = 'Douglas Peucker'


# If psyco is availiable, use it to speed up processing (2x)
try:
    import psyco
    psyco.full()
except ImportError:
    pass

class Statistics(GenStatistics):
    """Class that contains the statistics for the DP algorithm

    Attributes
        stat_names: Name of the statistics for the SpikeStatistics class. These name are
                    used by the Statistics class

    """

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        GenStatistics.__init__(self)
        self.stats_names = ((_ALGO, GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS))

    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"

        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
                DETAILED: Detailed information

        """

        str_out = []
        str_out.append( "%s algorithm Statistics" %(_ALGO))
        str_out.append( "-" * len(str_out[-1]))
        if (type == GenStatistics.DETAILED):
            for i in xrange((self.get_nbr_iteration())):
                str_out.append("Detailed statistics")
                str_out.append("Iteration # " + str(i))
                str_out.append("Vertices removed: " + str(self.get_stats_name_count_iter( _ALGO, i)))
                str_out.append( "--------" )
                str_out.append( "Conflicts:" )
                str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_iter( GenUtil.SIMPLE_LINE, i)))
                str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_iter( GenUtil.CROSSING_LINE, i)))
                str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_iter( GenUtil.SIDEDNESS, i)))
                str_out.append( "--------" )
        str_out.append( "Summary statistics" )
        str_out.append("Vertices removed: " + str(self.get_stats_name_count_total( _ALGO)))
        str_out.append( "--------" )
        str_out.append( "Conflicts:" )
        str_out.append( "    Simple Line   :  " + str(self.get_stats_name_count_total( GenUtil.SIMPLE_LINE)))
        str_out.append( "    Crossing Line :  " + str(self.get_stats_name_count_total( GenUtil.CROSSING_LINE)))
        str_out.append( "    Sidedness     :  " + str(self.get_stats_name_count_total( GenUtil.SIDEDNESS)))
        str_out.append( "--------" )
        str_out.append( "Number of iteration: " + str(self.get_nbr_iteration()) )

        return str_out




class AlgoDouglasPeucker(Algorithm):
    """
    This is the main class for the spike algorithm

    Attributes:
        - params: Parameters for the algorithm
    """

    def __init__(self, test_crossing_line=True,
                       test_simple_line=True, 
                       test_sidedness=True, 
                       keep_iter_results=False, 
                       debug=False):
        """Initialize the attributes of an object of the class DPAlgorithm

        Parameters:
            test_crossing_line: Flag to enable(TRUE)/disable(FLASE) checking the line crossing constraint
            test_simple_line: Flag to enable(TRUE)/disable(FLASE) checking the simple line constraint
            test_sidedness:  to enable(TRUE)/disable(FLASE) for checking the sidedness constraint
            keep_iter_results: Flag to enable(TRUE)/disable(FLASE) keeping the iterative results
            keep_iter_results: Flag for the iterative results
            debug: Flag to enable(TRUE)/disable(FLASE) for debug output 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()
        self.params.test_crossing_line=test_crossing_line
        self.params.test_simple_line=test_simple_line 
        self.params.test_sidedness=test_sidedness 
        self.params.keep_iter_results=keep_iter_results 
        self.params.debug=debug
        
        self.stats = Statistics()

    def _find_farthest_point(self, line, first, last):
        """
        Returns a tuple with the farthest point's index and it's distance of a line's section

        Parameters:
            - line: The line to process 
            - first: Index of the first point the line's section to test
            - last: Index of the last point the line's section to test

        """

        fartest_index = first
        fartest_dist = 0

        p_first = line.coords_dual[first]
        p_last = line.coords_dual[last]
        for i in range(first, last):
            dist = GenUtil.distance_line_point(p_first, p_last, line.coords_dual[i] )
            if  dist > fartest_dist:
                fartest_dist = dist
                fartest_index = i

        return (fartest_index, fartest_dist)

    def _process_closed_line (self, line):
        """
        Special simplification for closed line as we want to have at least 4 vertices to have a closed line with an area > 0

        Parameters:
            line: The line to process
            
        Return value:
            Set containing the index of the vertices to keep
            
        """
        
        first = 0
        last = len(line.coords_dual)-1
                
        # Calculate the farthest index from the first/last
        mid_index = None
        if ((first+1) <= last):
            mid_index, dummy = self._find_farthest_point(line, first, last)
            
        # Calculate the farthest point from the lower half of the closed area (if a lower half exist)
        lower_index = None
        lower_dist = None
        if ( mid_index is not None and ( (first+1) <= mid_index )):
            lower_index, lower_dist = self._find_farthest_point(line, first, mid_index)
                
        # Calculate the farthest point from the upper half of the closed area (if a upper half exist)
        upper_index = None
        upper_dist = None
        if ( mid_index is not None and ( (mid_index+1) <= last )):
            upper_index, upper_dist = self._find_farthest_point(line, mid_index, last)
                
        # From here determine the points to keep on the line
        index = set()                # Set of points to keep
        
        # Always keep the first and the last point
        index.update([first, last])  

        if (mid_index is not None):
            # Add the mid point
            index.update([mid_index])
            
            fourth_vertice = None
            if (lower_index is not None and upper_index is not None):
                if (lower_dist > upper_dist):
                    fourth_vertice = lower_index
                else:
                    fourth_vertice = upper_index
            else:
                if (lower_index is not None): fourth_vertice = lower_index
                if (upper_index is not None): fourth_vertice = upper_index
                
            if fourth_vertice is not None: index.update([fourth_vertice]) 
                        
        return (list(index))

    def reduce(self, line, pass_nbr):
        """
        This method is simplifying a line with the Douglas Peucker algorithm plus contraints checking if they are enabled.
        
        This method is checking the line differently the first time from the remaining time.  The idea behind it is only to 
        have a faster process. We assume that most of the lines will not have a problem so we check for problems (SIMPLE_LINE,
        CROSSING_LINE and SIDEDNESS) against the whole line for the first time. If there are some problems tthe next time we will
        check each sub portion of the line. This strategy is making a huge difference in time.
        
        Parameters:
            line: The line to process
            pass_nbr: The number of the pass. At their first pass we process the line as a whole. We do not
            check each sub portion of line simplification. For the other passes we check each sub portion of the line
            for constraints
            
        Return value:
            True: The line is simplified
            False: The line is not simplified

        """

        index = set()  # Set of the points to keep
        stack = []     # Stack to simulate the recursion
        line.simpliest = True
        first = 0
        last = len(line.coords_dual) - 1
        stack.append((first, last))

        while stack:
            (first, last) = stack.pop()
            if first + 1 < last:  # The segment to check has only 2 points
                add_point = True
                (farthest_index, farthest_dist) = self._find_farthest_point(line, first, last)
                if farthest_index != first:
                    if farthest_dist <= line.tolerance:
                    
                        if (pass_nbr != 0):
                            # For all the pass except the first one we check each sub line simplification individually
                            line_simple_line = LineString(line._coords_dual[:first+1] + line.coords_dual[last:])
                            new_segment_coords = [line.coords_dual[first], line.coords_dual[last]]
                            old_segment_coords = line.coords_dual[first:last+1]
                            line_crossing_line = LineString(new_segment_coords)
                            sidedness_polygon = GenUtil.calculate_sidedness_polygon ( LineString(old_segment_coords),
                                                                                      LineString(new_segment_coords) )
    
                            conflict_type = GenUtil.test_constraints (self, line._sif_id, line_simple_line, line_crossing_line,
                                                                      sidedness_polygon, line.s_container, line.s_container_key)
                        else:
                            # We check for conflict only at the end of the process so here we assume no conflict
                            conflict_type = None
                                                                   
                        if (conflict_type is not None):
                            line.simpliest = False # This line is not at it's simplest form since
                        else:                      # a constraint is blocking the simplification
                            index.update([first, last])
                            add_point = False

                    if add_point:
                        stack.append((first, farthest_index))
                        stack.append((farthest_index, last))
                else:
                    index.update([first, last])

            else:
                index.update([first, last])

        replacement_index = list(index)
        
        if (line.is_closed and (len(replacement_index) <= 3)):
       #   Closed line must have at least 4 vertices
            replacement_index = self._process_closed_line(line)

        # Check if the line has been simplified
        nbr_vertice_simplified =  len(line.coords_dual) - len(replacement_index) 
        if nbr_vertice_simplified==0:
            simplified =  False      # No change done (same quantity of coordinates)
            line.is_simplest = True  # The line is at its simplest form
        else:
            new_coords = [line.coords_dual[i] for i in sorted(replacement_index)]
            if (pass_nbr != 0):
                # If we process each sub modifification of the line inividually
                simplified =  True     # The line has been simplified
            else:
                # For the first iteration we process the line as a whole
                # Check for conglict
                line_simple_line = LineString(new_coords)
                new_segment_coords = new_coords
                old_segment_coords = line.coords_dual
                line_crossing_line = LineString(new_segment_coords)
                sidedness_polygon = GenUtil.calculate_sidedness_polygon ( LineString(old_segment_coords),
                                                                          LineString(new_segment_coords) )

                conflict_type = GenUtil.test_constraints (self, line._sif_id, line_simple_line, line_crossing_line,
                                                          sidedness_polygon, line.s_container, line.s_container_key)
                if (conflict_type is None):
                    simplified = True # The line was  simplified
                    line.is_simplest = True # If at the first pass the whole line as no conflict it is at its simplest form
                else:
                    simplified = False # The line was not simplified
                    
            if (simplified):
                for i in xrange(nbr_vertice_simplified): self.stats.add_stats(_ALGO)
                line.update_coords( new_coords )
                
        return simplified

    def _set_line_attributes (self):
        """This routine sets the attributes to the line

        The routine checks:
            - if the first/last vertices are the same the line is closed; otherwise it is open
            - if the line has less than 3 vertices it is at its simplest form; otherwise it is not

        Parameters: None

        Return value: None

        """

        # Select only the lines
        for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"):
            line.is_closed = False
            line.is_simplest = False
            if ( line.tolerance == 0. or len(line.coords_dual) <= 2):
                # A tolerance of 0 or a line with 2 vertices is at its simplest form
                line.is_simplest = True
                
            if (GenUtil.distance(line.coords_dual[0], line.coords_dual[-1]) <= GenUtil.ZERO):
                line.is_closed = True
                if (len(line.coords_dual) <= 4):
                    # A closed line with less than 4 points is at its simplest form
                    line.is_simplest = True
            else:
                line.is_closed = False

        return

    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the points
        class_type = MA_Point
        attributes_to_check = []
        for point in self.points:
            GenUtil.check_feature_integrity(point, class_type, attributes_to_check)
            point.attributes_to_clone = attributes_to_check
            
        # Check the line string
        class_type = MA_LineString
        attributes_to_check = ["tolerance"]
        for line_string in self.line_strings:
            GenUtil.check_feature_integrity(line_string, class_type, attributes_to_check)
            line_string.attributes_to_clone = attributes_to_check
    
    # def process(self):
    #
    #     GenUtil.print_debug (self.params, "Start of %s simplification algorithm)" %(_ALGO) )
    #     GenUtil.print_debug (self.params, "Parameter description:")
    #     GenUtil.print_debug (self.params, "  - Simple line constraint: %s" %(self.params.test_simple_line))
    #     GenUtil.print_debug (self.params, "  - Crossing line constraint: %s" %(self.params.test_crossing_line))
    #     GenUtil.print_debug (self.params, "  - Sidedness constraint: %s" %(self.params.test_sidedness))
    #     GenUtil.print_debug (self.params, "  - Keep iterative results: %s" %(self.params.keep_iter_results))
    #
    #     # Check the feature's class and attributes
    #     self.check_features()
    #
    #     # Load the shapely features into the spatial container
    #     self.s_container = self.load_features()
    #
    #     if (self.params.debug):
    #         #Test if print is needed before scanning the s_container for nothing... waste of time...
    #         nbr_lines = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
    #         nbr_points = len(self.s_container.get_features(filter="feature.feature_type==GenUtil.POINT"))
    #         GenUtil.print_debug (self.params, "Number of lines imported: %s"  %(nbr_lines) )
    #         GenUtil.print_debug (self.params, "Number of points imported: %s"  %(nbr_points) )
    #
    #     self._set_line_attributes()
    #
    #     line_simplified = True
    #     pass_nbr = 0
    #
    #     while (line_simplified or pass_nbr <= 1):
    #         # At each new iteration we create a new itearation
    #         self.stats.add_iteration()
    #         GenUtil.print_debug (self.params, 'Start of iteration # %s' %(self.stats.get_nbr_iteration()) )
    #
    #         line_simplified = False
    #
    #         for line in self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING and not feature.is_simplest"):
    #
    #             line_simplified = self.reduce(line, pass_nbr) or line_simplified
    #
    #         GenUtil.print_debug (self.params, 'Number of points removed %s: ' %(self.stats.get_stats_name_count_iter(_ALGO)))
    #         GenUtil.print_debug (self.params, 'End of iteration # %s' %(self.stats.get_nbr_iteration()) )
    #
    #         if (line_simplified):
    #             # If line_simplified is True there will be an other iteration and we reset the error counter
    #             self.stats.reset_stats_names([GenUtil.SIMPLE_LINE, GenUtil.CROSSING_LINE, GenUtil.SIDEDNESS]) # Reset error counter
    #             self.error_positions = []  # Reset error position we don't need the old one
    #
    #         if (self.params.keep_iter_results):
    #             self.iter_results.add_iteration()
    #             self.iter_results.add_features(self.s_container.get_features(filter="feature.feature_type==GenUtil.LINE_STRING"))
    #
    #         pass_nbr += 1
    #
    #     self.extract_features_out(self.s_container.get_features())
    #
    #     # Destruction of the cyclic references and garbage collection
    #     self.s_container.delete_container()
    #
    #     GenUtil.print_debug (self.params, "End of %s" % (_ALGO))

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
        self.load_features(self.geo_content.features)

        s_constraints = SpatialConstraints(s_container=self.s_container)

        self._manage_lines_simplification(s_constraints)

        out_features = []

        for feature in self.s_container.get_features():
            if feature.sb_geom_type == GenUtil.POINT:
                out_features.append(feature)
            elif feature.sb_geom_type == GenUtil.LINE_STRING:
                if feature.sb_original_type == GenUtil.LINE_STRING:
                    out_features.append(feature)
                else:
                    if feature.sb_original_type == GenUtil.POLYGON_EXTERIOR:
                        # The LineString was an exterior Polygon so reconstruct the original Polygon
                        interiors = [list(interior.coords) for interior in feature.sb_interiors]
                        polygon = Polygon(feature.coords, interiors)
                        polygon.sb_layer_name = feature.sb_layer_name
                        polygon.sb_properties = feature.sb_properties
                        out_features.append(polygon)
                    else:
                        pass  # Nothing to do with the holes here

        return out_features