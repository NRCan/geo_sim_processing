#! /usr/bin/env python
# -*- coding: UTF-8 -*-

#####################################################################################################################################

"""
    This algorithm calculates the bends along the line. A bend is defined as all the line segments that turn in the same sense either
    clockwise or counter clokwise.
    
    Usage:
        import algo_bend

    Limits and constraints
        A straight line or a line with only 2 points has no bends
          
"""

#####################################################################################################################################


import math

from shapely.geometry import Point, LineString, Polygon
from shapely.prepared import prep

from lib_geobato import MA_LineString, InvalidParameterError, GenUtil,\
                         GenStatistics, IterationResults,Algorithm
                         
########################################################

# Public constant
_ALGO = "Bend"
_BENDS_DETECTED = "Bend detected"
_UNKNOWN = "Unknown"
        
class _Bend(object):
    """
    This class defines attributes and operations for bends
    
    Attributes: None
    """

    def __init__(self):
        """
        Initialize attributes:
        """
        
        self.i = None
        self.j = None

class Statistics(GenStatistics):
    """Class that contains the statistics for the talweg cohenrence algorithm
    
    Attributes
        stat_names: Name of the statistics for the TalwegStatistics class. These name are
                    used by the Statistics class

    """
    
    def __init__(self):
        """Initialize the attributes of an object of the class"""
        
        GenStatistics.__init__(self)
        self.stats_names = ((_BENDS_DETECTED, ))
        self.add_iteration()
        
    def get_stats (self, type=GenStatistics.SUMMARY):
        """Extract the current statistics and build  a list of string that forms the statistical message"
        
        Parameters:
            type: Give the form of statistics to extract. Can take 2 values.
                SUMMARY: Summary information
        
        """
        
        str_out = []
        str_out.append( "%s algorithm Statistics" %(_ALGO) )
        
        str_out.append( "-------------------------------------" )
        str_out.append("Bends detected: " + str(self.get_stats_name_count_total( _BENDS_DETECTED)))
        str_out.append( "-------------------------------------" )
        
        return str_out
    
class AlgoBends(Algorithm):
    """This is the main class for the bend algorithm 
    
    Attributes:
        - params: Parameters of the algorithm
        
    """

    def __init__(self, debug=False):
        """Initialize the attributes of an object of the class

        Parameters:
            debug: Flag to enable(TRUE)/disable(FLASE) for debug output 

        Return value:
            None

        """

        Algorithm.__init__(self)
        
        self.params = Parameters()  
        self.params.debug=debug       
        self.stats = Statistics()      
        

    def locate_bends(self): 
        """Calculates the position of each individual bends in the line
        
        The position of the bends are calculated according to the definition of the bencds 
        in the orginal paper Wang 1998.
        
        Parameters: None
        
        Return value: None  
        """
        
        for line in (line for line in self.features):
        
            line.i = None
            line.j = None
            line.bends = []
            nbr_pnt = len(line.coords_dual)
        
            if (nbr_pnt >= 3):
                
                # A first loop to determine the rotation sense of the first 
                i = 1
                orientation = 0.0
                # We loop until it is not a straight line
                while  (orientation == 0.0 and i < nbr_pnt-1):
                    orientation = self.direction (line.coords_dual[i-1], line.coords_dual[i], line.coords_dual[i+1])
                    i += 1
                    
                if (orientation != 0.0):
                    i = 1
                    last_bend_last_angle = 0
                    last_orientation = orientation
                    # Detect all the bends of the line
                    while (i < nbr_pnt-1):
                        orientation = self.direction (line.coords_dual[i-1], line.coords_dual[i], line.coords_dual[i+1])
                        if (orientation > 0.0 and last_orientation > 0.0):
                            i_last_angle = i
                        elif (orientation == 0.0):
                            pass
                        elif (orientation < 0.0 and last_orientation < 0.0):
                            i_last_angle = i
                        else:
                            # A new bend is detected and loaded
                            line.bends.append(_Bend())
                            self.stats.add_stats(_BENDS_DETECTED)
                            line.bends[-1].i = last_bend_last_angle
                            line.bends[-1].j = i
                            last_bend_last_angle = i_last_angle
                            i_last_angle = i
                            last_orientation = orientation
                        i += 1
                    
                    # Manage the last bend of the line
                    line.bends.append(_Bend())
                    self.stats.add_stats(_BENDS_DETECTED)
                    line.bends[-1].i = last_bend_last_angle
                    line.bends[-1].j = i
                            
                else:
                    # A straight is detected no bends are created
                    line.bends = []
            else:
                # A line with only 2 points will never have a bend
                line.bends = []

    def direction(self, p0, p1, p2):
        """ Calculate the type angle (clockwise or anticlockwise) of a line formed by 3 vertices using the dot product
        
        Parameters:
            p0, p1, p2: Three (x,y) coordinates tuple 
            
        Return value
            float the direction of the line an
                0: Straight line
                <0: Counter clockwise angle
                >0: Clockwise angle
                
        """
        
        return ((p0[0] - p1[0]) * (p2[1] - p1[1])) - ((p2[0] - p1[0]) * (p0[1] - p1[1]))
    
    def check_features(self):
        """
        Check if the features passed in parameters are of the good class type and have the good attributes
        
        Parameters: None
        
        Return value: None
        
        """
        
        # Check the MA_LineString
        class_type = MA_LineString
        properties_to_check = []
        for feature in self.features:
            GenUtil.check_feature_integrity(feature, MA_LineString, properties_to_check)
               
    def process(self):
        """Main routine for the Bend algorithm
        
        This method will detect the bends in list of lines.  A list of object of type if
        is created on each line even if it has no bends. A line with no bend will have an empty list
        
        Parameters: None
            
        """
                
        # Check the feature's class and attributes
        self.check_features()
        
        self.locate_bends()
