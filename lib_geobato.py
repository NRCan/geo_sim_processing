#!/usr/bin/env python
# -=- encoding: utf-8 -=-

"""
General classes and utilities needed for the GENeralization MEta ALgorithm (GENMTEAL) tool

"""

from abc import ABCMeta
from abc import abstractmethod
import math
from collections import Iterable
from copy import deepcopy
from itertools import count

from shapely.geometry import Point, LineString, Polygon
from shapely.geometry.polygon import LinearRing
from shapely.ops import cascaded_union

from rtree import Rtree



class GenUtil:
    """This class defines a series of generic utility static method"""

    # Define some constants...
    POINT = 'Point'
    LINE_STRING = 'LineString'
    POLYGON = 'Polygon'
    POLYGON_EXTERIOR = 'PolygonExterior'
    POLYGON_INTERIOR = 'PolygonInterior'


    ANTI_CLOCKWISE = 1
    CLOCKWISE = -1
    CROSSING_LINE = 'Crossing line'
    SIDEDNESS = 'Sidedness'
    SIMPLE_LINE = 'Simple line'
    INVALID = 'Invalid'
    ZERO = 0.000001
    RADIAN = 'Radian'
    DEGREE = 'Degree'

    @staticmethod
    def distance (p1, p2):
        """Calculate the euclidian distance between 2 points

        *Parameters*:
            - p1: (x,y) tuple of the first coordinate
            - p2: (x,y) typle of the second coordinate

        *Returns*:
            - Distance between the 2 points (real)
            
        """

        return ( math.sqrt((p2[0]-p1[0])**2.0 + (p2[1]-p1[1])**2.0)   )
    
    @staticmethod
    def distance_line_point(p1, p2, p0):
        """Calculates the distance between a line (p1,p2) and a point p0
        
        Returns the shortest distance from the point to the line.
        If the intersecting point on the line is outside of the end points of the line segment, 
        it takes the distance to the nearest end point.
        
        *Parameters*:
            - p1: (x,y) tuple of the first coordinate of the line 
            - p2: (x,y) tuple of the second coordinate of the line
            - p0: (x,y) tuple of the point
    
        *Returns*:
            - Shortest distance from the point to the line
        """
        
        dist_p1_p2 = GenUtil.distance (p1, p2)
        if (dist_p1_p2 > GenUtil.ZERO):
            # p1 and p2 are not colinear
            x1,y1  = p1[0], p1[1]
            x2,y2  = p2[0], p2[1]
            x0, y0 = p0[0], p0[1]
            
            dot_product = (x0-x1)*(x2-x1) + (y0-y1)*(y2-y1)
            param = dot_product / dist_p1_p2**2.0
            
            if (param < 0):
                # The point is behind the first coordinate
                dist = GenUtil.distance(p1, p0)
            elif (param > 1):
                # The point is after the second coordinate
                dist = GenUtil.distance(p2, p0)
            else:
                # The point is between p1 and p2
                dist = abs((x2-x1)*(y1-y0) - (x1-x0)*(y2-y1)) / dist_p1_p2
        else:
            # p1 and p2 are colinear take the distance between  p0, p1
            dist = GenUtil.distance(p1, p0)
            
        return dist


    @staticmethod
    def print_debug (params, buffer):
        """
        This routine print the buffer id the params debug is True

        *Parameters*:
            - params: Parameters
            - buffer: String to print

        *Returns*: *None*

        """

        if params.debug:
            print (buffer)

        return


    @staticmethod
    def build_bounding_box(tolerance, coord):
        """Create and adjust a bounding box (xmin, ymin, xmax, ymax) 
        
        *Parameters*: 
            tolerance: A delta to add to the bounding box
            coord: (x,y) tuple
            
         *Returns*: 
            Tuple of the bounding box (xmin, ymin, xmax, ymax)  
            
        """
        
        xmin = coord[0] - tolerance
        ymin = coord[1] - tolerance
        xmax = coord[0] + tolerance
        ymax = coord[1] + tolerance
        
        return (xmin, ymin, xmax, ymax)    
    
    @staticmethod
    def make_iterable(iter_feature):
        """Test if the parameter is iterable; if not make it iterable by creating a tuple of one element

        *Parameters*:
            - iter_feature: Object to test if it's iterable

        *Returns*:
            - Iterable object

        """

        if not isinstance(iter_feature, Iterable):
            iter_feature = (iter_feature,)

        return iter_feature

    @staticmethod
    def create_LineString (lst_coords):
        """This method create a LineString from a list of coordinates

        If the number of coordinates is 1 or 0 it will create an Empty Geometry

        *Parameters*:
            - List of coordinate

        *Returns*
            - LineString object or Empty Geometry if no coordinates
        """

        if lst_coords is None or len(lst_coords) <= 1:
            lst_coords = []

        return LineString(lst_coords)

    @staticmethod
    def add_err_position (algorithm, id, line_error, err_type):
        """This method add an error in the error position list
        
         *Parameters*:
           - algorithm: Algorithm object to add the error
           - id: ID of the line being processed
           - line_error: LineString object containing the position of the error
           - err_type: Type of error to add in the statistics
           
         *Return*:
           - None 
        
        """
        
        # Increment statistics error count
        algorithm.stats.add_stats(err_type)
        
        # Log the position of the error 
        algorithm.error_positions.append(LineStringErrorPosition(id, err_type, line_error.coords))
    
    @staticmethod
    def test_constraints (algorithm, id, line_simple_line, line_crossing_line, polygon_sidedness, s_container, keys, 
                          add_err_position=True):
        """Test if the constraints are respected for a specific feature

        Each constraint can have a separate feature (line or polygon) to check.  If the feature corresponding to the feature
        to check is *None* the constraint is not checked. As soon as a constraint is violated the other constraints are not checked.

        *Parameters*:
           - algorithm: Algorithm object containing a params structure and information about the algorithm
           - id: ID of the line being processed
           - line_simple_line: The line feature to check against the simple line constraint
           - line_crossing_line: The line feature to check against the crossing line constraint
           - polygon_sideness: The polygon feature to check against the sidedness line constraint
           - s_container: the spatial container containing the feature and the index
           - key: The key of the feature being checked or the list of keys of features being checked.  When a search is done in the container
                  the keys are removed from the search
           - add_err_position: Flag to enable (True) or disable (False) the logging of the error in the stattistics and in the error position
        
        *Returns*:
           - conflict_type: Returns the *first detected* conflict type

        """
        conflict_type = None
        
        keys = GenUtil.make_iterable(keys)

        # Check if the algorithm have the "test_simple_line" attribute if no nothing to do
        if (hasattr(algorithm.params, 'test_simple_line')):
            if algorithm.params.test_simple_line and conflict_type is None:
                # Check the simple line constraint
                iter_line_simple_line = GenUtil.make_iterable(line_simple_line)
                for line_simple_line in iter_line_simple_line:
                    if (GenUtil.is_simple_line_constraint_violated (line_simple_line)):
                        conflict_type = GenUtil.SIMPLE_LINE
                        if (add_err_position):
                            algorithm.stats.add_stats(conflict_type)
                            algorithm.error_positions.append(LineStringErrorPosition(id, conflict_type, line_simple_line.coords))
                        break

        # Check if the algorithm have the "test_simple_line" attribute if no nothing to do
        if (hasattr(algorithm.params, 'test_crossing_line')):
            if (algorithm.params.test_crossing_line and conflict_type is None):
                # Check the crossing line constraint
                iter_line_crossing_line = GenUtil.make_iterable(line_crossing_line)
                for line_crossing_line in iter_line_crossing_line:
                    if ( GenUtil.is_crossing_line_constraint_violated (line_crossing_line, s_container, keys) ):
                        conflict_type = GenUtil.CROSSING_LINE
                        if (add_err_position):
                            algorithm.stats.add_stats(conflict_type)
                            algorithm.error_positions.append(LineStringErrorPosition(id, conflict_type, line_crossing_line.coords))
                        break

        # Check if the algorithm have the "test_sidedness" attribute if no nothing to do
        if (hasattr(algorithm.params, 'test_sidedness')):
            if (algorithm.params.test_sidedness and conflict_type is None):
                # check the sidedness constraint
                iter_polygon_sidedness = GenUtil.make_iterable(polygon_sidedness)
                for polygon_sidedness in iter_polygon_sidedness:
                    if ( GenUtil.is_sidedness_constraint_violated (polygon_sidedness, s_container, keys) ):
                        conflict_type = GenUtil.SIDEDNESS
                        if (add_err_position):
                            algorithm.stats.add_stats(conflict_type)
                            algorithm.error_positions.append(LineStringErrorPosition(id, conflict_type,
                                                                                     polygon_sidedness.exterior.coords))
                        break

        return conflict_type

    @staticmethod
    def is_line_closed (line):
        """Check if the first and last coordinate are the same or within a small distance
        
        *Parameters*: LineString or MA_LineString to verify for closeness
        
        *Returns:*
            Boolean: True the line is closed; False the line is open

        """
        if (isinstance(line,MA_LineString) and line.is_dual()):
            # Faster to access dual coordinate
            first = line.coords_dual[0]
            last =  line.coords_dual[-1]
        else:
            first = line.coords[0]
            last =  line.coords[-1]
            
        if (first[0] == last[0] and first[1] == last[1]):
            # First last coordinate are the same
            closed = True
        else:
            if (GenUtil.distance(first, last)  <GenUtil.ZERO):
                closed = True
            else:
                closed = False
                
        return closed

    @staticmethod
    def is_simple_line_constraint_violated (source_line):
        """
        Detect if the replacement line of the bend will cause the line to self intersect
        """

        # In some rare and complex situation Shapely can bug in that case wrap Shapely in a try catch
        # if there is a problem we suppose there in an error and we let the method going on...
        try:
            if (source_line is None  or source_line.is_simple):
                in_conflict = False
            else:
                in_conflict = True
        except:
            print ("genmetal_lib.py: Problem with shapely routine: is_simple in is_simple_line_constraint_violated")
            in_conflict = True
            
        return in_conflict

    @staticmethod
    def is_crossing_line_constraint_violated (source_line, s_container, keys):
        """Routine is checking line crossing constraint.

        The line crossing constraint is violated when the new line is intersecting any other points or line.
        The line crossing detection is done in 2 passes for speed. First a bounding box comparison and secondly
        a spatial filtering with Shapely.
        
        *Parameters*:
            - source_line: Line feature to check against the container for line crossing
            - s_container: Spatial container containing all the features
            - keys: the list of keys of features being check
            
        *Returns*:
            - Flag True/False indicating if the constraint is violated
             
        """

        if (source_line is None or source_line.is_empty):
            # Earlier processes might have created an empty line_string if so there is no line crossing problem
            target_features = []
        else:    
            target_features = s_container.get_features(bounds=source_line.bounds, remove_keys=keys)
    
        in_conflict = False
        if target_features:
            # In some rare and complex situation Shapely can bug in that case wrap Shapely in a try catch
            # if there is a problem we suppose there in an error and we let the method going on...
            try:
               features_intersects = filter(source_line.intersects, target_features)
               if (features_intersects):
                   # If all the features that intersects the line only touches the line
                   # (touches means connec at the extremity tha it's not in conflict
                   features_touches = filter(source_line.touches, features_intersects)
                   if (len(features_intersects) == len(features_touches)):
                       pass
                   else:
                       in_conflict = True
            except:
                print ("genmetal_lib.py: Problem with shapely routine: interscts in is_crossing_line_constraint_violated")
                in_conflict = True

        return in_conflict

    @staticmethod
    def is_sidedness_constraint_violated (source_polygon, s_container, keys):
        """
        This routine is checking sidedness constraint.
        The sidedness constraint is violated when a point or a line in completly inside the bend polygon
        we want to simplify.  The routine exit at the conflict encontered.
        
        *Parameters*:
            - source_polygon: Polygon feature to check against the container for line crossing
            - s_container: Spatial container containing all the features
            - keys: the list of keys of features being check
            
        *Returns*:
            - Flag True/False indicating if the constraint is violated 
            
        """
        
        if (source_polygon is None or source_polygon.is_empty):
            # Earlier processes might have created an empty polygon if so there is no sidedness problem
            target_features = []
        else:
            target_features = s_container.get_features(bounds=source_polygon.bounds, remove_keys=keys)
    
        in_conflict = False
        if target_features:
            # In some rare and complex situation Shapely can bug in that case wrap Shapely in a try catch
            # if there is a problem we suppose there in an error and we let the method going on...
            try:
                features_in_conflict = filter(source_polygon.contains, target_features)
                if (features_in_conflict):
                    in_conflict = True
            except:
                print ("genmetal_lib.py: Problem with shapely routine: contains in is_sidedness_constraint_violated")
                in_conflict = True


        return in_conflict

    @staticmethod
    def extract_bounds(lst_coords):

        list_x = [ x[0] for x in lst_coords ] # Create list of x coordinate
        list_y = [ y[1] for y in lst_coords ] # Create list of y coordinate

        return ( min(list_x), min(list_y), max(list_x), max(list_y) )

    @staticmethod
    def print_log (log, buffer):
        """
        This routine print the content to the log or the standard output
        """

        if log is None:
            print (buffer)
        else:
            log.log (buffer)

        return

    @staticmethod
    def rescale_vector(p1, p2, scale_factor):
        """Rescale the vector defined by the points P1 and P2 by a factor
        of SCALE_FACTOR

        Rescale the vector defined by the points P1 and P2 by a factor
        of SCALE_FACTOR

        *Parameters*:
            - P1: First coordinate of the vector. Tuple (x,y)
            - P2: Last  coordinate vector to rescale (second point)
            - scale_factor: factor to scale the vector (same for x and y)
            
        *Returns*: *TBA*

        """

        x1 = p1[0]
        y1 = p1[1]
        x2 = p2[0]
        y2 = p2[1]

        vec_x = x2 - x1
        vec_y = y2 - y1

        vec_x = vec_x * scale_factor
        vec_y = vec_y * scale_factor

        x_out = vec_x + x1
        y_out = vec_y + y1

        return (x_out, y_out)

    @staticmethod
    def angle_vector (p1, p2, p3, type=DEGREE):
        """Calculate the angle formed by the vector p1-p2 and p2-p3

        *Parameters*:
            - p1: (x,y) tuple of coordinates
            - p2: (x,y) tuple of coordinates
            - p3: (x,y) tuple of coordinates
            - type: Angle type DEGREE or ANGLE

        *Returns*:
            - The angle between the vector p1-p2 and p2-p3 (float)

        """

        a = (p2[0]-p1[0],p2[1]-p1[1])
        b = (p2[0]-p3[0],p2[1]-p3[1])
        len_a = (a[0]**2. + a[1]**2.)**.5
        len_b = (b[0]**2. + b[1]**2.)**.5

        dot_p = a[0]*b[0] + a[1]*b[1]

        # If P1 == P2 or P2 == P3 ===> angle is 180.
        if (len_a*len_b != 0.0):
            value = dot_p/(len_a*len_b)
            if value >= 1.0:  value = 1.0
            if value <= -1.0: value = -1.0
        else:
            value = -1.0

        theta = math.acos(value)
        
        if (type == GenUtil.DEGREE):
            theta = math.degrees(theta) 

        return theta 

    @staticmethod
    def compute_angle (p1, p2, p3, type=DEGREE):
        """
        Function to calculate angle between two vectors.
        """

        return (GenUtil.angle_vector(p1, p2, p3, type))

    @staticmethod
    def locate_bends1(lst_coords):

        bends = []

        if len(lst_coords) <= 2:
            # There is no bend in a line with two coordinates
            bend_direction = 0.0
        else:
            # check if the line is a straight line
            for i in range(1, len(lst_coords) - 1):
                bend_direction = GenUtil.direction(lst_coords[i - 1], lst_coords[i], lst_coords[i + 1])
                if bend_direction != 0.0:
                    break

        if bend_direction != 0.0:

            before_inflexion = 0  # Position of the index beofre the inflexion
            #            if GenUtil.is_line_closed(lst_coords):
            #                line_is_cloed = True
            #            else:
            #                line_is_closed = False

            last_bend_direction = GenUtil.direction(lst_coords[0], lst_coords[1], lst_coords[2])

            for i in range(1, len(lst_coords) - 1):
                bend_direction = GenUtil.direction(lst_coords[i - 1], lst_coords[i], lst_coords[i + 1])
                if bend_direction == 0.0:
                    # Nothing to do it's a straight line
                    pass
                elif last_bend_direction * bend_direction > 0.0:
                    # Nothing to do the bend is in the same orientation as the last bend
                    pass
                else:
                    # Change of bend direction; a bend is detected and created
                    bends.append((before_inflexion, i))
                    before_inflexion = i - 1
                    last_bend_direction = bend_direction

            # Manage the last bend of the line
            bends.append((before_inflexion, len(lst_coords) - 1))

        return bends


    @staticmethod
    def locate_bends_closed_line(lst_coords):

        bends = []

        if len(lst_coords) <= 3:
            # There is no bend in a closed line with 3 coordinates or less
            last_vertice = -1
        else:
            # Locate the position of the first change in inflexion in the line
            start_vertice = -1
            lst_tmp = list(lst_coords)
            del lst_tmp[-1]  # Delete the last coordinates because it's the same as the first one
            nbr_tmp = len(lst_tmp)
            last_bend_direction = GenUtil.direction(lst_tmp[-1%nbr_tmp], lst_tmp[0], lst_tmp[1])
            for i in range(1, nbr_tmp):
                bend_direction = GenUtil.direction(lst_tmp[(i-1)%nbr_tmp], lst_tmp[i], lst_tmp[(i+1)%nbr_tmp])
                if bend_direction == 0.0:
                    # It's a straight line pass to the next vertice
                    pass
                elif last_bend_direction*bend_direction >= 0.0:
                    # The bend is in the same direction pass to the next vertice
                    pass
                else:
                    # The is a change in bend direction
                    start_vertice = i-1
                    break

            if start_vertice == -1:
                # There is no bend in this closed line
                pass
            else:
                before_inflexion = start_vertice
                start_vertice += 1
                j = start_vertice
                last_bend_direction = GenUtil.direction(lst_tmp[(j-1)%nbr_tmp], lst_tmp[(j)%nbr_tmp], lst_tmp[(j+1)%nbr_tmp] )
                for i in range(nbr_tmp):
                    j = i+start_vertice
                    bend_direction = GenUtil.direction(lst_tmp[(j-1)%nbr_tmp], lst_tmp[(j)%nbr_tmp],
                                                       lst_tmp[(j+1)%nbr_tmp])
                    if last_bend_direction == 0.0:
                        # It's a staight line ; pass
                        pass
                    elif last_bend_direction*bend_direction >= 0:
                        # The bend is in the same directio; pass
                        pass
                    else:
                        # Change of bend direction; a bend is detected and created
                        bends.append((before_inflexion, j))
                        before_inflexion = j
                        last_bend_direction = bend_direction

                # Create a last bend because it' s a closed line
                bends.append((before_inflexion, start_vertice-1))

        return bends




    @staticmethod
    def locate_bends(lst_coords):
        """Calculates the position of each individual bends in a line

        The position of the bends are calculated according to the definition of the bencds
        in the orginal paper Wang 1998.

        Keyword definition
            lst_coords -- list of (x,y) tuple forming

        Return value: Bend
        """

        nbr_coords = len(lst_coords)
        bends = []
        if (nbr_coords >= 3):

            # A first loop to determine the rotation sense of the first
            i = 1
            orientation = 0.0
            # We loop until it is not a straight line
            while (orientation == 0.0 and i < nbr_coords - 1):
                orientation = GenUtil.direction(lst_coords[i - 1], lst_coords[i], lst_coords[i + 1])
                i += 1

            if (orientation != 0.0):
                i = 1
                last_bend_last_angle = 0
                last_orientation = orientation
                # Detect all the bends of the line
                while (i < nbr_coords - 1):
                    orientation = GenUtil.direction(lst_coords[i - 1], lst_coords[i], lst_coords[i + 1])
                    if (orientation > 0.0 and last_orientation > 0.0):
                        i_last_angle = i
                    elif (orientation == 0.0):
                        pass
                    elif (orientation < 0.0 and last_orientation < 0.0):
                        i_last_angle = i
                    else:
                        # A new bend is detected and created
                        bends.append((last_bend_last_angle, i))
                        last_bend_last_angle = i_last_angle
                        i_last_angle = i
                        last_orientation = orientation
                    i += 1

                # Manage the last bend of the line
                bends.append((last_bend_last_angle, i))

            else:
                # A straight is detected no bends are created
                pass
        else:
            # A line with only 2 points will never have a bend
            pass

        return bends

    @staticmethod
    def direction(p0, p1, p2):
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

    @staticmethod
    def compute_angles(coords_list):
        """
        Compute angle for each line vertice except first and last ones.
        Both ends receive a value of 360.0
        """

        angles = []
        angles.append(360.0)             # replace the uncomputable first vertice (index 0)

        cnt = len(coords_list)
        for v in range (1, cnt-1):       # '1' and 'cnt-1' to 'forget' first and last vertice
            angle = GenUtil.compute_angle(coords_list[v-1],coords_list[v],coords_list[v+1])
            angles.append(angle)

        angles.append(360.0)             # replace the uncomputable last vertice

        return angles

    @staticmethod
    def compute_flat_angles(coords_list):
        """Compute flat angle for each vertice in the list except first and last vertices."""

        flat_angles = []
        cnt = len(coords_list)

        # Get angle of each line vertice
        angles = GenUtil.compute_angles(coords_list)

        for i in range (0, cnt):
            flat_angles.append(abs(abs(angles[i])-180))

        return flat_angles

    @staticmethod
    def find_flatter_vertice(coords_list, use_group=True):
        """Find the vertice with the flatter angle."""

        cnt = len(coords_list)
        mean_angles = []

        values = GenUtil.compute_flat_angles(coords_list)

        if ((cnt < 4) or (use_group is False)):
            mean_angles = values

        if (cnt in [4,5]):
            u_function = "(values[i-1]*0.25)+(values[i]*0.5)+(values[i+1]*0.25)"
            mean_angles.append(180)                     # for first vertice
            for i in range(1,cnt-1):                    # '1' because 1 vertice on each side
                mean_angles.append(eval(u_function))
            mean_angles.append(180)                     # for penultimate vertice

        if (cnt in [6,7]):
            u_function = "(values[i-2]*0.05)+(values[i-1]*0.2)+(values[i]*0.5)+(values[i+1]*0.2)+(values[i+2]*0.05)"
            mean_angles.append(180)                     # for first vertice
            mean_angles.append(180)                     # for second vertice
            for i in range(2,cnt-2):                    # '2' because 2 vertices on each side
                mean_angles.append(eval(u_function))
            mean_angles.append(180)                     # for pre-penultimate vertice
            mean_angles.append(180)                     # for penultimate vertice

        if (cnt > 7):
            u_function = "(values[i-3]*0.05)+(values[i-2]*0.05)+(values[i-1]*0.2)+(values[i]*0.4)+(values[i+1]*0.2)+(values[i+2]*0.05)+(values[i+3]*0.05)"
            mean_angles.append(180)                     # for first vertice
            mean_angles.append(180)                     # for second vertice
            mean_angles.append(180)                     # for third vertice
            for i in range(3,cnt-3):                    # '3' because 3 vertices on each side
                mean_angles.append(eval(u_function))
            mean_angles.append(180)                     # for pre-penultimate vertice
            mean_angles.append(180)                     # for penultimate vertice
            mean_angles.append(180)                     # for last vertice

        # return the smallest flat angle and its list index
        return min(zip(mean_angles, count()))

    @staticmethod
    def find_bounds(coords_list, nbVertices):
        """
        Determine boundaries based on the required number of vertices.
        Bounds are 0 based. Both ends are included in result.
        """

        cnt = len(coords_list)
        i = 0
        bounds = []

        # Append first bound (vertice 0)
        bounds.append(0)

        # Iterate through vertices and calculate bounds
        inProgress = True
        while (inProgress):
            i += 1
            r=(i*nbVertices)-1

            # Check if end of vertices (end of line)
            if (r >= cnt-1):
                bounds.append(cnt-1)
                inProgress = False
            else:
                bounds.append(r)

        return bounds
    
    @staticmethod
    def translate_coords (coords, x_delta, y_delta):
        """Translate the list of tuple of coordinates
        
        *Parameters*:
            - coords: List of tuple to translate
            - x_delta: x translation factor
            - y_delta: y translation factor
            
        *Returns*: 
            - List of translated coordinates
        
        """
        
        new_coords = []
        for i in xrange(len(coords)):
            new_coords.append((coords[i][0]+x_delta, coords[i][1]+y_delta))
            
        return new_coords
    
    @staticmethod
    def cut_line_distance(line, distance, copy_att=False):
        """Cut the line at specific distance along the line
    
        *Parameters*:
            - line: LineString to cut
            - distance: distance in grounf unit fom the start of the line where to cut the line
            - copy_att: Flag to enable (True) or disable (False) the copy of the attributes in
                        in ma_properties for MA_LineString only
        
        *Returns*:
            - List of LineString or MA_LineString
    """
    
        if (isinstance(line, MA_LineString)):
            if line.is_dual():
                coords = line.coords_dual
            else:
                coords = list(line.coords)
        elif (isinstance(line, LineString)):
            coords = list(line.coords)
            if (copy_att):
                raise InternalError ("Cannot copy attributes on LineString only on MA_LineString")
        else:
            raise InternalError ("Cannot cut a %s feature" %(line.geom_type))
        
        if distance <= 0.0 or distance >= line.length:
            lst_coords_out = [coords]
        else:
            lst_coords_out = [coords]  # This line will prevent errors on non simple lines
            last_i = len(coords)-1
            for i, p in enumerate(coords):
                if (i < last_i):
                    pd = line.project(Point(p))
                else:
                    pd = line.length
                if pd == distance:
                    lst_coords_out =  [coords[:i+1], coords[i:] ]
                if pd > distance:
                    cp = line.interpolate(distance)
                    lst_coords_out =  [list(coords[:i]) + [(cp.x, cp.y)], [(cp.x, cp.y)] + list(coords[i:])]
                    break
        
        lines_out = []
        for coords_out in lst_coords_out:
            if isinstance(line, MA_LineString):
                line_out =  MA_LineString(coords_out)
            else:
                line_out =  LineString(coords_out)
            if (copy_att):
                # Copy the attributes
                line_out.ma_properties = deepcopy(line.ma_properties)
            lines_out.append(line_out)
                                    
        return lines_out
        
    
    @staticmethod        
    def validate(features, valid_types, raise_exception=False):
        """Validate a feature type against a list of expected valid types.
        
        *Parameters*:
            - features : single feature or list of features to validate
            - valid_types : (list) list of valid types
            - raise_exception : (boolean) flag to raise an exception on first invalid feature
            
        *Returns*:
            a list of boolean(s) giving, in order, the validity status of each feature
            
        *Examples*:
            - validate(LineString(), [Point, MA_Point]) returns [False]
            - validate(MA_Point(0,0), [Point, MA_Point]) returns [True]
            - validate(3, [int]) returns [True]
            - validate([3,'hello'], [Point, int]) returns [True, False]
        
        """
        
        valid_status = []
        
        if (type(features) is not list and type(features) is not tuple):
            features = list([features])
            
        for feature in features:
            feature_type = type(feature)     
            valid = (feature_type in (valid_types))
            
            if ((not valid) and raise_exception):
                raise Exception ("Invalid feature type : " + str(type(feature)))
            else:
                    valid_status.append(valid)
                    
        return valid_status
    

    @staticmethod        
    def rotate_coords (coords, theta):
        """Rotate the list of tuple of coordinates by theta
        
        *Parameters*:
            - coords: List of tuple to rotate
            - theta: Rotation factor/angle in raduans
            
        *Returns*: 
            - List of rotated coordinates
        
        """
        
        new_coords = []
        for i in xrange(len(coords)):
            x = coords[i][0]
            y = coords[i][1]
            x_angle = x*math.cos(theta) - y*math.sin(theta)
            y_angle = x*math.sin(theta) + y*math.cos(theta)
            new_coords.append((x_angle,y_angle))
        
        return new_coords
               
    @staticmethod
    def shear_coords (coords, shear_x, shear_y):
        """Shear the list of tuple of coodinates by a shear factor
        
        *Parameters*: 
            - coords: List of tuple to shear
            - shear_x: shear factor in x axis
            - shear_y: shear factir in y axis
            
        *Returns*: 
            - List of sheared coordinates
            
        """
        
        new_coords = []
        for i in xrange(len(coords)):
            x = coords[i][0]
            y = coords[i][1]
            x = x + y*shear_x
            y = y + x*shear_y
            new_coords.append((x,y))
                          
        return new_coords
                    
    @staticmethod
    def polygonize(feature):
        """ Transform the feature into a polygon feature
        
        If the input feature is a line and the line is not closed the line is closed automatically when the polygon is created.
        If the feature do not form a valid polygon a lien with 2 vertices or a polygon with an area of 0 than an empty polygon
        is created instead.
        If a non empty polygon is passed and it is not a simple polygon we try to clean it (pol.buffer(0))
        
        Extra tests have been inserted into the code for performance issue.  Polygon manipulation is time consuming
        
        *Parameters*:
            - feature: Feature to transform into a Polygon. The feature can be of type LineString or Polygon
            
        *Returns*:
            - The feature transformed into a polygon
            
        """
        
        if (feature.is_empty):
            # whatever the geometry it becomes empty geometry
            pol = Polygon()
        
        else:
            
            feature_type = feature.geom_type
            if (not feature_type in ("Point", "LineString", "Polygon", "MultiPolygon") ):
                raise InternalError ("Feature type not supported: " + feature.geom_type)        
            elif (feature_type == "Point"):
                    pol = Polygon()
            elif (feature.geom_type == "LineString"):
                try:
                    if (len(feature.coords) <= 2):
                        # If there is only 2 coordinates or less it is not forming a polygon for sure
                        pol = Polygon()
                    else:
                        # Check if the line is a straight by comparing the length of the line against the length of the first/last coordinate
                        if (feature.length - GenUtil.distance(feature.coords[0],feature.coords[-1]) > GenUtil.ZERO):
                            pol = Polygon(list(feature.coords))
                        else:
                            # If the line is almost a straight don't try to create a polygon it will create problem else where
                            pol = Polygon()
                except:
                    # If the line is composed of less than 3 coordinate it's a Shapely Exception we create an empty polygon 
                    pol = Polygon()
            else:
                # It's already a polygon
                pol = feature
                
        if (not pol.is_empty and not pol.is_valid):
            # If the polygon is not simple it is probably because there is a self intersection in the polygon
            # pol.buffer(0.0) will clean the self intersecting polygon
            try:
                pol = pol.buffer(0)
                # Check for almost 0 area polygon
                if (pol.area <= GenUtil.ZERO):
                    # The polygon is to small. Should be an empty polygon.. so create an empty polygon instead
                    pol = Polygon()
            except ValueError: 
                pol = Polygon()
        
        return pol
    
    @staticmethod
    def calculate_sidedness_polygon(feature_a, feature_b):
        """Calculate the sidedness region between 2 features.
        
        If the features A and/or B are LineString, the feature are transformed into Polygon.  The polygon are than cleaned
        to resolve self touching or self crossing.  The symmetric difference is applied on the feature in order to create 
        the sidedness regions
            
        *Parameters*:
            - feature_a: First feature
            - feature_b: Second feature
                
        *Returns*:
            - Polygon or multipolygon representing the symmetric difference between the 2 features
            - *None* if there is an error calculating the symmetric difference
                
            """
    
        
        # Transform the features into polygons
        pol_a = GenUtil.polygonize(feature_a)
        pol_b = GenUtil.polygonize(feature_b)
        
        try:
            sym_dif = pol_a.symmetric_difference(pol_b)
        except:
            # This exception is very rare and it usually caused by the use of a almost 0 or empty polygon so we create an empty polygon 
            sym_dif = Polygon()
            
        return sym_dif
    
    @staticmethod
    def mid_point(p1, p2):
        """Return a point in the middle of the 2 points
        
        *Parameters*: 
            - p1: x,y tuple for the first point
            - p2: x,y tuple for the second point
            
        *Returns*:
            - x,y tuple of the milddle point
        """

        x = (p1[0] + p2[0]) / 2.
        y = (p1[1] + p2[1]) / 2.
        
        return (x,y)
    
    @staticmethod
    def check_feature_integrity(feature, class_type, attributes):
        """Check if the feature are of the good class and have the good attribute in the property ma_properties
        
        *Parameters*: 
            - feature to check
            - class_type: Class type of the features; if None the class is not checked
            - attributes: List of attributes name that the feature must have
            
        *Returns*: *None*
        
        """
        
        # Check the class type
        if class_type is not None:
            if not isinstance(feature, class_type):
                raise InternalError ("The geometry instance of type: %s" %(str(class_type)))
            
        # Check the attributes
        for attribute in attributes:
            if not feature.ma_properties.has_key(attribute):
                raise InternalError("The geometry don't have the property/attribute: %s" %(attribute))
            
    @staticmethod
    def add_vertex(line, norm_lr_pos, tolerance=0., s_container=None):
        """This method insert a vertex on a line based on a specific normalized linear reference position along the line.
        
        If the vertex to add is within the tolerance of an existing vertex the vertex is not added.
        
        *Parameters*:
            - line: MA_LineString object inside which we have to insert a vertex. If the object is not 
                    of type MA_LineString an exception is raised.
            - norm_lr_pos: Normalized linear reference position to create a vertex along the line. The domain 
                           value is ]0..1[
            - tolerance: If the vertex to insert is within this tolerance of actual vertex along the line, the vertex is not added.
                         If tolerance is none a vertex is always added.
            - s_container: The container containing the line or None if the line is not contained in a container.
        
        *Returns*:
            - False: No vertex added
            - True: A vertex was added
        
        """
        
        vertex_added = False
        # It is mandatory to work on a MA_LineString object
        if isinstance(line, MA_LineString):

            if ( 0. < norm_lr_pos < 1. ):            
                for i in reversed(xrange(len(line.coords_dual)-1)):
                    # Starting from the penultimate vertex it's looping back until it finds a
                    # linear reference smaller than the one we want 
                    coord_lr =  line.project(Point(line.coords_dual[i]), normalized = True)
                    if (coord_lr <= norm_lr_pos):
                        break
                    
                i += 1 # Point to the next coordinate
                # The vertex to add is between i and i+1
                point = line.interpolate(norm_lr_pos, normalized=True)
                coords = list(point.coords)
                new_coord = coords[0]
                if (coord_lr != norm_lr_pos): 
                    if ( GenUtil.distance(new_coord, line.coords_dual[i]) <= tolerance or
                         GenUtil.distance(new_coord, line.coords_dual[i-1]) <= tolerance ):
                        # The vertex to add is within the tolerance of an actual vertex of the line... nothing to do
                        pass
                    else:
                        # Add a new vertex on the line
                        line_coords = list(line.coords_dual)
                        line_coords = line_coords[0:i] + [new_coord] + line_coords[i:]
                        # Update the coordinate of the line
                        line.update_coords(line_coords, s_container)
                        vertex_added = True
                else:
                    # there is already a vertice on the line where we want to add a vertice
                    pass
            else:
                # The vertice to add is either on the first vertice or on the last vertice... nothing to do
                pass
        else:
            # We can only add a vertice on a MA_LineString line
            raise InternalError ("We can only add a vertex on a MA_LineString object")
        
        return vertex_added
    
    @staticmethod
    def warp_line(line_to_warp, warp_infos, propagation=0):
        """This method will deform or warp a list of vertice according to the constraint contained into the warp information. 

        *Parameters*:
            - line_to_warp: MA_LineString object to warp
            - warp_infos: List of tuple containing the information to warp the LineString. 
                          Each tuple contains three information (index, delta_x, delta_y)
                            - index: Vertex index number in the LineString (>=0 and < the number of vertex in the line string)
                            - delta_x: Shift to apply in the x coordinate 
                            - delta_y: Shift to apply in the y coordinate
                          There is always at least two tuple of information. Only one tuple of information is a fatal error.
                          The warping is done from the first warp_information to the second warp information, than from the second warp
                          information to the third warp information and from the third to the fourth and so on until the last warp information
                          Example the following warp information [(3,0,0),(7,20,-15),(10,10,-30),(15,0,0)] has to be interpret has follow
                          First: Vertex:3 is not moved (0,0) and vertex 7 is moved 20 in x and -15 in y; vertex 4, 5 and 6  are moved
                                 and the warp propagation (from 0 to 20 in x and from 0 to -15 in y) is distributed between the vertice.
                          Second: Vertex 10 is moved (10 in x and -30 in y); vertex 8, 9 are moved
                                  and the warp propagation (from 20 to 10 in x and from -15 to -30 in y) is distributed between the vertice.
                          Third: Vertex 15 is not moved; vertex 11, 12, 13, 14 are moved and the warp propagation between these vertice
                                 (from 10 to 0 in x and from -30 to 0 in y) is distributed between these vertice.
                                       
            - propagation: Type of propagation between 2 warp information. The possible values are
                           0: Linear propagation
        
        *Returns*:
            - MA_LineString object warped
        
        """
        
        len_lst_coords = len(line_to_warp.coords_dual)
        lst_delta_xy = [(0.,0.) for i in range(len_lst_coords)]
        
        # Create the list of the cumulative distance of the coordinate
        lst_cumul_dist = [0.]
        coord_0 = line_to_warp.coords_dual[0]
        for i in xrange(1, len_lst_coords):
            distance = ((line_to_warp.coords_dual[i][0]-line_to_warp.coords_dual[i-1][0])**2.0 + 
                        (line_to_warp.coords_dual[i][1]-line_to_warp.coords_dual[i-1][1])**2.0)**.5 
            lst_cumul_dist.append(lst_cumul_dist[i-1]+distance)
            
        # Set the value for the first warp information
        start_i = warp_infos[0][0]
        start_delta_x = warp_infos[0][1]
        start_delta_y = warp_infos[0][2]
        start_dist =  lst_cumul_dist[start_i]
        lst_delta_xy[start_i] = (start_delta_x,start_delta_y)
        
        # Loop over each warp information
        for warp_info in warp_infos[1:]:
            # Initialize the end information
            end_i = warp_info[0]
            end_delta_x = warp_info[1]
            end_delta_y = warp_info[2]
            end_dist =  lst_cumul_dist[end_i]

            # Validate some information
            if (start_i >= end_i):
                raise InternalCorruptionError
            if (end_i > len_lst_coords-1):
                raise InternalCorruptionError 
            
            # Loop over each vertice between the start point and end point in order to propagate the displacement
            for i in range(start_i+1, end_i+1):
                current_dist = lst_cumul_dist[i]
                # A test to avoid division by zero...
                if (end_dist-start_dist >= GenUtil.ZERO):
                    percent_dist =  (current_dist-start_dist) / (end_dist-start_dist)
                else:
                    percent_dist = 1.0
                delta_x = start_delta_x +(end_delta_x - start_delta_x)*percent_dist
                delta_y = start_delta_y +(end_delta_y - start_delta_y)*percent_dist
                lst_delta_xy[i] = (delta_x,delta_y)
                
            # For the next loop the end info becomes the start info
            start_i = end_i
            start_delta_x = end_delta_x
            start_delta_y = end_delta_y
            start_i = end_i
            start_dist = end_dist
            
                
        # Apply the list of delta to each coordinates
        lst_coords = []
        for i in xrange(len_lst_coords):
            coords = line_to_warp.coords_dual[i]
            delta_xy = lst_delta_xy[i]
            new_coords = (coords[0]+delta_xy[0], coords[1]+delta_xy[1])
            lst_coords.append(new_coords)
            
        return MA_LineString(lst_coords)
    

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
            if GenUtil.distance(self.coords[0], self.coords[-1]) <= GenUtil.ZERO:
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
        print ("Need to update the spatial container...")
        LineString.coords.__set__(self, coords)
        if self.fast_access:
            self.__lst_coords = list(super().coords)

    def remove_colinear_vertex(self):
        """This method remove the colinear verxtex in the line string. Also handles closed line"""

        # Detect the position of the colinear vertex
        vertex_to_del = [i for i, angle in (enumerate(self.vertex_to_del)) if angle == 0.]
        if len(vertex_to_del) >= 1:
            # Delete the colinear vertex
            lst_coords = list(self.coords)
            for i in reverse(vertex_to_del):
                del(lst_coords[i])
            self.coords = lst_coords




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
        Poin.coords.__set__(self, coords)
        if self.fast_access:
            self.__lst_coords = list(super().coords)


class Holder(object):
    """Generic class creator
    
    Creates one Holder object with as many properties there are key words key
    
    *Parameters*: As many parameters as needed
    
    *Returns*:
        Holder object with as many properties as there was parameters
    
    """ 
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class GenException (Exception):
    """
    This is the base exception class for genmetal algorithms
    """

    def __init__(self, *arguments, **keywords):
        Exception.__init__(self, *arguments, **keywords)

class InvalidParameterError (GenException):
    """
    This exception is raised when an parameter is set to an invalid value (or not set)
    """

    def __init__(self, *param_names):
        """
        Initialise an Invalid Parameter Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GenException.__init__(self, *param_names)

class MetricTestError (GenException):
    """
    This exception is raised when a metric test failed
    """

    def __init__(self, *param_names):

        """
        Initialise an Invalid Parameter Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GenException.__init__(self, *param_names)

class InvalidGeometryError (GenException):
    """
    This exception is raised when an invalid geometry is passed
    """

    def __init__(self, *param_names):

        """
        Initialise an Invalid Geometry Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GenException.__init__(self, *param_names)

class InternalError (GenException):
    """
    This exception is raised when an internal error as occurred
    """

    def __init__(self, *param_names):

        """
        Initialise an Invalid Geometry Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GenException.__init__(self, *param_names)

#class Parameters(object):
#    """Class to hold parameters.
#
#    The parameters are defines on the fly by each algorithm that instantiate a
#    Parameters object
#
#    """
#
#    def __init__(self):
#        """Create an object of type Parameters
#
#        *Parameters*: *None*
#
#        *Returns*: *None*
#
#        """
#
#        pass

class MA_Point(Point):

    """
    The MA_Point class extends the Point of the shapely library by adding
    new methods and attributes that enables the creation of Meta algorithm.
    Here are some possibilities of this class:
     - MA_Point can be stored  in an instance of a SpatialContainer class
     - MA_Point can have than coordinates updated through the update_coords
       method which also updates the spatial index if used
     - MA_Point disable the capacity to update a coordinate through the property coords;
       you must always use the update_coords method to update coordinates
     - MA_Point can have a copy of the coordinates in an the attribute coords_dual;
       this copy of the attributes allows a faster access to the coordinates than
       when accessing them through Shapely coords attributes. This can be an advantage when
       you have to loop over and over each coordinates of your features.

    """

    def __init__(self, coords, dual=True):
        """Initialize the properties of a MA_Point object.  The properties of an
        instance are:
          -feature_type: The type of feature
          - _coords_dual: duplicate of the coordinates
          - ma_properties: Attributes of the spatial feature

        *Parameters*:
            - coords: Coordinates of the point
            - dual: Flag to enable (True) or disable (False) the dual property
        """

        # Convert the coords in float equivalent to c_double which is the internal type of the Geos library
        tuple_coords = self._coords_to_tuple(coords)
        if dual:
            Point.__init__(self, tuple_coords)
            self._coords_dual = tuple_coords
        else:
            Point.__init__(self, tuple_coords)

        self.feature_type = GenUtil.POINT
        self.ma_properties = {} # Dictionary for attributes

    def _coords_to_tuple (self, coords):
        """Convert a coords which can be either an iterable (x,y) or ((x,y)) into a tuple (x,y) containing floats"""

        if (len(coords) == 2):
            tuple_coords = ((float(coords[0]),float(coords[1])),) # Creation of a tuple
        else:
            tuple_coords = ((float(coords[0][0]),float(coords[0][1])),)

        return tuple_coords

    @property
    def coords_dual(self):
        # Manage the coords_dual property
        if self.is_dual():
            return self._coords_dual
        else:
            raise GenException ("Dual property no set... cannot access it")

    @coords_dual.setter
    def coords_dual (self, coords_dual):
        raise InternalError ("Cannot set coordinates using .coords. Use .coords_dual instead...")

    @property
    def coords(self):
        return Point.coords.__get__(self, Point)

    @coords.setter
    def coords (self, coords):
        raise InternalError ("Cannot set coordinates using .coords. Use .coords_dual instead...")

    def is_dual(self):
        """Return True if the feature is dual; false otherwise

        *Parameters*: *None*

        *Returns*: True: The feature is dual
                   False: The feature is not dual
        """

        if (hasattr(self, "_coords_dual")):
            return True
        else:
            return False

    def mid_point(self, point):
        """Return a point in the middle of the 2 points"""

        x = (point._coords_dual[0][0] + self._coords_dual[0][0]) / 2.
        y = (point._coords_dual[0][1] + self._coords_dual[0][1]) / 2.

        return MA_Point([x,y])

    def cloner (self):
        """This method clones the object

        *Parameters*:
            *None*

        *Returns*:
            - A shapely cloned object

        """

        # Clone the geometry
        clone = MA_Point(self.coords_dual)

        # Clone the attributes
        clone.ma_properties = self.ma_properties.copy()

        return clone

    def update_coords (self, coords, s_container=None):
        """Update the coordinates of the feature and update the bounding in the spatial container

        *Parameters*:
            - coords: List of the new coordinates to apply to the feature
            - s_container: SpatialContainer which contains the feature to update. If None the feature
                           is not contained in a spatial container.

        *Returns*:
            *None*

        """

        # Update the coordinates
        if (self.is_dual()):
            # Manage the coords dual property
            # Convert the coords in float equivalent to c_double which is the internal type of the Geos library
            tuple_coords = self._coords_to_tuple(coords)
            Point.coords.__set__(self,tuple_coords)
            self._coords_dual = tuple_coords
        else:
            Point.coords.__set__(self,coords)

        # If the feature is placed in a spatial container; once must update the spatial index
        if s_container is not None:
            s_container.update_spatial_index(self)
        else:
            if s_container is None and hasattr(self, "_sci_id"):
                raise InternalError("Feature is in a spatial container and no spatial container is supplied...")


class MA_LineString(LineString):

    """
    The MA_LineString class extends the LineString of the shapely library by adding
    new methods and attributes that enables the creation of Meta algorithm.
    Here are some possibilities of this class:
     - MA_LineString can be stored  in an instance of a SpatialContainer class
     - MA_LineString can have than coordinates updated through the update_coords
       method which also updates the spatial index if used
     - MA_LineString disable the capacity to update a coordinate through the property coords;
       you must always use the update_coords method to update coordinates
     - MA_LineString can have a copy of the coordinates in an the attribute coords_dual;
       this copy of the attributes allows a faster access to the coordinates than
       when accessing them through Shapely coords attributes. This can be an advantage when
       you have to loop over and over each coordinates of your features.

    """

    def __init__(self, coords, dual=True):
        """Initialize the properties of a MA_LineString object.  The properties of an
        instance are:
          -feature_type: The type of feature
          - _coords_dual: duplicate of the coordinates
          - ma_properties: Attributes of the spatial feature

        *Parameters*:
            coords: List or tuple of coordinates
            dual: Enable (true) or diable (False) the possibility of a dual property

        *Return*: None

        """

        # Convert the coords in float equivalent to c_double which is the internal type of the Geos library
        tuple_coords = tuple( (float(x),float(y)) for x,y in coords)

        if dual:
            self._coords_dual = tuple_coords
            LineString.__init__(self, tuple_coords)
        else:
            LineString.__init__(self, tuple_coords)

        self.feature_type = GenUtil.LINE_STRING
        self.ma_properties = {} # Dictionary for attributes

    @property
    def coords_dual(self):
        if (self.is_dual()):
            return self._coords_dual
        else:
            raise GenException ("Dual property no set... cannot acces it")

    @coords_dual.setter
    def coords_dual (self, lst_coords_dual):
        raise InternalError ("Cannot set coordinates using .coords. Use update_coords method instead...")

    @property
    def coords(self):
        return LineString.coords.__get__(self, LineString)
    @coords.setter
    def coords (self, coords):
        raise InternalError ("Cannot set coordinates using .coords. Use update_coords method instead...")

    def is_dual(self):
        """Return True if the feature is dual; false otherwise

        *Parameters*: *None*

        *Returns*: True: The feature is dual
                   False: The feature is not dual
        """

        if (hasattr(self, "_coords_dual")):
            return True
        else:
            return False

    def cloner (self):
        """
        This method clones the object

        *Parameters*:
            *None*

        *Returns*:
            - A shapely cloned object

        """

        # Clone the geometry
        clone = MA_LineString(self.coords_dual)

        # Clone the attributes
        clone.ma_properties = self.ma_properties.copy()

        return clone

    def update_coords (self, coords, s_container=None):
        """Update the coordinates of the feature and update the bounding in the spatial container

        *Parameters*:
            - coords: List of the new coordinates to apply to the feature
            - s_container: SpatialContainer which contains the feature to update. If None the feature
                           is not contained in a spatial container.

        *Returns*:
            *None*

        """

        # Update the coordinates
        if (self.is_dual()):
            # Convert the coords in float equivalent to c_double which is the internal type of the Geos library
            tuple_coords = tuple( (float(x),float(y)) for x,y in coords)
            self._coords_dual = tuple_coords
            LineString.coords.__set__(self,tuple_coords)
        else:
            LineString.coords.__set__(self,coords)

        # If the feature is placed in a spatial container; once must update the spatial index
        if s_container is not None:
            s_container.update_spatial_index(self)
        else:
            # The feature is in a spatial container and no spatial container is supplied...ERROR
            if s_container is None and hasattr(self, "_sci_id"):
                raise InternalError("Feature is in a spatial container and no spatial container is supplied...")

class MA_Polygon(Polygon):

    """
    The MA_Polygon class extends the Polygon of the shapely library by adding
    new methods and attributes that enables the creation of Meta algorithm.
    Here are some possibilities of this class:
     - MA_Polygon can be stored  in an instance of a SpatialContainer class
     - MA_Polygon can have than coordinates updated through the update_coords
       method which also updates the spatial index if used
     - MA_Polygon disable the capacity to update a coordinate through the property coords;
       you must always use the update_coords method to update coordinates

    """

    def __init__(self, exterior, interiors = None):
        """Initialize the properties of an area object"""

        Polygon.__init__(self, exterior, interiors)
        self.feature_type = GenUtil.POLYGON
        self.ma_properties = {} # Dictionary for attributes

    def cloner (self):
        """
        This method clones the object

        *Parameters*:
            *None*

        *Returns*:
            - A shapely cloned object

        """

        # Clone the geometry
        clone = MA_Polygon(self.exterior, self.interiors)

        # Clone the attributes
        clone.ma_properties = self.ma_properties.copy()

        return clone

    def update_coords (self, coords):
        """Update the coordinates of the feature and update the bounding in the spatial container

        *Parameters*:
            - coords: List of the new coordinates to apply to the feature
            - container: Spatial container into which we update the spatial index

        *Returns*:
            *None*

        """

        # Update the coordinates
        raise InternalError ("Method not yet implemented...")

class ErrorPosition(object):
    """ Abstract class used to manage the geometry (spatial position) of the errors

    Errors can be of type Point, LineString or Polygon.

    *Parameters*:
        - id: The id of the error feature.  Usually it's the ID of the input feature which made the error
        - error_type: Name of the error

    *Returns*:
        - object of class ErrorPosition
        
    """
    
    _err_id = -1

    def __init__(self, id, error_type ):
        """Creator of the abstract class ErrorPosition"""

        self.ma_properties = {}
        if id is None:
            ErrorPosition._err_id += 1
            id = ErrorPosition._err_id
        self.ma_properties['id'] = id
        self.error_type = error_type

class PointErrorPosition(ErrorPosition, Point):
    """Maintains the position of a point error

    *Parameters*:
        - id: The id of the error feature.  Usually it's the ID of the input feature which made the error
        - error_type: Name of the error
        - coord: (x,y) tuple position of the error

    *Returns*:
        - object of class PointErrorPosition

    """


    def __init__(self, id, error_type, coord):
        """Creator of the class PointErrorPosition"""

        ErrorPosition.__init__(self, id, error_type)
        Point.__init__(self, coord)

class LineStringErrorPosition(ErrorPosition, LineString):
    """Maintains the position of a line error

    *Parameters*:
        - id: The id of the error feature.  Usually it's the ID of the input feature which made the error
        - error_type: Name of the error
        - coords: list of (x,y) tuple position of the error

    *Returns*:
        - object of class LineStringErrorPosition

    """

    def __init__(self, id, error_type, coords):
        """Creator of the class LineStringErrorPosition"""

        ErrorPosition.__init__(self, id, error_type)
        LineString.__init__(self, coords)

class PolygonErrorPosition(ErrorPosition, Polygon):
    """Maintains the position of a polygon error

    *Parameters*:
        - id: The id of the error feature.  Usually it's the ID of the input feature which made the error
        - error_type: Name of the error
        - exterior: list of (x,y) tuple position of the exterior of the polygon
        - interiors: contains 0 to n lists of (x,y) tuple defining the interiors of the polygon

    *Returns*:
        - object of class PolygonErrorPosition

    """

    def __init__(self, id, error_type, exterior, interiors=None):
        """Creator of the class PolygonErrorPosition"""

        ErrorPosition.__init__(self, id, error_type)
        Polygon.__init__(exterior, interiors)

class IterationResults(object):
    """Contains the iteration results when iterative algorithm are involved

    *Parameters*:
        - iter: list containing the intermediate result for each iteration

    *Returns*:
        - object of class IterationResults
            
    """

    def __init__(self):
        """Creator of the class IntermediateResults"""

        self.iters = []

    def add_iteration(self):
        """Add a new iteration in order to add new feature"

        *Parameters*: *None*

        *Returns*: *None*

        """

        self.iters.append([])

    def add_features(self, features):
        """Add intermediate results for the current iteration

        *Parameters*:
            - features: list of features to add

        *Returns*: *None*

        """

        # Add the points
        for feature in features:
            self.iters[-1].append(feature.cloner())

class GenStatistics (object):
    """Abstract class that contains the statistics for a process
        
        *Parameters*: *None*

        *Returns*:
            - object of class GenStatistics

    """

    DETAILED = "Detailed"
    SUMMARY  = "Summary"

    def __init__(self):
        """Creator of the class Genstatistics"""
        self.stats_iters = []

    def print_stats (self, type=DETAILED):
        """Print the statistics with a simple print command"""

        strs = self.get_stats(type)

        for str in strs:
            print (str)

    def add_iteration(self):
        """Add one iteration to the statistics list

        *Parameters*: *None*

        *Returns*: *None*

        """

        stats_iter = self.stats_iters.append({})
        for stats_name in self.stats_names:
            self.stats_iters[-1][stats_name] = 0

    def add_stats(self, stats_name, value=1):
        """Add one to the statistics count for the specified statistics name

        *Parameters*:
            - stats_name: the name of the statistics to add
            - value: Number of increment  to add

        *Returns*: *None*

        """

        self._validate_parameters(stats_name)
        self.stats_iters[-1][stats_name] += value
        
    def reset_stats_names(self, lst_stats_names, reset_value=0):
        """Reset the value of the statistics names for the specified statistics value
        
        *Parameters*:
            - lst_stats_names: List of the statistics names to reset
            - reset_value: Value to reset
            
        *Returns*: *None*
        
        """
        
        for stats_name in lst_stats_names:
            self._validate_parameters(stats_name)
            self.stats_iters[-1][stats_name] = reset_value
        

    def get_stats_name_count_iter (self, stats_name, iter=None):
        """Get the count of statistics for a specific statistic name and iteration

        *Parameters*:
            - stats_name: Name of the statistics
            - iter: Iteration number to get the statistics.  If no value is given it will take the stats from the last iteration

        *Returns*:
            - Integer: statistic count for the specified statistics name

        """

        if iter == -1:
            iter = None

        self._validate_parameters(stats_name, iter)

        if iter is None:
            count =  self.stats_iters[-1][stats_name]
        else:
            count =  self.stats_iters[iter][stats_name]

        return count

    def get_stats_name_count_total (self, stats_name):
        """Get the total count of statistics for a specific statistic name and all iteration

        *Parameters*:
            - stats_name: Name of the statistics

        *Returns*:
            - Integer: statistic count for the specified statistics name

        """

        self._validate_parameters(stats_name)
        total = 0
        for stats_iter in self.stats_iters:
            total +=  stats_iter[stats_name]

        return total

    def get_nbr_iteration (self):
        """Get the total number of iteration

        *Parameters*: *None*

        *Returns*:
            - Integer: total number of iteration

        """

        return (len(self.stats_iters))

    def _validate_parameters(self, stats_name, iter=None):
        """Validate the name of the statistics and the iteration number

        *Parameters*:
            - stats_name: Name of the statistics
            - iter: iteration number

        *Returns*:
            - Bolean: True: the parameters are valid
            - Exception: the parameters are invalid

        """

        valid = True
        if (stats_name not in self.stats_names):
            raise InternalError('Unknown statistics name: ' + stats_name)

        if (iter is not None):
            if (str(iter).isdigit() and iter <= len(self.stats_iters)-1):
                pass
            else:
                raise InvalidIteration ('Invalid iteration number: ' + iter)

        return valid

class Algorithm(object):
    """This is the base class for the input / output of the genmetal algorithms

    Attributes:
        - points: input/output points
        - line_strings: input/output line strings
        - polygons input/output polygons
        - error_position: output spatial position of the errors
        - iter_results: output spatial position of the iterative intermediate results
        
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize the attributes of an object of the class"""

        self.features = []
        self.error_positions = []
        self.iter_results = IterationResults()

    def check_integrity (self, features):
            """
            Checks the integrity of the data structure by checking that the coords
            contained in line.coords and line.lst_coords are the same.

            This check is dome in development mode to verify that a programming error
            have corrupted the data structure.

            """

            for feature in features:
                if feature.feature_type == GenUtil.POINT:
                    if ( not feature.almost_equals(Point(feature.coords_dual),6) ):
                        raise InternalError ("Internal corruption of MA_Point coords structure")

                if feature.feature_type == GenUtil.LINE_STRING:
                    if ( not feature.almost_equals(LineString(feature.coords_dual),6) ):
                        raise InternalError ("Internal corruption of MA_LineString coords structure")

            return

    def load_features(self):
        """Load the points, line strings and polygons in the spatial container

        *Parameters*: *None*

        *Returns*: *None*
        
        """

        # Create the spatial container that will receive all the spatial features
        s_container = SpatialContainer()
        
        # Load all the features in the spatial container
        s_container.add_features(self.features)
        
        # Reset the feature container
        self.features = []

        return (s_container)

    def extract_features_out (self, lst_features):
        """Extract the features from the container and fill the output container

        *Parameters*:
            - *TBD*

        *Returns*: *None*

        """

        for feature in lst_features:
            self.features.append(feature)
                
        return 

    @abstractmethod
    def process(self):
        """Start the process of a generalization algorithm

        This method must be overriden

        """

class SpatialContainer(object):
    """This class manages the spatial features and a spatial index.

    This class enables the management of spatial features by incorporation 
    transparently a spatial index.  The spatial index is an implementation
    of the Rtree open source softawre.  The spatial container offers the following
    main options:
      - add features in the constainer and update the spatial index
      - delete features in the container and update the spatial index
      - update the coordinates of a feature and update the spatial index if needed
      - make spatial queries by bounding box
      - make attributes queries
      - delete the container

    """

    # Class variable that contains the Spatial Container Internal ID
    _sb_sc_id = 0
    
    def __init__(self, line_opt_value=0):
        """Create an object of type SpatialContainer

        The init will create one container for the feature a dictionary and one
        container for the spatial index (Rtree)

        *Parameters*: 
            - line_opt_value: Value of 0 or greater than 2.  The line optimizer value helps
                              to optimize the search in the Rtree.  When very long LineString are
                              used the bounding becomes very large and when a bounding box intersection
                              is done all the long LineString are intersected.  When this parameter is used
                              (value != 0), sub bounding boxes are calculated like if we were splitting
                              the line in smaller lines where each sublines is conraining a maximum of
                              line_opt_value coordinates.  So when this value is small more bounding boxes
                              are created and the Rtree searches are finer and when the value is high
                              less bounding boxes are created and the Rtree searches are coarser. This 
                              parameters applies to LineString features only.  

        *Returns*: *None*

        """

        self._r_tree = Rtree()  # Container for the Rtree
        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def _is_bbox_the_same(self, feature, old_lst_bbox, new_lst_bbox):
        """Checks if the bbox are the same
        
        This method checks if the bounding needs to be updated. When there is only
        one bbox in the old and new lst_bbox, we compare them and if they are the
        same, the bbox are the same.  If in the new_lst_bbox there is more than 
        one bbox than we check that all the coordinates of the feature are contained 
        in the old_lst_bbox; if they are all in the old_lst_bbox the bbox are the same 
        otherwise the bbox are different
        
        *Parameters:*
            - feature: MA_* spatial feature to check
            - old_lst_bbox: List of the bbox corresponding to the previous feature
            - new_lst_bbox: List of the bbox corresponding to the updated feature
        """ 
        
        if (len(old_lst_bbox) == 1 == len(new_lst_bbox)):
            # We just check that the bbox are the same
            old_bbox = old_lst_bbox[0]
            new_bbox = new_lst_bbox[0]
            if ( new_bbox[0] != old_bbox[0] or
                 new_bbox[1] != old_bbox[1] or
                 new_bbox[2] != old_bbox[2] or
                 new_bbox[3] != old_bbox[3] ):
                is_the_same = False
            else:
                is_the_same = True
        else:
            if (len(new_lst_bbox)==1):
                # Because the new list of box contains only one box we consider it is always better
                # to reduce the number of bbox so we consider them not the same
                is_the_same = False
            else:
                # This the more complex case wehre we loop over each coordinates in order to check if allt he
                # coordinate are contained in the old bbox if so we don't need to recreate the spatial index
                if (feature.is_dual()):
                    line_coords = feature.coords_dual    
                else:
                    line_coords = list(feature.coords)
                i_bbox = 0
                len_old_bbox = len(old_lst_bbox)
                xmin,ymin = old_lst_bbox[i_bbox][0], old_lst_bbox[i_bbox][1]
                xmax,ymax = old_lst_bbox[i_bbox][2], old_lst_bbox[i_bbox][3]
                try:
                    for coord in line_coords:
                        if ( not (xmin <= coord[0] <= xmax and ymin <= coord[1] <= ymax)  ):
                            # Try to find the next bbox that contains the coordinate
                            i_bbox += 1
                            if (i_bbox < len_old_bbox):
                                xmin,ymin = old_lst_bbox[i_bbox][0], old_lst_bbox[i_bbox][1]
                                xmax,ymax = old_lst_bbox[i_bbox][2], old_lst_bbox[i_bbox][3]
                                if ( not (xmin <= coord[0] <= xmax and ymin <= coord[1] <= ymax) ):
                                    # The coordinate is outside all bounding boxes
                                    raise Exception
                    is_the_same = True
                except Exception:
                    is_the_same = False
                except:
                    raise ("Unknown error...")  
                
            
        return is_the_same
    
    def _adjust_bounding_box (self, bounds):
        """Modify the bounds of a feature when the bounds of almost zero
        
        *Parameters*: 
            - bounds: Tuple of a bounding box (xmin, ymin, xmax, ymax)
            
        *Returns*: 
            - Tuple of a bounding box (xmin, ymin, xmax, ymax)
        
        """
        
        if (abs(bounds[2]-bounds[0]) < GenUtil.ZERO or abs(bounds[3]-bounds[1]) < GenUtil.ZERO):
            if abs(bounds[2]-bounds[0]) >= GenUtil.ZERO:
                xmin = bounds[0]
                xmax = bounds[2]
            else:
                xmin = bounds[0] - GenUtil.ZERO
                xmax = bounds[2] + GenUtil.ZERO
        if abs(bounds[3]-bounds[1]) >= GenUtil.ZERO:
            ymin = bounds[1]
            ymax = bounds[3]
        else:
            ymin = bounds[1] - GenUtil.ZERO
            ymax = bounds[3] + GenUtil.ZERO
        bounds = (xmin, ymin, xmax, ymax)
            
        return bounds
    
    def _extract_bounding_box(self, feature):
        """Extract the bounding box of a features
        
        If the feature is a LineString and the line optimizer is on; the line is splitted in 
        smaller fragment and a bounding box is computed for each small fragment.
        
        *Parameters*: 
            - feature: Shapely feature (Point, LineString or Polygon
            
        *Returns*: 
            - List of tuple of bounding boxes in the form of [(xmin,ymin,xmax,ymax),(...),...]
              It will always contains at least one bounding box
        
        """
        
        bounds = feature.bounds
        # if (isinstance(feature, MA_Point) or isinstance(feature, MA_Polygon)):
        #     lst_bounds = [feature.bounds]
        # else:
        #     if (self._line_opt_value == 0):
        #         # Line optimizer is diabled
        #         lst_bounds = [feature.bounds]
        #     else:
        #         if (feature.is_dual()):
        #             line_coords = feature.coords_dual
        #         else:
        #             line_coords = list(feature.coords)
        #
        #         max_coords = self._line_opt_value - 1
        #         # Split the list of coordinates into a list of list of coordinates where each list of coordinates contains a
        #         # maximum of max_coords. It also copy the last coordinate of one group as the first coordinate of the next group
        #         split_coords = [line_coords[i:i+max_coords+1] for i in range(0, len(line_coords), max_coords)]
        #         # If the last group contains only 1 coordinate (x,y) we delete it
        #         if (len(split_coords[-1]) == 1):
        #             del split_coords[-1]
        #
        #         # Calculates for each list of coordinates its bounding box in the form (xmin,ymin,xmax,ymax)
        #         lst_min_max   =  [ map((lambda lst: (min(lst),max(lst)) ),zip(*coords)) for coords in split_coords]
        #         lst_bounds    =  [ (min_max[0][0],min_max[1][0], min_max[0][1], min_max[1][1]) for min_max in lst_min_max ]
        #
        #if len(lst_bounds) == 1:
            #Presently there is a little bug in RTree when the xmin and xmax or ymin and ymax are the same value
            #The problem is resolved when we add a very small delta between the 2 values
            #This bug is supposed to be solved in the next release. We're running now on 0.6
        #    lst_bounds[0] = self._adjust_bounding_box (lst_bounds[0])
            
        return bounds

    # def _set_if_statement(self, filter, remove_keys):
    #     """Set the if statement needed to extract features from the container
    #
    #     *Parameters*:
    #         - filter: This parameter define the if statement to be used in order to filter
    #                 the features based on the value of some properties. If *None* no filter are applied.
    #                 Example of a filters to filter on feature type: "if feature.feature_type == 'LINE'"
    #                 The filters must be a valid python expression
    #         - remove_keys: List of keys to be removes from the selection
    #
    #     *Returns*:
    #         - String containing the if statement
    #
    #     """
    #
    #     if remove_keys is not None and len(remove_keys) != 0:
    #         str_remove_keys = "".join([str(key)+"," for key in remove_keys])
    #         str_remove_keys = "feature._sci_id not in [%s]" %(str_remove_keys)
    #     else:
    #         str_remove_keys = ""
    #
    #     str_filter = filter
    #     if (str_filter is not None):
    #         str_filter = str_filter.rstrip()
    #         str_filter = str_filter.lstrip()
    #         str_filter = "(" + str_filter + ")"
    #     else:
    #         str_filter = ""
    #
    #     str_if = ""
    #     if (str_filter != ""):
    #         str_if = str_filter
    #
    #     if (str_remove_keys != ""):
    #         if (str_filter != ""):
    #             str_if = "%s and %s" %(str_filter, str_remove_keys)
    #         else:
    #             str_if = str_remove_keys
    #
    #     return str_if
    
    def add_feature(self, feature):
        """Adds a feature in the container and update the spatial index with the feature's bound
        
        To be added in the container a spatial feature must be a MA_Point, MA_LineString or
        MA_Polygon. 

        *Parameters*:
            - feature: A spatial feature derives from MA_Point, MA_LineString or MA_Polygon

        *Returns*: *None*

        """
        
        # Check if the type is valid
        if feature.sb_geom_type == "Point" or feature.sb_geom_type == "LineString":
            pass
        else:
            raise GenException ('Unsupported feature type...')
        
        # Check if the feature is already in a spatial container
        if hasattr(feature, "_gbt_sc_id"):
            raise GenException ('Feature is already in a spatial container')
        
        bounds = self._extract_bounding_box(feature)

        # Container unique internal counter
        SpatialContainer._sb_sc_id += 1
        
        # Add the spatial id to the feature
        feature._sb_sc_id = SpatialContainer._sb_sc_id
        
        # Add the feature in the feature container 
        self._features[feature._sb_sc_id] = feature

        # Add the bounding box in the bbox_container
        self._bbox_features[feature._sb_sc_id] = bounds
        self._r_tree.add(feature._sb_sc_id, bounds)
        
 #       # Add the each boundfing box in the RTree
 #       for bounds in lst_bounds:
 #           self._r_tree.add(self._sci_id, bounds)

        return

    def add_features(self, features):
        """Adds a list of feature in the continer and update the spatial index with the feature's bound

        *Parameters*:
            - feature: A spatial feature derives from Point, LineString or Polygon

        *Returns*: *None*

        """

        for feature in features:
            self.add_feature(feature)

        return

    def del_feature(self, feature):
        """Delete the feature in the spatial container and in the RTree.
        
        If the feature is included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...

        *Parameters*:
            - feature: The feature to delete in the spatial container

        *Returns*:
            - 0 if feature is deleted from the patial container
            - 1 if feature was not included in the spatial container

        """
        
        ret_value = 0
        
        # Check if the feature has a container_key
        if hasattr(feature, "_gbt_sci_id"):
            
            if (feature._gbt_sci_id in self._features and
                feature._gbt_sci_id in self._bbox_features):

                try:
                    # Retreive the bounding boxes of this feature
                    lst_bbox = self._bbox_features[feature._sci_id]
                    # Delete the feature from the features and the bbox_features
                    del self._features[feature._sci_id]
                    del self._bbox_features[feature._sci_id]
                    # Delete the different bounds in RTree
                    self._r_tree.delete(feature._sci_id, bbox)
                    #Delete the property _sci_id
                    del feature._sci_id
                except:
                    raise InternalError ("Internal corruption, problem with the container and/or the RTree")
            else:
                raise InternalError ("Internal corruption, key {} has disappear...".format(feature._sci_id))
                
        return ret_value


    def del_features(self, features):
        """Delete a list of features in the spatial container
        
        If the features are included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...

        *Parameters*:
            - features: list of features to delete

        *Returns*:
            - List of value for each feature to delete.
            - 0 if feature is deleted from the patial container
            - 1 if feature was not included in the spatial container  
           
        Exception InternalError: If the key is not in one of the structure

        """

        ret_value = []
        
        for feature in features:
            ret_value.append(self.del_feature(feature))
            
        return ret_value
        
    def update_spatial_index (self, feature):
        """Update the bounds of the feature in the spatial index
        
        It will only modify the Rtree spatial index if the bounding
        box of the feature is changed in comparison with the old one.

        *Parameters*:
            - feature: Feature containing the bounds to update
            
        *Returns*: *None*
        
        """
                    
        old_bbox = self._bbox_features[feature._sci_id]
        new_bbox = self._extract_bounding_box(feature)
        
        if (self._is_bbox_the_same(feature, old_bbox, new_bbox)):
            # Nothing special to do
            pass
        else:
            # The bounding box has changed
            # Delete The old bounding box in Rtree
            self._r_tree.delete(feature._gbt_sci_id, old_bbox
                                   )
            #Add the new bounding boxes in Rtree
            self._r_tree.add(feature.__gbt_sci_id, new_bbox)
            
            #Save the bounding boxes
            self._bbox_features[feature._sci_id] = new_bbox

        return

    def get_keys_by_bounds(self, bounds, keys_to_remove=None):
        """Extract keys in the container based on the value of a bounding box

        *Parameters*:
            - bounds: Bounding box defined as a list: xmin, ymin, xmax, ymax
            - keys_to_remove: List of keys to remove  removed from the list of key features returned

        *Returns*:
            - List of keys contained in the bounding box

        """

        # Extract the keys from the RTree
        keys = list(self._r_tree.intersection(bounds))
        
        # Remove the keys which are in the keys to remove list
        keys = list(set(keys)- set(keys_to_remove))

        return keys
    
    def get_features(self, bounds=None, filter=True, remove_features=[]):
        """Extract the features from the spatial container.
        
        According to the parameters the extraction can manage the extraction based on a bounding box using 
        the spatial index RTree, some filters to manage extraction based on properties and the possibility
        to remove specific features based on a list of keys

        *Parameters*:
            - bounds: Bounding for the spatial extraction. *None* means all the features
            - filter: This parameter define the if statement to be used in order to filter
                      the features based on the value of some properties. If *None* no filter are applied.
                      The filters must be a valid python expression.
                      Example of a filter to filter on feature type: "if feature.feature_type == 'LINE'"
            - remove_keys: List of keys to be removed from the selection 

        *Returns*:
            - List of features extracted from spatial container

        """

        # Extract the features by bounds if requested
        if (bounds != None):
            # Extract features by bounds
            keys = self.get_keys_by_bounds(bounds, remove_features)
            features = (self._features[key] for key in keys if key in self._features)
        else:
            features = (feature for feature in self._features.values() if feature not in remove_features)

        # Filter the result if requested
        if filter:
            features = [feature for feature in features if filter]

        return features


class ChordalAxisTransformer(object):
    """This class is creating  a skeleton and identify bottleneck based on the Chordal Axis Transform CAT
    
       The CAT is very interesting as it simulate the skeleton created by the Medial Axis Transform (MAT).
        
       The CAT is created using the output of the Constrained Delanauy Triangulation (CAT).
       The class containes only 2 public method
          - get_skeleton: To extract the skeleton (line network)
          - get_triangle: To extract the original triangles. Each triangle containing a the property "type"
                          which can take 2 values:
                              - bottleneck: If the width of the triangle is below the minimal width
                              - other:  If the width of the trioangle is over the minimal width
                              
        This class also use the class _Triangle and _LineSegments
    
       For more information of the Chordal Axis Transform refer to:
          Rectification of the Chordal Axis Transform and a new criterion for shape decompositio,  Laksham Prasad,
          12 International Conference on Discrte Geometry for Computer Imagery, Poitiers, France, 13-15 April 2005
    """
    CODE = 'code'
    WIDTH = 'width'
    CENTER_LINE = 'center_line'
    BOTTLENECK = 'bottleneck'
    OTHER = 'other'
    
    def __init__(self, polygon, triangles, minimal_width, search_tolerance=GenUtil.ZERO):
        """Constructor of the class
        
        Parameters:
            - polygon: MA_Polygon to process
            - triangles: List of MA_LineString. The line represent the triangles of a Constriant Delanuay Triangulation (CDT)
                         Each triangle is composed of 4 non colinear vertice
            - minimal_width: Float used to prune the skeleton outputted by the Chordal Axis Transform and to identify 
                             bottleneck triangles
            - search_tolerance: Search tolerance can vary depending on the dynamic of the data set from lat-lon to Lambert conformal
                             
        Return value: None
        """
        self.s_cont_triangles = SpatialContainer()
        self._minimal_width = minimal_width
        self._search_tolerance = search_tolerance
        self._process_polygon(polygon)
        self._load_triangles(triangles)
        
        _Triangle.line_segments = self.line_segments
        _Triangle.perimeters = self.perimeter_distance
        
        self._build_skeleton()
        self._prune_skeleton()
          
        
    def _process_polygon(self, polygon):
        """Process a polygon to create the object property line_segments and perimeter_distance
        
        Parameters: 
          - polygon: MA_Polygon to process
          
        Return value: None 
        """
        
        # Create one list containing the exterior ring plus the interior rings (if any)
        rings = [polygon.exterior]
        for ring in polygon.interiors:
            rings.append(ring)
        
        # _line_Segments object are used to calculate the distance from a point to the different line segment of the polygon
        self.line_segments = _LineSegments(rings, self._search_tolerance)
        
        # Perimeter distance objects are used to calculate the shortest distance between 2 vertice on a ring
        self.perimeter_distance = PerimeterDistance(rings, self._search_tolerance)
        
    def _load_triangles(self, triangles):
        """Load the triangles 
        
        Parameters: 
            - List of MA_LineString.  Each MA_LineString is composed of 4 non colinear vertice forming a triangle
            
        Return value: None
        """
        
        for tri in triangles:
            triangle = _Triangle(list(tri.coords))
            self.s_cont_triangles.add_feature(triangle)
            
    def _build_skeleton(self):
        """Build the Chordal Axis Transform skeleton with the triangles
        
        Parameters: None
        
        Ruturn value: None
        
        """
        
        for triangle in self.s_cont_triangles.get_features():
            triangle.get_centre_line()
            
    def _prune_skeleton(self):
        """Prune the noise of Chordal Axis Transform skeleton to remove detail below the minimum width
        
        Remove the small lines of the skeleton.  In the chordal axis tranform each junction triangles
        is a bifurcation (skeleton splitting).  The is removing noisy skeleton arm that are 
        below a treshold.  We use the distance along the perimeter to determine if an arm must be pruned
        
        Refer to the original paper for more details
        
        Parameters: None
        
        Return value: None 
         
        """
        
        # Process all the triangles
        for triangle in self.s_cont_triangles.get_features():
            if triangle.get_nbr_internal() == 3:
                # Only process Junction triangle
                side_demoted = []  # Contains the side number of the triangle that will be pruned
                for i in xrange(3):
                    p0 = triangle.coords_dual[i]
                    p1 = triangle.coords_dual[(i+1)%3]
                    # Check if the perimeter distance is below the minimal width
                    extremity = self.perimeter_distance.is_extremity(p0,p1, self._minimal_width)
                    if extremity:
                        # The side/skeleton must be pruned
                        lst_sub_coords = self.perimeter_distance.get_sub_perimeter(p0,p1)
                        # Creation of a polygon with the vertice included in the perimeter
                        extremity_polygon= Polygon(lst_sub_coords)
                        centroid_coord = triangle.get_centroid()  # Centroid of the triangle in hand
                        centroid_point = Point(centroid_coord)
                        # Make sure that the perimeter polygon do not contain the triangle in hand
                        if (extremity_polygon.disjoint(centroid_point)):
                            side_demoted.append(i)
                            
                            # All the triangle located inside the polygon must have there skeleton removed
                            # Proceed in 2 phases
                            # First phase: Bounding box search
                            potential_triangles = self.s_cont_triangles.get_features(bounds=extremity_polygon.bounds)
                            for potential_triangle in potential_triangles:
                                centroid_coord = potential_triangle.get_centroid()
                                centroid_point = Point(centroid_coord)
                                # Second phase spatial search
                                if (extremity_polygon.contains(centroid_point)):
                                    # Triangle pruned from its triangle
                                    potential_triangle.del_centre_line()
                                    
                # If some side of the tirangle were demoted once have to recalculate the skeleton for this triangle
                nbr_side_demoted = len(side_demoted) 
                if nbr_side_demoted >= 2:
                    # If 2 or 3 side were demoted we also remove the skeleton from this triangle
                    triangle.del_centre_line()
                else:
                    if nbr_side_demoted == 1:
                        # The center_line of the triangle has changed
                        triangle.demote_junction(side_demoted[0])
                    else:
                        # Nothing has changed on the state of the triangle
                        pass
                                    
    def get_skeletton(self):
        """Extract the Chordal Axis Transform skeleton from a constrained Delanauy trianulation
        
        Parameters: None
        
        Return value: 
            - List of MA_LineString of the skeleton of the polygon
        """
        
        center_lines = []
        for triangle in self.s_cont_triangles.get_features():
            center_lines += triangle.get_center_line()
            
        return center_lines
    
    def get_triangles(self):
        """Extract the triangles from the constrained delanauy triangulation
        
        Parameters: None
        
        Return value:
            - List of the MA_LineString of 4 vertice each.
              Each triangle also have the following ma_properties
                  - type: The type of triangle can take 2 values:
                          - bottleneck: If the width of the triangle is below the minimal widt
                          - other: If the width of the triangle is over the minimal width
                  - width: The width of the triangle.  If the width is below the minimal width the
                           exact value of width is return; otherwise the value 1.0E+99 is output. The real width
                           of each triangle is not calculated for a performance reason. 
            
        """
        triangles = []
        for triangle in self.s_cont_triangles.get_features():
            (category, width) = triangle.get_category(self._minimal_width)
            tri = triangle.cloner() # Create a copy of the triangle not a reference
            tri.ma_properties[ChordalAxisTransformer.CODE] = category 
            tri.ma_properties[ChordalAxisTransformer.WIDTH] = width
            tri.ma_properties[ChordalAxisTransformer.CENTER_LINE] = triangle.get_centre_line()
            triangles.append(tri)
            
        return triangles

class _Triangle(MA_LineString):
    """Calculates the 
    """
    
    SUPERIMPOSED = 'superimposed'
    INTERNAL = 'internal'
    
    
    line_segments = None
    perimeters = None

    
    def __init__(self, lst_coords ):
        if (len(lst_coords) == 4):
            MA_LineString.__init__(self, lst_coords)
            self._nbr_internal = None  # 
            self._centre_lines = None
            self._side_type = None
            self._category = None
            self._mid_triangle = None
        else:
            raise Exception ("A triangle must have 4 and only 4 coordinates")
            
    def _get_obtuse_angle(self):
        """Return the vertice number (0, 1, 2) of the obtuse angle of the triangle
        
        Parameters: None
        
        Return value:
            - integer: Vertice number of the obtuse angle or None if there is not obtuse angle
            
        """
        
        acute_angle = None
        
        for i in xrange(3):
            p0 = self.coords_dual[(i-1)%3]
            p1 = self.coords_dual[(i)%3]
            p2 = self.coords_dual[(i+1)%3]
            angle = GenUtil.compute_angle(p0, p1, p2)
            if (angle > 90.):
                # There is only angle greater than 90 in a triangle so break after
                acute_angle = i 
                break
            
        return acute_angle
    
    def _is_acute_triangle(self):
        """Check for sleeve triangle, if the triangle is acute or obtuse
        
        To be considered acute, a sleeve triangle must have the 2 angles of the side superimposed with the polygon
        below 90 degrees.
        
        This method is also calculating the heigth (which is the width of the bottleneck) of the triangle using the
        Heron formula
        
        Parameters: None
        
        Return value:
            boolean: True: The sleeve polygon is acute
                     False: The sleeve polygon is not acute
        
        """
        
        for i, side_type in enumerate(self._side_type): 
            if (side_type == _Triangle.SUPERIMPOSED):
                superimposed = i
                        
        p0_base  = self.coords_dual[superimposed%3]
        p1_base  = self.coords_dual[(superimposed+1)%3] 
        p_summit = self.coords_dual[(superimposed+2)%3]
                    
        # Extract the side of the triangle
        a = GenUtil.distance(p0_base, p1_base)
        b = GenUtil.distance(p1_base, p_summit)
        c = GenUtil.distance(p0_base, p_summit)
        
        # Check if it's side a has an obtuse angle
        angle_p1 = GenUtil.compute_angle(p_summit, p0_base, p1_base)
        angle_p2 = GenUtil.compute_angle(p0_base, p1_base, p_summit)
        
        # Accute triangle
        # Take the height of the triangle as the bottleneck which is a better
        # evaluation of the width than the smallest side
        # Calculates the area of the triangle using the Heron formula
        p = (a+b+c)/2. # Calculates the half perimeter
        # Now use the Heron formula  A = sqrt(s(s-a)(s-b)(s-c)) where s = (a+b+c)/2
        area = (p*(p-a)*(p-b)*(p-c))**0.5
        base = a

        # Calculate the height using area = (base*height)/2 ===> height= (2*area)/base
        self._height = (2.*area) / base
        
        if (angle_p1 < 90. and angle_p2 < 90.):
            acute = True
        else:
            acute = False
            
        return acute


    def get_centroid (self):
        """Calculate the position of the baricenter of the triangle
        
        The gravity centre is always 2/3 the distance between the middle of one side and the opposite angle
        
        Parameters: None
            
        Return value
            Tuple of (x,y) float representing the position of the centroid
            
        """
        
        mid_point = GenUtil.mid_point(self.coords_dual[0], self.coords_dual[1])
        centroid = GenUtil.rescale_vector(mid_point, self.coords_dual[2], 1./3.)
        
        return centroid
    
    def demote_junction(self, side_number):
        """Demote a junction polygon if the skeleton of one of the side is considered as noise and calculate a new skeleton
        
        The new skeleton is calculated by joining to 2 remaining mid side
        
        Parameter:
            side_number: The side number of the triangle to demote. Value [0..2]
        """
        
        p0 = self.coords_dual[(side_number+1)%3]
        p1 = self.coords_dual[(side_number+2)%3]
        p2 = self.coords_dual[(side_number+3)%3]
        
        mid_p0_p1 = GenUtil.mid_point(p0, p1)
        mid_p1_p2 = GenUtil.mid_point(p1, p2)
        
        centre_line = MA_LineString([mid_p0_p1, mid_p1_p2])
        
        self._centre_lines = [centre_line]
        
    def del_centre_line(self):
        """Delete the center line for this triangle
        
        Parameters: None
        
        Return value: None
        
        """
        
        self._centre_lines = []
    
    def get_centre_line(self):
        """Calculates and extract the center line of the triangle
        
        The center line depends of the type of triangle
            Terminal triangle: No center line
            Sleeve triangle: Joining the mid side of the internal side
            Junction triangle:  If the triangle is obtuse: 
                                           - find the opposite side of the obtuse angle ()
                                           - creates 2 lines from the opposite side to the mid point of the 2 other sides
                                If the triangle is acute:
                                        - calculate the baricenter
                                        - creates 3 lines from the mid side to the baricenter
        """
        
        if self._centre_lines is None:
            # Implement late calculation
            self._centre_lines = []    # List of the centre lines
            mid_side_points = [] # List of the mid point on each side of the triangle
            internal_sides = []  # List of the number of the internal side
            external_sides = []  # List of the number of the external side
            
            nbr_internal = self.get_nbr_internal()
            
            for i in xrange(3):
                mid_side_points.append(GenUtil.mid_point(self.coords_dual[i], self.coords_dual[i+1]))
                if self._side_type[i] == _Triangle.INTERNAL:
                    internal_sides.append(i)
                else:
                    external_sides.append(i)
                    
            # Process each case depending on the number of internal side of the triangle
            if (nbr_internal == 0):
                # Degenerated polygon with one triangle no skeleton line added
                pass

                    
            if (nbr_internal == 1):
                # Terminal triangle no skeleton line added
                pass
    
            if (nbr_internal == 2):
                # Sleeve triangle skeleton added between the mid point of each chord
                internal_side0 = internal_sides[0]
                internal_side1 = internal_sides[1]
                self._centre_lines.append( MA_LineString([mid_side_points[internal_side0], mid_side_points[internal_side1]]))
                self._mid_triangle = GenUtil.mid_point(mid_side_points[internal_side0], mid_side_points[internal_side1])
    
            if (nbr_internal == 3):
                # Junction triangle skeleton added.
                obtuse_angle = self._get_obtuse_angle()
                if (obtuse_angle is None):
                    # With an acute triangle a mid point is calculated in the middle of the triangle
                    centroid = self.get_centroid()
                    self._mid_triangle = centroid
                    for mid_side_point in mid_side_points:
                        self._centre_lines.append(MA_LineString([centroid, mid_side_point]))
                else:
                    # With an obtuse triangle the mid point is placed on the sided opposite to the obtuse angle
                    opposite_side = (obtuse_angle+1)%3
                    left_side = (opposite_side+1)%3
                    right_side = (opposite_side-1)%3
                    self._mid_triangle = mid_side_points[opposite_side]
                    self._centre_lines.append(MA_LineString([mid_side_points[opposite_side], mid_side_points[left_side]]))
                    self._centre_lines.append(MA_LineString([mid_side_points[opposite_side], mid_side_points[right_side]]))
        else:
            # Centre line was already calculated... nothing to do
            pass
                
        return self._centre_lines
        
    def get_nbr_internal(self):
        """Extract the number of side of the triangle which are completely inside the polygon.
        
        Three scenarios are possible 0, 1, 2 ou 3 sides completely inside the polygon:
            0: The triangle is completely inside. This is called a Junction triangle
            1: The Triangle as on side that lies on the polygon. This is called a Sleeve triangle
            2: The Triangle has only one side completely inside the polygon. This is called a Terminal triangle.
            3: Special case where the polygon has only 3 sides
            
        Return value:
            - Number of side completely inside the polygon. Value between 0 and 3.
        """
        
        if self._nbr_internal is None:
            self._side_type = []
            self._nbr_internal = 0
            for i in range(3):
                p0 = self.coords_dual[i]
                p1 = self.coords_dual[i+1]
                if (_Triangle.line_segments.is_line_segment_present(p0, p1)):
                    self._side_type.append(_Triangle.SUPERIMPOSED)
                else:
                    self._side_type.append(_Triangle.INTERNAL)
                    self._nbr_internal += 1
        else:
            # nbr_internal has already been calculated... nothing to do
            pass
        
        return self._nbr_internal
        
    def get_category(self, minimal_width):
        """Determines if a triangle is a type bottleneck (to narrow) or a polygon other
        
        Parameters:
            - minimal_width: Float of the minimal width to check
            
        Return value:
            Tuple of 2 values
                - The category of the triangle: bottleneck or other
                - If the category is bottleneck the width of the bottle if it is below the minimal width
        
        """
        if self._category is None:
            
            if len(self._centre_lines) == 0:
                # If there is no center line the polygon is of type other
                self._category = ChordalAxisTransformer.OTHER
            else:
                nbr_internal = self.get_nbr_internal()
                if (nbr_internal == 2 and self._is_acute_triangle() ):
                    # The triangle is a sleeve triangle and is acute
                    # In this case the height of the triangle is also the width of the bottleneck
                    if self._height < minimal_width:
                        neighbour = True
                        self._width = self._height
                    else:
                        neighbour = False
                else:
                    # Check through a spatial search if there are any neighbous
                    (neighbour, self._width) = _Triangle.line_segments.check_chordal_axis(minimal_width/2., self._mid_triangle)
                if (neighbour):
                    extremity = False
                    for i, type in enumerate(self._side_type):
                        if (type == _Triangle.INTERNAL):
                            coord0 = self.coords_dual[i]
                            coord1 = self.coords_dual[i+1]
                            # Check if the triangle is located near the extremity of the polygon
                            # A triangle near an extremity of a polygon is not considered as a bottleneck
                            extremity = extremity or _Triangle.perimeters.is_extremity(coord0, coord1, minimal_width)
                    if extremity:
                        self._category = ChordalAxisTransformer.OTHER
                    else:
                        self._category = ChordalAxisTransformer.BOTTLENECK
                else:
                    self._category = ChordalAxisTransformer.OTHER
                
            if self._category != ChordalAxisTransformer.BOTTLENECK:
                self._width = None
                            
        return self._category, self._width
                    

class _LineSegments(object):
    """This class allows to make specific search in a list of line segment 
    
    """
    
    ID_RING = 'id_ring'
    
    def __init__(self, rings, search_tolerance):
        """Load the line segment in the spatial container in the form of 2 vertice line segment
        
        Parameters: 
            - rings: List of MA_LineString to store in the spatial container
            - search_tolerance: Search tolerance can vary depending on the dynamic of the data set from lat-lon to Lambert conformal
            
        Return value: None
        
        """
        
        self._search_tolerance = search_tolerance
        self._s_container = SpatialContainer()
        
        for id_ring, ring in enumerate(rings):
            lst_coords = list(ring.coords)
            nbr_coord_ring = len(lst_coords)
            # Split the line string into line segment of 2 vertice
            for j in xrange(nbr_coord_ring-1):
                p0 = lst_coords[j]
                p1 = lst_coords[j+1]
                line = MA_LineString([p0,p1], dual=None) # Dual is False as there are so many LineSegment
                del line.ma_properties  # Save some space as there are so many LineSegment
#                line.ma_properties[_LineSegments.ID_RING] = id_ring
                self._s_container.add_feature(line)
                        
    def is_line_segment_present(self, p0, p1):
        """Check if a line segment is located at a specific place
        
        Parameters: p0, p1: Tuple of x,y coordinates
        
        Return value
            - boolean: True: There is a line segement there
                       False: There is no line segment there
                       
        """
            
        # Find the mid point between p0 and p1
        mid_p0_p1 = GenUtil.mid_point(p0,p1)
        # Search for line segment there
        b_box = GenUtil.build_bounding_box(self._search_tolerance, mid_p0_p1)
        lines = self._s_container.get_features(bounds=b_box)
        
        present = False
        for line in lines:
            # Check if the first/vertice of the line correpond to the p1, p2 of the triangle if so 
            # there is exactly one line there
            line_coords = list(line.coords)
            if ( (GenUtil.distance(line_coords[0], p0) < GenUtil.ZERO and
                  GenUtil.distance(line_coords[1], p1) < GenUtil.ZERO ) or 
                 (GenUtil.distance(line_coords[0], p1) < GenUtil.ZERO and
                  GenUtil.distance(line_coords[1], p0) < GenUtil.ZERO ) ):  
                present = True
                break
            
        return present
    
    def check_chordal_axis(self, tolerance, target_coord):
        """Check at a specific location if a chordal axis is found within a certain tolerance 
        
        The chordal is calculated with the following steps:
            - Find all the line segment within a certain tolerance
            - Find the first closest line segment and the position p0 along the line
            - Find the second closest line segment and the position p1 along the line
            - p0 and p1 must form alomost a straight line
            - distance (p0,p1) must be below the tolerance
        
        Parameters:
            - tolerance: Search tolerance
            - target_coords: Search location, tuple of x,y coordinate
        """
        
        min_distance_0 = 1.0E+99
        min_distance_1 = 1.0E+99
        chord_distance = 1.0E+99
        
        # Extract the line segment near the search location
        b_box = GenUtil.build_bounding_box(tolerance, target_coord)
        lines = self._s_container.get_features(bounds=b_box)
        if (len(lines) == 1):
            # We need at least 2 lines otherwise there is no chrdal axes 
            lines = []
        ref_line = None
        for line in lines:
            # First pass find the closest line segment
            line_coords = list(line.coords)
            distance = GenUtil.distance_line_point(line_coords[0], line_coords[1], target_coord)
            if (distance < min_distance_0):
                min_distance_0 = distance
                ref_line = line
                    
        if ref_line is not None:
            # Extract the coordinate on the line
            target_point = Point(target_coord)
            distance_lr = ref_line.project(target_point)
            point_on_line_0 = ref_line.interpolate(distance_lr)
            coord_on_line_0 = point_on_line_0.coords[0]
            
            for line in lines:
                # Second pass find the second closest line 
                if line != ref_line:
                    line_coords = list(line.coords)
                    distance = GenUtil.distance_line_point(line_coords[0], line_coords[1], target_coord)
                    distance_lr = line.project(target_point)
                    point_on_line_1 = line.interpolate(distance_lr)
                    coord_on_line_1 = point_on_line_1.coords[0]
                    angle = GenUtil.compute_angle(coord_on_line_0, target_coord, coord_on_line_1)
                    # Check that p0 and p1 formed almost a straight line
                    if (angle > 120.):
                        distance = GenUtil.distance(coord_on_line_0, coord_on_line_1)
                        if (distance < chord_distance):
                            chord_distance = distance
                            
        # Check that the chord distance is within the tolerance
        if chord_distance < tolerance*2.:
            neighbours = True
        else:
            neighbours = False
            chord_distance = 1.0E+99
            
        return (neighbours,chord_distance)
    
    
class PolygonModifier(object):
    """Class used to modify polygon.  Polygon can be buffered eliminated and the topology is checked.
    
       Even if we manage polygon the reality is contour line string and they can overlap.  The 1000m contour can surround the
       1100m contour.  So for topology checking the polygons are broken into lines because contours can overlap. 
       If we were checking the topology with polygon there would be a topology error when a contour surrounds another contour  

    """
    
    # Define attribute name
    ACTION = 'action'
    BUFFER= 'buffer'
    ID_RING = 'id_ring'
    ID_POLYGON = 'id_polygon'
    INNER_RING = 'inner_ring'
    STATUS = 'status'
    
    # Define attribute value
    DONE = 'done'
    NOT_DONE = 'not done'
    TO_BE_DELETED = 'to_be_deleted'
    TO_BE_BUFFERED = 'to_be_buffered'
    
    def __init__(self, lst_tuple_id_polygon, test_sidedness=True, test_crossing_line=True):
        """Constructor to load the features
        
        Parameters: 
            - lst_tuple_id_polygon: List of tuple. the tuple contains
                                    - integer for the polygon id
                                    - MA_Polygon the polygon to modify
            - test_sidedness: Flag to enable/disable sidedness testing
            - test_crossing_line: Flag to enable/disable crossing line testing
        """
        
        # Public properties of the class
        self.nbr_widen = 0
        self.nbr_eliminated = 0
        self.nbr_conflict = 0
        
        # Private properties of the class
        self.params = Parameters()
        self.params.test_sidedness = test_sidedness
        self.params.test_crossing_line = test_crossing_line
        
        # Validating the list of tuple
        dict_ids = {}
        for tuple_id_polygon in lst_tuple_id_polygon:
            id = tuple_id_polygon[0]
            polygon = tuple_id_polygon[1]
            # Check the type of the polygon
            if not isinstance(polygon,MA_Polygon): 
                raise InternalError ("Not all polygons are of type MA_Polygon")
            else:
                # Check there are no duplicate id
                if dict_ids.has_key(id): 
                    raise InternalError ("Duplicate ID polygon")
                else:
                    dict_ids[id] = polygon
                    
        # Break the polygon into line strings in order to check the topology
        self._s_cont_pol_lines = SpatialContainer(line_opt_value=0)
        self._dict_pol_lines = {}
        self._dict_id_processed = {}
        self._dict_ma_properties = {}
        total_pol_lines = []
        for (id, polygon) in  dict_ids.iteritems():
            pol_lines = self._polygon_to_lines(id, polygon)
            total_pol_lines += pol_lines
            self._dict_pol_lines[id] = pol_lines[0]
            self._dict_ma_properties[id] = polygon.ma_properties
            polygon.ma_properties = []
        self._s_cont_pol_lines.add_features(total_pol_lines)
        
    def _polygon_to_lines(self, id_pol, polygon):
        """Transform a polygon into a list of line string
        
           The outer poygon always have a ring ID of 0.
           The outer ring always have a reference to the inner ring for performance issue
        
        Parameters:
            id_pol: ID of the polygon
            polygon: MA_Polygon to process
            
        Return value:
            List of MA_LineString
        """

        pol_lines = []
        pol_line = MA_LineString(list(polygon.exterior.coords), dual=False)
        pol_line.ma_properties[PolygonModifier.ID_RING] = 0
        pol_line.ma_properties[PolygonModifier.ID_POLYGON] = id_pol
        pol_lines.append(pol_line)            
        for i, pol_line in enumerate(polygon.interiors):
            pol_line = MA_LineString(list(pol_line.coords), dual=False)
            pol_line.ma_properties[PolygonModifier.ID_POLYGON] = id_pol
            pol_line.ma_properties[PolygonModifier.ID_RING] = i+1
            pol_lines.append(pol_line)
            
        pol_lines[0].ma_properties[PolygonModifier.INNER_RING] = pol_lines[1:]

        return pol_lines
    
    def _lines_to_polygon(self, outer_ring):
        """Transform a series of line into a polygon
        
        Parameter:
            - outer_ring: MA_LineString of the outer line of the polygon
            
        Return value:
            MA_Polygon built from the lines
            
        """
        
        inner_rings = outer_ring.ma_properties[PolygonModifier.INNER_RING]
        polygon = MA_Polygon(outer_ring, inner_rings)
        id_pol = outer_ring.ma_properties[PolygonModifier.ID_POLYGON]
        polygon.ma_properties[PolygonModifier.ID_POLYGON] = id_pol
        
        return polygon

    def _check_topology (self, polygon_modifier):
        """Check if the polygon to be modified will create a sidedness or a crossing line problem
        
        Paramaters: 
            - polygon_modifier: Polygon modification to applied to the original polygon 
        """
        
        line_simple_line = None # simple line constraint is not check as the convex hull is always simple
        line_crossing_line = LineString(list(polygon_modifier.exterior.coords))
        sidedness_polygon = None
        conflict_type = GenUtil.test_constraints (self, None, line_simple_line, line_crossing_line, 
                                                  sidedness_polygon, self._s_cont_pol_lines, [], False)
        
        return conflict_type
    
    def modify (self, id_polygon, lst_action_polygon):
        """Apply a list of modification (buffer or delete to one polygon)
        
        Because the modification of the polygon can lead to multipolygon, applying all the modification (add and substract)
        on the same polygon at the same time will make the code simpler. In consequence once a modification is done on a
        polygon, the polygon cannot be modify again.
        
        Parameters: 
            - id_polygon: ID of the polygon to modify
            - lst_action_polygon: List of polygon modifier to apply to the original polygon
        """
        
        if not self._dict_id_processed.has_key(id_polygon):
        
            if self._dict_pol_lines.has_key(id_polygon):
                
                # Rebuild the original polygon from the list of MA_LineString
                outer_ring = self._dict_pol_lines[id_polygon]
                inner_ring = outer_ring.ma_properties[PolygonModifier.INNER_RING]
                pol_lines = [outer_ring] + [inner_ring]
                edited_polygon = self._lines_to_polygon(self._dict_pol_lines[id_polygon])
                
                # Delete the features from the spatial container... they will be reinserted later
                self._s_cont_pol_lines.del_features(pol_lines)
                
                # Process each action modifier
                for action_polygon in lst_action_polygon:                    
                    if action_polygon.ma_properties[PolygonModifier.ACTION] == PolygonModifier.TO_BE_BUFFERED:
                        # The bottleneck needs to be buffered
                        buffer_pol = action_polygon.ma_properties[PolygonModifier.BUFFER]
                        # Checks if the buffer will create a topology error
                        conflict_type = self._check_topology(buffer_pol)
                        if (conflict_type is None):                            
                            edited_polygon = edited_polygon.union(buffer_pol)
                            self.nbr_widen += 1
                        else:
                            self.nbr_conflict += 1                                
                    else:
                        # The bottleneck needs to be deleted
                        # There is no topology to check when we delete a bottleneck
                        edited_polygon = edited_polygon.difference(action_polygon)
                        self.nbr_eliminated += 1
                        
                # Transform the multi polygon into a list of polygon 
                if isinstance(edited_polygon, Iterable):
                    edited_polygons = [ed_polygon for ed_polygon in edited_polygon]
                else:
                    edited_polygons = [edited_polygon]
                    
                for edited_polygon in edited_polygons:
                    ma_edited_polygon = MA_Polygon(edited_polygon.exterior, edited_polygon.interiors)
                    # Transform the MA_Polygon into MA_LineString
                    pol_lines = self._polygon_to_lines(id_polygon, ma_edited_polygon)
                    # Add the MA_LineString into the container
                    self._s_cont_pol_lines.add_features(pol_lines)
                    
                # Flag the polygon has processed... cannot be processed again
                self._dict_id_processed[id_polygon] = None
                
            else:
                # Unknown ID polygon
                pass
            
        else:
            # ID already processed we cannot apply 2 modification at the same polygon
            pass
        
    def get_polygons(self):
        """Extract the polygon features
        
        """
                
        polygons = []
        # Rebuild the polygons 
        for outer_ring in self._s_cont_pol_lines.get_features():
            # Only keep the outer ring
            if (outer_ring.ma_properties[PolygonModifier.ID_RING] == 0):
                polygon = self._lines_to_polygon(outer_ring)
                id_pol = polygon.ma_properties[PolygonModifier.ID_POLYGON]
                polygon.ma_properties = deepcopy(self._dict_ma_properties[id_pol])
                polygons.append(polygon)
                # Delete reference... to avoid possible circular references
                outer_ring.ma_properties[PolygonModifier.INNER_RING] = None
                
                # Get rid of the reference 
                outer_ring.ma_properties[PolygonModifier.INNER_RING] = None
            
        # Delete the content of the spatial container and reset everything so the object is uncallable...
        self._s_cont_pol_lines.del_features(self._s_cont_pol_lines.get_features())
        self._dict_pol_lines = {}
        self._dict_id_processed = {}
        self._dict_ma_properties = {}

        return polygons
    
class PerimeterDistance(object):
    """This class allows to calculate the distance between 2 coordinates on a line closed string

    
    Internal data structure are maintained to accelerate the computation 
    
    """
    
    _ID_RING = 'id_ring'
    _ID_COORD = 'id_coord'
    
    def __init__(self, rings, search_tolerance):
        """Load the internal structure with the ring information
        
        Parameters: 
           - rings: List of closed MA_LineString of a polygon
           - search_tolerance: Search tolerance can vary depending on the dynamic of the data set from lat-lon to Lambert conformal
        
        Return value: None
        """
        
        self._search_tolerance = search_tolerance
        
        for ring in rings:
            if ( isinstance(ring, LinearRing) ):
                pass
            else:
                raise Exception ("Can only work on LineString or MA_LineString")
            
        self.lst_cumm_distance = []
        self.lst_lst_coords = []
        self.s_cont_points = SpatialContainer()
        
        for id_ring, ring in enumerate(rings):
            lst_coords = list(ring.coords)
            self._init_load_points(id_ring, lst_coords)
            self._init_load_cum_distance(lst_coords)
        
    def _init_load_points(self, id_ring, lst_coords):
        """Load all the point in a spatial container for fast point search
        
        Parameters: 
            - id_rings: Ring ID 
            - lst_coords: List of coordinates to load in the spatial container
            
        Return value:
            - None
        
        """

        nbr_coords = len(lst_coords)
        for i in range(1,nbr_coords):
            point = MA_Point(lst_coords[i])
            point.ma_properties[PerimeterDistance._ID_RING] = id_ring
            point.ma_properties[PerimeterDistance._ID_COORD] = i
            self.s_cont_points.add_feature(point)

    def _init_load_cum_distance(self, lst_coords):
        """Create a list of cummulative distance for a list of coordinates of a closed line
        
        Parameters: 
            - lst_coords: List of tuple of x,y coordinate
            
        Return value: None
        """
                
        cum_distance = [0.] 
        nbr_coords = len(lst_coords)
        
        # Build the cumulative distance list
        for i, current_coord in enumerate(lst_coords):
            if (i==0): 
                # Do not add the first vertex in the container
                previous_coord = current_coord
            else:
                dist = GenUtil.distance(previous_coord, current_coord)
                last_dist = cum_distance[-1]
                cum_distance.append(last_dist+dist)
                previous_coord = current_coord
                
        self.lst_cumm_distance.append(cum_distance)
        self.lst_lst_coords.append(lst_coords)
        
    def _get_points_info(self, coord, raise_exception=True):
        """Extract the information for a specific coordinate
        
        Paramneter:
            - coord: Location to search as a x,y coordinate
            - raise_exception: Flag to enable (True) or disable (False) an exception if nothing is found
            
        Return value:
            - Tuple of 2 elements:
                - integer identifying the ring id of the ring
                - integer identifying the position (index) of the coordinate in the list of coordinate of the ring
        """
        
        b_box = GenUtil.build_bounding_box(self._search_tolerance, coord)
        points = self.s_cont_points.get_features(bounds=b_box)
        nbr_points = len(points)
        if (nbr_points == 0):
            # Nothing is found
            id_ring = -1
            id_coord = -1
            if (raise_exception):
                raise Exception ("Integrity problem at coordinate: (%f,%f)" %(coord[0],coord[1]) )
        else:
            if (nbr_points == 1):
                # There is only onpoint
                point = points[0]
            else:
                # Take the closest point
                min_dist = 1.0E+99
                for p in points:
                    dist =  GenUtil.distance(p.coords_dual[0], coord)
                    if (dist < min_dist):
                        point = p
                        dist = min_dist
            id_ring = points[0].ma_properties[PerimeterDistance._ID_RING]
            id_coord = points[0].ma_properties[PerimeterDistance._ID_COORD]

        return (id_ring, id_coord) 
    
    def get_sub_perimeter(self, coord0, coord1):
        """Extract the coordinate between coord0 and coord1
        
           If the 2 coordinates belongs to the same ring extract the smallest perimeter coordinates list 
           between the 2 coordinates.
        
        Parameters:
            - coord0: First coordinate as a x,y Tuple
            - coord1: Second coordinate as a x,y Tuple

        """
        
        id_ring0, id_coord0 = self._get_points_info(coord0)
        id_ring1, id_coord1 = self._get_points_info(coord1)
        
        if (id_ring0 != -1 and id_ring0 != -1):
            
            if (id_ring0 == id_ring1):
                # Extract the appropriate coordinate ring
                lst_coords = self.lst_lst_coords[id_ring0] 
                sub_coords1 = []
                sub_coords2 = []
                nbr_coords = len(lst_coords)
                
                # Loop to extract coordinate from first to last
                start,end = id_coord0, id_coord1
                i = start
                while (i != end):
                    sub_coords1.append(lst_coords[i])
                    i = (i+1)%nbr_coords
                sub_coords1.append(lst_coords[i])
                
                # Loop to extract coordinate from last to first
                start,end = id_coord1, id_coord0
                i = start
                while (i != end):
                    sub_coords2.append(lst_coords[i])
                    i = (i+1)%nbr_coords
                sub_coords2.append(lst_coords[i])
                
                # Take the smallest list
                if (len(sub_coords1) < len(sub_coords2)):
                    sub_lst_coords = sub_coords1
                else:
                    sub_lst_coords = sub_coords2
            
            else:
                # Coordinates belongs to different ring
                sub_lst_coords = None
                
        else:
            # One or both coordinate not found
            sub_lst_coords = None
            
        return sub_lst_coords
        
    def is_extremity(self, coord0, coord1, minimal_width):
        """Determine if the 2 coordinates is located on an extremity of the polygon or not on an extremity
        
           To be considered as an extremity the 2 points must meet the following 2 conditions:
               - It must be on the same ring
               - The smallest of the areas formed by coord0, coord1 cutting the polygon is smaller
                 than an area formed by a rectangle where one side is of length minimal_width and the other
                 side of length coord0,coord1  
           
           Parameters:
               - coord1,coord2: Tuple of x,y coordinate of the location to check on the perimeter
               - minimal_width: Tolerance to determine if it is an extremity
               
           Return value:
               - boolean: True it is an extremity; False otherwise
        """
        
        id_ring0, id_coord0 = self._get_points_info(coord0)
        id_ring1, id_coord1 = self._get_points_info(coord1)
        
        if (id_ring0 != -1 and id_ring0 != -1):
            
            if (id_ring0 == id_ring1):
                # Coordinates are on the same ring
                cumm_distance = self.lst_cumm_distance[id_ring0]
                if ( (id_coord0 < 0 or id_coord0 > len(cumm_distance))  or 
                     (id_coord1 < 0 or id_coord1 > len(cumm_distance)) ):
                    raise Exception ("Internal Error...")
                
                if (id_coord0 < id_coord1):
                    i,j = id_coord0, id_coord1
                else:
                    i,j = id_coord1, id_coord0
                    
                # Extracting the 2 perimetre formed by the line (p1,p2) cutting the polygon in two parts 
                peri_distance_1 = cumm_distance[j] - cumm_distance[i]
                peri_distance_2 = cumm_distance[-1] - peri_distance_1    
                peri_distance = min(peri_distance_1,peri_distance_2)
                
                # Check if the area formed is smaller than a minimal_width area 
                dist_0_1 = GenUtil.distance(coord0, coord1)
                if dist_0_1 + peri_distance < 2*dist_0_1 + 2*minimal_width:
                    extremity = True
                else:
                    extremity = False
            else:
                # Different rings
                extremity = False
        else:
            # Some coordinates not found
            extremity = False

        return extremity

class PolygonMerger(object):
    """Class that merged polygons together and add in the merged polygon a list properties of the original merged polygon    
    """

    _REF_POLYGON = 'ref_polygon'
    
    def __init__(self, lst_polygons):
        """Constructor of the PolygonMerger class
        
        Parameters: lst_polygons: List of MA_Polygon
        
        Return value: None
        
        """
        
        for polygon in lst_polygons:
            if not isinstance(polygon,MA_Polygon):
                raise InternalError ("Not all polygons are of type MA_Polygon")
        
        self.lst_polygon = lst_polygons
        
    def getMerged(self, lst_attributes=None):
        """Merged the polygon together 
     
           Using shpely cascadeunion all the polygons are merged.
           For all the attribute to keep on the merged polygon 
              If the orginal polygon is inside the union polygon than add the attribute to the resulting polygon
        
        Parameters: List of attributes to keep (as a list of attributes) on the resulting merged polygons
        
        Return values:
            List of MA_Polygon
        """
        
        if (len(self.lst_polygon) >= 1):
            # Merge the polygons
            merged_pol = cascaded_union(self.lst_polygon)
            
            if not (isinstance(merged_pol, Iterable)):
                merged_pol = [merged_pol]
            
            #Transoform the multipolygon into a list of MA_Polygon
            lst_merged_pol = [MA_Polygon(pol.exterior, pol.interiors) for pol in merged_pol]
            
            if lst_attributes is not None:
                # For performance reason creates a spatial container         
                s_cont_input_pol = SpatialContainer()
                # Add the original polygon in the spatial container
                for polygon in self.lst_polygon:
                    point_inside = polygon.representative_point()
                    ma_point = MA_Point([point_inside.coords[0]])
                    ma_point.ma_properties[PolygonMerger._REF_POLYGON] = polygon
                    s_cont_input_pol.add_feature(ma_point)
                
                # Process each resulting polygon
                for merged_pol in lst_merged_pol:
                    # Process each attribute
                    for attribute in lst_attributes:
                        merged_pol.ma_properties[attribute] = []
                    
                    # First pass bounding box search
                    potential_pols = s_cont_input_pol.get_features(bounds=merged_pol.bounds)
                    for potential_pol in potential_pols:
                        # Second pass spatial search
                        if (merged_pol.contains(potential_pol)):
                            for attribute in lst_attributes:
                                if (attribute in potential_pol.ma_properties[PolygonMerger._REF_POLYGON].ma_properties):
                                    # Add the attribute
                                    value = potential_pol.ma_properties[PolygonMerger._REF_POLYGON].ma_properties[attribute]
                                else:
                                    # Add empty if no attribute
                                    value = None
                                merged_pol.ma_properties[attribute].append(value)
                                
        else:
            # There is no polygon to process
            lst_merged_pol = []
            
        return lst_merged_pol