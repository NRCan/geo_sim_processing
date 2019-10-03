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
    def distance(p1, p2):
        """Calculate the euclidian distance between 2 points

        *Parameters*:
            - p1: (x,y) tuple of the first coordinate
            - p2: (x,y) typle of the second coordinate

        *Returns*:
            - Distance between the 2 points (real)

        """

        return (math.sqrt((p2[0] - p1[0]) ** 2.0 + (p2[1] - p1[1]) ** 2.0))


    @staticmethod
    def orientation(p0, p1, p2):
        """ Calculate the orientation (clockwise or anticlockwise) of a line formed by 3 vertices using the dot product

        Parameters:
            p0, p1, p2: Three (x,y) coordinates tuple

        Return value
            float the direction of the line an
                0 : Straight line
                1: Counter clockwise angle
                -1 : Clockwise angle

        """

        orient =  ((p0[0] - p1[0]) * (p2[1] - p1[1])) - ((p2[0] - p1[0]) * (p0[1] - p1[1]))

        if orient > 0.:
            orient = 1
        elif orient< 0.:
            orient = -1
        else:
            orient = 0

        return orient


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
            if (new_bbox[0] != old_bbox[0] or
                    new_bbox[1] != old_bbox[1] or
                    new_bbox[2] != old_bbox[2] or
                    new_bbox[3] != old_bbox[3]):
                is_the_same = False
            else:
                is_the_same = True
        else:
            if (len(new_lst_bbox) == 1):
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
                xmin, ymin = old_lst_bbox[i_bbox][0], old_lst_bbox[i_bbox][1]
                xmax, ymax = old_lst_bbox[i_bbox][2], old_lst_bbox[i_bbox][3]
                try:
                    for coord in line_coords:
                        if (not (xmin <= coord[0] <= xmax and ymin <= coord[1] <= ymax)):
                            # Try to find the next bbox that contains the coordinate
                            i_bbox += 1
                            if (i_bbox < len_old_bbox):
                                xmin, ymin = old_lst_bbox[i_bbox][0], old_lst_bbox[i_bbox][1]
                                xmax, ymax = old_lst_bbox[i_bbox][2], old_lst_bbox[i_bbox][3]
                                if (not (xmin <= coord[0] <= xmax and ymin <= coord[1] <= ymax)):
                                    # The coordinate is outside all bounding boxes
                                    raise Exception
                    is_the_same = True
                except Exception:
                    is_the_same = False
                except:
                    raise ("Unknown error...")

        return is_the_same

    def _adjust_bounding_box(self, bounds):
        """Modify the bounds of a feature when the bounds of almost zero

        *Parameters*:
            - bounds: Tuple of a bounding box (xmin, ymin, xmax, ymax)

        *Returns*:
            - Tuple of a bounding box (xmin, ymin, xmax, ymax)

        """

        if (abs(bounds[2] - bounds[0]) < GenUtil.ZERO or abs(bounds[3] - bounds[1]) < GenUtil.ZERO):
            if abs(bounds[2] - bounds[0]) >= GenUtil.ZERO:
                xmin = bounds[0]
                xmax = bounds[2]
            else:
                xmin = bounds[0] - GenUtil.ZERO
                xmax = bounds[2] + GenUtil.ZERO
        if abs(bounds[3] - bounds[1]) >= GenUtil.ZERO:
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
        # if len(lst_bounds) == 1:
        # Presently there is a little bug in RTree when the xmin and xmax or ymin and ymax are the same value
        # The problem is resolved when we add a very small delta between the 2 values
        # This bug is supposed to be solved in the next release. We're running now on 0.6
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
            raise GenException('Unsupported feature type...')

        # Check if the feature is already in a spatial container
        if hasattr(feature, "_gbt_sc_id"):
            raise GenException('Feature is already in a spatial container')

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
                    # Delete the property _sci_id
                    del feature._sci_id
                except:
                    raise InternalError("Internal corruption, problem with the container and/or the RTree")
            else:
                raise InternalError("Internal corruption, key {} has disappear...".format(feature._sci_id))

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

    def update_spatial_index(self, feature):
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
            # Add the new bounding boxes in Rtree
            self._r_tree.add(feature.__gbt_sci_id, new_bbox)

            # Save the bounding boxes
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
        keys = list(set(keys) - set(keys_to_remove))

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


class GenException (Exception):
    """
    This is the base exception class for genmetal algorithms
    """

    def __init__(self, *arguments, **keywords):
        Exception.__init__(self, *arguments, **keywords)


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