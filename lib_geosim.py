#!/usr/bin/env python
# -=- encoding: utf-8 -=-

"""
General classes and utilities needed for the GeoSim.

"""

import sys
from math import atan, degrees, sqrt, acos, pi
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import linemerge
from shapely.ops import unary_union
from collections.abc import Iterable
from collections import OrderedDict
import fiona

try:
    from rtree import Rtree

    lib_Rtree = True
except:
    # If Rtree is not installed
    lib_Rtree = False
    from shapely.strtree import STRtree


class LineStringSc(LineString):
    """LineString specialization that allow LincestringSc to be included  in the SpatialContainer"""

    def __init__(self, coords):
        """ Constructor for the LineStringSc

            Parameters
            ----------
            coords : tuple
                Tuple of x,y coordinate

            Returns
            -------
            None
        """
        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        """ Return the coordinate of the line

        Parameters
        ----------
        None

        Returns
        -------
        Tuple
            List of x,y coordinate

        """
        return super().coords

    @coords.setter
    def coords(self, coords):
        """Update the coordinate value anf the spatial container if there is a spatial container

        Parameters
        ----------
        Coords : tuple
            List of x,y coordinates

        Returns
        -------
        None
        """
        LineString.coords.__set__(self, coords)

        if self._sc_scontainer != None:  # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_scontainer.update_spatial_index(self)


class PointSc(Point):
    """Point specialization that allow PointSc features to be included  in the SpatialContainer"""

    def __init__(self, coords):
        """ Constructor of the PointSc class

        Parameters
        ----------
        coords : tuple
            x,y coordinates of the point

        Returns
        -------
        None

        """

        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        """ Return the coordinate of a line

        Parameters
        ----------
        None

        Returns
        -------
        Tuple
            x,y coordinate of the point

        """

        return super().coords

    @coords.setter
    def coords (self, coords):
        """Update the coordinate of the LineString and update the spatial container if the spatial container exists

        Parameters
        ----------
        coords : tuple
            x,y coordinates of the point

        Returns
        -------
        None

        """

        Point.coords.__set__(self, coords)

        if self._sc_scontainer is not None:  # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_container.update_bbox(self)


class PolygonSc(Polygon):
    """Polygon specialization to be included in the SpatialContainer"""

    def __init__(self, exterior, interiors=None):
        """Constructor of the PolygoncSc

        Parameters
        ----------
        exterior : list
            List of x,y coordinate forming a loop
        interiors : list
            List of loops

        Returns
        -------
        None

        """
        super().__init__(exterior, interiors)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def exterior(self):
        """Return the exterior description

        Parameters
        ----------
        None

        Returns
        -------
        list
            List of x,y coordinate forming a loop

        """
        return super().exterior

    @property
    def interiors(self):
        """Return the interior description

        Parameters
        ----------
        None

        Returns
        -------
        list
            List of loops formint the interior

        """

        return super().interiors

    @exterior.setter
    def exterior(self, exterior):
        """Cannot set the interior raise exception

        Parameters
        ----------
        exterior : list
            List of xy coordinate forming a loop

        Returns
        -------
        Exception

        """

        raise GeoSimException("Cannot update the exterior coordinates of a polygon")

    @interiors.setter
    def interiors(self, interiors):
        """Cannot set the exterior raise exception

        Parameters
        ----------
        interiors : list
            List of interiors forming loops

        Returns
        -------
        Exception

        """


        raise GeoSimException("Cannot update the interior coordinates of a polygon")


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
    def make_iterable(iter_feature):
        """Test if the parameter is iterable; if not make it iterable by creating a tuple of one element

        Parameters
        ----------
        iter_feature: Object
            Object to test for iterability

        Returns
        -------
        Iterable object
        """

        if not isinstance(iter_feature, Iterable):
            iter_feature = (iter_feature,)

        return iter_feature

    @staticmethod
    def difference_angle_vector(p0, p1, zero_tolerance):
        """Calculate the angle between the 2 vectors

        Parameters
        ----------
        p0 : tuple
            x,y coordinate of the first vector (center to 0,0)
        p1 : tuple
            x,y coordinate of the first vector (center to 0,0)
        zero_tolerance : float
            Value for zero approximation

        Returns
        -------
        float
            angle between [0..360]

        """

        x0, y0 = p0[0], p0[1]
        x1, y1 = p1[0], p1[1]

        delta_y = y1 - y0
        delta_x = x1 - x0
        if abs(delta_x) <= zero_tolerance:  # Avoid division by zero
            delta_x = zero_tolerance
        angle = degrees(atan(delta_y / delta_x))

        # Apply a correction for the quadrant; in order to obtain an angle betweer 0..360
        if delta_x >= 0 and delta_y >= 0:
            # Quadrant 1 nothing to do
            pass
        elif delta_x < 0 and delta_y >= 0:
            # Quadrant 2
            angle += 180.
        elif delta_x < 0 and delta_y < 0:
            # Quadrant 3
            angle += 180.
        else:
            # delta_x > 0 anf delta_y < 0
            # Quandrant 4
            angle += 360.

        return angle


    @staticmethod
    def distance(p1, p2):
        """Calculate the euclidean distance between 2 points

        Parameters
        ----------
        p0 : tuple
            x,y coordinate of the first vector (center to 0,0)
        p1 : tuple
            x,y coordinate of the second vector (center to 0,0)

        Returns
        -------
        float
            Distance between the 2 points (real)

        """

        return sqrt((p2[0] - p1[0]) ** 2.0 + (p2[1] - p1[1]) ** 2.0)

#    @staticmethod
#    def compute_angle(p1, p2, p3, type=DEGREE):
#        """
#        Function to calculate angle between two vectors.
#        """
#
#        return GenUtil.angle_vector(p1, p2, p3, type)

    @staticmethod
    def angle_vector(p1, p2, p3, type=DEGREE):
        """Calculate the angle formed by the vector p1-p2 and p2-p3

        Parameters
        ----------
        p1 : tuple
            x,y coordinate of the first coordinate
        p2 : tuple
            x,y coordinate of the second coordinate
        p3 : tuple
            x,y coordinate of the third coordinate

        Returns
        -------
        float
            Angle between the vector p1-p2 and p2-p3 [0..180]

        """

        a = (p2[0] - p1[0], p2[1] - p1[1])
        b = (p2[0] - p3[0], p2[1] - p3[1])
        len_a = (a[0] ** 2. + a[1] ** 2.) ** .5
        len_b = (b[0] ** 2. + b[1] ** 2.) ** .5

        dot_p = a[0] * b[0] + a[1] * b[1]

        # If P1 == P2 or P2 == P3 ===> angle is 180.
        if len_a * len_b != 0.0:
            value = dot_p / (len_a * len_b)
            if value >= 1.0:
                value = 1.0
            if value <= -1.0:
                value = -1.0
        else:
            value = -1.0

        theta = acos(value)

        if type == GenUtil.DEGREE:
            theta = degrees(theta)

        return theta

    @staticmethod
    def orientation(p0, p1, p2):
        """ Calculate the orientation (clockwise or anticlockwise) of a line formed by 3 vertices using the dot product

        Parameters
        ----------
        p1 : tuple
            x,y coordinate of the first coordinate
        p2 : tuple
            x,y coordinate of the second coordinate
        p3 : tuple
            x,y coordinate of the third coordinate

        Returns
        -------
        int
            The direction of the line:
                0 : Straight line
                1: Counter clockwise angle
                -1 : Clockwise angle

        """

        orient = ((p0[0] - p1[0]) * (p2[1] - p1[1])) - ((p2[0] - p1[0]) * (p0[1] - p1[1]))

        if orient > 0.:
            orient = 1
        elif orient < 0.:
            orient = -1
        else:
            orient = 0

        return orient

    @staticmethod
    def rescale_vector(p1, p2, scale_factor):
        """Rescale the vector defined by the points P1 and P2 by a factor

        Parameters
        ----------
        p1 : tuple
            x,y coordinate of the first coordinate
        p2 : tuple
            x,y coordinate of the second coordinate
        scale_factor : real
            Factor to scale the vector

        Returns
        -------
        tuple
            x,y coordinate of the rescale vector

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
    def mid_point(p1, p2):
        """Return a point in the middle of the 2 points

        Parameters
        ----------
        p1 : tuple
            x,y coordinate of the first coordinate
        p2 : tuple
            x,y coordinate of the second coordinate

        Returns
        -------
        tuple
            x,y coordinate of the milddle point
        """

        x = (p1[0] + p2[0]) / 2.
        y = (p1[1] + p2[1]) / 2.

        return (x, y)

    @staticmethod
    def calculate_compactness_index(area, perimeter):
        """Calculate the compactness index based of the perimeter and area

        Parameters
        area : float
            Area of the polygon
        perimeter : float
            Perimeter of the area

        Return:
        float
            Compactness index

        """

        return 4 * area * pi / (perimeter ** 2.0)

    @staticmethod
    def build_bounding_box(tolerance, coord):
        """Create and adjust a bounding box (xmin, ymin, xmax, ymax) with a small tolerance

        Parameters
        ----------
        tolerance : float
            Delta value to add to the bounding box
        coord : tuple
            Two x,y coordinate defining the bounding box

        Returns
        -------
        tuple
            Two x,y coordinate defining the bounding box

        """

        xmin = coord[0] - tolerance
        ymin = coord[1] - tolerance
        xmax = coord[0] + tolerance
        ymax = coord[1] + tolerance

        return (xmin, ymin, xmax, ymax)

    @staticmethod
    def calculate_adjusted_area(area, cmp_index):
        """Calculate the adjusted area from the area and compactness index

        Parameters
        ----------
        area : float
            Area of the polygon
        cmp_index : float
            Compactness index of the areea

        Return:
        float
            Adjusted area of the polygon

            """
        return area * (0.75 / cmp_index)

    @staticmethod
    def read_in_file(in_file, geo_content, layer_in=None):
        """Read and load the vectors features in the input file

        Parameters
        ----------
        in_file : lst
            List of name of the input file layers (geopackage)
        geo_content : dict
            Dictionary containing the information of the spatial database (including the features)
        layer_in : str
            Name of the layer to read

        Returns
        -------
        None

        """

        if layer_in is None:
            # Extract the name of the layers in the file
            geo_content.layer_names = fiona.listlayers(in_file)
        else:
            # Only extract specific layer
            geo_content.layer_names = layer_in

        # extract the spatial feature in the file
        for layer_name in geo_content.layer_names:
            with fiona.open(in_file, 'r', layer=layer_name) as src:
                geo_content.crs = src.crs
                geo_content.driver = src.driver
                geo_content.schemas[layer_name] = src.schema
                geo_content.bounds.append(src.bounds)

                for in_feature in src:
                    geom = in_feature['geometry']
                    if geom['type'] == 'Point':
                        feature = PointSc(geom['coordinates'])
                    elif geom['type'] == 'LineString':
                        feature = LineStringSc(geom['coordinates'])
                    elif geom['type'] == 'Polygon':
                        exterior = geom['coordinates'][0]
                        interiors = geom['coordinates'][1:]
                        feature = PolygonSc(exterior, interiors)
                    else:
                        print("The following geometry type is unsupported: {}".format(geom['type']))
                        feature = None
                    if feature is not None:
                        feature.sb_layer_name = layer_name  # Layer name is the key for the schema
                        feature.sb_properties = in_feature['properties']
                        geo_content.in_features.append(feature)
            src.close()

    @staticmethod
    def write_out_file(out_file, geo_content):
        """Write the vectors in the output file

        Parameters
        ----------
        out_file :str
            Name of the output file (geopackage)
        geo_content : DataClass
            Contains information needed to create the spatial database

        Returns
        -------
        None

        """

        # Loop over each layer and write the content of the file
        for layer_name in geo_content.layer_names:
            with fiona.open(out_file, 'w',
                            driver=geo_content.driver,
                            layer=layer_name,
                            crs=geo_content.crs,
                            schema=geo_content.schemas[layer_name]) as dest:
                out_features = []
                for feature in (feature for feature in geo_content.out_features
                                if feature.sb_layer_name == layer_name):
                    # Transform the Shapely features for fiona writing
                    if feature.geom_type == GenUtil.POINT:
                        coordinates = (feature.x, feature.y)
                        geo_content.out_nbr_points += 1
                    elif feature.geom_type == GenUtil.LINE_STRING:
                        coordinates = list(feature.coords)
                        geo_content.out_nbr_line_strings += 1
                    elif feature.geom_type == GenUtil.POLYGON:
                        exterior = list(feature.exterior.coords)
                        interiors = [list(interior.coords) for interior in feature.interiors]
                        coordinates = [exterior] + interiors
                        geo_content.out_nbr_polygons += 1
                        geo_content.out_nbr_holes += len(interiors)

                    out_feature = {'geometry': {'type': feature.geom_type,
                                                'coordinates': coordinates},
                                   'properties': feature.sb_properties}
                    out_features.append(out_feature)

                dest.writerecords(out_features)

            dest.close()

    @staticmethod
    def write_out_file_append(out_file, geo_content):
        """Write the vectors in the output file

        Parameters
        ----------
        out_file : str
            Name of the output file (geopackage)
        geo_content : DataClass
            Contains information to create the spatial database

        Returns
        -------
        None

        """

        line_schema = landmarks_schema = {'geometry': 'LineString',
                                          'properties': OrderedDict([])
                                          }

        # Loop over each layer and write the content of the file
        for layer_name in geo_content.layer_names:
            with fiona.open(out_file, 'w',
                            driver=geo_content.driver,
                            layer=layer_name,
                            crs=geo_content.crs,
                            schema=line_schema) as dest:
                out_features = []
                for feature in (feature for feature in geo_content.out_features
                                if feature.sb_layer_name == layer_name):
                    # Transform the Shapely features for fiona writing
                    if feature.geom_type == GenUtil.POINT:
                        coordinates = (feature.x, feature.y)
                        geo_content.out_nbr_points += 1
                    elif feature.geom_type == GenUtil.LINE_STRING:
                        coordinates = list(feature.coords)
                        geo_content.out_nbr_line_strings += 1
                    elif feature.geom_type == GenUtil.POLYGON:
                        exterior = list(feature.exterior.coords)
                        interiors = [list(interior.coords) for interior in feature.interiors]
                        coordinates = [exterior] + interiors
                        geo_content.out_nbr_polygons += 1
                        geo_content.out_nbr_holes += len(interiors)

                    out_feature = {'geometry': {'type': feature.geom_type,
                                                'coordinates': coordinates},
                                   'properties': feature.sb_properties}
                    out_features.append(out_feature)

                dest.writerecords(out_features)

            dest.close()


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
    _sc_id = 0

    def __init__(self):
        """Create an object of type SpatialContainer

        The init will create one container for the feature a dictionary and one
        container for the spatial index (Rtree)

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self._r_tree = Rtree()  # Container for the Rtree
        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def adjust_bbox(self, bounds, delta=GenUtil.ZERO):
        """Adjust the bounding box by increasing by a very small delta

        Parameters
        ----------
        bounds: Tuple
            Values forming the bounding box of the feature (xmin, ymin, wmax, ymax)

        Returns
        -------
        Tuple
            Altered bounding box (xmin, ymin, wmax, ymax)
        """

        xmin, ymin, xmax, ymax = bounds

        xmin -= delta
        ymin -= delta
        xmax += delta
        ymax += delta

        return (xmin, ymin, xmax, ymax)

    def add_feature(self, feature):
        """Adds a feature in the container and update the spatial index with the feature's bound

        To be added in the container a spatial feature must be a PointSc, LineStringSc or PolygonSc.

        Parameters
        ----------
        feature: PointSc, LineStringSc or PolygonSc
            A spatial feature to insert into the SpatialContainer

        Returns
        -------
        None
        """

        # Check if the type is valid
        if issubclass(feature.__class__, (PointSc, LineStringSc, PolygonSc)):
            pass
        else:
            raise GeoSimException('Unsupported class: {}'.format(str(feature.__class__)))

        bounds = feature.bounds

        # Adjust the bounding box
        bounds = self.adjust_bbox(bounds)

        # Container unique internal counter
        SpatialContainer._sc_id += 1

        # Add the spatial id to the feature
        feature._sc_id = SpatialContainer._sc_id
        feature._sc_scontainer = self

        # Add the feature in the feature container
        self._features[feature._sc_id] = feature

        # Add the bounding box in the bbox_container
        self._bbox_features[feature._sc_id] = bounds
        self._r_tree.add(feature._sc_id, bounds)

        return

    def add_features(self, features):
        """Adds a list of feature in the spatial container and

        Parameters
        ----------
        feature: List
            A list spatial feature derived from PointSc, LineStringSc or PolygonSc

        Returns
        -------
        None
        """

        for feature in features:
            self.add_feature(feature)

        return

    def del_feature(self, feature):
        """Delete the feature in the spatial container and in the RTree.

        If the feature is included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...

        Parameters
        ----------
        feature: PointSc, LineStringSc, PolygonSc
            The feature to delete in the spatial container

        Returns
        -------
        None
        """

        ret_value = 0

        # Check if the feature has a container_key
        if hasattr(feature, "_sc_id"):

            if (feature._sc_id in self._features and
                    feature._sc_id in self._bbox_features):

                try:
                    # Retrieve the bounding boxes of this feature
                    bbox = self._bbox_features[feature._sc_id]
                    # Delete the feature from the features and the bbox_features
                    del self._features[feature._sc_id]
                    del self._bbox_features[feature._sc_id]
                    # Delete the different bounds in RTree
                    self._r_tree.delete(feature._sc_id, bbox)
                    # Reset the property _sc_id and _sc_scontainer
                    feature._sc_id = None
                    feature._sc_scontainer = None
                except:
                    raise InternalError("Internal corruption, problem with the container and/or the RTree")
            else:
                raise InternalError("Internal corruption, key {} has disappear...".format(feature._sc_id))

        return ret_value

    def del_features(self, features):
        """Delete a list of features in the spatial container

        If the features are included in the spatial container the feature is deleted;
        if the feature is not included in the spatial container... nothing happen...
        Exception InternalError: If the key is not in one of the structure

        Parameters
        ----------
        features List
            The list of features to delete

        Returns
        -------
        List
            Value for each feature to delete. 0 : deleted; 1 : not delete
        """

        ret_value = []

        for feature in features:
            ret_value.append(self.del_feature(feature))

        return ret_value

    def update_spatial_index(self, feature):
        """Update the bounds of the feature in the spatial index

        It will only modify the Rtree spatial index if the bounding
        box of the feature is changed in comparison with the old one.

        Parameters
        ----------
        feature: PointSc or LineStringSc or PolygonSc
            Feature containing the bounds to update

        Returns
        -------
        None
        """

        old_bbox = self._bbox_features[feature._sc_id]
        new_bbox = feature.bounds
        old_x_min, old_y_min, old_x_max, old_y_max = old_bbox[0], old_bbox[1], old_bbox[2], old_bbox[3]
        new_x_min, new_y_min, new_x_max, new_y_max = new_bbox[0], new_bbox[1], new_bbox[2], new_bbox[3]

        if old_x_min <= new_x_min and \
                old_y_min <= new_y_min and \
                old_x_max >= new_x_max and \
                old_y_max >= new_y_max:
            # Nothing to do new bounding box is completely included into the old one
            pass
        else:
            # Adjust the bounding box
            new_bbox = self.adjust_bbox(new_bbox)
            # Delete The old bounding box in Rtree
            self._r_tree.delete(feature._sc_id, old_bbox)
            # Add the new bounding boxes in Rtree
            self._r_tree.add(feature._sc_id, new_bbox)

            # Save the bounding boxes
            self._bbox_features[feature._sc_id] = new_bbox

        return

    def get_features(self, bounds=None, remove_features=[]):
        """Extract the features from the spatial container.

        According to the parameters the extraction can manage the extraction based on a bounding box using
        the spatial index RTree, some filters to manage extraction based on properties and the possibility
        to remove specific features based on a list of keys

        Parameters
        ----------
        bounds: tuple
            Bounding for the spatial extraction. *None* means all the features
        remove_keys: List
            List of ID of the keys to be removed from the selection

        *Returns*:
            - List of features extracted from spatial container
        """

        tmp_remove_features = []
        for feature in remove_features:
            if isinstance(feature, int):
                tmp_remove_features.append(feature)
            else:
                tmp_remove_features.append(feature._sc_id)

        remove_features = tmp_remove_features

        # Extract the features by bounds if requested
        if bounds is not None:
            # Extract features by bounds
            keys = (key for key in self._r_tree.intersection(bounds) if key not in remove_features)
            features = [self._features[key] for key in keys if key in self._features]
        else:
            features = [feature for feature in self._features.values() if feature not in remove_features]

        # The following code allows to delete feature while iterating over a get_features request
        for feature in features:
            if feature._sc_id is not None:
                yield feature
            else:
                # If the feature has been deleted do not return it
                pass

        return


class SpatialContainerSTRtree(object):
    """This class manages the spatial features and a spatial index for the STRtree.

    The STRtree is using the STRtree of shapely and has the following limitations
    compared to RTree:
       - immutable rtree
       - only able load features one time (all at the same time)
       - no edit after
       - edition is needed an exception is thrown

    """

    # Class variable that contains the Spatial Container Internal ID
    _sc_id = 0

    def __init__(self):
        """Create an object of type SpatialContainer

        The init will create one container for the feature a dictionary and one
        container for the spatial index (Rtree)

        *Parameters*: None

        *Returns*: *None*

        """

        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def adjust_bbox(self, bounds, delta=GenUtil.ZERO):
        """Adjust the bounding box by increasing by a very small delta

        Parameters:
            bounds: Tuple forming the bounding box (xmin, ymin, wmax, ymax)

        return value:
            altered bounding box (xmin, ymin, wmax, ymax)"""

        xmin, ymin, xmax, ymax = bounds

        xmin -= delta
        ymin -= delta
        xmax += delta
        ymax += delta

        return (xmin, ymin, xmax, ymax)

    def add_features(self, features):
        """Adds all the features in the container and update the spatial index with the feature's bound


        *Parameters*:
            - feature: A spatial feature derives from PointSc, LineStringSc

        *Returns*: *None*

        """

        tmp_features = []
        for feature in features:
            # Check if the type is valid
            if issubclass(feature.__class__, (PointSc, LineStringSc, PolygonSc)):
                pass
            else:
                raise GeoSimException('Unsupported class: {}'.format(str(feature.__class__)))

            bounds = feature.bounds

            # Adjust the bounding box
            bounds = self.adjust_bbox(bounds)

            # Container unique internal counter
            SpatialContainer._sc_id += 1

            # Add the spatial id to the feature
            feature._sc_id = SpatialContainer._sc_id
            feature._sc_scontainer = self

            # Add the feature in the feature container
            self._features[feature._sc_id] = feature

            # Add the bounding box in the bbox_container
            self._bbox_features[feature._sc_id] = bounds

            # Transform the feature as its corresponding bounding box... to simulate the Rtree class
            xmin, ymin, xmax, ymax = bounds
            tmp_feature = LineString(((xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)))
            tmp_feature._sc_id = feature._sc_id
            tmp_features.append(tmp_feature)

        # Load all the features at the same time in the shapely rtree
        self._r_tree = STRtree(tmp_features)

        return

    def add_feature(self, features):
        """Adds a list of feature in the spatial container and

        *Parameters*:
            - feature: A spatial feature derived from Point, LineString or Polygon

        *Returns*: *None*

        """

        raise GeoSimException("Cannot add feature with shapely STRtree")

        return

    def del_feature(self, feature):
        """Delete the feature in the spatial container and in the RTree.

        Raise exception cannot delete feature in shapely RTree

        *Parameters*:
            - feature: The feature to delete in the spatial container

        *Returns*:
            None

        """

        raise GeoSimException("Cannot delete feature with shapely STRtree")

    def del_features(self, features):
        """Delete a list of features in the spatial container

        Raise exception cannot delete feature in shapely RTree

        *Parameters*:
            - features: list of features to delete

        *Returns*:
            Exception

        Exception InternalError: If the key is not in one of the structure

        """

        raise GeoSimException("Cannot delete features with shapely STRtree")

        return ret_value

    def update_spatial_index(self, feature):
        """Update the bounds of the feature in the spatial index

        It will only modify the Rtree spatial index if the bounding
        box of the feature is changed in comparison with the old one.

        *Parameters*:
            - feature: Feature containing the bounds to update

        *Returns*: *None*

        """

        old_bbox = self._bbox_features[feature._sc_id]
        new_bbox = feature.bounds
        old_x_min, old_y_min, old_x_max, old_y_max = old_bbox[0], old_bbox[1], old_bbox[2], old_bbox[3]
        new_x_min, new_y_min, new_x_max, new_y_max = new_bbox[0], new_bbox[1], new_bbox[2], new_bbox[3]

        if old_x_min <= new_x_min and \
                old_y_min <= new_y_min and \
                old_x_max >= new_x_max and \
                old_y_max >= new_y_max:
            # Nothing to do new bounding box is completely included into the old one
            pass
        else:
            # The bounding box has changed...
            raise GeoSimException("Cannot change the bounding box of a feature with shapely STRtree")

        return

    def get_features(self, bounds=None, remove_features=[]):
        """Extract the features from the spatial container.

        According to the parameters the extraction can manage the extraction based on a bounding box using
        the spatial index RTree, some filters to manage extraction based on properties and the possibility
        to remove specific features based on a list of keys

        *Parameters*:
            - bounds: Bounding for the spatial extraction. *None* means all the features
            - remove_keys: List of keys to be removed from the selection

        *Returns*:
            - List of features extracted from spatial container

        """

        tmp_remove_features = []
        for feature in remove_features:
            if isinstance(feature, int):
                tmp_remove_features.append(feature)
            else:
                tmp_remove_features.append(feature._sc_id)

        remove_features = tmp_remove_features

        # Extract the features by bounds if requested
        if (bounds != None):
            # Extract features by bounds
            keys = (key for key in self._r_tree.intersection(bounds) if key not in remove_features)
            features = [self._features[key] for key in keys if key in self._features]
        else:
            features = [feature for feature in self._features.values() if feature not in remove_features]

        # The following code allows to delete feature while iterating over a get_features request
        for feature in features:
            if feature._sc_id is not None:
                yield feature
            else:
                # If the feature has been deleted do not return it
                pass

        return


class ChordalAxis(object):

    # Define the type of Triangle
    ISOLATED = 0  # No side touch another triangl
    TERMINAL = 1  # One side touche another triangle
    SLEEVE = 2  # Two sides touch a triangle
    SLEEVE_X = 3  # Sleeve triangle part of a X junction
    JUNCTION = 4  # All sides touch a triangle
    JUNCTION_T = 5  # T Junction triangle where there is a straight between two of its branches
    JUNCTION_X_FIRST = 6  # First (primary) X Junction triangle
    JUNCTION_X_LAST = 7  # Last (secondary) X Junction triangle
    JUNCTION_X_LENGTH = .2  #  Distance between two X junction

    # Define the type of action(type_action) for creating the centre line
    #   NONE = 0  # No special action is needed when creating the center line
    #   T_EDIT_CENTRE_LINE = 1  # A correction is needed for a T Junction (only used when correcting the centre line)
    #   X_EDIT_CENTRE_LINE = 2  # A correction is needed for a X Junction (only used when correcting the centre line)
    #   NO_CENTRE_LINE = 3  # No centre line is needed for this triangle (only used when correcting the centre line)

    ANGLE_JUNCTION_T = 45.  # Delta used to test if 2 branches or contiguous
    SEARCH_TOLERANCE = None

    def __init__(self, lst_triangle, search_tolerance=GenUtil.ZERO):
        """Constructor

        Parameters
        ----------
        lst_triangle : list
            List of LineString triangle to process
        search_tolerance : float, optional
            Value used for estimating zero

        Return
        ------
        None
        """

        ChordalAxis.SEARCH_TOLERANCE = search_tolerance

        self._validate_triangles(lst_triangle)

        # Transform the Triangle LineString into _TriangleSc to be loaded in SpatialContainer
        for i, triangle in enumerate(lst_triangle):
            triangle = _TriangleSc(triangle.coords)
            lst_triangle[i] = triangle

        # Create spatial container
        self.s_container = SpatialContainer()

        # Load triangles
        self.s_container.add_features(lst_triangle)

        # Load some class variables
        _TriangleSc.s_container = self.s_container

        # Build the cluster (group of triangles part of a polygon)
        self.triangle_clusters = self._build_clusters()

        self.nbr_polygons = len(self.triangle_clusters)
        self.nbr_triangles = len(lst_triangle)

        # Initialise stats value
        self.nbr_lines_pruned = 0
        self.nbr_iteration = 0
        self.nbr_t_junction = 0
        self.nbr_x_junction = 0

        return

    def _validate_triangles(self, lst_triangle):
        """ Validate the each triangle

        Parameters
        ----------
        lst_triangle : lst
            List of LineString or Polygon triangles

        Return
        ------
        None
        """

        if len(lst_triangle) >= 1:
            triangle_type = lst_triangle[0].geom_type
        triangle_valid = True
        for i, triangle in enumerate(lst_triangle):

            # Triangle are LineString
            if triangle.geom_type == GenUtil.LINE_STRING:
                coords = list(triangle.coords)
                # Check number of vertice
                if len(coords) != 4:
                    print("Triangle does not contain exactly 4 vertices: {0}".format(coords[0]))
                    triangle_valid = False
                # Check if all geometry are identical
                if triangle.geom_type != triangle_type:
                    print("Triangle has mixed geometry type: {0}".format(coords[0]))
                    triangle_valid = False
                # Check if the triangle is closed
                if Point(coords[0]).distance(Point(coords[3])) >= ChordalAxis.SEARCH_TOLERANCE:
                    print("Triangle is not closed: {0}".format(coords[0]))
                    triangle_valid = False

            # Triangle are polygons
            if triangle.geom_type == GenUtil.POLYGON:
                coords = list(triangle.exterior.coords)
                # Validate triangle has 4 vertice
                if len(coords) != 4:
                    print("Triangle does not contain exactly 4 vertices: {0}".format(coords[0]))
                    triangle_valid = False
                    # Check if all geometry are identical
                if triangle.geom_type != triangle_type:
                    print("Triangle has mixed geometry type: {0}".format(coords[0]))
                    triangle_valid = False

                # Transform the polygon into a LineString
                if triangle_valid:
                    lst_triangle[i] = LineString(coords)

            # Raise exception if there is an error
            if not triangle_valid:
                # There are one or more errors in the triangles
                raise GeoSimException("Error in the triangles... cannot process them...")

        return

    def _build_clusters(self):
        """Build the clusters of triangle

        One cluster of triangles are all the triangles that have a common edge. One cluster of polygon is equivalent
        to the area of one polygon (including the holes). The set of clusters represent all the polygons

        Parameters
        ----------
        None

        Returns
        -------
        list
           List of clusters formint the different polygon
        """

        clusters = []
        dict_triangles = {}
        for triangle in self.s_container.get_features():
            dict_triangles[triangle.id] = triangle

        # Loop unitl all the triangles are processed
        while len(dict_triangles) >= 1:
            # Create one cluster (one polygon)
            seed_triangle = next(iter(dict_triangles.values()))
            cluster = self._build_one_cluster(dict_triangles, seed_triangle)

            cluster = {triangle.id: triangle for triangle in cluster}
            clusters.append(cluster)

        return clusters

    def _build_one_cluster(self, dict_triangles, seed_triangle):
        """Identify all the triangle that shared an edge (triangles part of one polygon)

        This method is simulating a recursive function using a stack

        Parameters
        ----------
        dict_triangle : dict
            All the triangles to process
        seed_triangle : LineStringSc
            Randomly chosen triangle from which we find all the adjacent triangles

        Return
        ------
        list
            a list of all the raingle forming a polygon
        """

        # Create cluster list to accumulate triangle in cluster
        cluster = []

        # Create the stack to simulate recursivity
        stack = []

        # Initialize the stack with the seed triangle
        stack.append(seed_triangle)

        # Loop over the stack until no more triangle to process
        while stack:
            # Fetch the next triangle
            triangle = stack.pop()
            # Add the triangle in the cluster
            cluster.append(triangle)

            if triangle.id in dict_triangles:
                # Remove from the triangle from the dictionary
                del dict_triangles[triangle.id]
                # Process the adjacent sides
                for adjacent_side_ref in (triangle.adjacent_sides_ref):
                    if adjacent_side_ref is not None:
                        if adjacent_side_ref.id in dict_triangles:
                            stack.append(adjacent_side_ref)  # Simulate recursivity
                        else:
                            # triangle already processed
                            pass
                    else:
                        # No triangle to process
                        pass
            else:
                # Triangle alrerady processed
                pass

        return cluster

    def get_skeleton(self):
        """extract the ceneter line of each triangle merged them and create a list of LineString feature

                This method is simulating recursivity using a stack

        *Parameters*:
            - None

        *Returns*:
            - None

        """
        """Extract the centre line of each triangle; merged them and create a list of LineString features

        Parameters
        ----------
        None

        Return
        ------
        list
            A list of LineString representing the centre line
        """

        merged_centre_lines = []

        # Process each cluster (polygon) independently
        for triangle_cluster in self.triangle_clusters:
            centre_lines = []
            # Process each triangle of one cluster
            for triangle in triangle_cluster.values():
                centre_lines += triangle.centre_line

            merge_centre_line = linemerge(centre_lines)
            merged_centre_line = GenUtil.make_iterable(merge_centre_line)
            merged_centre_lines += merged_centre_line

        return merged_centre_lines

    def correct_skeleton(self):
        """ Apply correction to the skeleton in order to remove unwanted artifact.

        This method is applying 3 corrections on the skeleton
          - Firstly, the skeleton is pruned; the small skeleton segment from a junction triangleare removed based
            are removed (deleted); we are not using an absolute tolerance for pruning but a tolerance relative to
            the size of the junction triangle
          - Secondly, T junction are corrected and the 2 branch that are forminf a natural continuity are joind
          - Thirdly, X junction are reconstructed bu joining together two adjacent or near adjacent T junctions

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        # Prune the small branch from the Junction triangle until there are no more small branch to prune (iterative process)
        nbr_iteration = 0
        for triangle_cluster in self.triangle_clusters:  # Loop over each polygon (cluster)
            while True:
                nbr_pruned = 0
                nbr_iteration += 1
                for triangle in list(triangle_cluster.values()):  # Loop over each triangle of one cluster (polygon)
                    if triangle.type == ChordalAxis.JUNCTION:
                        junction_pruned = self.prune_junction(triangle_cluster, triangle)
                        if junction_pruned >= 1:
                            nbr_pruned += junction_pruned
                self.nbr_lines_pruned += nbr_pruned
                if nbr_pruned == 0:
                    self.nbr_iteration = max(self.nbr_iteration, nbr_iteration)
                    break  # Nothing more to prune in this cluster... exit

        # Correct the T junction to form a straight line if 2 brancj have the same orientation
        for triangle_cluster in self.triangle_clusters:  # Loop over each polygon (cluster)
            for triangle in list(triangle_cluster.values()):  # Loop over each triangle of one polygon
                if triangle.type == ChordalAxis.JUNCTION:
                    sides_t_junction = self.adjust_t_junction(triangle)
                    if sides_t_junction is not None:
                        self.nbr_t_junction += 1
                        triangle.junction_side_a = sides_t_junction[0]
                        triangle.junction_side_b = sides_t_junction[1]
                        triangle.reset_attributes()

        # Correct the X junction to join adjacent junction and remove the line joining the 2 junctions
        total_x_junctions_infos = []
        for triangle_cluster in self.triangle_clusters:  # Loop over each polygon (cluster)
            for triangle in list(triangle_cluster.values()):  # Loop over each triangle of one polygon
                if triangle.type in (ChordalAxis.JUNCTION, ChordalAxis.JUNCTION_T):
                    junction_x_infos = self.adjust_x_junction(triangle_cluster, triangle)
                    if len(junction_x_infos) >= 1:
                        total_x_junctions_infos += [junction_x_infos]
                    else:
                        # Not a valid X junction. Nothing to do
                        pass

        # Remove from the X junction correction, the junction that are indirectly linked
        # If junction A and B are to be merged and B and C are to be merged this means
        # that A, B, C are closed.  Correcting three or more X junction creates bad artifacts
        # so do not correct junction A, B, C
        id_to_remove = []
        for junction_x_infos in total_x_junctions_infos:
            if len(junction_x_infos) >= 2:  # This junction is closed to two of more possible X junction
                for x_infos in junction_x_infos:
                    id_to_remove += [x_infos.first_junction.id, x_infos.last_junction.id]

        for junction_x_infos in total_x_junctions_infos:
            if len(junction_x_infos) == 1:
                junction_x_infos = junction_x_infos[0]
                first_junction = junction_x_infos.first_junction
                last_junction = junction_x_infos.last_junction
                if first_junction.id not in id_to_remove and last_junction.id not in id_to_remove:
                    # Valid X junction to process
                    possible_junction = [ChordalAxis.JUNCTION, ChordalAxis.JUNCTION_T]
                    if first_junction.type in possible_junction and last_junction.type in possible_junction:
                        self.nbr_x_junction += 1
                        # Set some attribute in the triangle
                        first_junction.type = ChordalAxis.JUNCTION_X_FIRST
                        first_junction.junction_x_mid_pnt_sides = junction_x_infos.mid_pnt_sides
                        first_junction.junction_x_centroid = junction_x_infos.x_centroid
                        last_junction.type = ChordalAxis.JUNCTION_X_LAST
                        for sleeve_triangle in junction_x_infos.sleeve_in_branch:
                            # Track the sleeve triangle between the 2 junctions
                            sleeve_triangle.type = ChordalAxis.SLEEVE_X
                    else:
                        # Junction has already been processed... nothing to do
                        pass

#       Important note: After the X junction correction the complex data structure is not maintained anymore...
#                       To keep in mind if more correction has to be done followinf this...

        return

    def adjust_t_junction(self, junction_triangle):
        """ Apply the correction for the T junction triangle

        A junction is a valid X junction when 2 T junction of near each other and should form instead a T junction

        Parameters
        ----------
        junction_triangle : Triangle
            The junction triangle to process

        Return
        ------
        list
            A list of 2 numbers between [0..2] that determines the 2 sides forming the T junction.  None if the junction is not a T junction
        """

        sides_t_junction = None
        for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
            # Special case when a junction triangle is adjacent to another junction triangle.
            # The junction triangl is automaticaly a T junction
            if junction_triangle.adjacent_sides_ref[i].type in (ChordalAxis.JUNCTION, ChordalAxis.JUNCTION_T) and \
               junction_triangle.adjacent_sides_ref[j].type == ChordalAxis.SLEEVE and \
               junction_triangle.adjacent_sides_ref[k].type == ChordalAxis.SLEEVE:
                sides_t_junction = [j, k]
                break

        if sides_t_junction is None:
            branches = []
            # Loop each branch of the junction triangle
            for next_triangle in junction_triangle.adjacent_sides_ref:
                if next_triangle.type == ChordalAxis.SLEEVE:
                    branch = Branch(junction_triangle, next_triangle)
                    branches.append(branch)

            # Evaluate if 2 branches is forming an almost straight line
            branch_angle = [branch.angle for branch in branches]  # Extract the angles of the branch
            if len(branches) == 3:
                angle_max = ChordalAxis.ANGLE_JUNCTION_T
                for i, j in [(0, 1), (1, 2), (2, 0)]:
                    delta_angle = abs(180. - abs(branch_angle[i] - branch_angle[j]))
                    if delta_angle < angle_max:
                        angle_max = delta_angle
                        sides_t_junction = [i, j]
            else:
                # No T junction to process
                pass

        return sides_t_junction

    def adjust_x_junction(self, cluster_triangle, current_junction):
        """ Apply the correction for the X junction triangle

        A junction is a valid T junction when the 3 branches of a triangle is composed of sleeve or terminal triangle;  and each branch is of a certain length
        and 2 of the branches form an almost straight line.

        Parameters
        ----------
        junction_triangle : Triangle
            The junction triangle to process

        Return
        ------
        object
            An object containing the information needed to create the X junction
        """

        junction_x_infos = []
        #  Loop over each branch of the junction triangle
        for adjacent_junction in current_junction.adjacent_sides_ref:
            branch = Branch(current_junction, adjacent_junction)
            last_triangle = branch.triangle_in_branch[-1]  # Extract the last triangle of the branch

            # If the last triangle in a branch is junctin and within a certain tolerance it's a candidate for X Junction
            if last_triangle.type in (ChordalAxis.JUNCTION, ChordalAxis.JUNCTION_T) and \
                    branch.length < min(current_junction.width, last_triangle.width) * ChordalAxis.JUNCTION_X_LENGTH:
                # Merge the triangle lin the branch to form only one polygon
                line_triangles = [current_junction] + branch.triangle_in_branch
                pol_triangles = [Polygon(line.coords) for line in line_triangles]
                merged_pol = unary_union(pol_triangles)
                if merged_pol.geom_type == GenUtil.POLYGON:  # Merged the triangle to form only one polygon
                    x_centroid = merged_pol.centroid
                    merged_line = LineString(merged_pol.exterior.coords)

                    # Detect which mid side point we must keep (we must keep four only)
                    mid_pnt_sides = current_junction.mid_pnt_sides + last_triangle.mid_pnt_sides
                    new_mid_pnt_sides = []
                    for mid_pnt_side in mid_pnt_sides:
                        if mid_pnt_side.distance(merged_line) < ChordalAxis.SEARCH_TOLERANCE:
                            new_mid_pnt_sides.append(mid_pnt_side)

                    # Validate the the center line
                    if self.validate_x_junction(merged_pol, x_centroid, new_mid_pnt_sides):
                        junction_x_info = Holder(first_junction=current_junction, last_junction=branch.triangle_in_branch[-1],
                                                 sleeve_in_branch=branch.triangle_in_branch[0:-1], mid_pnt_sides=new_mid_pnt_sides, x_centroid=x_centroid)
                        junction_x_infos.append(junction_x_info)
                else:
                    # Invalid merging nothing to do (should not happen)
                    pass
            else:
                # Not a triangle to process for X junction
                pass

        return junction_x_infos

    def validate_x_junction(self, merged_pol, centroid, new_mid_pnt_sides):
        """ Validate the X junction centre line

        The center line of a x junction is valid when each new center line are located inside the new polygon.
        In some special the centre line can be outside and the junction is not a candidate for a X junction

        Parameters
        ----------
        merged_pol : Polygon
            the polugon formed by the 2 junction and the possible sleeves polygon

        Return
        ------
        bool
            True: valid X junction; False: Invalid X junction
        """

        buf_merged_pol = merged_pol.buffer(.01)

        status = True
        for mid_pnt_side in new_mid_pnt_sides:
            line = LineString((centroid, mid_pnt_side))
            if line.crosses(buf_merged_pol):
                # Centre line not valid
                status = False
                break

        # Validate the angle between the line
        for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
            x0,y0 = new_mid_pnt_sides[i].x - centroid.x, new_mid_pnt_sides[i].y - centroid.y
            x1,y1 = new_mid_pnt_sides[j].x - centroid.x, new_mid_pnt_sides[j].y - centroid.y
            angle = GenUtil.difference_angle_vector((x0, y0), (x1, y1), ChordalAxis.SEARCH_TOLERANCE)

        return status

    def prune_junction(self, cluster_triangle, junction_triangle):
        """This function prune a junction triangle of the branches that are below a certain tolerance

        This function is looping over the 3 branches of the junction triangle.  If one or more branches
        are below the tolerance than the triangle part of the branches are deleted in order to prune
        the sjkeleton from unwanted small lines.

        Parameters
        ----------
        cluster_triangle : List
            List of the Triangle forming part of one cluster (polygon)
        junction_triangle : Triangle
            Junction triangle to process

        Return
        ------
        int
            Number of branches deleted
        """

        branches = []

        for next_triangle in junction_triangle.adjacent_sides_ref:
            branch = Branch(junction_triangle, next_triangle)
            # Only branches below tolrance and finishing with a Terminal triangle
            if branch.last_triangle_type == ChordalAxis.TERMINAL and branch.length <= junction_triangle.width:
                branches.append(branch)

        if len(branches) == 3:
            # If the three branches of the junction Triangle are below the tolerance
            # Remove the  branch with the smallest tolerance
            max_length = sys.float_info.max
            for branch in branches:
                if branch.length < max_length:
                    del_branches = [branch]
                    max_length = branch.length

        elif len(branches) == 2:
            # Two branches are below the tolerance
            if branches[0].length < branches[1].length:
                branch_0 = branches[0]
                branch_1 = branches[1]
            else:
                branch_0 = branches[1]
                branch_1 = branches[0]
            if branch_0.length < .3 * branch_1.length:
                # If one branch is very smaller than the other one delete it
                del_branches = [branch_0]
            else:
                # Otherwise delete both
                del_branches = [branch_0, branch_1]
        elif len(branches) == 1:
            # Only one branch is smaller than the tolerance ==> delete it
            del_branches = [branches[0]]
        else:
            del_branches = []

        if len(del_branches) >= 1:
            #  Process the triangle to delete
            triangles_to_reset = []
            triangles_to_isolate = []
            for branch in del_branches:
                # Accumulate each triangle of each branch
                for triangle in branch.triangle_in_branch:
                    # Special case for the triangle that are referenced (adjacent) by the junction triangle
                    for ref_triangle in triangle.adjacent_sides_ref:
                        if ref_triangle is not None:
                            triangles_to_reset.append(ref_triangle)
                    triangles_to_isolate.append(triangle)

            # Reset all the attributes of the referenced triangles
            for triangle in triangles_to_reset:
                triangle.reset_attributes()

            # Delete the triangle in the branch
            self.s_container.del_features(triangles_to_isolate)

            #  Delete the triangles in the cluster
            for triangle in triangles_to_isolate:
                del cluster_triangle[triangle.id]

        return len(del_branches)


class _TriangleSc(LineStringSc):
    """LineString specialization to be included in the SpatialContainer"""

    id = 0

    def __init__(self, coords):
        super().__init__(coords)

        # Add unique id to each Triangle
        self.id = _TriangleSc.id

        _TriangleSc.id += 1

    @property
    def mid_pnt_sides(self):
        """Attribute function to extract the mid point of each side of the triangle

        This attribute function is using on the fly calculation to extract theattribute

        A junction is a valid T junction when the 3 branches of a triangle is composed of sleeve or terminal triangle;  and each branch is of a certain length
        and 2 of the branches form an almost straight line.

        Parameters
        ----------
        None

        Return
        ------
        list
            List containing the coordinate (x,y) the mid point of each side of the triangle
        """

        try:
            # Implement attribute late calculation
            return self._mid_pnt_sides
        except AttributeError:
            if not hasattr(self, 'junction_x_mid_pnt_sides'):
                # Calculate the mid point of each side of the triangle
                coords = list(self.coords)
                mid_pnt_side_0 = LineString([coords[0], coords[1]]).interpolate(0.5, normalized=True)
                mid_pnt_side_1 = LineString((coords[1], coords[2])).interpolate(0.5, normalized=True)
                mid_pnt_side_2 = LineString((coords[2], coords[0])).interpolate(0.5, normalized=True)
                self._mid_pnt_sides = [mid_pnt_side_0, mid_pnt_side_1, mid_pnt_side_2]
            else:
                # The figure is not anymore a triangle
                self._mid_pnt_sides = self.junction_x_mid_pnt_sides
            return self._mid_pnt_sides

    @property
    def type(self):
        """Attribute function to extract the type of triangle

        This attribute function is using on the fly calculation to extract the attribute

        Parameters
        ----------
        None

        Return
        ------
        int
            The type of the triangle
        """

        try:
            # Attribute late calculation
            return self._type
        except AttributeError:
            nbr_side = 0
            for adjacent_side_ref in self.adjacent_sides_ref:
                if adjacent_side_ref != None:
                    nbr_side += 1

            if nbr_side == 0:
                self._type = ChordalAxis.ISOLATED
            elif nbr_side == 1:
                self._type = ChordalAxis.TERMINAL
            elif nbr_side == 2:
                self._type = ChordalAxis.SLEEVE
            elif nbr_side == 3:
                self._type = ChordalAxis.JUNCTION
                if hasattr(self, 'junction_side_a'):
                    self._type = ChordalAxis.JUNCTION_T  # Junction form a T Junction
            elif nbr_side >= 4:
                raise GeoSimException ('Internal error...!!!')

            return self._type

    @type.setter
    def type(self, value):
        """Attribute function to set the type of triangle

        There is a special case where for X junction we delete the centre line

        Parameters
        ----------
        Value : int
            The type of the triangle

        Return
        ------
        None
        """

        self._type = value

        if self._type in ( ChordalAxis.JUNCTION_X_FIRST, ChordalAxis.JUNCTION_X_LAST, ChordalAxis.SLEEVE_X):
            try:
                del self._centre_line
            except AttributeError:
                pass

    @property
    def width(self):
        """Attribute function to extract the width of a junction triangle

        This attribute function is using on the fly calculation to extract the attribute.
        the width of the triangle is 2 times the length of the longest centre line of a junction triangle

        Parameters
        ----------
        None

        Return
        ------
        real
            The width of the triangle
        """

        try:
            # On the fly calculation
            return self._width
        except AttributeError:
            lst_length = [line.length for line in self.centre_line]
            max_length = max(lst_length)
            self._width = max_length * 2.

            return self._width

    @property
    def adjacent_sides_ref(self):
        """Attribute function to extract the adjacent triangle

        This attribute function is using on the fly calculation to extract the attribute

        Parameters
        ----------
        None

        Return
        ------
        list
            List of 3 containing a reference to the adjacent triangle
        """

        try:
            # On the fly calculation
            return self._adjacent_sides_ref
        except AttributeError:

            self._adjacent_sides_ref = []
            # Loop on each side (each mid_pnt_side) to find adjacent triangles
            for mid_pnt_side in self.mid_pnt_sides:

                # Find all potential adjacent triangles
                potential_triangles = _TriangleSc.s_container.get_features(bounds=mid_pnt_side.bounds, remove_features=[self])

                # Find the closest triangle
                triangles = [(triangle, mid_pnt_side.distance(triangle)) for triangle in potential_triangles]
                sorted(triangles, key=lambda triangle: triangle[1])  # Sort by distance
                triangles = [triangle for (triangle, distance) in triangles if distance < ChordalAxis.SEARCH_TOLERANCE]
                if len(triangles) == 0:
                    self._adjacent_sides_ref.append(None)  # No  triangle
                if len(triangles) == 1:
                    self._adjacent_sides_ref.append(triangles[0])
                elif len(triangles) >= 1:
                    xy = (mid_pnt_side.x, mid_pnt_side.y)
                    print("***Warning*** Triangles may be to small: {0} Try to simplify them (e.g. Douglas Peucker)".format(xy))
                    self._adjacent_sides_ref.append(triangles[0])

            return self._adjacent_sides_ref

    @property
    def centre_line(self):
        """Attribute function to extract the centre line of a triangle

        This attribute function is using on the fly calculation to extract the attribute.
        The centre line depends on the type of triangle each type of triangle having a different representation of centre line

        Parameters
        ----------
        None

        Return
        ------
        List
            Zero to 4 LineString defining the centre line of the triangle
        """
        try:
            # On the fly attribute calculation
            return self._centre_line
        except AttributeError:

            centre_line = []
            coords = list(self.coords)

            # Process each case depending on the number of internal side of the triangle
            if self.type == ChordalAxis.ISOLATED:
                # Degenerated polygon with one triangle no skeleton line added
                pass

            elif self._type == ChordalAxis.TERMINAL:
                # Terminal triangle add line from the extremity of the triangle up to mid opposite side
                if self.adjacent_sides_ref[0] != None:
                    coords_line = [coords[2], self.mid_pnt_sides[0]]
                if self.adjacent_sides_ref[1] != None:
                    coords_line = [coords[0], self.mid_pnt_sides[1]]
                if self.adjacent_sides_ref[2] != None:
                    coords_line = [coords[1], self.mid_pnt_sides[2]]

                centre_line.append(LineString(coords_line))

            elif self.type == ChordalAxis.SLEEVE:
                # Sleeve triangle skeleton added between the mid point of side adjacent to another triangle
                mid_pnt = []
                for i, adjacent_side_ref in enumerate(self._adjacent_sides_ref):
                    if adjacent_side_ref != None:
                        mid_pnt.append(self.mid_pnt_sides[i])
                centre_line.append(LineString([mid_pnt[0], mid_pnt[1]]))

            elif self.type == ChordalAxis.SLEEVE_X:
                # No center line to create
                pass

            elif self.type == ChordalAxis.JUNCTION:
                # Regular triangle T type. Centroid is the centroid of the triangle
                pnt_x = (coords[0][0] + coords[1][0] + coords[2][0]) / 3.
                pnt_y = (coords[0][1] + coords[1][1] + coords[2][1]) / 3.
                centroid = [pnt_x, pnt_y]
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    centre_line.append(LineString([centroid, mid_side_pnt]))

            elif self.type == ChordalAxis.JUNCTION_T:
                # Corrected triangle T. Centroid is the middle point between the 2 aligned branches
                pnt0 = self.mid_pnt_sides[self.junction_side_a]
                pnt1 = self.mid_pnt_sides[self.junction_side_b]
                pnt = LineString([(pnt0.x, pnt0.y), (pnt1.x, pnt1.y)]).interpolate(0.5, normalized=True)
                centroid = [pnt.x, pnt.y]
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    centre_line.append(LineString([centroid, mid_side_pnt]))

            elif self.type == ChordalAxis.JUNCTION_X_FIRST:
                centroid = (self.junction_x_centroid.x, self.junction_x_centroid.y)
                #  create the center line
                for mid_pnt_side in self.junction_x_mid_pnt_sides:
                    centre_line.append(LineString([centroid, mid_pnt_side]))

            elif self.type == ChordalAxis.JUNCTION_X_LAST:
                # No center line to create
                pass

            else:
                raise GeoSimException ("Unknow triangle type: {1}".format(self.type))

            return centre_line

    def reset_attributes(self):
        """Function used to reset all the attributes of a triangle

        Because all the attributes are using on the fly calculation the attributes will be recalculated as needed

        Parameters
        ----------
        None

        Return
        ------
        None
        """

        if hasattr(self, '_adjacent_sides_ref'): del self._adjacent_sides_ref
        if hasattr(self, '_width'): del self._width
        if hasattr(self, '_centre_line'): del self._centre_line
        if hasattr(self, '_mid_pnt_sides'): del self._mid_pnt_sides
        if hasattr(self, '_type'): del self._type


class Branch:
    """Class used to build and manage a branch of a junction triangle"""

    def __init__(self, current_triangle, next_triangle):
        """Constructor of a branch of a junction triangle

        Because all the attributes are using on the fly calculation the attributes will be recalculated as needed

        Parameters
        ----------
        current_triangle : Triangle
            The starting junction triangle to calculate the branch
        next_triangle : Triangle
            The starting point (side) to evaluate the branch

        Return
        ------
        None
        """

        self.current_triangle = current_triangle
        self.triangle_in_branch = []
        self.length = 0.
        max_length = current_triangle.width * 3.  # Maximum length of a branch

        while True:
            self.triangle_in_branch.append(next_triangle)
            if next_triangle.type in (ChordalAxis.SLEEVE, ChordalAxis.TERMINAL):
                self.length += next_triangle.centre_line[0].length  # Add the length
                if next_triangle.type == ChordalAxis.TERMINAL:
                    # Terminal triangle is the end of the branch
                    break
            else:
                # Junction triangle is the end of the branch
                break

            # Loop for the next adjacent triangle
            if self.length < max_length:
                adjacents = [adjacent for adjacent in next_triangle.adjacent_sides_ref if adjacent is not None]
                if adjacents[0].id == current_triangle.id:
                    current_triangle, next_triangle = next_triangle, adjacents[1]
                else:
                    current_triangle, next_triangle = next_triangle, adjacents[0]
            else:
                # End of looping reached
                break

        # Extract the type of the last triangle in the branch
        self.last_triangle_type = self.triangle_in_branch[-1].type

        return

    @property
    def angle(self):
        """Attribute function to extract the angle of the branch

        The angle of the branch is calculated using the angle between the first and the last coordinate
        of the centre line of the branch

        Parameters
        ----------
        None

        Return
        ------
        float
            The angle between [¸0..360]
        """

        try:
            return self._angle
        except AttributeError:
            # Extract the skeleton line from the branch
            lines = []
            for triangle in self.triangle_in_branch:
                if triangle.type in [ChordalAxis.SLEEVE, ChordalAxis.TERMINAL]:
                    lines += triangle.centre_line

            # Merge the lines to form one line
            merged_line = linemerge(lines)
            if merged_line.geom_type == GenUtil.LINE_STRING:
                # Extract the angle formed by the first and last coordinate
                x0, y0 = merged_line.coords[0][0], merged_line.coords[0][-1]
                x1, y1 = merged_line.coords[-1][0], merged_line.coords[-1][-1]
            else:
                # There was an error in line merge so just take the first line (less good but should not happen often)
                line  = lines[0]
                x0, y0 = line.coords[0][0], line.coords[0][-1]
                x1, y1 = line.coords[-1][0], line.coords[-1][-1]

            # Checked that the first coordinate of the line is located on the triangle
            if self.current_triangle.distance(Point(x0, y0)) < ChordalAxis.SEARCH_TOLERANCE:
                # The line is well oriented
                pass
            else:
                # Invert the orientation of the line
                x0, y0, x1, y1 = x1, y1, x0, y0

            self._angle = GenUtil.difference_angle_vector((x0,y0), (x1,y1), ChordalAxis.SEARCH_TOLERANCE)

            return self._angle


class Holder(object):
    """Generic class creator"""

    def __init__(self, **kwds):
        """Construction function


        Parameters
        ----------
        **kwds
            Lis of keywords to transform into attribute

        Return
        ------
        None
        """
        self.__dict__.update(kwds)


class GeoSimException(Exception):
    """Generic GeoSimException.  Base class exception for other GoSimException"""

    def __init__(self, *arguments, **keywords):
        """Constructor"""
        Exception.__init__(self, *arguments, **keywords)


class InternalError(GeoSimException):
    """InternalError exception"""

    def __init__(self, *param_names):
        """Constructor"""

        GeoSimException.__init__(self, *param_names)
