#!/usr/bin/env python
# -=- encoding: utf-8 -=-

"""
General classes and utilities needed for the GENeralization MEta ALgorithm (GENMTEAL) tool

"""

import math
from shapely.geometry import Point, LineString, LinearRing, Polygon
from shapely.ops import linemerge
import fiona

try:
    from rtree import Rtree
    lib_Rtree = True
except :
    # If Rtree is not installed
    lib_Rtree = False
    from shapely.strtree import STRtree

class LineStringSc(LineString):

    """LineString specialization to be included in the SpatialContainer"""

    def __init__(self, coords):
        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        return super().coords

    @coords.setter
    def coords(self, coords):
        # Update the coord attribute in the parent class
        LineString.coords.__set__(self, coords)

        if self._sc_scontainer != None:    # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_scontainer.update_spatial_index(self)


class PointSc(Point):

    """LineString specialization to be included in the SpatialContainer"""

    def __init__(self, coords):
        super().__init__(coords)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def coords(self):
        return super().coords

    @coords.setter
    def coords(self, coords):
        Point.coords.__set__(self, coords)

        if self._sc_scontainer != None:  # Is the feature is a spatial container
            # The coordinate has changed so update the bounding box in the spatial container
            self._sc_container.update_bbox(self)


class PolygonSc(Polygon):

    """Polygon specialization to be included in the SpatialContainer"""

    def __init__(self, exterior, interiors=None):
        super().__init__(exterior, interiors)
        self._sc_id = None
        self._sc_scontainer = None

    @property
    def exterior(self):
        return super().exterior

    @property
    def interiors(self):
        return super().interiors

    @exterior.setter
    def exterior(self, exterior):
        raise GeoSimException ("Cannot update the exterior coordinates of a polygon")

    @interiors.setter
    def interiors(self, interiors):
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
    def distance(p1, p2):
        """Calculate the euclidean distance between 2 points

        *Parameters*:
            - p1: (x,y) tuple of the first coordinate
            - p2: (x,y) tuple of the second coordinate

        *Returns*:
            - Distance between the 2 points (real)

        """

        return math.sqrt((p2[0] - p1[0]) ** 2.0 + (p2[1] - p1[1]) ** 2.0)

    @staticmethod
    def compute_angle (p1, p2, p3, type=DEGREE):
        """
        Function to calculate angle between two vectors.
        """

        return GenUtil.angle_vector(p1, p2, p3, type)

    @staticmethod
    def angle_vector(p1, p2, p3, type=DEGREE):
        """Calculate the angle formed by the vector p1-p2 and p2-p3

        *Parameters*:
            - p1: (x,y) tuple of coordinates
            - p2: (x,y) tuple of coordinates
            - p3: (x,y) tuple of coordinates
            - type: Angle type DEGREE or ANGLE

        *Returns*:
            - The angle between the vector p1-p2 and p2-p3 (float)

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

        theta = math.acos(value)

        if type == GenUtil.DEGREE:
            theta = math.degrees(theta)

        return theta

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

        orient = ((p0[0] - p1[0]) * (p2[1] - p1[1])) - ((p2[0] - p1[0]) * (p0[1] - p1[1]))

        if orient > 0.:
            orient = 1
        elif orient< 0.:
            orient = -1
        else:
            orient = 0

        return orient

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

        return (x, y)

    @staticmethod
    def calculate_compactness_index(area, perimeter):
        """Calculate the compactness index based of the perimeter and area

        Args:
            area (float): Area of the polygon
            perimeter (float): Perimeter of the area

        Return:
            (float): Compactness index

        """

        return 4 * area * math.pi / (perimeter ** 2.0)

    @staticmethod
    def build_bounding_box(tolerance, coord):
        """Create and adjust a bounding box (xmin, ymin, xmax, ymax) with a small tolerance

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
    def calculate_adjusted_area(area, cmp_index):
        """Calculate the adjusted area from the area and compactness index

        Args:
            area (float): Area of the polygon
            cmp_index (float): Compactness index of the areea

        Return:
            flot: Adjusted area of the polygon

            """
        return area * (0.75 / cmp_index)

    @staticmethod
    def read_in_file(in_file, geo_content, layer_in=None):
        """
        Read and load the vectors in the input file

        Args:
            in_file (str): Name of the input file (geopackage)
            geo_content (dict): Dictionary containing information to create the spatial database
            layer_in (str): Layer name to read

        Return:
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
        """
        Write the vectors in the output file

        Args:
            out_file (str): Name of the output file (geopackage)
            geo_content (DataClass): Contains information to create the spatial database

        Return:
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
                                if feature.sb_layer_name==layer_name):
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
                        coordinates = [exterior]+interiors
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
        """
        Write the vectors in the output file

        Args:
            out_file (str): Name of the output file (geopackage)
            geo_content (DataClass): Contains information to create the spatial database

        Return:
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
                                if feature.sb_layer_name==layer_name):
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
                        coordinates = [exterior]+interiors
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

        *Parameters*: None

        *Returns*: *None*

        """

        self._r_tree = Rtree()  # Container for the Rtree
        self._features = {}  # Container to hold the features
        self._bbox_features = {}  # Container to hold the bounding boxes

    def adjust_bbox(self, bounds, delta = GenUtil.ZERO):
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

    def add_feature(self, feature):
        """Adds a feature in the container and update the spatial index with the feature's bound

        To be added in the container a spatial feature must be a MA_Point, MA_LineString or
        MA_Polygon.

        *Parameters*:
            - feature: A spatial feature derives from PointSc, LineStringSc

        *Returns*: *None*

        """

        # Check if the type is valid
        if issubclass(feature.__class__, (PointSc, LineStringSc, PolygonSc)):
            pass
        else:
            raise GenException('Unsupported class: {}'.format(str(feature.__class__)))

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

        *Parameters*:
            - feature: A spatial feature derived from Point, LineString or Polygon

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
            # The bounding box has changed
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

    def adjust_bbox(self, bounds, delta = GenUtil.ZERO):
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
            tmp_feature = LineString(((xmin,ymin),(xmin,ymax),(xmax,ymax), (xmax,ymin), (xmin,ymin)))
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

class ChordalAxis1(object):

    def __init__(self, lst_triangles, search_tolerance=GenUtil.ZERO):

        self._load_triangles(lst_triangles)

        self._build_triangle_clusters()


    def _load_triangles(self, lst_triangles):

        lst_triangles_sc = []
        for triangle in lst_triangles:
            triangle_sc = LineStringSc(triangle.coords)
            triangle_sc._sc_processed = False
            lst_triangles_sc.append(triangle_sc)

            


class ChordalAxis(object):
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

    def __init__(self, triangles, search_tolerance=GenUtil.ZERO):
        """Constructor of the class

        Parameters:
            - polygon: PolygonSc to process
            - triangles: List of LineString. The line represent the triangles of a Constriant Delanuay Triangulation (CDT)
                         Each triangle is composed of 4 non colinear vertice
            - minimal_width: Float used to prune the skeleton outputted by the Chordal Axis Transform and to identify
                             bottleneck triangles
            - search_tolerance: Search tolerance can vary depending on the dynamic of the data set from lat-lon to Lambert conformal

        Return value: None
        """
        self.s_cont_triangles = SpatialContainer()
        self._search_tolerance = search_tolerance
#        self._process_polygon(polygon)
        self._load_triangles(triangles)

#        _Triangle.line_segments = self.line_segments
#        _Triangle.perimeters = self.perimeter_distance

        self._build_skeleton()
#        self._prune_skeleton()

    def _process_polygon(self, polygon):
        """Process a polygon to create the object property line_segments and perimeter_distance

        Parameters:
          - PolygonSc to process

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
            - List of LineString.  Each LineString is composed of 4 non colinear vertice forming a triangle

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
        is a bifurcation (skeleton splitting).  The method is removing noisy skeleton arm that are
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
                for i in range(3):
                    p0 = triangle.coords[i]
                    p1 = triangle.coords[(i + 1) % 3]
                    # Check if the perimeter distance is below the minimal width
                    extremity = self.perimeter_distance.is_extremity(p0, p1, self._minimal_width)
                    if extremity:
                        # The side/skeleton must be pruned
                        lst_sub_coords = self.perimeter_distance.get_sub_perimeter(p0, p1)
                        # Creation of a polygon with the vertice included in the perimeter
                        extremity_polygon = PolygonSc(lst_sub_coords)
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

    def get_skeletton(self, noise=0.):
        """Extract the Chordal Axis Transform skeleton from a constrained Delanauy trianulation

        Parameters: None

        Return value:
            - List of LineString of the skeleton of the polygon
        """

        center_lines = []
        for triangle in self.s_cont_triangles.get_features():
            center_lines += triangle.get_centre_line()

        # Merge the center line
        merged_center_lines = linemerge(center_lines)
        if merged_center_lines.geom_type == GenUtil.LINE_STRING:
            center_lines = [merged_center_lines]
        else:
            center_lines = [center_line for center_line in merged_center_lines]

        # Transform LineString into LineStringSc
        tmp_center_lines = []
        for center_line in center_lines:
            tmp_center_lines.append(LineStringSc(center_line.coords))
        center_lines = tmp_center_lines

        if noise > 0.:
            center_lines = self._remove_noise(center_lines, noise)

        return center_lines

    def get_triangles(self):
        """Extract the triangles from the constrained delanauy triangulation

        Parameters: None

        Return value:
            - List of the LineString of 4 vertice each.
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
            tri = triangle.cloner()  # Create a copy of the triangle not a reference
            tri.ma_properties[ChordalAxis.CODE] = category
            tri.ma_properties[ChordalAxis.WIDTH] = width
            tri.ma_properties[ChordalAxis.CENTER_LINE] = triangle.get_centre_line()
            triangles.append(tri)

        return triangles


    def  _remove_noise(self, center_lines, noise):
        """remove the noise (small lines) in the skeleton

        Parameters:
            center_lines: List of LineString

        Return:
            List of LineString
        """

#        a = LineString(((0,0),(15,0)))
#        b = LineString(((15,0),(30,0)))
#        c = LineString(((30,0),(30,5)))
#        d = LineString(((30,0),(45,0)))
#        e = LineString(((45,0),(60,0)))
#        center_lines = [a,b,c,d,e]

        # Load the features in the spatial container (accelerate the search)
        s_container = SpatialContainer()
        for center_line in center_lines:
            center_line.sb_geom_type = GenUtil.LINE_STRING
            s_container.add_feature(center_line)

        # Build topology. For each create list of connecting lines
        for line in center_lines:
            p_start = line.coords[0]
            p_end = line.coords[-1]
            b_box = GenUtil.build_bounding_box(self._search_tolerance, p_start)
            lines_b_box = s_container.get_features(b_box, remove_features=[line])
            line.start_lines = []
            for line_b_box in lines_b_box:
                if line.touches(line_b_box):
                    line.start_lines.append(line)
            b_box = GenUtil.build_bounding_box(self._search_tolerance, p_end)
            lines_b_box = s_container.get_features(b_box, remove_features=[line])
            line.end_lines = []
            for line_b_box in lines_b_box:
                if line.touches(line_b_box):
                    line.end_lines.append(line)

        lines = s_container.get_features()
        center_lines = []
        for line in lines:
            keep_line = True
            if line.length <= noise:
                # Only line below noise length are candidate to be removed
                if len(line.start_lines) == 0 and len(line.end_lines) == 0:
                    # Isolated line . Nothing to do
                    pass
                else:
                    if len(line.start_lines) != 0:
                        tmp_lines = line.start_lines
                    else:
                        tmp_lines = line.end_lines
                    for tmp_line in tmp_lines:
                        if len(tmp_line.start_lines) == 0 or len(tmp_line.end_lines) == 0:
                            keep_line = False
                            break
            if keep_line:
                center_lines.append(line)
            else:
                print ("Line deleted")

        # Merge the center line
        merged_center_lines = linemerge(center_lines)
        if merged_center_lines.geom_type == GenUtil.LINE_STRING:
            center_lines = [merged_center_lines]
        else:
            center_lines = [center_line for center_line in merged_center_lines]


        return center_lines


class _Triangle(LineStringSc):
    """Calculates the
    """

    SUPERIMPOSED = 'superimposed'
    INTERNAL = 'internal'

    line_segments = None
    perimeters = None

    def __init__(self, lst_coords):
        if len(lst_coords) == 4:
            super().__init__(lst_coords)
            self.sb_geom_type = GenUtil.LINE_STRING
            self._nbr_internal = None  #
            self._centre_lines = None
            self._side_type = None
            self._category = None
            self._mid_triangle = None
        else:
            raise Exception("A triangle must have 4 and only 4 coordinates")

    def _get_obtuse_angle(self):
        """Return the vertice number (0, 1, 2) of the obtuse angle of the triangle

        Parameters: None

        Return value:
            - integer: Vertice number of the obtuse angle or None if there is not obtuse angle

        """

        acute_angle = None

        for i in range(3):
            coords = list(self.coords)
            p0 = coords[(i - 1) % 3]
            p1 = coords[(i) % 3]
            p2 = coords[(i + 1) % 3]
            angle = GenUtil.compute_angle(p0, p1, p2)
            if angle > 90.:
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
            if side_type == _Triangle.SUPERIMPOSED:
                superimposed = i

        p0_base = self.coords[superimposed % 3]
        p1_base = self.coords[(superimposed + 1) % 3]
        p_summit = self.coords[(superimposed + 2) % 3]

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
        p = (a + b + c) / 2.  # Calculates the half perimeter
        # Now use the Heron formula  A = sqrt(s(s-a)(s-b)(s-c)) where s = (a+b+c)/2
        area = (p * (p - a) * (p - b) * (p - c)) ** 0.5
        base = a

        # Calculate the height using area = (base*height)/2 ===> height= (2*area)/base
        self._height = (2. * area) / base

        if (angle_p1 < 90. and angle_p2 < 90.):
            acute = True
        else:
            acute = False

        return acute

    def get_centroid(self):
        """Calculate the position of the baricenter of the triangle

        The gravity centre is always 2/3 the distance between the middle of one side and the opposite angle

        Parameters: None

        Return value
            Tuple of (x,y) float representing the position of the centroid

        """

        coords = list(self.coords)
        mid_point = GenUtil.mid_point(coords[0], coords[1])
        centroid = GenUtil.rescale_vector(mid_point, coords[2], 1. / 3.)

        return centroid

    def demote_junction(self, side_number):
        """Demote a junction polygon if the skeleton of one of the side is considered as noise and calculate a new skeleton

        The new skeleton is calculated by joining to 2 remaining mid side

        Parameter:
            side_number: The side number of the triangle to demote. Value [0..2]
        """

        coords = list(self.coords)
        p0 = coords[(side_number + 1) % 3]
        p1 = coords[(side_number + 2) % 3]
        p2 = coords[(side_number + 3) % 3]

        mid_p0_p1 = GenUtil.mid_point(p0, p1)
        mid_p1_p2 = GenUtil.mid_point(p1, p2)

        centre_line = LineStringSc([mid_p0_p1, mid_p1_p2])

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
            self._centre_lines = []  # List of the centre lines
            mid_side_points = []  # List of the mid point on each side of the triangle
            internal_sides = []  # List of the number of the internal side
            external_sides = []  # List of the number of the external side

            nbr_adjacent = self.get_nbr_adjacent()

            coords = list(self.coords)
            for i in range(3):
                mid_side_points.append(GenUtil.mid_point(coords[i], coords[i + 1]))
                if self._side_type[i] == _Triangle.INTERNAL:
                    internal_sides.append(i)
                else:
                    external_sides.append(i)

            # Process each case depending on the number of internal side of the triangle
            if nbr_adjacent == 0:
                # Degenerated polygon with one triangle no skeleton line added
                pass

            if nbr_adjacent == 1:
                # Terminal triangle add line from the extremity of the triangle up to mid opposite side
                if internal_sides[0] == 0:
                    coords_line = [coords[2], mid_side_points[0]]
                if internal_sides[0] == 1:
                    coords_line = [coords[0], mid_side_points[1]]
                if internal_sides[0] == 2:
                    coords_line = [coords[1], mid_side_points[2]]

                self._centre_lines.append(LineStringSc(coords_line))

            if nbr_adjacent == 2:
                # Sleeve triangle skeleton added between the mid point of each chord
                internal_side0 = internal_sides[0]
                internal_side1 = internal_sides[1]
                self._centre_lines.append(LineStringSc([mid_side_points[internal_side0], mid_side_points[internal_side1]]))
                self._mid_triangle = GenUtil.mid_point(mid_side_points[internal_side0], mid_side_points[internal_side1])

            if nbr_adjacent == 3:
                # Junction triangle skeleton added.
                obtuse_angle = self._get_obtuse_angle()
                if (obtuse_angle is None):
                    # With an acute triangle a mid point is calculated in the middle of the triangle
                    centroid = self.get_centroid()
                    self._mid_triangle = centroid
                    for mid_side_point in mid_side_points:
                        self._centre_lines.append(LineStringSc([centroid, mid_side_point]))
                else:
                    # With an obtuse triangle the mid point is placed on the sided opposite to the obtuse angle
                    opposite_side = (obtuse_angle + 1) % 3
                    left_side = (opposite_side + 1) % 3
                    right_side = (opposite_side - 1) % 3
                    self._mid_triangle = mid_side_points[opposite_side]
                    self._centre_lines.append(LineStringSc([mid_side_points[opposite_side], mid_side_points[left_side]]))
                    self._centre_lines.append(LineStringSc([mid_side_points[opposite_side], mid_side_points[right_side]]))
        else:
            # Centre line was already calculated... nothing to do
            pass

        return self._centre_lines

    def get_nbr_adjacent(self):
        """Extract the number of side of the triangle which are adjacent to another polygon.

        Three scenarios are possible 0, 1, 2 ou 3 sides completely inside the polygon:
            3: The triangle is completely inside. This is called a Junction triangle
            2: The Triangle as on side that lies on the polygon. This is called a Sleeve triangle
            1: The Triangle has only one side completely inside the polygon. This is called a Terminal triangle.
            0: Special case where the polygon has only 3 sides

        Return value:
            - Number of side completely inside the polygon. Value between 0 and 3.
        """

        if self._nbr_internal is None:
            self._side_type = []
            self._nbr_adjacent = 0
            coords = list(self.coords)
            for i in range(3):
                tmp_line = LineString((coords[i], coords[i+1]))
                mid_point = tmp_line.interpolate(.5, normalized=True)
                bbox = GenUtil.build_bounding_box(GenUtil.ZERO, mid_point.bounds)
                adjacents = list(self._sc_scontainer.get_features(bounds=bbox, remove_features=[self]))
                tmp_superimposed = True
                for adjacent in adjacents:
                    if mid_point.distance(adjacent) <= GenUtil.ZERO:
                        tmp_superimposed = False
                        break
                if tmp_superimposed:
                    self._side_type.append(_Triangle.SUPERIMPOSED)
                else:
                    self._side_type.append(_Triangle.INTERNAL)
                    self._nbr_adjacent += 1
        else:
            # nbr_internal has already been calculated... nothing to do
            pass

        return self._nbr_adjacent

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
                self._category = ChordalAxis.OTHER
            else:
                nbr_internal = self.get_nbr_internal()
                if nbr_internal == 2 and self._is_acute_triangle():
                    # The triangle is a sleeve triangle and is acute
                    # In this case the height of the triangle is also the width of the bottleneck
                    if self._height < minimal_width:
                        neighbour = True
                        self._width = self._height
                    else:
                        neighbour = False
                else:
                    # Check through a spatial search if there are any neighbous
                    (neighbour, self._width) = _Triangle.line_segments.check_chordal_axis(minimal_width / 2.,
                                                                                          self._mid_triangle)
                if (neighbour):
                    extremity = False
                    for i, type in enumerate(self._side_type):
                        if type == _Triangle.INTERNAL:
                            coord0 = self.coords[i]
                            coord1 = self.coords[i + 1]
                            # Check if the triangle is located near the extremity of the polygon
                            # A triangle near an extremity of a polygon is not considered as a bottleneck
                            extremity = extremity or _Triangle.perimeters.is_extremity(coord0, coord1, minimal_width)
                    if extremity:
                        self._category = ChordalAxis.OTHER
                    else:
                        self._category = ChordalAxis.BOTTLENECK
                else:
                    self._category = ChordalAxis.OTHER

            if self._category != ChordalAxis.BOTTLENECK:
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
            for j in range(nbr_coord_ring - 1):
                p0 = lst_coords[j]
                p1 = lst_coords[j + 1]
                line = LineStringSc((lst_coords[j], lst_coords[j+1]))
#                line.sb_geom_type = GenUtil.LINE_STRING
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
        mid_p0_p1 = GenUtil.mid_point(p0, p1)
        # Search for line segment there
        b_box = GenUtil.build_bounding_box(self._search_tolerance, mid_p0_p1)
        lines = self._s_container.get_features(bounds=b_box)

        present = False
        for line in lines:
            # Check if the first/vertice of the line correpond to the p1, p2 of the triangle if so
            # there is exactly one line there
            line_coords = list(line.coords)
            if ((GenUtil.distance(line_coords[0], p0) < GenUtil.ZERO and
                 GenUtil.distance(line_coords[1], p1) < GenUtil.ZERO) or
                    (GenUtil.distance(line_coords[0], p1) < GenUtil.ZERO and
                     GenUtil.distance(line_coords[1], p0) < GenUtil.ZERO)):
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
            if distance < min_distance_0:
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
        if chord_distance < tolerance * 2.:
            neighbours = True
        else:
            neighbours = False
            chord_distance = 1.0E+99

        return (neighbours, chord_distance)


class PerimeterDistance(object):
    """This class allows to calculate the distance between 2 coordinates on a line closed string


    Internal data structure are maintained to accelerate the computation

    """

    _ID_RING = 'id_ring'
    _ID_COORD = 'id_coord'

    def __init__(self, rings, search_tolerance):
        """Load the internal structure with the ring information

        Parameters:
           - rings: List of closed LineStringSc of a polygon
           - search_tolerance: Search tolerance can vary depending on the dynamic of the data set from lat-lon to Lambert conformal

        Return value: None
        """

        self._search_tolerance = search_tolerance

        for ring in rings:
            if (isinstance(ring, LinearRing)):
                pass
            else:
                raise Exception("Can only work on LineStringSc")

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
        for i in range(1, nbr_coords):
            point = PointSc(lst_coords[i])
            point.sb_geom_type = GenUtil.POINT
            point.sb_id_ring = id_ring
            point.sb_id_coord = i
#            point.ma_properties[PerimeterDistance._ID_RING] = id_ring
#            point.ma_properties[PerimeterDistance._ID_COORD] = i
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
            if (i == 0):
                # Do not add the first vertex in the container
                previous_coord = current_coord
            else:
                dist = GenUtil.distance(previous_coord, current_coord)
                last_dist = cum_distance[-1]
                cum_distance.append(last_dist + dist)
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
        points = list(self.s_cont_points.get_features(bounds=b_box))
        nbr_points = len(points)
        if nbr_points == 0:
            # Nothing is found
            id_ring = -1
            id_coord = -1
            if (raise_exception):
                raise Exception("Integrity problem at coordinate: (%f,%f)" % (coord[0], coord[1]))
        else:
            if nbr_points == 1:
                # There is only one point
                point = points[0]
            else:
                # Take the closest point
                min_dist = 1.0E+99
                for p in points:
                    dist = GenUtil.distance(p.coords[0], coord)
                    if dist < min_dist:
                        point = p
                        dist = min_dist
            id_ring = points[0].sb_id_ring
            id_coord = points[0].sb_id_coord

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
                start, end = id_coord0, id_coord1
                i = start
                while (i != end):
                    sub_coords1.append(lst_coords[i])
                    i = (i + 1) % nbr_coords
                sub_coords1.append(lst_coords[i])

                # Loop to extract coordinate from last to first
                start, end = id_coord1, id_coord0
                i = start
                while i != end:
                    sub_coords2.append(lst_coords[i])
                    i = (i + 1) % nbr_coords
                sub_coords2.append(lst_coords[i])

                # Take the smallest list
                if len(sub_coords1) < len(sub_coords2):
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

        if id_ring0 != -1 and id_ring0 != -1:

            if id_ring0 == id_ring1:
                # Coordinates are on the same ring
                cumm_distance = self.lst_cumm_distance[id_ring0]
                if ((id_coord0 < 0 or id_coord0 > len(cumm_distance)) or
                        (id_coord1 < 0 or id_coord1 > len(cumm_distance))):
                    raise Exception("Internal Error...")

                if id_coord0 < id_coord1:
                    i, j = id_coord0, id_coord1
                else:
                    i, j = id_coord1, id_coord0

                # Extracting the 2 perimetre formed by the line (p1,p2) cutting the polygon in two parts
                peri_distance_1 = cumm_distance[j] - cumm_distance[i]
                peri_distance_2 = cumm_distance[-1] - peri_distance_1
                peri_distance = min(peri_distance_1, peri_distance_2)

                # Check if the area formed is smaller than a minimal_width area
                dist_0_1 = GenUtil.distance(coord0, coord1)
                if dist_0_1 + peri_distance < 2 * dist_0_1 + 2 * minimal_width:
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


class GeoSimException (Exception):
    """
    This is the base exception class for genmetal algorithms
    """

    def __init__(self, *arguments, **keywords):
        Exception.__init__(self, *arguments, **keywords)


class InternalError (GeoSimException):
    """
    This exception is raised when an internal error as occurred
    """

    def __init__(self, *param_names):

        """
        Initialise an Invalid Geometry Error

        *Parameters*:
            - param_names: one or more names of invalid parameters

        """

        GeoSimException.__init__(self, *param_names)
