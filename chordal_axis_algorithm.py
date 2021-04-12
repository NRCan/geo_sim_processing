# -*- coding: utf-8 -*-

# /***************************************************************************
# chordal_axis_algorithm.py
# ----------
# Date                 : January 2021
# copyright            : (C) 2020 by Natural Resources Canada
# email                : daniel.pilon@canada.ca
#
#  ***************************************************************************/
#
# /***************************************************************************
#  *                                                                         *
#  *   This program is free software; you can redistribute it and/or modify  *
#  *   it under the terms of the GNU General Public License as published by  *
#  *   the Free Software Foundation; either version 2 of the License, or     *
#  *   (at your option) any later version.                                   *
#  *                                                                         *
#  ***************************************************************************/


import sys
from math import atan, degrees
from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsFields,
                       QgsFeature,
                       QgsFeatureRequest,
                       QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingException,
                       QgsProcessingParameterFeatureSource,
                       QgsWkbTypes,
                       QgsGeometry,
                       QgsFeatureSink,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean)
from qgis.core import QgsPoint, QgsLineString, \
                      QgsSpatialIndex, QgsPolygon, QgsRectangle, \
                      QgsMultiLineString, QgsVectorLayer, QgsApplication
import os
import inspect
import processing
from qgis._3d import Qgs3DAlgorithms
from qgis.PyQt.QtGui import QIcon


class ChordalAxisAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer,
    creates some new layers and returns some results.
    """

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        """
        Returns a new copy of your algorithm.
        """
        return ChordalAxisAlgorithm()

    def name(self):
        """
        Returns the unique algorithm name.
        """
        return 'chordalaxis'

    def displayName(self):
        """
        Returns the translated algorithm name.
        """
        return self.tr('Chordal axis')

    def group(self):
        """
        Returns the name of the group this algorithm belongs to.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs
        to.
        """
        return ''

    def shortHelpString(self):
        """
        Returns a localised short help string for the algorithm.
        """
        help_str = """
    <b>Chordal Axis</b>
    ChordalAxis is a geospatial tool that creates a skeleton (the center line). ChordalAxis creates a \
    triangulation (using QGIS Tessellate tools) and use it to extract the chordal axis (center line).

    <b>Usage</b>
    <u>Input</u>: A polygon layer to extract the chordal axis.

    <u>Correct skeleton</u>:  Correct the skeleton for small centre line, T junction and X junction. Useful in the case \
    of long any narrow polygon (ex.: polygonized road network)
    
    <b>Output</b>
    <u>Chordal axis</u>: The center line of the polygons.
    <u>Tessellate</u>: The result of the triangulation.

    """
        return self.tr(help_str)

    def icon(self):
        cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]
        icon = QIcon(os.path.join(os.path.join(cmd_folder, 'logo.png')))
        return icon

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and outputs of the algorithm.
        """
        # 'INPUT' is the recommended name for the main input
        # parameter.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                'INPUT',
                self.tr('Input layer'),
                types=[QgsProcessing.TypeVectorAnyGeometry]
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                'CORRECTION',
                self.tr('Correction'),
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                'OUTPUT',
                self.tr('Chordal axis')
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                'TRIANGLES',
                self.tr('Tessellation')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        in_source = self.parameterAsSource(parameters, "INPUT", context)
        correction = self.parameterAsBool(parameters, "CORRECTION", context)

        if in_source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        # Transform the in source into a vector layer
        in_vector_layer = in_source.materialize(QgsFeatureRequest(), feedback)

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                               QgsFields(),
                                               QgsWkbTypes.LineString,
                                               in_vector_layer.sourceCrs() )

        (sink_t, dest_id_t) = self.parameterAsSink(parameters, "TRIANGLES", context,
                                                   QgsFields(),
                                                   QgsWkbTypes.MultiPolygon,
                                                   in_vector_layer.sourceCrs())

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))
        if sink_t is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "TRIANGLES"))

        nbr_polygon = in_vector_layer.featureCount()
        nbr_centre_line = 0
        total = 100.0 / in_vector_layer.featureCount() if in_vector_layer.featureCount() else 0
        qgs_multi_triangles = ChordalAxis.tessellate_polygon(in_vector_layer, feedback)  # Tessellate polygons
        for i, qgs_multi_triangle in enumerate(qgs_multi_triangles):
            sink_t.addFeature(qgs_multi_triangle, QgsFeatureSink.FastInsert)  # Add to Triangle sink
            # Call the chordal axis
            try:
                ca = ChordalAxis(qgs_multi_triangle, GenUtil.ZERO)
                if correction:
                    ca.correct_skeleton()
                centre_lines = ca.get_skeleton()
            except Exception:
                import traceback
                traceback.print_exc()

            # Load the centre line (skeleton) in the sink
            for line in centre_lines:
                out_feature = QgsFeature()
                geom_feature = QgsGeometry(line.clone())
                out_feature.setGeometry(geom_feature)
                sink.addFeature(out_feature, QgsFeatureSink.FastInsert)  # Add to skeleton sink
                nbr_centre_line += 1

            if feedback.isCanceled():
                break

            feedback.setProgress(int(i * total))

        # Push some output statistics
        feedback.pushInfo("Number of input polygons: {0}".format(nbr_polygon))
        feedback.pushInfo("Number of centre lines: {0}".format(nbr_centre_line))

        return {"OUTPUT": dest_id, "TRIANGLES": dest_id_t}





#Chordal Axis Algorithm
#----------------------

"""
General classes and utilities needed for ChordalAxis.
"""
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
    def merge_line_string(lst_line_string):
        """Merge a list of QgsLineString and return a list of QgsLineString merged

        Note If all the QgsLineString are merged the lise will contain only one element

        Parameters
        ----------
        lst_line_string: List
            List of QgsLineString to merge
        Returns
        -------
        List of merged LineString
        """

        # Accumulate all the QgsLineString into one QgsMultiLineString
        qgs_mline = QgsMultiLineString()
        for line_string in lst_line_string:
            qgs_mline.addGeometry(line_string.clone())

        # Merge the QgsLineString
        geom_mline = QgsGeometry(qgs_mline)
        merged_line_string = geom_mline.mergeLines()

        # Extract the QgsLineString from the qgsMultiLineString
        out_line_string = []
        for qgs_line in merged_line_string.constParts():
            out_line_string.append(qgs_line.clone())

        return out_line_string


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

class SpatialContainer(object):
    """This class manages the spatial features and a spatial index using the QgsSpatialIndex.

    """

    # Class variable that contains the Spatial Container Internal ID
    _sc_id = 1

    def __init__(self):
        """Create an object of type SpatialContainer
        The init will create one container for the feature a dictionary and one
        container for the spatial index (QgsSpatialIndex
        *Parameters*: None
        *Returns*: *None*
        """

        self._features = {}  # Container to hold the QgsFeature

        self._index = QgsSpatialIndex()

    def add_features(self, features):
        """Adds all the features in the container and update the spatial index with the feature's bound
        *Parameters*:
            - feature: A spatial feature derives from PointSc, LineStringSc
        *Returns*: *None*
        """

        for feature in features:
            # Check if the type is valid
            if issubclass(feature.__class__, (_TriangleSc)):
                pass
            else:
                raise QgsProcessingException("Cannot add feature in the QgsSpatialIndex")

            qgs_rect = QgsRectangle(feature.rectangle)
            qgs_rect.grow(GenUtil.ZERO)

            # Container unique internal counter
            SpatialContainer._sc_id += 1

            # Add the spatial id to the feature
            feature._sc_id = SpatialContainer._sc_id
            feature._sc_scontainer = self

            # Add the feature in the feature container
#            self._features[feature._sc_id] = feature

            self._features[SpatialContainer._sc_id] = feature

            # Load all the features at the spatial index
            self._index.addFeature(SpatialContainer._sc_id, qgs_rect)

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

        Because it's impossible the STRtree structure is an immutable structure we mark the feature to
        delete by changing the ID (_sc_id) to a negative value.

        *Parameters*:
            - feature: The feature to delete in the spatial container
        *Returns*:
            None
        """

#        print ("Implement real del_feature")
        if feature._sc_id >= 1:
            feature._sc_id = -(feature._sc_id)
        else:
            raise GeoSimException("Cannot delete a feature already deleted")

        return

    def del_features(self, features):
        """Delete a list of features in the spatial container

        *Parameters*:
            - features: list of features to delete
        *Returns*:
            None
        """

        for feature in features:
            self.del_feature(feature)

        return

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

    def get_features(self, qgs_rectangle=None, remove_features=[]):
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
        if (qgs_rectangle != None):
            # Extract features by bounds
            qgs_rectangle.grow(GenUtil.ZERO)
            features_id = self._index.intersects(qgs_rectangle)
        else:
            features_id = [self._features.keys()]

        # Extract the feature from the features ID
        features = [self._features[feature_id] for feature_id in features_id if feature_id not in remove_features]

#        for feature in features:
#            yield feature
        features = [feature for feature in features if feature._sc_id >= 0]

        return features

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

    ANGLE_JUNCTION_T = 45.  # Delta used to test if 2 branches or contiguous
    SEARCH_TOLERANCE = None

    @staticmethod
    def tessellate_polygon(source, feedback):
        """Method to tesellate a polygon

        This method pre process the input source by doing the following task a
         - execute "Multi to single part" processing in order to accept even multipolygon
         - excecute "Drop  Z and M values" processing as they are not usefull
         - Validate if the resulting layer is Polygom
         - execute "Tessellate" processing to create the triangles
         - excecute "Drop  Z and M values" processing as they are not usefull

        :param source:
        :return:
        """

        # Execute MultiToSinglePart processing
        params = {'INPUT': source,
                  'OUTPUT': 'memory:'}
        result_ms = processing.run("native:multiparttosingleparts", params, feedback=feedback)
        ms_part_layer = result_ms['OUTPUT']

        # Execute Drop Z M processing
        params = {'INPUT': ms_part_layer,
                  'DROP_M_VALUES': True,
                  'DROP_Z_VALUES': True,
                  'OUTPUT': 'memory:'}
        result_drop_zm = processing.run("native:dropmzvalues", params, feedback=feedback)
        drop_zm_layer = result_drop_zm['OUTPUT']

        # Validate the resulting geometry is Polygon
        if drop_zm_layer.wkbType() == QgsWkbTypes.Polygon:
            source_crs = source.sourceCrs()
            epsg_id = source_crs.authid()
            if epsg_id is None or epsg_id == "":
                uri_p = "Polygon"  # Should not happen... just in case...
                uri_mp = "MultiPolygon"
            else:
                uri_p = "Polygon?" + epsg_id
                uri_mp = "MultiPolygon?" + epsg_id
            QgsApplication.processingRegistry().addProvider(Qgs3DAlgorithms(QgsApplication.processingRegistry()))

            qgs_multi_triangles = []  # Contains the tessellated features
           # Loop over each Polygon to create the tesselation (Triangles)
            qgs_features = drop_zm_layer.getFeatures()
            for qgs_feature in qgs_features:
                # Process one feature at the time so it's possible to apply some correction if the tessellation crashes
                single_feature_layer = QgsVectorLayer(uri_p, "temporary_polygon", "memory")
                data_provider = single_feature_layer.dataProvider()
                data_provider.addFeatures([qgs_feature])
                single_feature_layer.updateExtents()
                try:
                    params = {'INPUT': single_feature_layer,
                              'OUTPUT': 'memory:'}
                    # Execute Tessellate processing to create the triangles
                    result_tessellate = processing.run("3d:tessellate", params, feedback=feedback)
                except Exception:
                    # Code to be place here to try cleaning the feature when it crashes
                    raise QgsProcessingException("Unable to tessellate feature: {}".format(qgs_feature.id()))

                tessellate_layer = result_tessellate['OUTPUT']

                qgs_features = tessellate_layer.getFeatures()
                for qgs_feature in qgs_features:
                    qgs_geom = qgs_feature.geometry()
                    qgs_geoms = qgs_geom.coerceToType(QgsWkbTypes.MultiPolygon)  # Drop Z and M
                    qgs_feature.setGeometry(qgs_geoms[0])
                    qgs_multi_triangles.append(qgs_feature)
        else:
            # Unable to process this geometry
            geom_type_str = QgsWkbTypes.geometryDisplayString(drop_zm_layer.geometryType())
            raise QgsProcessingException ("Unable to process geometry: {}". format(geom_type_str))

        return (qgs_multi_triangles)









    def __init__(self, in_feature, search_tolerance=GenUtil.ZERO):
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

        abs_geom = in_feature.geometry()
        multi_geom = abs_geom.constGet().clone()

        lst_line_qgs_point = self._validate_triangles(multi_geom)

        # Transform the polygons (triangles) into _TriangleSc to be loaded in SpatialContainer
        lst_triangle = []
        for line_qgs_point in lst_line_qgs_point:
            triangle = _TriangleSc(line_qgs_point)
            lst_triangle.append(triangle)

        # Create spatial container
        self.s_container = SpatialContainer()

        # Load triangles
        self.s_container.add_features(lst_triangle)

        # Load some class variables
        _TriangleSc.s_container = self.s_container

        # Build the cluster (group of triangles part of a polygon)

#        cluster = {triangle.id: triangle for triangle in lst_triangle}
#        self.triangle_clusters = [cluster]
        self.triangle_clusters = [lst_triangle]

#        self.triangle_clusters = self._build_clusters()

        self.nbr_polygons = len(self.triangle_clusters)
        self.nbr_triangles = len(lst_triangle)

        # Initialise stats value
        self.nbr_lines_pruned = 0
        self.nbr_iteration = 0
        self.nbr_t_junction = 0
        self.nbr_x_junction = 0

        return

    def _validate_triangles(self, multi_geom):
        """ Validate the each triangle
        Parameters
        ----------
        lst_triangle : QGIS MultiGeometry
            MultiGeometry contaiing the triangles of one polygon
        Return
        ------
        None

        Exception if the geometry is invalid
        """

        #
        multi_geom = QgsGeometry(multi_geom.clone())  # create a geometry in order to iterate over it
        lst_line_qgs_point = []
        if multi_geom.wkbType() in [QgsWkbTypes.MultiPolygon, QgsWkbTypes.MultiPolygonZ]:
            lst_pol_qgs_point = multi_geom.asMultiPolygon()  # Tranform the MultiPolygon into a list of QgsPoint
            for pol_qgs_point in lst_pol_qgs_point:
                if len(pol_qgs_point) == 1:  # Validate that the polygon has only an exterior
                    if len(pol_qgs_point[0]) == 4:  # Validate that the polygon has 4 QgsPoint (a triangle)
                        #for qgs_point in pol_qgs_point[0]:
                        # Transform the QgsPointXY ==> QgsPoint of the exterior ring
                        tmp_qgs_point = [QgsPoint(qgs_point.x(), qgs_point.y()) for qgs_point in pol_qgs_point[0]]
                        lst_line_qgs_point.append(tmp_qgs_point)  # Only keep the exterior ring
                    else:
                        QgsProcessingException("Polygon must have 4 vertice")
                else:
                    raise QgsProcessingException("Polygon cannot have interior")
            return lst_line_qgs_point
        else:
            raise QgsProcessingException("Can only process MultiPolygon")

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
           List of clusters forming the different polygon
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
            for triangle in triangle_cluster:
                centre_lines += triangle.centre_line

        merged_centre_lines = GenUtil.merge_line_string(centre_lines)
#        # Accumulate all the QgsLineString into one QgsMultiLineString
#        qgs_mline = QgsMultiLineString()
#        for centre_line in centre_lines:
#            qgs_mline.addGeometry(centre_line)
#
#        # Merge the QgsLineString
#        geom_mline = QgsGeometry(qgs_mline)
#        merged_lines = geom_mline.mergeLines()
#
#        # Extract the QgsLineString from the qgsMultiLineString
#        qgs_centre_lines = []
#        for qgs_line in merged_lines.constParts():
#            qgs_centre_lines.append(qgs_line.clone())#
#
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
                for triangle in triangle_cluster:  # Loop over each triangle
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
            for triangle in triangle_cluster:  # Loop over each triangle of one polygon
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
            for triangle in triangle_cluster:  # Loop over each triangle of one polygon
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
            # The junction triangle is automaticaly a T junction
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

            # If the last triangle in a branch is junction and within a certain tolerance it's a candidate for X Junction
            if last_triangle.type in (ChordalAxis.JUNCTION, ChordalAxis.JUNCTION_T) and \
                    branch.length < min(current_junction.width, last_triangle.width) * ChordalAxis.JUNCTION_X_LENGTH:
                # Merge the triangle in the branch to form only one polygon
                triangles = [current_junction] + branch.triangle_in_branch
                pol_triangles = [QgsGeometry(QgsPolygon(triangle.qgs_line_string.clone())) for triangle in triangles]
                merged_pol = QgsGeometry.unaryUnion(pol_triangles)
                if merged_pol.wkbType() == QgsWkbTypes.Polygon: ## == GenUtil.POLYGON:  # Merged the triangle to form only one polygon
                    merged_line = QgsGeometry(merged_pol.constGet().clone().boundary()) #LineString(merged_pol.exterior.coords)

                    # Detect which mid side point we must keep (we must keep four only)
                    mid_pnt_sides = current_junction.mid_pnt_sides + last_triangle.mid_pnt_sides
                    new_mid_pnt_sides = []
                    for mid_pnt_side in mid_pnt_sides:
                        if mid_pnt_side.distance(merged_line) < ChordalAxis.SEARCH_TOLERANCE:
                            new_mid_pnt_sides.append(mid_pnt_side)

                    # Validate the center line
                    qgs_geom_x_centroid = merged_pol.centroid()
                    qgs_pnt_x_centroid = qgs_geom_x_centroid.constGet().clone()
                    if self.validate_x_junction(merged_pol, qgs_pnt_x_centroid, new_mid_pnt_sides):
                        junction_x_info = Holder(first_junction=current_junction, last_junction=branch.triangle_in_branch[-1],
                                                 sleeve_in_branch=branch.triangle_in_branch[0:-1], mid_pnt_sides=new_mid_pnt_sides,
                                                 x_centroid=qgs_pnt_x_centroid)
                        junction_x_infos.append(junction_x_info)
                else:
                    # Invalid merging nothing to do (should not happen)
                    pass
            else:
                # Not a triangle to process for X junction
                pass

        return junction_x_infos

    def validate_x_junction(self, merged_pol, x_centroid, new_mid_pnt_sides):
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

        buf_merged_pol = merged_pol.buffer(.01, 3)

        status = True
        for mid_pnt_side in new_mid_pnt_sides:
            line = QgsGeometry(QgsLineString((x_centroid, mid_pnt_side.constGet().clone())))
            if line.crosses(buf_merged_pol):
                # Centre line not valid
                status = False
                break

        # Validate the angle between the line
        # Validation not yes implemented...
        for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)):
            x0,y0 = new_mid_pnt_sides[i].constGet().x() - x_centroid.x(), new_mid_pnt_sides[i].constGet().y() - x_centroid.y()
            x1,y1 = new_mid_pnt_sides[j].constGet().x() - x_centroid.x(), new_mid_pnt_sides[j].constGet().y() - x_centroid.y()
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
#                print ("à corriger par un dictionnaire...")
                for i, to_del in enumerate(cluster_triangle):
                    if triangle.id == to_del.id:
                        del cluster_triangle[i]
                        break
#                del cluster_triangle[triangle.id]

        return len(del_branches)


class _TriangleSc():
    """LineString1 specialization to be included in the SpatialContainer"""

    id = 0

    def __init__(self, line_qgs_point):
#        super().__init__(coords)

        self.qgs_pnt_p0 = QgsPoint(line_qgs_point[0].x(), line_qgs_point[0].y())
        self.qgs_pnt_p1 = QgsPoint(line_qgs_point[1].x(), line_qgs_point[1].y())
        self.qgs_pnt_p2 = QgsPoint(line_qgs_point[2].x(), line_qgs_point[2].y())
        self.qgs_line_string = QgsLineString(line_qgs_point)
        self.geom_line_string = QgsGeometry(self.qgs_line_string.clone())
        self.rectangle = self.qgs_line_string.calculateBoundingBox()
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
            List containing the mid QgsPoint of each side of the triangle
        """

        try:
            # Implement attribute late calculation
            return self._mid_pnt_sides
        except AttributeError:
            if not hasattr(self, 'junction_x_mid_pnt_sides'):
                # Calculate the mid point of each side of the triangle
 #               coords = list(self.coords)
                line = QgsLineString(self.qgs_pnt_p0, self.qgs_pnt_p1)
                mid_pnt_side_0 = line.interpolatePoint(line.length()/2.)
#                mid_pnt_side_0 = LineString([coords[0], coords[1]]).interpolate(0.5, normalized=True)
                line = QgsLineString(self.qgs_pnt_p1, self.qgs_pnt_p2)
                mid_pnt_side_1 = line.interpolatePoint(line.length() / 2.)
#                mid_pnt_side_1 = LineString((coords[1], coords[2])).interpolate(0.5, normalized=True)
                line = QgsLineString(self.qgs_pnt_p2, self.qgs_pnt_p0)
                mid_pnt_side_2 = line.interpolatePoint(line.length() / 2.)
#                mid_pnt_side_2 = LineString((coords[2], coords[0])).interpolate(0.5, normalized=True)
                self._mid_pnt_sides = [QgsGeometry(mid_pnt_side_0),
                                       QgsGeometry(mid_pnt_side_1),
                                       QgsGeometry(mid_pnt_side_2)]
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
            lst_length = [qgs_line.length() for qgs_line in self.centre_line]
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
                potential_triangles = _TriangleSc.s_container.get_features(qgs_rectangle=mid_pnt_side.boundingBox(), remove_features=[self])

                # Find the closest triangle
                triangles = [(triangle, mid_pnt_side.distance(triangle.geom_line_string)) for triangle in potential_triangles]
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

            self._centre_line = []

            # Process each case depending on the number of internal side of the triangle
            if self.type == ChordalAxis.ISOLATED:
                # Degenerated polygon with one triangle no skeleton line added
                pass

            elif self._type == ChordalAxis.TERMINAL:
                # Terminal triangle add line from the extremity of the triangle up to mid opposite side
                if self.adjacent_sides_ref[0] != None:
                    coords_line = [self.qgs_pnt_p2, self.mid_pnt_sides[0]]
                if self.adjacent_sides_ref[1] != None:
                    coords_line = [self.qgs_pnt_p0, self.mid_pnt_sides[1]]
                if self.adjacent_sides_ref[2] != None:
                    coords_line = [self.qgs_pnt_p1, self.mid_pnt_sides[2]]

                coords_line[1] = coords_line[1].constGet()  # Return the QgsPoint
                self._centre_line.append(QgsLineString(coords_line))

            elif self.type == ChordalAxis.SLEEVE:
                # Sleeve triangle skeleton added between the mid point of side adjacent to another triangle
                mid_pnt = []
                for i, adjacent_side_ref in enumerate(self._adjacent_sides_ref):
                    if adjacent_side_ref != None:
                        mid_pnt.append(self.mid_pnt_sides[i].constGet())
                self._centre_line.append(QgsLineString([mid_pnt[0], mid_pnt[1]]))

            elif self.type == ChordalAxis.SLEEVE_X:
                # No center line to create
                pass

            elif self.type == ChordalAxis.JUNCTION:
                # Regular triangle T type. Centroid is the centroid of the triangle
                pnt_x = (self.qgs_pnt_p0.x() + self.qgs_pnt_p1.x() + self.qgs_pnt_p2.x()) / 3.
                pnt_y = (self.qgs_pnt_p0.y() + self.qgs_pnt_p1.y() + self.qgs_pnt_p2.y()) / 3.
                centroid = QgsPoint(pnt_x, pnt_y)
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    self._centre_line.append(QgsLineString(centroid, mid_side_pnt.constGet()))

            elif self.type == ChordalAxis.JUNCTION_T:
                # Corrected triangle T. Centroid is the middle point between the 2 aligned branches
#######                voir https://brilliant.org/wiki/triangles-centroid/ finding the centroid
                pnt0 = self.mid_pnt_sides[self.junction_side_a].constGet().clone()
                pnt1 = self.mid_pnt_sides[self.junction_side_b].constGet().clone()
                line = QgsLineString(pnt0.clone(), pnt1.clone())
                centroid = line.interpolatePoint(line.length() / 2.)
#
#                line = QgsLineString()
#               mid_pnt_side_1 = line.interpolatePoint(line.length() / 2.)
#
#
#                centroid = [pnt.x, pnt.y]
                # Create the centre line
                for mid_side_pnt in self.mid_pnt_sides:
                    self._centre_line.append(QgsLineString(centroid, mid_side_pnt.constGet().clone()))

            elif self.type == ChordalAxis.JUNCTION_X_FIRST:
                #  create the center line
                for mid_pnt_side in self.junction_x_mid_pnt_sides:
                    self._centre_line.append(QgsLineString(self.junction_x_centroid, mid_pnt_side.constGet().clone()))

            elif self.type == ChordalAxis.JUNCTION_X_LAST:
                # No center line to create
                pass

            else:
                raise GeoSimException ("Unknow triangle type: {1}".format(self.type))

            return self._centre_line

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
                self.length += next_triangle.centre_line[0].length()  # Add the length
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
            merged_line = GenUtil.merge_line_string(lines)
            line = merged_line[0]
            start_pnt = line.startPoint()
            end_pnt = line.endPoint()
            x0, y0 = start_pnt.x(), start_pnt.y()
            x1, y1 = end_pnt.x(), end_pnt.y()
#            if len(merged_line) == 1:
#                # Extract the angle formed by the first and last coordinate
#                x0, y0 = merged_line.coords[0][0], merged_line.coords[0][-1]
#                x1, y1 = merged_line.coords[-1][0], merged_line.coords[-1][-1]
#            else:
#                # There was an error in line merge so just take the first line (less good but should not happen often)
#                line  = lines[0]
#                x0, y0 = line.coords[0][0], line.coords[0][-1]
#                x1, y1 = line.coords[-1][0], line.coords[-1][-1]

            # Checked that the first coordinate of the line is located on the triangle
            geom_p0 = QgsGeometry(QgsPoint(x0,y0))
            if self.current_triangle.geom_line_string.distance(geom_p0) < ChordalAxis.SEARCH_TOLERANCE:
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
