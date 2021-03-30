# -*- coding: utf-8 -*-
# pylint: disable=no-name-in-module
# pylint: disable=too-many-lines
# pylint: disable=useless-return
# pylint: disable=too-few-public-methods

# /***************************************************************************
# reduce_bend_algorithm.py
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

"""
QGIS Plugin for Bend reduction
"""


import os
import inspect
import sys
import math
from abc import ABC, abstractmethod
from qgis.PyQt.QtCore import QCoreApplication

from qgis.PyQt.QtGui import QIcon

from qgis.core import (QgsFeature, QgsProcessing, QgsProcessingAlgorithm, QgsProcessingParameterDistance,
                       QgsProcessingParameterFeatureSource, QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean, QgsFeatureSink, QgsFeatureRequest, QgsPoint,
                       QgsPointXY, QgsLineString, QgsPolygon, QgsWkbTypes, QgsSpatialIndex, QgsGeometry,
                       QgsGeometryUtils, QgsRectangle, QgsProcessingException, QgsMultiPolygon)

import processing


class ReduceBendAlgorithm(QgsProcessingAlgorithm):
    """Main class defining the Reduce Bend as a QGIS processing algorithm.
    """

    def tr(self, string):  # pylint: disable=no-self-use
        """Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):  # pylint: disable=no-self-use
        """Returns a new copy of the algorithm.
        """
        return ReduceBendAlgorithm()

    def name(self):  # pylint: disable=no-self-use
        """Returns the unique algorithm name.
        """
        return 'reducebend'

    def displayName(self):  # pylint: disable=no-self-use
        """Returns the translated algorithm name.
        """
        return self.tr('Reduce bend')

    def group(self):
        """Returns the name of the group this algorithm belongs to.
        """
        return self.tr(self.groupId())

    def groupId(self):  # pylint: disable=no-self-use
        """Returns the unique ID of the group this algorithm belongs to.
        """
        return ''

    def shortHelpString(self):
        """Returns a localised short help string for the algorithm.
        """
        help_str = """
    Reduce bend is a geospatial simplification and generalization tool for lines and polygons. The \
    particularity of this algorithm is that for each line or polygon it analyzes its bends (curves) and \
    decides which one to reduce, trying to emulate what a cartographer would do manually \
    to simplify or generalize a line. Reduce bend will accept lines and polygons as input.  Reduce bend will \
    preserve the topology (spatial relations) within and between the features during the bend reduction. \
    Reduce bend also accept multi lines and multi polygons but will output lines and polygons.

    <b>Usage</b>
    <u>Input layer</u> : Any LineString or Polygon layer.  Multi geometry are transformed into single geometry.
    <u>Diameter tolerance</u>: Theoretical diameter of a bend to remove.
    <u>Smooth line</u>: If you want to smooth the reduced bends (when possible).
    <u>Exclude hole</u>: If you want to exclude (delete )holes below the diameter of the bend.
    <u>Exclude polygon</u>: If you want to exclude (delete) polygon below the diameter of the bend.
    <u>Reduced bend</u> : Output layer of the algorithm.

    <b>Rule of thumb for the diameter tolerance</b>
    Reduce bend can be used for line simplifying in the context of line generalization. The big \
    question will often be what diameter should we use? A good starting point is the cartographic rule of \
    thumb -- the .5mm on the map -- which says that the minimum distance between two lines should be \
    greater than 0.5mm on a paper map. So to simplify (generalize) a line for representation at a scale of \
    1:50 000 for example a diameter of 25m should be a good starting point

    """

        return self.tr(help_str)

    def icon(self):  # pylint: disable=no-self-use
        """Define the logo of the algorithm.
        """

        cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]
        icon = QIcon(os.path.join(os.path.join(cmd_folder, 'logo.png')))
        return icon

    def initAlgorithm(self, config=None):  # pylint: disable=unused-argument
        """Define the inputs and outputs of the algorithm.
        """

        # 'INPUT' is the recommended name for the main input parameter.
        self.addParameter(QgsProcessingParameterFeatureSource(
                          'INPUT',
                          self.tr('Input layer'),
                          types=[QgsProcessing.TypeVectorAnyGeometry]))

        # 'TOLERANCE' to be used bor bend reduction
        self.addParameter(QgsProcessingParameterDistance(
                          'TOLERANCE',
                          self.tr('Diameter tolerance'),
                          defaultValue=0.0,
                          parentParameterName='INPUT'))  # Make distance units match the INPUT layer units

        # 'SMOOTH' to create rounded line instead of straight line (when possible)
        self.addParameter(QgsProcessingParameterBoolean(
            'SMOOTH',
            self.tr('Smooth line'),
            defaultValue=False))

        # 'EXCLUDE_POLYGON' to delete polygon (including holes) below the tolerance
        self.addParameter(QgsProcessingParameterBoolean(
                          'EXCLUDE_POLYGON',
                          self.tr('Exclude polygon'),
                          defaultValue=True))

        # 'EXCLUDE_POLYGON' to delete polygon holes below the tolerance
        self.addParameter(QgsProcessingParameterBoolean(
                          'EXCLUDE_HOLE',
                          self.tr('Exclude hole'),
                          defaultValue=True))

        # 'VERBOSE' mode for more output information
        self.addParameter(QgsProcessingParameterBoolean(
            'VERBOSE',
            self.tr('Verbose'),
            defaultValue=False))

        # 'OUTPUT' for the results
        self.addParameter(QgsProcessingParameterFeatureSink(
                          'OUTPUT',
                          self.tr('Reduced bend')))

    def processAlgorithm(self, parameters, context, feedback):
        """Main method that extract parameters and call ReduceBend algorithm.
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        # Extract parameter
        source_in = self.parameterAsSource(parameters, "INPUT", context)
        diameter_tol = self.parameterAsDouble(parameters, "TOLERANCE", context)
        smooth_line = self.parameterAsDouble(parameters, "SMOOTH", context)
        exclude_hole = self.parameterAsBool(parameters, "EXCLUDE_HOLE", context)
        exclude_polygon = self.parameterAsBool(parameters, "EXCLUDE_POLYGON", context)
        validate_structure = self.parameterAsBool(parameters, "VALIDATE_STRUCTURE", context)
        verbose = self.parameterAsBool(parameters, "VERBOSE", context)

        if source_in is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        # Transform the in source into a vector layer
        vector_layer_in = source_in.materialize(QgsFeatureRequest(), feedback)

        # Normalize and extract QGS input features
        qgs_features_in, geom_type = ReduceBend.normalize_in_vector_layer(vector_layer_in, feedback)

        # Validate input geometry type
        if geom_type not in (QgsWkbTypes.LineString, QgsWkbTypes.Polygon):
            raise QgsProcessingException("Can only process: (Multi)LineString or (Multi)Polygon vector layers")

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                               vector_layer_in.fields(),
                                               geom_type,
                                               vector_layer_in.sourceCrs())

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))

        # Set progress bar to 1%
        feedback.setProgress(1)

        # Call ReduceBend algorithm
        rb_return = ReduceBend.reduce(qgs_features_in, diameter_tol, smooth_line, exclude_polygon,
                                      exclude_hole, validate_structure, feedback)

        for qgs_feature_out in rb_return.qgs_features_out:
            sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

        # Push some output statistics
        feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
        feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
        feedback.pushInfo("Number of iteration needed: {0}".format(rb_return.nbr_pass))
        feedback.pushInfo("Number of bends detected: {0}".format(rb_return.nbr_bend_detected[0]))
        feedback.pushInfo("Number of bends reduced: {0}".format(sum(rb_return.nbr_bend_reduced)))
        feedback.pushInfo("Number of deleted polygons: {0}".format(rb_return.nbr_pol_del))
        feedback.pushInfo("Number of deleted polygon holes: {0}".format(rb_return.nbr_hole_del))
        feedback.pushInfo("Number of line smoothed: {0}".format(rb_return.nbr_line_smooth))
        if validate_structure:
            if rb_return.is_structure_valid:
                status = "Valid"
            else:
                status = "Invalid"
            feedback.pushInfo("Debug - State of the internal data structure: {0}".format(status))
        if verbose:
            for i in range(rb_return.nbr_pass):
                str_value = "Iteration: {}; Bends detected: {}; Bend reduced: {}" \
                            .format(i, rb_return.nbr_bend_detected[i], rb_return.nbr_bend_reduced[i])
                feedback.pushInfo("Verbose - {0}".format(str_value))

        return {"OUTPUT": dest_id}


# --------------------------------------------------------
# Start of the algorithm
# --------------------------------------------------------

# Define global constant
ANTI_CLOCK_WISE = -1
CLOCK_WISE = 0


class RbFeature(ABC):
    """Contain one QgsFeature

    Abstract class specialized into processing specific geometries
    """

    _id_counter = 0  # Counter of feature

    @staticmethod
    def is_point(feature_type):
        """Static method which determine if a QgsFeature is any kind of Point.

        :param feature_type: Feature type to validate.
        :return: True if a point False otherwise
        :rtype: bool
        """

        val = feature_type in [QgsWkbTypes.Point, QgsWkbTypes.Point25D, QgsWkbTypes.PointM, QgsWkbTypes.PointZ,
                               QgsWkbTypes.PointZM]

        return val

    @staticmethod
    def is_line_string(feature_type):
        """Static method which determine if a QgsFeature is any kind of LineString.

        :param feature_type: Feature type to validate.
        :return: True if it's a LineString; False otherwise
        :rtype: bool
        """

        val = feature_type in [QgsWkbTypes.LineString, QgsWkbTypes.LineString25D, QgsWkbTypes.LineStringZ,
                               QgsWkbTypes.LineStringM, QgsWkbTypes.LineStringZM]

        return val

    @staticmethod
    def is_polygon(feature_type):
        """Static method which determine if a QgsFeature is any kind of Polygon.

        :param feature_type: Feature type to validate.
        :return: True if a Polygon False otherwise
        :rtype: bool
        """
        val = feature_type in [QgsWkbTypes.Polygon, QgsWkbTypes.Polygon25D, QgsWkbTypes.PolygonZ, QgsWkbTypes.PolygonM,
                               QgsWkbTypes.PolygonZM]

        return val

    def __init__(self, qgs_feature):
        """Constructor of the RbFeature class.

        :param qgs_feature: QgsFeature to process.
        """

        self.qgs_feature = qgs_feature
        self.id = RbFeature._id_counter
        RbFeature._id_counter += 1
        abs_geom = qgs_feature.geometry().constGet()
        self.qgs_geom = QgsGeometry(abs_geom.clone())
        self.qgs_feature.clearGeometry()  # Empty the geometry.  Geometry to be recreated at the end

    @abstractmethod
    def get_rb_geom(self):
        """Define an abstract method.
        """

    @abstractmethod
    def get_qgs_feature(self):
        """Define an abstract method.
        """


class RbPolygon(RbFeature):
    """Class description for RbPolygon"""

    def __init__(self, qgs_feature):
        """Constructor that breaks the Polygon into a list of closed LineString (RbGeom).

        :param qgs_feature: QgsFeature polygon to process.
        """

        super().__init__(qgs_feature)
        if self.qgs_geom.wkbType() != QgsWkbTypes.Polygon:
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.Polygon)  # Force geometry to be a QgsPolygon
        # Transform geometry into a list a LineString first ring being outer ring
        self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.LineString)
        # Breaks the rings into a list of closed RbGeom (LineString). The first one being the outer ring
        self.rb_geom = [RbGeom(qgs_geom, QgsWkbTypes.Polygon) for qgs_geom in self.qgs_geom]
        self.qgs_geom = None

    def get_rb_geom(self):
        """Return the RbGeom.

        :return: The RbGeom of the instance
        :rtype: List of RbGeom
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the new geometry.

        :return: The new QgsFeature
        :rtype: QgsFeature
        """

        qgs_pol = QgsPolygon()
        qgs_pol.setExteriorRing(self.rb_geom[0].qgs_geom.constGet().clone())
        for rb_geom in self.rb_geom[1:]:
            qgs_pol.addInteriorRing(rb_geom.qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_pol)

        return self.qgs_feature


class RbLineString(RbFeature):
    """Class managing a RbLineString.
    """

    def __init__(self, qgs_feature):
        """Constructor that breaks the LineString into a list of LineString (RbGeom).

        :param qgs_feature: QgsFeature LineString to process.
        """
        super().__init__(qgs_feature)
        if self.qgs_geom.wkbType() != QgsWkbTypes.LineString:
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.LineString)  # Force geometry to a QgsPoint
        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.LineString)]
        self.qgs_geom = None

    def get_rb_geom(self):
        """Return the RbGeom.

        :return: The RbGeom of the instance
        :rtype: List of RbGeom
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the new geometry.

        :return: The new QgsFeature
        :rtype: QgsFeature
        """

        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbPoint(RbFeature):
    """Class managing a RbPoint
    """

    def __init__(self, qgs_feature):
        """Constructor that breaks the Point into a list of Point (RbGeom).

        :param: qgs_feature: QgsFeature Point to process.
        """

        super().__init__(qgs_feature)
        if self.qgs_geom.wkbType() != QgsWkbTypes.Point:
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.Point)  # Force geometry to QgsPoint
        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.Point)]
        self.rb_geom[0].is_simplest = True  # A point cannot be reduced
        self.qgs_geom = None

    def get_rb_geom(self):
        """Return the RbGeom.

        :return: The RbGeom of the instance.
        :rtype: List of RbGeom.
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the original geometry.

        A Point cannot be reduced but is needed for the spatial constraints validation

        :return: The new QgsFeature
        :rtype: QgsFeature
        """

        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbCollection:
    """Class used for managing the feature spatially.

    QgsSpatialIndex class is used to store and retrieve the features.
    """

    __slots__ = ('_spatial_index', '_dict_qgs_rb_geom', '_dict_qgs_segment', '_id_qgs_segment')

    def __init__(self, rb_results):
        """Constructor that initialize the RbCollection.

        """

        self._spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
        self._dict_qgs_segment = {}  # Contains a reference to the original geometry
        self._id_qgs_segment = 0

    def _get_next_id_segment(self):
        """Increment the id of the segment.

        :return: Value of the next ID
        :rtype: int
        """

        self._id_qgs_segment += 1

        return self._id_qgs_segment

    def _create_feature_segment(self, origin_id, qgs_geom):
        """Creates a new QgsFeature to load in the QgsSpatialIndex.

        :param int origin_id: ID of the source feature (for reference purpose)
        :return: The feature created
        :rtype: QgsFeature
        """

        id_segment = self._get_next_id_segment()
        self._dict_qgs_segment[id_segment] = origin_id  # Creates a reference for to the original RbGeom
        qgs_feature = QgsFeature(id=id_segment)
        qgs_feature.setGeometry(qgs_geom)

        return qgs_feature

    def add_features(self, rb_geoms):
        """Add a RbGeom object in the spatial index.

        For the LineString geometries. The geometry is broken into each line segment that are individually
        loaded in the QgsSpatialIndex.  This strategy takes longer to load than if the feature was loaded as a whole
        but is better for much of the cases.

        :param [RbGeom] rb_geoms: List of RbGeom to load
        """

        for rb_geom in rb_geoms:
            qgs_features = []
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.Point:
                qgs_features.append(self._create_feature_segment(rb_geom.id, rb_geom.qgs_geom))
            else:
                qgs_points = rb_geom.qgs_geom.constGet().points()
                for i in range(0, (len(qgs_points)-1)):
                    qgs_geom = QgsGeometry(QgsLineString(qgs_points[i], qgs_points[i+1]))
                    qgs_feature = self._create_feature_segment(rb_geom.id, qgs_geom)
                    qgs_features.append(qgs_feature)

            self._spatial_index.addFeatures(qgs_features)  # Load all the segment of one RbGeom at the same time

        return

    def get_segment_intersect(self, qgs_geom_id, qgs_rectangle, qgs_geom_subline):
        """Find the feature that intersects the bounding box.

        Once the line string intersecting the bounding box are found. They are separated into 2 lists.
        The first one being the line string with the same id (same line) the second one all the others line string.

        :param qgs_geom_id: ID of the line string that is being simplified
        :param qgs_rectangle: QgsRectangle used for feature intersection
        :param qgs_geom_subline: LineString used to remove line segment superimposed to this line string
        :return: Two lists  of line string segment. First: Line string with same id; Second all the others
        :rtype: tuple of 2 lists
        """

        qgs_geoms_with_itself = []
        qgs_geoms_with_others = []
        qgs_rectangle.grow(Epsilon.ZERO_RELATIVE*100.)  # Always increase the b_box to avoid degenerated b_box
        qgs_segment_ids = self._spatial_index.intersects(qgs_rectangle)
        for qgs_segment_id in qgs_segment_ids:
            qgs_geom_segment = self._spatial_index.geometry(qgs_segment_id)
            if qgs_geom_segment.wkbType() == QgsWkbTypes.Point:
                qgs_geoms_with_others.append(qgs_geom_segment)
            else:
                if self._dict_qgs_segment[qgs_segment_id] == qgs_geom_id:
                    if not qgs_geom_segment.within(qgs_geom_subline):
                        qgs_geoms_with_itself.append(qgs_geom_segment)
                else:
                    qgs_geoms_with_others.append(qgs_geom_segment)

        return qgs_geoms_with_itself, qgs_geoms_with_others

    def _delete_segment(self, qgs_pnt0, qgs_pnt1):
        """Delete a line segment in the spatial index based on start/end points.

        To minimise the number of feature returned we search for a very small bounding box located in the middle
        of the line segment.  Usually only one line segment is returned.

        :param qgs_pnt0 : QgsPoint start point of the target line segment.
        :param qgs_pnt1 : QgsPoint end point of the target line segment.
        """

        qgs_geom_target = QgsGeometry(QgsLineString(qgs_pnt0, qgs_pnt1))
        qgs_mid_point = QgsGeometryUtils.midpoint(qgs_pnt0, qgs_pnt1)
        qgs_rectangle = qgs_mid_point.boundingBox()
        qgs_rectangle.grow(Epsilon.ZERO_RELATIVE*100)
        feat_ids = self._spatial_index.intersects(qgs_rectangle)
        deleted = True
        for feat_id in feat_ids:
            qgs_geom_line = self._spatial_index.geometry(feat_id)  # Extract geometry
            #  Check if it's the target geometry
            if qgs_geom_target.hausdorffDistance(qgs_geom_line) <= Epsilon.ZERO_RELATIVE:
                feature = QgsFeature(id=feat_id)
                feature.setGeometry(qgs_geom_target)
                if self._spatial_index.deleteFeature(feature):  # Delete the line segment
                    deleted = True
                    break
                else:
                    raise Exception(QgsProcessingException("Unable to delete entry in QgsSpatialIndex..."))
            else:
                deleted = False

        if not deleted:
            raise Exception(QgsProcessingException("Internal structure corruption..."))

        return

    def delete_vertex(self, rb_geom, v_id_start, v_id_end):
        """Delete a vertex in the line and update the spatial index.

        When a vertex in a line string is deleted.  Two line segments are deleted and one line segment is
        created in the spatial index.  Cannot delete the first/last vertex of a line string

        :param rb_geom: LineString object to update.
        :param v_id_start: start of the vertex to delete.
        :param v_id_end: end of the vertex to delete.
        """

        is_closed = rb_geom.qgs_geom.constGet().isClosed()
        v_ids_to_del = list(range(v_id_start, v_id_end+1))
        if is_closed and v_id_start == 0:
            # Special case for closed line where we simulate a circular array
            nbr_vertice = rb_geom.qgs_geom.constGet().numPoints()
            v_ids_to_del.insert(0, nbr_vertice - 2)
        else:
            v_ids_to_del.insert(0, v_ids_to_del[0]-1)
        v_ids_to_del.append(v_ids_to_del[-1]+1)

        # Delete the line segment in the spatial index
        for i in range(len(v_ids_to_del)-1):
            qgs_pnt0 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[i])
            qgs_pnt1 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[i+1])
            self._delete_segment(qgs_pnt0, qgs_pnt1)

        # Add the new line segment in the spatial index
        qgs_pnt0 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[0])
        qgs_pnt1 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[-1])
        qgs_geom_segment = QgsGeometry(QgsLineString(qgs_pnt0, qgs_pnt1))
        qgs_feature = self._create_feature_segment(rb_geom.id, qgs_geom_segment)
        self._spatial_index.addFeature(qgs_feature)

        # Delete the vertex in the line string geometry
        for v_id_to_del in reversed(range(v_id_start, v_id_end+1)):
            rb_geom.qgs_geom.deleteVertex(v_id_to_del)
            if is_closed and v_id_start == 0:
                # Special case for closed line where we simulate a circular array
                nbr_vertice = rb_geom.qgs_geom.constGet().numPoints()
                qgs_pnt_first = rb_geom.qgs_geom.vertexAt(0)
                rb_geom.qgs_geom.insertVertex(qgs_pnt_first, nbr_vertice-1)
                rb_geom.qgs_geom.deleteVertex(nbr_vertice)

        return

    def add_vertex(self, rb_geom, bend_i, bend_j, qgs_geom_new_subline):
        """Update the line segment in the spatial index

        :param rb_geom: RbGeom line to update
        :param bend_i: Start of the bend to delete
        :param bend_j: End of the bend to delete (always bend_i + 1)
        :param qgs_geom_new_subline: New sub line string to add in the spatial index
        :return:
        """

        # Delete the base of the bend
        qgs_pnt0 = rb_geom.qgs_geom.vertexAt(bend_i)
        qgs_pnt1 = rb_geom.qgs_geom.vertexAt(bend_j)
        self._delete_segment(qgs_pnt0, qgs_pnt1)

        qgs_points = qgs_geom_new_subline.constGet().points()
        tmp_qgs_points = qgs_points[1:-1]  # Drop first/last item
        # Insert the new vertex in the QgsGeometry. Work reversely to facilitate insertion
        for qgs_point in reversed(tmp_qgs_points):
            rb_geom.qgs_geom.insertVertex(qgs_point, bend_j)

        # Add the new segment in the spatial container
        for i in range(len(qgs_points)-1):
            qgs_pnt_a = qgs_points[i]
            qgs_pnt_b = qgs_points[i+1]
            qgs_geom_segment = QgsGeometry(QgsLineString(qgs_pnt_a, qgs_pnt_b))
            qgs_feature = self._create_feature_segment(rb_geom.id, qgs_geom_segment)
            self._spatial_index.addFeature(qgs_feature)

        return

    def validate_integrity(self, rb_geoms):
        """This method is used to validate the data structure at the end of the process

        This method is executed only when requested and for debug purpose only.  It's validating the data structure
        by removing element from it the data structure is unusable after. Validate integrity must be the last
        operation before ending the program as it destroy the data structure...

        :param rb_geoms: Geometry contained in the spatial container
        :return: Flag indicating if the structure is valid. True: is valid; False: is not valid
        :rtype: Boolean
        """

        is_structure_valid = True
        # from the geometry remove all the segment in the spatial index.
        for rb_geom in rb_geoms:
            qgs_line_string = rb_geom.qgs_geom.constGet()
            if qgs_line_string.wkbType() == QgsWkbTypes.LineString:
                qgs_points = qgs_line_string.points()
                for i in range(len(qgs_points)-1):
                    self._delete_segment(qgs_points[i], qgs_points[i+1])

        if is_structure_valid:
            # Verify that there are no other feature in the spatial index; except for QgsPoint
            qgs_rectangle = QgsRectangle(-sys.float_info.max, -sys.float_info.max,
                                         sys.float_info.max, sys.float_info.max)
            feat_ids = self._spatial_index.intersects(qgs_rectangle)
            for feat_id in feat_ids:
                qgs_geom = self._spatial_index.geometry(feat_id)
                if qgs_geom.wkbType() == QgsWkbTypes.Point:
                    pass
                else:
                    # Error
                    is_structure_valid = False

        return is_structure_valid


class RbGeom:
    """Class defining the line string used for the bend reduction"""

    __slots__ = ('id', 'original_geom_type', 'is_simplest', 'qgs_geom', 'qgs_rectangle', 'bends', 'nbr_bend_reduced',
                 'need_pivot')

    _id_counter = 0  # Unique ID counter

    @staticmethod
    def next_id():
        """Get the next counterID.

        :param: QgsMultiLineString qgs_multi_line_string: Multi line string to merge together
        :return: ID of the RbGeom object
        :rtype: int
        """

        RbGeom._id_counter += 1

        return RbGeom._id_counter

    def __init__(self, qgs_abs_geom, original_geom_type):
        """Constructor that initialize a RbGeom object.

        :param: qgs_abs_geom: QgsAbstractGeometry to process
        :param: original_geom_type: Original type of the geometry

        """

        self.id = RbGeom.next_id()
        self.original_geom_type = original_geom_type
        qgs_geometry = qgs_abs_geom.constGet()
        self.qgs_geom = QgsGeometry(qgs_geometry.clone())
        self.is_simplest = False
        self.need_pivot = False
        self.bends = []
        self.nbr_bend_reduced = 0
        # Set some variable depending on the attribute of the feature
        if self.original_geom_type == QgsWkbTypes.Point:
            self.is_simplest = True  # A point cannot be simplified
        elif self.original_geom_type == QgsWkbTypes.LineString:
            if qgs_geometry.length() >= Epsilon.ZERO_RELATIVE:
                if qgs_geometry.isClosed():  # Closed LineString
                    if abs(qgs_geometry.sumUpArea()) > Epsilon.ZERO_RELATIVE:
                        self.need_pivot = True
                    else:
                        self.is_simplest = True  # Zero area polygon (degenerated).  Do not try to simplify
            else:
                self.is_simplest = True  # Zero length line (degenerated). Do not try to simplify
        elif self.original_geom_type == QgsWkbTypes.Polygon:
            qgs_polygon = QgsPolygon(qgs_geometry.clone())  # Create QgsPolygon to calculate area
            if qgs_polygon.area() > Epsilon.ZERO_RELATIVE:
                self.need_pivot = True
            else:
                self.is_simplest = True  # Zero area polygon. Do not simplify the closed line


class Bend:
    """Define a Bend object which is the reduction goal of this algorithm"""

    __slots__ = ('i', 'j', 'area', 'perimeter', 'adj_area', 'to_reduce', '_qgs_geom_new_subline',
                 '_qgs_geom_old_subline', '_qgs_points', 'qgs_geom_bend')

    def __init__(self, i, j, qgs_points):
        """Constructor that initialize a Bend object.

        :param: int i: start position of the vertice in the LineString to reduce
        :param: int j: end position of the vertice in the LineString to reduce
        :param: qgs_points: List of QgsPoint defining the bend
        :return: None
        :rtype: None
        """

        self.i = i
        self.j = j
        self._qgs_points = qgs_points
        self._qgs_geom_new_subline = None
        self._qgs_geom_old_subline = None
        self.qgs_geom_bend = QgsGeometry(QgsPolygon(QgsLineString(qgs_points)))  # QgsPolygon will close the polygon
        self.area = self.qgs_geom_bend.area()
        self.perimeter = self.qgs_geom_bend.length()
        self.adj_area = ReduceBend.calculate_adj_area(self.area, self.perimeter)
        self.to_reduce = False

    @property
    def qgs_geom_new_subline(self):
        """Late attribute evaluation as this attribute is costly to evaluate"""
        if self._qgs_geom_new_subline is None:
            self._qgs_geom_new_subline = QgsGeometry(QgsLineString(self._qgs_points[0], self._qgs_points[-1]))
        return self._qgs_geom_new_subline

    @property
    def qgs_geom_old_subline(self):
        """Late attribute evaluation as this attribute is costly to evaluate"""
        if self._qgs_geom_old_subline is None:
            self._qgs_geom_old_subline = QgsGeometry(QgsLineString(self._qgs_points))
        return self._qgs_geom_old_subline


class BendReduced:
    """Class containing the information of a bend that has been reduced. """

    CASE_1 = "Case1"
    CASE_2 = "Case2"
    CASE_3 = "Case3"

    @staticmethod
    def _calculate_angle(angle_i, angle_j, smooth_case):
        """

        :param: angle_j: Angle in degrees at vertex i
        :param: angle_j: Angle in degrees at vertex j
        :param: smooth_case: Type of bend smoothing value between [1..3]
        :return:
        """
        if angle_i > math.pi:
            angle_i = (2 * math.pi) - angle_i  # Normalize angle between [0..180]
        if angle_j > math.pi:
            angle_j = (2 * math.pi) - angle_j  # Normalize angle between [0..180]
        angle_smooth = max(angle_i, angle_j)
        angle_smooth = (math.pi - angle_smooth)
        if smooth_case == BendReduced.CASE_1:
            angle_smooth /= 1.5
            if math.degrees(angle_smooth) > 30.:
                angle_smooth = math.radians(30.)
        elif smooth_case == BendReduced.CASE_2:
            angle_smooth /= 2.5
            if math.degrees(angle_smooth) > 20.:
                angle_smooth = math.radians(20.)
        elif smooth_case == BendReduced.CASE_3:
            angle_smooth /= 3
            if math.degrees(angle_smooth) > 20.:
                angle_smooth = math.radians(20.)

        return angle_smooth

    __slots__ = ('rb_geom', 'qgs_point_start', 'qgs_point_end', 'qgs_geom_bend', 'qgs_geom_old_subline',
                 'i', 'j', 'is_line_smoothable', 'qgs_geom_smooth_line', 'qgs_geom_smooth_polygon')

    def __init__(self, rb_geom, qgs_point_start, qgs_point_end, qgs_geom_bend):
        """Constructor that initialize a BendReduced object

        :param: rb_geom: RbGeom object containing the reduced bend
        :param: qgs_point_start: QgsPoint of the start of the bend reduced
        :param: qgs_point_end: QgsPoint of the end of the bend reduced
        :param: qgs_geom_bend: QgsPolygon representing the bend reduced
        """

        self.rb_geom = rb_geom
        self.qgs_point_start = qgs_point_start
        self.qgs_point_end = qgs_point_end
        self.qgs_geom_bend = qgs_geom_bend
        self.qgs_geom_old_subline = QgsGeometry(QgsLineString([qgs_point_start, qgs_point_end]))
        self.qgs_geom_smooth_line = None
        self.qgs_geom_smooth_polygon = None
        self.i = None
        self.j = None
        self.is_line_smoothable = None

    def _resolve_non_valid_polygon(self, epsilon):
        """This method split a self intersecting line (non valid line) into a non intersecting multipolygon

        The method closes the input line string and if it is self intersecting. It uses unary union and
        polygonize to split the line string into a list of non self intersecting polygon (multipolygon)

        :param: epsilon: Float value representing near zero value
        :return: Non self intersecting multipolygon
        :rtype: QgsGeometry
        """

        qgs_multi_pol = QgsMultiPolygon()
        qgs_line_string = self.qgs_geom_smooth_line.constGet().clone()
        qgs_geom_pol = QgsGeometry(QgsPolygon(qgs_line_string.clone()))
        if qgs_geom_pol.isGeosValid():
            # Polygon is valid nothing to do
            qgs_multi_pol.addGeometry(qgs_geom_pol.constGet().clone())
        else:
            # Polygon is invalid... let's try to repair it
            qgs_pnt_start = qgs_line_string.startPoint()
            qgs_line_string.addVertex(qgs_pnt_start)  # Close the line with the startpoint
            qgs_geom_close_line = QgsGeometry(qgs_line_string)
            qgs_geom_unary = QgsGeometry.unaryUnion([qgs_geom_close_line])  # Create node at each overlap
            qgs_geom_polygonize = QgsGeometry.polygonize([qgs_geom_unary])  # Create multi polygon
            for qgs_part in qgs_geom_polygonize.parts():
                # Validate that each part is valid
                qgs_geom_part = QgsGeometry(qgs_part.clone())
                if qgs_geom_part.isGeosValid() and qgs_geom_part.area() > epsilon:
                    qgs_multi_pol.addGeometry(qgs_geom_part.constGet().clone())
                else:
                    self.is_line_smoothable = False  # Something went wrong... nothing else to try... quit
                    break

        self.qgs_geom_smooth_polygon = QgsGeometry(qgs_multi_pol.clone())

        return

    def _extract_vertex_ind(self, qgs_ls, qgs_point):
        """Extract the vertice number of a QgsPoint in a QgsLineString

        :param qgs_ls: Input QgsLineString
        :param qgs_point: QgsPoint for which we want the vertice number in the line
        :return: Vertice number of the point in the QgsLineString
        :rtype: int
        """

        vertex_info = QgsGeometryUtils.closestVertex(qgs_ls, qgs_point)
        vertex_id = vertex_info[1].vertex
        qgs_point_target = self.rb_geom.qgs_geom.vertexAt(vertex_id)
        if qgs_point_target.distance(qgs_point) < Epsilon.ZERO_RELATIVE:
            ind = vertex_id
        else:
            ind = None

        return ind

    def _calculate_smooth_line(self):
        """This method will transform a straight subline into a "smooth" sub line

        To ease the smoothing of the sub line, the subline is translate so that the start of the bend (bend.i)
        is at (0,0).  The sub line is than rotate to align on the x axis. If the sub line meet the requirements
        2 points are added on the line. The line is than rotated and translated to there original position.
        """

        qgs_geom_bend_centroid = self.qgs_geom_bend.centroid()  # Centroid of the bend
        qgs_points_subline = []
        for ind in [self.i-1, self.i, self.j, self.j+1]:
            qgs_points_subline.append(self.rb_geom.qgs_geom.vertexAt(ind))  # Information needed to smooth the line

        qgs_point_translate = qgs_points_subline[1].clone()  # Bend i is the point used for translation and rotation
        qgs_geom_smooth = QgsGeometry(QgsLineString(qgs_points_subline).clone())
        qgs_geom_smooth.translate(-qgs_point_translate.x(), -qgs_point_translate.y())
        qgs_geom_bend_centroid.translate(-qgs_point_translate.x(), -qgs_point_translate.y())

        qgs_points_tr = qgs_geom_smooth.constGet().points()
        x_axis_length = QgsLineString([qgs_points_tr[1], qgs_points_tr[2]]).length()
        qgs_point_x_axis = QgsPoint(x_axis_length, 0.)
        # Calculate the angle between the translated sub line and the x axis
        angle_x_axis = QgsGeometryUtils.angleBetweenThreePoints(qgs_point_x_axis.x(), 0.,
                                                                qgs_points_tr[1].x(), qgs_points_tr[1].y(),
                                                                qgs_points_tr[2].x(), qgs_points_tr[2].y())
        angle_x_axis_degree = math.degrees(angle_x_axis)
        qgs_geom_smooth.rotate(-angle_x_axis_degree, QgsPointXY(0, 0))  # Rotate the subline on the x axis
        qgs_geom_bend_centroid.rotate(-angle_x_axis_degree, QgsPointXY(0, 0))  # Rotate the centroid of the bend
        qgs_points_ro = qgs_geom_smooth.constGet().points()

        base_length = qgs_points_ro[2].x()  # The length of the base of the bend on the x axis
        p0_x = base_length * (1. / 3.)  # x position of the first point
        p1_x = base_length * (2. / 3.)  # x position of the second point

        # Set the different smoothing cases of bend to smooth
        if qgs_points_ro[0].y() * qgs_points_ro[3].y() > 0:
            if qgs_points_ro[0].y() * qgs_geom_bend_centroid.constGet().y() < 0:
                # The bend and the segment previous and after the bend are located on the opposite side of the x axis
                smooth_case = BendReduced.CASE_1
            else:
                # The bend and the segment previous and after the bend are located on the same side of the x axis
                smooth_case = BendReduced.CASE_2
        else:
            # The previous and after the bend are located on different side of the x axis
            smooth_case = BendReduced.CASE_3

        angle_i = QgsGeometryUtils.angleBetweenThreePoints(qgs_points_ro[0].x(), qgs_points_ro[0].y(),
                                                           qgs_points_ro[1].x(), qgs_points_ro[1].y(),
                                                           qgs_points_ro[2].x(), qgs_points_ro[2].y(),)
        angle_j = QgsGeometryUtils.angleBetweenThreePoints(qgs_points_ro[1].x(), qgs_points_ro[1].y(),
                                                           qgs_points_ro[2].x(), qgs_points_ro[2].y(),
                                                           qgs_points_ro[3].x(), qgs_points_ro[3].y(),)

        angle_smooth = BendReduced._calculate_angle(angle_i, angle_j, smooth_case)

        p0_y = math.tan(angle_smooth) * p0_x  # Trigonometric formula to find length of opposite side (value of y)

        if smooth_case in [BendReduced.CASE_1, BendReduced.CASE_2]:
            if qgs_points_ro[0].y() > 0.:
                p0_y *= -1
            qgs_point_smooth_0 = QgsPoint(p0_x, p0_y)
            qgs_point_smooth_1 = QgsPoint(p1_x, p0_y)
        else:  # CASE_3
            if qgs_points_ro[0].y() > 0.:
                p0_y *= -1
            qgs_point_smooth_0 = QgsPoint(p0_x, p0_y)
            p0_y *= -1
            qgs_point_smooth_1 = QgsPoint(p1_x, p0_y)

        qgs_geom_smooth = QgsGeometry(QgsLineString([qgs_points_ro[1], qgs_point_smooth_0,
                                                     qgs_point_smooth_1, qgs_points_ro[2]]).clone())
        qgs_geom_smooth.rotate(angle_x_axis_degree, QgsPointXY(0, 0))
        qgs_geom_smooth.translate(qgs_point_translate.x(), qgs_point_translate.y())

        self.qgs_geom_smooth_line = qgs_geom_smooth

        return

    def set_values(self, diameter_tol):
        """Set different values required for the smoothing od the line

        This method also validate if this reduced bend is a candidate for line smoothing

        :param diameter_tol: Diameter tolerance used for the bend reduction
        """

        if self.qgs_geom_old_subline.constGet().length() > diameter_tol * (2. / 3.):  # Do not smooth short bend base
            qgs_line_string = self.rb_geom.qgs_geom.constGet()
            self.i = self._extract_vertex_ind(qgs_line_string, self.qgs_point_start)
            self.j = self._extract_vertex_ind(qgs_line_string, self.qgs_point_end)

            if self.i is not None and self.j is not None:
                if self.i + 1 == self.j:
                    if self.i >= 1 and self.j <= self.rb_geom.qgs_geom.constGet().numPoints() - 2:
                        self.is_line_smoothable = True  # Good candidate for line smoothing
                    else:
                        self.is_line_smoothable = False  # QgsPoint must not be first or last point in the line
                else:
                    self.is_line_smoothable = False  # The QgsPoint must be consecutive
            else:
                self.is_line_smoothable = False  # One or both QgsPoint are not there any more

        if self.is_line_smoothable:
            self._calculate_smooth_line()
            self._resolve_non_valid_polygon(Epsilon.ZERO_RELATIVE)

        return


class RbResults:
    """Class defining the stats and result"""

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'nbr_bend_reduced', 'nbr_bend_detected',
                 'qgs_features_out', 'nbr_hole_del', 'nbr_pol_del', 'nbr_pass', 'is_structure_valid',
                 'nbr_line_smooth')

    def __init__(self):
        """Constructor that initialize a RbResult object.

        :param: None
        :return: None
        :rtype: None
        """

        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_bend_reduced = []  # One value per iteration
        self.nbr_bend_detected = []  # One value per iteration
        self.qgs_features_out = None
        self.nbr_hole_del = 0
        self.nbr_pol_del = 0
        self.nbr_pass = 0
        self.nbr_line_smooth = 0
        self.is_structure_valid = None


class Epsilon:
    """Class defining the value of the zero"""

    ZERO_RELATIVE = None
    ZERO_ABSOLUTE = None
    ZERO_ANGLE = None

    __slots__ = '_zero_relative', '_zero_absolute', '_zero_angle', '_map_range'

    def __init__(self, features):
        """Constructor that initialize the Epsilon (near zero) object.

        The dynamic (range) of the feature can vary a lot. We calculate the dynamic of the bounding box of all
        the features and we use it to estimate an epsilon (zero).  when the range of the bounding box is very small
        the epsilon can be very small and the opposite when the bigger the bounding box is.

        :param: [QgsFeatures] features: List of QgsFeature to process.
        :return: None
        :rtype: None
        """

        if len(features) >= 1:
            b_box = features[0].geometry().boundingBox()  # Initialize the bounding box
        else:
            b_box = QgsRectangle(0, 0, 1, 1)  # Manage empty list of feature

        for feature in features:
            b_box.combineExtentWith(feature.geometry().boundingBox())  # Update the bbox

        delta_x = abs(b_box.xMinimum()) + abs(b_box.xMaximum())
        delta_y = abs(b_box.yMinimum()) + abs(b_box.yMaximum())
        dynamic_xy = max(delta_x, delta_y)  # Dynamic of the bounding box
        log_loss = int(math.log(dynamic_xy, 10)+1)
        max_digit = 15  # Number of significative digits for real number
        security = 2  # Keep 2 order of magnitude of security
        abs_digit = max_digit - security
        rel_digit = max_digit - log_loss - security
        self._zero_relative = (1. / (10**rel_digit))
        self._zero_absolute = (1. / (10**abs_digit))
        self._zero_angle = math.radians(.0001)  # Angle used to decide a flat angle

    def set_class_variables(self):
        """Set the different epsilon values.

        :return: None
        :rtype: None
        """

        Epsilon.ZERO_RELATIVE = self._zero_relative
        Epsilon.ZERO_ABSOLUTE = self._zero_absolute
        Epsilon.ZERO_ANGLE = self._zero_angle

        return


class ReduceBend:
    """Main class for bend reduction"""

    @staticmethod
    def normalize_in_vector_layer(in_vector_layer, feedback):
        """Method used to normalize the input vector layer

        Two processing are used to normalized the input vector layer
         - execute "Multi to single part" processing in order to accept even multipolygon
         - execute "Drop  Z and M values" processing as they are not useful
         - Validate if the resulting layer is Point LineString or Polygon

        :param in_vector_layer:  Input vector layer to normalize
        :param feedback: QgsFeedback handle used to communicate with QGIS
        :return Output vector layer and Output geometry type
        :rtype Tuple of 2 values
        """

        # Execute MultiToSinglePart processing
        feedback.pushInfo("Start normalizing input layer")
        params = {'INPUT': in_vector_layer,
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

        # Extract the QgsFeature from the vector layer
        qgs_in_features = []
        qgs_features = drop_zm_layer.getFeatures()
        for qgs_feature in qgs_features:
            qgs_in_features.append(qgs_feature)
        if len(qgs_in_features) > 1:
            geom_type = qgs_in_features[0].geometry().wkbType()
        else:
            geom_type = drop_zm_layer.wkbType()  # In case of empty layer
        feedback.pushInfo("End normalizing input layer")

        return qgs_in_features, geom_type

    @staticmethod
    def get_angles(qgs_line_string):
        """Extract the list of angles of a LineString

        The angle of a vertice is the angle formed by the preceding, the current and the next vertice.
        For an open LineString the start/end vertice do not have angles

        :param: qgs_line_string: QgsLineString to extract angles
        :return: Angle of each vertice of the LineString
        :rtype: [real]
        """

        qgs_pnts = qgs_line_string.points()
        xy = [(qgs_pnt.x(), qgs_pnt.y()) for qgs_pnt in qgs_pnts]
        if len(xy) >= 3:
            if qgs_line_string.isClosed():
                # Add two vertice at the start/end for the circularity of a closed line
                end_xy = xy[-2]  # Not the last vertice because it's the same position as the first vertice
                xy.insert(0, end_xy)

            angles = [QgsGeometryUtils.angleBetweenThreePoints(xy[i-1][0], xy[i-1][1], xy[i][0], xy[i][1],
                                                               xy[i+1][0], xy[i+1][1]) for i in range(1, len(xy)-1)]
        else:
            # A line string must have 3 vertice in order to calculate an angle
            angles = []

        return angles

    @staticmethod
    def reduce(qgs_in_features, diameter_tol, smooth_line=False, flag_del_outer=False,
               flag_del_inner=False, validate_structure=False, feedback=None):
        """Main static method used to launch the bend reduction.

        :param: [QgsFeatures] qgs_features: List of features to process.
        :param: real diameter_tol: Tolerance of the diameter of the bend to reduce.
        :param: Bool smooth_line: Smooth line after bend reduction if possible
        :param: Bool flag_del_outer: Delete polygon if area below the diameter tolerance.
        :param: Bool flag_del_inner: Delete polygon holes if area below the diameter tolerance.
        :param: Bool validate_structure: Validate internal data structure after processing (for debugging)
        :param: QgsFeedback feedback: Handle for interaction with QGIS.
        :return: Statistics and result object.
        :rtype: RbResult
        """

        rb = ReduceBend(qgs_in_features, diameter_tol, smooth_line, flag_del_outer, flag_del_inner, validate_structure,
                        feedback)
        results = rb.reduce_bends()

        return results

    @staticmethod
    def _extract_polygon_attributes(qgs_geom):
        """Static method to calculate the area and perimeter of a LineString.

       :param: QgsGeometry qgs_geom: Geometry to process.
       :return: Area and perimeter of the geometry.
       :rtype: Tuple
       """

        qgs_line_string = qgs_geom.constGet()
        qgs_pol = QgsPolygon(qgs_line_string.clone())
        area = qgs_pol.area()
        perimeter = qgs_pol.perimeter()

        return area, perimeter

    @staticmethod
    def calculate_adj_area(area, perimeter):
        """Static method to calculate the adjusted area.

        The adjusted area is used to determine if a bend must be reduce.

       :param: real area: area of a polygon.
       :param: real perimeter: perimeter of a polygon.
       :return: Adjusted area of a polygon
       :rtype: Real
       """

        compactness_index = 4 * area * math.pi / perimeter ** 2
        adj_area = area * (.75 / compactness_index)

        return adj_area

    @staticmethod
    def flag_bend_to_reduce(rb_geom, diameter_tol):
        """This method identifies the bend that need to be reduced

        The method starts by finding all the bends that are below the minimum adjusted area. It starts with the
        smallest area bend and if the adjacent area is flag to be simplified, it flags that bend to be simplified;
        it goes after to next smaller bend and so on for each bend that are below the minimum adjusted area

        Note: For closed line it deleted the first and last bend because the line was pivoted on a bend that do not
        need to be reduced

        :param: rb_geom: Geometry to delete co-linear vertices
        :param: diameter_tol: float tolerance for bend reduction
        """

        min_adj_area = ReduceBend.calculate_min_adj_area(diameter_tol)
        if rb_geom.qgs_geom.constGet().isClosed() and len(rb_geom.bends) >= 3:
            # The closed line start/end point lie on a bend that do not need to be reduced
            del rb_geom.bends[0]  # Remove the first bend
            del rb_geom.bends[-1]  # Remove the last bend

        lst_bends = [(bend.adj_area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
        lst_bends.sort(key=lambda item: item[0])

        start = 0
        end = len(rb_geom.bends) - 1

        for (adj_area, i) in lst_bends:
            if adj_area <= min_adj_area:
                if len(lst_bends) == 1:
                    rb_geom.bends[i].to_reduce = True  # Only one bend process it...
                else:
                    if i == start:
                        if rb_geom.bends[i + 1].to_reduce:
                            pass  # Cannot reduce two bend adjacent
                        else:
                            rb_geom.bends[i].to_reduce = True
                    elif i == end:
                        if rb_geom.bends[i - 1].to_reduce:
                            pass  # Cannot reduce two bend adjacent
                        else:
                            rb_geom.bends[i].to_reduce = True
                    elif rb_geom.bends[i - 1].to_reduce or rb_geom.bends[i + 1].to_reduce:
                        pass  # Cannot reduce two bend adjacent
                    else:
                        rb_geom.bends[i].to_reduce = True
            else:
                # Over minimum adjusted area
                break

        if len(rb_geom.bends) == 0:
            # No more bends to reduce
            rb_geom.is_simplest = True

        return

    @staticmethod
    def create_polygon(i, j, qgs_points):
        """This method create a polygon from a subset of point (vertice)

        Note: The first vertice is also the last vertice

        :param: i: Start of point in the list
        :param: j: End of point in the list
        :param: qgs_points: List of QgsPoint
        :return: A polygon formed by a subset of the list of QgsPoint
        :rtype: QgsPolygon
        """

        # Create the list of point to create
        if i < j:
            index = list(range(i, j+1)) + [i]
        else:
            index = list(range(i, len(qgs_points))) + list(range(0, j+1)) + [i]  # Manage circular array

        qgs_sub_points = [qgs_points[k] for k in index]
        qgs_polygon = QgsPolygon(QgsLineString(qgs_sub_points))

        return qgs_polygon

    @staticmethod
    def pivot_closed_line(rb_geom, diameter_tol):
        """For closed LineString, this method will move the start/end vertice.

        The start/end vertice of a closed LineString is moved over a bend that does not need to be simplified.
        By moving the first/start vertice at a position that we know does not need bend reduction, this algorithm
        does not have to deal with the complexity of circular array

        :param: rb_geom: Geometry to set bend direction
        :param: diameter_tol: float tolerance used for determining if a bend must be reduced
        """

        if rb_geom.need_pivot:
            bend_location = None
            bend_area = 0.0
            for bend in rb_geom.bends:
                if bend.area > bend_area:
                    bend_location = bend
                    bend_area = bend.area
                if bend.j - bend.i >= 4:  # Ideal bend for pivot. The bend  has 4 vertices
                    if bend.area >= ReduceBend.calculate_min_adj_area(diameter_tol):
                        bend_location = bend
                        rb_geom.need_pivot = False  # Optimal bend found
                        break

            if bend_location is not None:
                # There is bend candidate for a line rotation
                # Move the start/end of the line in the middle of the bend candidate
                qgs_points = rb_geom.qgs_geom.constGet().points()
                new_start_end = (bend_location.j + bend_location.i) // 2
                new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end + 1]
                rb_geom.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))

        return

    @staticmethod
    def detect_bends(rb_geom):
        """This method detect the bends in a LineString

        In order to detect bends the algorithm is doing two distinct operation:
         - First: for each vertice in the LineString it calculates if the angle is going clockwise or anticlockwise
         - Second: The algorithm is grouping together all the clockwise and all the anticlockwise vertice; when
           there is a change in direction there is a new bend created. A convex figure will only have one bend.

        Process is also slightly different for open and closed LineString

        :param rb_geom: The geometry LineString for which to detect the bends
        :return: Number of bend created on the LineString
        :rtype: int
        """

        rb_geom.bends = []  # Reset the list of bends
        angles = ReduceBend.get_angles(rb_geom.qgs_geom.constGet())
        # Modify the angle to binary orientation: clockwise or anti clockwise
        orientation = [CLOCK_WISE if angle >= math.pi else ANTI_CLOCK_WISE for angle in angles]
        if rb_geom.qgs_geom.constGet().isClosed():
            if len(set(orientation)) == 1:
                orientation = []  # All the angles have the same orientation.  No bend to reduce
            else:
                del orientation[0]  # Do not process the first angle as it is the angle of the start/end

        if len(orientation) >= 1:
            if orientation[0] == CLOCK_WISE:
                orientation.insert(0, ANTI_CLOCK_WISE)
            else:
                orientation.insert(0, CLOCK_WISE)

            if orientation[-1] == CLOCK_WISE:
                orientation.append(ANTI_CLOCK_WISE)
            else:
                orientation.append(CLOCK_WISE)

        # Find the inflexion points in the line.
        inflexion = [i for i in range(0, len(orientation)-1) if orientation[i] != orientation[(i + 1)]]
        qgs_points = rb_geom.qgs_geom.constGet().points()
        if len(inflexion) != 0:
            for k in range(len(inflexion)-1):
                i = inflexion[k]
                j = (inflexion[(k+1)]+1)
                rb_geom.bends.append(Bend(i, j, qgs_points[i:j+1]))

        else:
            # If there is no inflexion the line cannot be simplified
            rb_geom.is_simplest = True

        return len(rb_geom.bends)

    @staticmethod
    def calculate_min_adj_area(diameter_tol):
        """Static method to calculate the adjusted area of the maximum diameter tolerance.

       :param: diameter_tol: float diameter tolerance to used for bend reduction
       :return: Minimum adjusted area of a polygon to reduce
       :rtype: Real
       """

        min_adj_area = .75 * math.pi * (diameter_tol / 2.) ** 2

        return min_adj_area

    @staticmethod
    def validate_simplicity(qgs_geoms_with_itself, qgs_geom_new_subline):
        """Validate the simplictity constraint

        This constraint assure that the new sub line is not intersecting with any other segment of the same line

        :param: qgs_geoms_with_itself: List of QgsLineString segment to verify for self intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        geom_engine_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet().clone())
        for qgs_geom_potential in qgs_geoms_with_itself:
            de_9IM_pattern = geom_engine_subline.relate(qgs_geom_potential.constGet().clone())
            # de_9IM_pattern[0] == '0' means that their interiors intersect (crosses)
            # de_9IM_pattern[1] == '0' means that one extremity is touching the interior of the other (touches)
            if de_9IM_pattern[0] == '0' or de_9IM_pattern[1] == '0':
                # The new sub line intersect or touch with itself. The result would create a non OGC simple line
                constraints_valid = False
                break

        return constraints_valid

    @staticmethod
    def validate_intersection(qgs_geom_with_others, qgs_geom_new_subline):
        """Validate the intersection constraint

        This constraint assure that the new sub line is not intersecting with any other lines (not itself)

        :param: qgs_geoms_with_others: List of QgsLineString segment to verify for intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        for qgs_geom_potential in qgs_geom_with_others:
            if not qgs_geom_potential.disjoint(qgs_geom_new_subline):
                # The bend area intersects with a point
                constraints_valid = False
                break

        return constraints_valid

    @staticmethod
    def validate_sidedness(qgs_geom_with_others, qgs_geom_bend):
        """Validate the sidedness constraint

        This constraint assure that the new sub line will not change the relative position of an object compared to
        the polygon formed by the bend to reduce. ex.: an interior ring of a polygon going outside of the exterior ring.

        :param: qgs_geoms_with_others: List of QgsLineString segment to verify for intersection
        :param: qgs_geom_bend: QgsPolygon formed by the bend to reduce
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        for qgs_geom_potential in qgs_geom_with_others:
            if qgs_geom_bend.contains(qgs_geom_potential):
                # A feature is totally located inside
                constraints_valid = False
                break

        return constraints_valid

    @staticmethod
    def find_alternate_bends(ind, rb_geom):
        """This method fins alternate possible bend

        This method is called when a bend reduction is causing self intersection.  In most cases this self intersection
        is caused by a bend with wave or spiral shape that cause the bend reduction to self intersect the line.
        This method finds a list of possible alternate bends.  Usually one of these alternate bend will allow
        to resolve the wave or spiral shape of the bend.

        :param: ind: Index number of the bend to process
        :param: rb_geom: Geometry used to validate constraints
        :return: List of potential bends
        :rtype: [Bend]
        """

        bend = rb_geom.bends[ind]
        alternate_bends = []
        j = bend.j
        qgs_points = rb_geom.qgs_geom.constGet().points()
        while j - 1 >= 2:
            i = bend.i
            while j - i >= 2:
                alternate_bend = Bend(i, j, qgs_points[i:j+1])
                alternate_bends.append((alternate_bend.area, alternate_bend))
                i += 1
            j -= 1
        # Sort the bends from the biggest area to the smallest (bigger area are better candidate)
        alternate_bends.sort(key=lambda item: item[0], reverse=True)
        # Just keep the Bend from the list
        alternate_bends = [bend for (area, bend) in alternate_bends]

        return alternate_bends

    __slots__ = ('qgs_in_features', 'diameter_tol', 'smooth_line', 'flag_del_outer', 'flag_del_inner',
                 'validate_structure', 'feedback', 'rb_collection', 'eps', 'rb_results', 'rb_features', 'rb_geoms',
                 'bends_reduced')

    def __init__(self, qgs_in_features, diameter_tol, smooth_line, flag_del_outer, flag_del_inner, validate_structure,
                 feedback):
        """Constructor for the bend reduction.

       :param: qgs_in_features: List of features to process.
       :param: diameter_tol: Float tolerance of the diameter of the bend to reduce.
       :param: smooth_line: Flag to smooth line after bend reduction if possible
       :param: flag_del_outer: Flag to delete polygon if area below the diameter tolerance.
       :param: flag_del_inner: Flag to delete polygon holes if area below the diameter tolerance.
       :param: validate_structure: flag to validate internal data structure after processing (for debugging)
       :param: feedback: QgsFeedback handle for interaction with QGIS.
       """

        self.qgs_in_features = qgs_in_features
        self.diameter_tol = diameter_tol
        self.smooth_line = smooth_line
        self.flag_del_outer = flag_del_outer
        self.flag_del_inner = flag_del_inner
        self.validate_structure = validate_structure
        self.feedback = feedback
        self.bends_reduced = []  # List containing the reduced bend
        self.eps = None
        self.rb_results = None
        self.rb_features = None
        self.rb_geoms = None
        self.rb_collection = None

    def reduce_bends(self):
        """Main method to manage bend reduction.

        :return: Statistics and result object.
        :rtype: RbResult
        """

        #  Code used for the profiler (uncomment if needed)
#        import cProfile, pstats, io
#        from pstats import SortKey
#        pr = cProfile.Profile()
#        pr.enable()

        # Calculates the epsilon and initialize some stats and results value
        self.eps = Epsilon(self.qgs_in_features)
        self.eps.set_class_variables()
        self.rb_results = RbResults()

        # Create the list of RbPolygon, RbLineString and RbPoint to process
        self.rb_features = self.create_rb_feature()
        self.rb_results.in_nbr_features = len(self.qgs_in_features)

        # Pre process the LineString: remove to close point and co-linear points
        self.rb_geoms = self.pre_reduction_process()

        # Create the RbCollection a spatial index to accelerate search
        self.rb_collection = RbCollection(self.rb_results)
        self.rb_collection.add_features(self.rb_geoms)

        # Execute the bend reduction for each LineString
        self._manage_reduce_bend()

        # Manage the line smoothing if needed
        if self.smooth_line:
            self.manage_smooth_line()

        # Recreate the QgsFeature
        qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in self.rb_features]

        # Set return values
        self.rb_results.out_nbr_features = len(qgs_features_out)
        self.rb_results.qgs_features_out = qgs_features_out

        # Validate inner spatial structure. For debug purpose only
        if self.rb_results.is_structure_valid:
            self.rb_collection.validate_integrity(self.rb_geoms)

        #  Code used for the profiler (uncomment if needed)
 #       pr.disable()
 #       s = io.StringIO()
 #       sortby = SortKey.CUMULATIVE
 #       ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
 #       ps.print_stats()
 #       print(s.getvalue())

        return self.rb_results

    def create_rb_feature(self):
        """Create the different RbFeatures from the QgsFeatures.

        :return: List of rb_features
        :rtype: [RbFeatures]
        """

        rb_features = []

        for qgs_feature in self.qgs_in_features:
            qgs_geom = qgs_feature.geometry()  # extract the Geometry

            if RbFeature.is_polygon(qgs_geom.wkbType()):
                rb_features.append(RbPolygon(qgs_feature))
            elif RbFeature.is_line_string(qgs_geom.wkbType()):
                rb_features.append(RbLineString(qgs_feature))
            elif RbFeature.is_point(qgs_geom.wkbType()):
                rb_features.append(RbPoint(qgs_feature))
            else:
                raise QgsProcessingException("Internal geometry error")

        return rb_features

    def pre_reduction_process(self):
        """This method execute the pre reduction process

        Pre reduction process includes remove small polygon or polygon hole simplify line when vertice
        are really too close and transform the RbFeature into LineString

        :return: List of rb_geom
        :rtype: [RbGeom]
        """

        # Delete the outer or inner ring below the diameter tolerance
        if self.flag_del_outer or self.flag_del_inner:
            self.del_outer_inner_ring()

        # Create the list of RbGeom ==> List of geometry to reduce the bend
        rb_geoms = []
        for rb_feature in self.rb_features:
            rb_geoms += rb_feature.get_rb_geom()

        # Remove duplicate nodes
        for rb_geom in rb_geoms:
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
                if not rb_geom.is_simplest:
                    rb_geom.qgs_geom.removeDuplicateNodes(epsilon=Epsilon.ZERO_RELATIVE)

        return rb_geoms

    def del_outer_inner_ring(self):
        """This method deletes the polygons and polygon holes below the diameter tolerance

        """

        # Loop over each rb_features
        for i in reversed(range(len(self.rb_features))):  # List visited in reverse order for easier deletion of entry
            if isinstance(self.rb_features[i], RbPolygon):  # Only process Polygon
                min_adj_area = ReduceBend.calculate_min_adj_area(self.diameter_tol)
                for j in reversed(range(len(self.rb_features[i].rb_geom))):  # List visited in reverse order (delete)
                    area, perimeter = ReduceBend._extract_polygon_attributes(self.rb_features[i].rb_geom[j].qgs_geom)
                    adj_area = self.calculate_adj_area(area, perimeter)
                    if j == 0:
                        # Process the exterior ring (outer ring always at position 0)
                        if self.flag_del_outer and adj_area < min_adj_area:
                            del self.rb_features[i]  # Delete the rb_feature
                            self.rb_results.nbr_pol_del += 1
                            break
                    else:
                        # Process an interior ring
                        if self.flag_del_inner and adj_area < min_adj_area:
                            del self.rb_features[i].rb_geom[j]  # Delete the ring
                            self.rb_results.nbr_hole_del += 1

        return

    def _manage_reduce_bend(self):
        """Loop over the geometry until there is no more bend to reduce

        An iterative process for bend reduction is needed in order to maximise the bend reduction.  The process
        will always stabilize and exit when there are no more bends to reduce.

        """

        min_nbr_pass = 2
        nbr_geoms = 100.0 / len(self.rb_geoms) if len(self.rb_geoms) >= 1 else 0
        while True:
            self.feedback.setProgress(max(1, int(self.count_rb_geoms_done() * nbr_geoms)))
            nbr_bend_reduced = 0
            nbr_bend_detected = 0
            for rb_geom in self.rb_geoms:
                if self.feedback.isCanceled():
                    break
                if not rb_geom.is_simplest:  # Only process geometry that are not at simplest form
                    self.delete_co_linear(rb_geom)
                    nbr_bend_detected = ReduceBend.detect_bends(rb_geom)
                    if rb_geom.need_pivot:
                        #  Pivoting a closed line moves the first/last vertice on a bend that do not need simplification
                        ReduceBend.pivot_closed_line(rb_geom, self.diameter_tol)
                        nbr_bend_detected = ReduceBend.detect_bends(rb_geom)  # Bend detection needed after pivot
                    ReduceBend.flag_bend_to_reduce(rb_geom, self.diameter_tol)
                    nbr_bend_reduced += self.process_bends(rb_geom)

            self.rb_results.nbr_bend_reduced.append(nbr_bend_reduced)
            self.rb_results.nbr_bend_detected.append(nbr_bend_detected)

            # While loop breaking condition
            if self.rb_results.nbr_pass > min_nbr_pass and nbr_bend_reduced == 0:
                break
            self.rb_results.nbr_pass += 1

        return

    def count_rb_geoms_done(self):
        """Count the number of geometry  that are at there simplest form

        """

        nbr_done = 0
        for rb_geom in self.rb_geoms:
            if rb_geom.is_simplest:
                nbr_done += 1

        return nbr_done

    def delete_co_linear(self, rb_geom):
        """Delete co-linear vertice on a LineString

        This method delete co-linear and near co-linear vertice because they are unnecessary but moreover in certain
        condition when they are forming near 0 (empty) area these vertices are creating spatial calculus errors

        :param: RbGeom rb_geom: Geometry to delete co-linear vertices
        """

        # Build the list of angles for each vertice
        vertex_ids_to_del = []
        angles = ReduceBend.get_angles(rb_geom.qgs_geom.constGet())
        if rb_geom.qgs_geom.constGet().isClosed() and len(angles) >= 1:
            del angles[0]  # Do not process the start/end vertice (even if co-linear)
        for i, angle in enumerate(angles):
            if abs(angle - math.pi) <= Epsilon.ZERO_ANGLE or abs(angle) <= Epsilon.ZERO_ANGLE:
                # Co-linear point or flat angle delete the current point
                vertex_ids_to_del.append(i+1)

        # Delete co-linear vertex
        for vertex_id_to_del in reversed(vertex_ids_to_del):
            self.rb_collection.delete_vertex(rb_geom, vertex_id_to_del, vertex_id_to_del)

        # Special case to process closed line string to find ans delete co-linear points at the first/last vertice
        if rb_geom.qgs_geom.constGet().isClosed():
            num_points = rb_geom.qgs_geom.constGet().numPoints()
            if num_points >= 5:  # Minimum of 5 vertices are needed to have co-linear vertices in closed line
                qgs_ls = QgsLineString([rb_geom.qgs_geom.vertexAt(num_points-2), \
                                        rb_geom.qgs_geom.vertexAt(0), \
                                        rb_geom.qgs_geom.vertexAt(1)])
                angles = ReduceBend.get_angles(qgs_ls)
                angle = angles[0]
                if abs(angle - math.pi) <= Epsilon.ZERO_ANGLE or abs(angle) <= Epsilon.ZERO_ANGLE:
                    self.rb_collection.delete_vertex(rb_geom, 0, 0)

        if rb_geom.qgs_geom.length() <= Epsilon.ZERO_RELATIVE:
            # Something wrong.  do not try to simplify the LineString
            rb_geom.is_simplest = True

        return

    def validate_alternate_bend(self, alternate_bends, ind, rb_geom):
        """Validate alternate bends

        For each alternate bend this method validate if the alternate bend is valid until an alternate bend is found
        or the list of alternate bends is all validated.

        :param: alternate_bends: List of alterate Bend
        :param: ind: ndex of the position of the bend in the list of bend
        :param: rb_geom: RbGeom geometry containing the bend to reduce
        :return: Flag indicating if one alternate bend is valid
        :rtype: boolean
        """

        constraints_valid = False
        for alternate_bend in alternate_bends:
            b_box = alternate_bend.qgs_geom_bend.boundingBox()
            qgs_geoms_with_itself, dummy = self.rb_collection.get_segment_intersect(rb_geom.id, b_box,
                                                                                    alternate_bend.qgs_geom_old_subline)

            new_bend_ok = True

            gqs_ls_new_subline = alternate_bend.qgs_geom_new_subline.constGet().clone()
            geom_engine_subline = QgsGeometry.createGeometryEngine(gqs_ls_new_subline)
            for qgs_geom_potential in qgs_geoms_with_itself:
                de_9IM_pattern = geom_engine_subline.relate(qgs_geom_potential.constGet().clone())
                # de_9IM_pattern[0] == '0' means that their interiors intersect (crosses)
                # de_9IM_pattern[1] == '0' means that one extremity is touching the interior of the other (touches)
                if de_9IM_pattern[0] == '0' or de_9IM_pattern[1] == '0':
                    # The new sub line intersect or touch with itself. The result would create a non OGC simple line
                    new_bend_ok = False
                    break

            if new_bend_ok:
                rb_geom.bends[ind] = alternate_bend  # Reset the bend with the alternate bend
                constraints_valid = True
                break  # No need to process the next alternate bend

        return constraints_valid

    def validate_constraints(self, ind, rb_geom):
        """Validate the spatial relationship in order maintain topological structure

        Three distinct spatial relation are tested in order to assure that each bend reduce will continue to maintain
        the topological structure in a feature between the features:
         - Simplicity: Adequate validation is done to make sure that the bend reduction will not cause the feature
                       to cross  itself.
         - Intersection : Adequate validation is done to make sure that a line from other features will not intersect
                          the bend being reduced
         - Sidedness: Adequate validation is done to make sure that a line is not completely contained in the bend.
                      This situation can happen when a ring in a polygon complete;y lie in a bend ans after bend
                      reduction, the the ring falls outside the polygon which make it invalid.

        Note if the topological structure is wrong before the bend correction no correction will be done on these
        errors.

        :param: ind: Index number of the bend to process
        :param: rb_geom: Geometry used to validate constraints
        :param: detect_alternate_bend: Indicates if alternate bend can be find when self intersection is detected
        :return: Flag indicating if the spatial constraints are valid for this bend reduction
        :rtype: Bool
        """

        constraints_valid = True
        bend = rb_geom.bends[ind]
        b_box = bend.qgs_geom_bend.boundingBox()
        qgs_geoms_with_itself, qgs_geoms_with_others = \
            self.rb_collection.get_segment_intersect(rb_geom.id, b_box, bend.qgs_geom_old_subline)

        # First: check if the bend reduce line string is an OGC simple line
        # We test with a tiny smaller line to ease the testing and false positive error
        if bend.qgs_geom_new_subline.length() >= Epsilon.ZERO_RELATIVE:
            constraints_valid = ReduceBend.validate_simplicity(qgs_geoms_with_itself, bend.qgs_geom_new_subline)
            if not constraints_valid:
                # The bend reduction caused self intersection; try to find an alternate bend
                alternate_bends = ReduceBend.find_alternate_bends(ind, rb_geom)
                constraints_valid = self.validate_alternate_bend(alternate_bends, ind, rb_geom)
        else:
            # Error in the input file
            qgs_line_string = bend.qgs_geom_new_subline.constGet()
            x = qgs_line_string.startPoint().x()
            y = qgs_line_string.startPoint().y()
            text = "Possibly non OGC simple feature at {},{} use Fix Geometries".format(x, y)
            self.feedback.pushInfo(text)

        # Second: check that the new line does not intersect any other line or points
        if constraints_valid:
            constraints_valid = ReduceBend.validate_intersection(qgs_geoms_with_others, bend.qgs_geom_new_subline)

        # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
        # sidedness or relative position error
        if constraints_valid:
            constraints_valid = ReduceBend.validate_sidedness(qgs_geoms_with_others, bend.qgs_geom_bend)

        return constraints_valid

    def validate_constraints_smooth(self, reduced_bend):
        """Validate the spatial relationship in order maintain topological structure for a smoothed line

        Three distinct spatial relation are tested in order to assure that each bend reduce will continue to maintain
        the topological structure between the feature:
         - Simplicity: Adequate validation is done to make sure that the bend reduction will not cause the feature
                       to cross  itself.
         - Intersection : Adequate validation is done to make sure that a line from other features will not intersect
                          the bend being reduced
         - Sidedness: Adequate validation is done to make sure that a line is not completely contained in the bend.
                      This situation can happen when a ring in a polygon complete;y lie in a bend ans after bend
                      reduction, the the ring falls outside the polygon which make it invalid.

        Note if the topological structure is wrong before the bend correction no correction will be done on these
        errors.

        :param: reduce_bend: Contains a ReducedBend to smooth
        :return: Flag indicating if the spatial constraints are valid for this bend reduction
        :rtype: Bool
        """

        rb_geom = reduced_bend.rb_geom
        b_box = reduced_bend.qgs_geom_smooth_polygon.boundingBox()
        qgs_geoms_with_itself, qgs_geoms_with_others = \
            self.rb_collection.get_segment_intersect(rb_geom.id, b_box, reduced_bend.qgs_geom_old_subline)

        # First: check if the bend reduce line string is an OGC simple line
        # We test with a tiny smaller line to ease the testing and false positive error
        constraints_valid = ReduceBend.validate_simplicity(qgs_geoms_with_itself,
                                                           reduced_bend.qgs_geom_smooth_line)

        # Second: check that the new line does not intersect any other line or points
        if constraints_valid:
            constraints_valid = ReduceBend.validate_intersection(qgs_geoms_with_others,
                                                                 reduced_bend.qgs_geom_smooth_line)

        # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
        # sidedness or relative position error
        if constraints_valid:
            constraints_valid = ReduceBend.validate_sidedness(qgs_geoms_with_others,
                                                              reduced_bend.qgs_geom_smooth_polygon)

        return constraints_valid

    def process_bends(self, rb_geom):
        """This method manages the reduction of the bends of one lineString

        :param: RbGeom rb_geom: Geometry
        :param: int ind: index number of the bend to process
        :return: Flag indicating if the spatial constraints are valid for this bend reduction
        :rtype: Bool
        """

        nbr_bend_reduced = 0
        for ind in reversed(range(len(rb_geom.bends))):  # Work reversely to be able to delete easily in the list
            if rb_geom.bends[ind].to_reduce:
                # Check spatial constraints
                spatial_constraints = self.validate_constraints(ind, rb_geom)
                if spatial_constraints:
                    bend = rb_geom.bends[ind]
                    if self.smooth_line:
                        qgs_pnt_i = rb_geom.qgs_geom.vertexAt(bend.i)
                        qgs_pnt_j = rb_geom.qgs_geom.vertexAt(bend.j)
                        self.bends_reduced.append(BendReduced(rb_geom, qgs_pnt_i, qgs_pnt_j, bend.qgs_geom_bend))
                    self.rb_collection.delete_vertex(rb_geom, bend.i+1, bend.j-1)
                    nbr_bend_reduced += 1

        return nbr_bend_reduced

    def manage_smooth_line(self):
        """Manage the line smoothing of all the bend reduced to a straight line

        The smoothing of the line is done after all bend reduction is done.  The spatial constraints are validated
        ob the smoothed line to make sure the smooth line respect spatial constraints
        """

        for bend_reduced in self.bends_reduced:
            bend_reduced.set_values(self.diameter_tol)
            if bend_reduced.is_line_smoothable:
                if self.validate_constraints_smooth(bend_reduced):
                    self.rb_collection.add_vertex(bend_reduced.rb_geom, bend_reduced.i, bend_reduced.j,
                                                  bend_reduced.qgs_geom_smooth_line)
                    self.rb_results.nbr_line_smooth += 1
                else:
                    # Smooth line validation of spatial constraint did not pass stay with the straight line
                    pass
            else:
                # Smoothing of the line impossible
                pass
