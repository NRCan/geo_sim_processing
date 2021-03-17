from qgis.PyQt.QtCore import QCoreApplication
from PyQt5.QtGui import QTransform
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterDistance,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean)

from abc import ABC, abstractmethod
import sys
import math
from qgis.core import QgsFeatureSink, QgsFeatureRequest, QgsFeature, QgsPoint, QgsPointXY, QgsLineString, QgsPolygon, \
    QgsWkbTypes, QgsSpatialIndex, QgsGeometry, QgsGeometryUtils, QgsRectangle, \
    QgsProcessingException, QgsMultiLineString, QgsMultiPolygon
import os
import inspect
import processing
from qgis.PyQt.QtGui import QIcon, QPainter
from qgis.utils import iface
from qgis.gui import QgsRubberBand


class BuildingPatternAlgorithm(QgsProcessingAlgorithm):
    """Main class defining the Reduce Bend as a QGIS processing algorithm.
    """

    def tr(self, string):
        """Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        """Returns a new copy of the algorithm.
        """
        return BuildingPatternAlgorithm()

    def name(self):
        """Returns the unique algorithm name.
        """
        return 'buildingpattern'

    def displayName(self):
        """Returns the translated algorithm name.
        """
        return 'Building Pattern'

    def group(self):
        """Returns the name of the group this algorithm belongs to.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """Returns the unique ID of the group this algorithm belongs to.
        """
        return ''

    def shortHelpString(self):
        """Returns a localised short help string for the algorithm.
        """
        help_str = """
    Building pattern is a geospatial simplification and generalization tool for lines and polygons. The \
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


    """

        return self.tr(help_str)

    def icon(self):
        """Define the logo of the algorithm.
        """

        cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]
        icon = QIcon(os.path.join(os.path.join(cmd_folder, 'logo.png')))
        return icon

    def initAlgorithm(self, config=None):
        """Define the inputs and outputs of the algorithm.
        """

        # 'INPUT' is the recommended name for the main input parameter.
        self.addParameter(QgsProcessingParameterFeatureSource(
            'INPUT',
            self.tr('Input layer'),
            types=[QgsProcessing.TypeVectorAnyGeometry]))


        # 'RECTANGULARITY' to be used bor bend reduction
        self.addParameter(QgsProcessingParameterDistance(
            'RECTANGULARITY',
            self.tr('Rectangularity tolerance [0..1]'),
            defaultValue=0.0,
            parentParameterName='INPUT'))  # Make distance units match the INPUT layer units

        # 'COMPACTNESS' to be used bor bend reduction
        self.addParameter(QgsProcessingParameterDistance(
            'COMPACTNESS',
            self.tr('Compactness tolerance [0..1]'),
            defaultValue=0.0,
            parentParameterName='INPUT'))  # Make distance units match the INPUT layer units

        # 'VERBOSE' mode for more output information
        self.addParameter(QgsProcessingParameterBoolean(
            'VERBOSE',
            self.tr('Verbose'),
            defaultValue=False))

        # 'OUTPUT' for the results
        self.addParameter(QgsProcessingParameterFeatureSink(
            'OUTPUT',
            self.tr('Building pattern')))

    def processAlgorithm(self, parameters, context, feedback):
        """Main method that extract parameters and call ReduceBend algorithm.
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        # Extract parameter
        source_in = self.parameterAsSource(parameters, "INPUT", context)
        rectangularity_tol = self.parameterAsDouble(parameters, "RECTANGULARITY", context)
        compactness_tol = self.parameterAsDouble(parameters, "COMPACTNESS", context)
        verbose = self.parameterAsBool(parameters, "VERBOSE", context)

        if source_in is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        # Transform the in source into a vector layer
        vector_layer_in = source_in.materialize(QgsFeatureRequest(), feedback)

        # Normalize and extract QGS input features
        qgs_features_in, geom_type = BuildingPattern.normalize_in_vector_layer(vector_layer_in, feedback)

        # Validate input geometry type
        if geom_type not in [QgsWkbTypes.Polygon]:
            raise QgsProcessingException("Can only process: (Multi)Polygon vector layers")

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                               vector_layer_in.fields(),
                                               geom_type,
                                               vector_layer_in.sourceCrs())

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))

        # Set progress bar to 1%
        feedback.setProgress(1)

        # Call BuildingPattern algorithm
        bp_return = BuildingPattern.match(qgs_features_in, rectangularity_tol, compactness_tol, feedback)

        for qgs_feature_out in bp_return.qgs_features_out:
            sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

        # Push some output statistics
        feedback.pushInfo("Number of features in: {0}".format(bp_return.in_nbr_features))
        feedback.pushInfo("Number of features out: {0}".format(bp_return.out_nbr_features))
        feedback.pushInfo("Number of rectangularity corrections: {0}".format(bp_return.nbr_rectangularity))
        feedback.pushInfo("Number of compactness corrections: {0}".format(bp_return.nbr_compactness))
        feedback.pushInfo("Number of pattern corrections: {0}".format(bp_return.nbr_pattern))

        return {"OUTPUT": dest_id}


# --------------------------------------------------------
# Start of the algorithm
# --------------------------------------------------------

# Define global constant


class BpResults:
    """Class defining the stats and result"""

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'qgs_features_out', 'nbr_rectangularity', 'nbr_compactness',\
                 'nbr_pattern')

    def __init__(self):
        """Constructor that initialize a RbResult object.

        :param: None
        :return: None
        :rtype: None
        """

        self.in_nbr_features = 0
        self.out_nbr_features = 0
        self.qgs_features_out = []
        self.nbr_rectangularity = 0
        self.nbr_compactness = 0
        self.nbr_pattern = 0


class Epsilon:
    """Class defining the value of the zero"""

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
        log_loss = int(math.log(dynamic_xy, 10) + 1)
        max_digit = 15  # Number of significative digits for real number
        security = 2  # Keep 2 order of magnitude of security
        abs_digit = max_digit - security
        rel_digit = max_digit - log_loss - security
        self._zero_relative = (1. / (10 ** rel_digit))
        self._zero_absolute = (1. / (10 ** abs_digit))
        self._zero_angle = math.radians(.0001)  # Angle used to decide a flat angle

    def set_class_variables(self):
        """Set the different epsilon values.

        :return: None
        :rtype: None
        """

        BuildingPattern.ZERO_RELATIVE = self._zero_relative
        BuildingPattern.ZERO_ABSOLUTE = self._zero_absolute
        BuildingPattern.ZERO_ANGLE = self._zero_angle

        return


class BuildingPattern:
    """Main class for bend reduction"""

    ZERO_ANGLE = None
    ZERO_RELATIVE = None
    ZERO_ABSOLUTE = None

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
    def match(qgs_in_features, rectangularity_tol, compactness_tol, feedback=None):
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

        pb = BuildingPattern(qgs_in_features, rectangularity_tol, compactness_tol, feedback)
        results = pb.building_pattern()

        return results



    __slots__ = ('qgs_in_features', 'rectangularity_tol', 'compactness_tol', 'feedback', 'eps', 'bp_results')

    def __init__(self, qgs_in_features, rectangularity_tol, compactness_tol, feedback):
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
        self.rectangularity_tol = rectangularity_tol
        self.compactness_tol = compactness_tol
        self.feedback = feedback
        self.eps = None
        self.bp_results = None

    def create_geom(self, points):

        qgs_pnts = []
        for point in points:
            qgs_pnts.append(point[0], point[1])
        qgs_geom_pol = QgsGeometry(QgsLineString(qgs_pnts))

        return qgs_geom_pol

    def building_pattern(self):

        self.bp_results = BpResults()

        patterns = []
        # Pattern in L form
        qgs_geom = QgsLineString([QgsPoint(0,0), QgsPoint(0,1), QgsPoint(1,1), QgsPoint(1,.8), QgsPoint(.8,.8), QgsPoint(.8,0), QgsPoint(0,0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .8), QgsPoint(.6, .8), QgsPoint(.6, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .8), QgsPoint(.5, .8), QgsPoint(.5, 0),QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .8), QgsPoint(.4, .8), QgsPoint(.4, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .6), QgsPoint(.6, .6), QgsPoint(.6, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .5), QgsPoint(.5, .5), QgsPoint(.5, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(1, 1), QgsPoint(1, .4), QgsPoint(.4, .4), QgsPoint(.4, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        # Pattern in T form
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(.8, 1), QgsPoint(.8, .8), QgsPoint(1, .8), QgsPoint(1, .6), QgsPoint(.8, .6), QgsPoint(.8, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(.8, 1), QgsPoint(.8, .8), QgsPoint(1, .8), QgsPoint(1, .4), QgsPoint(.8, .4), QgsPoint(.8, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)
        qgs_geom = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 1), QgsPoint(.6, 1), QgsPoint(.6, .8), QgsPoint(1, .8), QgsPoint(1, .4), QgsPoint(.6, .4), QgsPoint(.6, 0), QgsPoint(0, 0)])
        qgs_geom = QgsGeometry(QgsPolygon(qgs_geom.clone()))
        patterns += Pattern.build_pattern(qgs_geom)

#        y_delta = 0
#        for i in range(len(patterns)):
#            qgs_geom = patterns[i].qgs_geom_in_pattern
#            qgs_geom.translate(0, y_delta)
#            y_delta += 2
#            feat = QgsFeature()
#            feat.setGeometry(qgs_geom)#
#            self.bp_results.qgs_features_out.append(feat)
#        return self.bp_results

        for qgs_feature in self.qgs_in_features:
            self.bp_results.in_nbr_features += 1
            tf = TargetFeature(qgs_feature)
            replaced = False
            if tf.replace_bbox_rectangularity(self.rectangularity_tol):
                replaced = True
                self.bp_results.nbr_rectangularity += 1
            elif tf.replace_bbox_compactness(self.compactness_tol):
                replaced = True
                self.bp_results.nbr_compactness += 1
            elif tf.replace_pattern(patterns, .7):
                replaced = True
                self.bp_results.nbr_pattern += 1

            self.bp_results.qgs_features_out.append(tf.qgs_feature)

        return self.bp_results

class TargetFeature:

    __slots__ = ('qgs_feature', 'qgs_geom', 'qgs_pnt_centroid', 'area', 'perimeter', 'area_bbox', 'angle',\
                 'width', 'height', 'rectangularity_ratio', 'compactness_ratio')

    def __init__(self, qgs_feature):

        self.qgs_feature = qgs_feature
        self.qgs_geom = qgs_feature.geometry()
        qgs_geom_centroid = self.qgs_geom.centroid()
        self.qgs_pnt_centroid = qgs_geom_centroid.constGet().clone()
        self.area = self.qgs_geom.area()
        self.perimeter = self.qgs_geom.constGet().perimeter()
        tuple_results = self.qgs_geom.orientedMinimumBoundingBox()
        self.area_bbox = tuple_results[1]
        self.angle = tuple_results[2]
        self.width = tuple_results[3]
        self.height = tuple_results[4]
        self.rectangularity_ratio = self.area/(self.width*self.height)
        self.compactness_ratio = 4.* self.area * math.pi / (self.perimeter**2)

    def replace_bbox_rectangularity(self, rectangularity_tol):

        feature_replaced = False
        if self.rectangularity_ratio > rectangularity_tol:
            self._replace_bbox()
            feature_replaced = True

        return feature_replaced

    def replace_bbox_compactness(self, compactness_tol):

        feature_replaced = False
        if self.compactness_ratio > compactness_tol:
            self._replace_bbox()
            feature_replaced = True

        return feature_replaced

    def replace_pattern(self, patterns, tolerance_ratio):

        max_tolerance_ratio = 0.0
        for pattern in patterns:
            qgs_geom_pattern = QgsGeometry(pattern.qgs_geom_in_pattern)
            qgs_geom_pattern.transform(QTransform().scale(self.width, self.height))
            qgs_geom_pattern.rotate(self.angle, QgsPointXY(0,0))
            qgs_centroid_pattern = qgs_geom_pattern.centroid()
            qgs_pnt_centroid_rect = qgs_centroid_pattern.constGet().clone()
            delta_x = self.qgs_pnt_centroid.x() - qgs_pnt_centroid_rect.x()
            delta_y = self.qgs_pnt_centroid.y() - qgs_pnt_centroid_rect.y()
            qgs_geom_pattern.translate(delta_x, delta_y)
            area_pattern = qgs_geom_pattern.area()

            qgs_a_diff_b = qgs_geom_pattern.difference(self.qgs_geom)
            qgs_b_diff_a = self.qgs_geom.difference(qgs_geom_pattern)
            qgs_symmetric_diff = QgsGeometry.unaryUnion([qgs_a_diff_b, qgs_b_diff_a])

            area_ratio = 1 - (qgs_symmetric_diff.area() / ((area_pattern+self.area)/2.))
            if area_ratio > tolerance_ratio:
                if area_ratio > max_tolerance_ratio:
                    max_tolerance_ratio = area_ratio
                    self.qgs_feature.setGeometry(qgs_geom_pattern)
        if max_tolerance_ratio > tolerance_ratio:
            replaced = True
        else:
            replaced = False

        return replaced


    def _replace_bbox(self):
        shrink_bbox = self.rectangularity_ratio ** .5
        xmin = shrink_bbox
        ymin = shrink_bbox
        xmax = self.width - shrink_bbox
        ymax = self.height - shrink_bbox
        qgs_rect_feature = QgsLineString([xmin, xmin, xmax, xmax], [ymin, ymax, ymax,ymin])
        qgs_rect_target = QgsGeometry(QgsPolygon(qgs_rect_feature.clone()))
        qgs_rect_target.rotate(self.angle, QgsPointXY(0., 0.))
        qgs_centroid_rect = qgs_rect_target.centroid()
        qgs_pnt_centroid_rect = qgs_centroid_rect.constGet().clone()
        delta_x = self.qgs_pnt_centroid.x() - qgs_pnt_centroid_rect.x()
        delta_y = self.qgs_pnt_centroid.y() - qgs_pnt_centroid_rect.y()
        qgs_rect_target.translate(delta_x, delta_y)
        self.qgs_feature.setGeometry(qgs_rect_target)


class Pattern:

    @staticmethod
    def scale_geom(qgs_geom_pattern, x_scale, y_scale):

        qgs_geom = QgsGeometry(qgs_geom_pattern)
        qgs_geom.transform(QTransform().scale(x_scale, y_scale))
        qgs_geom_bbox = qgs_geom.boundingBox()
        xmin = qgs_geom_bbox.xMinimum()
        ymin = qgs_geom_bbox.yMinimum()
        qgs_geom.translate(-xmin, -ymin)

        return qgs_geom

    @staticmethod
    def build_pattern(qgs_geom_in_pattern, qgs_geom_out_pattern=None):

        if qgs_geom_out_pattern is None:
            qgs_geom_out_pattern = QgsGeometry(qgs_geom_in_pattern)

        patterns = []

        for x_scale, y_scale in ((1.,1.), (-1.,1),(1,-1.),(-1., -1.)):
            qgs_geom_in = Pattern.scale_geom(qgs_geom_in_pattern, x_scale, y_scale)
            qgs_geom_out = Pattern.scale_geom(qgs_geom_out_pattern, x_scale, y_scale)
            patterns.append(Pattern(qgs_geom_in, qgs_geom_out))

        return patterns

    def __init__(self, qgs_geom_in_pattern, qgs_geom_out_pattern):

        self.qgs_geom_in_pattern = qgs_geom_in_pattern
        self.qgs_geom_out_pattern = qgs_geom_out_pattern






