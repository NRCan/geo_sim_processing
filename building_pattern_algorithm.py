from qgis.PyQt.QtCore import QCoreApplication
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
from qgis.PyQt.QtGui import QIcon


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

        # 'TOLERANCE' to be used bor bend reduction
        self.addParameter(QgsProcessingParameterDistance(
            'TOLERANCE',
            self.tr('Diameter tolerance'),
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
            self.tr('Reduced bend')))

    def processAlgorithm(self, parameters, context, feedback):
        """Main method that extract parameters and call ReduceBend algorithm.
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        # Extract parameter
        source_in = self.parameterAsSource(parameters, "INPUT", context)
        rectangularity_tol = self.parameterAsDouble(parameters, "TOLERANCE", context)
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
        rb_return = BuildingPattern.match(qgs_features_in, rectangularity_tol, feedback)

        for qgs_feature_out in rb_return.qgs_features_out:
            sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

        # Push some output statistics
#        feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
#        feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
#        feedback.pushInfo("Number of iteration needed: {0}".format(rb_return.nbr_pass))
#        feedback.pushInfo("Number of bends detected: {0}".format(rb_return.nbr_bend_detected))
#        feedback.pushInfo("Number of bends reduced: {0}".format(rb_return.nbr_bend_reduced))
#        feedback.pushInfo("Number of deleted polygons: {0}".format(rb_return.nbr_pol_del))
#        feedback.pushInfo("Number of deleted polygon holes: {0}".format(rb_return.nbr_hole_del))
#        feedback.pushInfo("Number of line smoothed: {0}".format(rb_return.nbr_line_smooth))

        return {"OUTPUT": dest_id}


# --------------------------------------------------------
# Start of the algorithm
# --------------------------------------------------------

# Define global constant
ANTI_CLOCK_WISE = -1
CLOCK_WISE = 0


class BpResults:
    """Class defining the stats and result"""

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'nbr_bend_reduced', 'nbr_bend_detected',
                 'qgs_features_out', 'nbr_hole_del', 'nbr_pol_del', 'nbr_pass', 'is_structure_valid',
                 'lines_log_info', 'nbr_line_smooth')

    def __init__(self):
        """Constructor that initialize a RbResult object.

        :param: None
        :return: None
        :rtype: None
        """

        self.in_nbr_features = None
        self.out_nbr_features = None
        self.qgs_features_out = []
        self.nbr_line_smooth = 0
        self.lines_log_info = []


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
    def match(qgs_in_features, rectangularity_tol, feedback=None):
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

        pb = BuildingPattern(qgs_in_features, rectangularity_tol, feedback)
        results = pb.building_pattern()

        return results



    __slots__ = ('qgs_in_features', 'rectangularity_tol', 'feedback', 'eps', 'bp_results')

    def __init__(self, qgs_in_features, rectangularity_tol, feedback):
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
        self.feedback = feedback
        self.eps = None
        self.bp_results = None

    def building_pattern(self):
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
#        self.eps = Epsilon(self.qgs_in_features)
#        self.eps.set_class_variables()
#        self.rb_results = RbResults()

        self.bp_results = BpResults()

        for qgs_feature in self.qgs_in_features:
            qgs_geom_target = qgs_feature.geometry()
            area_geom_target = qgs_geom_target.area()
            qgs_geom_centroid = qgs_geom_target.centroid()
            qgs_pnt_centroid_target = qgs_geom_centroid.constGet().clone()
            tuple_results = qgs_geom_target.orientedMinimumBoundingBox()
            qgs_geom_orient_bbox = tuple_results[0]
            area_orient_bbox = tuple_results[1]
            angle_orient_bbox = tuple_results[2]
            width_orient_bbox = tuple_results[3]
            height_orient_bbox = tuple_results[4]
            ratio_rectangularity = area_geom_target / area_orient_bbox
            if ratio_rectangularity > self.rectangularity_tol:
                growth = ratio_rectangularity**.5
                xmin = growth
                ymin = growth
                xmax = width_orient_bbox - growth
                ymax = height_orient_bbox - growth
                qgs_pnt0 = QgsPoint(xmin,ymin)
                qgs_pnt1 = QgsPoint(xmin, ymax)
                qgs_pnt2 = QgsPoint(xmax, ymax)
                qgs_pnt3 = QgsPoint(xmax, ymin)
                qgs_pnts = [qgs_pnt0, qgs_pnt1, qgs_pnt2, qgs_pnt3]
                qgs_rect_target = QgsGeometry(QgsPolygon(QgsLineString(qgs_pnts)))
                qgs_rect_target.rotate(angle_orient_bbox, QgsPointXY(0.,0.))
                qgs_centroid_rect = qgs_rect_target.centroid()
                qgs_pnt_centroid_rect = qgs_centroid_rect.constGet().clone()
                delta_x = qgs_pnt_centroid_target.x() - qgs_pnt_centroid_rect.x()
                delta_y = qgs_pnt_centroid_target.y() - qgs_pnt_centroid_rect.y()
                qgs_rect_target.translate(delta_x, delta_y)
                qgs_feature.setGeometry(qgs_rect_target)

            self.bp_results.qgs_features_out.append(qgs_feature)

#            area_orient_bbox = qgs_geom_orient_bbox.area()
#            ratio_rectangularity = area_geom_target / area_orient_bbox
#            if ratio_rectangularity > self.rectangularity_tol:
#                qgs_rect_ext_ring = qgs_geom_orient_bbox.constGet().exteriorRing()
#                qgs_points = qgs_rect_ext_ring.points()
#                qgs_l0 = QgsLineString(qgs_points[0], qgs_points[1])
#                qgs_l1 = QgsLineString(qgs_points[1], qgs_points[2])
#                qgs_l2 = QgsLineString(qgs_points[2], qgs_points[3])
#                qgs_l3 = QgsLineString(qgs_points[3], qgs_points[4])
#                ratio_rect_square = ratio_rectangularity**0.5
#                qgs_pnt0 = qgs_l0.interpolatePoint(ratio_rect_square*qgs_l0.length())
#                qgs_pnt1 = qgs_l1.interpolatePoint(ratio_rect_square*qgs_l1.length())
#                qgs_pnt2 = qgs_l2.interpolatePoint(ratio_rect_square*qgs_l2.length())
#                qgs_pnt3 = qgs_l3.interpolatePoint(ratio_rect_square*qgs_l3.length())
#                print (qgs_l0.length(), qgs_l1.length(), qgs_l2.length(), qgs_l3.length())
#                qgs_geom_rect = QgsGeometry(QgsPolygon(QgsLineString([qgs_pnt0, qgs_pnt1, qgs_pnt2, qgs_pnt3])))
#                qgs_geom_centroid_rect = qgs_geom_rect.centroid()
#                qgs_pnt_centroid_rect = qgs_geom_centroid_rect.constGet().clone()
#                delta_x = qgs_pnt_centroid_target.x() - qgs_pnt_centroid_rect.x()
#                delta_y = qgs_pnt_centroid_target.y() - qgs_pnt_centroid_rect.y()
#                qgs_geom_rect.translate(delta_x,delta_y)
#                qgs_feature.setGeometry(qgs_geom_rect)
#
#            self.bp_results.qgs_features_out.append(qgs_feature)

        return self.bp_results


