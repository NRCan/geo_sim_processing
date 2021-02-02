from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingException,
                       QgsProcessingParameterDistance,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterVectorDestination,
                       QgsWkbTypes,
                       QgsTopologyPreservingSimplifier,
                       QgsGeometry,
                       QgsFeatureSink,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean)

from abc import ABC, abstractmethod
import sys, math
from qgis.core import QgsFeatureSink, QgsFeatureRequest, QgsFeature, QgsLineString, QgsPolygon, QgsWkbTypes, \
                      QgsSpatialIndex,  QgsGeometry, QgsGeometryUtils, QgsRectangle, QgsProcessingException, \
                      QgsMultiLineString
import os
import inspect
from qgis.PyQt.QtGui import QIcon

class ReduceBendAlgorithm(QgsProcessingAlgorithm):
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
        return ReduceBendAlgorithm()

    def name(self):
        """
        Returns the unique algorithm name.
        """
        return 'reducebend'

    def displayName(self):
        """
        Returns the translated algorithm name.
        """
        return 'Reduce Bend'

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
    Reduce bend is a geospatial simplification and generalization tool for lines and polygons. \
    Reduce bend is an implementation and an improvement of the algorithm described in the paper \
    "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and \
    Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm". The \
    particularity of this algorithm is that for each line or polygon it analyzes its bends (curves) and \
    decides which one to simplify, trying to emulate what a cartographer would do manually \
    to simplify or generalize a line. Reduce bend will accept lines and polygons as input. \

    <b>Usage</b>

    <u>Input layer</u> : Any Line string or Polygon layer

    <u>Diameter tolerance</u>: Theoritical diameter of a bend to remove

    <u>Exclude hole</u>: If you want to exclude holes below the diameter of the bend

    <u>Exclude polygon</u>: If you want to exclude polygon below the diameter of the bend
    
    <u>Reduced bend</u> : Output layer of the algorithm

    <b>Rule of thumb for the diameter</b>
    Reduce bend can be used for line simplifying in the context of line generalization. The big \
    question will often be what diameter should we use? A good starting point is the cartographic rule of \
    thumb -- the .5mm on the map -- which says that the minimum distance between two lines should be \
    greater than 0.5mm on a paper map. So to simplify (generalize) a line for representation at a scale of \
    1:50 000 for example a diameter of 25m should be a good starting point.

    for more information: https:...

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
            QgsProcessingParameterDistance(
                'TOLERANCE',
                self.tr('Diameter tolerance'),
                defaultValue = 0.0,
                # Make distance units match the INPUT layer units:
                parentParameterName='INPUT'
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                'VERBOSE',
                self.tr('Verbose'),
                defaultValue=False
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                'EXCLUDE_POLYGON',
                self.tr('Exclude polygon'),
                defaultValue=True
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                'EXCLUDE_HOLE',
                self.tr('Exclude hole'),
                defaultValue=True
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                'OUTPUT',
                self.tr('Reduced bend')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
       
        """

        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

        source = self.parameterAsSource(parameters, "INPUT", context )
        exclude_hole = self.parameterAsBool(parameters, "EXCLUDE_HOLE", context)
        exclude_polygon = self.parameterAsBool(parameters, "EXCLUDE_POLYGON", context)
        diameter_tol = self.parameterAsDouble(parameters, "TOLERANCE", context)
        validate_structure = self.parameterAsBool(parameters, "VALIDATE_STRUCTURE", context)
        verbose = self.parameterAsBool(parameters, "VERBOSE", context)

        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        # Validate input source type
        if RbFeature.is_polygon(source.wkbType()):
            type = QgsWkbTypes.Polygon
        elif RbFeature.is_line_string(source.wkbType()):
            type = QgsWkbTypes.LineString
        else:
            #  Cannot process this feature type
            raise QgsProcessingException("Can only process: LineString or Polygon layers")

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                                   source.fields(),
                                                   type,
                                                   source.sourceCrs())

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))

        features = source.getFeatures()
        qgs_features_in = []
        for qgs_feature_in in features:
            # Load all the QgsFeature
            qgs_features_in.append(qgs_feature_in)
        try:
            # Call the bend reduction
            rb_return = ReduceBend.reduce(qgs_features_in, diameter_tol, feedback, exclude_polygon,
                                          exclude_hole, validate_structure)
        except Exception:
            import traceback
            traceback.print_exc()

        for qgs_feature_out in rb_return.qgs_features_out:
            sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

        # Push some output statistics
        feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
        feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
        feedback.pushInfo("Number of iteration needed: {0}".format(rb_return.nbr_pass))
        feedback.pushInfo("Number of bends detected: {0}".format(rb_return.nbr_bend_detected))
        feedback.pushInfo("Number of bends reduced: {0}".format(rb_return.nbr_bend_reduced))
        feedback.pushInfo("Number of deleted polygons: {0}".format(rb_return.nbr_pol_del))
        feedback.pushInfo("Number of deleted polygon holes: {0}".format(rb_return.nbr_hole_del))
        if validate_structure:
            if rb_return.is_structure_valid:
                status = "Valid"
            else:
                status = "Invalid"
            feedback.pushInfo("Debug - State of the internal data structure: {0}".format(status))
        if verbose:
            for line_log_info in rb_return.lines_log_info:
                feedback.pushInfo("Verbose - {0}".format(line_log_info))

        return {"OUTPUT": dest_id}







#Remaining modifications to do (to place as issues):
#  - put comment in the code
#  - verify if progrssive erosion bring a plus (min_pass)
#  - Add some comprehension list instead of for loop
#  - edit line with a smooth line instead of a straight line (new parameter)
#  - edit when create an acute angle could mesaure the average angle and use it at a floor (new parameter)
#  - in the bend flagging process prioritize the bend that goes outside the polygon first (for polygon)
#  - put some code to correct J-Bend
#  - Do not delete the line segment in the spatial index as it is very expensive mark them as delete (feature id -1)
#  - Put __slotss in RbFeature


# Define global constant
ANTI_CLOCK_WISE = -1
CLOCK_WISE = 0


class RbFeature(ABC):
    """Contain one QgsFeature

    Abstract class specialized into scpecific geometries
    """

    _id_counter = 0  # Counter of feature

    @staticmethod
    def is_point(feature_type):
        """Static method which determine if a QgsFeature is any kind of Point.

        :param int feature_type: Feature type to validate.
        :return: True if a point False otherwise
        :rtype: bool
        """

        if feature_type in [QgsWkbTypes.Point, QgsWkbTypes.Point25D, QgsWkbTypes.PointM, QgsWkbTypes.PointZ,
                            QgsWkbTypes.PointZM]:
            val = True
        else:
            val = False

        return val

    @staticmethod
    def is_line_string(feature_type):
        """Static method which determine if a QgsFeature is any kind of LineString.

        :param int feature_type: Feature type to validate.
        :return: True if a LineString False otherwise
        :rtype: bool
        """

        if feature_type in [QgsWkbTypes.LineString, QgsWkbTypes.LineString25D, QgsWkbTypes.LineStringZ,
                            QgsWkbTypes.LineStringM, QgsWkbTypes.LineStringZM]:
            val = True
        else:
            val = False

        return val

    @staticmethod
    def is_polygon(feature_type):
        """Static method which determine if a QgsFeature is any kind of Polygon.

        :param int feature_type: Feature type to validate.
        :return: True if a Polygon False otherwise
        :rtype: bool
        """
        if feature_type in [QgsWkbTypes.Polygon, QgsWkbTypes.Polygon25D, QgsWkbTypes.PolygonZ, QgsWkbTypes.PolygonM,
                            QgsWkbTypes.PolygonZM]:
            val = True
        else:
            val = False

        return val

    def __init__(self, qgs_feature):
        """Constructor of the RbFeature class.

        :param QgsFeature feature_type: QgsFeature to process.
        :return: None
        :rtype: None
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

        pass

    @abstractmethod
    def get_qgs_feature(self):
        """Define an abstract method.
         """

        pass


class RbPolygon(RbFeature):

    def __init__(self, qgs_feature):
        """Constructor that breaks the Polygon into a list of closed LineString (RbGeom).

        :param QgsFeature feature_type: QgsFeature polygon to process.
        :return: None
        :rtype: None
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

        :param None.
        :return: The RbGeom of the instance
        :rtype: List of RbGeom
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the new geometry.

        :param None.
        :return: The new Qgsfeature
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

        :param QgsFeature feature_type: QgsFeature LineString to process.
        :return: None
        :rtype: None
        """
        super().__init__(qgs_feature)
        if self.qgs_geom.wkbType() != QgsWkbTypes.LineString:
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.LineString)  # Force geometry to a QgsPoint
        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.LineString)]
        self.qgs_geom = None

    def get_rb_geom(self):
        """Return the RbGeom.

        :param None.
        :return: The RbGeom of the instance
        :rtype: List of RbGeom
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the new geometry.

        :param None
        :return: The new Qgsfeature
        :rtype: QgsFeature
        """

        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbPoint(RbFeature):
    """Class managin a RbPoint
    """

    def __init__(self, qgs_feature):
        """Constructor that breaks the Point into a list of Point (RbGeom).

        :param QgsFeature feature_type: QgsFeature Point to process.
        :return: None
        :rtype: None
        """

        super().__init__(qgs_feature)
        if self.qgs_geom.wkbType() != QgsWkbTypes.Point:
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.Point)  # Force geometry to QgsPoint
        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.Point)]
        self.rb_geom[0].is_simplest = True  # A point cannot be reduced
        self.qgs_geom = None

    def get_rb_geom(self):
        """Return the RbGeom.

        :param None.
        :return: The RbGeom of the instance.
        :rtype: List of RbGeom.
        """

        return self.rb_geom

    def get_qgs_feature(self):
        """Reconstruct the original QgsFeature with the original geometry.

        A Point cannot be reduced but is needed for the spatial constraints

        :param None
        :return: The new Qgsfeature
        :rtype: QgsFeature
        """

        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbCollection(object):
    """Class used for managing the feature spatially.

    QgsSpatialIndex class is used to store and retreive the features.
    """

    __slots__ = ('_spatial_index', '_dict_qgs_rb_geom', '_dict_qgs_segment', 'rb_results', '_id_qgs_segment')

    def __init__(self, rb_results):
        """Constructor that initialize the RbCollection.

        :param RbResults rb_results: Object containing the results and statistics of the execution
        :return: None
        :rtype: None
        """

        self._spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
        self._dict_qgs_segment = {}  # Contains a reference to the original geometry
        self._id_qgs_segment = 0
        self.rb_results = rb_results

        return

    def _get_next_id_segment(self):
        """Increment the id of the segment.

        :param None.
        :return: Value of the next ID
        :rtype: Nint
        """

        self._id_qgs_segment += 1

        return self._id_qgs_segment

    def _create_feature_segment(self, origin_id, qgs_geom):
        """Creates a new QgsFeature to load in the QgsSpatialIndex.

        :param int origin_id: ID of the source feature (for reference purpose)
        :return: The feature created
        :rtype: Qgsfeature
        """

        id_segment = self._get_next_id_segment()
        self._dict_qgs_segment[id_segment] = origin_id  #  Creates a reference for to the original RbGeom
        qgs_feature = QgsFeature(id=id_segment)
        qgs_feature.setGeometry(qgs_geom)

        return qgs_feature

    def add_features(self, rb_geoms):
        """Add a RbGeom object in the spatial index.

        For the LineString geometries. The geometry is broken into each line segment that are individually
        loaded in the QgsSpatialIndex.  This strategy takes longer to load than if the feature was loaded as a whole
        but is better for much of the cases.

        :param [RbGeom] rb_geoms: List of RbGeom to load
        :return: None
        :rtype: None
        """

        for rb_geom in rb_geoms:
            qgs_features = []
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.Point:
                qgs_features.append(self._create_feature_segment(rb_geom.id, rb_geom.qgs_geom))
            else:
                qgs_points = rb_geom.qgs_geom.constGet().points()
                for i in range(0,(len(qgs_points)-1)):
                    qgs_geom = QgsGeometry(QgsLineString(qgs_points[i], qgs_points[i+1]))
                    qgs_feature = self._create_feature_segment(rb_geom.id, qgs_geom)
                    qgs_features.append(qgs_feature)

            self._spatial_index.addFeatures(qgs_features)  #  Load all the segment of one RbGeom at the same time

        return

    def _merge_lines(self, qgs_multi_line_string):
        """Merge the line segment to form a continous line.

        In some cases, we can end up with multiple merged line

        :param QgsMultiLineString qgs_multi_line_string: Multi line string to merge together
        :return: The merged line string
        :rtype: List of QgsLineString
        """

        qgs_geom_tmp = QgsGeometry(qgs_multi_line_string.clone())
        if qgs_geom_tmp.isEmpty():
            qgs_geoms_ls = []
        else:
            qgs_geom_merged = qgs_geom_tmp.mergeLines()
            qgs_geoms_ls = qgs_geom_merged.coerceToType(QgsWkbTypes.LineString)  # Creates the list of QgsLineString

        return qgs_geoms_ls

    def get_segment_intersect(self, qgs_geom_id, qgs_rectangle ):
        """Find the feature that intersects the bounding box.

        Once the line string intersecting the bounding box are found are found. They are separated into 2 lists.
        The first one being the line string with the same id (same line) the second one all the others

        :param qgs_geom_id int: ID of the line string that is being simplified
        :return: Two lists  of line string segment. First: Line string with same id; Second all the others
        :rtype: tuple of 2 lists
        """

        qgs_multi_ls_with_itself = QgsMultiLineString()
        qgs_multi_ls_with_others = QgsMultiLineString()
        qgs_points = []
        qgs_rectangle.grow(ReduceBend.ZERO_RELATIVE*100)  # Always increase the b_box to avoid degenerated b_box
        qgs_segment_ids = self._spatial_index.intersects(qgs_rectangle)
        for qgs_segment_id in qgs_segment_ids:
            qgs_geom = self._spatial_index.geometry(qgs_segment_id)
            if qgs_geom.wkbType() == QgsWkbTypes.Point:
                qgs_points.append(qgs_geom)
            else:
                qgs_line_string = qgs_geom.constGet().clone()
                if self._dict_qgs_segment[qgs_segment_id] == qgs_geom_id:
                    qgs_multi_ls_with_itself.addGeometry(qgs_line_string)
                else:
                    qgs_multi_ls_with_others.addGeometry(qgs_line_string)

        # Merge the line_segment
        qgs_geoms_with_itself = self._merge_lines(qgs_multi_ls_with_itself)
        qgs_geoms_with_others = self._merge_lines(qgs_multi_ls_with_others) + qgs_points

        return qgs_geoms_with_itself, qgs_geoms_with_others

    def _delete_segment(self, qgs_pnt0, qgs_pnt1):
        """Delete a line segment in the spatial index based on start/end points.

        To minimise the number of feature returned we search for a very small bounding box located in the middle
        of the line segment.  Usually only one line segment is returned.

        :param qgs_pnt0 QgsPoint: start point of the line segment.
        :param qgs_pnt1 QgsPoint: end point of the line segment.
        :return: None
        """

        qgs_mid_point = QgsGeometryUtils.midpoint(qgs_pnt0, qgs_pnt1)
        qgs_rectangle = qgs_mid_point.boundingBox()
        qgs_rectangle.grow(ReduceBend.ZERO_RELATIVE*100)
        ids = self._spatial_index.intersects(qgs_rectangle)
        for id in ids:
            qgs_geom_line = self._spatial_index.geometry(id)  # Extract geometry
            qgs_pnt_start = qgs_geom_line.vertexAt(0)
            qgs_pnt_end = qgs_geom_line.vertexAt(1)
            #  Check if it's the right geometry
            if qgs_pnt_start.distance(qgs_pnt0) <= ReduceBend.ZERO_RELATIVE and \
               qgs_pnt_end.distance(qgs_pnt1) <= ReduceBend.ZERO_RELATIVE:
                feature = QgsFeature(id=id)
                feature.setGeometry(QgsLineString([qgs_pnt_start, qgs_pnt_end]))
                if (self._spatial_index.deleteFeature(feature)):  # Delete the line segment
                    deleted = True
                    break
                else:
                    raise Exception (QgsProcessingException("Unable to delete entry in QgsSpatialIndex..."))
            else:
                deleted = False

        if not deleted:
            raise Exception(QgsProcessingException("Internal structure corruption..."))

        return

    def delete_vertex(self, rb_geom, v_id_start, v_id_end):
        """Delete a vertex in the line and update the spatial index.

        When a vertex in a line string is deleted.  Two line segments are deleted and one line segment is
        created in the spatial index.  cannot delete the first/last vertex of a line string

        :param rb_geom RbGeom: LineString object to update.
        :param v_id_start int: start of the vertex to delete.
        :param v_id_end int: end of the vertex to delete.
        :return: None
        """

        v_ids_to_del = list(range(v_id_start, v_id_end+1))
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

        return


class RbGeom:
    """Class defining the line used for the bend reduction"""

    __slots__ = ('id', 'original_geom_type', 'is_simplest', 'qgs_geom', 'qgs_rectangle', 'bends', 'nbr_bend_reduced',
                 'is_closed', 'need_pivot')

    _id_counter = 0  # Unique ID counter

    def __init__(self, qgs_abs_geom, original_geom_type):
        """Constructor that initialize a RbGeom object.

        :param RbResults rb_results: Object containing the results and statistics of the execution
        :return: None
        :rtype: None
        """

        self.id = self.next_id()
        self.original_geom_type = original_geom_type
        qgs_geometry = qgs_abs_geom.constGet()
        self.qgs_geom = QgsGeometry(qgs_geometry.clone())
        self.is_simplest = False
        self.is_closed = False
        self.need_pivot = False
        self.bends = []
        self.nbr_bend_reduced = 0
        # Set some variable depending on the attribute of the feature
        if self.original_geom_type == QgsWkbTypes.Point:
            self.is_simplest = True  # A point cannot be simplified
        elif self.original_geom_type == QgsWkbTypes.LineString:
            if qgs_geometry.length() >= ReduceBend.ZERO_RELATIVE:
                if qgs_geometry.isClosed():  # Closed LineString
                    if abs(qgs_geometry.sumUpArea()) > ReduceBend.ZERO_RELATIVE:
                        self.is_closed = True
                        self.need_pivot = True
                    else:
                        self.is_simplest = True # Zero area polygon (degenerated).  Do not try to simplify
            else:
                self.is_simplest = True  # Zero length line (degenerated). Do not try to simplify
        else:
            # It's a polygon
            if abs(qgs_geometry.sumUpArea()) > ReduceBend.ZERO_RELATIVE:
                self.is_closed = True
                self.need_pivot = True
            else:
                self.is_simplest = True  # Zero area polygon. Do not simplify the closed line

    def next_id(self):
        """Get the next counterID.

        :param: QgsMultiLineString qgs_multi_line_string: Multi line string to merge togethernONE
        :return: ID of the RbGeom object
        :rtype: int
        """

        RbGeom._id_counter += 1

        return RbGeom._id_counter

    def get_angles(self):
        """Extract the list of angles of the LineString

        The angle od a vertice is the angle formed by the preceding, the current and the next vertice.
        For an open LineString the start/end vertice do not have angles

        :param: None
        :return: Angle of each vertice of the LineString
        :rtype: [real]
        """

        if self.original_geom_type != QgsWkbTypes.Point:
            qgs_line_string = self.qgs_geom.constGet()
            num_xy = qgs_line_string.numPoints()
            xy = [(qgs_line_string.xAt(i), qgs_line_string.yAt(i)) for i in range(num_xy)]
            if self.is_closed:
                # Add two vertice at the start/end for the circularity of a closed line
                end_xy = xy[-2]  # Not that last because it's the same position as the first vertice
                xy.insert(0, end_xy)

            angles = [QgsGeometryUtils.angleBetweenThreePoints(xy[i-1][0], xy[i-1][1], xy[i][0], xy[i][1],
                                                               xy[i+1][0], xy[i+1][1]) for i in range(1, len(xy)-1)]
        else:
            angles = []

        return angles


class Bend:
    """Define a Bend object which is the reduction goal of this algorithm"""

    __slots__ = ('i', 'j', 'area', 'perimeter', 'adj_area', 'to_reduce', 'qgs_geom_new_subline',
                 'qgs_geom_new_subline_trimmed', 'orientation', 'qgs_geom_bend', 'direction')

    def __init__(self, i, j, qgs_polygon):
        """Constructor that initialize a Bend object.

        :param: int i: start osition of the vertice in the LineString to reduce
        :param: int j: end position of the vertice in the LineString to reduce
        :param: QgsPolygon qgs_polygon: Polygon formed by the bend to reduce
        :return: None
        :rtype: None
        """

        self.i = i
        self.j = j
        self.qgs_geom_bend = QgsGeometry(qgs_polygon.clone())
        self.area = self.qgs_geom_bend.area()
        self.perimeter = self.qgs_geom_bend.length()
        self.adj_area = ReduceBend.calculate_adj_area(self.area, self.perimeter)
        self.to_reduce = False
        self.orientation = None
        self.qgs_geom_new_subline = None
        self.qgs_geom_new_subline_trimmed = None

    def get_new_subline(self, rb_geom):
        """Create the new line that will replace the bend.

       :param: RbGeom rb_geom: Geometry containing the lien to reduce
       :return: The new line to close the polygon
       :rtype: QgsGeometry
       """

        if self.qgs_geom_new_subline is None:
            qgs_pnt_i = rb_geom.qgs_geom.vertexAt(self.i)
            qgs_pnt_j = rb_geom.qgs_geom.vertexAt(self.j)
            qgs_ls_new_subline = QgsLineString([qgs_pnt_i, qgs_pnt_j])
            self.qgs_geom_new_subline = QgsGeometry(qgs_ls_new_subline)

        return self.qgs_geom_new_subline

    def get_new_subline_trimmed(self, rb_geom):
        """Create the new line that will replace the bend but just a little bit shorter.

       :param: RbGeom rb_geom: Geometry containing the lien to reduce
       :return: The new line to close the polygon
       :rtype: QgsGeometry
       """

        if self.qgs_geom_new_subline_trimmed is None:
            qgs_ls_new_line = self.get_new_subline(rb_geom).constGet()
            if qgs_ls_new_line.length() >= ReduceBend.ZERO_RELATIVE*100.:
                qgs_pnt_i_trimmed = qgs_ls_new_line.interpolatePoint(ReduceBend.ZERO_RELATIVE)
                qgs_pnt_j_trimmed = qgs_ls_new_line.interpolatePoint(qgs_ls_new_line.length() - ReduceBend.ZERO_RELATIVE)
                qgs_ls_new_subline_trimmed = QgsLineString([qgs_pnt_i_trimmed, qgs_pnt_j_trimmed])
            else:
                qgs_ls_new_subline_trimmed = qgs_ls_new_line.clone()
            self.qgs_geom_new_subline_trimmed = QgsGeometry(qgs_ls_new_subline_trimmed)


        return self.qgs_geom_new_subline_trimmed


class RbResults:
    """Class defining the stats and result"""

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'nbr_bend_reduced', 'nbr_bend_detected', \
                 'qgs_features_out', 'nbr_hole_del', 'nbr_pol_del', 'nbr_pass', 'is_structure_valid', \
                 'lines_log_info')

    def __init__(self):
        """Constructor that initialize a RbResult object.

        :param: None
        :return: None
        :rtype: None
        """

        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_bend_reduced = 0
        self.nbr_bend_detected = 0
        self.qgs_features_out = None
        self.nbr_hole_del = 0
        self.nbr_pol_del = 0
        self.nbr_pass = 0
        self.is_structure_valid = None
        self.lines_log_info = []


class Epsilon:
    """Class defining the value of the"""

    __slots__ = '_zero_relative', '_zero_absolute', '_zero_angle', '_map_range'

    def __init__(self, features):
        """Constructor that initialize the Epsilon (near zero) object.

        The dynamic (range) of the feature can vary a lot. We calculate the dynamic of the bounding box of all
        the features and we use it to estimate an epsilon (zero).  when the range of the bounding box is very small
        the epsilon can be very small and the opposite when the bigger the boundoing box is.

        :param: [QgsFeatures] feaures: List of QgsFeature to process.
        :return: None
        :rtype: None
        """

        if len(features) >= 1:
            b_box = features[0].geometry().boundingBox()  #  Initialize the bounding box
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
        abs_digit = max_digit  - security
        rel_digit = max_digit - log_loss - security
        self._zero_relative = (1. / (10**(rel_digit)))
        self._zero_absolute = 1. / (10**(abs_digit))
        self._zero_angle= math.radians(.0001)  # Angle used to decide a flat angle


    def set_class_variables(self):
        """Set the different epsilon values.

       :param: None
       :return: None
       :rtype: None
       """

        ReduceBend.ZERO_RELATIVE =  self._zero_relative
        ReduceBend.ZERO_ABSOLUTE = self._zero_absolute
        ReduceBend.ZERO_ANGLE = self._zero_angle

        return


class ReduceBend():
    """Main class for bend reduction"""

    @staticmethod
    def reduce(qgs_in_features, diameter_tol, feedback=None, flag_del_outer=False,
               flag_del_inner=False, validate_structure=False):
        """Main static method used to launch the bend reduction.

       :params: [QgsFeatures] qgs_features: List of geatures to process.
       :params: real diameter_tol: Tolerance of the diameter of the bend to reduce.
       :params: QgsFeedback feedback: Handle for interaction with QGIS.
       :params: Bool flag_del_outer: Delete polygon if area below the diameter tolerance.
       :params: Bool flag_del_inner: Delete polygon holes if area below the diameter tolerance.
       :params: Bool validate_structure: Validate internal data structure after processing (for debugging)
       :return: Statistics and result object.
       :rtype: RbResult
       """

        rb = ReduceBend(qgs_in_features, diameter_tol, feedback, flag_del_outer, flag_del_inner, validate_structure)
        results = rb.reduce_bends()

        return results

    @staticmethod
    def extract_polygon_attributes(qgs_geom):
        """Static method to calculate the area and perimeter of a LineString.

       :params: QgsGeometry qgs_geom: Geometry to process.
       :return: Area and perimeter of the geometry.
       :rtype: Tuple
       """

        qgs_line_string = qgs_geom.constGet()
        qgs_pol = QgsPolygon(qgs_line_string.clone())
        area = qgs_pol.area()
        perimeter = qgs_pol.perimeter()

        return (area, perimeter)

    @staticmethod
    def calculate_adj_area(area, perimeter):
        """Static method to calculate the adjusted area.

        The adjusted area is used to determine if a bend must be reduce.

       :params: real area: area of a polygon.
       :params: real perimeter: perimeter of a polygon.
       :return: Adjusted area of a polygon
       :rtype: Real
       """

        compactness_index = 4 * area * math.pi / perimeter ** 2
        adj_area = area * (.75 / compactness_index)

        return adj_area

    @staticmethod
    def calculate_min_adj_area(diameter_tol):
        """Static method to calculate the adjusted area of the maximum diameter tolerance.

       :params: real diameter: Diameter tolerance to used for bend reduction
       :return: Minimum adjusted area of a polygon to reduce
       :rtype: Real
       """

        min_adj_area = .75 * math.pi * (diameter_tol / 2.) ** 2

        return min_adj_area

    __slots__ = ('qgs_in_features', 'diameter_tol', 'feedback', 'flag_del_outer', 'flag_del_inner',
                 'validate_structure', 'rb_collection', 'eps', 'rb_results', 'rb_features', 'rb_geoms')

    def __init__(self, qgs_in_features, diameter_tol, feedback, flag_del_outer, flag_del_inner, validate_structure):
        """Constructor for the bend reduction.

       :params: [QgsFeatures] qgs_features: List of geatures to process.
       :params: real diameter_tol: Tolerance of the diameter of the bend to reduce.
       :params: QgsFeedback feedback: Handle for interaction with QGIS.
       :params: Bool flag_del_outer: Delete polygon if area below the diameter tolerance.
       :params: Bool flag_del_inner: Delete polygon holes if area below the diameter tolerance.
       :params: Bool validate_structure: Validate internal data structure after processing (for debugging)
       :return: None
       :rtype: None
       """

        self.qgs_in_features = qgs_in_features
        self.diameter_tol = diameter_tol
        self.feedback = feedback
        self.flag_del_outer = flag_del_outer
        self.flag_del_inner = flag_del_inner
        self.validate_structure = validate_structure
        self.rb_collection = None

    def reduce_bends(self):
        """Main method to manage bend reduction.

        :params: None
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

        # Recreate the QgsFeature
        qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in self.rb_features]

        # Set return values
        self.rb_results.out_nbr_features = len(qgs_features_out)
        self.rb_results.qgs_features_out = qgs_features_out

        # Validate inner spatial structure. For debug purpose only
        if self.validate_structure:
            self.rb_results.is_structure_valid = self.validate_integrity()

        #  Code used for the profiler (uncomment if needed)
#        pr.disable()
#        s = io.StringIO()
#        sortby = SortKey.CUMULATIVE
#        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#        ps.print_stats()
#        print(s.getvalue())

        return self.rb_results

    def create_rb_feature(self):
        """Create the different RbFeatures from the QgsFeatures.

        :params: None
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
                # Manage usupported geometry type
                raise Exception("Unsupported GeometryType: {}".format(qgs_geom.wkbType()))

        return rb_features

    def pre_reduction_process(self):
        """This method execute the pre reduction process

        Pre reduction process includes remove small polygon or polygon hole simplify line when vertice
        are really too close and transform the RbFeature into LineString

        :params: None
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


        # Remove co-linear and duplicate nodes
        for rb_geom in rb_geoms:
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
                if not rb_geom.is_simplest:
                    rb_geom.qgs_geom.removeDuplicateNodes(epsilon=ReduceBend.ZERO_RELATIVE)
                    self.delete_co_linear(rb_geom)

        return rb_geoms


    def del_outer_inner_ring(self):
        """This method deletes the polygons and polygon holes below the diameter tolerance

        :params: None
        :return: None
        :rtype: None
        """

        # Loop over each rb_features
        for i in reversed(range(len(self.rb_features))):  # List visited in reverse order for easier deletion of entry
            if isinstance(self.rb_features[i], RbPolygon):  # Only process Polygon
                min_adj_area = ReduceBend.calculate_min_adj_area(self.diameter_tol)
                for j in reversed(range(len(self.rb_features[i].rb_geom))):  # List also visited in reverse order
                    area, perimeter = ReduceBend.extract_polygon_attributes(self.rb_features[i].rb_geom[j].qgs_geom)
                    adj_area = self.calculate_adj_area(area, perimeter)
                    if j == 0:
                        # Process the exterior ring (outer ring always at position 0)
                        if self.flag_del_outer:
                            if adj_area < min_adj_area:
                                del self.rb_features[i]  # Delete the rb_feature
                                self.rb_results.nbr_pol_del += 1
                                break
                    else:
                        # Process an interior ring
                        if self.flag_del_inner:
                            if adj_area < min_adj_area:
                                del self.rb_features[i].rb_geom[j]  # Delete the ring
                                self.rb_results.nbr_hole_del += 1

        return

    def set_bend_direction(self, rb_geom):
        """For closed LineString, this method set the direction of each bend.

        For a closed LineString each bend is going either inside or outside.  Once a bend is determined to go
        inside or outside all the other bend alternate from inside to outside or the opposite

        :params: RbGeom rb_geom: Geometry to set bend direction
        :return: None
        :rtype: None
        """

        if rb_geom.is_closed and len(rb_geom.bends) >= 1:
            # Detemine if the first bend is going inside or outside the polygon
            first_bend = rb_geom.bends[0]
            qgs_geom_centroid = first_bend.qgs_geom_bend.centroid()
            qgs_geom_pol = QgsGeometry(QgsPolygon(rb_geom.qgs_geom.constGet().clone(), []))
            if qgs_geom_centroid.within(qgs_geom_pol):
                rb_geom.bends[0].direction = 'IN'
            else:
                rb_geom.bends[0].direction = 'OUT'

            #  All the other bend alternate [in, out, in, out] or [out, in, out, in...]
            for i in range(1, len(rb_geom.bends)):
                if rb_geom.bends[i - 1].direction == 'IN':
                    rb_geom.bends[i].direction = 'OUT'
                else:
                    rb_geom.bends[i].direction = 'IN'


    def pivot_closed_line(self, rb_geom, diameter_tol):

        if rb_geom.is_closed and rb_geom.need_pivot:
            bend_location = None
            bend_area = 0.0
            for bend in rb_geom.bends:
                if bend.area > bend_area:
                    bend_location = bend
                    bend_area = bend.area
                if bend.j - bend.i >= 4:
                    if bend.area >= ReduceBend.calculate_min_adj_area(diameter_tol):
                        bend_location = bend
                        if bend.direction == 'OUT':
                            rb_geom.need_pivot = False  # Optimal bend found
                            break

            if bend_location is not None:
                # There is bend candidate for a rotation
                qgs_points = rb_geom.qgs_geom.constGet().points()
                new_start_end = (bend_location.j + bend_location.i) // 2
                new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end + 1]
                rb_geom.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))

        return


    def _manage_reduce_bend(self):

        rb_geoms_done = []
        nbr_pass = 0
        min_pass = 1
        nbr_geoms = 100.0 / len(self.rb_geoms) if len(self.rb_geoms) >= 1 else 0
        while True:
            current_diameter_tol = self.diameter_tol * min(nbr_pass+1, min_pass)/float(min_pass)
            self.remove_rb_geoms_done(rb_geoms_done)  # Remove feature done to accelerate process
            if nbr_pass == 0:
                self.feedback.setProgress(1) # Set the progress bar
            else:
                self.feedback.setProgress(max(1, int(len(rb_geoms_done)) * nbr_geoms))
            nbr_bend_reduced = 0
            nbr_bend_detected = 0
            for rb_geom in self.rb_geoms:
                if self.feedback.isCanceled(): break
                nbr_bend_detected += self.manage_bend_creation(nbr_pass, rb_geom, self.diameter_tol)
                self.flag_bend_to_reduce(rb_geom, current_diameter_tol)
                nbr_bend_reduced += self.process_bends(rb_geom)

            str = "Iteration: {}; Bends detected: {}; Bend reduced: {}; Tolerance used: {}"\
                  .format(nbr_pass, nbr_bend_detected, nbr_bend_reduced, current_diameter_tol)
            self.rb_results.lines_log_info.append(str)
            if nbr_pass == 0: self.rb_results.nbr_bend_detected = nbr_bend_detected

            # While loop breaking condition
            if nbr_pass > min_pass and nbr_bend_reduced == 0:
                break
            nbr_pass += 1

        # Reset the rb_geoms list
        self.rb_geoms += rb_geoms_done
        self.rb_results.nbr_pass = nbr_pass

        return

    def remove_rb_geoms_done(self, rb_geoms_done):

        for i in reversed(range(len(self.rb_geoms))):
            if self.rb_geoms[i].is_simplest:
                rb_geoms_done.append(self.rb_geoms[i])
                del self.rb_geoms[i]

        return


    def delete_co_linear(self, rb_geom):
        # Delete co-linear vertice with angle of near 0 or 180 degrees
#        qgs_line_string = rb_geom.qgs_geom.constGet()
#        qgs_points = qgs_line_string.points()
#        num_points = len(qgs_points)
#        num_points_remaining = num_points
#        if rb_geom.original_geom_type == QgsWkbTypes.LineString:
#            min_remaining_points = 2  # Minimum number of vertice for a LineString
#        else:
#            min_remaining_points = 4  # Minimum number of vertice for a Polygon
        vertex_ids_to_del = []
        angles = rb_geom.get_angles()
        if rb_geom.is_closed and len(angles) >= 1:
            del angles[0]  # Do not process the start/end points (even if co-linear)
        for i, angle in enumerate(angles):
            if abs(angle - math.pi) <= ReduceBend.ZERO_ANGLE or abs(angle) <= ReduceBend.ZERO_ANGLE:
                # Co-linear point or flat angle delete the current point
                vertex_ids_to_del.append(i+1)
#                num_points_remaining -= 1
#        i = num_points - 2
#        p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
#        p2_x, p2_y = qgs_points[i + 1].x(), qgs_points[i + 1].y()
#        while i >= 1 and num_points_remaining >= min_remaining_points:
#            p0_x, p0_y = qgs_points[i - 1].x(), qgs_points[i - 1].y()
#            angle = QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
#            if abs(angle - math.pi) <= ReduceBend.ZERO_ANGLE or abs(angle) <= ReduceBend.ZERO_ANGLE:
#                # Co-linear point or flat angle delete the current point
#                vertex_ids_to_del.append(i)
#                num_points_remaining -= 1
#
#            i -= 1
#            p2_x, p2_y = p1_x, p1_y
#            p1_x, p1_y = p0_x, p0_y

        # Delete co-linerar vertex
        for vertex_id_to_del in reversed(vertex_ids_to_del):
            if self.rb_collection is None:
                rb_geom.qgs_geom.deleteVertex(vertex_id_to_del)
            else:
                self.rb_collection.delete_vertex(rb_geom, vertex_id_to_del, vertex_id_to_del)

        if rb_geom.qgs_geom.length() <= ReduceBend.ZERO_RELATIVE:
            # Something wrong.  do not try to simplify the LineString
            rb_geom.is_simplest = True

        return

    def flag_bend_to_reduce(self, rb_geom, diameter_tol):
        # Minimum adjusted area used to find bend to reduce
        min_adj_area = ReduceBend.calculate_min_adj_area(diameter_tol)

        if rb_geom.is_closed and len(rb_geom.bends) >= 3:
            # The closed line start/end point lie on a bend that do not need to be reduced
            del rb_geom.bends[0]  # Remove the first bend
            del rb_geom.bends[-1]  # Remove the last bend

#        if rb_geom.is_closed:
#            for bend in rb_geom.bends:
#                if bend.direction == "OUT": bend.area *= .66
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


    def validate_spatial_constraints(self, ind, rb_geom):

        check_constraints = True
        bend = rb_geom.bends[ind]
        qgs_geom_new_subline = bend.get_new_subline(rb_geom)

        b_box = bend.qgs_geom_bend.boundingBox()
        qgs_geoms_with_itself, qgs_geoms_with_others = self.rb_collection.get_segment_intersect(rb_geom.id, b_box)

        # First: check if the bend reduce line string is an OGC simple line
        # We test with a tiny smaller line to ease the testing and false positive error
        if check_constraints:
            if qgs_geom_new_subline.length() >= ReduceBend.ZERO_RELATIVE:
                qgs_geom_new_subline_trimmed = bend.get_new_subline_trimmed(rb_geom)
#                qgs_geom_potentials = self.rb_collection.get_line_segments(rb_geom.id,
#                                                                      qgs_geom_new_subline_trimmed.boundingBox())
                for qgs_geom_potential in qgs_geoms_with_itself:
                    if qgs_geom_new_subline_trimmed.disjoint(qgs_geom_potential):
                        #            if not qgs_geom_new_sub_trim_engine.disjoint(qgs_geom_potential.constGet()):
                        # Everything is OK
                        pass
                    else:
                        # The new sub line intersect the line itself. The result would create a non OGC simple line
                        check_constraints = False
                        break
            else:
                qgs_line_string = qgs_geom_new_subline.constGet()
                x = qgs_line_string.startPoint().x()
                y = qgs_line_string.startPoint().y()
                text = "Possibly non OGC simple feature at {},{} use Fix geometries".format(x,y)
                self.feedback.pushInfo(text)


        # Second: check that the new line does not intersect any other line or points
        if check_constraints:
            qgs_rectangle = bend.qgs_geom_bend.boundingBox()
            #        qgs_geom_engine_new_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet())
#            qgs_geom_potentials = self.rb_collection.get_features(qgs_rectangle, [rb_geom.id])
            for qgs_geom_potential in qgs_geoms_with_others:
                if not qgs_geom_potential.disjoint(qgs_geom_new_subline):
                    #            if not qgs_geom_new_subline.disjoint(qgs_geom_potential):
                    #            if not qgs_geom_engine_new_subline.disjoint(qgs_geom_potential.constGet()):
                    # The bend area intersects with a point
                    check_constraints = False
                    break

        # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
        # sidedness or relative position error
        if check_constraints:
            #        qgs_geom_engine_bend_area = QgsGeometry.createGeometryEngine(bend.qgs_geom_bend)
            #        qgs_geom_bend = QgsGeometry(bend.qgs_geom_bend.clone())
            for qgs_geom_potential in qgs_geoms_with_others:
                #            if qgs_geom_bend_area.contains(qgs_geom_potential.constGet()):
                if bend.qgs_geom_bend.contains(qgs_geom_potential):
                    # A feature is totally located inside
                    check_constraints = False
                    break

        return check_constraints

    def process_bends(self, rb_geom):

        nbr_bend_reduced = 0
        for ind in reversed(range(len(rb_geom.bends))):
            bend = rb_geom.bends[ind]
            if bend.to_reduce:
                # Check spatial constraints
                spatial_constraints = self.validate_spatial_constraints(ind, rb_geom)
                if spatial_constraints:
#                    v_ids_to_del = list(range(bend.i + 1, bend.j))  # List of the vertex id to delete
                    self.rb_collection.delete_vertex(rb_geom, bend.i+1, bend.j-1)
                    #                bend.reduce(rb_geom)
                    self.rb_results.nbr_bend_reduced += 1  # Global counter of bend reduced
                    nbr_bend_reduced += 1  # Local counter of bend reduced

        return nbr_bend_reduced

    def create_polygon (self, i, j, qgs_points):

        # Create the list of point to create
        if i < j:
            index = list(range(i, j+1)) + [i]
        else:
            index = list(range(i, len(qgs_points))) + list(range(0,j+1)) + [i]  # Manage circular array

        qgs_sub_points = [qgs_points[k] for k in index]
        qgs_polygon = QgsPolygon(QgsLineString(qgs_sub_points))

        return qgs_polygon


    def manage_bend_creation(self, nbr_pass, rb_geom, diameter_tol):

        self.delete_co_linear(rb_geom)
        nbr_bends_detected = self.detect_bends(rb_geom)
        if rb_geom.is_closed:
            if nbr_pass >= 2 or len(rb_geom.bends) < 10:
                self.pivot_closed_line(rb_geom, diameter_tol)
                self.delete_co_linear(rb_geom)
                nbr_bends_detected = self.detect_bends(rb_geom)

        return nbr_bends_detected

    def detect_bends(self, rb_geom):

        rb_geom.bends = []  # Reset the list of bends
        angles = rb_geom.get_angles()
        # Modify the angle to binary orientation: clockwise or anti clockwise
        orientation = [CLOCK_WISE if angle >= math.pi else ANTI_CLOCK_WISE for angle in angles]
        if rb_geom.is_closed:
            if len(set(orientation)) == 1:
                orientation = []  # All the angles have the same orientation.  No bend to reduce
            else:
                del orientation[0]  # Do not process the fist angle as it is the angle of the start/end

        if len(orientation) >= 1:
            orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
            orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

        #            if rb_geom.is_closed:
#                # Copy the first orientation to the last position for the circularity
#                orientation.append(orientation[0])
#            else:
#                orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
#                orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

        # Find the inflexion points in the line.
        inflexion = [i for i in range(0, len(orientation)-1) if orientation[i] != orientation[(i + 1)]]
        if len(inflexion) != 0:
#            if rb_geom.is_closed:
#                inflexion.append(inflexion[0])

            qgs_points = rb_geom.qgs_geom.constGet().points()
            for k in range(len(inflexion)-1):
                i = inflexion[k]
                j = (inflexion[(k+1)]+1)
                qgs_polygon = self.create_polygon(i, j, qgs_points)
                rb_geom.bends.append(Bend(i, j, qgs_polygon))

        else:
            # If there is no inflexion the line cannot be simplified
            rb_geom.is_simplest = True

        # Set the direction of the bend inside or outside the polygon
        if rb_geom.is_closed:
            self.set_bend_direction(rb_geom)
        #        if len(rb_geom.bends) == 0:
#            # A line with no inflexion cannot be simplified more
#            rb_geom.is_simplest = True

#        if rb_geom.is_closed:
#            # Process the line string as a closed line
#            angle = QgsGeometryUtils.angleBetweenThreePoints(xy[-1][0],xy[-1][1],xy[0][0],xy[0][1],xy[1][0],xy[1][1])
#            # Add a first orientation linking the start and end of the line (circular array)
#            orientation.insert(0, CLOCK_WISE) if angle >= math.pi else orientation.insert(0, ANTI_CLOCK_WISE)
#            # Find the inflexion points in the line.  Managing the circular array
#            inflexion = [i for i in range(0, len(orientation)) if orientation[i] != orientation[(i+1)%len(orientation)]]
#            for k in range(len(inflexion)):
#                i = inflexion[k]
#                j = inflexion[(k+1)%len(inflexion)]+1
#                if i != j:
#                    qgs_polygon = self.create_polygon(i, j, xy)
#                    rb_geom.bends.append(Bend(i, j, qgs_polygon))
#                else:
#                    # Extemely rare case. Happen with special case polygon
#                    pass

#            if 2 <= len(rb_geom.bends) <= 4:
#                del rb_geom.bends[-1]
#                del rb_geom.bends[-1]
#        else:
#            # For open lines add a point at the start and one at the end (to facilitate )
#            if len(orientation) >= 1:
#                orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
#                orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

#            # Find the inflexion points in the line
#            inflexions = [i for i in range(0,len(orientation)-1) if orientation[i] != orientation[i+1]]

#            for k in range(len(inflexions)-1):
#                i = inflexions[k]
#                j = inflexions[k+1]+1
#                qgs_polygon = self.create_polygon(i, j, xy)
#                rb_geom.bends.append(Bend(i, j, qgs_polygon))

        return len(rb_geom.bends)


    def validate_integrity(self):

        is_structure_valid = True
        for rb_geom in self.rb_geoms:
            qgs_line_string = rb_geom.qgs_geom.constGet()
            if qgs_line_string.wkbType() == QgsWkbTypes.LineString:
                qgs_points = qgs_line_string.points()
                for i in range(len(qgs_points)-1):
                    self.rb_collection._delete_segment(qgs_points[i], qgs_points[i+1])

        if is_structure_valid:
            qgs_rectangle = QgsRectangle(-sys.float_info.max, -sys.float_info.max,
                                         sys.float_info.max,sys.float_info.max)
            ids = self.rb_collection._spatial_index.intersects(qgs_rectangle)
            for id in ids:
                qgs_geom = self.rb_collection._spatial_index.geometry(id)
                if qgs_geom.wkbType() == QgsWkbTypes.Point:
                    pass
                else:
                    # Error
                    is_structure_valid = False

        return is_structure_valid
