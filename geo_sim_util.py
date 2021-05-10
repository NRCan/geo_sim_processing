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

import math
import sys
from abc import ABC, abstractmethod
from qgis.core import (QgsLineString, QgsWkbTypes, QgsSpatialIndex, QgsGeometry, QgsPolygon,
                       QgsGeometryUtils, QgsRectangle, QgsProcessingException)


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
        max_digit = 15  # Number of meaningful digits for real number
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


class ProgressBar:
    """Class used for managing the progress bar in the QGIS desktop

    """

    def __init__(self, feedback, max_value, message=None):
        """Constructor of the ProgressBar class

        :param: feedback: feedback handle for interaction with the QGIS desktop
        :param: max_value: Integer of the maximum value """

        self.feedback = feedback
        self.max_value  = max_value
        self.progress_bar_value = 0
        self.feedback.setProgress(self.progress_bar_value)
        if message is not None or message != "":
            self.feedback.pushInfo(message)

    def set_value(self, value):
        """Set the value of the progress bar

        :param: value: Integer value to use to set the progress bar

        """

        percent_value = int(value/self.max_value*100.)
        if percent_value != self.progress_bar_value:
            self.progress_bar_value = percent_value
            self.feedback.setProgress(self.progress_bar_value)




class GsCollection:
    """Class used for managing the QgsFeature spatially.

    QgsSpatialIndex class is used to store and retrieve the features.
    """

    __slots__ = ('_spatial_index', '_dict_qgs_segment', '_id_qgs_segment')

    def __init__(self):
        """Constructor that initialize the GsCollection.

        """

        self._spatial_index = QgsSpatialIndex()
        self._dict_qgs_segment = {}  # Contains a reference to the original geometry
        self._id_qgs_segment = 0

    def _get_next_id_segment(self):
        """Increment the id of the segment.

        :return: Value of the next ID
        :rtype: int
        """

        self._id_qgs_segment += 1

        return self._id_qgs_segment

    def _create_rectangle(self, geom_id, qgs_geom):
        """Creates a new QgsRectangle to load in the QgsSpatialIndex.

        :param: geom_id: Integer ID of the geometry
        :param: qgs_geom: QgsGeometry to use for bounding box extraction
        :return: The feature created
        :rtype: QgsFeature
        """

        id_segment = self._get_next_id_segment()
        self._dict_qgs_segment[id_segment] = (geom_id, qgs_geom)  # Reference to the RbGeom ID and geometry

        return id_segment, qgs_geom.boundingBox()

    def add_features(self, rb_geoms, feedback):
        """Add a RbGeom object in the spatial index.

        For the LineString geometries. The geometry is broken into each line segment that are individually
        loaded in the QgsSpatialIndex.  This strategy accelerate the validation of the spatial constraints.

        :param: rb_geoms: List of RbGeom to load in the QgsSpatialIndex
        :feedback: QgsFeedback handle used to update the progress bar
        """

        progress_bar = ProgressBar(feedback, len(rb_geoms), "Building internal structure...")
        for val, rb_geom in enumerate(rb_geoms):
            progress_bar.set_value(val)
            qgs_rectangles = []
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.Point:
                qgs_rectangles.append(self._create_rectangle(rb_geom.id, rb_geom.qgs_geom))
            else:
                qgs_points = rb_geom.qgs_geom.constGet().points()
                for i in range(0, (len(qgs_points)-1)):
                    qgs_geom = QgsGeometry(QgsLineString(qgs_points[i], qgs_points[i+1]))
                    qgs_rectangles.append(self._create_rectangle(rb_geom.id, qgs_geom))

            for geom_id, qgs_rectangle in qgs_rectangles:
                self._spatial_index.addFeature(geom_id, qgs_rectangle)

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
        ids = self._spatial_index.intersects(qgs_rectangle)
        for geom_id in ids:
            target_qgs_geom_id, target_qgs_geom = self._dict_qgs_segment[geom_id]
            if target_qgs_geom_id is None:
                # Nothing to do; segment was deleted
                pass
            else:
                if target_qgs_geom_id == qgs_geom_id:
                    # Test that the segment is not part of qgs_subline
                    if not target_qgs_geom.within(qgs_geom_subline):
                        qgs_geoms_with_itself.append(target_qgs_geom)
                else:
                    qgs_geoms_with_others.append(target_qgs_geom)

        return qgs_geoms_with_itself, qgs_geoms_with_others

    def _delete_segment(self, qgs_geom_id, qgs_pnt0, qgs_pnt1):
        """Delete a line segment in the spatial index based on start/end points.

        To minimise the number of feature returned we search for a very small bounding box located in the middle
        of the line segment.  Usually only one line segment is returned.

        :param qgs_geom_id: Integer ID of the geometry
        :param qgs_pnt0 : QgsPoint start point of the target line segment.
        :param qgs_pnt1 : QgsPoint end point of the target line segment.
        """

        qgs_geom_to_delete = QgsGeometry(QgsLineString(qgs_pnt0, qgs_pnt1))
        qgs_mid_point = QgsGeometryUtils.midpoint(qgs_pnt0, qgs_pnt1)
        qgs_rectangle = qgs_mid_point.boundingBox()
        qgs_rectangle.grow(Epsilon.ZERO_RELATIVE*100)
        deleted = False
        ids = self._spatial_index.intersects(qgs_rectangle)
        for geom_id in ids:
            target_qgs_geom_id, target_qgs_geom = self._dict_qgs_segment[geom_id]  # Extract id and geometry
            if qgs_geom_id == target_qgs_geom_id:
                # Only check for the same ID
                if target_qgs_geom.equals(qgs_geom_to_delete):  # Check if it's the same geometry
                    deleted = True
                    self._dict_qgs_segment[geom_id] = (None, None)  # Delete from the internal structure
                    break

        if not deleted:
            raise Exception(QgsProcessingException("Internal structure corruption..."))

        return

    def _delete_vertex(self, rb_geom, v_id_start, v_id_end):
        """Delete consecutive vertex in the line and update the spatial index.

        When a vertex in a line string is deleted.  Two line segments are deleted and one line segment is
        created in the spatial index.  Cannot delete the first/last vertex of a line string

        :param rb_geom: LineString object to update.
        :param v_id_start: start of the vertex to delete.
        :param v_id_end: end of the vertex to delete.
        """

        is_closed = rb_geom.qgs_geom.constGet().isClosed()
        v_ids_to_del = list(range(v_id_start, v_id_end+1))
        if v_id_start == 0 and is_closed:
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
            self._delete_segment(rb_geom.id, qgs_pnt0, qgs_pnt1)

        # Add the new line segment in the spatial index
        qgs_pnt0 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[0])
        qgs_pnt1 = rb_geom.qgs_geom.vertexAt(v_ids_to_del[-1])
        qgs_geom_segment = QgsGeometry(QgsLineString(qgs_pnt0, qgs_pnt1))
        geom_id, qgs_rectangle = self._create_rectangle(rb_geom.id, qgs_geom_segment)
        self._spatial_index.addFeature(geom_id, qgs_rectangle)

        # Delete the vertex in the line string geometry
        for v_id_to_del in reversed(range(v_id_start, v_id_end+1)):
            rb_geom.qgs_geom.deleteVertex(v_id_to_del)
            if v_id_start == 0 and is_closed:
                # Special case for closed line where we simulate a circular array
                nbr_vertice = rb_geom.qgs_geom.constGet().numPoints()
                qgs_pnt_first = rb_geom.qgs_geom.vertexAt(0)
                rb_geom.qgs_geom.insertVertex(qgs_pnt_first, nbr_vertice-1)
                rb_geom.qgs_geom.deleteVertex(nbr_vertice)

        return

    def delete_vertex(self, rb_geom, v_id_start, v_id_end):
        """Manage deletion of consecutives vertex.

        If v_id_start is greater than v_id_end the delete is broken into up to 3 calls

        :param rb_geom: LineString object to update.
        :param v_id_start: start of the vertex to delete.
        :param v_id_end: end of the vertex to delete.
        """

        num_points = rb_geom.qgs_geom.constGet().numPoints()
        # Manage closes line where first/last vertice are the same
        if v_id_start == num_points-1:
            v_id_start = 0  # Last point is the same as the first vertice
        if v_id_end == -1:
            v_id_end = num_points -2  # Preceding point the first/last vertice

        if v_id_start <= v_id_end:
            self._delete_vertex(rb_geom, v_id_start, v_id_end)
        else:
            self._delete_vertex(rb_geom, v_id_start, num_points-2)
            self._delete_vertex(rb_geom, 0, 0)
            if v_id_end > 0:
                self._delete_vertex(rb_geom, 1, v_id_end)
#            lst_vertex_to_del = list(range(v_id_start, num_points)) + list(range(0, v_id_end+1))
#            for vertex_to_del in lst_vertex_to_del:
#                self._delete_vertex(rb_geom, vertex_to_del, vertex_to_del)

#        num_points = rb_geom.qgs_geom.constGet().numPoints()
#        lst_vertex_to_del = list(range(v_id_start, num_points)) + list(range(0, v_id_end + 1))
#        for vertex_to_del in lst_vertex_to_del:
#            self._delete_vertex(rb_geom, vertex_to_del, vertex_to_del)

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
        self._delete_segment(rb_geom.id, qgs_pnt0, qgs_pnt1)

        qgs_points = qgs_geom_new_subline.constGet().points()
        tmp_qgs_points = qgs_points[1:-1]  # Drop first/last item
        # Insert the new vertex in the QgsGeometry. Work reversely to facilitate insertion
        for qgs_point in reversed(tmp_qgs_points):
            rb_geom.qgs_geom.insertVertex(qgs_point, bend_j)

        # Add the new segment in the spatial container
        for i in range(len(qgs_points)-1):
            qgs_geom_segment = QgsGeometry(QgsLineString(qgs_points[i], qgs_points[i+1]))
            geom_id, qgs_rectangle = self._create_rectangle(rb_geom.id, qgs_geom_segment)
            self._spatial_index.addFeature(geom_id, qgs_rectangle)

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
                    self._delete_segment(rb_geom.id, qgs_points[i], qgs_points[i+1])

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


class GsFeature(ABC):
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

    @staticmethod
    def create_gs_feature(qgs_in_features):
        """Create the different GsFeatures from the QgsFeatures.

        :param: qgs_in_features: List of QgsFeature to process
        :return: List of rb_features
        :rtype: [GsFeature]
        """

        rb_features = []

        for qgs_feature in qgs_in_features:
            qgs_geom = qgs_feature.geometry()  # extract the Geometry

            if GsFeature.is_polygon(qgs_geom.wkbType()):
                rb_features.append(GsPolygon(qgs_feature))
            elif GsFeature.is_line_string(qgs_geom.wkbType()):
                rb_features.append(GsLineString(qgs_feature))
            elif GsFeature.is_point(qgs_geom.wkbType()):
                rb_features.append(GsPoint(qgs_feature))
            else:
                raise QgsProcessingException("Internal geometry error")

        return rb_features

    def __init__(self, qgs_feature):
        """Constructor of the GsFeature class.

        :param qgs_feature: QgsFeature to process.
        """

        self.qgs_feature = qgs_feature
        self.id = GsFeature._id_counter
        GsFeature._id_counter += 1
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


class GsPolygon(GsFeature):
    """Class description for GsPolygon"""

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


class GsLineString(GsFeature):
    """Class managing a GsLineString.
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


class GsPoint(GsFeature):
    """Class managing a GsPoint
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


class RbGeom:
    """Class defining the line string used for the bend reduction"""

    __slots__ = ('id', 'original_geom_type', 'is_simplest', 'qgs_geom', 'bends', 'need_pivot')

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
        self.bends = None
        # Set some variable depending on the geometry of the feature
        if self.original_geom_type == QgsWkbTypes.Point:
            self.is_simplest = True  # A point cannot be simplified
        else:
            # Attribute setting for LineString and Polygon
            if qgs_geometry.length() >= Epsilon.ZERO_RELATIVE:
                if qgs_geometry.isClosed():
                    self.need_pivot = True  # A closed lined string can be pivoted
            else:
                self.is_simplest = True  # Degenerated LineString... Do not try to simplify...


class SimGeom:
    """Class defining the line string used for the douglas peucker simplification"""

    __slots__ = ('id', 'original_geom_type', 'is_simplest', 'qgs_geom', 'furthest_index')

    _id_counter = 0  # Unique ID counter

    @staticmethod
    def next_id():
        """Get the next counterID.

        :param: QgsMultiLineString qgs_multi_line_string: Multi line string to merge together
        :return: ID of the SimGeom object
        :rtype: int
        """

        SimGeom._id_counter += 1

        return SimGeom._id_counter

    def __init__(self, qgs_abs_geom, original_geom_type):
        """Constructor that initialize a SimGeom object.

        :param: qgs_abs_geom: QgsAbstractGeometry to process
        :param: original_geom_type: Original type of the geometry

        """

        self.id = SimGeom.next_id()
        self.original_geom_type = original_geom_type
        qgs_geometry = qgs_abs_geom.constGet()
        self.qgs_geom = QgsGeometry(qgs_geometry.clone())
        self.is_simplest = False
        self.furthest_index = None
        # Set some variable depending on the attribute of the feature
        if self.original_geom_type == QgsWkbTypes.Point:
            self.is_simplest = True  # A point cannot be simplified
        else:
            # Original geometry is LineString or Polygon
            if qgs_geometry.length() >= Epsilon.ZERO_RELATIVE:
                if qgs_geometry.isClosed():  # Closed LineString
                    qgs_polygon = QgsPolygon(qgs_geometry.clone())  # Create QgsPolygon to calculate area
                    if qgs_polygon.area() > Epsilon.ZERO_RELATIVE:
                        if qgs_geometry.numPoints() <= 4:
                            self.is_simplest = True  # Cannot simplify a closed line with less than 4 vertices
                    else:
                        self.is_simplest = True  # Degenerated area cannot simplify
                else:
                    if qgs_geometry.numPoints() <= 2:
                        self.is_simplest = True  # Cannot simplify a line with less than 2 vertice
            else:
                self.is_simplest = True  # Degenerated area cannot simplify


class Bend:
    """Define a Bend object which is the reduction goal of this algorithm"""

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
        self.adj_area = Bend.calculate_adj_area(self.area, self.perimeter)
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


class GeoSimUtil:
    """Class containing a list general static method"""

    @staticmethod
    def validate_simplicity(qgs_geoms_with_itself, qgs_geom_new_subline):
        """Validate the simplicitity constraint

        This constraint assure that the new sub line is not intersecting with any other segment of the same line

        :param: qgs_geoms_with_itself: List of QgsLineString segment to verify for self intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        if qgs_geom_new_subline.length() > Epsilon.ZERO_RELATIVE:
            geom_engine_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet().clone())
            for qgs_geom_potential in qgs_geoms_with_itself:
                de_9im_pattern = geom_engine_subline.relate(qgs_geom_potential.constGet().clone())
                # de_9im_pattern[0] == '0' means that their interiors intersect (crosses)
                # de_9im_pattern[1] == '0' means that one extremity is touching the interior of the other (touches)
                if de_9im_pattern[0] == '0' or de_9im_pattern[1] == '0':
                    constraints_valid = False
                    break
        else:
            # Special case do not validate simplicity for almost zero length line (equivalent to a point)
            pass

        return constraints_valid

    @staticmethod
    def validate_intersection(qgs_geoms_with_others, qgs_geom_new_subline):
        """Validate the intersection constraint

        This constraint assure that the new sub line is not intersecting with any other lines (not itself)

        :param: qgs_geoms_with_others: List of QgsLineString segment to verify for intersection
        :param: qgs_geom_new_subline: New QgsLineString replacement sub line.
        :return: Flag indicating if the spatial constraint is valid
        :rtype: Bool
        """

        constraints_valid = True
        if len(qgs_geoms_with_others) >= 1:
            geom_engine_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet().clone())
            for qgs_geom_potential in qgs_geoms_with_others:
                de_9im_pattern = geom_engine_subline.relate(qgs_geom_potential.constGet().clone())
                # de_9im_pattern[0] == '0' means that their interiors intersect (crosses)
                if de_9im_pattern[0] == '0':
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
