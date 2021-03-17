"""
Unit test for reduce_bend algorithm
"""

import unittest
from qgis.core import QgsApplication
from building_pattern_algorithm import BuildingPattern
from qgis.core import QgsPoint, QgsLineString, QgsPolygon, QgsFeature, QgsGeometry, QgsProcessingFeedback, \
                      QgsVectorLayer, QgsWkbTypes, QgsPointXY
from qgis.analysis import QgsNativeAlgorithms

def qgs_line_string_to_xy(qgs_line_string):

    qgs_points = qgs_line_string.points()
    lst_x = []
    lst_y = []
    for qgs_point in qgs_points:
        lst_x.append(qgs_point.x())
        lst_y.append(qgs_point.y())

    return (lst_x, lst_y)

def plot_lines(qgs_line_string, qgs_new_line):

    line0_lst_x, line0_lst_y = qgs_line_string_to_xy(qgs_line_string)
#    line1_lst_x, line1_lst_y = qgs_line_string_to_xy(qgs_new_line)

    import matplotlib.pyplot as plt
    plt.plot(line0_lst_x, line0_lst_y, 'b')
#    plt.plot(line1_lst_x, line1_lst_y, 'r')
    plt.show()


def build_and_launch(title, qgs_geoms, rectangularity_tol, compactness_tol):

    print(title)
    qgs_features = []
    feedback = QgsProcessingFeedback()
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    bp_results = BuildingPattern.match(qgs_features, rectangularity_tol, compactness_tol, feedback)
    log = feedback.textLog()
    print (log)
    qgs_features_out = bp_results.qgs_features_out

    qgs_geoms_out = []
    for qgs_feature_out in qgs_features_out:
        qgs_geoms_out.append(qgs_feature_out.geometry())

    return qgs_geoms_out

def create_line(coords, ret_geom=True):

    qgs_points = []
    for coord in coords:
        qgs_points.append(create_point(coord, False))

    if ret_geom:
        ret_val = QgsGeometry(QgsLineString(qgs_points))
    else:
        ret_val = QgsLineString(qgs_points).clone()

    return ret_val

def create_point(coord, ret_geom=True):

    qgs_point = QgsPoint(coord[0], coord[1])
    if ret_geom:
        ret_val = QgsGeometry(qgs_point)
    else:
        ret_val = qgs_point.clone()

    return ret_val

def create_polygon(outer, inners):

    outer_line = create_line(outer, False)
    qgs_pol = QgsPolygon()
    qgs_pol.setExteriorRing(outer_line)
    for inner in inners:
        inner_line = create_line(inner, False)
        qgs_pol.addInteriorRing(inner_line)
    qgs_geom = QgsGeometry(qgs_pol)

    return qgs_geom


class Test(unittest.TestCase):
    """
    Class allowing to test the algorithm
    """

    def test_case01(self):
        title = "Test 01: Empty file"


        qgs_feature_out = build_and_launch(title, [], 5, True, True)
        out_qgs_geom0 = create_polygon([(0, 10), (10, 10), (10, 0), (0, 0), (0, 10)], [])
        if len(qgs_feature_out) == 0:
            val0 = True
        else:
            val0 = False
        self.assertTrue(val0, title)

    def test_case02(self):
        title = "Test 02: Polygon with start/end point colinear"
        qgs_geom0 = create_polygon([(0,0), (0,10), (10,10), (8,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], .85, .9)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case03(self):
        title = "Test 03: Pattern in L form"
        qgs_geom0 = create_polygon([(0,0), (0,10), (10,10), (10,8), (8,8), (8,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], .85, .9)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", True)

# profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
#profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], True)

# Load providers and init QGIS
app.initQgis()
from processing.core.Processing import Processing
Processing.initialize()
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
