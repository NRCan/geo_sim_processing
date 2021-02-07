"""
Unit test for reduce_bend algorithm
"""

import unittest
from qgis.core import QgsApplication
from reduce_bend_algorithm import ReduceBend
from qgis.core import QgsPoint, QgsLineString, QgsPolygon, QgsFeature, QgsGeometry, QgsProcessingFeedback

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


def build_and_launch(title, qgs_geoms, diameter_tol, del_pol=False, del_hole=False):

    print(title)
    qgs_features = []
    feedback = QgsProcessingFeedback()
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    rb_results = ReduceBend.reduce(qgs_features, diameter_tol, feedback, del_pol, del_hole, True)
    log = feedback.textLog()
    print (log)
    qgs_features_out = rb_results.qgs_features_out

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
        qgs_geom0 = create_polygon([(0,10), (5,10), (10,10), (10,0), (0,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 300)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case03(self):
        title = "Test 03: Polygon with one bend the first/end vertice located on the bend to reduce"
        qgs_geom0 = create_polygon([(5,10), (5,11), (6,11), (6,10), (10,10), (10,0), (0,0), (0,10), (5,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_polygon([(10,0), (0,0), (0,10), (10,10), (10,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        print (qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case04(self):
        title = "Test 04: Square polygon with one bend"
        qgs_geom0 = create_polygon([(0,10), (5,9), (10,10), (10,0), (0,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,0), (0,0), (0,10), (10,10), (10,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        print (qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case05(self):
        title = "Test 05: triangle polygon with one bend"
        qgs_geom0 = create_polygon([(0,10), (5,9), (10,10), (5,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3000)
        out_qgs_geom0 = create_polygon([(10,10), (5,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case06(self):
        title = "Test 06: A polygon with no bend.  A line with no bend"
        qgs_geom0 = create_polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)], [])
        qgs_geom1 = create_line([(10, 0), (20, 0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case07(self):
        title = "Test 07: 1 polygon with no bend to reduce"
        coords0 = [(0, 0), (0, 5), (2.5,4), (5, 5), (5, 0), (0,0)]
        qgs_geom0 = create_polygon(coords0, [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        qgs_geom0_out = create_polygon([(5, 5), (5, 0), (0, 0), (0, 5), (2.5, 4), (5, 5)], [])
        val0 = qgs_geom0_out.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


    def test_case08(self):
        title = "Test 08: 1 point and 3 line string to validate bounding box sub dividing"
        qgs_geom0 = create_line([(0, 0), (1, 1), (2,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_line([(0, 0), (2, 0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case09(self):
        title = "Test 09: 1 point and 3 line string to validate bounding box sub dividing"
        qgs_geom0 = create_point((0,0))
        qgs_geom1 = create_line([(0, 0), (100, 0)])
        qgs_geom2 = create_line([(0, 0), (0, 100)])
        qgs_geom3 = create_line([(0, 0), (100, 100)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1, qgs_geom2, qgs_geom3], 30)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        val2 = qgs_geom2.equals(qgs_feature_out[2])
        val3 = qgs_geom3.equals(qgs_feature_out[3])
        self.assertTrue (val0 and val1 and val2, title)

    def test_case10(self):
        title = "Test 10: 2 simple line segment, simple triangle  and one point"
        qgs_geom0 = create_line([(0,0),(30,0)])
        qgs_geom1 = create_polygon([(10,10),(15,20), (20,10), (10,10)], [])
        qgs_geom2 = create_point((0,100))
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1, qgs_geom2], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        val2 = qgs_geom2.equals(qgs_feature_out[2])
        self.assertTrue (val0 and val1 and val2, title)

    def test_case11(self):
        title = "Test 11: Zero length line"
        qgs_geom0 = create_line([(10, 10), (10, 10)])
        qgs_geom1 = create_line([(20, 20), (20, 20), (20,20)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case12(self):
        title = "Test 12: Degenerated line"
        qgs_geom0 = create_line([(10, 10), (10, 20), (10,10)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3)
        out_qgs_geom0 = create_line([(10, 10), (10,20), (10,10)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case13(self):
        title = "Test 13: Line with segment parrallel to itself"
        qgs_geom0 = create_line([(0,0),(30,0), (20,0)])
        qgs_geom1 = create_line([(0, 10), (-5,10), (30, 10)])
        qgs_geom2 = create_line([(0, 20), (-5, 20), (30,20), (20, 20)])
        out_qgs_geom0 = create_line([(0,0), (20,0)])
        out_qgs_geom1 = create_line([(0, 10), (30, 10)])
        out_qgs_geom2 = create_line([(0, 20), (20, 20)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1, qgs_geom2], 3)
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        val2 = out_qgs_geom2.equals(qgs_feature_out[2])
        self.assertTrue (val0 and val1 and val2, title)

    def test_case14(self):

        title = "Test 14: Co-linear and alomost co-linear point"
        in_geom0 = create_line([(0, 0), (20, 0), (25.000000000000001, 0.0000000000001), (30, 0)])
        in_geom1 = create_line([(0, 10), (30, 10), (35.000000000001, 10.00000000000001), (40, 10)])
        in_geom2 = create_point((0, 100))
        out_geom0 = create_line([(0, 0), (30, 0)])
        out_geom1 = create_line([(0, 10), (40, 10)])
        qgs_feature_out = build_and_launch(title, [in_geom0, in_geom1, in_geom2], 3)
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        val2 = in_geom2.equals(qgs_feature_out[2])
        self.assertTrue(val0 and val1 and val2, title)

    def test_case15(self):
        title = "Test 15: Small bend"
        in_geom0 = create_line([(0, 0), (30, 0)])
        in_geom1 = create_line([(0, 10), (30, 10), (30, 11), (31, 11), (31, 10), (40, 10), (50, 10), (50, 11), (51, 10), (60, 10)])
        in_geom2 = create_point((0, 100))
        out_geom0 = create_line([(0, 0), (30, 0)])
        out_geom1 = create_line([(0, 10), (60, 10)])
        qgs_feature_out = build_and_launch(title, [in_geom0, in_geom1, in_geom2], 3)
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        val2 = in_geom2.equals(qgs_feature_out[2])
        self.assertTrue(val0 and val1 and val2, title)

    def test_case16(self):
        title = "Test 16: Polygon with bend"
        outer = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0), (0, 0)]
        inner = [(5, 5), (5, 6), (6, 6), (6, 5)]
        in_geom0 = create_polygon(outer, [inner])
        qgs_feature_out = build_and_launch(title, [in_geom0], 300)
        outer = [(20, 20), (20, 0), (0, 0), (0,20), (20,20)]
        inner = [(5, 5), (5, 6), (6, 6), (6, 5), (5, 5)]
        out_geom0 = create_polygon(outer, [inner])
        val0 = out_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case17(self):
        title = "Test 17: Polygon with line in bend"
        coord = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0)]
        qgs_geom0 = create_polygon(coord, [])
        qgs_geom1 = create_line([(10.1, 20.5), (10.2, 20.6), (10.3, 20.5)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        coord = [(20, 20), (20, 0), (0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        out_geom1 = create_line([(10.1, 20.5), (10.3, 20.5)])
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        self.assertTrue(val0 and val1, title)

    def test_case18(self):
        title = "Test 18: Polygon with point in bend"
        coord = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0), (0,0)]
        qgs_geom0 = create_polygon(coord, [])
        qgs_geom1 = create_point((10.1,20.5))
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        coord = [(20, 20), (20, 0), (0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue(val0 and val1, title)

    def test_case19(self):
        title = "Test 19: Line String self intersecting after bend reduction"
        coord = [(0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (30, 20), (30,0), (10.5,0), (10.5,20.5)]
        qgs_geom0 = create_line(coord)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case20(self):
        title = "Test 20: Polygon with a bend containing a hole to delete"
        coord0 = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0)]
        coord1 = [(10.1, 20.1), (10.1, 20.2), (10.2, 20.2), (10.2, 20.1), (10.1, 20.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        coord = [(20,20), (20,0), (0,0), (0,20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        val0 = out_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case21(self):
        title = "Test 21: Small polygon no bend with one small hole. Hole is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=False, del_hole=True)
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        qgs_geom0 = create_polygon(coord0, [])
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case22(self):
        title = "Test 22: Small polygon no bend with one small hole. Feature is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=False)
        self.assertEqual(len(qgs_feature_out), 0, title)

    def test_case23(self):
        title = "Test 23: Small polygon no bend with one small hole. Feature is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        self.assertEqual(len(qgs_feature_out), 0, title)

    def test_case24(self):
        title = "Test 24: A line with a bend where the length of the base is zero"
        coord0 = [(0, 0), (50, 0), (49, 1), (51, 1), (50, 0), (100, 0)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        coord0 = [(0, 0), (100, 0)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)


# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", True)

# profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
#profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False)

# Load providers
app.initQgis()
