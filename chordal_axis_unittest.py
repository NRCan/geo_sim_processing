"""
Unit test for ChordalAxis algorithm
"""

import unittest
from chordal_axis_algorithm import ChordalAxis, GenUtil
from qgis.core import QgsPoint, QgsLineString, QgsPolygon, QgsFeature, QgsGeometry, QgsProcessingFeedback, \
                      QgsVectorLayer
from qgis.analysis import QgsNativeAlgorithms
import processing
from qgis.core import QgsApplication
from qgis._3d import Qgs3DAlgorithms



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


def build_and_launch(title, qgs_geom_pol, correction=False):

    print(title)
    vl = QgsVectorLayer("Polygon", "temporary_polygon", "memory")
    pr = vl.dataProvider()
    fet = QgsFeature()
    fet.setId(1)
    fet.setGeometry(qgs_geom_pol)
    pr.addFeatures([fet])

    # update layer's extent when new features have been added
    # because change of extent in provider is not propagated to the layer
    vl.updateExtents()
    feedback = QgsProcessingFeedback()
    qgs_multi_triangles = ChordalAxis.tessellate_polygon(vl, feedback)

    ca = ChordalAxis(qgs_multi_triangles[0], GenUtil.ZERO)
    if correction:
        ca.correct_skeleton()
    centre_lines = ca.get_skeleton()

    return centre_lines

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

def coords_shift(pos, coords):

    new_coords = coords[pos:] + coords[0:pos-1] + [coords[pos]]

    return new_coords

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
        title = "Test 01 - Triangle - No skeleton produced"
        qgs_geom_pol = create_polygon([(0,0), (10,10), (20,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        self.assertTrue(len(centre_lines)==0, title)

    def test_case02(self):
        title = "Test 02 Small square - 2 triangles terminal without correction"
        qgs_geom_pol = create_polygon([(0,0), (10,0), (10,10), (0,10), (0,0)], [])

        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(0, 0), (5,5), (10,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        self.assertTrue(val0, title)

    def test_case03(self):
        title = "Test 03 Small square - 2 triangles terminal with correction"
        qgs_geom_pol = create_polygon([(0,0), (10,0), (10,10), (0,10), (0,0)], [])

        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(0, 0), (5,5), (10,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        self.assertTrue(val0, title)

    def test_case04(self):
        title = "Test 04 Long rectangle - 2 triangles terminal and 2 sleeves without correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (10,10), (20,10), (20,0), (10,0), (0,0)], [])

        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(0,0), (5,5), (10,5), (15,5), (20,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        self.assertTrue(val0, title)

    def test_case05(self):
        title = "Test 05 Long rectangle - 2 triangles terminal and 2 sleeves with correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (10,10), (20,10), (20,0), (10,0), (0,0)], [])

        centre_lines = build_and_launch(title, qgs_geom_pol, correction=True)
        qgs_geom0 = create_line([(0,0), (5,5), (10,5), (15,5), (20,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        self.assertTrue(val0, title)

    def test_case06(self):
        title = "Test 06 Long rectangle - with junction terminal without correction"
        qgs_geom_pol = create_polygon( [(0, 0), (0, 10), (9, 10), (10, 11), (11, 10), (20, 10), (20, 0), (10, 0), (0, 0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(10,6.66666666666666696), (9.5,5), (4.5,5), (0,10)])
        qgs_geom1 = create_line([(10,6.66666666666666696), (10,10), (10,11)])
        qgs_geom2 = create_line([(10,6.66666666666666696), (10.5,5), (15.5,5), (20,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        val1 = qgs_geom1.equals(QgsGeometry(centre_lines[1].clone()))
        val2 = qgs_geom2.equals(QgsGeometry(centre_lines[2].clone()))
        self.assertTrue(val0 and val1 and val2, title)

    def test_case07(self):
        title = "Test 07 Long rectangle - with junction terminal with correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (9,10), (10,11), (11,10), (20,10), (20,0), (10,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=True)
        qgs_geom0 = create_line([(0,10), (4.5,5), (9.5,5), (10.5,5), (15.5,5), (20,10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        self.assertTrue(val0, title)

    def test_case08(self):
        title = "Test 08 Narrow T Junction  without correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (25,10), (50,10), (50,0), (30,0), (30,-30), (20,-30), (20,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(0, 0), (10, 5), (22.5, 5), (25, 3.33333333333333348)])
        qgs_geom1 = create_line([(20, -30), (25, -15), (25, 0), (25, 3.33333333333333348)])
        qgs_geom2 = create_line([(25, 3.33333333333333348), (27.5, 5), (40, 5), (50, 0)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        val1 = qgs_geom1.equals(QgsGeometry(centre_lines[1].clone()))
        val2 = qgs_geom2.equals(QgsGeometry(centre_lines[2].clone()))
        self.assertTrue(val0 and val1 and val2, title)

    def test_case09(self):
        title = "Test 09 Narrow T Junction  with correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (25,10), (50,10), (50,0), (30,0), (30,-30), (20,-30), (20,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=True)
        qgs_geom0 = create_line([(0, 0), (10, 5), (22.5, 5), (25, 5)])
        qgs_geom1 = create_line([(20, -30), (25, -15), (25, 0), (25, 5)])
        qgs_geom2 = create_line([(25, 5), (27.5, 5), (40, 5), (50, 0)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        val1 = qgs_geom1.equals(QgsGeometry(centre_lines[1].clone()))
        val2 = qgs_geom2.equals(QgsGeometry(centre_lines[2].clone()))
        self.assertTrue(val0 and val1 and val2, title)

    def test_case10(self):
        title = "Test 10 Narrow X Junction  without correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (20,10), (20,40),(30,40),(30,10), (50,10), (50,0), (30,0), (30,-30), (20,-30), (20,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=False)
        qgs_geom0 = create_line([(0, 0), (10, 5), (20, 5), (23.33333333333333215, 3.33333333333333348)])
        qgs_geom1 = create_line([(20, -30), (25, -15), (25, 0), (23.33333333333333215, 3.33333333333333348)])
        qgs_geom2 = create_line([(23.33333333333333215, 3.33333333333333348), (25, 5), (26.66666666666666785, 6.66666666666666696)])
        qgs_geom3 = create_line([(26.66666666666666785, 6.66666666666666696), (25, 10), (25, 25), (30, 40)])
        qgs_geom4 = create_line([(26.66666666666666785, 6.66666666666666696), (30, 5), (40, 5), (50, 10)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        val1 = qgs_geom1.equals(QgsGeometry(centre_lines[1].clone()))
        val2 = qgs_geom2.equals(QgsGeometry(centre_lines[2].clone()))
        val3 = qgs_geom3.equals(QgsGeometry(centre_lines[3].clone()))
        val4 = qgs_geom4.equals(QgsGeometry(centre_lines[4].clone()))
        self.assertTrue(val0 and val1 and val2 and val3 and val4, title)

    def test_case11(self):
        title = "Test 11 Narrow X Junction  with correction"
        qgs_geom_pol = create_polygon([(0,0), (0,10), (20,10), (20,40),(30,40),(30,10), (50,10), (50,0), (30,0), (30,-30), (20,-30), (20,0), (0,0)], [])
        centre_lines = build_and_launch(title, qgs_geom_pol, correction=True)
        qgs_geom0 = create_line([(0, 0), (10, 5), (20, 5), (25, 5)])
        qgs_geom1 = create_line([(20, -30), (25, -15), (25, 0), (25, 5)])
        qgs_geom2 = create_line([(25, 5), (30, 5), (40, 5), (50, 10)])
        qgs_geom3 = create_line([(25, 5), (25, 10), (25, 25), (30, 40)])
        val0 = qgs_geom0.equals(QgsGeometry(centre_lines[0].clone()))
        val1 = qgs_geom1.equals(QgsGeometry(centre_lines[1].clone()))
        val2 = qgs_geom2.equals(QgsGeometry(centre_lines[2].clone()))
        val3 = qgs_geom3.equals(QgsGeometry(centre_lines[3].clone()))
        self.assertTrue(val0 and val1 and val2 and val3, title)


QgsApplication.processingRegistry().addProvider(Qgs3DAlgorithms(QgsApplication.processingRegistry()))
QgsApplication.setPrefixPath("/usr/bin/qgis", False)
profile_folder = '.'

# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers and init QGIS
app.initQgis()
QgsApplication.processingRegistry().addProvider(Qgs3DAlgorithms(QgsApplication.processingRegistry()))
from processing.core.Processing import Processing
Processing.initialize()
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())