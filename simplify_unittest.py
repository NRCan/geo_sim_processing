# -*- coding: utf-8 -*-

# /***************************************************************************
# simplify_unittest.py
# ----------
# Date                 : april 2021
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
Unit test for simplify algorithm
"""

import unittest
from qgis.core import QgsApplication
from .simplify_algorithm import Simplify
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


def build_and_launch(title, qgs_geoms, tolerance):

    print(title)
    qgs_features = []
    feedback = QgsProcessingFeedback()
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    rb_results = Simplify.douglas_peucker(qgs_features, tolerance, True, feedback)
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
        qgs_feature_out = build_and_launch(title, [], 2)
        if len(qgs_feature_out) == 0:
            val0 = True
        else:
            val0 = False
        self.assertTrue(val0, title)

    def test_case02(self):
        title = "Test 02: Open line with 2 vertice"
        qgs_geom0 = create_line([(0, 0), (10,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 5)
        out_qgs_geom0 = create_line([(0, 0), (10,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case03(self):
        title = "Test 03: Open line with 3 vertice"
        qgs_geom0 = create_line([(0, 0), (5, 1), (10,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 5)
        out_qgs_geom0 = create_line([(0, 0), (10,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case04(self):
        title = "Test 04: Open line with 4 vertice, farthest point is index: 2"
        qgs_geom0 = create_line([(0, 0), (5, 3), (10,4), (20,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 5)
        out_qgs_geom0 = create_line([(0, 0), (20,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case05(self):
        title = "Test 05: Open line with 4 vertice, farthest point is index: 1"
        qgs_geom0 = create_line([(0, 0), (5, 4), (10,3), (20,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 5)
        out_qgs_geom0 = create_line([(0, 0), (20,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case06(self):
        title = "Test 06: Open line with 5 vertice, farthest point is in the middle"
        qgs_geom0 = create_line([(0, 0), (5, 3), (10,4), (15,3), (20,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 5)
        out_qgs_geom0 = create_line([(0, 0), (20,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case07(self):
        title = "Test 07: Open line with 5 vertice, only vertice 1 and 4 are simplified"
        qgs_geom0 = create_line([(0, 0), (5, 3), (10,4), (15,3), (20,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3.5)
        out_qgs_geom0 = create_line([(0, 0), (10,4), (20,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case08(self):
        title = "Test 08: Open line with 5 vertice, no simplification; below tolrance"
        qgs_geom0 = create_line([(0, 0), (5, 3), (10,4), (15,3), (20,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], .1)
        out_qgs_geom0 = create_line([(0, 0), (5, 3), (10,4), (15,3), (20,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case09(self):
        title = "Test 09: Closed line in form of triangle, no simplification"
        qgs_geom0 = create_polygon([(0, 0), (5, 5), (10,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 10)
        out_qgs_geom0 = create_polygon([(0, 0), (5, 5), (10,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10(self):
        title = "Test 10: Closed line in form of a square"
        qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 10)
        out_qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case11(self):
        title = "Test 11: Closed line in form of a square"
        qgs_geom0 = create_polygon([(0, 0), (0,5), (3,6), (5,5), (5,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 10)
        out_qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case12(self):
        title = "Test 12: Closed line in form of a square"
        qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (2,1), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 10)
        out_qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case13(self):
        title = "Test 13: Closed line in form of a square"
        qgs_geom0 = create_polygon([(0, 0), (0, 5), (3, 5), (5, 5), (5, 0), (0, 0)], [])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 2)
        out_qgs_geom0 = create_polygon([(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case14(self):
        title = "Test 14: Closed line in form of a square"
        qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (2,1), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 2)
        out_qgs_geom0 = create_polygon([(0, 0), (0,5), (5,5), (5,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case15(self):
        title = "Test 15: Open line self intersecting"
        qgs_geom0 = create_line([(0, 0), (5,0), (5,2), (10,2), (10,0), (50,0), (50, -5), (7,-5), (7,1)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_line([(0, 0), (5,2), (50,0), (50, -5), (7,-5), (7,1)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case16(self):
        title = "Test 16: Open line intersecting another line (no simplification done)"
        qgs_geom0 = create_line([(0, 0), (2,2), (4,0)])
        qgs_geom1 = create_line([(2, -1), (2, 1)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (2,2), (4,0)])
        out_qgs_geom1 = create_line([(2, -1), (2, 1)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case17(self):
        title = "Test 17: Open line intersecting another line: simplification done"
        qgs_geom0 = create_line([(0, 1), (3,3), (6,1)])
        qgs_geom1 = create_line([(0, 0), (3,1.5), (6,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 1), (6, 1)])
        out_qgs_geom1 = create_line([(0, 0), (6, 0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case18(self):
        title = "Test 18: Open line intersecting another line: simplification done"
        qgs_geom0 = create_line([(0, 0), (3, 1.5), (6, 0)])
        qgs_geom1 = create_line([(0, 1), (3,3), (6,1)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (6, 0)])
        out_qgs_geom1 = create_line([(0, 1), (6, 1)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case19(self):
        title = "Test 19: Open line intersecting another line (partial line simplification done)"
        qgs_geom0 = create_line([(0, 0), (2,2), (4,0), (6,2), (8,0)])
        qgs_geom1 = create_line([(2, -1), (2, .5)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (2,2), (8,0)])
        out_qgs_geom1 = create_line([(2, -1), (2, .5)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case20(self):
        title = "Test 20: Open line intersecting another line (partial line simplification done)"
        qgs_geom0 = create_line([(0, 0), (2,2), (4,0), (6,2.5), (8,0)])
        qgs_geom1 = create_line([(6, -1), (6, .5)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (6,2.5), (8,0)])
        out_qgs_geom1 = create_line([(6, -1), (6, .5)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case21(self):
        title = "Test 21: Open line with sidedness problem"
        qgs_geom0 = create_line([(0, 0), (2,2), (4,0)])
        qgs_geom1 = create_line([(2,.1), (2,.2)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (2,2), (4,0)])
        out_qgs_geom1 = create_line([(2,.1), (2,.2)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case22(self):
        title = "Test 22: Open line with sidedness problem"
        qgs_geom0 = create_line([(0, 0), (2,2), (4,0), (6,2.5), (8,0)])
        qgs_geom1 = create_line([(2,.1), (2,.2)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (6,2.5), (8,0)])
        out_qgs_geom1 = create_line([(2,.1), (2,.2)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case23(self):
        title = "Test 23: Two disjoint open line string with extremity touching ==> simplified"
        qgs_geom0 = create_line([(0,0), (2,2), (4,0)])
        qgs_geom1 = create_line([(6,0), (8,2), (10,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        out_qgs_geom1 = create_line([(6, 0), (10,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case24(self):
        title = "Test 24: Two open line string with extremity touching ==> simplified"
        qgs_geom0 = create_line([(0,0), (2,2), (4,0)])
        qgs_geom1 = create_line([(4,0), (6,2), (8,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        out_qgs_geom1 = create_line([(4, 0), (8,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case25(self):
        title = "Test 25: Two open line string with one extremity touching the middle of the other line: simplified"
        qgs_geom0 = create_line([(0,0), (2,2), (4,0)])
        qgs_geom1 = create_line([(-2,0), (2,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        out_qgs_geom1 = create_line([(-2, 0), (2,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case26(self):
        title = "Test 26: Two open line with one extremity superimposed in the middle of the other line: simplified"
        qgs_geom0 = create_line([(0,0), (2,2), (4,0)])
        qgs_geom1 = create_line([(1,0), (3,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        out_qgs_geom1 = create_line([(1, 0), (3,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case27(self):
        title = "Test 27: One open line with with duplicate points: simplified"
        qgs_geom0 = create_line([(0,0), (2,2), (2,2), (4,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case28(self):
        title = "Test 28: One open line with with duplicate points: simplified"
        qgs_geom0 = create_line([(0,0), (0,0), (2,2), (2,2), (2,2), (4,0), (4,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_line([(0, 0), (4,0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case29(self):
        title = "Test 29: Different degenerated line string (points identical)"
        qgs_geom0 = create_line([(0,0), (0,0)])
        qgs_geom1 = create_line([(10, 10), (10, 10), (10,10)])
        qgs_geom2 = create_line([(20, 20), (20, 20), (20, 20), (20,20)])
        qgs_geom3 = create_line([(30, 30), (30, 30), (30, 30), (30, 30), (30, 30)])
        qgs_geom4 = create_line([(40, 40), (40, 40), (40, 40), (40, 40), (40, 40), (40,40)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1, qgs_geom2, qgs_geom3, qgs_geom4,],15)
        out_qgs_geom0 = create_line([(0, 0), (0, 0)])
        out_qgs_geom1 = create_line([(10, 10), (10, 10), (10, 10)])
        out_qgs_geom2 = create_line([(20, 20), (20, 20), (20, 20), (20, 20)])
        out_qgs_geom3 = create_line([(30, 30), (30, 30), (30, 30), (30, 30), (30, 30)])
        out_qgs_geom4 = create_line([(40, 40), (40, 40), (40, 40), (40, 40), (40, 40), (40, 40)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        val1 = out_qgs_geom1.equals(qgs_feature_out[1])
        val2 = out_qgs_geom2.equals(qgs_feature_out[2])
        val3 = out_qgs_geom3.equals(qgs_feature_out[3])
        val4 = out_qgs_geom4.equals(qgs_feature_out[4])
        self.assertTrue (val0 and val1 and val2 and val3 and val4, title)




# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", True)

# profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
#profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False)

# Load providers and init QGIS
app.initQgis()
from processing.core.Processing import Processing
Processing.initialize()
QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
