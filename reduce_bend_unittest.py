# -*- coding: utf-8 -*-

# /***************************************************************************
# reduce_bend_unittest.py
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
Unit test for reduce_bend algorithm
"""

import unittest
from qgis.core import QgsApplication
from .reduce_bend_algorithm import ReduceBend
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


def build_and_launch(title, qgs_geoms, diameter_tol, del_pol=False, del_hole=False, smooth_line=False):

    print(title)
    qgs_features = []
    feedback = QgsProcessingFeedback()
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    rb_results = ReduceBend.reduce(qgs_features, diameter_tol, smooth_line, del_pol, del_hole, True, feedback)
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
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case04(self):
        title = "Test 04: Square polygon with one bend"
        qgs_geom0 = create_polygon([(0,10), (5,9), (10,10), (10,0), (0,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case05(self):
        title = "Test 05: triangle polygon with one bend"
        qgs_geom0 = create_polygon([(0,10), (5,9), (10,10), (5,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3000)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (5,0), (0,10)], [])
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
        qgs_geom0_out = create_polygon(coords0, [])
        val0 = qgs_geom0_out.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


    def test_case08(self):
        title = "Test 08: 1 line string with one bend to simplify"
        qgs_geom0 = create_line([(0, 0), (1, 1), (2,0)])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 3)
        out_qgs_geom0 = create_line([(0, 0), (2, 0)])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case09(self):
        title = "Test 09: 1 point and 3 line string no simplification"
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

    def test_case10_1_1(self):
        title = "Test 10_1_1: Triangle with one bend (1/4)"
        qgs_geom0 = create_polygon([(10,10),(15,20), (20,10), (15,11), (10,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10),(15,20), (20,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_1_2(self):
        title = "Test 10_1_2: Triangle with one bend (2/4 pivot first vertice)"
        qgs_geom0 = create_polygon([(15,20), (20,10), (15,11), (10,10), (15,20)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(15,20), (20,10), (10,10), (15,20)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_1_3(self):
        title = "Test 10_1_3: Triangle with one bend (3/4 pivot first vertice)"
        qgs_geom0 = create_polygon([(20,10), (15,11), (10,10), (15,20), (20,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(20,10), (10,10), (15,20), (20,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_1_4(self):
        title = "Test 10_1_4: Triangle with one bend (4/4 pivot first vertice)"
        qgs_geom0 = create_polygon([(15,11), (10,10), (15,20), (20,10), (15,11)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (15,20), (20,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_2_1(self):
        title = "Test 10_2_1: Square with 2 bends (1/4)"
        qgs_geom0 = create_polygon([(0,0),(1,1),(0,2),(10,2), (9,1),(10,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,0),(0,2),(10,2),(10,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_2_2(self):
        title = "Test 10_2_2: Square with 2 bends (2/4) pivot first vertice"
        qgs_geom0 = create_polygon([(1,1),(0,2),(10,2), (9,1),(10,0), (0,0), (1,1)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,2),(10,2),(10,0), (0,0), (0,2)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_2_3(self):
        title = "Test 10_7: Square with 2 bends (3/4) pivot first vertice"
        qgs_geom0 = create_polygon([(0,2),(10,2), (9,1),(10,0), (0,0), (1,1), (0,2)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,2),(10,2),(10,0), (0,0), (0,2)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_2_4(self):
        title = "Test 10_8: Square with 2 bends (4/4) pivot first vertice"
        qgs_geom0 = create_polygon([(10,2), (9,1),(10,0), (0,0), (1,1), (0,2), (10,2)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,2),(10,0), (0,0), (0,2), (10,2)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_1(self):
        title = "Test 10_3_1: Square with 2 bends each bend with 2 vertices (1/7)"
        qgs_geom0 = create_polygon([(0,0), (0,10), (4,10), (4,9), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,0), (0,10), (10,10), (10,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_2(self):
        title = "Test 10_3_2: Square with 2 bends each bend with 2 vertices (2/7) pivot first vertice"
        qgs_geom0 = create_polygon([(0,10), (4,10), (4,9), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_3(self):
        title = "Test 10_9_2: Square with 2 bends each bend with 2 vertices (3/7) pivot first vertice"
        qgs_geom0 = create_polygon([(4,10), (4,9), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10), (4,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_4(self):
        title = "Test 10_3_4: Square with 2 bends each bend with 2 vertices (4/7) pivot first vertice"
        qgs_geom0 = create_polygon([(4,9), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10), (4,10), (4,9)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_5(self):
        title = "Test 10_3_5: Square with 2 bends each bend with 2 vertices (5/7) pivot first vertice"
        qgs_geom0 = create_polygon([(6,9), (6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (6,9)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_6(self):
        title = "Test 10_3_6: Square with 2 bends each bend with 2 vertices (6/7) pivot first vertice"
        qgs_geom0 = create_polygon([(6,10), (10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (6,9), (6,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_3_7(self):
        title = "Test 10_3_7: Square with 2 bends each bend with 2 vertices (7/7) pivot first vertice"
        qgs_geom0 = create_polygon([(10,10), (10,0), (6,0), (6,1),(4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (6,9), (6,10), (10,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


    def test_case10_4_1(self):
        title = "Test 10_4_1: Square with 2 bends each bend with 3 vertices (1/7)"
        qgs_geom0 = create_polygon([(0,0), (0,10), (4,10), (4,9), (5,9.5), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,0), (0,10), (10,10), (10,0), (0,0)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_4_2(self):
        title = "Test 10_4_2: Square with 2 bends each bend with 3 vertices (2/7)"
        qgs_geom0 = create_polygon([(0,10), (4,10), (4,9), (5,9.5), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(0,10), (10,10), (10,0), (0,0), (0,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


    def test_case10_4_3(self):
        title = "Test 10_4_3: Square with 2 bends each bend with 3 vertices (3/7)"
        qgs_geom0 = create_polygon([(4,10), (4,9), (5,9.5), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10), (4,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_4_4(self):
        title = "Test 10_4_4: Square with 2 bends each bend with 3 vertices (4/7)"
        qgs_geom0 = create_polygon([(4,9), (5,9.5), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10), (4,10), (4,9)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_4_5(self):
        title = "Test 10_4_5: Square with 2 bends each bend with 3 vertices (5/7)"
        qgs_geom0 = create_polygon([(5,9.5), (6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (5,9.5)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_4_6(self):
        title = "Test 10_4_6: Square with 2 bends each bend with 3 vertices (6/7)"
        qgs_geom0 = create_polygon([(6,9), (6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (5,9.5), (6,9)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)

    def test_case10_4_7(self):
        title = "Test 10_4_7: Square with 2 bends each bend with 3 vertices (7/7)"
        qgs_geom0 = create_polygon([(6,10), (10,10), (10,0), (6,0), (6,1), (5,1.5), (4,1), (4,0), (0,0), (0,10), (4,10), (4,9), (5,9.5), (6,9), (6,10)], [])
        qgs_feature_out = build_and_launch(title,[qgs_geom0], 30)
        out_qgs_geom0 = create_polygon([(10,10), (10,0), (0,0), (0,10), (10,10)], [])
        val0 = out_qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue (val0, title)


    def test_case11(self):
        title = "Test 11: Zero length line"
        qgs_geom0 = create_line([(10, 10), (10, 10)])
        qgs_geom1 = create_line([(20, 20), (20, 20), (20,20)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue (val0 and val1, title)

    def test_case12(self):
        title = "Test 12: Degenerated polygon"
        qgs_geom0 = create_line([(10, 10), (10, 20), (10,10)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3)
        out_qgs_geom0 = create_line([(10, 10), (10, 20), (10,10)])
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
        outer = [(0, 0), (0,20), (20,20), (20,0), (0,0)]
        inner = [(5, 5), (5, 6), (6, 6), (6, 5), (5, 5)]
        out_geom0 = create_polygon(outer, [inner])
        val0 = out_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case17(self):
        title = "Test 17: Polygon with line in bend"
        coord = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0), (0,0)]
        qgs_geom0 = create_polygon(coord, [])
        qgs_geom1 = create_line([(10.1, 20.5), (10.2, 20.6), (10.3, 20.5)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
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
        coord = [(0,0), (0,20), (20,20), (20,0), (0,0)]
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
        title = "Test 24: A line with a bend where the length of the base is zero (non simple line)"
        coord0 = [(0, 0), (50, 0), (49, 1), (51, 1), (50, 0), (100, 0)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        coord0 = [(0, 0), (100, 0)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case25(self):
        title = "Test 25: A line with a bend with the form of a wave"
        coord0 = [(0, 0), (50, 0), (50,2), (49,2), (49,1), (48,1), (48,3), (51,3), (51,0), (100,0)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 10, del_pol=True, del_hole=True)
        coord0 = [(0, 0), (100, 0)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case26(self):
        title = "Test 26: A line with a smooth line replacing the bend"
        tmp_coord = [(0, -25), (25,0), (25,1), (29,1), (29,0), (50,-25)]
        tmp_coord_out = [(0, -25), (25, 0), (26.33333333333333215, 0.76980035891950094),
                         (27.66666666666666785, 0.76980035891950094), (29, 0), (50, -25)]
        coord0 = list(tmp_coord)
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3.9, del_pol=True, del_hole=True, smooth_line=True)
        coord0_out = list(tmp_coord_out)
        qgs_geom0 = create_line(coord0_out)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)
        #
        coord0 = list(tmp_coord)
        coord0_out = list(tmp_coord_out)
        coord0.reverse()
        coord0_out.reverse()
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3.9, del_pol=True, del_hole=True, smooth_line=True)
        qgs_geom0 = create_line(coord0_out)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

        for angle in [45., 90, 135, 180, 225, 270, 300]:
            coord0 = list(tmp_coord)
            qgs_geom0 = create_line(coord0)
            qgs_geom0.rotate(angle, QgsPointXY(0,0))
            qgs_geom0.translate(25,25)
            qgs_feature_out = build_and_launch(title, [qgs_geom0], 3.9, del_pol=True, del_hole=True, smooth_line=True)
            qgs_geom_out = qgs_feature_out[0]
            qgs_geom_out.translate(-25, -25)
            qgs_geom_out.rotate(-angle, QgsPointXY(0, 0))
            grid = .0000000001
            coord_ref = list(tmp_coord_out)
            qgs_geom0 = create_line(coord_ref)
            qgs_geom_grid_out = qgs_geom_out.snappedToGrid (grid,grid)
            qgs_geom_grid_ref = qgs_geom0.snappedToGrid(grid,grid)
            val0 = qgs_geom_grid_out.equals(qgs_geom_grid_ref)
            self.assertTrue(val0, title)


    def test_case27(self):
        title = "Test 27: Bend simplification with line smoothing but smoothing breaks spatial constraints"
        coord0 = [(-50,-25), (0,0), (0,-1), (3,-1), (3,0), (50,-25)]
        coord1 = [(1.5, .1), (1.5,3)]
        qgs_geom0 = create_line(coord0)
        qgs_geom1 = create_line(coord1)
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3, del_pol=True, del_hole=True, smooth_line=True)
        coord0 = [(-50,-25), (0,0), (3,0), (50,-25)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case28(self):
        title = "Test 28: Bend simplification with line smoothing but with start/end segment going in opposite direction"
        coord0 = [(-50,-25), (0,0), (0,-1), (3,-1), (3,0), (50,25)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True, smooth_line=True)
        coord0 = [(-50, -25), (0, 0), (1, 0.15579156685976017), (2, -0.15579156685976017), (3, 0), (50, 25)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case29(self):
        title = "Test 29: Bend simplification with line smoothing but with self intersection"
        coord0 = [(-50,-25), (0,0), (0,-1), (3,-1), (3,0), (50,25), (50,0.05), (-50,0.05)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True, smooth_line=True)
        coord0 = [(-50, -25), (0, 0), (3, 0), (50, 25), (50,0.05), (-50,0.05)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case30(self):
        title = "Test 30: Bend simplification with line smoothing but with start/end segment going in opposite direction"
        coord0 = [(-50, -25), (0, 0), (0, -1), (3, -1), (3, 0), (50, 25)]
        coord1 = [(.9, .1), (1.1, .1)]
        qgs_geom0 = create_line(coord0)
        qgs_geom1 = create_line(coord1)
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3, del_pol=True, del_hole=True, smooth_line=True)
        coord0 = [(-50, -25), (0, 0), (3, 0), (50, 25)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case31(self):
        title = "Test 31: Co-linear vertices at first/last vertice"
        coord0 = [(5,0), (0,0), (0,10), (5, 10), (10,10), (10,0), (5,0)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True, smooth_line=True)
        coord0 = [(0,0), (0,10), (10,10), (10,0), (0,0)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case32(self):
        title = "Test 31: Non simple line"
        coord0 = [(0,0), (5,0), (4,1), (6,1), (5,0), (10,0)]
        qgs_geom0 = create_line(coord0)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3)
        coord0 = [(0,0), (10,0)]
        qgs_geom0 = create_line(coord0)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case33(self):
        title = "Test 33: Normalization of in vector layer"
        print (title)
        vl = QgsVectorLayer("LineString", "temporary_polygon", "memory")
        pr = vl.dataProvider()
        fet = QgsFeature()
        fet.setId(1)
        qgs_line = QgsLineString((QgsPoint(0,0,0),QgsPoint(10,10,0),QgsPoint(20,20,0)))
        qgs_geom = QgsGeometry(qgs_line.clone())
        fet.setGeometry(qgs_geom)
        pr.addFeatures([fet])
        vl.updateExtents()
        feedback = QgsProcessingFeedback()
        qgs_features, geom_type = ReduceBend.normalize_in_vector_layer(vl, feedback)
        val0 = len(qgs_features) == 1
        qgs_geom = qgs_features[0].geometry()
        val1 = qgs_geom.wkbType() == QgsWkbTypes.LineString
        self.assertTrue(val0 and val1, title)


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
