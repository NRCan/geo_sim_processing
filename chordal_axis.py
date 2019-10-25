#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys, os
from dataclasses import dataclass
from typing import List
from shapely.geometry import Point, LineString, Polygon


from lib_geosim import GenUtil, ChordalAxis


def managae_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("out_file", help="output vector file simplified")
    parser.add_argument("-t", "--triangle", action='store_true', help="input file containing the result of the triangulation")
    parser.add_argument("-p", "--polygon", action='store_true', help="input file containing the polygon")
    parser.add_argument("-a", "--attribute", action='store_true', help="name of attribute linking the triangle to the polygon")


    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not os.path.isfile(command.triangle):
        raise Exception('Input triangle file is missing: {}'.format(command.in_file))

    # Check that the polygon input file exist. Exit if missing
    if not os.path.isfile(command.triangle):
        raise Exception('Input triangle file is missing: {}'.format(command.in_file))

    # Check that the output file exist. Exit if present
    if os.path.isfile(command.out_file):
        raise Exception('Output file is present: {}'.format(command.out_file))

    return command

"""
a = LinearRing(((0,0),(1,1),(2,0)))
b = a.is_ccw
a = LinearRing(((2,0),(1,1),(0,0)))
b = a.is_ccw


a = Polygon((((1,1), (2,2), (2,0), (0,0), (0,2), (1,1))))
a = orient(Polygon(a.exterior.coords), GenUtil.ANTI_CLOCKWISE) # Orient line clockwiswe
a = LineStringSb(a.exterior.coords)
a.simplify(15)

a = LineStringSb(((0,0), (0,3), (1.5,2.5), (3,3), (3,0), (0,0) ))
a.simplify(5)

a = Polygon (( (1647625.889999593, 195454.0860011009),\
(1647630.371999593, 195435.4470011005),\
(1647640.775999592, 195439.1370011028),\
(1647649.498999593, 195447.547001102),\
(1647644.202999593, 195459.5080011021),\
(1647638.619999593, 195469.246001103),\
(1647618.486999592, 195492.9860011013),\
(1647623.151999593, 195464.8150011022),\
(1647625.889999593, 195454.0860011009)))
a.simplify(1.5)

a = LineStringSb(((0,0), (2,2)))
a.simplify(5)

a = LineStringSb(((0,0), (1,1), (2,2)))
a.simplify(5)

a = LineStringSb(((0,0), (1,1), (2,0)))
a.simplify(5)


a = LineStringSb(((0,0), (1,1), (2,0), (3,2), (4,0), (5,3), (6,0), (7,.1), (8,0)))
a.simplify(5)

# Closed star
a = LineStringSb(((0,0), (0,3), (1.5,2.5), (3,3), (2.5,1.5), (3,0), (0,0) ))
a.simplify(5)

# Closed star
a = LineStringSb(((0,0), (0,3), (1.5,2.5), (3,3), (3,0), (0,0) ))
a.simplify(5)




a = LineStringSb(((0,0), (1,1), (2,1), (3,0)))
a.simplify(5)

a = LineStringSb(((0,0), (1,1), (2,0), (3,1)))
a.simplify(5)

a = LineStringSb(((0,0), (1,1), (2,0), (3,1), (4,0)))
a.simplify(5)

a = LineStringSb(((0,0), (1,1), (2,0), (0,0)))
a.simplify(5)

a = LineStringSb(((0,0), (0,2), (1,1), (2,2), (2,0), (0,0)))
a.simplify(5)

a = LineStringSb((((0,2), (1,1), (2,2), (2,0), (0,0), (0,2))))
a.simplify(5)

a = LineStringSb((((1,1), (2,2), (2,0), (0,0), (0,2), (1,1))))
a.simplify(5)

a = LineStringSb((((2,2), (2,0), (0,0), (0,2), (1,1), (2,2))))
a.simplify(5)

a = LineStringSb((((2,0), (0,0), (0,2), (1,1), (2,2), (2,0))))
a.simplify(5)

a = LineStringSb((( (0,0),(0,3),(1,2),(3,3),(3,0),(1,1),(0,0)) ))
a.simplify(5)
"""


@dataclass
class Command:
    """Contains the parameters of the command."""


@dataclass
class GeoContent:
    """Contains the geographical content of the file.

        Keyword arguments:
        crs -- coordinate reference system
        driver -- name of the drive
        schemas -- dictionary of schema with "layer name" as key
        in_features -- input list of geographic features in shapely structure
        layer_names -- Name of the layers in the spatial file

    """
    crs: None
    driver: None
    schemas: dict
    in_nbr_triangles: 0
    in_nbr_polygons: 0
    in_triangles: List[object] = None
    in_polygons: List[object] = None
    out_features: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, in_triangles=[], in_polygons = [], out_features=[],
                         in_nbr_triangles=0, in_nbr_polygons=0)



# Read the command line arguments
###command = managae_arguments()

# Extract and load the layers of the input file
###GenUtil.read_in_file (command.in_triangle, geo_content)


coords = ((0,0),(0,3),(1.5,2),(3,3),(3,0),(1.5,1),(0,0))
pol = Polygon(coords)
l1 = ((0,0),(0,3),(1.5,2),(0,0))
t1 = LineString(l1)
l2 = ((1.5,2),(1.5,1),(0,0),(1.5,2))
t2 = LineString(l2)
l3 = ((1.5,2),(1.5,1),(3,0),(1.5,2))
t3 = LineString(l3)
l4 = ((1.5,2),(3,3),(3,0),(1.5,2))
t4 = LineString(l4)

a = ChordalAxis(pol, [t1,t2,t3,t4], 0.0)
center_line = a.get_skeletton()


print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of output file: {}".format(command.out_file))
print ("Number of layers read: {}".format(len(geo_content.schemas)))
print ("Number of features read: {}".format(len(geo_content.in_features)))
print ("-----")

# Execute the Sherbend algorithm on the feature read
sherbend = AlgoSherbend(command, geo_content)
sherbend.process()

# Copy the results in the output file
GenUtil.write_out_file (command.out_file, geo_content)

print ("Number of polygons excluded: {}".format(geo_content.nbr_del_polygons))
print ("Number of holes excluded: {}".format(geo_content.nbr_del_holes))
print ("Number of point features written: {}".format(geo_content.out_nbr_points))
print ("Number of line string features written: {}".format(geo_content.out_nbr_line_strings))
print ("Number of polygon written: {}".format(geo_content.out_nbr_polygons))
print ("Number of holes written: {}".format(geo_content.out_nbr_holes))
