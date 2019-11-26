#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys, os
from dataclasses import dataclass
from typing import List
from shapely.geometry import Point, LineString, Polygon


from lib_geosim import GenUtil, ChordalAxis, LineStringSc


def managae_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input vector file")
    parser.add_argument("-p", "--polygon", type=str, help="input layer name containing the polygon")
    parser.add_argument("-t", "--tesselation", type=str, help="input layer name containing the result of the tesselation (triangle)")
    parser.add_argument("-a", "--attribute", type=str, help="name of attribute linking the triangle to the polygon")
    parser.add_argument("-s", "--skeleton", type=str, help="name of output skeleton layer (centre line)")

    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))


    return command

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
    in_features: List[object]
    out_features: List[object]
    out_nbr_points: 0
    out_nbr_line_strings: 0
    out_nbr_polygons: 0
    bounds: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, in_features=[], out_features=[],
                         in_nbr_triangles=0, in_nbr_polygons=0, bounds=[], out_nbr_points=0,
                         out_nbr_line_strings=0, out_nbr_polygons=0)



# Read the command line arguments
command = managae_arguments()

# Extract and load the layers of the input file
layers = [command.polygon, command.tesselation]
GenUtil.read_in_file (command.in_file, geo_content, layers)

polygon_dict = {}
triangle_dict = {}
for in_feature in geo_content.in_features:
    key = in_feature.sb_properties[command.attribute]
    if in_feature.sb_layer_name == command.polygon:
        polygon_dict[key] = in_feature
    else:
        if key in triangle_dict.keys():
            triangle_dict[key].append(in_feature)
        else:
            triangle_dict[key] = [in_feature]

geo_content.in_features = None


a = LineStringSc([(0,0), (1,1), (2,2), (0,0)])
case0 = [a]

a = LineStringSc([(0,0), (1,1), (2,0), (0,0)])
b = LineStringSc([(1,1), (3,1), (2,0), (1,1)])
case1 = [a,b]

a = LineStringSc([(0,0), (1,1), (2,0), (0,0)])
b = LineStringSc([(1,1), (3,1), (2,0), (1,1)])
c = LineStringSc([(2,0), (3,1), (4,0), (2,0)])
case2 = [a,b,c]

a = LineStringSc([(0,0), (1,1), (2,0), (0,0)])
b = LineStringSc([(1,1), (3,1), (2,0), (1,1)])
c = LineStringSc([(2,0), (3,1), (4,0), (2,0)])
d = LineStringSc([(1,1), (2,2), (3,1), (1,1)])
case2 = [a,b,c,d]

#triangle_dict = {1:case2}

for key in triangle_dict.keys():
    ca = ChordalAxis(triangle_dict[key], GenUtil.ZERO)
    centre_lines = ca.get_skeletton()
    # Store the chordal axis in the output
    for centre_line in centre_lines:
        centre_line.sb_layer_name = command.skeleton
        centre_line.sb_properties={}
        geo_content.out_features.append(centre_line)

print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of input polygon layer: {}".format(command.polygon))
print ("Name of input tesselation layer: {}".format(command.tesselation))
print ("Nampe of output skeleton layer: {}".format(command.skeleton))
print ("-----")


# Copy the results in the output file
geo_content.layer_names = [command.skeleton]
GenUtil.write_out_file_append (command.in_file, geo_content)

print ("Number of point features written: {}".format(geo_content.out_nbr_points))
print ("Number of line string features written: {}".format(geo_content.out_nbr_line_strings))
print ("Number of polygon written: {}".format(geo_content.out_nbr_polygons))
