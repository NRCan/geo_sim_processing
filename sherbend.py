#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys, os
from dataclasses import dataclass
from typing import List
from algo_sherbend import AlgoSherbend


from lib_geosim import GenUtil


def managae_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input vector file to simplify")
    parser.add_argument("out_file", help="output vector file simplified")
    parser.add_argument("-eh", "--exclude_hole", action='store_true', help="exclude holes (interior) below minimum adjusted area")
    parser.add_argument("-ep", "--exclude_polygon", action='store_true', help="exclude polygons below minimum adjusted area")

    # Set exclusively mutual parameters
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-d", "--diameter", type=float, help="diameter of the minimum adjusted area bend to simplify")
    group.add_argument("-dl", "--dlayer", type=str, help="diameter of the minimum adjusted area bend to simplify per layer name (ex: -dl Road=5,Hydro=7.5)")

    # Read the command line parameter
    command = parser.parse_args()

    if command.dlayer is not None:
        # extract the diameter for each layer
        command.layers = {}
        # Split on the comas
        layers = command.dlayer.split(',')
        for layer in layers:
            # Split the layer name and the tolerance
            layer_tol = layer.split('=')
            try:
                command.layers[layer_tol[0]] = float(layer_tol[1])
            except:
                print ('Error in the definition of the diameter per layer: "{}"'.format(command.dlayer))
                parser.print_help()
                sys.exit(1)

    # Check that the input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))

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
    """Contains the parameters of the command.

        Keyword arguments:
        in_file -- name of the input file
        out_file -- name of the output file
        diameter -- diameter of the bend to simplify
        rotate_coord -- flag to enable/disable the rotation of closed line
        simplicity -- flag to enable/disable the test for OGC simple line constraint
        adjacency -- flag to enable/disable the test for adjacency constraint
        intersection -- flag to enable/disable the test for connection constraint
        add_vertex -- flag to enable/disable to add new vertex during bend simplification
        multi_bend -- flag to enable/disable the simplification of multi bends (more than one bend)
        verbose -- flag to enable/disable the verbose mode

        """
#    in_file: str
#    out_file: str
#    diameter: float
#    verbose: bool


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
    in_nbr_points: 0
    in_nbr_line_strings: 0
    in_nbr_polygons: 0
    in_nbr_holes: 0
    out_nbr_points: 0
    out_nbr_line_strings: 0
    out_nbr_polygons: 0
    out_nbr_holes: 0
    nbr_del_polygons: 0
    nbr_del_holes: 0
    bounds: List[object] = None
    layer_names: List[object] = None
    in_features: List[object] = None
    out_features: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, bounds=[], layer_names=[], in_features=[], out_features=[],
                         in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, in_nbr_holes=0,
                         out_nbr_points=0, out_nbr_line_strings=0, out_nbr_polygons=0, out_nbr_holes=0,
                         nbr_del_polygons=0, nbr_del_holes=0 )



# Read the command line arguments
command = managae_arguments()

# Extract and load the layers of the input file
GenUtil.read_in_file (command.in_file, geo_content)

# Set the diameter per layer
geo_content.layrer_diameter = {}
for layer_name in geo_content.layer_names:
    if command.layer_dimater.exists(layer_name):
        # set the value as defines in the command line
        geo_content.layrer_diameter[layer_name] = command.layer_dimater.exists(layer_name)
    else:



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
