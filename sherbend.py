#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys, os
from dataclasses import dataclass
from typing import List
from algo_sherbend import AlgoSherbend
from lib_geosim import GenUtil


def manage_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="input vector file to simplify")
    parser.add_argument("out_file", help="output vector file simplified")
    parser.add_argument("-eh", "--exclude_hole", action='store_true', help="exclude holes (interior) below minimum adjusted area")
    parser.add_argument("-ep", "--exclude_polygon", action='store_true', help="exclude polygons below minimum adjusted area")
    parser.add_argument("-pl", "--per_layer", action='store_true', help="evaluate topology per layer only (feature from different layers can overlap " +
                                                                                "after simplification)")

    # Set exclusively mutual parameters
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-d", "--diameter", type=float, help="diameter of the minimum adjusted area bend to simplify")
    group.add_argument("-dl", "--dlayer", type=str, help="diameter of the minimum adjusted area bend to simplify per layer name (ex: -dl Road=5,Hydro=7.5")

    # Read the command line parameter
    command = parser.parse_args()

    if command.dlayer is not None:
        # extract the diameter for each layer
        command.dlayer_dict = {}
        # Split on the comas
        layers = command.dlayer.split(',')
        for layer in layers:
            # Split the layer name and the tolerance
            layer_tol = layer.split('=')
            try:
                command.dlayer_dict[layer_tol[0]] = float(layer_tol[1])
            except Exception:
                print('Error in the definition of the diameter per layer: "{}"'.format(command.dlayer))
                parser.print_help()
                sys.exit(1)

    # Check that the input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))

    # Check that the output file exist. Exit if present
    if os.path.isfile(command.out_file):
        raise Exception('Output file is present: {}'.format(command.out_file))

    return command


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
                         nbr_del_polygons=0, nbr_del_holes=0)

# Read the command line arguments
command = manage_arguments()

if command.dlayer:
    # Extract the list of layers
    in_layer_names = [layer_name for layer_name in command.dlayer_dict.keys()]
else:
    # Read all the layers in the file
    in_layer_names = None

# Extract and load the layers of the input file
GenUtil.read_in_file(command.in_file, geo_content, in_layer_names)

# Set the diameter for each layer
tmp_dlayer_dict = {}
for layer_name in geo_content.layer_names:
    if command.diameter is not None:
        # The same diameter value is applied to all the layers
        tmp_dlayer_dict[layer_name] = command.diameter
    else:
        # There is a specific diameter for each layer
        tmp_dlayer_dict[layer_name] = command.dlayer_dict[layer_name]

# Reset the value of dlayer
command.dlayer_dict = tmp_dlayer_dict

# Only keep in the in features the ones where the diameter is not -1
tmp_in_features = []
for feature in geo_content.in_features:
    if command.dlayer_dict[feature.sb_layer_name] != -1:
        tmp_in_features.append(feature)
geo_content.in_features = tmp_in_features

print("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of output file: {}".format(command.out_file))
print("Number of layers read: {}".format(len(geo_content.schemas)))
for layer_name, diameter in command.dlayer_dict.items():
    print("   - Layer name: {} with diamater: {}".format(layer_name, diameter))
print("Exclude polygon below minimum diameter: {}".format(str(command.exclude_polygon)))
print("Exclude polygon interior (hole) below minimum diameter: {}".format(str(command.exclude_hole)))
print("Number of features read: {}".format(len(geo_content.in_features)))
print("-----")

# Execute the Sherbend algorithm on the feature read
if command.per_layer:
    # Process each layer independently (do not check topology between layer)
    tmp_in_features = geo_content.in_features
    tmp_out_features = []
    for layer_name in command.dlayer_dict.keys():
        geo_content.in_features = []
        geo_content.out_features = []
        for feature in tmp_in_features:
            if feature.sb_layer_name == layer_name:
                geo_content.in_features.append(feature)

        # Process feature per layer in one pprocess
        sherbend = AlgoSherbend(command, geo_content)
        sherbend.process()
        # Copy the feature in the temporary outut
        for feature in geo_content.out_features:
            tmp_out_features.append(feature)

    geo_content.out_features = tmp_out_features
else:
    # Process all layer in one pprocess (check topology between layer)
    sherbend = AlgoSherbend(command, geo_content)
    sherbend.process()

# Print some processing stats
print("Number of polygons excluded: {}".format(geo_content.nbr_del_polygons))
print("Number of holes excluded: {}".format(geo_content.nbr_del_holes))
print("Number of point features written: {}".format(geo_content.out_nbr_points))
print("Number of line string features written: {}".format(geo_content.out_nbr_line_strings))
print("Number of polygon written: {}".format(geo_content.out_nbr_polygons))
print("Number of holes written: {}".format(geo_content.out_nbr_holes))

# Copy the results in the output file
GenUtil.write_out_file(command.out_file, geo_content)

