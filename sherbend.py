#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List

import fiona

from shapely.geometry import LineString


@dataclass
class Command:
    in_file: str # Name of the input file
    out_file: str # Name of the output file
    first_last: bool # Flag to simplify first last vertice
    tolerance: float # Tolerance of simplification
    simplicity: bool # Flag to validate the simple feature
    sidedness: bool # Flag to validate the sidedness
    crossing: bool # Flag to validate the crossing
    connection: bool # Flag to verify the connection


@dataclass
class Layer:
    name: str # Name of the layer
    type: str # Type of feature
    crs: dict # Dictionary of the coordinate reference system
    schema: dict # Schema dictionary of the layer
    features: List[object]=None # List of features


@dataclass
class Params:
    command: Command
    layers: List[Layer]=None

command = Command(in_file='', out_file='', first_last=True, tolerance=10., simplicity=True,
                  sidedness=True, crossing=True, connection=True)
params = Params(command=command, layers=[])

params.command.in_file = r'data\simple_file.gpkg'
params.command.in_file = r'data\simple_file_out.gpkg'

in_file = params.command.in_file


# Extract and load the layers of the file
layers = fiona.listlayers(in_file)
shapely_features = []
for layer_name in layers:
    with fiona.open(in_file, 'r', layer=layer_name) as source:
        crs = source.crs
        driver = source.driver
        schema = source.schema

        for feature in source:
            try:
                geom = feature['geometry']
                properties = feature['properties']
                type = geom['type']
                coords = geom['coordinates']
                if type == 'LineString':
                    shapely_feature = LineString(coords)
                    shapely_feature.properties = properties
                    shapely_features.append(shapely_feature)

            except:
                print ("Error processing feature) {0}".format(feature['id']))
        layer = Layer(name=layer_name, type=type, crs=source.crs, schema=source.schema, features=shapely_features)
    params.layers.append(layer)
    source.close()

print (1)



