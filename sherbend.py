#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" """
from dataclasses import dataclass
from typing import List

import fiona

#from lib_geobato import Holder
from shapely.geometry import LineString

@dataclass
class Command:
    first_last: str
    tolerance: float
    features: List

command = Command("True", 25., [])
command.features.append('coco')
print (command)


in_file = r'data\simple_file.gpkg'

# Extract list of layers
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
    source.close()



