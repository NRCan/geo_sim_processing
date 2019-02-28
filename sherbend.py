#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List

import fiona

from shapely.geometry import LineString


@dataclass
class Command:
    """Contains the parameters of the command.

        Keyword arguments:
        in_file -- name of the input file
        out_file -- name of the output file
        first_last -- flag to simplify the first/last bend
        simplicity -- flag to test for OGC simple line constraint
        adjacency/sidedness/consistency -- flag to test for consistency constraint
        connection -- flag to test connection constraint

        """
    in_file: str
    out_file: str
    first_last: bool
    tolerance: float
    simplicity: bool
    sidedness: bool
    crossing: bool
    connection: bool


@dataclass
class Layer:
    """Contains the parameters of a layer.

        Keyword arguments:
        name -- name of the layer
        type -- geometry type
        schema -- layer schema
        features -- list of geographic features

        """
    name: str
    type: str
    schema: dict
    features: List[object]=None


@dataclass
class Params:
    """Contains the parameters .

            Keyword arguments:
            command -- parameter of the command line
            crs -- coordinate reference system
            driver -- layer schema
            features -- list of geographic features

            """
    command: Command
    crs: str
    driver: str


placer le nom du layer dans le feature shapely  et créer un dict de schema
à l'écriture extraire tous les feature avec schemas identique et faire l'écriture dans le fichier par couche

command = Command (in_file='', out_file='', first_last=True, tolerance=10., simplicity=True,
                  sidedness=True, crossing=True, connection=True)
params = Params (command=command, crs=None, driver=None)

params.command.in_file = r'data\simple_file.gpkg'
params.command.out_file = r'data\simple_file_out.gpkg'

in_file = params.command.in_file


# Extract and load the layers of the file
layer_names = fiona.listlayers(in_file)
features = []
for layer_name in layer_names:
    with fiona.open(in_file, 'r', layer=layer_name) as source:
        params.crs = source.crs
        params.driver = source.driver
        schema = source.schema

        for in_feature in source:
            try:
                geom = in_feature['geometry']
                if geom['type'] == 'LineString':
                    feature = LineString(geom['coordinates'])
                    feature.properties = in_feature['properties']
                    features.append(feature)

            except:
                print ("Error processing feature) {0}".format(feature['id']))
        layer = Layer(name=layer_name, type=type, schema=source.schema, shapely_features=features)
    params.layers.append(layer)
    source.close()

for layer in params.layers:
    with fiona.open(params.command.out_file, 'w',
                    driver=params.driver,
                    layer=layer.name,
                    crs=params.crs,
                    schema=schema) as dest:
        for feature in layer.features:
            out_feature = {'geometry': {'type': layer.type,
                                        'coordinates': list(feature.coords) },
                            'properties': feature.properties }
            dest.write(out_feature)

        dest.close()



