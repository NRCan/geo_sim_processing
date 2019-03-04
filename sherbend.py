#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List

import fiona
import gdal

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
class GeoContent:
    """Contains the geographical content of the file.

        Keyword arguments:
        crs -- coordinate reference system
        driver -- name of the drive
        schemas -- dictionary of schema with "layer name" as key
        features -- list of geographic features in shapely structure

    """
    crs: str
    driver: str
    schemas: dict
    features: List[object]=None


@dataclass
class Params:
    """Contains the parameters .

            Keyword arguments:
            command -- parameter of the command line
            geo_content -- geographic content of the file

    """
    command: Command
    geo_content: GeoContent


#placer le nom du layer dans le feature shapely  et créer un dict de schema
#à l'écriture extraire tous les feature avec schemas identique et faire l'écriture dans le fichier par couche

command = Command (in_file='', out_file='', first_last=True, tolerance=10., simplicity=True,
                  sidedness=True, crossing=True, connection=True)
geo_content = GeoContent(crs=None, driver=None, schemas={}, features=[])
params = Params (command=command, geo_content=geo_content)

params.command.in_file = r'data\simple_file.gpkg'
params.command.out_file = r'data\simple_file_out.gpkg'

in_file = params.command.in_file


# Extract and load the layers of the file
layer_names = fiona.listlayers(in_file)
for layer_name in layer_names:
    with fiona.open(in_file, 'r', layer=layer_name) as src:
        params.geo_content.crs = src.crs
        params.geo_content.driver = src.driver
        params.geo_content.schemas[layer_name] = src.schema

        for in_feature in src:
            try:
                geom = in_feature['geometry']
                if geom['type'] == 'LineString':
                    feature = LineString(geom['coordinates'])
                    feature.layer_name = layer_name # Layer name is the key for the schema
                    feature.properties = in_feature['properties']
                    params.geo_content.features.append(feature)

            except:
                print ("Error processing feature) {0}".format(feature['id']))

    src.close()

layer_names = []
layer_names = [feature.layer_name for feature in params.geo_content.features if feature.layer_name not in layer_names]

for layer in params.layers:
    with fiona.open(params.command.out_file, 'w',
                    driver=params.geo_content.driver,
                    layer=layer.name,
                    crs=params.crs,
                    schema=schema) as dest:
        for feature in layer.features:
            out_feature = {'geometry': {'type': layer.type,
                                        'coordinates': list(feature.coords) },
                            'properties': feature.properties }
            dest.write(out_feature)

        dest.close()



