#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List
from algo_sherbend import AlgoSherbend

import fiona

from shapely.geometry import Point, LineString, Polygon


@dataclass
class Command:
    """Contains the parameters of the command.

        Keyword arguments:
        in_file -- name of the input file
        out_file -- name of the output file
        diameter -- diameter of the bend to simplify
        first_last -- flag to enable/disable the simplification the first/last bend
        simplicity -- flag to enable/disable the test for OGC simple line constraint
        adjacency -- flag to enable/disable the test for adjacency constraint
        connection -- flag to enable/disable the test for connection constraint
        add_vertex -- flag to enable/disable to add new vertex during bend simplification
        multi_bend -- flag to enable/disable the simplification of multi bends (more than one bend)
        verbose -- flag to enable/disable the verbose mode

        """
    in_file: str
    out_file: str
    simplify_first_last: bool
    diameter: float
    simplicity: bool
    adjacency: bool
    crossing: bool
    connection: bool
    add_vertex: bool
    multi_bend: bool
    verbose: bool


@dataclass
class GeoContent:
    """Contains the geographical content of the file.

        Keyword arguments:
        crs -- coordinate reference system
        driver -- name of the drive
        schemas -- dictionary of schema with "layer name" as key
        features -- list of geographic features in shapely structure

    """
    crs: None
    driver: None
    schemas: dict
    bounds: List[object] = None
    features: List[object] = None


command = Command (in_file='', out_file='', simplify_first_last=True, diameter=2., simplicity=True,
                   adjacency=True, crossing=True, connection=True, add_vertex=True, multi_bend=False, verbose=True)

geo_content = GeoContent(crs=None, driver=None, schemas={}, bounds=[], features=[])


command.in_file = r'data\simple_file.gpkg'
command.out_file = r'data\simple_file_out.gpkg'

# Extract and load the layers of the file
layer_names = fiona.listlayers(command.in_file)
for layer_name in layer_names:
    with fiona.open(command.in_file, 'r', layer=layer_name) as src:
        geo_content.crs = src.crs
        geo_content.driver = src.driver
        geo_content.schemas[layer_name] = src.schema
        geo_content.bounds.append(src.bounds)

        for in_feature in src:
            geom = in_feature['geometry']
            if geom['type'] == 'Point':
                feature = Point(geom['coordinates'])
            elif geom['type'] == 'LineString':
                feature = LineString(geom['coordinates'])
            elif geom['type'] == 'Polygon':
                print (geom['coordinates'])
                exterior = geom['coordinates'][0]
                interiors = geom['coordinates'][1:]
                feature = Polygon(exterior, interiors)
            else:
                print ("The following geometry type is unsupported: {}".format(geom['type']))
            feature._gbt_layer_name = layer_name  # Layer name is the key for the schema
            feature._gbt_properties = in_feature['properties']
            geo_content.features.append(feature)
    src.close()


print ("-----")
print("Name of input file: {}".format(command.in_file))
print ("Number of layers read: {}".format(len(geo_content.schemas)))
print ("Number of features read: {}".format(len(geo_content.features)))

# Execute the Sherbend algorithm on the feature read
sherbend = AlgoSherbend(command, geo_content)
results = sherbend.process()

# Extract the name of each layer
layer_names = set()
for feature in results:
    layer_names.add(feature._gbt_layer_name)

# Loop over each layer and write the content of the file
for layer_name in layer_names:
    with fiona.open(command.out_file, 'w',
                    driver=geo_content.driver,
                    layer=layer_name,
                    crs=geo_content.crs,
                    schema=geo_content.schemas[layer_name]) as dest:
        for feature in (feature for feature in results if feature._gbt_layer_name==layer_name):
            # Transform the Shapely features for fiona writing
            if feature.geom_type == 'Point' or feature.geom_type == 'LineString':
                coordinates = list(feature.coords)
            elif feature.geom_type == 'Polygon':
                exterior = list(feature.exterior.coords)
                interior = [list(interior.coords) for interior in feature.interiors]
                coordinates = [exterior]+interior


            out_feature = {'geometry': {'type': feature.geom_type,
                                        'coordinates': coordinates},
                            'properties': feature._gbt_properties}
            dest.write(out_feature)

        dest.close()
