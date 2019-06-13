#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import List
from algo_sherbend import AlgoSherbend
from lib_geobato import PointSc, LineStringSc

import fiona

from shapely.geometry import Polygon


from lib_geobato import GenUtil
a = GenUtil.locate_bends1([(0,0),(1,1)])
b = GenUtil.locate_bends1([(0,0),(1,1),(2,2)])
c = GenUtil.locate_bends1([(0,0),(1,1),(2,0)])
d = GenUtil.locate_bends1([(0,0),(1,1),(2,1),(3,0)])
e = GenUtil.locate_bends1([(0,0),(1,1),(2,0),(3,0),(4,1),(5,0)])


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
    in_file: str
    out_file: str
    diameter: float
    rotate_coord: bool
    simplicity: bool
    sidedness: bool
    crossing: bool
    intersection: bool
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


command = Command (in_file='', out_file='', diameter=1.5, rotate_coord=True, simplicity=True,
                   sidedness=True, crossing=True, intersection=True, add_vertex=True, multi_bend=False, verbose=True)

geo_content = GeoContent(crs=None, driver=None, schemas={}, bounds=[], features=[])


#command.in_file = r'data\hydro_pol.shp'
command.in_file = r'data\simple_file3.gpkg'
command.out_file = r'data\simple_file_out3.gpkg'

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
                feature = PointSc(geom['coordinates'])
            elif geom['type'] == 'LineString':
                feature = LineStringSc(geom['coordinates'])
            elif geom['type'] == 'Polygon':
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

import time
start = time. time()
# Execute the Sherbend algorithm on the feature read
sherbend = AlgoSherbend(command, geo_content)
results = sherbend.process()
end = time. time()
print("Le temp est:{}: ".format(end - start))



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
