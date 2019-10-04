#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from dataclasses import dataclass
from typing import List
from algo_sherbend import AlgoSherbend, LineStringSb, PointSb

import fiona

from shapely.geometry import Polygon

from lib_geosim import GenUtil

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
        features -- list of geographic features in shapely structure

    """
    crs: None
    driver: None
    schemas: dict
    bounds: List[object] = None
    features: List[object] = None


#command = Command (in_file='', out_file='', diameter=50, rotate_coord=True, simplicity=True,
#                   sidedness=True, crossing=True, intersection=True, add_vertex=True, multi_bend=False, verbose=True)

geo_content = GeoContent(crs=None, driver=None, schemas={}, bounds=[], features=[])


# Reading the parameter on the command line
parser = argparse.ArgumentParser()
parser.add_argument("in_file", help="input vector file to simplify")
parser.add_argument("out_file", help="output vector file simplified")
parser.add_argument("-d", "--diameter", type=float, help="diameter of the bend to simplify")
command = parser.parse_args()
print (command.in_file)

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
                feature = PointSb(geom['coordinates'])
            elif geom['type'] == 'LineString':
                feature = LineStringSb(geom['coordinates'])
            elif geom['type'] == 'Polygon':
                exterior = geom['coordinates'][0]
                interiors = geom['coordinates'][1:]
                feature = Polygon(exterior, interiors)
            else:
                print ("The following geometry type is unsupported: {}".format(geom['type']))
                feature = None
            if feature is not None:
                feature.sb_layer_name = layer_name  # Layer name is the key for the schema
                feature.sb_properties = in_feature['properties']
                geo_content.features.append(feature)
    src.close()


print ("-----")
print("Name of input file: {}".format(command.in_file))
print("Name of output file: {}".format(command.out_file))
print ("Number of layers read: {}".format(len(geo_content.schemas)))
print ("Number of features read: {}".format(len(geo_content.features)))
print ("-----")

# Execute the Sherbend algorithm on the feature read
sherbend = AlgoSherbend(command, geo_content)
results = sherbend.process()




# Extract the unique name of each layer
layer_names = set()
for feature in results:
    layer_names.add(feature.sb_layer_name)

# Loop over each layer and write the content of the file
for layer_name in layer_names:
    with fiona.open(command.out_file, 'w',
                    driver=geo_content.driver,
                    layer=layer_name,
                    crs=geo_content.crs,
                    schema=geo_content.schemas[layer_name]) as dest:
        for feature in (feature for feature in results if feature.sb_layer_name==layer_name):
            # Transform the Shapely features for fiona writing
            if feature.geom_type == 'Point' or feature.geom_type == 'LineString':
                coordinates = list(feature.coords)
            elif feature.geom_type == 'Polygon':
                exterior = list(feature.exterior.coords)
                interior = [list(interior.coords) for interior in feature.interiors]
                coordinates = [exterior]+interior

            out_feature = {'geometry': {'type': feature.geom_type,
                                        'coordinates': coordinates},
                           'properties': feature.sb_properties}
            dest.write(out_feature)

    dest.close()

print ("Number of features written: {}".format(len(results)))
