"""Topological feature simplifier

This algorithm is patching the known bug (1) of the shapely(2) simplify method

note 1: Known bug: https://locationtech.github.io/jts/javadoc/org/locationtech/jts/simplify/TopologyPreservingSimplifier.html

In order to patch the "known bug" once simplified if a feature is invalid, the algorithm is looking for the interiors located outside
of the exterior and delete them.

Note 2: Shapely is a python wrapper of GEOS which is a translation fo the Java Topological Suite (JTS)

"""

from argparse import ArgumentParser
from os import path
from dataclasses import dataclass
from typing import List
from shapely.geometry import LineString, Polygon, MultiPolygon
from shapely.validation import explain_validity
from lib_geosim import GenUtil, Holder


def manage_arguments():
    """Extract the parameters of the line of command

    Parameters
    ----------
    None

    Returns
    -------
    command
        Parameters of the command line
    """

    # Setting the parameters of the command line
    parser = ArgumentParser()
    parser.add_argument("in_file", help="input vector file to simplify")
    parser.add_argument("out_file", help="output vector file simplified")

    # Set exclusively mutual parameters
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--tolerance", type=float, help="tolerance for the simplification for all the layers in the file")
    group.add_argument("-tl", "--tlayer", type=str, help="folerance for the simplification per layer name (ex: -tl Road=5,Hydro=7.5")

    # Read the command line parameter
    command = parser.parse_args()

    if command.tlayer is not None:
        # extract the diameter for each layer
        command.dlayer_dict = {}
        # Split on the comas
        layers = command.tlayer.split(',')
        for layer in layers:
            # Split the layer name and the tolerance
            layer_tol = layer.split('=')
            try:
                command.dlayer_dict[layer_tol[0]] = float(layer_tol[1])
            except Exception:
                print('Error in the definition of the tolerance per layer (-tl or --tlayer): "{}"'.format(command.dlayer))
                parser.print_help()
                exit(1)

    # Check that the input file exist. Exit if missing
    if not path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))

    # Check that the output file exist. Exit if present
    if path.isfile(command.out_file):
        raise Exception('Output file is present: {}'.format(command.out_file))

    return command

def polygon_area(exterior):
    return Polygon(exterior).area


def topo_simplifier(tolerance, in_feature, geo_content):
    """Simplify each line with the shapely simplify method

    If after simplification the feature is invalid it's because a hole of the simplied polygon is now located outside the polygon.
    This method is used to locate and delete these holes located outside the simplified feature.
    Should only happen on Polygon feature.

    Parameters
    ----------
    tolerance : float
        The tolerance for the feature simplification
    in_feature : LineString or Polygon
        Feature to simplify
    geo_content: data class
        Statistic of the line simplification

    Returns
    -------
    LineString or Polygon
        Simplified feature
    """

    out_feature = in_feature.simplify(tolerance, preserve_topology=True)
    if out_feature.is_valid and out_feature.is_simple:
        # Simplification OK
        pass
    else:
        if out_feature.geom_type == GenUtil.POLYGON:
            # Extract outer and inner rings
            pol_rings = [out_feature.exterior] + list(out_feature.interiors)
            pol_rings.sort(key=polygon_area)  # Sort by area size so process them by size order
            polygons = []  # Output polygons
            pol_ring_ext = pol_rings.pop()
            polygon = Holder(exterior=pol_ring_ext, interiors=[])
            polygons.append(polygon)
            exterior = pol_ring_ext
            interiors = []
            # Check if each polygon ring is located in the polygon by adding a smaller ring at each loop
            while pol_rings:  # Loop until no more rings
                ring_interior = pol_rings.pop()  # get the next ring to add
                polygon = Polygon(exterior, interiors)  # Recreate the polygon
                # Test if the ring is located inside the polygon
                if ring_interior.within(polygon):
                    # Add this ring
                    interiors.append(ring_interior)
                else:
                    # Ring located outside the polygon... delete it
                    geo_content.nbr_del_holes += 1

            # Recreate the polygon
            out_feature = Polygon(exterior, interiors)
            # check if the simplification is OK
            if out_feature.is_valid and out_feature.is_simple:
                # Everything is OK now
                pass
            else:
                # should not happen... but
                error_txt = explain_validity(out_feature)
                print ('Error: {0}'.format(error_txt))
        else:
            # Should not happen... but
            error_txt = explain_validity(out_feature)
            print('Error: {0}'.format(error_txt))

    return out_feature


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
    nbr_del_holes: 0
    out_nbr_polygons : 0
    out_nbr_holes : 0
    out_nbr_points : 0
    out_nbr_line_strings : 0
    bounds: List[object] = None
    layer_names: List[object] = None
    in_features: List[object] = None
    out_features: List[object] = None


geo_content = GeoContent(crs=None, driver=None, schemas={}, bounds=[], layer_names=[], in_features=[], out_features=[],
                         in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, nbr_del_holes=0, out_nbr_polygons=0,
                         out_nbr_holes=0, out_nbr_points=0, out_nbr_line_strings=0)

# Read the command line arguments
command = manage_arguments()

if command.tlayer:
    # Extract the list of layers to read
    in_layer_names = [layer_name for layer_name in command.dlayer_dict.keys()]
else:
    # Read all the layers in the input file
    in_layer_names = None

# Read and load the layers of the input file
GenUtil.read_in_file(command.in_file, geo_content, in_layer_names)

# Set the tolerance for each layer
tmp_tlayer_dict = {}
for layer_name in geo_content.layer_names:
    if command.tolerance is not None:
        # The same tolerance value is applied to all the layers
        tmp_tlayer_dict[layer_name] = command.tolerance
    else:
        # There is a specific tolerance for each layer
        tmp_tlayer_dict[layer_name] = command.tlayer_dict[layer_name]

# Reset the value of tlayer
command.tlayer_dict = tmp_tlayer_dict

# Only keep in the in features the ones where the diameter is not -1
tmp_in_features = []
for feature in geo_content.in_features:
    if command.tlayer_dict[feature.sb_layer_name] != -1:
        tmp_in_features.append(feature)
geo_content.in_features = tmp_in_features


# Print statistics of the process
print("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of output file: {}".format(command.out_file))
print("Number of layers read: {}".format(len(geo_content.schemas)))
for layer_name, tolerance in command.tlayer_dict.items():
    print("   - Layer name: {} with tolerance: {}".format(layer_name, tolerance))
print("Number of features read: {}".format(len(geo_content.in_features)))


geo_content.out_features = []
for in_feature in geo_content.in_features:
    if in_feature.geom_type == GenUtil.POINT:
        out_feature = in_feature
    elif in_feature.geom_type in [GenUtil.LINE_STRING, GenUtil.POLYGON]:
        out_feature = topo_simplifier(command.tolerance, in_feature, geo_content)
    else:
        raise GenUtil.GeoSimException("Cannot simplify type: {0}".format(in_feature.geom_type))
    out_feature.sb_properties = in_feature.sb_properties
    out_feature.sb_layer_name = in_feature.sb_layer_name
    geo_content.out_features.append(out_feature)

# Print some processing stats
print("Number of holes excluded: {}".format(geo_content.nbr_del_holes))
print("Number of features written: {}".format(len(geo_content.out_features)))
print("-----")

# Copy the results in the output file
GenUtil.write_out_file(command.out_file, geo_content)

