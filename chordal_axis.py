"""Algorithm to extract the ChordalAxis of a set of triangles"""

from argparse import ArgumentParser
from os import path
from dataclasses import dataclass
from typing import List
from shapely.geometry import LineString
from lib_geosim import GenUtil, ChordalAxis


def managae_arguments():
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
    parser.add_argument("in_file", help="input vector file")
    parser.add_argument("-t", "--triangle", type=str, help="input layer name containing the result of the triangulation (triangles)")
    parser.add_argument("-s", "--skeleton", type=str, help="name of output skeleton layer (centre line)")
    parser.add_argument("-c", "--correction", action='store_true', help="name of output skeleton layer (centre line)")

    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))

    return command


@dataclass
class GeoContent:
    """Data class containing the geographical content of the file.

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
    in_nbr_triangles: 0
    in_features: List[object]
    out_features: List[object]
    out_nbr_points: 0
    out_nbr_line_strings: 0
    out_nbr_polygons: 0
    bounds: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, in_features=[], out_features=[],
                         in_nbr_triangles=0, bounds=[], out_nbr_points=0,
                         out_nbr_line_strings=0, out_nbr_polygons=0)

# Read the command line arguments
command = managae_arguments()

# Extract and load the layers of the input file
layers = [command.triangle]
GenUtil.read_in_file (command.in_file, geo_content, layers)

lst_triangles = []
for in_feature in geo_content.in_features:
    lst_triangles.append(in_feature)

geo_content.in_features = None  # Reset in_features

ca = ChordalAxis(lst_triangles, GenUtil.ZERO)
if command.correction:
    # Correct the skeleton
    ca.correct_skeleton()
centre_lines = ca.get_skeleton()

# Store the chordal axis in the output file
for centre_line in centre_lines:
    centre_line.sb_layer_name = command.skeleton
    centre_line.sb_properties={}
    geo_content.out_features.append(centre_line)

# Print the stats
print ("-------")
print("Name of input file: {}".format(command.in_file))
print ("Name of input triangle layer: {}".format(command.triangle))
print ("Nampe of output skeleton layer: {}".format(command.skeleton))
print ("Number of polygons: {}".format(ca.nbr_polygons))
print ("Number of triangles: {}".format(ca.nbr_triangles))
print ("Number of line pruned during skeleton correction: {}".format(ca.nbr_lines_pruned))
print ("Number of iteration done for skeleton correction: {}".format(ca.nbr_iteration))
print ("Number of T junction corrected: {}".format(ca.nbr_t_junction))
print ("Number of X junction corrected: {}".format(ca.nbr_x_junction))
print ("-----")

# Copy the results in the output file
geo_content.layer_names = [command.skeleton]
GenUtil.write_out_file_append (command.in_file, geo_content)
