#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, os
from dataclasses import dataclass
from typing import List
from lib_geosim import SpatialContainer, GenUtil
from shapely.geometry import Point, LineString, Polygon



def read_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="geopackage vector file to clean")
    parser.add_argument("-y", "--yjunction", type=str, default=0.0, help="max tolerance for cleaning Y junction ")
    parser.add_argument("-x", "--xjunction", type=str, help="max tolerance for cleaning X junction")
    parser.add_argument("-j", "--join", type=str, help="max tolerance for joining 2 roads")
    parser.add_argument("-r", "--remove", type=str, help="max tolerance for removing noise")
    parser.add_argument("-il", "--input", type=str, help="name input layer containing the roads")
    parser.add_argument("-ol", "--output", type=str, help="name output layer containing the roads")

    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))


    return command


def build_topology(center_lines, search_tolerance):


    # Load the features in the spatial container (accelerate the search)
    s_container = SpatialContainer()
    for center_line in center_lines:
        center_line.sb_geom_type = GenUtil.LINE_STRING
        s_container.add_feature(center_line)

    # Build topology. For each create list of connecting lines
    for line in center_lines:
        p_start = line.coords[0]
        p_end = line.coords[-1]
        b_box = GenUtil.build_bounding_box(search_tolerance, p_start)
        lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
        line.start_lines = []
        for line_b_box in lines_b_box:
            if line.touches(line_b_box):
                line.start_lines.append(line_b_box)
        b_box = GenUtil.build_bounding_box(search_tolerance, p_end)
        lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
        line.end_lines = []
        for line_b_box in lines_b_box:
            if line.touches(line_b_box):
                line.end_lines.append(line_b_box)

    return s_container

def clean_x_junction(s_container, command, geo_content):

    # Loop over each line in the container
    for x_line in s_container.get_features():
        # Only select line below the tolerance and linked to 2 other lines at each extremity
        if x_line.length <= command.xjunction and \
           len(x_line.start_lines) == 2 and \
           len(x_line.end_lines) == 2:
            mid_point = x_line.interpolate(.5, normalized=True) # Find the new intersection point
            p_start = Point(x_line.coords[0])
            p_end = Point(x_line.coords[-1])
            lst_lines = x_line.start_lines+x_line.end_lines
            for line in lst_lines:
                if mid_point.distance(p_start) <= GenUtil.ZERO or \
                   mid_point.distance(p_end) <= GenUtil.ZERO:
                    # Move the first coordinate to the new mid point
                    i = 0
                else:
                    # Move the last coordinate to the new mid point
                    i = len(line.coords)-1
                # reset the coordinate
                lst_coords = list(line.coords)
                lst_coords[i] = (mid_point.x,mid_point.y)
                line.coords = lst_coords
                geo_content.nbr_xjunction += 1 # Add stats counter

            # Delete de x_line from the spatial container
            s_container.del_feature(x_line)

        a = list(s_container.get_features())
        print (a)


def manage_cleaning(command, geo_content):
    """Manage the cleaning of the road features

    Parameters

    Return value
    """

    s_container = build_topology(geo_content.in_features, GenUtil.ZERO)

    # Clean junction in Y form
#    if command.yjunction >= 0:
#        clean_y_junction(command, geo_content)

    # Clean junction in X form
    if command.xjunction >= 0:
        clean_x_junction(s_container, command, geo_content)

    # Extend line
#    if command.extend_line >= 0:
#        extend_lines(command, geo_content)

#    if command.remove_noise >= 0.0
#        remove_noises(command, geo_content)

#    centre_lines = extract_lines(command, geo_content)

#    return centre_lines


@dataclass
class Command:
    """Contains the parameters of the command."""


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
    nbr_xjunction: 0
    nbr_yjunction: 0
    nbr_join: 0
    nbr_noise: 0
    in_features: List[object]
    out_features: List[object]
    nbr_in_line_strings: 0
    nbr_out_line_strings: 0
    bounds: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, in_features=[], out_features=[],
                         nbr_xjunction=0, nbr_yjunction=0, nbr_join=0, nbr_noise=0,
                         nbr_in_line_strings=0, nbr_out_line_strings=0, bounds=[])



# Read the command line arguments
#command = read_arguments()

# Extract and load the layers of the input file
#layers = [command.polygon, command.tesselation]
#GenUtil.read_in_file (command.in_file, geo_content, layers)


a = LineString(((0,0),(5,0)))
b = LineString(((5,0),(5,-5)))
c = LineString(((5,0),(10,5)))
d = LineString(((10,10),(10,5)))
e = LineString(((10,5),(15,5)))
geo_content.in_features = [a,b,c,d,e]

command = SpatialContainer()
command.xjunction = 10

manage_cleaning(command, geo_content)




print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of input polygon layer: {}".format(command.polygon))
print ("Name of input tesselation layer: {}".format(command.tesselation))
print ("Nampe of output skeleton layer: {}".format(command.skeleton))
print ("-----")


# Copy the results in the output file
geo_content.layer_names = [command.skeleton]
GenUtil.write_out_file_append (command.in_file, geo_content)

print ("Number of point features written: {}".format(geo_content.out_nbr_points))
print ("Number of line string features written: {}".format(geo_content.out_nbr_line_strings))
print ("Number of polygon written: {}".format(geo_content.out_nbr_polygons))
