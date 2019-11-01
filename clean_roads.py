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
    parser.add_argument("-y", "--yjunction", type=float, default=0.0, help="max tolerance for cleaning Y junction ")
    parser.add_argument("-x", "--xjunction", type=float, help="max tolerance for cleaning X junction")
    parser.add_argument("-j", "--join", type=float, help="max tolerance for joining 2 roads")
    parser.add_argument("-r", "--remove", type=float, help="max tolerance for removing noise")
    parser.add_argument("-il", "--input_layer", type=str, help="name input layer containing the roads")
    parser.add_argument("-ol", "--output_layer", type=str, help="name output layer containing the roads")

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


def clean_noise(s_container, command, geo_content):

    for noise_line in s_container.get_features():
        delete_line = False
        if noise_line.length <= command.noise:
            # Line is below tolerance... possible edition
            if len(noise_line.start_lines) != 0 or len(noise_line.end_lines) != 0:
                # Line is connected. No edition
                pass
                if len(noise_line.start_lines) == 0:
                    linked_lines = noise_line.start_lines
                else:
                    linked_lines = noise_line.end_lines
                if linked_line[0].length >



def clean_y_junction(s_container, command, geo_content):
    """Clean road junction that form forms a configuration in Y """

    # Loop over each line in the container
    for y_line in s_container.get_features():

        # Process both the start and the end of the line
        start_end_lines = [y_line.start_lines] + [y_line.end_lines]
        for i, linked_lines in enumerate(start_end_lines):
            if  len(linked_lines) == 2 and \
                len(linked_lines[0].coords) >= 3 and \
                len(linked_lines[1].coords) >= 3:
                # The current line must be link to 2 and exactly 2 lines
                line_to_edit = []  # List containing a reference to the lines to edit; if the Y crossing is edited
                coord_to_edit = []  # List containing tfor each the coord index to edit
                line_to_edit.append(y_line) # append the current line
                # Extract the position of the P0
                if i == 0:
                    # Testing Y crossing on the first vertice of the line
                    coord_to_edit.append(0)
                    p0 = y_line.coords[0]
                else:
                    # Testing Y crossing on the last vertice of the line
                    coord_to_edit.append(-1)
                    p0 = y_line.coords[-1]

                for num_line, linked_line in enumerate(linked_lines):
                    line_to_edit.append(linked_line)
                    if Point(p0).distance(Point(linked_line.coords[0])) <= GenUtil.ZERO:
                        # Extract the second and third vertices
                        i1 = 1
                        i2 = 2
                        coord_to_edit.append(0)
                    else:
                        # Extract before last and before before last vertice
                        i1 = -2
                        i2 = -3
                        coord_to_edit.append(-1)
                        # P0 is located at the start of the linked line
                    if num_line==0:
                        # Process of the first line
                        p11 = linked_line.coords[i1]
                        p12 = linked_line.coords[i2]
                    else:
                        # Process of the second line
                        p21 = linked_line.coords[i1]
                        p22 = linked_line.coords[i2]

                angle_p11_p0_21 = GenUtil.compute_angle(p11, p0, p21, type=GenUtil.DEGREE)
                angle_p11_p21_p22 = GenUtil.compute_angle(p11, p21, p22, type=GenUtil.DEGREE)
                angle_p12_p11_p21 = GenUtil.compute_angle(p11, p21, p22, type=GenUtil.DEGREE)
                if angle_p11_p21_p22 > 100. and angle_p12_p11_p21 > 100. and \
                   angle_p11_p0_21 < 165.:
                    # The base of the triangle is pseudo alignes
                    point_mid_p11_p21 = LineString((p11, p21)).interpolate(.5, normalized=True)
                    coord_mid_p11_21 = (point_mid_p11_p21.x,point_mid_p11_p21.y)
                    line_y_crossing = LineString((p0, coord_mid_p11_21))
                    if (line_y_crossing.length <= command.yjunction):
                        # The Y crossing requirements are all met ===> edit the line now
                        for i in range(3):
                            # Loop over the 3 line to edit
                            line = line_to_edit[i]
                            ind = coord_to_edit[i]
                            lst_coord = list(line.coords)
                            lst_coord[ind] = coord_mid_p11_21
                            line.coords = lst_coord
                            print ("Adjust update spatial index!!!")
                            #s_container.update_spatial_index(line)
                    else:
                        # Over the tolerance do not edit the line
                        pass
                else:
                    # the base of the triangle is not flat enough do not edit the line
                    pass
            else:
                # the basic requirment are not met do not edit the line
                pass


def adjust_coords (target_point, new_point, line, s_container):

    if target_point.distance(Point(line.coords[0])) <= GenUtil.ZERO:
        # Move/edit the first coordinate to the new mid point
        i = 0
    else:
        # Move/edit the last coordinate to the new mid point
        i = len(line.coords) - 1

    lst_coords = list(line.coords)
    lst_coords[i] = (new_point.x, new_point.y)
    line.coords = lst_coords
    #s_container.update_spatial_index(line)





def clean_x_junction(s_container, command, geo_content):
    """Clean a road junction that should form 4 line crossing"""

    x_lines = list(s_container.get_features())
    # Loop over each line in the container
    for x_line in x_lines:
        # Only select line below the tolerance and linked to 2 other lines at each extremity
        if x_line.length <= command.xjunction and \
           len(x_line.start_lines) == 2 and \
           len(x_line.end_lines) == 2:
            mid_point = x_line.interpolate(.5, normalized=True) # Find the new intersection point

            start_point = Point(x_line.coords[0])
            adjust_coords(start_point, mid_point, x_line.start_lines[0], s_container)
            adjust_coords(start_point, mid_point, x_line.start_lines[1], s_container)
            end_point = Point(x_line.coords[-1])
            adjust_coords(end_point, mid_point, x_line.end_lines[0], s_container)
            adjust_coords(end_point, mid_point, x_line.end_lines[1], s_container)

            # Delete de x_line from the spatial container
            s_container.del_feature(x_line)
            geo_content.nbr_xjunction += 1  # Add stats counter



def manage_cleaning(command, geo_content):
    """Manage the cleaning of the road features

    Parameters

    Return value
    """

    s_container = build_topology(geo_content.in_features, GenUtil.ZERO)

    # Clean the noise
    if commend.noise:
        clean_noise((s_container, command, geo_content))

    # Clean junction in X form
    if command.xjunction >= 0:
        clean_x_junction(s_container, command, geo_content)

    # Clean junction in Y form
#    if command.yjunction >= 0:
#        clean_y_junction(s_container, command, geo_content)



    # Extend line
#    if command.extend_line >= 0:
#        extend_lines(command, geo_content)

#    if command.remove_noise >= 0.0
#        remove_noises(command, geo_content)

#    centre_lines = extract_lines(command, geo_content)

#    return centre_lines

    geo_content.out_features = list(s_container.get_features())

    return


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
    in_nbr_line_strings: 0
    out_nbr_line_strings: 0
    bounds: List[object] = None

geo_content = GeoContent(crs=None, driver=None, schemas={}, in_features=[], out_features=[],
                         nbr_xjunction=0, nbr_yjunction=0, nbr_join=0, nbr_noise=0,
                         in_nbr_line_strings=0, out_nbr_line_strings=0, bounds=[])



# Read the command line arguments
command = read_arguments()


GenUtil.read_in_file (command.in_file, geo_content, [command.input_layer])


a = LineString(((0,0),(5,0)))
b = LineString(((5,0),(5,-5)))
c = LineString(((5,0),(10,5)))
d = LineString(((10,10),(10,5)))
e = LineString(((10,5),(15,5)))

f = LineString(((10,10),(10,20)))
g = LineString (((6,21),(8,21),(10,20)))
h = LineString (((10,20),(12,21),(14,21)))
i = LineString (((10,10), (8,9), (6,9)))
j = LineString (((14,9),(12,9),(10,10)))

#geo_content.in_features = [a,b,c,d,e]

#command = SpatialContainer()
#command.xjunction = 10
#command.yjunction=10

manage_cleaning(command, geo_content)




print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of input layer: {}".format(command.input_layer))
print("Name of output_layer: {}".format(command.output_layer))
print("Tolerance for Y junction: {}".format(command.yjunction))
print("Tolerance for X junction: {}".format(command.xjunction))
print ("-----")


# Copy the results in the output file
geo_content.layer_names = [command.output_layer]
for feature in geo_content.out_features:
    feature.sb_layer_name = command.output_layer

GenUtil.write_out_file_append (command.in_file, geo_content)
