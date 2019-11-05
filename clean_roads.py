#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, os
from dataclasses import dataclass
from typing import List
from lib_geosim import SpatialContainer, GenUtil
from shapely.geometry import Point, LineString
from shapely.ops import linemerge

FIRST = 0
LAST = -1
BOTH = "Both"


class Topology(object)

    def __init__(self, s_container, line, search_tolerance=GenUtil.ZERO):

        self.line = line
        self.s_container = s_container
        self.lst_coord = list(line.coords)
        self.start_linked_lines = []
        self.end_linked_lines = []

        # Find the linked lines of the start of the line
        p_start = self.lst_coord[0]
        b_box = GenUtil.build_bounding_box(search_tolerance, p_start)
        linked_lines = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
        for linked_line in linked_lines:
            if self.line.touches(linked_line):
                self.start_linked_lines.append(linked_line)

        # Find the linked lines of the end of the line
        p_end = self.line.coord[-1]
        b_box = GenUtil.build_bounding_box(search_tolerance, p_end)
        linked_lines = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
        for linked_line in linked_lines:
            if line.touches(linked_line):
                self.end_linked_lines.append(linked_line)

        return


    def orient_linked_lines(self, target=BOTH):

        if target == BOTH:
            for (ind,linked_lines) in ((0,self.start_linked_lines), (self.end_linked_lines, -1)):
                for linked_line in linked_lines:
                    lst_linked_line_coord = list(linked_line.coords)
                    if Point(self.lst_coord[ind]).distance(Point(lst_linked_line_coord[0])) > GenUtil.distance:
                        lst_linked_line_coord.reverse()
                        linked_line.coords = lst_linked_line_coord
        else:
            raise Exception "Only value BOTH is accepted for 'target' parameter"


    def reverse(self):

        self.lst_coord.reverse()
        self.line.coords = self.lst_coord
        self.lst_coord = list(line.coords)
        self.start_linked_lines = []
        self.end_linked_lines = []


def read_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="geopackage vector file to clean")
    parser.add_argument("-y", "--yjunction", type=float, default=0.0, help="max tolerance for cleaning Y junction ")
    parser.add_argument("-x", "--xjunction", type=float, help="max tolerance for cleaning X junction")
    parser.add_argument("-j", "--join", type=float, help="max tolerance for joining 2 roads")
    parser.add_argument("-n", "--noise", type=float, help="max tolerance for removing noise")
    parser.add_argument("-e", "--extend", type=float, help="max tolerance for extending line")
    parser.add_argument("-i", "--iteration", type=int, default=5, help="Number of iteration for extending line")
    parser.add_argument("-il", "--input_layer", type=str, help="name input layer containing the roads")
    parser.add_argument("-ol", "--output_layer", type=str, help="name output layer containing the roads")

    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))


    return command


def build_topology1(center_lines, search_tolerance):


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





def build_topology(s_container, line):
    """Extract the topology for a line and a position (first or last)"""

    line.start_lines=[]
    line.end_lines=[]
    for start_end,coord in enumerate((line.coords[0],line.coords[-1])):
        b_box = GenUtil.build_bounding_box(GenUtil.ZERO, coord)
        lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
        for line_b_box in lines_b_box:
            if Point(coord).distance(Point(line_b_box.coords[0])) <= GenUtil.ZERO:
                # The line is connected to the start of the other line
                ind = 0
            else:
                if Point(coord).distance(Point(line_b_box.coords[-1])) <= GenUtil.ZERO:
                    # The line is connected to the end of the other line
                    ind = -1
                else:
                    break
            if start_end == 0:
                # Append to the start of the line
                line.start_lines.append((line_b_box,ind))
            else:
                # Append to the end of the line
                line.end_lines.append((line_b_box, ind))





#def clean_noise_one_extrimity(line, extremity):


#    if are_open_arm(target_lines)



def clean_noise_lines(s_container, command, geo_content):

     delete_features = []
     lines = list(s_container.get_features())
     for line in lines:
         if id(line) in delete_features:
             print ("Trouve...")

         start_lines = point_topolgy(s_container, line, FIRST)
         end_lines = point_topolgy(s_container, line, LAST)
         if id(line) not in delete_features and \
            len(start_lines) >= 2 and len(end_lines) >= 2:
             for i, target_lines in enumerate([start_lines]+[line.end_lines]):
                 if are_open_arms(target_lines):
                     if i == 0:
                         ind_coord = 0
                     else:
                         ind_coord = -1

                     if len(target_lines) == 2 and \
                         target_lines[0].length <= command.noise and \
                         target_lines[1].length <= command.noise:
                         p = []
                         for target_line in target_lines:
                             if line.distance(Point(target_line.coords[0]))<= GenUtil.ZERO:
                                 p.append(target_line.coords[-1])
                             else:
                                 p.append(target_line.coords[0])
                         p_line = LineString((p[0], p[1]))
                         if p_line.disjoint(line):
                             new_point = p_line.interpolate(.5, normalized=True)
                             lst_coords = list(line.coords)
                             new_coord = [(new_point.x, new_point.y)]
                             if ind_coord == -1:
                                 lst_coords = lst_coords + new_coord
                             else:
                                 lst_coords = new_coord + lst_coords
                             line.coords = lst_coords
                             s_container.update_spatial_index(line)
                             delete_features.append(id(target_lines[0]))
                             delete_features.append(id(target_lines[1]))
                             s_container.del_features(target_lines)
                             command.noise += 1



def is_edition_needed(lst_coord_0, lst_coord_1, lst_coord_2):

    p0 = lst_coord_0[0]
    l1 = lst_coord_0[1]
    l2 = lst_coord_1[1]
    l3 = lst_coord_2[1]

    angle_1 = GenUtil.compute_angle(l1, p0, l2, type=GenUtil.DEGREE)
    angle_2 = GenUtil.compute_angle(l1, p0, l3, type=GenUtil.DEGREE)
    angle_3 = GenUtil.compute_angle(l2, p0, l3, type=GenUtil.DEGREE)

    if angle_1 >= 165. or angle_2 >= 165. or angle_3 >= 165.:
        edition_needed = False
        print ("Position {}".format(str(p0)))
    else:
        edition_needed = True

    return edition_needed


def metric_y_junction (lst_coord_0, lst_coord_1, command):

    if len(lst_coord_0) >=3 and len(lst_coord_1) >=3:
        p0 = lst_coord_0[0]
        p01 = lst_coord_0[1]
        p02 = lst_coord_0[2]

        p11 = lst_coord_1[1]
        p12 = lst_coord_1[2]

        tmp_mid_p01_p11 = LineString((p01,p11)).interpolate(.5, normalized=True)
        mid_p01_p11 = (tmp_mid_p01_p11.x, tmp_mid_p01_p11.y)
        line_mid_p0 = LineString((mid_p01_p11, p0))
        length_y_junction = line_mid_p0.length

        angle_mid_p01_p02 = GenUtil.compute_angle(mid_p01_p11, p01, p02, type=GenUtil.DEGREE)
        angle_mid_p11_p12 = GenUtil.compute_angle(mid_p01_p11, p11, p12, type=GenUtil.DEGREE)

        if (length_y_junction <= command.yjunction):
            sum_angles = angle_mid_p01_p02 + angle_mid_p11_p12
        else:
            sum_angles = None

    else:
        print ("Add a pseudo mid point...")
        sum_angles = None
        mid_p01_p11 = None

    return (sum_angles, mid_p01_p11)


def clean_y_junction_dummy(s_container, command, geo_content):
    """Clean road junction that form forms a configuration in Y """

    processed_y_junction = []  # create an empty set
    # Loop over each line in the container
    for y0_line in list(s_container.get_features()):

        build_topology(s_container, y0_line)
        y0_lst_coord = list(y0_line.coords)

        # Process both the start and the end of the line
        start_end_lines = [y0_line.start_lines] + [y0_line.end_lines]
        for num, linked_lines in enumerate(start_end_lines):
            if len(linked_lines) == 2:

                y1_line = linked_lines[0][0]
                y1_lst_coord = list(y1_line.coords)
                y1_order = linked_lines[0][1]

                y2_line = linked_lines[1][0]
                y2_lst_coord = list(y2_line.coords)
                y2_order = linked_lines[1][1]

                current_y_junction = sorted((id(y0_line), id(y1_line), id(y2_line))) #
                if current_y_junction not in processed_y_junction:

                    # Add the current junction
                    processed_y_junction.append(current_y_junction)

                    if num == 1:
                        y0_lst_coord.reverse()

                    if y1_order == -1: # Line is connect by the last coordinate
                        y1_lst_coord.reverse()

                    if y2_order == -1: # Line is connect by the last coordinate
                        y2_lst_coord.reverse()

                    if is_edition_needed (y0_lst_coord, y1_lst_coord, y2_lst_coord):

                        metric_0 = metric_y_junction(y0_lst_coord, y1_lst_coord, command)
                        metric_1 = metric_y_junction(y1_lst_coord, y2_lst_coord, command)
                        metric_2 = metric_y_junction(y0_lst_coord, y2_lst_coord, command)

                        max_sum_angles = -1.
                        new_mid_coord = None
                        for (sum_angles, mid_coord) in (metric_0, metric_1, metric_2):
                            if sum_angles is not None:
                                if sum_angles >= max_sum_angles:
                                    max_sum_angles = sum_angles
                                    new_mid_coord = mid_coord

                        if new_mid_coord is not None:

                            # The Y crossing requirements are all met ===> edit the line now
                            # Update the lines forming the Y crossing
                            y0_lst_coord[0] = new_mid_coord
                            y0_line.coords = y0_lst_coord
                            s_container.update_spatial_index(y0_line)

                            y1_lst_coord[0] = new_mid_coord
                            y1_line.coords = y1_lst_coord
                            s_container.update_spatial_index(y1_line)

                            y2_lst_coord[0] = new_mid_coord
                            y2_line.coords = y2_lst_coord
                            s_container.update_spatial_index(y2_line)
                            geo_content.nbr_yjunction += 1

                        else:
                            # No line met the basic requirements
                            pass
                    else:
                        # No edition is needed this is not a Y crossing
                        pass
                else:
                    # Junction already processed
                    pass
            else:
                # Junction must be formes by exactly 3 lines
                pass



def clean_y_junction(s_container, command, geo_content):
    """Clean road junction that form forms a configuration in Y """

    processed_y_junction = []  # create an empty set
    # Loop over each line in the container
    for y0_line in list(s_container.get_features()):

        y0_topology = Topology(s_container, y0_line)
        y0_lst_coord = list(y0_line.coords)

        # Process both the start and the end of the line
        start_end_lines = [(FIRST, y0_topology.start_lines)] + [(LAST, y0_topology.end_lines)]
        for (ind, linked_lines) in enumerate(start_end_lines):
            if len(linked_lines) == 2:
                # Line needs to be linked to 2 and only 2 lines to form a valid Y junction
                y0_topology.orient_linked_lines()
                y1_line = linked_lines[0]
                y1_lst_coord = list(y1_line.coords)

                y2_line = linked_lines[1]
                y2_lst_coord = list(y2_line.coords)

                current_y_junction = sorted((id(y0_line), id(y1_line), id(y2_line)))
                if current_y_junction not in processed_y_junction:

                    # Add the current junction to the list of process jnction
                    processed_y_junction.append(current_y_junction)

                    if is_edition_needed (y0_lst_coord, y1_lst_coord, y2_lst_coord):

                        metric_0 = metric_y_junction(y0_lst_coord, y1_lst_coord, command)
                        metric_1 = metric_y_junction(y1_lst_coord, y2_lst_coord, command)
                        metric_2 = metric_y_junction(y0_lst_coord, y2_lst_coord, command)

                        max_sum_angles = -1.
                        new_mid_coord = None
                        for (sum_angles, mid_coord) in (metric_0, metric_1, metric_2):
                            if sum_angles is not None:
                                if sum_angles >= max_sum_angles:
                                    max_sum_angles = sum_angles
                                    new_mid_coord = mid_coord

                        if new_mid_coord is not None:

                            # The Y crossing requirements are all met ===> edit the line now
                            # Update the line formiog the Y crossing
                            y0_lst_coord[0] = new_mid_coord
                            y0_line.coords = y0_lst_coord
                            s_container.update_spatial_index(y0_line)

                            y1_lst_coord[0] = new_mid_coord
                            y1_line.coords = y1_lst_coord
                            s_container.update_spatial_index(y1_line)

                            y2_lst_coord[0] = new_mid_coord
                            y2_line.coords = y2_lst_coord
                            s_container.update_spatial_index(y2_line)
                            geo_content.nbr_yjunction += 1

                        else:
                            # No line met the basic requirements
                            pass
                    else:
                        # No edition is needed this is not a Y crossing
                        pass
                else:
                    # Junction already processed
                    pass
            else:
                # Junction must be formes by exactly 3 lines
                pass


def is_open_arm(s_container, line):

    open_arm = False
    build_topology(s_container, line)
    if len(line.start_lines) == 0 or len(line.end_lines) == 0:
        open_arm = True

    return open_arm

def clean_dual_fork(s_container, command, geo_content):


    fork_lines = list(s_container.get_features())
    # Loop over each line in the list

    processed_lines = []
    for fork_line in fork_lines:

        if id(fork_line) not in processed_lines:

            # Build the topology for the line
            build_topology(s_container, fork_line)

            start_end_lines = [fork_line.start_lines]+[fork_line.end_lines]
            for num, linked_lines in enumerate(start_end_lines):
                if len(linked_lines) == 2:
                    line0 = linked_lines[0][0]
                    line1 = linked_lines[1][0]
                    if line0.length <= command.noise and \
                       line1.length <= command.noise and \
                       is_open_arm(s_container, line0) and \
                       is_open_arm(s_container, line1):
                        lst_coord0 = list(line0.coords)
                        if linked_lines[0][1] == -1:
                            lst_coord0.reverse()
                        lst_coord1 = list(line1.coords)
                        if linked_lines[1][1] == -1:
                            lst_coord1.reverse()

                        line = LineString((lst_coord0[-1], lst_coord1[-1]))
                        mid_point = line.interpolate(.5, normalized=True)
                        mid_coord = (mid_point.x, mid_point.y)
                        fork_line_coord = list(fork_line.coords)
                        if num == 0:
                            fork_line_coord = [mid_coord] + fork_line_coord
                        else:
                            fork_line_coord = fork_line_coord + [mid_coord]

                        fork_line.coords = fork_line_coord
                        s_container.del_features ([line0, line1])
                        s_container.update_spatial_index(fork_line)
                        processed_lines.append(id(line0))
                        processed_lines.append(id(line1))
                        geo_content.nbr_noise += 2

                    else:
                        # Line does not pass the requirements
                        pass
                else:
                    # Line does not pass the requirements
                    pass
        else:
            # Line already processed
            pass


def clean_x_junction(s_container, command, geo_content):
    """Clean a road junction that should form 4 line crossing"""

    x_lines = list(s_container.get_features())
    # Loop over each line in the container
    for x_line in x_lines:
        # Only select line below the tolerance and linked to 2 other lines at each extremity
        if x_line.length <= command.xjunction:
            # Build the topology for the line
            build_topology(s_container, x_line)
            if len(x_line.start_lines) == 2 and \
               len(x_line.end_lines) == 2:
                mid_point = x_line.interpolate(.5, normalized=True) # Find the new intersection point
                coord_point = [mid_point.x, mid_point.y]
                for tuple_line_ind in (x_line.start_lines+x_line.end_lines):
                    line = tuple_line_ind[0]
                    ind = tuple_line_ind[1]
                    lst_coord = list(line.coords)
                    lst_coord[ind] = coord_point
                    line.coords = lst_coord
                    #s_container.update_spatial_index(line)
                    print ("Edit Update spatial index...")

                # Delete de x_line from the spatial container
                s_container.del_feature(x_line)
                geo_content.nbr_xjunction += 1  # Add stats counter


def extend_line(s_container, command, geo_content):

    lines = list(s_container.get_features())

    # Loop over each line in the list
    processed_lines = []
    for line in lines:

        if id(line) not in processed_lines:

            # Loop over the first and last vertice of the line
            for ind in [0, -1]:
                lst_coord_line = list(line.coords)
                coord_line = lst_coord_line[ind]
                build_topology(s_container, line)
                if (ind == 0 and len(line.start_lines)) == 0 or \
                   (ind == -1 and len(line.end_lines)) == 0 :
                    #The line is open on the requested side
                    # Loop to with increasing tolerance to find target line to merge
                    for i in range(command.iteration):
                        tolerance = command.extend * (float(i+1)/command.iteration)
                        b_box = GenUtil.build_bounding_box(tolerance, coord_line)
                        lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
                        if lines_b_box:
                            target_line = None
                            for line_b_box in lines_b_box:
                                lst_coord_target = list(line_b_box.coords)
                                coord_tl_0 = lst_coord_target[0]
                                coord_tl_n = lst_coord_target[-1]
                                if Point(coord_line).distance(Point(coord_tl_0)) <= tolerance:
                                    target_line = line_b_box
                                    coord_target = coord_tl_0
                                elif Point(coord_line).distance(Point(coord_tl_n)) <= tolerance:
                                    target_line = line_b_box
                                    coord_target = coord_tl_n

                            if target_line:
                                #Merge the lines
                                new_line = LineString((coord_line, coord_target))
                                new_line.sb_layer_name = line.sb_layer_name
                                new_line.sb_properties = line.sb_properties
                                s_container.add_feature(new_line)
                                break

#                                merged_line = linemerge([line, new_line, target_line])
#                                if merged_line.geom_type == GenUtil.LINE_STRING:
#                                    lst_merged_coord = list(merged_line.coords)
#                                    if Point(coord_line).distance(Point(lst_merged_coord[ind])) >= GenUtil.ZERO:
#                                        # Line was swapped during the merge process. Reverse it
#                                        lst_merged_coord.reverse()
#                                    line.coords = lst_merged_coord
#                                    s_container.del_feature(target_line)
#                                    s_container.update_spatial_index(line)
#                                    processed_lines.append(id(target_line))
#                                    # Merged done...
#                                    break
#                                else:
#                                    # Something weird happen... pass to the next line
#                                    pass
#                                    break
                            else:
                                # Nothing found. Increase tolerance
                                pass

                else:
                    # Line has no open segment (both sides have lines)
                    pass

        else:
            # Line already processed
            pass



def manage_cleaning(command, geo_content):
    """Manage the cleaning of the road features

    Parameters

    Return value
    """

    # Load the features in the spatial container
    s_container = SpatialContainer()
    s_container.add_features(geo_content.in_features)
    geo_content.in_nbr_line_strings = len(geo_content.in_features)
    geo_content.in_features=[]

    # Clean the noise
    if command.noise:
        clean_dual_fork(s_container, command, geo_content)
        clean_dual_fork(s_container, command, geo_content)


    # Clean junction in X form
    if command.xjunction >= 0:
        clean_x_junction(s_container, command, geo_content)

    # Clean junction in Y form
    if command.yjunction >= 0:
        clean_y_junction(s_container, command, geo_content)

    # Extend line
    if command.extend >= 0:
        extend_line(s_container, command, geo_content)

    geo_content.out_features = list(s_container.get_features())
    geo_content.out_nbr_line_strings = len(geo_content.out_features)

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
    nbr_noise: 0
    nbr_join: 0
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


# test for Ycrossing
a = LineString(((0,0),(5,0)))
b = LineString(((5,0),(5,-5)))
c = LineString(((5,0),(10,5)))
d = LineString(((10,10),(10,5)))
e = LineString(((10,5),(15,5)))

# Test for X crossing
f = LineString(((10,10), (10,15),(10,20)))
g = LineString (((6,21),(8,21),(10,20)))
h = LineString (((10,20),(12,21),(14,21)))
i = LineString (((10,10), (8,9), (6,9)))
j = LineString (((14,9),(12,9),(10,10)))

# Test for dual fork
k = LineString (((10,0), (30,0)))
l = LineString (((30,0),(35,3)))
m = LineString (((30,0),(35,-3)))
n = LineString (((5,3),(10,0)))
o = LineString (((10,0),(5,-3)))

# Test for extend

p = LineString(((20,0),(30,0)))
q = LineString(((32,1),(35,0)))
r = LineString(((19,1),(10,0)))
s = LineString(((10,0),(0,0)))
t = LineString(((10,0),(10,10)))


#geo_content.in_features = [p,q,r]

#command = SpatialContainer()
#command.extend = 3
#command.iteration=5
#command.noise = 7

manage_cleaning(command, geo_content)




print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of input layer: {}".format(command.input_layer))
print("Name of output_layer: {}".format(command.output_layer))
print("Tolerance for Y junction: {}".format(command.yjunction))
print("Tolerance for X junction: {}".format(command.xjunction))
print("Tolerance for noise detection: {}".format(command.xjunction))
print("Tolerance for extending lines: {}".format(command.extending))
print ("-----")
print("Number of features (line) read: {}".format(geo_content.in_nbr_line_strings.extending))
print("Number of features (line) written: {}".forma(geo_content.out_nbr_line_strings.extending))
print("Number of Y junction corrected: {}".format(geo_content.nbr_y_junction))
print("Number of X junction corrected: {}".format(geo_content.nbr_x_junction))
print("Number of line (noise) deleted: {}".format(geo_content.nbr_noise))
print("Number of line extended: {}".format(geo_content.nbr_extended))


# Copy the results in the output file
geo_content.layer_names = [command.output_layer]
for feature in geo_content.out_features:
    feature.sb_layer_name = command.output_layer

GenUtil.write_out_file_append (command.in_file, geo_content)
