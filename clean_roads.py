#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse, os
from dataclasses import dataclass
from typing import List
from lib_geosim import SpatialContainer, GenUtil, LineStringSc
from shapely.geometry import Point, LineString
from shapely.ops import linemerge

FIRST = 0
LAST = -1
BOTH = "Both"


class Topology(object):

    def __init__(self, s_container, line, search_tolerance=GenUtil.ZERO):

        self.line = line
        self.search_tolerance = search_tolerance
        self.s_container = s_container
        self.start_linked_lines = []
        self.end_linked_lines = []

        # Find the linked lines of the start of the line
        self.start_linked_lines = self._extract_linked_lines(self.line.coords[0])

        # Find the linked lines of the end of the line
        self.end_linked_lines = self._extract_linked_lines(self.line.coords[-1])

        return

    def _extract_linked_lines(self, coord):

        linked_lines = []
        b_box = GenUtil.build_bounding_box(self. search_tolerance, coord)
        potential_lines = self.s_container.get_features(b_box, remove_features=[self.line])
        for potential_line in potential_lines:
            if potential_line.distance(Point(coord)) <= GenUtil.ZERO:
                linked_lines.append(potential_line)

        return linked_lines

    def orient_linked_lines(self, target=BOTH):

        line_lst_coord = list(self.line.coords)
        if target == BOTH:
            for (ind, linked_lines) in ((0, self.start_linked_lines), (-1, self.end_linked_lines)):
                for linked_line in linked_lines:
                    lst_linked_line_coord = list(linked_line.coords)
                    if Point(line_lst_coord[ind]).distance(Point(lst_linked_line_coord[0])) > GenUtil.ZERO:
                        lst_linked_line_coord.reverse()
                        linked_line.coords = lst_linked_line_coord
        else:
            raise GeoSimException("Only value BOTH is accepted for 'target' parameter")

    def reverse(self):

        lst_coord = list(self.line.coords)
        lst_coord.reverse()
        self.line.coords = lst_coord
        tmp = self.start_linked_lines
        self.start_linked_lines = self.end_linked_lines
        self.end_linked_lines = tmp


def read_arguments():
    """Read and manage the input arguments in the command line"""

    # Setting the parameters of the command line
    parser = argparse.ArgumentParser()
    parser.add_argument("in_file", help="geopackage vector file to clean")
    parser.add_argument("-y", "--yjunction", type=float, default=0.0, help="max tolerance for cleaning Y junction ")
    parser.add_argument("-x", "--xjunction", type=float, default=0.0, help="max tolerance for cleaning X junction")
    parser.add_argument("-j", "--join", type=float, default=0.0, help="max tolerance for joining 2 lines extremity")
    parser.add_argument("-n", "--noise", type=float, default=0.0, help="max tolerance for removing noise")
    parser.add_argument("-e", "--extend", type=float, default=0.0, help="max tolerance for extending a line to an extremity")
    parser.add_argument("-o", "--orphan", type=float, default=0.0, help="max tolerance for deleting orphan line")
    parser.add_argument("-i", "--iteration", type=int, default=5,  help="Number of iteration for extending line")
    parser.add_argument("-il", "--input_layer", type=str, help="name input layer containing the roads")
    parser.add_argument("-ol", "--output_layer", type=str, help="name output layer containing the roads")

    # Read the command line parameter
    command = parser.parse_args()

    # Check that the triangle input file exist. Exit if missing
    if not os.path.isfile(command.in_file):
        raise Exception('Input file is missing: {}'.format(command.in_file))

    return command

# def build_topology1(center_lines, search_tolerance):
#
#
#     # Load the features in the spatial container (accelerate the search)
#     s_container = SpatialContainer()
#     for center_line in center_lines:
#         center_line.sb_geom_type = GenUtil.LINE_STRING
#         s_container.add_feature(center_line)
#
#     # Build topology. For each create list of connecting lines
#     for line in center_lines:
#         p_start = line.coords[0]
#         p_end = line.coords[-1]
#         b_box = GenUtil.build_bounding_box(search_tolerance, p_start)
#         lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
#         line.start_lines = []
#         for line_b_box in lines_b_box:
#             if line.touches(line_b_box):
#                 line.start_lines.append(line_b_box)
#         b_box = GenUtil.build_bounding_box(search_tolerance, p_end)
#         lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
#         line.end_lines = []
#         for line_b_box in lines_b_box:
#             if line.touches(line_b_box):
#                 line.end_lines.append(line_b_box)
#
#     return s_container
#
# def build_topology(s_container, line):
#     """Extract the topology for a line and a position (first or last)"""
#
#     line.start_lines=[]
#     line.end_lines=[]
#     for start_end,coord in enumerate((line.coords[0], line.coords[-1])):
#         b_box = GenUtil.build_bounding_box(GenUtil.ZERO, coord)
#         lines_b_box = s_container.get_features(b_box, remove_features=[line._sb_sc_id])
#         for line_b_box in lines_b_box:
#             if Point(coord).distance(Point(line_b_box.coords[0])) <= GenUtil.ZERO:
#                 # The line is connected to the start of the other line
#                 ind = 0
#             else:
#                 if Point(coord).distance(Point(line_b_box.coords[-1])) <= GenUtil.ZERO:
#                     # The line is connected to the end of the other line
#                     ind = -1
#                 else:
#                     break
#             if start_end == 0:
#                 # Append to the start of the line
#                 line.start_lines.append((line_b_box,ind))
#             else:
#                 # Append to the end of the line
#                 line.end_lines.append((line_b_box, ind))

def clean_noise_type_2(s_container, command, geo_content):

    delete_features = []
    lines = list(s_container.get_features())
    for line in lines:

        if line.length <= command.noise:

            if id(line) not in delete_features:

                line_topology = Topology(s_container, line)
                line_topology.orient_linked_lines()

                for ind in [FIRST, LAST]:
                    if ind == LAST:
                        line_topology.reverse()

                    if len(line_topology.start_linked_lines) == 2 and len(line_topology.end_linked_lines) == 0:

                        line_a = line_topology.start_linked_lines[0]
                        line_b = line_topology.start_linked_lines[1]

                        line_a_topology = Topology(s_container, line_a)
                        line_b_topology = Topology(s_container, line_b)

                        if (len(line_a_topology.end_linked_lines) >= 2 or line_a.length > command.noise) and \
                           (len(line_b_topology.end_linked_lines) >= 2 or line_b.length > command.noise):
                            # the requirements are met to delete the line
                            merged_line = linemerge((line_a, line_b))
                            if merged_line.geom_type == GenUtil.LINE_STRING:
                                # Merge the line and house keeping in the spatial container
                                lst_coord = list(merged_line.coords)
                                line_a.coords = lst_coord
                                delete_features.append(id(line_b))
                                delete_features.append(id(line))
                                s_container.del_feature(line)
                                s_container.del_feature(line_b)
                                s_container.update_spatial_index(line_a)
                                break
                            else:
                                # Merging did not worked as planned... skip...
                                    pass
                        else:
                            # requirement not met for delete line
                            pass
                    else:
                        # The line is not open... skip...
                        pass
            else:
                # Feature already processed
                pass
        else:
            # Features over noise tolerace
            pass


def clean_noise(s_container, command, geo_content):

    lines = list(s_container.get_features())
    for line in lines:
        line_topology = Topology(s_container, line)
        line_topology.orient_linked_lines()

        if line.length > command.noise or len(line_topology.start_linked_lines) >= 2:
            for ind in [FIRST, LAST]:
                if ind == LAST:
                    line_topology.reverse()
                linked_lines = line_topology.start_linked_lines
                if len(linked_lines) == 2 and \
                   is_open_arm(s_container, linked_lines[0]) and \
                   is_open_arm(s_container, linked_lines[1]):

                    if linked_lines[0].length <= command.noise and linked_lines[1].length <= command.noise:
                        p = []
                        for linked_line in linked_lines:
                            p.append(linked_line.coords[-1])
                        p_line = LineString((p[0], p[1]))
                        if p_line.disjoint(line):
                            new_point = p_line.interpolate(.5, normalized=True)
                            lst_coords = list(line.coords)
                            new_coord = [(new_point.x, new_point.y)]
                            lst_coords = new_coord + lst_coords
                            line.coords = lst_coords
                            
                            geo_content.nbr_noise += 2
                        else:
                            # The new line is touching something ===> No edtion
                            pass
                    else:
                        # There must be 2 lines under the noise tolerance
                        pass
                else:
                    # The line to edit must have an open arm (no links)
                    pass
        else:
            # The line is not linked to other line
            pass


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
    else:
        edition_needed = True

    return edition_needed


def calculate_mid_point(lst_coord_0, lst_coord_1, lst_coord_2):

    lst_coord = [[[None], [None], [None]], [[None], [None], [None]], [[None], [None], [None]]]
    for coords in ([lst_coord_0, lst_coord_1, lst_coord_2]):
        if len(coords) <= 2:
            lst_coord[i][0] = coords[0]
            mid_point = LineString((coords[0], coords[1])).interpolate(.5, normalized=True)
            lst_coord[i][1] = (mid_point.x, mid_point.y)
            lst_coord[i][2] = coords[1]
        else:
            for j in range(3):
                lst_coord[i][j] = coords[j]

    try_0 = (lst_coord[0], lst_coord[1])
    try_1 = (lst_coord[1], lst_coord[2])
    try_2 = (lst_coord[0], lst_coord[2])
    max_sum_angles = -1.
    max_mid_p01_p11 = None

    for (lst_coord_a, lst_coord_b) in (try_0, try_1, try_2 ):

        p0 = lst_coord_a[0]
        p01 = lst_coord_a[1]
        p02 = lst_coord_a[2]

        p11 = lst_coord_b[1]
        p12 = lst_coord_b[2]

        tmp_mid_p01_p11 = LineString((p01, p11)).interpolate(.5, normalized=True)
        mid_p01_p11 = (tmp_mid_p01_p11.x, tmp_mid_p01_p11.y)
        line_mid_p0 = LineString((mid_p01_p11, p0))
        length_y_junction = line_mid_p0.length

        angle_mid_p01_p02 = GenUtil.compute_angle(mid_p01_p11, p01, p02, type=GenUtil.DEGREE)
        angle_mid_p11_p12 = GenUtil.compute_angle(mid_p01_p11, p11, p12, type=GenUtil.DEGREE)

        if length_y_junction <= command.yjunction:
            sum_angles = angle_mid_p01_p02 + angle_mid_p11_p12
            if sum_angles > max_sum_angles:
                max_sum_angles = sum_angles
                max_mid_p01_p11 = mid_p01_p11

    return max_mid_p01_p11

# def clean_y_junction_dummy(s_container, command, geo_content):
#     """Clean road junction that form forms a configuration in Y """
#
#     processed_y_junction = []  # create an empty set
#     # Loop over each line in the container
#     for y0_line in list(s_container.get_features()):
#
#         build_topology(s_container, y0_line)
#         y0_lst_coord = list(y0_line.coords)
#
#         # Process both the start and the end of the line
#         start_end_lines = [y0_line.start_lines] + [y0_line.end_lines]
#         for num, linked_lines in enumerate(start_end_lines):
#             if len(linked_lines) == 2:
#
#                 y1_line = linked_lines[0][0]
#                 y1_lst_coord = list(y1_line.coords)
#                 y1_order = linked_lines[0][1]
#
#                 y2_line = linked_lines[1][0]
#                 y2_lst_coord = list(y2_line.coords)
#                 y2_order = linked_lines[1][1]
#
#                 current_y_junction = sorted((id(y0_line), id(y1_line), id(y2_line))) #
#                 if current_y_junction not in processed_y_junction:
#
#                     processed_y_junction.append(current_y_junction)
#
#                     if is_edition_needed (y0_lst_coord, y1_lst_coord, y2_lst_coord):
#
#                         new_mid_coord = calculate_mid_point(y0_lst_coord, y1_lst_coord, y1_lst_coord, command)
#
#                         if new_mid_coord is not None:
#
#                             # The Y crossing requirements are all met ===> edit the line now
#                             # Update the lines forming the Y crossing
#                             y0_lst_coord[0] = new_mid_coord
#                             y0_line.coords = y0_lst_coord
#                             s_container.update_spatial_index(y0_line)
#
#                             y1_lst_coord[0] = new_mid_coord
#                             y1_line.coords = y1_lst_coord
#                             s_container.update_spatial_index(y1_line)
#
#                             y2_lst_coord[0] = new_mid_coord
#                             y2_line.coords = y2_lst_coord
#                             s_container.update_spatial_index(y2_line)
#                             geo_content.nbr_yjunction += 1
#
#                         else:
#                             # No line met the basic requirements
#                             pass
#                     else:
#                         # No edition is needed this is not a Y crossing
#                         pass
#                 else:
#                     # Junction already processed
#                     pass
#             else:
#                 # Junction must be formes by exactly 3 lines
#                 pass
#

def clean_y_junction(s_container, geo_content):
    """Clean road junction that form forms a configuration in Y """

    processed_y_junction = []  # create an empty set
    # Loop over each line in the container
    for y0_line in s_container.get_features():

        y0_topology = Topology(s_container, y0_line)
        y0_topology.orient_linked_lines()

        # Process both the start and the end of the line
        for ind in [FIRST,LAST]:
            if ind == LAST:
                y0_topology.reverse()
            y0_lst_coord = list(y0_line.coords)
            linked_lines = y0_topology.start_linked_lines
            if len(linked_lines) == 2:
                # Line needs to be linked to 2 and only 2 lines to form a valid Y junction
                y1_line = linked_lines[0]
                y1_lst_coord = list(y1_line.coords)

                y2_line = linked_lines[1]
                y2_lst_coord = list(y2_line.coords)

                current_y_junction = sorted((id(y0_line), id(y1_line), id(y2_line)))
                if current_y_junction not in processed_y_junction:

                    # Add the current junction to the list of process jnction
                    processed_y_junction.append(current_y_junction)

                    if is_edition_needed(y0_lst_coord, y1_lst_coord, y2_lst_coord):

                        new_mid_coord = calculate_mid_point(y0_lst_coord, y1_lst_coord, y2_lst_coord)

                        if new_mid_coord is not None:

                            # The Y crossing requirements are all met ===> edit the line now
                            # Update the line forming the Y crossing
                            y0_lst_coord[0] = new_mid_coord
                            y0_line.coords = y0_lst_coord

                            y1_lst_coord[0] = new_mid_coord
                            y1_line.coords = y1_lst_coord

                            y2_lst_coord[0] = new_mid_coord
                            y2_line.coords = y2_lst_coord
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
    line_topology = Topology(s_container, line)
    if len(line_topology.start_linked_lines) == 0 or len(line_topology.end_linked_lines) == 0:
        open_arm = True

    return open_arm

# def clean_dual_fork(s_container, command, geo_content):
#
#     fork_lines = list(s_container.get_features())
#     # Loop over each line in the list
#
#     processed_lines = []
#     for fork_line in fork_lines:
#
#         if id(fork_line) not in processed_lines:
#
#             # Build the topology for the line
#             build_topology(s_container, fork_line)
#
#             start_end_lines = [fork_line.start_lines]+[fork_line.end_lines]
#             for num, linked_lines in enumerate(start_end_lines):
#                 if len(linked_lines) == 2:
#                     line0 = linked_lines[0][0]
#                     line1 = linked_lines[1][0]
#                     if line0.length <= command.noise and \
#                        line1.length <= command.noise and \
#                        is_open_arm(s_container, line0) and \
#                        is_open_arm(s_container, line1):
#                         lst_coord0 = list(line0.coords)
#                         if linked_lines[0][1] == -1:
#                             lst_coord0.reverse()
#                         lst_coord1 = list(line1.coords)
#                         if linked_lines[1][1] == -1:
#                             lst_coord1.reverse()
#
#                         line = LineString((lst_coord0[-1], lst_coord1[-1]))
#                         mid_point = line.interpolate(.5, normalized=True)
#                         mid_coord = (mid_point.x, mid_point.y)
#                         fork_line_coord = list(fork_line.coords)
#                         if num == 0:
#                             fork_line_coord = [mid_coord] + fork_line_coord
#                         else:
#                             fork_line_coord = fork_line_coord + [mid_coord]
#
#                         fork_line.coords = fork_line_coord
#                         s_container.del_features ([line0, line1])
#                         s_container.update_spatial_index(fork_line)
#                         processed_lines.append(id(line0))
#                         processed_lines.append(id(line1))
#                         geo_content.nbr_noise += 2
#
#                     else:
#                         # Line does not pass the requirements
#                         pass
#                 else:
#                     # Line does not pass the requirements
#                     pass
#         else:
#             # Line already processed
#             pass


def clean_x_junction(s_container, command, geo_content):
    """Clean a road junction that should form 4 line crossing"""

    # Loop over each line in the container
    for x_line in s_container.get_features():
        # Only select line below the tolerance and linked to 2 other lines at each extremity
        if x_line.length <= command.xjunction:
            # Build the topology for the line
            x_line_topology = Topology(s_container, x_line)
            if len(x_line_topology.start_linked_lines) == 2 and \
               len(x_line_topology.end_linked_lines) == 2:
                x_line_topology.orient_linked_lines()
                mid_point = x_line.interpolate(.5, normalized=True) # Find the new intersection point
                coord_point = [mid_point.x, mid_point.y]
                for linked_line in (x_line_topology.start_linked_lines+x_line_topology.end_linked_lines):
                    lst_coord = list(linked_line.coords)
                    lst_coord[0] = coord_point
                    linked_line.coords = lst_coord

                # Delete de x_line from the spatial container
                s_container.del_feature(x_line)
                geo_content.nbr_xjunction += 1  # Add stats counter


def join_lines(s_container, command, geo_content, tolerance):


    # Loop over each line in the list
    for line in s_container.get_features():
        # Loop over the first and last vertice of the line
        for ind in [FIRST, LAST]:
            line_topology = Topology(s_container, line)
            if ind == LAST:
                line_topology.reverse()
            lst_coord_line = list(line.coords)
            if len(line_topology.start_linked_lines) == 0:
                b_box = GenUtil.build_bounding_box(tolerance, lst_coord_line[0])
                potential_lines = s_container.get_features(b_box, remove_features=[line])
                target_line = None
                for potential_line in potential_lines:
                    for i in [FIRST, LAST]:
                        lst_coord_potential = list(potential_line.coords)
                        if Point(lst_coord_line[0]).distance(Point(lst_coord_potential[i])) <= tolerance:
                            target_line = potential_line
                            target_coord = lst_coord_potential[i]

                    if target_line:
                        # We have a line to join with
                        new_line = LineString((lst_coord_line[0], target_coord))
                        merged_line = linemerge([line, new_line, target_line])
                        if merged_line.geom_type == GenUtil.LINE_STRING:
                            lst_merged_line_coord = list(merged_line.coords)
                            line.coords = lst_merged_line_coord
                            s_container.del_feature(target_line)
                            geo_content.nbr_join += 1
                        else:
                            # Possible problem with the merged line go to next line
                            pass
                    else:
                        # No line to merged with
                        pass
            else:
                # Line is not open. Go to next line
                pass

    return

def extend_line(s_container, command, geo_content, tolerance):

    # Loop over each line in the list
    for line in s_container.get_features():
        # Loop over the first and last vertice of the line
        for ind in [FIRST, LAST]:
            line_topology = Topology(s_container, line)
            if ind == LAST:
                line_topology.reverse()
            lst_coord_line = list(line.coords)
            if len(line_topology.start_linked_lines) == 0:
                b_box = GenUtil.build_bounding_box(tolerance, lst_coord_line[0])
                potential_lines = s_container.get_features(b_box, remove_features=[line])
                new_coord = None
                distance_min = tolerance + 1.
                for potential_line in potential_lines:
                    linear_ref = potential_line.project(Point(lst_coord_line[0]))
                    point_on_line = potential_line.interpolate(linear_ref)
                    distance = Point(lst_coord_line[0]).distance(point_on_line)
                    if distance <= tolerance:
                        if distance_min > distance:
                            distance_min = distance
                            new_coord = (point_on_line.x, point_on_line.y)

                    if new_coord:
                        # We have a line to join with
                        lst_new_coord_line = [new_coord] + lst_coord_line
                        line.coords = lst_new_coord_line
                        geo_content.nbr_extend += 1
                    else:
                        # No line to merged with
                        pass
            else:
                # Line is not open. Go to next line
                pass

    return


def manage_cleaning(command, geo_content):
    """Manage the cleaning of the road features

    Parameters

    Return value
    """

    # Load the features in the spatial container
    s_container = SpatialContainer()
    s_container.add_features(geo_content.in_features)
    geo_content.in_nbr_line_strings = len(geo_content.in_features)
    geo_content.in_features = []

    # Clean the noise type 2
#    if command.noise:
#        clean_noise_type_2(s_container, command, geo_content)

    # Clean the noise type 1
    if command.noise:
        clean_noise(s_container, command, geo_content)
        clean_noise(s_container, command, geo_content)

    # Clean junction in X form
#    if command.xjunction >= 0:
#        clean_x_junction(s_container, command, geo_content)

    # Clean junction in Y form
#    if command.yjunction >= 0:
#        clean_y_junction(s_container, geo_content)

    # Join line
#    if command.join >= 0:
#        for i in range(command.iteration):
#            tolerance = command.join * (float(i + 1) / command.iteration)
#            join_lines(s_container, command, geo_content, tolerance)

    # Extend line
    if command.extend >= 0:
        for i in range(command.iteration):
            tolerance = command.extend * (float(i + 1) / command.iteration)
            extend_line(s_container, command, geo_content, tolerance)

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
    nbr_extend: 0
    in_features: List[object]
    out_features: List[object]
    in_nbr_line_strings: 0
    out_nbr_line_strings: 0
    bounds: List[object] = None


geo_content = GeoContent(crs=None, driver=None, schemas={}, in_features=[], out_features=[],
                         nbr_xjunction=0, nbr_yjunction=0, nbr_join=0, nbr_noise=0, nbr_extend=0,
                         in_nbr_line_strings=0, out_nbr_line_strings=0, bounds=[])



# Read the command line arguments
#command = read_arguments()


#GenUtil.read_in_file (command.in_file, geo_content, [command.input_layer])


# test for Ycrossing
a = LineStringSc(((5,5),(5,10)))
b = LineStringSc(((0,4),(3,4),(5,5)))
c = LineStringSc(((5,5),(7,4),(10,4)))
d = LineStringSc(((5,10),(7,11), (10,11)))
e = LineStringSc(((0,11),(3,11),(5,10)))
lst_y_junction = [a,b,c,d,e]

# Test for X crossing
f = LineStringSc(((10,10), (10,15),(10,20)))
g = LineStringSc(((6,21),(8,21),(10,20)))
h = LineStringSc(((10,20),(12,21),(14,21)))
i = LineStringSc(((10,10), (8,9), (6,9)))
j = LineStringSc(((14,9),(12,9),(10,10)))
lst_x_junction = [f,g,h,i,j]

# Test for noise type 1
k = LineStringSc(((10,0), (30,0)))
l = LineStringSc(((30,0),(35,3)))
m = LineStringSc(((30,0),(35,-3)))
n = LineStringSc(((5,3),(10,0)))
o = LineStringSc(((10,0),(5,-3)))
lst_noise1 = [k,l,m,n,o]

# Test for join

p = LineStringSc(((20,0),(30,0)))
q = LineStringSc(((32,1),(35,0)))
r = LineStringSc(((19,1),(10,0)))
s = LineStringSc(((10,0),(0,0)))
t = LineStringSc(((10,0),(10,10)))
lst_join = [p,q,r,s,t]

# Test for extend
p = LineStringSc(((5,0),(20,0)))
q = LineStringSc(((4,-4),(4,0),(4,4)))
r = LineStringSc(((21,-4),(21,0),(21,4)))
lst_extend = [p,q,r]

# Test noise type  II

u = LineStringSc(((50,0),(0,0)))
v = LineStringSc(((50,0),(50,5)))
w = LineStringSc(((50,0),(55,0)))
x = LineStringSc(((55,-5),(55,0)))
y = LineStringSc(((60,0),(55,0)))
z = LineStringSc(((60,0),(65,-2)))
zz = LineStringSc(((65,2),(60,0)))

geo_content.in_features = lst_noise1
command = SpatialContainer()
command.yjunction = 5
command.xjunction = 10
command.join = 40
command.extend = 3
command.iteration=5
command.noise = 7

manage_cleaning(command, geo_content)




print ("-------")
print("Name of input file: {}".format(command.in_file))
print("Name of input layer: {}".format(command.input_layer))
print("Name of output_layer: {}".format(command.output_layer))
print("Tolerance for Y junction: {}".format(command.yjunction))
print("Tolerance for X junction: {}".format(command.xjunction))
print("Tolerance for noise detection: {}".format(command.xjunction))
print("Tolerance for extending lines: {}".format(command.extend))
print ("-----")
print("Number of features (line) read: {}".format(geo_content.in_nbr_line_strings))
print("Number of features (line) written: {}".format(geo_content.out_nbr_line_strings))
print("Number of Y junction corrected: {}".format(geo_content.nbr_yjunction))
print("Number of X junction corrected: {}".format(geo_content.nbr_xjunction))
print("Number of end line joined: {}".format(geo_content.nbr_join))
print("Number of line extended: {}".format(geo_content.nbr_extend))
print("Number of line (noise) deleted: {}".format(geo_content.nbr_noise))


# Copy the results in the output file
geo_content.layer_names = [command.output_layer]
for feature in geo_content.out_features:
    feature.sb_layer_name = command.output_layer

GenUtil.write_out_file_append (command.in_file, geo_content)
