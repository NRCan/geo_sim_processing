#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Use to test if the geopackages passed in parameters are equal.
Two geopackages are equal if each feature has an equivalent in the other geopackage.
Features are tested for equality with the shapely method object.almost_equals
Order of the feature in each geopackage does need not to be preserved
"""

import sys
from dataclasses import dataclass
from typing import List
from lib_geosim import GenUtil

# Extract name of the file from the parameter line
primary = sys.argv.pop()
secondary= sys.argv.pop()

# It is important to import the unit test after we extract the arguments otherwise it will not work
import unittest

class TestFeatures(unittest.TestCase):

    def test_per_feature(self):
        """Test if Primary and Secondary are identical

        Validate if each feature of each layer of the primary geopackage as an equivalent in the in the
        secondary geopackage

        """

        # Loop over each feature in the primary geopackage
        nbr_err = 0
        for feature_pr in geo_content_pr.in_features:
            identical = False
            # Loop to find equivalent feature in the secondary geopackage
            for feature_sc in geo_content_sc.in_features:
                if feature_pr.almost_equals(feature_sc):
                    identical = True
                    break

            if not identical:
                nbr_err += 1

        assert nbr_err == 0, "Number of feature in error: {}".format(nbr_err)


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
    in_nbr_holes: 0
    bounds: List[object] = None
    layer_names: List[object] = None
    in_features: List[object] = None
    out_features: List[object] = None


geo_content_pr = GeoContent(crs=None, driver=None, schemas={}, bounds=[], layer_names=[], in_features=[],
                            out_features=[],
                            in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, in_nbr_holes=0)

geo_content_sc = GeoContent(crs=None, driver=None, schemas={}, bounds=[], layer_names=[], in_features=[],
                            out_features=[],
                            in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, in_nbr_holes=0)

# Read and load the layers of the primary geopackage
GenUtil.read_in_file(primary, geo_content_pr, None)

# Read and load the layers of the secondary geopackage
GenUtil.read_in_file(secondary, geo_content_sc, None)

if __name__ == '__main__':
    unittest.main()



