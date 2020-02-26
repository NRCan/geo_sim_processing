# GeoSim

## Introduction

Sherbend is a geospatial simplification and generalization tool for lines and polygons.  Sherbend is an implementation and an improvement of the algorithm described in the paper "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this algorithm is that for each line it analyzes its bends (curves) and decides which one to simplify, trying to emulate what a cartographer would do manually to simplify or generalize a line.  Sherbend will accept points, lines and polygons as input.  Even though points cannot be simplified, they are used for topological relationship validations. Sherbend can accept GeoPackage and Esri Shape file as input/ouput but not a mixed of both.

## Requirements  
- Python 3.7 with the following libraries:
    - Shapely
    - Rtree
    - Fiona

## Installation on your workstation
Using conda, you can set and activate your python environment with the following commands:   
```
conda create --name YOUR_ENV python=3.7 shapely rtree fiona
source activate YOUR_ENV   (for Linux and macos)
activate YOUR_ENV          (for Windows)
```
Note on the installation:
  - For Windows users, it you are not using conda, do not forget that Shapely, Rtree and Fiona are all python wrapper of C libraries and need DLLs so use the appropriate installer (not just pip). This [site](https://www.lfd.uci.edu/~gohlke/pythonlibs/) contains a long list of Windows installers.

## Usage

usage: python sherbend.py \[-h] \[-eh] \[-ep] \[-pl] \[-d diameter | -dl dlayer] in_file out_file

positional arguments:
    
    in_file               Input Geopackage vector file to simplify (GPKG)
    out_file              Output Geopackage vector file simplified (GPKG)

optional arguments:

     -d , --diameter          Diameter of the minimum adjusted area bend to simplify (to remove)     
     -h, --help               Show this help message and exit
     -eh, --exclude_hole      Exclude (delete) polygon rings (interior holes) below the minimum adjusted area
     -ep, --exclude_polygon   Exclude (delete) polygons exteriors below the minimum adjusted area (delete also any interior holes if present)
     -pl, --per_layer         Analyze topology per layer only; this means features from different layers can overlap after simplification
     -dl, --dlayer            Specify the diameter of the minimum adjusted area bend to simplify per layer name (ex: -dl Road=5,Hydro=7.5)
     
Some example:

python sherbend.py -d 3 in_file.gpkg out\_file.gpkh
   
   - Simplify each feature of each layer of the input file (in_file.gpkg) with a bend diameter below 3 (in map unit) and create the output file out_file.gpkg
   
python sherbend.py -d 3 -pl in\_file.gpkg out_file.gpkh
   
   - Simplify each feature of each layer of the input file with a bend diameter below 3 and create the output file with each layer processed independently
   
python sherbend.py -d 3 -ep -eh in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file with a bend diameter below 3 and create the output file; delete the polygons including all their interiors if the exterior is below a bend diameter of 3; also delete the polygon interiors if the interior is below a bend diameter of 3
   
python sherbend.py -dl Road=3,Lake=5,River=0 in_file.gpkg out_file.gpkh

   - Simplify each feature of the Road, Lake and River layers of the input file with a bend diameter below 3 for the Road layer, 5 for the Lake layer and do not simplify the River layer features but use them for analysing the topology; finally create the output file

## Line Simplification versus Line Generalization

*Line Simplification* is the process of removing vertices in a line while trying to keep the maximum number of details within the line whereas *Line Generalization* is the process of removing meaningless (unwanted) details in a line usually for scaling down.  The well known Douglas-Peucker algorithm is a very good example of line simplification tool and Sherbend falls more in the category of line generalization tools. Keep in mind thay both algorithms can be complementary because Sherbend will not remove unnecessary vertices in the case of very high densities of vertices.  It may be a good idea to use Douglass Peucker before Sherbend in the case of very densed geometries.

## How it works

Sherbend will simplify (generalize) lines as well as polygons.  It will also take into account points, which are unsimplifiable, for analysis of topological relationships. Sherbend consists of three main steps: detect bends, determine which bends to simplify and preserve the topological (spatial) relationships.  These 3 steps are detailed below.

* __Detecting bends__
For each line and ring composing polygon features, Sherbend will detect the position of each bend.  Wang and Müller defined a bend as being the part of a line which contains a number of subsequent vertices with the inflection angles on all vertices being in opposite sign.
Figure 1a shows a line.  Figure 1b depicts the same line with inflexion signs on ech vertice.  Figure 1c shows the position of the 3 bends each forming an area.

* __Determining the bends to simplify__
For each bend of a line or polygon ring, Sherbend calculates an adjusted area value using the following formula: *\.75\*A/cmpi* where *A* is the area of the bend *(1)* and *cmpi* the compactness index of the bend.  The compactness index is computed using the following formula: *4\*π\*A/p\*\*2* where *A* is the area of the bend and *p* is the perimeter of the bend. The compactness index varies between \[0..1].  The more circular the bend, the closer the index to 1.  Conversely, the flatter the bend, the closer the index to 0.  The Sherbend parameter -d (ex.: -d 4) represents the diameter of a theoretical circle to define the minimum adjusted area value using *\.75\*2\*π\*r\*\*2/cmpi* where *r* is d/2.  Finally, each bend of a line that is below the minimum adjusted area value is replaced by a straight line.  Figure 1d shows the result with the middle bend of the line removed (simplified).

*(1)* The computations are always done in map unit: meters, feet, degrees...

![Figure1](/image/figure1.png)

* __Preserving topological relationship__
Before any bend simplifcation is applied, Sherbend will always analyze the following 3 topological relationships to ensure they are not affected by the simplification operation: simplicity, intersection and sidedness.  If simplification alters any of those relationships, then it is not performed.  Thereby Sherbend preserves the existing relative topology between the geospatial features to simplify.  

### Simplicity
Sherbend will not simplify a bend, if the simplified bend (dashed line in figure 2a) creates a self intersection.  

### Intersection
Sherbend will not simplify a bend, if the simplified bend creates an intersection between 2 existing features (figure 2b).  Conflicting features can be a line with another line or a line with a polygon ring.

### Sidedness
Sherbend will not simplify a bend, if simplifying the bend creates a sidedness or relative position error between 2 features. Two examples of sidedness issues are shown in figures 2c and 2d.  The preservation of the sidedness topological relationship is particulary important when it comes to simplifying polygon rings.  In figure 2c, simplifying the polygon as shown (dashed line) would make what was an inner hole "pop out" and become external to the new polygon.  In figure 2d, simplifying the bend in the line segment (dashed line) would result in the point feature changing its location relative to the original line.  If for example the original line represents a river and the point represents a building, it would mean the building would find itself on the other side of the river after simplification.   Conflicting features can be a line with a point or a line with a line or a line with a polygon ring.

Note: For any given line or polygon ring, only those bends the simplification of which do not cause any topological issues as expressed above will be simplified.

![Figure2](/image/figure2.png)

### Rule of thumb for the diameter
Sherbend can be used for line simplifying often in the context of line generalization. The big question will often be what diameter should we use?  A good starting point is the cartographic rule of thumb -- the *.5mm on the map* -- which says that the minimumm distance between two lines should be greater than 0.5mm on a paper map. So to simplify (generalize) a line for representation at a scale of 1:50 000 for example a diameter of 25m should be a good starting point... 

## Known issue with GeoPackage

The following problem can occur when using fiona libraries when creating GeoPackage.  It's a known issue, where the spatial index is not created for a specific layer.  The program still terminates with Exit Code 0 (meaning "success").  You can create the spatial index after in QGIS.

```
ERROR 1: sqlite3_exec(CREATE VIRTUAL TABLE "rtree_line_geom" USING rtree(id, minx, maxx, miny, maxy)) failed: table "rtree_line_geom" already exists
Traceback (most recent call last):
  File "fiona/_err.pyx", line 201, in fiona._err.GDALErrCtxManager.__exit__
fiona._err.CPLE_AppDefinedError: b'sqlite3_exec(CREATE VIRTUAL TABLE "rtree_line_geom" USING rtree(id, minx, maxx, miny, maxy)) failed: table "rtree_line_geom" already exists'
Exception ignored in: 'fiona._shim.gdal_flush_cache'
Traceback (most recent call last):
  File "fiona/_err.pyx", line 201, in fiona._err.GDALErrCtxManager.__exit__
fiona._err.CPLE_AppDefinedError: b'sqlite3_exec(CREATE VIRTUAL TABLE "rtree_line_geom" USING rtree(id, minx, maxx, miny, maxy)) failed: table "rtree_line_geom" already exists'

Process finished with exit code 0
```
