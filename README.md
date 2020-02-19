# GeoSim

## Introduction

Sherbend is a geospatial simplification and generalization tool for lines and polygons.  Sherbend is an implementation and an improvement of the algorithm described in the paper "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this algorithm is that for each line it analyses its bends (curves) and decide which one to simplify, trying to emulate what a cartographer would do manually to simplify or generalize a line.  Sherbend will accept as input point, line and polygon (points are unsimplifiable but are used for topological relationship validation). Sherbend can accept GeoPackage and Esri Shape file as input/ouput but not a mixed of both.

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
     -eh, --exclude_hole      Exclude (delete) polygon ring (interior holes) below the minimum adjusted area
     -ep, --exclude_polygon   Exclude (delete) polygons exterior below the minimum adjusted area (delete also any interior holes if present)
     -pl, --per_layer         Analyse topology per layer only which mean features from different layers can overlap after simplification
     -dl, --dlayer            Specify the diameter of the minimum adjusted area bend to simplify per layer name (ex: -dl Road=5,Hydro=7.5)
     
Some example:

python sherbend.py -d 3 in_file.gpkg out\_file.gpkh
   
   - Simplify each feature of each layer of the input file (in_file.gpkg) with a bend diameter of 3 (in map unit) and create the output file out_file.gpkg
   
python sherbend.py -d 3 -pl in\_file.gpkg out_file.gpkh
   
   - Simplify each feature of each layer of the input file with a bend diameter of 3 and create the output file but each layer are processed independently
   
python sherbend.py -d 3 -ep -eh in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file with a bend diameter of 3 and create the output file; delete the polygon including all the interiors if the exterior is below the minimum adjusted area also delete the polygon interiors if the interior is below the minimum adjusted area
   
python sherbend.py -dl Road=3,Lake=5,River=0 in_file.gpkg out_file.gpkh

   - Simplify each feature of the Road, Lake and River layers of the input file with a bend diameter of 3 for the Road layer, 5 for the Lake layer  and do not simplify the River layer features but use them for analysing the topology; finally create the output file out_file.gpkg

## Comparison with other simplication tool

Compared to the well known Douglas-Peucker algorithm, Sherbend algorithm will always try to remove unnecessary bends (line details) based on a bend diameter.  Douglas-Peucker on the other hand will always try to preserve the maximum number of line details (line definition) with the minimum number of vertices.  both algorithm are complimentary because Sherbend will not remove unnecessary vertice.

## How it works

Sherbend will simplify (generalize) line and polygon it also take into account point which are unsimplifiable and only used when analysing topological relationships. The 3 main steps of the algorithm are: Detecting the bends, Determining the bends to simplify and Preserving the topological (spatial) relationships.  These 3 steps are detailed below.

* __Detecting bends__
For each line and rings composing polygon Sherbend will detect the position of each bend.  Wang and Müller defined a bend as being the part of a line which contains a number of susequent vertices, with the inflections angles on all vertices being in opposite sign.
Figure 1a show a line, figure 1b the same line with inflexion sign on ech vertice, figure 1 c the same line with the position of the 3 bends forming each an area.

* __Determining the bends to simplify__
For each bend of a line or polygon ring Sherbend calculates an adjusted area value using the following formula: *\.75\*A/cmpi* where *A* is the area of the bend in map uni and *cmpi* the compactness index of the bend.  The compactness index is calculate using the following area: *4\*π\*A/p\*\*2* where *A* is the area of the bend and *p* is the perimeter of the bend. The compactness index vary between \[0..1] with a circular bend will value near 1 and an almost flat bend having a value near 0.  The Sherbend parameter -d (ex.: -d 4) represent the diameter of a theoritical circle that permit to define the minimum adjusted area value using *\.75\*2\*π\*r\*\*2/cmpi* where *r* is d/2.  Finally, each bend of a line that are below the minimum adjusted area value are replaced by a straight line.  Figure 1d represent the result with the middle bend of the line simplified.

![Figure1](/image/figure1.png)

* __Preserving topological relationship__
Before any bend simplifcation, Sherbend will analayse the following 3 topological relationship: simplicity, intersection and sidedness; if one of the topological relationship is not valid that particulatbend is not simplified.  This process preserve the existing relative topology between the geospatial features to simplify.  

### Simplicity
Sherbend will not simplify a bend, if the simplified bend creates a self intersection in the line (figure 2a).  

### Intersection
Sherbend will not permit bend simplification if the simplified bend creates an intersection between 2 features (figure 2b).  The features in conflict can be a line with a line or a line with a polygon ring.

### Sidedness
Sherbend will not permit bend simplification if the simplified bend creates a sidedness or relative position error between 2 features. Like a building that change side in regards with a river after simplification (figure 2c and 2d).  The features in conflict can be a line with a point or a line with line or a line with a polygon ring.  The analysis of this topological relationship is particulary important when it comes to simplify polygon ring.  In order to prevent interior rings to "pop out" outside its exterior ring afiter the bend simplification (figure x).

Note: For all 3 topological relationships, for any given a line or polygon ring if one or more of its bend simplification create topological error these bend will not be simplified but all the bend that do not create topological errors will be simplified.

![Figure2](/image/figure2.png)

### Rule of thumb for the diameter
Shebend will be used for line simplifying often in the context of map generalization. The big question will often be what diameter should we use?  A good starting point is the cartogrphic rule of thumb of the *.5mm on the map* which say that the minimumm distance between two lines should be greater than 0.5mm on a paper map. So to simplify (generalize) a line in order to acheive 1:50 000 on the map a diameter of 25 should be a good starting point... 

## Known issue with GeoPackage

The following problem can occured when creating a GeoPackage, it's a know issue, the spatial index is not created for a specific layer.  The program still terminate with Exit Code 0.  You can create the spatial index after in QGIS.

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
