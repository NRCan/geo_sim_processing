# GeoSim

GeoSim is a set of tools that aims to simplify/generalize line and polygon features. It is composed of 3 tools: [Sherbend](#Sherbend), [Chordal Axis](#Chordal-Axis) and [TopoSim](#TopoSim)

## Requirements  
- Python 3.7 with the following libraries:
    - [Shapely](https://pypi.org/project/Shapely/)
    - [Rtree](https://pypi.org/project/Rtree/)
    - [Fiona](https://pypi.org/project/Fiona/)

## Installation on your workstation
Using conda, you can set and activate your python environment with the following commands:   
```
conda create --name YOUR_ENV python=3.7 shapely rtree fiona
source activate YOUR_ENV   (for Linux and macos)
activate YOUR_ENV          (for Windows)
```
Note on the installation:
  - For Windows users, it you are not using conda, do not forget that Shapely, Rtree and Fiona are all python wrapper of C libraries and need DLLs so use the appropriate installer (not just pip). This [site](https://www.lfd.uci.edu/~gohlke/pythonlibs/) contains a long list of Windows installers.
  - For TopoSim, the rtree library is not needed

## Known issue with GeoPackage

The following problem can occur when using fiona libraries when creating GeoPackage or layers in GeoPackage when using GeoSim tools.  It's a known issue, where the spatial index is not created for a specific layer.  The program still terminates with Exit Code 0 (meaning "success").  You can create the spatial index after in QGIS.

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

# Sherbend

Sherbend is a geospatial simplification and generalization tool for lines and polygons.  Sherbend is an implementation and an improvement of the algorithm described in the paper "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this algorithm is that for each line it analyzes its bends (curves) and decides which one to simplify, trying to emulate what a cartographer would do manually to simplify or generalize a line.  Sherbend will accept points, lines and polygons as input.  Even though points cannot be simplified, they are used for topological relationship validations. Sherbend can accept GeoPackage and Esri Shape file as input/ouput but not a mixed of both.

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

python sherbend.py -d 3 in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file (in_file.gpkg) with a bend diameter below 3 (in map unit) and create the output file out_file.gpkg

python sherbend.py -d 3 -pl in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file with a bend diameter below 3 and create the output file with each layer processed independently

python sherbend.py -d 3 -ep -eh in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file with a bend diameter below 3 and create the output file; delete the polygons including all their interiors if the exterior is below a bend diameter of 3; also delete the polygon interiors if the interior is below a bend diameter of 3

python sherbend.py -dl Road=3,Lake=5,River=0 in_file.gpkg out_file.gpkh

   - Simplify each feature of the Road, Lake and River layers of the input file with a bend diameter below 3 for the Road layer, 5 for the Lake layer and do not simplify the River layer features but use them for analysing the topology; finally create the output file

## Line Simplification versus Line Generalization

*Line Simplification* is the process of removing vertices in a line while trying to keep the maximum number of details within the line whereas *Line Generalization* is the process of removing meaningless (unwanted) details in a line usually for scaling down.  The well known Douglas-Peucker algorithm is a very good example of line simplification tool and Sherbend falls more in the category of line generalization tools. Keep in mind thay both algorithms can be complementary because Sherbend will not remove unnecessary vertices in the case of very high densities of vertices.  It may be a good idea to use Douglass Peucker before Sherbend in the case of very densed geometries.

## How it works

Sherbend will simplify (generalize) lines as well as polygons.  It will also take into account points, which are unsimplifiable, for analysis of topological relationships. Sherbend consists of three main steps: detect bends, determine which bends to simplify and preserve the topological (spatial) relationships.  These 3 steps are detailed below.

* __Detecting bends__ -
For each line and ring composing polygon features, Sherbend will detect the position of each bend.  Wang and Müller defined a bend as being the part of a line which contains a number of subsequent vertices with the inflection angles on all vertices being in opposite sign.
Figure 1a shows a line.  Figure 1b depicts the same line with inflexion signs on ech vertice.  Figure 1c shows the position of the 3 bends each forming an area.

* __Determining the bends to simplify__ -
For each bend of a line or polygon ring, Sherbend calculates an adjusted area value using the following formula: *\.75\*A/cmpi* where *A* is the area of the bend *(1)* and *cmpi* the compactness index of the bend.  The compactness index is computed using the following formula: *4\*π\*A/p\*\*2* where *A* is the area of the bend and *p* is the perimeter of the bend. The compactness index varies between \[0..1].  The more circular the bend, the closer the index to 1.  Conversely, the flatter the bend, the closer the index to 0.  The Sherbend parameter -d (ex.: -d 4) represents the diameter of a theoretical circle to define the minimum adjusted area value using *\.75\*2\*π\*r\*\*2/cmpi* where *r* is d/2.  Finally, each bend of a line that is below the minimum adjusted area value is replaced by a straight line.  Figure 1d shows the result with the middle bend of the line removed (simplified).

*(1)* The computations are always done in map unit: meters, feet, degrees...

![Figure1](/image/figure1.png)

* __Preserving topological relationship__ -
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

## Known issue with GeoPackage format
# Chordal Axis

ChordalAxis is a geospatial tool that takes triangles, usually the result of a constraint Delauny trianglulation and creates a skeleton (the center line).  ChordalAxis is an improvement of the algorithm based of the paper "Rectification of the Chordal Axis Transform and a New Criterion for Shape
Decomposition", Lakshman Prasad, 2005".

## Medial Axis Versus Chordal Axis

The skeleton (center line) is a linear feature representation of a polygonized feature. In computational geometry, it is known as the medial axis and many algorithms are approximating it very well.  A major issue with those algorithms is the possible instability for very irregular complex polygons such as dense river or road network polygons. (Figure 4).  The Chordal Axis has shown excellent stability in very very polygons while extracting a very representative skeleton.

## Usage

usage: python Chordal_Axis.py \[-h] \[-t] \[-s] \ file

positional arguments:

    file                  Input/Output Geopackage vector file to extract skeleton

optional arguments:

     -h, --help          Show this help message and exit
     -t, --triangle      Name of the layer in the graphic file containing the triangle (Line string)
     -s, --skeleton      Name of the layer to create that will contain the skeleton
     -c, --correction    Correct the skeleton for small centre line, T junction and X junction


Example:

python chordal_axis.py -t tesselation -s skeleton road.gpkg

   - Load the triangle in the layer named tesselation; create the centre line using the chordal axis and create the layer skeleton in the file road.gpkg

 ## How it works

A user will probably create the triangulation from a set of polygons using a constraints Delaunay triangulation tool.  Delaunay triangulation is known to describe polygons well and to  be very robust and stable.  The resulting triangles are the input for the Chordal Axis program.  The Chordal Axis alogorithm will analyze each triangle, determine its type based on the number of adjacent triangles and build the appropriate skeleton (centre line).  All triangles fall within one of the following four types: 1)  _isolated triangle_, when a triangle has no adjacent triangle; 2) _terminal triangle_, when a triangle has only one adjacent triangle; 3) _sleeve triangle_, when a triangle has 2 adjacent triangles; 4) _junction triangle_, when a triangle has 3 adjacent triangles.  Each of the four triangle types will produce a specific centre line.  For the _isolated triangle_, (Figure 3a) no center line (degenerated case) is created; for the _terminal triangle_ (Figure 3b) the mid point of the adjecent side is connected with the opposite angle; for _sleeve triangle_ (Figure 3c) the mid point of the two adajcent sides are connected; for the _junction triangle_ (Figure 3d) the mid points of each side are connected to the centre point of the triangle.  After centre line creation all the centre lines are merged together.  The Chordal Axis transform will preserve [Simplicity](#Simplicity) and [Intersection](#Intersection) topological relationships between the lines forming the skeleton and the outer and inner boundaries of the polygon defined by the Delaynay triangulation.

![figure3](/image/figure3.png)

## Correction
The Chordal Axis algorithm gives a very good approximation of the true medial axis of a polygon but it produces unwanted artifacts when it creates the skeleton especially in the case of long and narrow polygons (figure 5) . The main artifact types are: meaningless small centre line (figure 4a); wrongly formed "T junctions" (figure 4c) and "X junctions" (crossing junction) (figure 4e).  When the correction parameter  is used (-c or --correction), the skeleton will be pruned of the meaningless small centre line (4b); it will correct "T junctions" and rectify the normal direction of the line (figure 4d); and, it will rectify the "X crossing" by merging two T junctions that are adjacent (figure 4f).

![figure5a](/image/figure5a.png "Figure 5a") ![figure5b](/image/figure5b.png "Figure 5b") ![figure5c](/image/figure5c.png "Figure 5c") ![figure5d](/image/figure5d.png "Figure 5d") ![figure5e](/image/figure5e.png "Figure 5e") ![figure5f](/image/figure5f.png "Figure 5f")

  Figure 4a    Figure 4b     Figure 4c    Figure 4d    Figure 4e    Figure 4f

## Rule of thumb for the use of Chordal Axis
Chordal Axis can be used for skeleton extraction and polygon to line transformation in the context of polygon generalization. Often the quality of the skeleton produced will depend on the density of polygon vertices and therefore overall quantity of triangles ingested by Chordal Axis : the more vertices, the higher the number of generated triangles and the better the skeleton (at the price of increased computation time).  Equilateral triangles produce the best skeleton while highly obtuse and/or acute triangles will produce a jagged line that can then be simplified.  The vertex density should not result in either over- or under-simplified features.  Delaunay triangulation and Chordal Axis will give excellent results in very complex situations like a densely polygonized road network such as the one shown in Figure 4 in which all road segments belong to the same polygon!

![figure4](/image/figure4.png)

Figure 5


# TopoSim

TopoSim is a geospatial simplification tool for lines and polygons.  TopoSim implements [Shapely](https://pypi.org/project/Shapely/)'s *simplify* tool with parameter  *preserve_topology=True*. For line and polygon simplification Shapely implements an algorithm similar to the [Douglas Peucker algorithm](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm).  The implementation preserves the topology within one feature but not between features of the same layer or from different layers.  There is also a [known bug](https://locationtech.github.io/jts/javadoc/org/locationtech/jts/simplify/TopologyPreservingSimplifier.html) where the algorithm may create invalid topologies if there are components which are small relative to the tolerance value.   In particular, if a small interior hole is very close to an edge, simplification may result in the hole being moved outside the polygon (figure 6a). Similarly, a small polygon close to a larger one may end up being swallowed into the larger polygon.  Toposim will detect situations like figure 6b where one or more rings (interior parts) fall outside the polygon after being simplified and make the polygon invalid. The algoritm will remove (delete) these ring(s) so the feature remains valid after simplification.

Note: While most GIS tools will handle and display invalid geometries like figure 6b, some spatial operation will not be allowed and this is why it's important to keep validity of the geometry after a spatial operation.

![figure6a](/image/figure6a.png "Figure 6a") ![figure6b](/image/figure6b.png "Figure 6b")

        Figure 6a                    Figure 6b


## Usage

usage: python toposim.py \[-h] \[-t tolerance | -tl tlayer] in_file out_file

positional arguments:

    in_file               Input Geopackage vector file to simplify (GPKG)
    out_file              Output Geopackage vector file simplified (GPKG)

optional arguments:

     -t , --tolerance         Tolerance for the line simplification (usage similar to Douglas Peucker)     
     -h, --help               Show this help message and exit
     -tl, --tlayer            Specify the tolerance for the line simplification per layer name (ex: -tl Road=5,Hydro=7.5)

Examples:

```python
python toposim.py -t 3 in_file.gpkg out_file.gpkh
```

Simplify each feature in *in_file.gpkg* with a tolerance of 3 and create *out_file.gpkg*

```python
python toposim.py -tl Road=3,Lake=5 in_file.gpkg out_file.gpkh
```

Simplify each feature of the Road layer with a tolerance of 3 and Lake layers with a tolerance of 5.

## How it works

Toposim is an excellent tool to remove vertices on features with high vertex densities.  Try it with small tolerance value and then use [Sherbend](#Sherbend) to [generalize features](##Line Simplification versus Line Generalization).
