# GeoSim
Line simplification tools for python using shapely and fiona libraries

Introduction

Sherbend is a geospatial line and polygon simplification tool.  Sherbend is the implementation of the algorithm from the papaer "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wand and Jean-Clsaude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this line simplification algorithm is that it analyses for a line each of it's curves and decide which one to simplify as a cartographer would do to manually simplify a line.  Compared to Douglas-Peucker this algorithm tries to preserve the maximum number of curves or bends (line definition) with the minimum number of vertices, Sherbend algorithm tries to remove unnecessar curves based on a tolerance (curve diameter

## Requirements  
- Python 3.7 with the following libraries:
    - Shapely
    - Rtree
    - Fiona

## Installation on your workstation
Using conda, you can set and activate your python environment with the following commands: 
    ```shell
    conda create -p YOUR_PATH python=3.7 shapely rtree
    source activate YOUR_ENV
    pip install fiona
    ```
    
  Note on the installation:
  - fiona needs to be installed separatly has there is a problem (wtih conda?) when you try to installes shapely, rtree, fiona at the same time
  - for Windos users, do not forget that shapely, rtree and fiona are all python wrapper of C libraries and need DLLs so use the appropriate installer (not just pip)

##Usage

usage: sherbend.py \[-h] \[-eh] \[-ep] \[-pl] \[-d DIAMETER | -dl DLAYER] in_file out_file

positional arguments:
  in_file               input vector file to simplif
  out_file              output vector file simplified

optional arguments:
     -d , --diameter          diameter of the minimum adjusted area bend to simplify (to remove)
     
     -h, --help               show this help message and exit
  
     -eh, --exclude_hole      for polygon exclude (delete) polygon holes (interior) below the minimum adjusted area
  
     -ep, --exclude_polygon   for polygon exclude polygons exterior below the minimum adjusted area (delete also any interior if present
  
     -pl, --per_layer         evaluate topology per layer only (features from different layers can overlap after simplification)
                        
     -dl, --dlayer            specify the diameter of the minimum adjusted area bend to simplify per layer name (ex: -dl Road=5,Hydro=7.5)
     
Some example:

python sherbend.py -d 3 in_file.gpkg out\_file.gpkh
   
   - Simplify each feature of each layer of the input file (in_file.gpkg) with a diameter of 3 (in map unit) and create the output file out_file.gpkg
   
python sherbend.py -d 3 -pl in\_file.gpkg out_file.gpkh
   
   - Simplify each feature of each layer of the input file with a diameter of 3 and create the output file out_file.gpkg but each layer are processed independently
   
python sherbend.py -d 3 -ep -eh in_file.gpkg out_file.gpkh

   - Simplify each feature of each layer of the input file with a diameter of 3 and create the output file out_file.gpkg delete the polygon including all the interiors if the exterior is below the minimum adjusted area also delete the polygon interiors if the interior is below the minimum adjusted area
   
python sherbend.py -dl Road=3,Lake=5,River=0 in_file.gpkg out_file.gpkh

   - Simplify each feature of the Road, Lake and River layers of the input file with a diameter of 3 for the Road layer, 5 for the Lake layer  and do no simplify for the River layer features but use them for topology constraints; finally create the output file out_file.gpkg

##How it works

Sherbend will simplify line and polygon it also take into account point which are unsimplifiable but used when analysing and validating topological relationships.

* __Detecting bends__
For each line and rings composing polygon Sherbend will detect the position of each bend.  Wang and Müller defined a bend as being as the part of a line which contains a number of susequent vertices, with the inflections angles on all vertices being in opposite sign.
Figure 1 a show a line, figure 1b the same line with inflexion sign on ech vertice, figure 1 c the same line with the position of the bends.

* __Calculating adjusted area__
For each bend Sherbend calculates the adjusted area of each bend with the following formula: *\.75\*A/cmpi* where *A* is the area in map unit of the bend and *cmpi* the compactness index of the bend.  The copactness index is calculate with *4\*π\*A/(p\*p)* where *A* is the area and *p* is the perimeter of the bend in map unit. The compactness vary between \[0..1] with a circle having a value of 1 and an almost flat bend having a value of 0.  The parameter -d (ex.: -d 4) represent the diameter of a therotical circle that help  

* __Validating topological relationship__


Bend simplificatin

Topological 

##Topological relationships

ffff

###Simplicity

ddd

###Intersection

ddd

###Sidedness

ddd

### Rule of thimb
