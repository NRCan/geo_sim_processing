# GeoSim
Line simplification tools for python using shapely and fiona libraries

Introduction

Sherbend is a geospatial line simplification tool (another...).  It's the implementation of the algorithm from the papaer "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wand and Jean-Clsaude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm".  The particularity of this line simplification algorithm is that it analyses for a line each of it's curves and decide which one to simplify as a cartographer would do to manually simplify a line.  Compared to Douglas-Peucker this algorithm tries to preserve the maximum number of curves or bends (line definition) with the minimum number of vertices, Sherbend algorithm tries to remove unnecessar curves based on a tolerance (curve diameter

## Requirements  
- Python 3.7 with the following libraries:
    - Shapely
    - Rtree
    - fiona

## Installation on your workstation
Using conda, you can set and activate your python environment with the following commands: 
    ```shell
    conda create -p YOUR_PATH python=3.7 shapely rtree
    source activate YOUR_ENV
    pip install fiona
    ```
  Note on the installation
  - fiona needs to be installed separatly has there is a problem (wtih conda?) when you try to installes shapely, rtree, fiona at the same time
  - for Windos users, do not forget that shapely, rtree and fiona are all python wrapper of C libraries and need DLLs so use the appropriate installer

Exucution and options

How it works (Rule of thumb)

Topological relationships

Simplicity

Intersection

Sidedness
