# GeoSim
Line simplification tools for python using shapely and fiona libraries

Introduction

Sherbend is a geospatial line simplification tool (another...).  It's the implementation of the algorithm from the papaer "Line 
Generalization Based on Analysis of Shape Characteristics, Zeshen Wand and Jean-Clsaude Müller, 1998" often known as "Bend Simplify" or 
"Wang Algorithm".  The particularity of this line simplification algorithm is that it analyses for a line each of it's curves and decide
which one to simplify as a cartographer would do to manually simplify a line.  Compared to Douglas-Peucker this algorithm tries to preserve
the maximum number of curves or bends (line definition) with the minimum number of vertices, Sherbend algorithm tries to remove unnecessary 
curves based on a tolerance (curve diameter

Requirements

Installation on your workstation

Exucution and options

How it works (Rule of thumb)

Topological relationships

Simplicity

Intersection

Sidedness
