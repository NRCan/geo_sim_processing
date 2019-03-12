from shapely.geometry import Point, LineString
from time import time
import timeit

diag = LineString(((1,1),(2,2), (2,10)))
diag1 = LineString(((1,1),(2,2), (2,10)))
print (diag1)
diag1.coords = ((0,0),(10,10))
print (diag1)
diag1.coords[0] = (3,3)
a = diag.minimum_rotated_rectangle


loop = 100000
dist = 0.
start_time = time()
for i in range(loop):
    x1,y1,x2,y2 = i, i+1, i+2, i+3
    dist += ((x2-x1)**2 + (y2-y1)**2)**.5
print ("Distance: {}".format(dist))
print ("Le temps 1: {}".format( time() - start_time) )

dist = 0.
start_time = time()
for i in range(loop):
    x1,y1,x2,y2 = i, i+1, i+2, i+3
    dist += Point(x1,y1).distance(Point(x2,y2))
print("Distance: {}".format(dist))
print ("Le temps 2: {}".format( time() - start_time) )

0/0

lst = []
coords = []
x,y = 0,0
for i in range(100):
    coords.append((x,y))
    x+=1
    y+=1

loop=10000
start_time = time()
for i in range(loop):
    line = LineString(coords)
    lst.append(line)
print ("Le temps 1: {}".format( time() - start_time))

cpt=0
cpt_coord=0
start_time = time()
for line in lst:
    if line.geom_type=='LineString':
        cpt+=1
    coords = list(line.coords)
    for i in range(100):
        if coords[0] != None:
           cpt_coord+=1

print(cpt)
print (cpt_coord)

print ("Le temps 2: {}".format( time() - start_time) )

lst=[]
for i in range(loop):
    line = LineString(coords)
    line._gbt_type = 'LineString'
    lst.append(line)
print ("Le temps 3: {}".format( time() - start_time) )

cpt=0
cpt_coord=0
start_time = time()
for line in lst:
    if line._gbt_type=='LineString':
        cpt+=1
    for i in range(100):
        if line.coords[0] != None:
           cpt_coord+=1
print(cpt)
print(cpt_coord)

print ("Le temps 4: {}".format( time() - start_time) )