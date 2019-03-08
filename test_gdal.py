from shapely.geometry import LineString
from time import time
import timeit

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
print ("Le temps 1: {}".format( time() - start_time) )

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