from shapely.geometry import Point, LineString, Polygon
from time import time
from lib_geosim import GenUtil, PolygonSc
import math

from shapely import affinity
from shapely.ops import snap

from shapely.geometry import Point
from shapely.strtree import STRtree
points = []
for i in range(10):
    point = Point((i,i))
    point._t = i*100
    points.append(point)
tree = STRtree(points)
result = tree.query(Point(2,2).buffer(0.99))

0/0


class Foo(object):
    pass

class SubFoo(Foo):
    pass

a = SubFoo()

print (type(a))
print (a.__class__)
print (issubclass(a.__class__, SubFoo))
print (issubclass(a.__class__, Foo))

0/0


b = Point(0,0)
c = LineString(((2,-2),(2,2)))
a = c.project(b)


p1 = Point(0,0)
l1 =  LineString (((2,-5),(3,5) ))
result = snap(p1, l1, 15)

line = LineString([(0,0), (0.8, 0.8), (1.8, 0.95), (2.6, 0.5)])
start_time = time()
p1 = Point((5,5))
p2 = Point((3,6))
line1 = LineString(((0,3),(10,4)))
line2 = LineString(((5,5),(3,6)))

total = 0
for i in range(1):
    dist1 = line1.distance(p1)
    dist2 = line1.distance(p2)
    if dist1>dist2:
        total+=dist1
    else:
        total +=dist2

print("Distance: {}".format(total))
print ("Le temps 2: {}".format( time() - start_time) )

total = 0
start_time = time()
for i in range(1):
    dist = line2.hausdorff_distance(line1)
    total += dist
print("Distance: {}".format(total))
print ("Le temps 2: {}".format( time() - start_time) )

0/0

bend_i = 1
bend_j = 4
#line_ori = [(5,12),(10,10),(10,15),(15,15),(15,10),(20,12)]
line_ori = [(5,8),(10,10),(10,5),(15,5),(15,10),(20,12)]
for order, angle in ((0,0.), (0,45.), (0,90.), (0,135.), (0,180.), (0,225.), (0,270.), (0,315), (1,0.), (1,45.), (1,90.), (1,135.), (1,180.), (1,225.), (1,270.), (1,315)):
    if order == 0:
        line=LineString(line_ori)
    else:
        line_tmp = list(line_ori)
        line_tmp.reverse()
        line=LineString(line_tmp)
    a = math.cos(angle)
    b = -math.sin(angle)
    d = math.sin(angle)
    e = math.cos(angle)

    line1 = affinity.rotate(line, angle, line.coords[1])
    x, y = zip(*list(line1.coords))
    plt.plot(x, y)

    # Translate the line
    xoff, yoff = line1.coords[bend_i][0], line1.coords[bend_i][1]
    line2 = affinity.affine_transform(line1, [1, 0, 0, 1, -xoff, -yoff])

    x,y = zip(*list(line2.coords))
    plt.plot(x,y)
    x,y = zip(*list(line1.coords))
    plt.plot(x,y)

    line1_coord = list(line2.coords)
    p0_x = line1_coord[bend_j][0]
    p0_y = line1_coord[bend_j][1]
    p1_x = abs(p0_x) + 1.  # In case x == 0
    p1_y = 0.

    dot = p0_x*p1_x + p0_y*p1_y
    len_a = (p0_x**2+p0_y**2)**.5
    len_b = (p1_x**2+p1_y**2)**.5

    angle = math.acos(dot/(len_a*len_b))
    angle =  (angle*180/math.pi)
    print (angle)

    if p0_y >= 0.:
        angle = -angle
    a = math.cos(angle)
    b = -math.sin(angle)
    d = math.sin(angle)
    e = math.cos(angle)

    line4 = affinity.rotate(line2, angle, origin=(0,0))

    x,y = zip(*list(line4.coords))
    plt.plot(x,y)

    lst_coords = list(line4.coords)
    line_i = LineString(lst_coords[0:3])
    line_j = LineString(lst_coords[-2:])
    theta_i = lib_geobato.GenUtil.compute_angle(lst_coords[0], lst_coords[1], lst_coords[bend_j])
    theta_j = lib_geobato.GenUtil.compute_angle(lst_coords[bend_j], lst_coords[-2], lst_coords[-1])

    (minx, miny, maxx, maxy) = line4.bounds
    y_dynamic = (abs(miny) + abs(maxy))* 10.
    x_middle = (lst_coords[bend_i][0] + lst_coords[bend_j][0]) / 2.
    line_y_positive = LineString(((x_middle,0),(x_middle, y_dynamic)))
    line_y_negative = LineString(((x_middle, 0), (x_middle, -y_dynamic)))
    if line4.crosses(line_y_positive):
        bend_side = +1
    else:
        if line4.crosses(line_y_negative):
            bend_side = -1

    if lst_coords[0][1] >= 0.:
        start_line_side = 1
    else:
        start_line_side = -1

    if lst_coords[-1][1] >= 0.:
        end_line_side = 1
    else:
        end_line_side = -1

    if (start_line_side*end_line_side == -1):
        print ("Nothing to do....")
        line5 = LineString(lst_coords[0:bend_i+1] + lst_coords[bend_j:])
    else:
        # Both line are on the same side
        if start_line_side == 1 and end_line_side == 1:
            if bend_side == -1:
                angle_bias = 2.
                y_offset = -1
            else:
                angle_bias = 3.
                y_offset = -1
        if start_line_side == -1 and end_line_side == -1:
            if bend_side == 1:
                angle_bias = 2.
                y_offset = 1
            else:
                angle_bias = 3.
                y_offset = 1

        theta_i = (180. - theta_i) / angle_bias
        if theta_i >= 5.:
            hypothenus = x_middle / math.cos(theta_i*math.pi/180.)
            y_height = math.sqrt(hypothenus**2 - x_middle**2)
            if bend_side == -1:
                y_height *= y_offset
            new_coord = (x_middle, y_height)
            line5 = LineString(lst_coords[0:bend_i+1] + [new_coord] + lst_coords[bend_j:])
        else:
            print("Nothing to do....")
            line5 = LineString(lst_coords[0:bend_i+1] + lst_coords[bend_j:])

    x, y = zip(*list(line5.coords))
    plt.plot(x, y)
    plt.xlim(-20, 20)
    plt.ylim(-20, 20)
    plt.show()




pol = Polygon([(0,0),(0,10), (10,10),(10,0)], [[(5,5),(5,15),(8,15),(8,5)]])
print (pol.is_simple)
print (pol.area)
line = LineString( ((-1,-1), (11,11) ) )
print (pol.crosses(line))
coord = list(pol.coords)
int = list(pol.interiors)
pol.interiors = []
int_coords = pol.interiors.coords


print (pol)
0/0



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