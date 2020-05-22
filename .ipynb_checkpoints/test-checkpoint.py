from shapely.geometry import LineString, Point, Polygon, MultiLineString
from shapely.strtree import STRtree
#import random, time
#from rtree import index
from cmath import rect, phase
from math import radians, degrees
from shapely.ops import polygonize, polygonize_full, split


line = LineString(((0,0),(5,0),(10,0)))
a = line.project(Point(5,10))

ext = [(0,0), (0,10),(3,10),(3,14),(7,14),(7,10),(10,10), (10,0),(0,0)]
int_1 = [[(4,11), (4,12),(6,12),(6,11), (4,11)]]
int_2 = [[(4,9.5), (4,10.5),(6,10.5),(6,9.5), (4,9.5)]]

pol = Polygon(ext, int_1)
val = pol.is_valid
pol_s = pol.simplify(5, preserve_topology=True)
val = pol.is_simple



pol = Polygon(ext, int)
val = pol.is_valid
pol_s = pol.simplify(5, preserve_topology=True)
pol_s1 = pol.simplify(5, preserve_topology=False)
val2 = pol_s1.is_valid
pol = Polygon([(((0,0), 0,10),(10,10), (10,12), (12,12), (12,10), (20,10)), ((11,0), (11,4),(12,4),(12,5),(11,5), (11,7), (11,11))])
mline1 = mline.simplify(2, preserve_topology=True)


mline = MultiLineString([((0,10),(10,10), (10,12), (12,12), (12,10), (20,10)), ((11,0), (11,4),(12,4),(12,5),(11,5), (11,7), (11,11))])
mline1 = mline.simplify(2, preserve_topology=True)

line = LineString([(0,0),(10,0), (11,0), (11,1), (10,1), (20,0), (20,1), (22,1), (22,0), (30, 0), (30,-5), (21,-5), (21,.5)])
line1 = line.simplify(1, preserve_topology=True)



line = LineString([(0,0),(2,2), (4,-2), (6,2), (7,0), (8,0)])
splitter = LineString([(0,0),(8,0)])

line_a = LineString([(0,0),(10,0)])
line_b = LineString([(10,110),(20,20)])
line_c = LineString([(7,0),(9,0)])

line_split = split(line, splitter)
splitter_split = split(splitter, line)

pol = polygonize_full([line_split, splitter_split])

val_a = line_a.intersects(pol[0])
val_b = line_b.intersects(pol[0])
val_c = line_c.intersects(pol[0])


coords = [(0, 0), (0, 2), (1, 1), (2, 2), (2, 0), (1, 1), (0, 0)]
coords = [(0,0),(5,0), (10,0), (10,10), (5,10), (5,0), (5,-10),(0,-10), (0,0)]
bowtie = Polygon(coords)
va11 =  bowtie.is_valid
clean = bowtie.buffer(0)
val2 = clean.is_valid

l_a = [(0,0),(1,3),(2,-3),(10,10), (0,0)]
pol = Polygon (l_a)
pol1 = pol.buffer(0)
l_b = [(0,0),(10,10)]



line_a = LineString(l_a)
line_b = LineString(l_b)

p = polygonize_full([line_a,line_b])

p0 = [(0,0),(10,0),(10,10),(0,10),(0,0)]
p1 = [(0,10),(10,10),(10,20),(0,20),(0,10)]
p2 = [(0,0),(20,20),(-5,20),(-5,17),(5,17),(5,12),(-5,12),(0,0)]

pol1 = Polygon(p1)
pol2 = Polygon(p2)

pols = pol1.symmetric_difference(pol2)

p1 = [(0,0),(10,0),(10,10),(0,10),(0,0)]
p2 = [(2,2),(8,2),(8,8),(2,8),(2,2)]

pol1 = Polygon(p1)
pol2 = Polygon(p2)

pols1 = pol1.symmetric_difference(pol2)
pols2 = pol2.symmetric_difference(pol1)

p1 = [(0,0),(7,0),(7,7),(0,7),(0,0)]
p2 = p1 = [(0,10),(10,10),(10,20),(0,20),(0,10)]

pol1 = Polygon(p1)
pol2 = Polygon(p2)

pols1 = pol1.symmetric_difference(pol2)
pols2 = pol2.symmetric_difference(pol1)
0/0
0/0

def mean_angle(deg):
    a = None
#    a = sum(rect(1, radians(d)) for d,l in deg)
    for d,ll in deg:
        if a is None:
            a = rect(ll, radians(d))
        else:
            a +=  rect(ll, radians(d))

    b = phase(a)
    c = degrees(b)

    d =  degrees(phase(sum(rect(l, radians(d)) for d,l in deg)/len(deg)))

    return (c,d)

def mean_angle2(degrees):

    angle_sum = 0.
    tot_len = 0

    a = sum(rect(l, radians(d)) for d,l in degrees)
    a_phase = phase(a)
    a_degree = degrees(a_phase)


    for deg, len in degrees:
        angle_sum += rect(1, radians(deg))
        tot_len += len

    average_sum = degrees(phase(angle_sum))
    average_sum = average_sum / tot_len

    d = degrees(angle_sum)

    return d



    return degrees(phase(sum(rect(1, radians(d)*l) for d,l in deg)/sum([c[1] for c in deg])))

for angles in [[(350,1000), (10,1)], [(90,1), (180,1), (270,1), (360,1)], [(10,10), (20,1), (30,1)]]:
    print('The mean angle of', angles, 'is:', round(mean_angle(angles)[0], 12), 'degrees')

#for angles in [[(350,2), (10,4)], [(90,2), (180,2), (270,2), (360,2)], [(10,1), (20,2), (30,3)]]:
#    print('The mean angle of', angles, 'is:', round(mean_angle2(angles), 12), 'degrees')

0/0

for xy in [(1,.1),(1,1),(0.1,1),(-0.1,1),(-1,1),(-1,-1),(1,-1)]:
    line0 = LineString([(0,0), xy])
    line1 = LineString([xy, (0,0)])
    for line in (line0,line1):
        x0, y0 = line.coords[0][0], line.coords[0][1]
        x1, y1 = line.coords[1][0], line.coords[1][1]
        delta_y = (y1 - y0)
        delta_x = (x1 - x0)
        angle = math.atan(delta_y / delta_x)
        angle = math.degrees(angle)
        print (x0, y0, x1, y1, angle)

0/0

# Create the triangles
for i in range(250000):
    x = random.random() * 10000.
    y = random.random() * 10000.
    coords = [(x,y),(x+5, y+5),(x,y+10),(x,y)]
    lst_lines.append(LineString(coords))

# Create the bounding boxes
for i in range(10000):
    x = random.random() * 10000.
    y = random.random() * 10000.
    coords = [(x,y),(x+15,y),(x+15,y+15),(x,y+15),(x,y)]
    lst_intersects.append(LineString(coords))

# Create shapely STRtree
tree = STRtree(lst_lines)

# Create RTree
idx = index.Index()
for i, line in enumerate(lst_lines):
    idx.insert(i, line.bounds)
print (time.time())

sec1 = time.time()

# finf the intersection with STRtree
str_tree_nbr = 0
for intersect in lst_intersects:
    str_tree = tree.query(intersect)
    str_tree_nbr += len(str_tree)

sec2 = time.time()
print("Seconds for STRtree =", sec2-sec1)
print ("Str tree number: ", str_tree_nbr)

# Find the intersections with RTree
rtree_nbr = 0
for intersect in lst_intersects:
    rtree = idx.intersection(intersect.bounds)
    rtree_nbr += len(list(rtree))

sec3 = time.time()
print("Seconds for RTree =", sec3-sec2)
print ("Rtree number: ", rtree_nbr)
