from shapely.ops import linemerge
from shapely.geometry import LineString
l0 = LineString(((0,0),(1,1)))
l0i = LineString(((1,1),(0,0)))
l1 = LineString(((1,1),(2,2)))
l1i = LineString(((2,2),(1,1)))
l2 = LineString(((2,2),(3,3)))
l2i = LineString(((3,3),(2,2)))
lm0 = linemerge([l0,l1,l2])
lm1i = linemerge([l0i,l1i,l2i])
lm2i = linemerge([l0,l2,l1])
lm3i = linemerge([l2,l1,l0])
lm4i = linemerge([l2,l1i,l0])
lm5i = linemerge([l2,l1i,l0])
print
