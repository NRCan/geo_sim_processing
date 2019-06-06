from shapely.geometry import LineString

class LineStringSc(LineString):

    def __init__(self, coords, fast_access=True):
        super().__init__(coords)
        self.fast_access = fast_access
        if self.fast_access:
            self.__lst_coords = list(super().coords)

    @property
    def coords(self):
        if self.fast_access:
            return self.__lst_coords
        else:
            return super().coords

    @coords.setter
    def coords(self, coords):
        print ("Need to update the spatial container...")
        LineString.coords.__set__(self, coords)
        if self.fast_access:
            self.__lst_coords = list(super().coords)


xy1 = [(10,10),(20,20)]
xy2 = [(11,11),(21,21)]
line_B = Line(xy1)
line_C = Line(xy2)

print  ("Len B: ", len(line_B.coords))
print  ("Len C: ", len(line_C.coords))
print ("List B: ", line_B.coords[0])
print ("List C: ", line_C.coords[0])

line_B.coords =  [(10,10),(20,20),(30,30)]
line_C.coords = [(1,1),(2,2),(3,3),(4,4)]
print (len(line_B.coords), list(line_B.coords))
print (len(line_C.coords), list(line_C.coords))
#t = a.xy[0:2]
#print (t)
#print (len(xy))
print ("Fin...")
