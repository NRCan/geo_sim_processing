v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori2
v.info map=hydro_pol_ori2@PERMANENT
v.out.ogr input=hydro_pol_ori2@PERMANENT output=C:\Users\berge\PycharmProjects\geobato\data\test\ori\coco format=GPKG
v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori3 thin=1
v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori4 thin=20
v.import input=D:\OneDrive\Personnel\Daniel\QGIS\Kingston\Kingston_dp_Road_To_Skel.gpkg layer=Road output=Road
v.voronoi -s input=Road@PERMANENT output=Road
v.voronoi -s input=Road@PERMANENT output=Road_Linear
v.voronoi -s input=Road@PERMANENT output=Road_Linear1
v.voronoi -s -t input=Road@PERMANENT output=Road_Linear1
v.voronoi -s input=Road@PERMANENT output=Road_Linear
v.voronoi -s input=Road@PERMANENT output=Road_Linear
g.remove type=vector name=Road@PERMANENT
g.remove -f type=vector name=Road@PERMANENT
v.import input=D:\OneDrive\Personnel\Daniel\QGIS\Kingston\Kingston_dp_Road_To_Skel.gpkg layer=Road output=Road
v.voronoi -s input=Road@PERMANENT output=RoadL
