v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori2
v.info map=hydro_pol_ori2@PERMANENT
v.out.ogr input=hydro_pol_ori2@PERMANENT output=C:\Users\berge\PycharmProjects\geobato\data\test\ori\coco format=GPKG
v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori3 thin=1
v.voronoi -s input=hydro_pol_ori@PERMANENT output=hydro_pol_ori4 thin=20
