from qgis.core import QgsApplication, QgsPolygon, QgsGeometry, QgsLineString, QgsPoint, QgsFeature, QgsVectorLayer
from processing.core.Processing import Processing
import processing

def create_point(coord, ret_geom=True):

    qgs_point = QgsPoint(coord[0], coord[1])
    if ret_geom:
        ret_val = QgsGeometry(qgs_point)
    else:
        ret_val = qgs_point.clone()

    return ret_val

def create_line(coords, ret_geom=True):

    qgs_points = []
    for coord in coords:
        qgs_points.append(create_point(coord, False))

    if ret_geom:
        ret_val = QgsGeometry(QgsLineString(qgs_points))
    else:
        ret_val = QgsLineString(qgs_points).clone()

    return ret_val

def create_polygon(outer, inners):

    outer_line = create_line(outer, False)
    qgs_pol = QgsPolygon()
    qgs_pol.setExteriorRing(outer_line)
    for inner in inners:
        inner_line = create_line(inner, False)
        qgs_pol.addInteriorRing(inner_line)
    qgs_geom = QgsGeometry(qgs_pol)

    return qgs_geom

QgsApplication.setPrefixPath("/usr/bin/qgis", False)

profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)
from qgis._3d import Qgs3DAlgorithms
QgsApplication.processingRegistry().addProvider(Qgs3DAlgorithms(QgsApplication.processingRegistry()))

# Load providers
app.initQgis()
Processing.initialize()
for alg in QgsApplication.processingRegistry().algorithms():
    print("{}:{} --> {}".format(alg.provider().name(), alg.name(), alg.displayName()))



vl = QgsVectorLayer("Polygon", "temporary_polygon", "memory")
pr = vl.dataProvider()



# add a feature
fet = QgsFeature()
fet.setId(1)
geom = create_polygon([(0,0),(0,10),(10,10),(10,0),(0,0)],[])
geom = create_polygon([(0,0),(10,10),(0,10),(10,0),(0,0)],[])
fet.setGeometry(geom)
pr.addFeatures([fet])

# update layer's extent when new features have been added
# because change of extent in provider is not propagated to the layer
vl.updateExtents()

try:
    output = processing.run("3d:tessellate",
                          {'INPUT': vl, 'OUTPUT': 'TEMPORARY_OUTPUT'})
except Exception:
    print ("erreur: ", fet.id())

source = output['OUTPUT']
features = source.getFeatures()
for qgs_feature_in in features:
    print (qgs_feature_in)
    geom = qgs_feature_in.geometry()
    print (geom)