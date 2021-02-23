from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingException,
                       QgsProcessingParameterDistance,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterVectorDestination,
                       QgsWkbTypes,
                       QgsTopologyPreservingSimplifier,
                       QgsGeometry,
                       QgsFeatureSink,
                       QgsProcessingParameterFeatureSink,
                       QgsProcessingParameterBoolean)
import os
import inspect
import processing
from qgis.PyQt.QtGui import QIcon

class SimplifyAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer,
    creates some new layers and returns some results.
    """

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        # Must return a new copy of your algorithm.
        return SimplifyAlgorithm()

    def name(self):
        """
        Returns the unique algorithm name.
        """
        return 'simplify'

    def displayName(self):
        """
        Returns the translated algorithm name.
        """
        return 'Simplify'

    def group(self):
        """
        Returns the name of the group this algorithm belongs to.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs
        to.
        """
        return ''

    def shortHelpString(self):
        """
        Returns a localised short help string for the algorithm.
        """
        help_str = """<b>Simplify</b> \
    Simplify is a geospatial simplification tool for lines and polygons. Simplify implements \
    QGIS's QgsTopologyPreservingSimplifier tool. For line and polygon simplification \
    that tool implements an algorithm similar to the Douglas Peucker algorithm. The implementation \
    preserves the topology within one vector feature but not between vector features. \
    There is also a known bug where the algorithm may create invalid topologies \
    if there are components which are small relative to the tolerance value. In particular, if a \
    small interior hole is very close to an edge, the resulting simplification may result in the hole being \
    moved outside the polygon. This algorithm will detect these situations where one or more rings (interior \
    parts) fall outside the polygon after being simplified and make the polygon invalid. \
    The algoritm will remove (delete) these ring(s) so the feature remains valid after simplification.

    <b>Usage</b>
    
    <u>Input layer</u>: The Line String or Polygon layer to simplify
    
    <u>Tolerance</u>: The tolerance in ground unit used to simplify the line
    
    <u>Simplified</u>: The simplified Line String or Polygon Layer

    For more information: https:...
    """
        return self.tr(help_str)

    def icon(self):
        cmd_folder = os.path.split(inspect.getfile(inspect.currentframe()))[0]
        icon = QIcon(os.path.join(os.path.join(cmd_folder, 'logo.png')))
        return icon
    
    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and outputs of the algorithm.
        """
        # 'INPUT' is the recommended name for the main input
        # parameter.
        self.addParameter(QgsProcessingParameterFeatureSource('INPUT',
                                                              self.tr('Input vector layer'),
                                                              types=[QgsProcessing.TypeVectorAnyGeometry]))

        self.addParameter(QgsProcessingParameterDistance('TOLERANCE',
                                                         self.tr('Tolerance'),
                                                         defaultValue = 1.0,
                                                         # Make distance units match the INPUT layer units:
                                                         parentParameterName='INPUT'))

        self.addParameter(QgsProcessingParameterBoolean('VERBOSE',
                                                        self.tr('Verbose'),
                                                        defaultValue=1.0))

        self.addParameter(QgsProcessingParameterFeatureSink('OUTPUT',
                                                            self.tr('Simplified')))





#from qgis.processing import alg
#from qgis.core import QgsFeature, QgsFeatureSink, QgsTopologyPreservingSimplifier, QgsProcessingException,  QgsWkbTypes, QgsGeometry

#@alg(name="topological_simplifier", label=alg.tr("Topological simplifier"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
#@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
#@alg.input(type=alg.DISTANCE, name="TOLERANCE", label="Tolerance", default=1.0)
#@alg.input(type=alg.BOOL, name="VERBOSE", label="Verbose mode", default=False)
#@alg.input(type=alg.SINK, name="OUTPUT", label="Output layer")

    def processAlgorithm(self, parameters, context, feedback):
        """
        <b>Topological simplifier</b>
        TopoSim is a geospatial simplification tool for lines and polygons. TopoSim implements \
        QGIS's QgsTopologyPreservingSimplifier tool. For line and polygon simplification \
        that tool implements an algorithm similar to the Douglas Peucker algorithm. The implementation \
        preserves the topology within one feature but not between features of the same layer or from \
        different layers. There is also a known bug where the algorithm may create invalid topologies \
        if there are components which are small relative to the tolerance value. In particular, if a \
        small interior hole is very close to an edge, simplification may result in the hole being moved \
        outside the polygon. Toposim will detect these situations where one or more rings (interior \
        parts) fall outside the polygon after being simplified and make the polygon invalid. \
        The algoritm will remove (delete) these ring(s) so the feature remains valid after simplification.

        <b>Usage</b>
        <u>Input</u>: A Line String or Polygon layer
        <u>Tolerance</u>: The tolerance in ground unit used to simplify the line

        For more information: https://github.com/Dan-Eli/GeoSim
        """
        source = self.parameterAsVectorLayer(parameters, "INPUT", context)
        tolerance = self.parameterAsDouble(parameters,"TOLERANCE", context)
        verbose = self.parameterAsBool(parameters, "VERBOSE", context)

        if source is None:
            raise QgsProcessingException(self.invalidSourceError(parameters, "INPUT"))

        (sink, dest_id) = self.parameterAsSink(parameters, "OUTPUT", context,
                                                   source.fields(),
                                                   source.wkbType(),
                                                   source.sourceCrs() )

        # Execute MultiToSinglePart processing
        feedback.pushInfo("Start normalizing input layer")
        params = {'INPUT': source,
                  'OUTPUT': 'memory:'}
        result_ms = processing.run("native:multiparttosingleparts", params, feedback=feedback)
        vector_layer_in = result_ms['OUTPUT']
        feedback.pushInfo("End normalizing input layer")

        # Validate input vector layer in type
        if vector_layer_in.wkbType() not in [QgsWkbTypes.Polygon, QgsWkbTypes.LineString]:
            raise QgsProcessingException("Can only process: LineString and Polygon type layer")

        # Validate sink
        if sink is None:
            raise QgsProcessingException(self.invalidSinkError(parameters, "OUTPUT"))

        nbr_feature = 0
        nbr_hole_del = 0
        nbr_feat_inv_correct = 0
        nbr_feat_inv_uncorrect = 0
        total = 100.0 / vector_layer_in.featureCount() if vector_layer_in.featureCount() else 0
        simplifier = QgsTopologyPreservingSimplifier(tolerance)
        features = vector_layer_in.getFeatures()
        print ("Nbr features: ", vector_layer_in.featureCount())
        for i, feature in enumerate(features):
            geom = feature.geometry()
            s_geom = simplifier.simplify(geom)
            if s_geom.isGeosValid() and s_geom.isSimple():
                # Simplification OK
                pass
            else:
                if s_geom.wkbType() == QgsWkbTypes.Polygon:
                    # Process the invalid geometry Polygon
                    s_geom_parts = s_geom.asPolygon() # Extract the outer and inner rings in a list
                    # Recreate independant Polygon for each part
                    s_geom_pols = [QgsGeometry.fromPolygonXY([s_geom_part]) for s_geom_part in s_geom_parts ]
                    s_geom_pols.sort(key=polygon_area())  # Sort polygon by ascending area size
                    s_geom_outer = s_geom_pols.pop()  # extract the outer ring
                    while s_geom_pols:
                        # Extract each inner ring and test if located inside the polygon in construction
                        s_geom_inner = s_geom_pols.pop()
                        if s_geom_inner.within(s_geom_outer):
                            s_geom_outer.addRing(s_geom_inner.asPolygon()[0])  # Add a new inner ring
                        else:
                            nbr_hole_del += 1
                            if verbose:
                                inner_hole = s_geom_inner.asPolygon()
                                inner_hole_xy = inner_hole[0][0]  # Extract the first coordinate of the ring
                                xy = str(inner_hole_xy.x()) + ", " + str(inner_hole_xy.y())
                                feedback.pushInfo("Inner hole deleted: {0}".format(xy))

                    if s_geom.isGeosValid() and s_geom.isSimple():
                        # Feature corrected... OK
                        nbr_feat_inv_correct += 1
                    else:
                        # Should not happen
                        nbr_feat_inv_uncorrect += 1
                else:
                    # Should not happen...
                    nbr_feat_inv_uncorrect += 1

            # Add the feature in the sink
            feature.setGeometry(s_geom)
            sink.addFeature(feature, QgsFeatureSink.FastInsert)
            nbr_feature += 1

            if feedback.isCanceled():
                break
            feedback.setProgress(int(i * total))

        # Push some output statistics
        feedback.pushInfo("Number of features simplified: {0}".format(nbr_feature))
        feedback.pushInfo("Number of invalid features corrected: {0}".format(nbr_feat_inv_correct))
        feedback.pushInfo("Number of holes deleted: {0}".format(nbr_hole_del))
        feedback.pushInfo("Number of features left invalid: {0}".format(nbr_feat_inv_uncorrect))

        return {"OUTPUT": dest_id}

def polygon_area(pol):
    return pol.area()
