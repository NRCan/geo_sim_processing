import numpy, sys
from osgeo import gdal
from osgeo.gdalconst import *

# register all of the GDAL drivers
gdal.AllRegister()

# open the image
inDs = gdal.Open("mono_lake_rgb.tif")
if inDs is None:
  print ("couldn't open input dataset")
  sys.exit(1)
else:
  print ("opening was successful!")
cols = inDs.RasterXSize
rows = inDs.RasterYSize
bands =  inDs.RasterCount
driver = inDs.GetDriver()
#driver.Create("newfile.tif",cols,rows,3)
#outDs = gdal.Open("newfile.tif")

#if outDs is None:
#  print ("failure to create new file")
#  sys.exit(1)


# Create a new raster data source
outDs = driver.Create("newfile.tif", cols, rows, 3, gdal.GDT_UInt16)

# Write metadata
outDs.SetGeoTransform(inDs.GetGeoTransform())
outDs.SetProjection(inDs.GetProjection())

# Write raster data sets
for i in range(3):
    outBand = outDs.GetRasterBand(i + 1)
    data = inDs.GetRasterBand(i+1).ReadAsArray()
    print (data[0:5,0:5])
    outBand.WriteArray(data)

# Close raster file
print ("Fin")
outDs = None

#outBand1 = outDs.GetRasterBand(1)
#outBand2 = outDs.GetRasterBand(2)
#outBand3 = outDs.GetRasterBand(3)
#data1 = inDs.GetRasterBand(1).ReadAsArray()
#data2 = inDs.GetRasterBand(2).ReadAsArray()
#data3 = inDs.GetRasterBand(3).ReadAsArray()

#outBand1.WriteArray(data1,0,0)
#outBand2.WriteArray(data2,0,0)
#outBand3.WriteArray(data3,0,0)
#outDs.SetProjection(inDs.GetProjection())
#outDs.SetGeoTransform(inDs.GetGeoTransform())
#outDs = None


outDs = gdal.Open("newfile.tif")
print ("after reopening")
print (outDs.GetRasterBand(1).ReadAsArray())
print (outDs.GetRasterBand(2).ReadAsArray())
print (outDs.GetRasterBand(3).ReadAsArray())