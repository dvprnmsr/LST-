# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 15:48:49 2023

@author: purnamas
"""

import ee, folium, geemap
from GEE_EB import GEEImage

# Initialize the Earth Engine API
try:
    ee.Initialize()
    
except Exception as error:
    ee.Authenticate()
    ee.Initialize()
    
#Define your initial date and final date 
initial_date = ee.Date("2019-11-30")
final_date = ee.Date("2019-12-31")

#Define the geometry of the "rhine_basin" feature collection
roi = ee.FeatureCollection("projects/ee-devpurnamasari/assets/rhine_basin")
bbox = roi.geometry()
# metadata = roi.getInfo()
# points = ee.Geometry.Point([5.2913, 52.1326])
# bbox = points.buffer(20000) # var bbox is used as region interest filter

def clip(image):
  return image.clip(bbox)

#Get filtered image collection
collection = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_RAW")
collection = collection.filterDate(initial_date, final_date)
collection = collection.filterBounds(bbox)

# Create a GEEImage instance from the selected bands
sb_collection = GEEImage(collection).from_bands()
preprocessed_images = sb_collection.preprocessing()

#Get clipped filtered image collection (filter: date, region, bands, additional bands)
data_roi = preprocessed_images.map(clip);





###Prepare wflow_input
wflow_input = ee.Image("projects/ee-devpurnamasari/assets/rhine")
bandNames = wflow_input.bandNames();
time_start = data_roi.aggregate_array('system:time_start'); #taking properties of 'system: time_start from GEE collection  
time_end = data_roi.aggregate_array('system:time_end'); #taking properties of 'system: time_start from GEE collection


def addBandToCollection(bandName, initial):
    extracted = wflow_input.select([bandName]) 
    newImage = extracted.rename(['actevap']) #renaming wflow bands from 1 to ... to actevap
    return ee.ImageCollection(initial).merge([newImage])


#new image collection from wflow
wflow_actevap = ee.ImageCollection(bandNames.iterate(addBandToCollection, ee.ImageCollection([])))

# Change system index
# make a list of length of features
idList = ee.List.sequence(0,wflow_actevap.size().subtract(1));

# featureCollection to a List
list = wflow_actevap.toList(wflow_actevap.size());

# set the system:index
actevap = ee.FeatureCollection(idList.map(lambda newSysIndex: \
    ee.Feature(list.get(newSysIndex)) \
        .set({
            'system:index': ee.Number(newSysIndex).format('%3d'),
            'system:time_start': time_start.get(newSysIndex),
            'system:time_end': time_end.get(newSysIndex),
        })
))

filter = ee.Filter.equals(leftField='system:time_start', rightField='system:time_start')

# Define an inner join object
inner_join = ee.Join.inner()

# Apply the inner join to the two collections
joined_collection = ee.ImageCollection(inner_join.apply(actevap, data_roi, filter))

# Map over the joined collection to create a new collection
# with the primary and secondary images concatenated
joined = joined_collection.map(lambda feature: ee.Image.cat(feature.get('primary'), feature.get('secondary')))

LE = joined.first().select('actevap').multiply(1000).multiply(2.5e6).divide(1000).divide(3600).rename('LE') #convert from mm/day to W/m2
joined.first().addBands(LE)


test = GEEImage(joined).postprocessing()





#Get the information needed to display the image in a folium map
b = test.first().select('Ts')
mapid = b.getMapId({'min':0, 'max':310, 'palette': '0000FF, FFFFFF, FF0000'})

#Create a folium map centered on the region
map = folium.Map(location=bbox.centroid().getInfo()['coordinates'][::-1], zoom_start=6)

#Display the image on the folium map
folium.TileLayer(
    tiles=mapid['tile_fetcher'].url_format,
    attr='Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
    overlay=True,
    name='image',
  ).add_to(map)

#Save the folium map to an HTML file
map.save("map.html")

map = geemap.Map(location=bbox.centroid().getInfo()['coordinates'][::-1], zoom_start=6)




