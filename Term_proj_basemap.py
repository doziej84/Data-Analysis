#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 17:09:12 2020

@author: user
"""
# T plot the locations of the earthquakes (Solomon Islands) and Depth on Basemap for 
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv

# setting parameters for title and axes
font = {'family' : 'verdana',
        'size'   : 16}
matplotlib.rc('font', **font)

# Grabbing the .csv data
lats,lons,depth, magnitude = [],[],[],[]
temp_dat = []
with open('solomon.csv') as csvfile:
    reader = csv.DictReader(csvfile,delimiter=',')
    for data in reader:
        #if float(data['UTC-6'])>-5. or float(data['UTC-6'])<-8. or float(data['depth'])<0.0:
            #continue
        magnitude.append(float(data['mag']))
        lats.append(float(data['latitude']))
        lons.append(float(data['longitude']))
        depth.append(float(data['depth']))
        
# How much to zoom from coordinates (in degrees)
zoom_scale = 0

# Setup the bounding box for the zoom and bounds of the map
bbox = [np.min(lats)-zoom_scale,np.max(lats)+zoom_scale,\
        np.min(lons)-zoom_scale,np.max(lons)+zoom_scale]

fig, ax = plt.subplots(figsize=(12,7))
plt.title("Earthquake, depth and Subduction zone")
# Define the projection, scale, the corners of the map, and the resolution.
m = Basemap(projection='merc',llcrnrlat=bbox[0],urcrnrlat=bbox[1],\
            llcrnrlon=bbox[2],urcrnrlon=bbox[3],lat_ts=10,resolution='f')

# Draw coastlines and fill continents and water with color
m.drawcoastlines()
m.fillcontinents(color='peru',lake_color='lightblue')

# draw parallels, meridians, and color boundaries
m.drawparallels(np.arange(bbox[0],bbox[1],(bbox[1]-bbox[0])/5),labels=[1,0,0,0])
m.drawmeridians(np.arange(bbox[2],bbox[3],(bbox[3]-bbox[2])/5),labels=[0,0,0,1],rotation=15)
m.drawmapboundary(fill_color='lightblue')

# format colors for elevation range
depth_min = np.min(depth)
depth_max = np.max(depth)
cmap = plt.get_cmap('gist_earth')
normalize = matplotlib.colors.Normalize(vmin=depth_min, vmax=depth_max)

# plot elevations with different colors using the numpy interpolation mapping tool
# the range [10,400] can be changed to create different colors and ranges
for ii in range(0,len(depth)):
    x,y = m(lons[ii],lats[ii])
    color_interp = np.interp(depth[ii],[depth_min,depth_max],[10,400])
    plt.plot(x,y,marker='o',markersize=4,color=cmap(int(color_interp)))

# format the colorbar 
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,norm=normalize,label='Depth')

# save the figure and show it
#plt.savefig('depth_subductionzone.png', format='png', dpi=500,transparent=True)
plt.show()