# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:11:26 2023

@author: Ryan.Larson
"""

from datetime import datetime
from geopy import Nominatim
from tzwhere import tzwhere
from pytz import timezone, utc

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos, stellarium
from skyfield.projections import build_stereographic_projection

# load celestial data

# de421 shows position of earth and sun in space
eph = load('de421.bsp')

# hipparcos dataset contains star location data
with load.open(hipparcos.URL) as f:
    stars = hipparcos.load_dataframe(f)
location = 'Payson, UT'
# location = 'Times Square, New York, NY'
when = '2015-08-26 22:00'
# get latitude and longitude of our location 
locator = Nominatim(user_agent='myGeocoder')
location = locator.geocode(location)
lat, long = location.latitude, location.longitude
# convert date string into datetime object
dt = datetime.strptime(when, '%Y-%m-%d %H:%M')

# define datetime and convert to utc based on our timezone
timezone_str = tzwhere.tzwhere().tzNameAt(lat, long)
local = timezone(timezone_str)

# get UTC from local timezone and datetime
local_dt = local.localize(dt, is_dst=None)
utc_dt = local_dt.astimezone(utc)

# find location of earth and sun and set the observer position
sun = eph['sun']
earth = eph['earth']

# define observation time from our UTC datetime
ts = load.timescale()
t = ts.from_datetime(utc_dt)

# define an observer using the world geodetic system data
observer = wgs84.latlon(latitude_degrees=lat, longitude_degrees=long).at(t)

# define the position in the sky where we will be looking
position = observer.from_altaz(alt_degrees=90, az_degrees=0)
# center the observation point in the middle of the sky
ra, dec, distance = observer.radec()
center_object = Star(ra=ra, dec=dec)

# Get constellation lines from Stellarium
url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
        '/skycultures/western_SnT/constellationship.fab')

with load.open(url) as f:
    constellations = stellarium.parse_constellations(f)

edges = [edge for name, edges in constellations for edge in edges]
edges_star1 = [star1 for star1, star2 in edges]
edges_star2 = [star2 for star1, star2 in edges]

# find where our center object is relative to earth and build a projection with 180 degree view
center = earth.at(t).observe(center_object)
projection = build_stereographic_projection(center)
field_of_view_degrees = 180.0

# calculate star positions and project them onto a plain space
star_positions = earth.at(t).observe(Star.from_dataframe(stars))
stars['x'], stars['y'] = projection(star_positions)

chart_size = 10
max_star_size = 200
limiting_magnitude = 10

bright_stars = (stars.magnitude <= limiting_magnitude)
magnitude = stars['magnitude'][bright_stars]
fig, ax = plt.subplots(figsize=(chart_size, chart_size), dpi=300)
    
border = plt.Circle((0, 0), 1, color='navy', fill=True)
ax.add_patch(border)

marker_size = max_star_size * 10 ** (magnitude / -2.5)

# The constellation lines will each begin at the x,y of one star and end
# at the x,y of another.  We have to "rollaxis" the resulting coordinate
# array into the shape that matplotlib expects.

xy1 = stars[['x', 'y']].loc[edges_star1].values
xy2 = stars[['x', 'y']].loc[edges_star2].values
lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)

# Draw the constellation lines.
ax.add_collection(LineCollection(lines_xy, colors='white', linewidths=(0.5)))

ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
           s=marker_size, color='white', marker='.', linewidths=0, 
           zorder=2)

horizon = Circle((0, 0), radius=1, transform=ax.transData)
for col in ax.collections:
    col.set_clip_path(horizon)


# other settings
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
plt.axis('off')

plt.show()