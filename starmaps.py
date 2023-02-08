# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 09:53:54 2022

@author: Ryan.Larson
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from datetime import datetime
import pytz
from tzwhere import tzwhere

from skyfield.api import Star, load, N, W, E, S, wgs84
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from skyfield.positionlib import Astrometric, position_of_radec
from skyfield.vectorlib import VectorFunction, VectorSum

def generate_all_sky(dt_str, lat, lon, filename):
    date_format  = "%m/%d/%Y %H:%M:%S"
    
    # Create datetime object in local timezone
    tz = tzwhere.tzwhere()
    timezone_str = tz.tzNameAt(lat,lon)
    tzone = pytz.timezone(timezone_str)        
    local_dt = tzone.localize(datetime.strptime(dt_str, date_format))
    dt_utc = local_dt.astimezone(pytz.UTC)
    
    ts = load.timescale()
    t = ts.utc(dt_utc.year,dt_utc.month,dt_utc.hour,dt_utc.minute,dt_utc.second)    # This is UTC time as written and must be converted to local time
    if lat >= 0.0:
        lat *= N
    else:
        lat = np.abs(lat) * S
    if lon < 0.0:
        lon = np.abs(lon) * W
    else:
        lon *= E
    geographic = wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon)
    observer = geographic.at(t)
    zenith = observer.from_altaz(alt_degrees=90, az_degrees=0)
    zenith_ra, zenith_dec, dist = zenith.radec()
    
    eph = load('de421.bsp')
    sun = eph['sun']
    earth = eph['earth']
    viewloc = earth + wgs84.latlon(lat*N, lon*W)
    
    
    # The Hipparcos mission provides our star catalog.
    
    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)
    
    # And the constellation outlines come from Stellarium.  We make a list
    # of the stars at which each edge stars, and the star at which each edge
    # ends.
    
    url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
            '/skycultures/western_SnT/constellationship.fab')
    
    with load.open(url) as f:
        constellations = stellarium.parse_constellations(f)
    
    edges = [edge for name, edges in constellations for edge in edges]
    edges_star1 = [star1 for star1, star2 in edges]
    edges_star2 = [star2 for star1, star2 in edges]
    
    projection = build_stereographic_projection(zenith)
    field_of_view_degrees = 180.0
    limiting_magnitude = 6.0
    
    # Now that we have constructed our projection, compute the x and y
    # coordinates that each star and the comet will have on the plot.
    
    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    stars['x'], stars['y'] = projection(star_positions)
    
    # Create a True/False mask marking the stars bright enough to be
    # included in our plot.  And go ahead and compute how large their
    # markers will be on the plot.
    
    bright_stars = (stars.magnitude <= limiting_magnitude)
    magnitude = stars['magnitude'][bright_stars]
    marker_size = (1.0 + limiting_magnitude - magnitude) ** 2.0
    
    # The constellation lines will each begin at the x,y of one star and end
    # at the x,y of another.  We have to "rollaxis" the resulting coordinate
    # array into the shape that matplotlib expects.
    
    xy1 = stars[['x', 'y']].loc[edges_star1].values
    xy2 = stars[['x', 'y']].loc[edges_star2].values
    lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)
    
    # Time to build the figure!
    
    inches = 12
    fig, ax = plt.subplots(figsize=[inches, inches], dpi=300)
    
    # Draw the constellation lines.
    
    ax.add_collection(LineCollection(lines_xy, colors='b'))
    
    # Draw the stars.
    
    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
                s=marker_size, color='k')
    
    # Finally, title the plot and set some final parameters.
    angle = np.pi - field_of_view_degrees / 360.0 * np.pi
    limit = np.sin(angle) / (1.0 - np.cos(angle))
    
    # Plot a circle with the field of view limit
    circle = plt.Circle((0,0), limit, color='k', fill=False)
    ax.add_patch(circle)
    
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_aspect(1.0)
    ax.set_title("Star chart - time: {}".format(t))
    
    # Save.
    fig.savefig(filename, bbox_inches='tight')
    

if __name__ == "__main__":
    dt_str  = "08/26/2015 22:00:00"
    lat = 40.0196
    lon = -111.7495
    
    filename = "test.png"
    
    generate_all_sky(dt_str, lat, lon, filename)