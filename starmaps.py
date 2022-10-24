# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 09:53:54 2022

@author: Ryan.Larson
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection

from skyfield.api import Star, load, N, W, wgs84, Topos
from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN
from skyfield.data import hipparcos, mpc, stellarium
from skyfield.projections import build_stereographic_projection
from skyfield.positionlib import Astrometric, position_of_radec
from skyfield.vectorlib import VectorFunction, VectorSum

ts = load.timescale()
t = ts.utc(2022,10,24,7)    # This is UTC time as written and must be converted to local time
lat = 41.0
lon = 111.0
geographic = wgs84.latlon(latitude_degrees=lat*N, longitude_degrees=lon*W)
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


# We will center the chart on the comet's middle position.
# comet = sun + mpc.comet_orbit(row, ts, GM_SUN)
# barycentric = earth.at(t)
# center_position = position_of_radec(zenith_ra.hours, zenith_dec.degrees)
# # center_position = VectorSum(center_position.position.au)
# astrometric = Astrometric(center_position.position.au) # Not sure what is needed as an argument here
# center = barycentric.observe(astrometric)
# center = earth.at(t).observe(astrometric)
projection = build_stereographic_projection(zenith)
# projection = build_stereographic_projection(center)
field_of_view_degrees = 180.0
limiting_magnitude = 7.0

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

fig, ax = plt.subplots(figsize=[15, 15], dpi=300)

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

fig.savefig('star_chart.png', bbox_inches='tight')


# # The comet is plotted on several dates `t_comet`.  But the stars only
# # need to be drawn once, so we take the middle comet date as the single
# # time `t` we use for everything else.

# ts = load.timescale()
# t_comet = ts.utc(2020, 7, range(17, 27))
# t = t_comet[len(t_comet) // 2]  # middle date

# # An ephemeris from the JPL provides Sun and Earth positions.

# eph = load('de421.bsp')
# sun = eph['sun']
# earth = eph['earth']
# lat = 41.0
# lon = 111.0
# viewloc = earth + wgs84.latlon(lat*N, lon*W)




# # The Minor Planet Center data file provides the comet orbit.

# with load.open(mpc.COMET_URL) as f:
#     comets = mpc.load_comets_dataframe(f)

# comets = (comets.sort_values('reference')
#           .groupby('designation', as_index=False).last()
#           .set_index('designation', drop=False))

# row = comets.loc['C/2020 F3 (NEOWISE)']
# comet = sun + mpc.comet_orbit(row, ts, GM_SUN)

# # The Hipparcos mission provides our star catalog.

# with load.open(hipparcos.URL) as f:
#     stars = hipparcos.load_dataframe(f)

# # And the constellation outlines come from Stellarium.  We make a list
# # of the stars at which each edge stars, and the star at which each edge
# # ends.

# url = ('https://raw.githubusercontent.com/Stellarium/stellarium/master'
#         '/skycultures/western_SnT/constellationship.fab')

# with load.open(url) as f:
#     constellations = stellarium.parse_constellations(f)

# edges = [edge for name, edges in constellations for edge in edges]
# edges_star1 = [star1 for star1, star2 in edges]
# edges_star2 = [star2 for star1, star2 in edges]

# # We will center the chart on the comet's middle position.

# center = earth.at(t).observe(comet)
# projection = build_stereographic_projection(center)
# field_of_view_degrees = 45.0
# limiting_magnitude = 7.0

# # Now that we have constructed our projection, compute the x and y
# # coordinates that each star and the comet will have on the plot.

# star_positions = earth.at(t).observe(Star.from_dataframe(stars))
# stars['x'], stars['y'] = projection(star_positions)

# comet_x, comet_y = projection(earth.at(t_comet).observe(comet))

# # Create a True/False mask marking the stars bright enough to be
# # included in our plot.  And go ahead and compute how large their
# # markers will be on the plot.

# bright_stars = (stars.magnitude <= limiting_magnitude)
# magnitude = stars['magnitude'][bright_stars]
# marker_size = (0.5 + limiting_magnitude - magnitude) ** 2.0

# # The constellation lines will each begin at the x,y of one star and end
# # at the x,y of another.  We have to "rollaxis" the resulting coordinate
# # array into the shape that matplotlib expects.

# xy1 = stars[['x', 'y']].loc[edges_star1].values
# xy2 = stars[['x', 'y']].loc[edges_star2].values
# lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)

# # Time to build the figure!

# fig, ax = plt.subplots(figsize=[9, 9])

# # Draw the constellation lines.

# ax.add_collection(LineCollection(lines_xy, colors='#00f2'))

# # Draw the stars.

# ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars],
#             s=marker_size, color='k')

# # Draw the comet positions, and label them with dates.

# comet_color = '#f00'
# offset = 0.002

# ax.plot(comet_x, comet_y, '+', c=comet_color, zorder=3)

# for xi, yi, tstr in zip(comet_x, comet_y, t_comet.utc_strftime('%m/%d')):
#     tstr = tstr.lstrip('0')
#     text = ax.text(xi + offset, yi - offset, tstr, color=comet_color,
#                     ha='left', va='top', fontsize=9, weight='bold', zorder=-1)
#     text.set_alpha(0.5)

# # Finally, title the plot and set some final parameters.

# angle = np.pi - field_of_view_degrees / 360.0 * np.pi
# limit = np.sin(angle) / (1.0 - np.cos(angle))

# ax.set_xlim(-limit, limit)
# ax.set_ylim(-limit, limit)
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# ax.set_aspect(1.0)
# ax.set_title('Comet NEOWISE {} through {}'.format(
#     t_comet[0].utc_strftime('%Y %B %d'),
#     t_comet[-1].utc_strftime('%Y %B %d'),
# ))

# # Save.

# fig.savefig('neowise-finder-chart.png', bbox_inches='tight')