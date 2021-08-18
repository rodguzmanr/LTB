#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Variables needed by the LTB_simple_model.py script to compute
   the main features of an indealized Light Tower of Babel.

This file contains all the atmospheric and LTB-related variables
needed in order to compute all the other variables (module size,
mass, ...) describing the LTB with such hypotheses/configuration.

May 2021, v0 by Rodrigo Guzman

"""

import numpy as np


### Physical variables to be set by the user (I) in order ###
### to perform calculations (II) related to the main      ###
### features of LTB.                                      ###
### All units are in International System of Units (SI)   ###


#############################################################
#####  I - Physical variables values given by the user  #####
#############################################################

#############  a) Static atmosphere hypotheses  #############

# Vertical and constant temperature gradient [K/m]
temp_grad = -6.5e-3
# Upper altitude limit temperature
min_temp = 200.

# Surface conditions
# Temperature in [K]
surf_temp = 290.
# Pressure in [Pa]
surf_pres = 102000.
# Surface elevation in [m] above the sea level
surf_alt = 0.

# Location
# Latitude where the LTB is located, in [degrees North]
lat = 50.
# Longitude where the LTB is located, in [degrees East]
lon = 100.
# Days of interest to plot the daily shadows at the surface,
# by default the 2020 solstices and the spring equinox
n_days = ['20200321','20200621','20201221']


###########  b) Light Tower of Babel hypotheses  ############

# Part 1, the cubic casing
# All cubic casing have the same shape and dimensions
# Length [m]
x_cube = 20.
# Safety distance between cubic casings and the rest of the structure,
# in percent with respect to length considered [%]
x_cube_margin = 10.
# Safety volume between theoretical maximum volume for the cubic
# casings and the maximum volume accepted, in percent with respect
# to the theoretical maximum volume [%]
max_vol_margin = 10.

# Part 2, the tower quasi-cylindrical modules
# Safety distance between cubic casings and the top of each module,
# in percent with respect to the maximum height of the casing [%]
z_mod_margin = 25.

# Part 3, the tower
# Height of the tower, from the surface elevation
# to the tropopause (lower stratosphere), in [m]
z_tower = 15000.
# Surface density in [kg/m^2], corresponding to a "volume density"
# for the equivalent of a solid angle having à 1 m^2 surface at
# the edge of the tower section.
# This makes the mass of each module proportional to its radius.
dens_surf_mod = 20.

# Part 4, PV panel at the top of the tower
# Length, this dimension always stays horizontal
x_pv = 2000.
# Width, this dimension will go from horizontal to vertical
# on a daily basis once deployed
y_pv = 500.
# Effective conversion rate of solar radiative flux to electricity, in [%]
conv_rate_pv = 20.
# Average hours per day of full exposure of the PV panel to solar flux
hours_day = 11


####################################################################
##  II - Physical variables derived from information given above  ##
####################################################################

##############  b) Light Tower of Babel computation  ###############

# Part 1, the cubic casing
# One casing side surface [m^2]
xy_cube = x_cube**2
# Maximum volume for the casing [m^3]
max_vol_cube = x_cube**3

# Part 2, the tower quasi-cylindrical modules
# Module altitude considering z_mod_margin [m]
z_mod = x_cube*(1+z_mod_margin/100.)

# Part 3, the tower
# Minimum diameter of a module considering x_cube_margin
min_d_tower = x_cube*np.sqrt(2)*(1+x_cube_margin/100.)

# Part 4, PV panel at the top of the tower
# Total surface of the PV panel
surf_pv = x_pv*y_pv
# Number of cubic casings attached the PV panel, x dimension
nb_cube_x_pv = int(x_pv/(x_cube*np.sqrt(2)*(1+x_cube_margin/100.)))
# Number of cubic casings attached to the PV panel, y dimension
nb_cube_y_pv = int(y_pv/(x_cube*(1+x_cube_margin/100.)))
# Total number of cubic casings attached to the PV panel
nb_cube_pv = nb_cube_x_pv*nb_cube_y_pv
# Maximum distance between the center of the PV panel
# and the corners of it
max_r_pv = np.sqrt(x_pv**2+y_pv**2)/2
# Maximum horizontal speed for the corners of the PV panel
max_speed_pv = 2.*np.pi*max_r_pv/(60*60*24)
