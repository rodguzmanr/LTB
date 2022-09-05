#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions needed by the LTB_simple_model.py script to compute or
   plot the main features of an indealized Light Tower of Babel (LTB).

This file contains all the thermodynamic formulas needed in order to
compute all the atmospheric and LTB-related variables (pressure, mass,
...) describing the LTB with such hypotheses/configuration. It also
contains the functions to plot and save graphically some features of
the LTB.

May 2021, v0 by Rodrigo Guzman

"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pysolar.solar as pysol

import constants as cst


def compute_temp(temp, surf_temp, temp_grad, min_temp, alt, surf_alt):
    """Function to compute the idealized atmosphere temperature profile.

    Parameters
    ----------   
    temp : array (Numpy object)
        Initialized temperature profile of the atmosphere, in [K].

    surf_temp : float
        Surface temperature of the atmosphere, first value of the profile.

    temp_grad : float
        Constant profile temperature gradient, until min_temp is reached.

    min_temp : float
        Absolute minimum temperature of the profile, reached at the
        tropopause and considered as constant upwards.

    alt : array (Numpy object)
        Altitudes of the profile (from the surface to the tropopause), in [m].

    surf_alt : float
        Surface elevation altitude, first value of the altitude profile.

    Returns
    -------
    temp : array (Numpy object)
        Updated temperature profile of the atmosphere, in [K].

    """

    # > < : ces signes ne marchent pas sur mon clavier...
    # Assuming a constant decrease until the tropopause
    temp = surf_temp + (alt-surf_alt)*temp_grad
    # temp is constant from the tropopause upwards
    temp[temp < min_temp] = min_temp 

    return temp


def compute_pres(pres, surf_pres, alt, surf_alt, temp_grad, surf_temp):
    """Function to compute the idealized atmosphere pressure profile.

    Parameters
    ----------   
    pres : array (Numpy object)
        Initialized pressure profile of the atmosphere.

    surf_pres : float
        Surface pressure of the atmosphere, first value of the profile.

    alt : array (Numpy object)
        Altitudes of the profile (from the surface to the tropopause), in [m].

    surf_alt : float
        Surface elevation altitude, first value of the altitude profile.

    temp_grad : float
        Constant profile temperature gradient, until min_temp is reached.

    surf_temp : float
        Surface temperature of the atmosphere, first value of the profile.

    Returns
    -------
    pres : array (Numpy object)
        Updated pressure profile of the atmosphere, in [Pa].

    """

    # Formula assuming hydrostatic equilibrium,
    # constant temp_grad and ideal gas
    pres = surf_pres/(1+(temp_grad*alt)/surf_temp)**(cst.gravity/(cst.R_air*temp_grad))
    
    return pres


def mass_from_ideal_gas(pres, vol, mol_mass, temp):
    """Function to compute the mass of an ideal gas.

    Parameters
    ----------   
    pres : array (Numpy object)
        Pressure of the ideal gas, in [Pa].

    vol : array (Numpy object)
        Volume of the ideal gas, in [m^3].

    mol_mass : float
        Molar mass of the ideal gas.

    temp : array (Numpy object)
        Temperature of the ideal gas, in [K].

    Returns
    -------
    mass : array (Numpy object)
        Total mass of the ideal gas considered.

    """

    # Mass from the ideal gas law
    mass = pres*vol*mol_mass/(cst.R*temp)

    return mass


def module_mass(dens_surf_mod, z_mod, d_mod):
    """Function to compute the mass of a cylindrical module.

    Parameters
    ----------   
    dens_surf_mod : float
        Module surface density used as a "volume density", in [kg/m^2].

    z_mod : float
        Module altitude considering z_mod_margin, in [m].

    d_mod : float
        Module diameter, in [m]

    Returns
    -------
    mass_mod : float
        Total mass of the cylindrical module considered.

    """

    # Mass derived from the module variables
    mass_mod = dens_surf_mod*z_mod*d_mod*np.pi

    return mass_mod


def plot_profiles(pres, alt, nb_cube_mod, temp, d_mod, d_mod_cone, nb_cube_tower):
    """Function to plot the atmospheric and LTB profiles.

    Parameters
    ----------   
    pres : array (Numpy object)
        Pressure of the ideal atmosphere (from the surface to the tropopause), in [Pa].

    alt : array (Numpy object)
        Altitudes of the profile (from the surface to the tropopause), in [m].

    nb_cube_mod : array (Numpy object)
        Number of casings per module (from the surface to the tropopause), no unit.

    temp : array (Numpy object)
        Temperature of the ideal atmosphere (from the surface to the tropopause), in [K].

    d_mod : array (Numpy object)
        Diameter of the LTB modules (from the surface to the tropopause), in [m].

    d_mod_cone : array (Numpy object)
        Modules' diameter for a (upside down) conical LTB (from the surface to the tropopause), no unit.

    nb_cube_tower : float
        Total number of cubic casings within the tower, no unit.

    """

    fig = plt.figure(figsize=(8, 6))
    plt.suptitle('Figure 1: Idealized atmosphere and LTB profiles', fontsize=14)
    
    # Top left subplot 
    plt.subplot(2, 2, 1)
    plt.plot(pres, alt, color='blue')
    plt.title('a) Atmospheric pressure')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Pressure [Pa]')

    # Top right subplot 
    plt.subplot(2, 2, 2)
    plt.plot(nb_cube_mod, alt, color='black')
    plt.title('b) Number of casings per module')
    plt.text(nb_cube_mod[-1]*0.5, alt[0]+1000, 'Total number of\ncasings within the\ntower : '+str(nb_cube_tower), fontsize=10)
    plt.ylabel('Altitude [m]')
    plt.xlabel('Nb casings [no unit]')
    
    # Bottom left subplot
    plt.subplot(2, 2, 3)
    plt.plot(temp, alt, color='red')
    plt.title('c) Atmospheric temperature')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Temperature [K]')

    # Bottom right subplot 
    plt.subplot(2, 2, 4)
    plt.plot(d_mod, alt, color='black', label='cylindrical modules tower')
    plt.plot(d_mod_cone, alt, color='green', label='conical tower')
    plt.title('d) Tower diameter')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Diameter [m]')
    plt.legend(fontsize=9)
    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig1_profiles_LTB_simple_model.png')


def plot_sections(alt, nb_cube_mod, d_mod, x_cube, x_cube_margin):
    """Function to plot the casings in three sections (modules) of the LTB.

    Parameters
    ----------   
    alt : array (Numpy object)
        Altitudes of the profile (from the surface to the tropopause), in [m].

    nb_cube_mod : array (Numpy object)
        Number of casings per module (from the surface to the tropopause), no unit.

    d_mod : array (Numpy object)
        Diameter of the LTB modules (from the surface to the tropopause), in [m].

    x_cube : float
        Side length of all cubic casings, in [m].

    x_cube_margin : float
        Safety distance between cubic casings and the rest of the structure, in [%].

    """

    # Three sections will be plotted: top, middle and bottom of the LTB
    top, mid, bot = -1, int(len(alt)/2), 0
    # Information needed for the figure
    sub_letter = ["a", "b", "c"]
    sub_position = ["Top", "Middle", "Bottom"]
    # The features of the three sections are selected
    radius_mod = [d_mod[top]/2, d_mod[mid]/2, d_mod[bot]/2]
    cubes_per_mod = [nb_cube_mod[top], nb_cube_mod[mid], nb_cube_mod[bot]]
    alt_mod = [int(alt[top]), int(alt[mid]), int(alt[bot])]
    
    # Three subplots figure
    fig, axes = plt.subplots(3, 1, figsize=(5, 14))
    plt.suptitle('Figure 2: Three sections (modules) of the LTB\n', fontsize=14)
    
    # For each module plotted
    for k in np.arange(3):
    
        # Drawing the circle representing the module
    	module = plt.Circle((0, 0), radius_mod[k], fill=False)
    	axes[k].add_patch(module)
    	# Distance between the center of the module and the center of the
    	# most distant cubic casing
    	c1 = radius_mod[k] - x_cube*np.sqrt(2)*(1+x_cube_margin/100)/2
    	# Compute initial position within the circle to draw all the casings    
    	x_c1_ll = -(np.sqrt(2)*c1/2 + x_cube/2)
    	y_c1_ll = -(np.sqrt(2)*c1/2 + x_cube/2)
    	# From this position, we draw all the casings within the module
    	for i in np.arange(np.sqrt(cubes_per_mod[k])):
    	    for j in np.arange(np.sqrt(cubes_per_mod[k])):
            	casing = plt.Rectangle((x_c1_ll+i*x_cube*(1+x_cube_margin/100), y_c1_ll+j*x_cube*(1+x_cube_margin/100)), x_cube, x_cube, facecolor="grey", alpha=0.5)
            	axes[k].add_patch(casing)
        # Subplot text
    	axes[k].set_title(sub_letter[k]+') '+sub_position[k]+' module, altitude = '+str(alt_mod[k])+' m')
    	axes[k].set_xlim((-radius_mod[0], radius_mod[0]))
    	axes[k].set_ylim((-radius_mod[0], radius_mod[0]))
    	axes[k].set_ylabel('Distance [m]')
    	axes[k].set_xlabel('Distance [m]')

    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig2_sections_LTB_simple_model.png')


def plot_shadows(lat, lon, n_days, x_pv, y_pv, z_tower):
    """Function to plot the surface shadows caused by the PV panel.

    Parameters
    ----------   
    lat : float
        Latitude where the LTB is located, in [degrees North].

    lon : float
        Longitude where the LTB is located, in [degrees East].

    n_days : list (Python object)
        List of days for which we wish to plot the surface shadows at every hour of the day, expected string format for each day [yyyymmdd].

    x_pv : float
        Length of the PV panel at the top of the tower, this dimension always stays horizontal, in [m].

    y_pv : float
        Width of the PV panel at the top of the tower, this dimension will go from horizontal to vertical, in [m].

    z_tower : float
        Height of the tower, from the surface elevation to the tropopause (lower stratosphere), in [m].

    """
    
    # Defining coordinate tables for the days to be plotted
    x_box_centers = np.ones((len(n_days),24))*np.nan
    y_box_centers = np.ones((len(n_days),24))*np.nan
    x_box_lls = np.ones((len(n_days),24))*np.nan
    y_box_lls = np.ones((len(n_days),24))*np.nan
    x_shadows = np.ones((len(n_days),24))*np.nan
    y_shadows = np.ones((len(n_days),24))*np.nan
    saas = np.ones((len(n_days),24))*np.nan
    # Maximum value for plots
    max_val = 999999
    # Defining variable to switch from one hemisphere to the other
    lat_sign = 1.
    if lat < 0:
        lat_sign = -1.
    
    # Creating figure
    fig, ax = plt.subplots(figsize=(6, 6))
    plt.suptitle('Figure 3: Surface shadows in the surrounding LTB area\ncaused by the PV panel array in clear-sky conditions', fontsize=14)
    
    # Loop over the number of days in n_days
    for i in np.arange(len(n_days)):
        # Initializing daily counter
        k = 0
        # Loop over the 24 hours of the day      
        for j in np.arange(24):
            # Defining the datetime object at local hour j at the
            # geographical location corresponding to the given longitude
            date = dt.datetime(int(n_days[i][0:4]), int(n_days[i][4:6]), int(n_days[i][6:8]), j, tzinfo=dt.timezone(offset=dt.timedelta(hours=int(lon/15))))
            # Getting the Solar Zenithal Angle
            sza = pysol.get_altitude(lat, lon, date)
            # Computation performed if the sun is above the horizon
            if sza > 0:
                # Getting the Solar Azimuth Angle
                saa = pysol.get_azimuth(lat, lon, date)
                # The x component of the surface shadow is x_pv
                # because there is no angle involved, the PV panel
                # is always perpendicular to the solar flux
                x_shadow = x_pv
                # Computing the y component of the surface shadow
                # which depends on the altitude of the sun in the sky
                y_shadow = y_pv/np.sin((2*np.pi/360.)*sza)
                # Computing the distance between the center of the
                # surface shadow and the bottom of the tower
                l_shadow = z_tower/np.tan((2*np.pi/360.)*sza)

                # Center of the box coordinates
                x_box_center = l_shadow*(-np.sin((2*np.pi/360.)*saa))
                y_box_center = l_shadow*(-np.cos((2*np.pi/360.)*saa))
                
                # Lower left corner of the box coordinates
                x_box_ll = x_box_center + (x_shadow/2.)*(np.cos((2*np.pi/360.)*saa)) + (y_shadow/2.)*(np.sin((2*np.pi/360.)*saa))
                y_box_ll = y_box_center + (x_shadow/2.)*(-np.sin((2*np.pi/360.)*saa)) + (y_shadow/2.)*(np.cos((2*np.pi/360.)*saa))

                # We save the box coordinates to plot the shadows
                # created by the PV panel array at the ground
                if k == 0:
                    x_box_centers[i, :] = x_box_center
                    y_box_centers[i, :] = y_box_center
                    x_box_lls[i, :] = x_box_ll
                    y_box_lls[i, :] = y_box_ll
                    x_shadows[i, :] = x_shadow
                    y_shadows[i, :] = y_shadow
                    saas[i, :] = saa
                    k = 1
                else:
                    x_box_centers[i, j:] = x_box_center
                    y_box_centers[i, j:] = y_box_center
                    x_box_lls[i, j:] = x_box_ll
                    y_box_lls[i, j:] = y_box_ll
                    x_shadows[i, j:] = x_shadow
                    y_shadows[i, j:] = y_shadow
                    saas[i, j:] = saa

    # Plotting each day's surface shadow trajectory
    for i in np.arange(len(n_days)):
        ax.plot(x_box_centers[i, :], y_box_centers[i, :], label=n_days[i], lw=0.6)

    # Filling the yearly shadow zone only when abs(lat) less than 66 deg
    if np.abs(lat) < 66:
        ax.fill_between(x_box_centers[0, :], y_box_centers[0, :], np.ones(len(x_box_centers[0,:]))*max_val, color="grey", alpha=0.6, label="Yearly surface\nshadow zone")
        # Limiting the shadow zone northwards
        ax.fill_between(x_box_centers[-1, :], y_box_centers[-1, :], np.ones(len(x_box_centers[-1,:]))*max_val, color="white", alpha=1)
        
    # Plotting hourly shadows
    for i in np.arange(len(n_days)):
        for j in np.arange(24):
            if i==0 and j==0:
                rect = plt.Rectangle((x_box_lls[i, j], y_box_lls[i, j]), x_shadows[i, j], y_shadows[i, j], angle=(180-saas[i, j]), facecolor="black", alpha=0.8, label="Instantaneous hourly\nsurface shadows")
                ax.add_patch(rect)
            else:
                rect = plt.Rectangle((x_box_lls[i, j], y_box_lls[i, j]), x_shadows[i, j], y_shadows[i, j], angle=(180-saas[i, j]), facecolor="black", alpha=0.8)
                ax.add_patch(rect)
	    
    # Plotting LTB position at the ground
    ax.scatter(np.zeros(1), np.zeros(1), color='red', marker=".", label='LTB', lw=0.2)
    
    # Reference surface at the ground
    rect = plt.Rectangle((-1000, -20000*lat_sign), 2000, 500, facecolor="black")
    ax.add_patch(rect)
    ax.text(-3000, -23000*lat_sign, '2000 m', fontsize=6)
    ax.text(-9000, -20000*lat_sign, '500 m', fontsize=6)

    ax.set_title('lat = '+str(lat)+' N' )

    # Orthonormal axes
    plt.xlim(-50000,50000)
    if lat_sign < 0:
    	plt.ylim(-60000,40000)
    else:
    	plt.ylim(-40000,60000)
    plt.xlabel("Eastward distance with respect to LTB [m]")
    plt.ylabel("Northward distance with respect to LTB [m]")
    plt.legend(fontsize=8)
    plt.tight_layout()
    
    # Save figure
    plt.savefig('Fig3_shadows_LTB_simple_model.png')


def plot_daily_prod(lat, lon, n_days, max_power):
    """Function to plot the idealized daily power production.

    Parameters
    ----------   
    lat : float
        Latitude where the LTB is located, in [degrees North].

    lon : float
        Longitude where the LTB is located, in [degrees East].

    n_days : list (Python object)
        List of days for which we wish to plot the surface shadows at every hour of the day, expected string format for each day [yyyymmdd].

    max_power : float
        Maximum power produced by the structure, in [W].

    """

    # Defining hourly 2-day matrices
    daytime_flag = np.zeros((2,24))
    hour_power = np.zeros((2,24))

    # Only the 2 solstices are ploted, first and last elements of n_days
    for i in [0, -1]:
        # Loop over the 24 hours of the day
        for j in np.arange(24):
            # Defining the datetime object at local hour j at the
            # geographical location corresponding to the given longitude
            date = dt.datetime(int(n_days[i][0:4]), int(n_days[i][4:6]), int(n_days[i][6:8]), j, tzinfo=dt.timezone(offset=dt.timedelta(hours=int(lon/15))))
            # Getting the Solar Zenithal Angle
            sza = pysol.get_altitude(lat, lon, date)
            # If the sun is above the horizon, daytime_flag = 1
            if sza > 0:
                daytime_flag[i, j] = 1
    # Computing idealized daily power for the 2 days, in MW
    hour_power = daytime_flag*max_power/1e6

    fig = plt.figure(figsize=(8, 4))
    plt.suptitle('Figure 4: Maximum and minimum daily energy production, lat = '+str(lat)+' N', fontsize=14)

    # Left subplot 
    plt.subplot(1, 2, 1)
    plt.plot(np.arange(24), hour_power[0, :], color='red', label='Total daily energy\nproduction = '+str(np.sum(daytime_flag[0, :])*max_power/1e9)+' GWh')
    plt.title('a) Summer solstice power production')
    plt.ylabel('Power [MW]')
    plt.xlabel('Solar hour [no unit]')
    plt.legend(loc='center left', fontsize=9)

    # Right subplot
    plt.subplot(1, 2, 2)
    plt.plot(np.arange(24), hour_power[-1, :], color='blue', label='Total daily energy\nproduction = '+str(np.sum(daytime_flag[-1, :])*max_power/1e9)+' GWh')
    plt.title('b) Winter solstice power production')
    plt.ylabel('Power [MW]')
    plt.xlabel('Solar hour [no unit]')
    plt.legend(loc='center left', fontsize=9)
    plt.tight_layout()

    # Save figure
    plt.savefig('Fig4_daily_prod_LTB_simple_model.png')
