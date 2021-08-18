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
import numpy.ma as ma
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
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


def plot_profiles(pres, alt, nb_cube_mod, temp, d_mod, d_mod_cone):
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
        Diameter of the LTB modules (from the surface to the tropopause), no unit.

    d_mod_cone : array (Numpy object)
        Modules' diameter for a (upside down) conical LTB (from the surface to the tropopause), no unit.

    """

    fig = plt.figure(figsize=(8, 6))

    # Top left subplot 
    plt.subplot(2, 2, 1)
    plt.plot(pres, alt, color='blue')
    plt.title('Atmospheric pressure')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Pressure [Pa]')

    # Top right subplot 
    plt.subplot(2, 2, 2)
    plt.plot(nb_cube_mod, alt, color='black')
    plt.title('Number of casings per module')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Nb casings [no unit]')
    
    # Bottom left subplot
    plt.subplot(2, 2, 3)
    plt.plot(temp, alt, color='red')
    plt.title('Atmospheric temperature')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Temperature [K]')

    # Bottom right subplot 
    plt.subplot(2, 2, 4)
    plt.plot(d_mod, alt, color='black', label='cylindrical modules tower')
    plt.plot(d_mod_cone, alt, color='green', label='cone tower')
    plt.title('Tower diameter')
    plt.ylabel('Altitude [m]')
    plt.xlabel('Diameter [m]')
    plt.legend(fontsize=9)

    plt.tight_layout()
    # Save figure
    plt.savefig('profiles_LTB_simple_model.png')


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
    
    # Defining coordinate matrices for the days to be plotted
#    x_box_centers = np.ones((len(n_days),24))*-9999
#    y_box_centers = np.ones((len(n_days),24))*-9999
    
    # Creating figure
    fig, ax = plt.subplots()

    # Loop over the number of days in n_days
    for i in np.arange(len(n_days)):
        # Loop over the 24 hours of the day
        for j in np.arange(24):
            # Defining the datetime object at local hour j at the
            # location corresponding to the given longitude
            date = dt.datetime(int(n_days[i][0:4]), int(n_days[i][4:6]), int(n_days[i][6:8]), j, tzinfo=dt.timezone(offset=dt.timedelta(hours=int(lon/15))))
            # Getting the Solar Zenithal Angle
            sza = pysol.get_altitude(lat, lon, date)
            # Computation keeps going if the sun is above the horizon
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

                # We save the box center coordinates to plot them
#                x_box_centers[i, j] = x_box_center
#                y_box_centers[i, j] = y_box_center
                
                rect = plt.Rectangle((x_box_ll, y_box_ll), x_shadow, y_shadow, angle=(180-saa), facecolor="black", alpha=0.8)
                ax.add_patch(rect)

    # Masking useless elements of the coordinate matrix
#    x_box_centers = ma.masked_where(x_box_centers <= -9999., x_box_centers)
#    y_box_centers = ma.masked_where(y_box_centers <= -9999., y_box_centers)

#    for i in np.arange(len(n_days)):
        # Plotting each day's surface shadow trajectory
#        ax.plot(x_box_centers[i, :], y_box_centers[i, :], label=n_days[i], lw=1.5)

    # LTB position at the ground, try scatterplot
    ax.plot(np.zeros(2), np.zeros(2), '-o', color='blue', label='LTB', lw=3)
    
    # Reference surface at the ground
    rect = plt.Rectangle((-1000, -20000), 2000, 500, facecolor="black")
    ax.add_patch(rect)
    ax.text(-3000, -24000, '2000 m', fontsize=6)
    ax.text(-8000, -20000, '500 m', fontsize=6)

    ax.set_title('Hourly surface shadows')
#    ax.plot([-500, 500], [-20000, -20000], color="black", lw=3)

    plt.xlim(-70000,70000)
    plt.ylim(-30000,100000)
    plt.legend()
    plt.tight_layout()
    # Save figure
    plt.savefig('shadows_LTB_simple_model.png')
