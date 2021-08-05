#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions needed by the LTB_simple_model.py script to compute
   the main features of an indealized Light Tower of Babel.

This file contains all the thermodynamic formulas needed in order to
compute all the atmospheric and LTB-related variables (pressure, mass,
...) describing the LTB with such hypotheses/configuration.

May 2021, v0 by Rodrigo Guzman

"""

import numpy as np

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
