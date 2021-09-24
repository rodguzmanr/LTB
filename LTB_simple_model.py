#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Script to compute the main features of an indealized Light Tower of Babel

This script needs a few variables to be defined in the main part of the
program in order to compute all the other variables (height, mass, ...)
describing the LTB with such hypotheses/configuration.
Some plots to illustrate the output are also created at the end of it.

May 2021, v0 by Rodrigo Guzman

"""

import numpy as np
import numpy.ma as ma

import constants as cst
import variables as var
import functions as func


if __name__ == "__main__":
    """Script to compute main features of the idealized LTB.

    """

    ######################################################################
    # I - Atmosphere variables derived from variables.py and function.py #
    ######################################################################

    # Number of modules needed to build the tower, also corresponds to the
    # number of elements for the vertically extended (profile) variables 
    nb_mod = int(var.z_tower/var.z_mod)

    # Initializaing the atmosphere variables profiles (nb_mod elements)
    # Temperature profile
    temp = np.zeros(nb_mod)
    # Pressure profile
    pres = np.zeros(nb_mod)
    # Altitude profile
    alt = np.arange(nb_mod)*var.z_mod + var.surf_alt
    # Computing idealized profiles from functions
    # Temperature
    temp = func.compute_temp(temp, var.surf_temp, var.temp_grad, var.min_temp, alt, var.surf_alt)
    # Pressure
    pres = func.compute_pres(pres, var.surf_pres, alt, var.surf_alt, var.temp_grad, var.surf_temp)

    ######################################################################
    ##  II - Tower variables derived from variables.py and function.py  ##
    ######################################################################
    
    # > < : ces signes ne marchent pas sur mon clavier...
   
    # Number of casings per module (vertical profile), initialized with
    # 1 casing per module 
    nb_cube_mod = np.ones(nb_mod)
    # Initializing number of casings per module side
    nb_cube_mod_per_side = 1
    # Diameter of each module (vertical profile), initialized with the
    # smallest diameter allowed by x_cube and x_cube_margin
    d_mod = np.ones(nb_mod)*var.min_d_tower
    # Mass of each module (vertical profile), in [kg]
    mass_mod = np.zeros(nb_mod)
    # Mass of each module (vertical profile), in [kg]
    max_vol_mod = np.zeros(nb_mod)
    # Mass of the air volume displaced for each module (vertical profile), in [kg]
    mass_air_mod = np.zeros(nb_mod)
    # Mass of the h2 volume needed for each module (vertical profile), in [kg]
    mass_h2_mod = np.zeros(nb_mod)

    # Computing each module features from surface to top of the tower
    for i in np.arange(nb_mod):
        # Computing mass of the module
        mass_mod[i] = var.dens_surf_mod*var.z_mod*d_mod[i]
        # Computing maximum gas volume allowed within the module
        max_vol_mod[i] = var.max_vol_cube*(1-var.max_vol_margin/100.)*nb_cube_mod[i]
        # Air mass associated to max_vol_mod 
        mass_air_mod[i] = func.mass_from_ideal_gas(pres[i], max_vol_mod[i], cst.mol_mass_air, temp[i])
        # H2 mass associated to max_vol_mod 
        mass_h2_mod[i] = func.mass_from_ideal_gas(pres[i], max_vol_mod[i], cst.mol_mass_h2, temp[i])
        # While H2 volume within the module is not enough to lift the module
        while (mass_air_mod[i]-mass_h2_mod[i]) < mass_mod[i]:
            # Add 1 casing per module side
            nb_cube_mod_per_side = nb_cube_mod_per_side+1
            # Resulting in an increase of ((nb+1)^2 - nb^2) number of casings
            nb_cube_mod[i] = nb_cube_mod_per_side**2
            # Updating the minimum module diameter associated to new nb_cube_mod
            d_mod[i] = var.x_cube*np.sqrt(2)*(1+var.x_cube_margin/100.)*nb_cube_mod_per_side
            # Computing maximum gas volume allowed within the module
            max_vol_mod[i] = var.max_vol_cube*(1-var.max_vol_margin/100.)*nb_cube_mod[i]
            # Air mass associated to max_vol_mod 
            mass_air_mod[i] = func.mass_from_ideal_gas(pres[i], max_vol_mod[i], cst.mol_mass_air, temp[i])
            # H2 mass associated to max_vol_mod 
            mass_h2_mod[i] = func.mass_from_ideal_gas(pres[i], max_vol_mod[i], cst.mol_mass_h2, temp[i])
            # Computing mass of the module
            mass_mod[i] = var.dens_surf_mod*var.z_mod*d_mod[i]
        # Updating the number of casings needed for the modules above
        nb_cube_mod[i:] = nb_cube_mod[i]
        # Updating the minimum module diameter for the modules above
        d_mod[i:] = d_mod[i]

    # Computing total number of casings within the tower
    nb_cube_tower = int(np.sum(nb_cube_mod))
    # Total mass of the tower
    max_mass_tower = np.sum(mass_mod)
    # Total h2 mass needed to lift the tower
    mass_h2_tower = np.sum(mass_h2_mod)
    # If the tower diameter has to increase constantly with altitude,
    # tower becomes an inverted cone
    d_mod_max = d_mod[-1]
    d_mod_min = d_mod[0]
    d_grad = (d_mod_max-d_mod_min)/nb_mod
    d_mod_cone = np.arange(nb_mod)*d_grad + d_mod_min
    

    ######################################################################
    # III - PV panel variables derived from variables.py and function.py #
    ######################################################################

    # Maximum gas volume within the cubic casings attached to the PV panel,
    # with the max_vol_margin safety margin
    max_vol_pv = var.nb_cube_pv*var.max_vol_cube*(1-var.max_vol_margin/100.)
    # Air mass associated to max_vol_pv 
    mass_air_pv = func.mass_from_ideal_gas(pres[-1], max_vol_pv, cst.mol_mass_air, temp[-1])
    # H2 mass associated to max_vol_pv 
    mass_h2_pv = func.mass_from_ideal_gas(pres[-1], max_vol_pv, cst.mol_mass_h2, temp[-1])
    # Maximum total mass of the PV panel
    max_mass_pv = mass_air_pv-mass_h2_pv
    # Maximum surface density of the PV panel
    max_dens_pv = max_mass_pv/var.surf_pv

    
    ######################################################################
    #####  IV - Total structure features derived from results above  #####
    ######################################################################

    # Computing total structure mass
    total_mass_structure = max_mass_tower + max_mass_pv
    # Computing total H2 mass needed to lift all the structure
    total_mass_h2 = mass_h2_tower + mass_h2_pv
    # Computing ratio between total H2 mass needed to lift the structure
    # and the total mass of the structure
    total_mass_ratio = total_mass_h2/total_mass_structure

    
    #######################################################################
    ## V - Electricity production calculations derived from variables.py ##
    #######################################################################

    # Computing maximum power produced by the structure
    max_power = var.surf_pv*cst.S_flux*(var.conv_rate_pv/100.)
    # Computing average daily energy production
    day_energy_prod = max_power*var.hours_day
    # Computing energy production for a full year
    year_energy_prod = day_energy_prod*365

    
    #######################################################################
    #########  VI - Displaying main results from the simulation  ##########
    #######################################################################

    print('\n*** Main results from this simulation ***\n')
    print('Main input variables:')
    print('x_cube = '+str(var.x_cube)+' [m]')
    print('dens_surf_mod = '+str(var.dens_surf_mod)+' [kg.m^-2]')
    print('x_pv = '+str(var.x_pv)+' [m]')
    print('y_pv = '+str(var.y_pv)+' [m]')
    print('\n***********\n')    
    print('Main structure output variables:')
    print('max_mass_pv = '+str(np.round(max_mass_pv/1000.))+' [ton]')
    print('max_dens_pv = '+str(np.round(max_dens_pv, decimals=3))+' [kg.m^-2]')
    print('max_mass_tower + max_mass_pv = '+str(np.round(max_mass_tower/1000.))+' + '+str(np.round(max_mass_pv/1000.))+' [ton] = ')
    print('total_mass_structure = '+str(np.round(total_mass_structure/1000.))+' [ton]')
    print('mass_h2_tower + mass_h2_pv = '+str(np.round(mass_h2_tower/1000.))+' + '+str(np.round(mass_h2_pv/1000.))+' [ton] =')
    print('total_mass_h2 = '+str(np.round(total_mass_h2/1000.))+' [ton]')
    print('total_mass_h2/total_mass_structure =')
    print('total_mass_ratio = '+str(np.round(total_mass_ratio, decimals=5))+' [no unit]')
    print(str(np.round(total_mass_ratio*100., decimals=2))+' % of the mass of the structure is H2')
    print('nb_cube_tower + nb_cube_pv = '+str(nb_cube_tower)+' + '+str(var.nb_cube_pv)+' =')
    print('total_nb_cube_structure = '+str(var.nb_cube_pv+nb_cube_tower))
    print('total_nb_mod_tower = '+str(nb_mod))
    print('\n***********\n')    
    print('Main production output variables:')
    print('max_power = '+str(max_power/1e6)+' [MW]')
    print('day_energy_prod = '+str(day_energy_prod/1e9)+' [GWh/day]')
    print('year_energy_prod = '+str(year_energy_prod/1e12)+' [TWh/year]\n')   


    ########################################
    ############  VII - Plots  #############
    ########################################
    
    # Ploting Atmospheric and LTB profiles
    func.plot_profiles(pres, alt, nb_cube_mod, temp, d_mod, d_mod_cone)

    # Ploting surface shadows caused by the LTB
    func.plot_shadows(var.lat, var.lon, var.n_days, var.x_pv, var.y_pv, var.z_tower)
    
    # Ploting idealized daily production for the 2 solstices
    func.plot_daily_prod(var.lat, var.lon, var.n_days, max_power)
