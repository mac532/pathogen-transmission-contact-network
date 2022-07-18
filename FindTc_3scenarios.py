#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:19:19 2020

@author: melissacollier
"""


#############################################################
import numpy as np
import os
import csv
import important_functions as inf
###############################################################

os.chdir("mac532/pathogen-transmission-contact-network")

Fluid_nets = inf.get_networks_with_directory("mac532/pathogen-transmission-contact-network/232_Network_Subsample/FluidExchange")
print("Got Fluid Networks")
Nonphys_nets = inf.get_networks_with_directory("mac532/pathogen-transmission-contact-network/232_Network_Subsample/Nonphysical")
print("Got Nonphysical Networks")
Phys_nets = inf.get_networks_with_directory("mac532/pathogen-transmission-contact-network/232_Network_Subsample/Physical")
print("Got Physical Networks")
Ind_nets = inf.get_networks_with_directory("mac532/pathogen-transmission-contact-network/232_Network_Subsample/Indirect")
print("Got Indirect Networks")

T_range = np.arange(0.01, 0.8, 0.01)
R0 = 1


os.chdir("mac532/pathogen-transmission-contact-network")
with open('Tc_3Scenarios.csv', 'w', encoding='utf-8-sig') as mycsvfile:
    writer = csv.writer(mycsvfile, lineterminator = '\n')
 
#writer = csv.writer(open('Simulation_Results_GLMM.csv','w', encoding='utf-8-sig'))
    writer.writerow([ "NetID","transmission.route", "T_ad", "T_dh", "T_sim", "Prob"
                    ])
    i=1
    for net in Fluid_nets:
        NetName = i
        G = net
        xdeg_ad = inf.get_avg_excess_degree(G)   
        xdeg_dh = inf.get_xdeg_with_deghet(G)
        
        if xdeg_ad <= 1 :
            T_ad = "NA"
        else:
            T_ad = R0/xdeg_ad

        if xdeg_dh <= 1:
            #dhIP = "NA"
            T_dh = "NA"
        else:
            T_dh = R0/xdeg_dh
   
        for t in T_range:
            print("f",i, t)
            R0, prob = inf.many_simulations(250,G, t)
            if R0 > 1:
                T_sim = t
                break
            else:
                T_sim = "NA"

        elements = [ NetName, "Fluid", T_ad, T_dh, T_sim, prob
                        ]
        writer.writerow(elements)
            
        i=i+1  

    i=1
    for net in Phys_nets:
        NetName = i
        G = net
        xdeg_ad = inf.get_avg_excess_degree(G)   
        xdeg_dh = inf.get_xdeg_with_deghet(G)
        
        if xdeg_ad <= 1 :
            T_ad = "NA"
        else:
            T_ad = R0/xdeg_ad

        if xdeg_dh <= 1:
            T_dh = "NA"
        else:
            T_dh = R0/xdeg_dh
     
        for t in T_range:
            print("p",i, t)
            R0, prob = inf.many_simulations(250,G, t)
            if R0 > 1:
                T_sim = t
                break
            else:
                T_sim = "NA"

        elements = [ NetName, "Phys", T_ad, T_dh, T_sim, prob
                        ]
        writer.writerow(elements)
            
        i=i+1     

    i=1
    for net in Nonphys_nets:
        NetName = i
        G = net
        xdeg_ad = inf.get_avg_excess_degree(G)   
        xdeg_dh = inf.get_xdeg_with_deghet(G)
        
        if xdeg_ad <= 1 :
            T_ad = "NA"
        else:
            T_ad = R0/xdeg_ad

        if xdeg_dh <= 1:
            T_dh = "NA"
        else:
            T_dh = R0/xdeg_dh
    
        for t in T_range:
            print("np",i, t)
            R0, prob = inf.many_simulations(250,G, t)
            if R0 > 1:
                T_sim = t
                break
            else:
                T_sim = "NA"

        elements = [ NetName, "Nonphys", T_ad, T_dh, T_sim, prob
                        ]
        writer.writerow(elements)
            
        i=i+1   
        
    i=1
    for net in Ind_nets:
        NetName = i
        G = net
        xdeg_ad = inf.get_avg_excess_degree(G)   
        xdeg_dh = inf.get_xdeg_with_deghet(G)
        
        if xdeg_ad <= 1 :
            T_ad = "NA"
        else:
            T_ad = R0/xdeg_ad

        if xdeg_dh <= 1:
            T_dh = "NA"
        else:
            T_dh = R0/xdeg_dh
   
        for t in T_range:
            print("i",i, t)
            R0, prob = inf.many_simulations(250,G, t)
            if R0 > 1:
                T_sim = t
                break
            else:
                T_sim = "NA"

        elements = [ NetName, "Ind", T_ad, T_dh, T_sim, prob
                        ]
        writer.writerow(elements)
            
        i=i+1   


             