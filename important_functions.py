#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:06:41 2020

@author: mcollier
"""

#############################################################
import networkx as nx
import random as rnd
import numpy as np
import os
import statistics as st
###############################################################

def get_networks_with_directory(directory):
    """
    This function returns a list of networks from a directory with multiple .graphml files
    """
    
    #first get files
    filelist = [filename for filename in sorted(os.listdir(os.path.abspath(directory))) if filename.endswith(".graphml")]
    networks = []
    os.chdir(directory)
    for file in filelist:
        network = nx.read_graphml(file)
        if len(network.nodes()) > 0:
            
            networks.append(network)
    return(networks)

def get_avg_excess_degree(G):
    """
    This function gets the avg excess degree of the network using only avg degree
    (homogenous network)
    """
    print(len(G.nodes()))
    d = [G.degree(node) for node in list(G.nodes())]
    avg_deg = round(np.mean(d),3)

    xdeg = avg_deg - 1

    return xdeg

def get_xdeg_with_deghet(G):
    """
    This function get the excess degree of a network incorporating the 
    deg het of that network    
    """
    deg_list = [G.degree(node) for node in list(G.nodes())]
    avg_deg = round(np.mean(deg_list),3)
    var_deg = round(st.variance(deg_list),3)
    numer = var_deg + (avg_deg ** 2) - avg_deg
    xdeg = numer/avg_deg
    
    return xdeg

def infected_neighbors(G, node, infected_list):
    """Calculates the number of infected neighbors for every susceptible node"""
    infected_deg = [x for x in list(G.neighbors(node)) if x in infected_list]   
    length = len(infected_deg)
    return length, infected_deg

########################
def simulation(Network,T):
# This function implements a percolation simulation to find R0
       
    net_size = Network.number_of_nodes()
    # Initialize variables for the list of infected and recovered individuals
    infected = []
    recovered = []

    ##################
    # Choose one node to infect (patient zero) so that outbreak can be seeded
    p_zero = rnd.choice(list(Network.nodes())) # Randomly choose one node from the network
    infected = [p_zero]                  # The node p_zero is now infected
    infected_count = 1
    R0_count = []  
  
    gen1_infected =[]
    while infected:
        
        infector = infected[0]      
        for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 

            if neigh not in infected and neigh not in recovered: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                if rnd.random() < T:                     # if infector does infect neigh
                    gen1_infected.append(neigh) #append them to the gen1 infected list
                    infected_count = infected_count +1
        
        infected.remove(infector)
        recovered.append(infector)

######### For generation 2 ###############    
    infected = gen1_infected
    gen2_infected = []
    
    if len(infected) > 0:
        
        while infected:       
            num_infected = 0
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
                
                if neigh not in infected and neigh not in recovered and neigh not in gen2_infected: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                    if rnd.random() < T:                     # if infector does infect neigh
                        gen2_infected.append(neigh)
                        infected_count = infected_count +1
                        num_infected = num_infected +1
       
            R0_count.append(num_infected) # add the number they infected to the R0 list
            infected.remove(infector)
            recovered.append(infector)
        
        infected = gen2_infected
    else: infected = []
    
##### For gen 3 ################## 
    gen3_infected = []
    if len(infected) > 0:

        while infected:
            num_infected = 0
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
               
                if neigh not in infected and neigh not in recovered and neigh not in gen3_infected: # check if this neighbor is susceptible

                # figure out if infector is successful at infecting neighbor "neigh"
                    if rnd.random() < T:                     # if infector does infect neigh
                        gen3_infected.append(neigh)
                        infected_count = infected_count +1
                        num_infected = num_infected +1
            
            R0_count.append(num_infected) # add the number they infected to the R0 list
            infected.remove(infector)
            recovered.append(infector)

    
            
        infected = gen3_infected
    else: infected = []

    #FOR THE REMAINIG SIMULATION
    if len(infected) > 0:
        while infected:        
            infector = infected[0]
            for neigh in list(Network.neighbors(infector)): # for all the nodes connected to (i.e. neighbors of) the infector 
                
                if neigh not in infected and neigh not in recovered: # check if this neighbor is susceptible

                    if rnd.random() < T:                     
                        infected.append(neigh)
                        infected_count = infected_count +1        

            infected.remove(infector)
            recovered.append(infector)

    if len(R0_count) > 0:
        R0 = round(np.mean(R0_count),1)
    else: R0 = 0
    epi_size = infected_count/net_size
    return R0, T, epi_size
    


def many_simulations(num_sims,network, T):
    """
    This function finds an average R0 value for a certain T value on a network
    It will only calculate R0 if the epidemic probability is 10% and the epi size is 10%
    """   
    R0_list = []
    
    for x in range(num_sims):
        #print(x)
        R0, T, size = simulation(network,T)
        
        if size >= 0.1: #only want epidemics so 10% or more of the network infected
                R0_list.append(R0)   
    
    epi_prob = len(R0_list)/num_sims
    if epi_prob >= 0.1: #basically epidemic probability is 10% or higher
        R0 = round(np.mean(R0_list),1)   
    else: R0 = 0

        

    return R0, epi_prob



