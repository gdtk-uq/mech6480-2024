#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MECH6480 - Week 3
Study steady state diffusion and convection-diffusion

- demo 'wiggles'

@author: uqtmitc3
"""

__author__ = "Travis Mitchell"
__version__ = "1.0.0"
__date__ = "2024-08-07"

import numpy as np
import matplotlib.pyplot as plt
import time as time

def steady_diffusion(k, dx, num_cells, T_A, T_B, q=0):
    """
        solve steady diffusion term with dirichlet bounds.
    """
    # internal nodes
    a_E = k / dx
    a_W = k / dx
    a_P = a_E + a_W
    b_P = -q*dx
    # boundary conditions
    S_u = -2.0 * k / dx
    a_P_star = 3.0 * k / dx
    
    # Set up system of equations
    A = (np.diag(np.repeat(a_W,num_cells-1),-1) + 
         np.diag(np.repeat(-a_P,num_cells),0) + 
         np.diag(np.repeat(a_E,num_cells-1),1) )
    # Overwrite boundary rows
    A[0,0] = -a_P_star
    A[-1,-1] = -a_P_star
    # set up RHS
    b = np.ones(num_cells) * b_P
    b[0] += S_u * T_A
    b[-1] += S_u * T_B
    return A, b.transpose()

def steady_convection(F, num_cells, phi_A, phi_B):
    """
        Solve convection term with Dirichlet bounds.
    """
        # Interior
    a_E = F/2
    a_W = -F/2
    a_P = a_E + a_W
    b_P = 0

    # Boundaries
    S_u = F
    a_P_star = F/2

    # Set up the linear system of equations
    A = (np.diag(np.repeat(a_W,num_cells-1),-1) + 
         np.diag(np.repeat(a_P,num_cells),0) +     
         np.diag(np.repeat(a_E,num_cells-1),1) )  
    
    # Overwrite the first and last rows for boundaries
    A[0,0] = a_P_star
    A[-1,-1] = -a_P_star

    # Set up the RHS
    b = np.ones(num_cells) * b_P
    b[0]  += S_u * phi_A
    b[-1] += -S_u * phi_B
    return A, b.transpose()

if __name__ == "__main__":
    # Problem Settings
    gamma = 0.1     # kg/m.s
    length = 1      # m
    rho = 1.2       # kg/m3
    u = 0.1         # m/s
    phi_0 = 1       
    phi_L = 0      
    q = 0           
    
    # Discretisation
    num_cells = 5
    dx = length / num_cells
    x = np.linspace(0.5*dx,(num_cells-0.5)*dx,num_cells)
    
    D, b_D = steady_diffusion(gamma, dx, num_cells, phi_0, phi_L)
    A, b_A = steady_convection(rho*u, num_cells, phi_0, phi_L)
    T = np.linalg.solve(D-A, b_D-b_A)
    plt.plot(x, T, 'bo')
    
    def analytic_soln(x):
        return phi_0 + (phi_L - phi_0) * (np.exp(rho*u*x/gamma)-1) \
            / (np.exp(rho*u*length/gamma)-1)
   
    x = np.linspace(0, length, 100)
    soln = analytic_soln(x)
    plt.plot(x,soln,'r--')