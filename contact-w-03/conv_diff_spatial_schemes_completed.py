#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MECH6480 - Week 3
Study steady state convection-diffusion with different spatial schemes.
Implementation of to do prints completed in lecture/contact session.

- central difference
- upwind
- hybrid
"""

__author__ = "Travis Mitchell"
__version__ = "1.0.0"
__date__ = "2024-08-07"

import numpy as np
import matplotlib.pyplot as plt

def conv_diff_dirichlet(n, dx, F, D, boundary_left, boundary_right,method="centered"):
    """ Generate solution matrix for the FVM """
    if method == "centered":
        a_W = D + 0.5*F
        a_E = D - 0.5*F    
        S_p_left  = (2.0*D + F)
        S_p_right = (2.0*D - F)
    elif method == "upwind":
        a_W = D + F
        a_E = D
        S_p_left = (2.0*D + F)
        S_p_right = 2.0*D
    elif method == "hybrid":
        a_W = np.max([F,D+0.5*F,0])
        a_E = np.max([-F,D-0.5*F,0])
        S_p_left = (2.0*D + F)
        S_p_right = 2.0*D
    else:
        raise NotImplementedError("Method not implemented. Pick one of the following: 'centered', 'upwind', 'hybrid'")
    a_P = a_W + a_E

    matrix_A = np.diag(np.repeat(-a_W,n-1),-1) + \
               np.diag(np.repeat(a_P,n),0) + \
               np.diag(np.repeat(-a_E,n-1),1)
    matrix_A[0,0]   = a_E + S_p_left
    matrix_A[-1,-1] = a_W + S_p_right

    # Make source vector
    vector_b = np.zeros(n)
    vector_b[0]  = S_p_left * boundary_left
    vector_b[-1] = S_p_right * boundary_right
    return (matrix_A, vector_b.T)

def generate_matrix_and_solve(num_cells, dx, F, D, phi_0, phi_L, method = "centered"):
    solution_matrix, source_matrix = conv_diff_dirichlet(num_cells, dx, F, D, phi_0, phi_L, method)
    concentration = np.linalg.solve(solution_matrix, source_matrix)
    return concentration

def L2_norm(num_cells, x_location, concentration, function):
    L2 = 0
    for cell in range(num_cells):
        L2 += (concentration[cell] - function(x_location[cell]))**2
    return np.sqrt(L2/num_cells)

#def analytic_solution(x, phi_0, phi_L, rho, u, Gamma):
#    return phi_0 + (phi_L - phi_0) * (np.exp(rho*u*x/Gamma)-1) / (np.exp(rho*u*length/Gamma) -1)

if __name__ == "__main__":
    # SYSTEM PARAMETERS:
    Gamma = 0.1             # kg/m.s      
    u = 2.5                 # m/s
    rho = 1.0               # kg/m3
    length = 1              # m
    phi_0 = 1
    phi_L = 0
    def analytic_solution(x):
        return phi_0 + (phi_L - phi_0) * (np.exp(rho*u*x/Gamma)-1) / (np.exp(rho*u*length/Gamma) -1)
   
    # GRID GENERATION
    num_cells = 5              #[-]
    dx = length / (num_cells)  #[m]
    x_locations = np.linspace(0.5*dx,(num_cells-0.5)*dx,num_cells)
    
    # CONVECTIVE AND DIFFUSIVE TERMS
    F = rho * u
    D = Gamma / dx
    
    # SYSTEM OF EQUATIONS
    concentration = generate_matrix_and_solve(num_cells, dx, F, D, phi_0, phi_L, method = "hybrid")
    
    test_resolutions =[5, 8, 16, 32, 64, 128, 256, 512]
    errors = np.zeros((len(test_resolutions),3)) # Update: create a row for each method
    dx_store = np.zeros(len(test_resolutions))
    for indx,res in zip(range(len(test_resolutions)),test_resolutions):
        dx = length / (res)  #[m]
        dx_store[indx] = dx    
        
        # Update: Sorry forgot to update diffusive conductance in the lecture!
        F = rho * u
        D = Gamma / dx
        
        
        x_locations = np.linspace(0.5*dx,(res-0.5)*dx,res)
        concentration = generate_matrix_and_solve(res, dx, F, D, phi_0, phi_L, method = "upwind")
        errors[indx][0] = L2_norm(res, x_locations, concentration, analytic_solution)
        
        # Update: calc errors for convergence of other discretisations
        concentration = generate_matrix_and_solve(res, dx, F, D, phi_0, phi_L, method = "centered")
        errors[indx][1] = L2_norm(res, x_locations, concentration, analytic_solution)
        concentration = generate_matrix_and_solve(res, dx, F, D, phi_0, phi_L, method = "hybrid")
        errors[indx][2] = L2_norm(res, x_locations, concentration, analytic_solution)
    
    x = np.linspace(0,length,100)
    solution = analytic_solution(x)
    plt.plot(x, solution, 'r--')
            
    plt.plot(x_locations, concentration, 'b-o')
    plt.xlabel('Distance along hallway (m)')
    plt.ylabel('concentration')
    
    fig,ax = plt.subplots(1)
    ax.plot(dx_store, errors[:,0],'ro')
    
    # Update: add in extra plots and indicating gradient lines
    ax.plot(dx_store, errors[:,1], 'bo')
    ax.plot(dx_store, errors[:,2], 'go')
    y_tmp = dx_store**2
    ax.plot(dx_store, y_tmp, 'k-')
    ax.set_xlabel("dx [m]")
    ax.set_ylabel("L2 norm")
    
    ax.set_xscale("log")
    ax.set_yscale("log")
    
