""" 
    Solve steady state problems using the finite volume method.
"""

__author__ = "Travis Mitchell"
__version__ = "1.0.0"
__date__ = "2024-07-31"

import numpy as np
import matplotlib.pyplot as plt
import time as time

def steady_diffusion(k, dx, num_cells, T_A, T_B, q=0):
    """
        Solve d(k.dT/dx)/dx + q = 0 for T using the finite volume method. With Dirichlet T_A and T_B boundary conditions.   
    """
    # Interior
    a_E = k / dx
    a_W = k / dx
    a_P = a_E + a_W
    b_P = -q * dx

    # Boundaries
    S_u = -2 * k / dx
    a_P_star = 3 * k / dx

    # Set up the linear system of equations
    A = (np.diag(np.repeat(a_W,num_cells-1),-1) + 
         np.diag(np.repeat(-a_P,num_cells),0) +     
         np.diag(np.repeat(a_E,num_cells-1),1) )  
    
    # Overwrite the first and last rows for boundaries
    A[0,0] = -a_P_star
    A[-1,-1] = -a_P_star

    # Set up the RHS
    b = np.ones(num_cells) * b_P
    b[0]  += S_u * T_A
    b[-1] += S_u * T_B

    return A, b.transpose()

def steady_convection(F, num_cells, phi_A, phi_B):
    """
        Solve d(rho u phi)/dx = 0 for T using the finite volume method. With Dirichlet T_A and T_B boundary conditions.
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
    problem = 3
    # SYSTEM PARAMETERS:
    if problem == 1:
        k = 1000                # W/m.K       
        length = 0.5            # m
        q = 0

        # BOUNDARIES:
        T_A = 100
        T_B = 500

        def analytic_solution(x):
            return (T_B - T_A)/length * x + T_A
    elif problem == 2:
        k = 0.5                # W/m.K
        q = 1000e3             # W/m^3
        length = 0.02          # m

        # BOUNDARIES:
        T_A = 100
        T_B = 200

        def analytic_solution(x):
            return ( (T_B - T_A)/length + q*(length - x)/(2*k) )*x + T_A
    elif problem == 3:
        gamma = 0.1
        q = 0
        rho = 1
        length = 1
        u = 2.5

        phi_0 = 1
        phi_L = 0
        def analytic_solution(x):
            return phi_0 + (phi_L - phi_0) * (np.exp(rho*u*x/gamma)-1) / (np.exp(rho*u*length/gamma) -1)
    else:
        print('Invalid problem number')

    # GRID GENERATION
    num_cells = 20              #[-]
    dx = length / (num_cells)  #[m]
    x_locations = np.linspace(0.5*dx,(num_cells-0.5)*dx,num_cells)

    if problem < 3:
        # Diffusion only problem examples
        A, b = steady_diffusion(k, dx, num_cells, T_A, T_B, q)
        tic = time.time()
        T = np.linalg.solve(A,b)
        toc = time.time()
        print('Time to solve linalg: ', toc-tic)

        # PLOT
        x_analytic = np.linspace(0,length,100)
        plt.plot(x_analytic, analytic_solution(x_analytic), 'r--')
        plt.plot(x_locations, T, 'bo')
        plt.xlabel('x [m]')
        plt.ylabel('T [K]')
        plt.show()
    else:
        # convection-diffusion problem example
        D, b_D = steady_diffusion(gamma, dx, num_cells, phi_0, phi_L, q)
        A, b_A = steady_convection(rho*u, num_cells, phi_0, phi_L)
        tic = time.time()
        phi = np.linalg.solve(D-A,b_D-b_A)
        toc = time.time()
        print('Time to solve linalg: ', toc-tic)

        # PLOT
        x_analytic = np.linspace(0,length,100)
        plt.plot(x_analytic, analytic_solution(x_analytic), 'r--')
        plt.plot(x_locations, phi, 'bo')
        plt.xlabel('x [m]')
        plt.ylabel('phi')
        plt.show()
