# # MECH6480 - Finite volume methods
# ## Week 03 
# 
# ### Tutorial helper
#
# Author: @TravisMitchell - t.mitchell@uq.edu.au
# 

import matplotlib.pyplot as plt
import numpy as np

def analytic_solution(r, R, mu, dp):
    "this should match the analytical solution you derived last week"
    return (R**2/(4.*mu)) * (-dp) * (1-(r/R)**2)

def myL2_norm(velocity, x_locations, function, N, R, mu, dp):
    "L2 norm from lecture but with tailored function input"
    L2_norm  = 0
    for cell in range(N):
        L2_norm  += (velocity[cell] - function(x_locations[cell], R, mu, dp))**2
    return np.sqrt(L2_norm/N)

def generate_matrix_dirichlet_neumann(n, dr, r_locations, mu, dp):
    """ Generate solution matrix for the FVM """
    raise NotImplementedError("Implement this based on your discretisation")
    r_n = 
    r_s = 
    a_N = 
    a_S = 
    a_P = a_N + a_S

    # if your coefficients are arrays then we can generate with the below
    matrix_A = np.diag(-a_S[1:],-1) + \
               np.diag(a_P,0) + \
               np.diag(-a_N[:-1],1)
    # apply wall boundary condition
    matrix_A[-1,-1] = 

    # determine b vector
    vector_b = 

    return matrix_A, vector_b.transpose()

if __name__ == "__main__":
    # SYSTEM PARAMETERS:
    mu = 0.001              
    dp = -0.001
    R  = 1

    # BOUNDARIES:
    phi_0 = 1
    phi_L = 0

    # GRID GENERATION
    num_cells = 16        #[-]
    dr = R / (num_cells)  #[m]
    r_locations = np.linspace(0.5*dr,(num_cells-0.5)*dr,num_cells)

    # SYSTEM OF EQUATIONS
    solution_matrix, source_matrix = generate_matrix_dirichlet_neumann(num_cells, dr, r_locations, mu, dp)
    velocity = np.linalg.solve(solution_matrix, source_matrix)

    # ANALYTICAL SOLUTION
    r = np.linspace(0,R,100)
    solution = analytic_solution(r, R, mu, dp)
    plt.plot(solution, r, 'r--')
            
    # PLOTTING
    plt.plot(velocity, r_locations, 'bo')
    plt.xlabel('velocity (m/s)')
    plt.ylabel('radial location')
    plt.legend(["exact soln","fvm soln"])
    plt.show()

    # TODO: Calculate the L2 norm for various resolutions 
    #       and compare to expected rate of convergence
    resolutions = [8, 16, 32, 64, 128, 256, 512]
