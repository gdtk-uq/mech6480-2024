# Libraries
import numpy as np
import matplotlib.pyplot as plt
from time import time

# Lattice constants
DISCRETE_VELOCITIES = 9
DIRECTIONS = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),\
                       (1, 1), (-1, 1), (-1, -1), (1, -1)]).T
BOUNCED = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

WEIGHTS = np.array([ 4/9, \
                     1/9,  1/9,  1/9,  1/9, \
                    1/36, 1/36, 1/36, 1/36], dtype=np.float64)

def main():
    # TODO: Set up simulation params

    # TODO: Create the lattice and solid flags

    # TODO: Initialize the lattice with RHO0, U0

    def f_equilibrium(rho, u_x, u_y):
        #TODO: Implement the equilibrium distribution function
        feq = np.zeros((DISCRETE_VELOCITIES, NX, NY), dtype=float)
        return feq
    
    def collisionMRT(f, rho, Jx, Jy, S2 = 1-omega, S3=0, S4=0):
            rho = f[8,:,:] + f[7,:,:] + f[6,:,:] + f[5,:,:] + f[4,:,:] + f[3,:,:] + f[2,:,:] + f[1,:,:] + f[0,:,:];                         
            Jx = f[8,:,:] - f[7,:,:] - f[6,:,:] + f[5,:,:] - f[3,:,:] + f[1,:,:];  
            Jy = -f[8,:,:] - f[7,:,:] + f[6,:,:] + f[5,:,:] - f[4,:,:] + f[2,:,:];
            R3 = -f[4,:,:] - f[2,:,:] - f[0,:,:] + ( f[8,:,:] + f[7,:,:] + f[6,:,:] + f[5,:,:] + f[3,:,:] + f[1,:,:] )*2.;    
            R4 = -f[3,:,:] - f[1,:,:] - f[0,:,:] + ( f[8,:,:] + f[7,:,:] + f[6,:,:] + f[5,:,:] + f[4,:,:] + f[2,:,:] )*2.;
            R5 = ( -f[8,:,:] + f[7,:,:] - f[6,:,:] + f[5,:,:] )*3.;
            R6 = f[4,:,:] - f[2,:,:] + ( -f[8,:,:] - f[7,:,:] + f[6,:,:] + f[5,:,:] )*2.;
            R7 = f[3,:,:] - f[1,:,:] + ( f[8,:,:] - f[7,:,:] - f[6,:,:] + f[5,:,:] )*2.;
            R8 = f[0,:,:] + ( -f[4,:,:] - f[3,:,:] - f[2,:,:] - f[1,:,:] + ( f[8,:,:] + f[7,:,:] + f[6,:,:] + f[5,:,:] )*2. )*2.;   

            R3 = R3 - Jx*Jx/rho*3.;
            R4 = R4 - Jy*Jy/rho*3.;
            R5 = R5 - Jy*Jx/rho*3.;
            R6 = R6;
            R7 = R7;
            R8 = R8;            
            R3 = S2*R3;
            R4 = S2*R4;        
            R5 = S2*R5;  
            R6 = S3*R6;  
            R7 = S3*R7;  
            R8 = S4*R8;
            #Jx = Jx + GravitationX*rho;
            #Jy = Jy + GravitationY*rho;
            R3 = R3 + Jx*Jx/rho*3.;
            R4 = R4 + Jy*Jy/rho*3.;
            R5 = R5 + Jy*Jx/rho*3.;
            R6 = R6;     
            R7 = R7;    
            R8 = R8; 
            f[0,:,:] = ( R8 + ( -R4 - R3 + rho*2 )*2 )/9.;
            f[1,:,:] = ( -R8 - R4 - R7*3 + ( R3 + rho + Jx*3 )*2 )/18.;
            f[2,:,:] = ( -R8 - R3 - R6*3 + ( R4 + rho + Jy*3 )*2 )/18.;
            f[3,:,:] = ( -R8 - R4 + R7*3 + ( R3 + rho - Jx*3 )*2 )/18.;
            f[4,:,:] = ( -R8 - R3 + R6*3 + ( R4 + rho - Jy*3 )*2 )/18.;
            f[5,:,:] = ( R8 + R4 + R3 + rho + ( R7 + R6 + R5 + Jy + Jx )*3 )/36.;
            f[6,:,:] = ( R8 + R4 + R3 + rho + ( -R7 + R6 - R5 + Jy - Jx )*3 )/36.;
            f[7,:,:] = ( R8 + R4 + R3 + rho + ( -R7 - R6 + R5 - Jy - Jx )*3 )/36.;
            f[8,:,:] = ( R8 + R4 + R3 + rho + ( R7 - R6 - R5 - Jy + Jx )*3 )/36.; 
            return f
    
    f_prec = f_equilibrium(rho, u_x, u_y)
    f_post = np.zeros_like(f_prec)
    # Main loop
    fig, ax = plt.subplots(figsize=(10, 10*NY/NX))
    tic = time()
    for iter in range(N_ITER):
        # 1) TODO: Macroscopic quantities

        # 1a) OUTPUT HERE
        if iter % plot_every == 0:
            toc = time()
            print(f"Iteration {iter}: max velocity = {np.max(np.sqrt(u_x**2 + u_y**2)):.4f}.....speed = {NX*NY*plot_every/(toc-tic)/1e6:.4f} MLUPs")
            ax.clear()
            ax.imshow(np.sqrt(u_x**2 + u_y**2).T, cmap='viridis', origin='lower')
            ax.plot( NX//2 + 100*u_x[NX//2,:], yy[NX//2,:], 'r')
            circle = plt.Circle((cyl_cx, cyl_cy), cyl_r, color='red')
            ax.add_patch(circle)
            plt.pause(0.001)
            tic = time()
            pass

        # 2) TODO: Calculate equilibrium

        # 3) TODO: Collision step

        # 4) Streaming step
        
        # 5) Boundary conditions (bounce-back)
        
        # 5) Boundary conditions (wet nodes)
        # 5a) Inlet - Zou/He
        
                                
        # # 5b) Outlet - zero-gradient for unknown pops
        f_prec[3,-1,:] = f_prec[3,-2,:]
        f_prec[7,-1,:] = f_prec[7,-2,:]
        f_prec[6,-1,:] = f_prec[6,-2,:]

        # 5c) Top/Bottom walls (could also do these as bounceback, but boundary would then be between wall node and fluid node)
        rho_w = -(f_prec[0,:,0] + f_prec[1,:,0] + f_prec[3,:,0] + 2*f_prec[4,:,0] + 2*f_prec[7,:,0] + 2*f_prec[8,:,0])/(0.0 - 1)
        f_prec[2,:,0] = 1.0*f_prec[4,:,0] + 0.66666666666666663*rho_w*0.0
        f_prec[5,:,0] = -0.5*f_prec[1,:,0] + 0.5*f_prec[3,:,0] + 1.0*f_prec[7,:,0] + 0.5*rho_w*0.0 + 0.16666666666666666*rho_w*0.0
        f_prec[6,:,0] = 0.5*f_prec[1,:,0] - 0.5*f_prec[3,:,0] + 1.0*f_prec[8,:,0] - 0.5*rho_w*0.0 + 0.16666666666666666*rho_w*0.0
        
        rho_w = (f_prec[0,:,-1] + f_prec[1,:,-1] + 2*f_prec[2,:,-1] + f_prec[3,:,-1] + 2*f_prec[5,:,-1] + 2*f_prec[6,:,-1])/(0.0 + 1)
        f_prec[4,:,-1] = 1.0*f_prec[2,:,-1] - 0.66666666666666663*rho_w*0.0
        f_prec[7,:,-1] = 0.5*f_prec[1,:,-1] - 0.5*f_prec[3,:,-1] + 1.0*f_prec[5,:,-1] - 0.5*rho_w*0.0 - 0.16666666666666666*rho_w*0.0
        f_prec[8,:,-1] = -0.5*f_prec[1,:,-1] + 0.5*f_prec[3,:,-1] + 1.0*f_prec[6,:,-1] + 0.5*rho_w*0.0 - 0.16666666666666666*rho_w*0.0
        
    return 0

main()