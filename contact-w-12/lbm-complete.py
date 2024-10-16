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
    NX = 300; NY = 80
    N_ITER = 50_000
    plot_every = N_ITER//100
    cyl_cx = NX//4; cyl_cy = NY//2+2; 
    cyl_r = NY//10

    RE = 100
    RHO0 = 1.0
    U0 = 0.05
    nu = U0 * 2 * cyl_r / RE
    tau = 3*nu + 0.5
    omega = 1.0/tau
    collision = 'MRT'

    # TODO: Create the lattice and solid flags
    xx, yy = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    solid = np.zeros((NX, NY), dtype=bool)
    solid[(xx-cyl_cx)**2 + (yy-cyl_cy)**2 <= cyl_r**2] = True

    # TODO: Initialize the lattice with RHO0, U0
    rho = np.ones((NX, NY), dtype=float) * RHO0
    u_x = np.ones((NX, NY), dtype=float) * U0
    u_x[solid] = 0.0
    u_y = np.zeros((NX, NY), dtype=float)

    def f_equilibrium(rho, u_x, u_y):
        #TODO: Implement the equilibrium distribution function
        feq = np.zeros((DISCRETE_VELOCITIES, NX, NY), dtype=float)
        feq[0,:,:] = 0.444444444444444*rho*(-1.5*u_x**2 - 1.5*u_y**2 + 1)
        feq[1,:,:] = 0.111111111111111*rho*(3.0*u_x**2 + 3*u_x - 1.5*u_y**2 + 1)
        feq[2,:,:] = 0.111111111111111*rho*(-1.5*u_x**2 + 3.0*u_y**2 + 3*u_y + 1)
        feq[3,:,:] = 0.111111111111111*rho*(3.0*u_x**2 - 3*u_x - 1.5*u_y**2 + 1)
        feq[4,:,:] = 0.111111111111111*rho*(-1.5*u_x**2 + 3.0*u_y**2 - 3*u_y + 1)
        feq[5,:,:] = 0.0277777777777778*rho*(-1.5*u_x**2 + 3*u_x - 1.5*u_y**2 + 3*u_y + 4.5*(u_x + u_y)**2 + 1)
        feq[6,:,:] = 0.0277777777777778*rho*(-1.5*u_x**2 - 3*u_x - 1.5*u_y**2 + 3*u_y + 4.5*(-u_x + u_y)**2 + 1)
        feq[7,:,:] = 0.0277777777777778*rho*(-1.5*u_x**2 - 3*u_x - 1.5*u_y**2 - 3*u_y + 4.5*(-u_x - u_y)**2 + 1)
        feq[8,:,:] = 0.0277777777777778*rho*(-1.5*u_x**2 + 3*u_x - 1.5*u_y**2 - 3*u_y + 4.5*(u_x - u_y)**2 + 1)
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
        rho = np.sum(f_prec, axis=0)
        u_x = np.einsum('ijk, i -> jk', f_prec, DIRECTIONS[0,:]) / rho
        u_y = np.einsum('ijk, i -> jk', f_prec, DIRECTIONS[1,:]) / rho

        # 1a) OUTPUT HERE
        if iter % plot_every == 0:
            toc = time()
            if iter>0: print(f"Iteration {iter}: max velocity = {np.max(np.sqrt(u_x**2 + u_y**2)):.4f}.....speed = {NX*NY*plot_every/(toc-tic)/1e6:.4f} MLUPs")
            ax.clear()
            ax.imshow(np.sqrt(u_x**2 + u_y**2).T, cmap='viridis', origin='lower')
            ax.plot( NX//2 + 100*u_x[NX//2,:], yy[NX//2,:], 'r')
            circle = plt.Circle((cyl_cx, cyl_cy), cyl_r, color='red')
            ax.add_patch(circle)
            plt.pause(0.001)
            tic = time()
            pass

        # 2) TODO: Calculate equilibrium
        f_eq = f_equilibrium(rho, u_x, u_y)

        # 3) TODO: Collision step
        if collision == 'BGK':
            f_post = f_prec - omega*(f_prec - f_eq)
        elif collision == 'MRT':
            f_post = collisionMRT(f_prec, rho, u_x, u_y)

        # 4) Streaming step
        # for k in range(DISCRETE_VELOCITIES):
        #     for i in range(NX):
        #          for j in range(NY):
        #                 i2 = (i+DIRECTIONS[0,k]) % NX
        #                 j2 = (j+DIRECTIONS[1,k]) % NY
        #                 f_prec[k,i2,j2] = f_post[k,i,j]
        for i in range(DISCRETE_VELOCITIES):
             f_prec[i,:,:] = np.roll(f_post[i,:,:], DIRECTIONS[:,i], axis=(0,1))

        # 5) Boundary conditions (bounce-back)
        f_tmp = f_prec.copy()
        for i in range(DISCRETE_VELOCITIES):
            f_prec[i,solid] = f_tmp[BOUNCED[i],solid]

        # 5) Boundary conditions (wet nodes)
        # 5a) Inlet - Zou/He
        rho_w = (-f_prec[0,0,:] - f_prec[2,0,:] - 2*f_prec[3,0,:] - f_prec[4,0,:] - 2*f_prec[6,0,:] - 2*f_prec[7,0,:])/(U0 - 1)
        f_prec[1,0,:] = 1.0*f_prec[3,0,:] + 0.66666666666666663*rho_w*U0
        f_prec[5,0,:] = -0.5*f_prec[2,0,:] + 0.5*f_prec[4,0,:] + 1.0*f_prec[7,0,:] + 0.16666666666666666*rho_w*U0
        f_prec[8,0,:] = 0.5*f_prec[2,0,:] - 0.5*f_prec[4,0,:] + 1.0*f_prec[6,0,:] + 0.16666666666666666*rho_w*U0

        #f1 == 1.0*f3 + 0.66666666666666663*rho*u_w
        #f5 == -0.5*f2 + 0.5*f4 + 1.0*f7 - 0.16666666666666666*rho*u_w
        #f8 == 0.5*f2 - 0.5*f4 + 1.0*f6 - 0.16666666666666666*rho*u_w

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