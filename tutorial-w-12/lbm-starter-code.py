# Libraries
import numpy as np
import matplotlib.pyplot as plt
try:
    import tqdm
    progress_bar = True
except ImportError:
    progress_bar = False
    pass

# Lattice constants
DISCRETE_VELOCITIES = 9
DIRECTIONS = np.array([(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1),\
                       (1, 1), (-1, 1), (-1, -1), (1, -1)]).T
BOUNCED = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

WEIGHTS = np.array([ 4/9, \
                     1/9,  1/9,  1/9,  1/9, \
                    1/36, 1/36, 1/36, 1/36], dtype=np.float64)

def main(Re = 100, stop = 50_000):
    cyl_r  = 5
    NX = cyl_r*2 * 25
    NY = cyl_r*2 * 20
    N_ITER = stop
    plot_every = N_ITER // 100
    cyl_cx = NX // 4; cyl_cy = NY // 2
    
    RE   = Re
    RHO0 = 1.0
    U0   = 0.05
    nu   = U0 * 2 * cyl_r / RE
    tau  =  3 * nu + 0.5
    omega = 1 / tau
    print(f"omega = {omega}")
    collision = 'MRT'

    # Create the lattice and solid flags
    xx, yy = np.meshgrid(np.arange(NX), np.arange(NY), indexing='ij')
    solid = np.zeros((NX, NY), dtype=bool)
    solid[(xx - cyl_cx)**2 + (yy - cyl_cy)**2 <= cyl_r**2] = True
    
    # Initialize the lattice with RHO0, U0
    rho = RHO0 * np.ones((NX, NY), dtype=float)
    u_x = U0*np.ones((NX, NY), dtype=float)
    u_x[solid] = 0
    u_y = np.zeros((NX, NY), dtype=float)

    def f_equilibrium(rho, u_x, u_y):
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
    fig, (ax) = plt.subplots(1, 1, figsize=(10, 10*NY/NX))
    
    if progress_bar: iterations = tqdm.tqdm(range(N_ITER))
    else: iterations = range(N_ITER)

    for iter in iterations:
        # 1) Macroscopic quantities
        rho = np.sum(f_prec, axis=0)
        u_x = np.einsum("ijk, i -> jk",f_prec,DIRECTIONS[0,:]) / rho
        u_y = np.einsum("ijk, i -> jk",f_prec,DIRECTIONS[1,:]) / rho

        # 1a) OUTPUT HERE
        if iter % plot_every == 0:
            if progress_bar != True: print(f"Iteration {iter}")
            ax.clear()
            ax.imshow(np.sqrt(u_x**2 + u_y**2).T, cmap='viridis', origin='lower')
            ax.plot( NX//2 + 100*u_x[NX//2,:], yy[NX//2,:], 'r')
            circle = plt.Circle((cyl_cx, cyl_cy), cyl_r, color='red')
            ax.add_patch(circle)
            plt.pause(0.001)
            pass

        # 2) Calculate equilibrium
        f_eq = f_equilibrium(rho, u_x, u_y)

        # 3) Collision step
        if collision == 'BGK':
            f_post = f_prec - omega*(f_prec - f_eq)
        elif collision == 'MRT':
            f_post = collisionMRT(f_prec, rho, u_x, u_y) # much more stable than BGK

        # 4) Streaming step
        for i in range(DISCRETE_VELOCITIES):
            if iter == 0: print(f"Streaming {i}, direction {DIRECTIONS[:,i]}, shape {f_post[i,:,:].shape}")
            f_prec[i,:,:] = np.roll(f_post[i,:,:], DIRECTIONS[:,i], axis=(0,1)) # input 1: array, input 2: shift, axis: tuple of ints

        # 5a) Boundary conditions (bounce-back)
        f_tmp = f_prec.copy()
        for i in range(DISCRETE_VELOCITIES):
            f_prec[i,solid] = f_tmp[BOUNCED[i],solid]

        # 5b) slip boundary on top and bottom using bounce-back like approach
        #     for these we only 'bounce' the normal component of ci
        f_prec[4,:,-1]   = f_post[2,:,-1]
        f_prec[8,1:,-1]  = f_post[5,:-1,-1]
        f_prec[8,0,-1]   = f_post[5,-1,-1]
        f_prec[7,:-1,-1] = f_post[6,1:,-1]
        f_prec[7,-1,-1]  = f_post[6,0,-1]
 
        f_prec[2,:,0]   = f_post[4,:,-1]
        f_prec[6,1:,0]  = f_post[7,:-1,-1] 
        f_prec[6,0,0]   = f_post[7,-1,-1]
        f_prec[5,:-1,0] = f_post[8,1:,-1]
        f_prec[5,-1,0]  = f_post[8,0,-1]

        # 6) Boundary conditions (wet nodes)
        # 6a) Inlet - Zou/He
        rho_w = 1/(1-U0) * (f_prec[0,0,:] + f_prec[2,0,:] + f_prec[4,0,:] + 2*(f_prec[3,0,:] + f_prec[6,0,:] + f_prec[7,0,:]))
        f_prec[1,0,:] =  1.0*f_prec[3,0,:] + 0.66666666666666663*rho_w*U0
        f_prec[5,0,:] = -0.5*f_prec[2,0,:] + 0.5*f_prec[4,0,:] + 1.0*f_prec[7,0,:] + 0.16666666666666666*rho_w*U0
        f_prec[8,0,:] =  0.5*f_prec[2,0,:] - 0.5*f_prec[4,0,:] + 1.0*f_prec[6,0,:] + 0.16666666666666666*rho_w*U0
                                
        # # 6b) Outlet - zero-gradient for unknown pops
        f_prec[3,-1,:] = f_prec[3,-2,:]
        f_prec[7,-1,:] = f_prec[7,-2,:]
        f_prec[6,-1,:] = f_prec[6,-2,:]  
    return 0

Re = 60
main(Re)