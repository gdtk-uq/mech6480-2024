import matplotlib.pyplot as plt
import numpy as np
import time

# Parameters
LX = 10 # m
LY = 10 # m
ALPHA = 4
R = 2; r2 = R**2
CX = 5
CY = 5

# Grid
NX = 100
NY = 100
dx = LX / NX
dy = LY / NY
x = np.linspace(dx/2., LX-dx/2.0, NX) # C.V. centres
y = np.linspace(dy/2., LY-dy/2.0, NY)
xx, yy = np.meshgrid(x, y, indexing='ij') # NOTE: flip indices
                                          #       to match discretisation
print(f"dx={dx:.4f}, dy={dy:.4f}")

# Time
simulated_time = 0
iteration = 0
T_END = 62.5e-3 
DT = 0.0000625
PLOT_EVERY = 100

# Initial and boundary conditions
T_COOL = 300.0
T_HOT  = 700.0
T_WALL = 300.0

T = np.ones((NX,NY))*T_COOL

for i in range(NX):
    for j in range(NY):
        if ((xx[i,j]-CX)**2 + (yy[i,j]-CY)**2) < r2:
            T[i,j] = T_HOT

# plot I.C.
#plt.contourf(xx, yy, T, cmap='hot')
#plt.colorbar()
#plt.show()

# allocate memory for fluxes
x_flux = np.zeros((NX+1,NY))
y_flux = np.zeros((NX,NY+1))

# Plot setup
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(1, 1)
ax0 = fig.add_subplot(gs[0, 0])  

tic = time.time()
while simulated_time < T_END:
    # Calculate fluxes - interior
    x_flux[1:-1,:] = ALPHA*(T[1:,:]-T[:-1,:])/dx
    y_flux[:,1:-1] = ALPHA*(T[:,1:]-T[:,:-1])/dy

    # # Calculate fluxes - Boundaries
    x_flux[0,:]  = ALPHA*(T[0,:]-T_WALL)/(dx/2.)
    x_flux[-1,:] = ALPHA*(T_WALL-T[-1,:])/(dx/2.)
    y_flux[:,0]  = ALPHA*(T[:,0]-T_WALL)/(dy/2.)
    y_flux[:,-1] = ALPHA*(T_WALL-T[:,-1])/(dy/2.)

    # # Update T
    T =  T + DT*(dy*(x_flux[1:,:]-x_flux[:-1,:])\
                +dx*(y_flux[:,1:]-y_flux[:,:-1]))/(dx*dy)
    # Update time    
    simulated_time += DT
    iteration += 1
    
    # Plotting
    if np.mod(iteration, PLOT_EVERY) == 0:
        ax0.cla()
        ax0.contourf(xx, yy, T, vmin=T_COOL, vmax=T_HOT, cmap='hot')
        ax0.set_xlabel('x');  ax0.set_ylabel('y')
        ax0.set_title('Temperature ' + str(round(simulated_time,5)) + ' s')
        ax0.set_aspect('equal')
        #plt.pause(0.1)
        fig.savefig(f'heatmap_{iteration:04d}.png')

# plot final T distribution
ax0.cla()
ax0.contourf(xx, yy, T, vmin=T_COOL, vmax=T_HOT, cmap='hot')
ax0.set_xlabel('x');  ax0.set_ylabel('y')
ax0.set_title('Temperature ' + str(round(simulated_time,5)) + ' s')
ax0.set_aspect('equal')
plt.show()
# plot profile through the center
plt.figure()
plt.plot(x, T[:, NY//2])
plt.xlabel('x')
plt.ylabel('T')
plt.title('Temperature profile y=L/2')

print("Total elapsed time:", round(time.time()-tic,2), "s (", round((time.time()-tic)/60.0,2), "min)")
plt.show()

