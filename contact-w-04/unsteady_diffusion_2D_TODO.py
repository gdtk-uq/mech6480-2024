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
xx, yy = np.meshgrid(x, y, indexing='ij')
print(f"dx={dx:.4f}, dy={dy:.4f}")

# Time
simulated_time = 0
iteration = 0
T_END = 62.5e-3 
DT = 0.0000625
PLOT_EVERY = 10

# Initial and boundary conditions
T_COOL = 300.0
T_HOT  = 700.0
T_WALL = 300.0

T = 
T[:,0] = 
T[:,-1] = 
T[0,:] = 
T[-1,:] = 

for i in range(NX):
    for j in range(NY):
        if ((xx[i,j]-CX)**2 + (yy[i,j]-CY)**2) < r2:
            T[i,j] = T_HOT

# plot I.C.
plt.contourf(xx, yy, T, cmap='hot')
plt.colorbar()
plt.show()

# allocate memory for fluxes
x_flux = 
y_flux = 

# Plot setup
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(1, 1)
ax0 = fig.add_subplot(gs[0, 0])  

tic = time.time()
while simulated_time < T_END:
    # Calculate fluxes - interior
    x_flux[1:-1,:] = 
    y_flux[:,1:-1] = 

    # Calculate fluxes - Boundaries
    x_flux[0,:]  = 
    x_flux[-1,:] = 
    y_flux[:,0]  = 
    y_flux[:,-1] = 

    # Update T
    T =  T + 

    # Plotting
    if np.mod(iteration, PLOT_EVERY) == 0:
        ax0.cla()
        ax0.contourf(xx, yy, T, vmin=T_COOL, vmax=T_HOT, cmap='hot')
        ax0.set_xlabel('x');  ax0.set_ylabel('y')
        ax0.set_title('Temperature ' + str(round(simulated_time,5)) + ' s')
        ax0.set_aspect('equal')
        plt.pause(0.1)

    # Update time    
    simulated_time += DT
    iteration += 1

# plot final T distribution
ax0.cla()
ax0.contourf(xx, yy, T, vmin=T_COOL, vmax=T_HOT, cmap='hot')
ax0.set_xlabel('x');  ax0.set_ylabel('y')
ax0.set_title('Temperature ' + str(round(simulated_time,5)) + ' s')
ax0.set_aspect('equal')
# plot profile through the center
plt.figure()
plt.plot(x, T[:, NY//2])
plt.xlabel('x')
plt.ylabel('T')
plt.title('Temperature profile y=L/2')

print("Total elapsed time:", round(time.time()-tic,2), "s (", round((time.time()-tic)/60.0,2), "min)")
plt.show()

