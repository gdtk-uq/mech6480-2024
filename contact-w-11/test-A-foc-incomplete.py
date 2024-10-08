# test-A-foc-incomplete.py

from math import sin, pi
from numpy import zeros

def test_f(x):
    # -- 1. --
    # fill in test function for sine wave

def update_in_time(f, a, L, ncells, dt, nsteps):
    dx = float(L)/float(ncells)
    xs = zeros(ncells)
    us = zeros(ncells)
    fs = zeros(ncells+1)

    # Initialise positions and starting values
    for i in range(ncells):
        xs[i] = # -- 2. -- Fill in cell centre locations
        us[i] = # -- 3. -- Fill in starting values at cell centres

    # Now begin timestepping
    for n in range(nsteps):
        # 1. Compute fluxes.
        # 1a. At all INTERNAL interfaces
        for i in range(1,ncells):
            fs[i] = # -- 4. -- Fill in flux function based on first-order centred scheme
        # 1b. At the BOUNDARY interfaces
        fs[0] = 0.0; fs[-1] = a*us[-1]

        # 2. Compute new u value
        for i in range(ncells):
            us[i] = # -- 5. -- Fill in time update
        
    return xs, us

dt = 0.01
nsteps = 50
L = 5.0
a = 1.0
ncells = 100
dx = L/ncells

# generate analytical solution
# -- 6. -- Generate analytical solution
# Save in xa, ua

# -- 7. -- Generate numerical solution
# Call update_in_time
# Save in xn, un



import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(xa, ua, "-", xn, un, "o")
plt.show()


