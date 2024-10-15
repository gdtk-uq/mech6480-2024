import sympy as sp
import numpy as np

# Create symbols to generate lattice boltzmann kernels
DISCRETE_VELOCITIES = 9
DIRECTIONS = sp.Matrix([[0, 1, 0, -1, 0, 1, -1, -1, 1],\
                        [0, 0, 1, 0, -1, 1, 1, -1, -1]])
WEIGHTS = sp.Matrix([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])

# calc equilibrium expressions
