import sympy as sp

# Define symbols
rho, u_w, u_y, Ny = sp.symbols('rho u_w u_y, Ny')
f0, f1, f2, f3, f4, f5, f6, f7, f8 = sp.symbols('f0 f1 f2 f3 f4 f5 f6 f7 f8')
f0_eq, f1_eq, f2_eq, f3_eq, f4_eq, f5_eq, f6_eq, f7_eq, f8_eq = sp.symbols('f0_eq f1_eq f2_eq f3_eq f4_eq f5_eq f6_eq f7_eq f8_eq')
cix = sp.Matrix([0, 1, 0, -1,  0, 1, -1, -1,  1])
ciy = sp.Matrix([0, 0, 1,  0, -1, 1,  1, -1, -1])
wi = sp.Matrix([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])
fs = sp.Matrix([f0, f1, f2, f3, f4, f5, f6, f7, f8])
feqs = sp.Matrix([f0_eq, f1_eq, f2_eq, f3_eq, f4_eq, f5_eq, f6_eq, f7_eq, f8_eq])
unknowns = [1, 5, 8]