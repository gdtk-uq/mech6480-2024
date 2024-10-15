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

# Moments
density_eq = sp.Eq(sum(fs), rho)
x_velocity_eq = sp.Eq(sum([cix[i]*fs[i] for i in range(9)]), rho*u_w)
y_velocity_eq = sp.Eq(sum([ciy[i]*fs[i] for i in range(9)]), 0) # u_y = 0

x_solve = sp.solve(x_velocity_eq, sum([fs[i] for i in unknowns]))
density_eq = density_eq.subs(sum([fs[i] for i in unknowns]), x_solve[0])
rho_eq = sp.Eq(rho,sp.solve(density_eq, rho)[0])
print(sp.ccode(rho_eq))

feq_equations = []
for i in range(len(fs)):
    feq_equations.append(sp.Eq(feqs[i], wi[i]*rho*(1 + 3*(cix[i]*u_w ) \
                                                    + 9/2*(cix[i]*u_w )**2 \
                                                    - 3/2*(u_w**2 ))))

#NEBB
eq1 = sp.Eq(f3 + f3_eq - f1_eq, f1)
eq1 = eq1.subs(f3_eq, feq_equations[3].rhs)
eq1 = eq1.subs(f1_eq, feq_equations[1].rhs).simplify()
print(sp.ccode(eq1))
eq2 = sp.Eq(f7 + f7_eq - f5_eq + Ny, f5)
eq2 = eq2.subs(f7_eq, feq_equations[7].rhs)
eq2 = eq2.subs(f5_eq, feq_equations[5].rhs).simplify()
print(sp.ccode(eq2))
eq3 = sp.Eq(f6 + f6_eq - f8_eq - Ny, f8)
eq3 = eq3.subs(f6_eq, feq_equations[6].rhs)
eq3 = eq3.subs(f8_eq, feq_equations[8].rhs).simplify()
print(sp.ccode(eq3))

# Solve for tangential momentum
y_velocity_eq = y_velocity_eq.subs(f1, eq1.rhs)
y_velocity_eq = y_velocity_eq.subs(f5, eq2.rhs)
y_velocity_eq = y_velocity_eq.subs(f8, eq3.rhs)
y_velocity_eq = sp.solve(y_velocity_eq, Ny)[0].simplify()
y_velocity_eq = sp.Eq(Ny, y_velocity_eq)

eq1 = eq1.subs(Ny, y_velocity_eq.rhs)
eq2 = eq2.subs(Ny, y_velocity_eq.rhs)
eq3 = eq3.subs(Ny, y_velocity_eq.rhs)
print("final expressions")
print(sp.ccode(eq1))
print(sp.ccode(eq2))
print(sp.ccode(eq3))