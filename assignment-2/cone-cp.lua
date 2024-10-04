-- Compute Cp on cone
-- Invoke at command line:
-- $ e4shared --custom-script --script-file=cone-cp.lua

U_inf = 2720; p_inf = 100.0; T_inf = 300.0; theta=math.rad(25.0)
R_N2 = 296.8
rho_inf = p_inf/(R_N2*T_inf)
q_inf = 0.5*rho_inf*U_inf^2

-- First, compute beta (stream deflection angle)
-- since it is an input to the cone properties
-- calculator
beta = idealgasflow.beta_cone(U_inf, p_inf, T_inf, theta)

-- Second, compute properties at cone surface
theta_c, V_c, p_c, T_c = idealgasflow.theta_cone(U_inf, p_inf, T_inf, beta)

-- Third, compute C_p using definition
C_p = (p_c - p_inf)/q_inf
print("------------------------------------------\n")
print("  C_p at cone surface= ", C_p)
print("")
print("------------------------------------------")


