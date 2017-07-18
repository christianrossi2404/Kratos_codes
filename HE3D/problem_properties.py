
# copper properties
# SPECIFIC_HEAT  0.39 //[KJ/(kg K)]
# DENSITY   8940 // [kg/m³]
# CONDUCTIVITY  401 //[W/(m K)]

# HeatEquation3D
# HeatEquationNeumannCondition3D

# air properties
# DENSITY 1.205 // [kg/m³]
# SPECIFIC_HEAT 1.005 //[KJ/(kg K)]
# CONDUCTIVITY 0.0257 //[W/(m K)]
# VISCOSITY   1.47760E-05
import numpy as np
const = np.arcsin(1)
print(const)

DENSITY = 1.205
VISCOSITY =  1.47760E-05
u1 = 5
u2 = u1* 0.2/0.0625
print(u2)
Re1 = 0.2 *  u1 / VISCOSITY
print(Re1)
Re2 = 0.0625*u2   /VISCOSITY
print(Re2)
