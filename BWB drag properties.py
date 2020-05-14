import Aero_properties_concepts as ap
from math import degrees, radians

#planform parameters estimates
AR_BWB = 5.75
l_BWB = 15
b_BWB = 27
S_BWB = 128.1
W_start_cruise_BWB = 31537.95 * 9.81
W_end_cruise_BWB = 28337.66 * 9.81
taper_BWB = 0.3
altitude = 11000
M = 0.8
quarter_sweep_BWB = radians(50)

q = ap.dynamic_pressure(M, altitude)

CL_BWB = ap.CLdes(q, W_start_cruise_BWB, W_end_cruise_BWB, S_BWB, quarter_sweep_BWB)[0]
CL2_BWB = W_start_cruise_BWB / (q * W_start_cruise_BWB)
LEsweep_BWB =(ap.LE_sweep(quarter_sweep_BWB, AR_BWB, taper_BWB))

e_BWB = ap.oswald_factor(quarter_sweep_BWB, AR_BWB)
CD_induced_BWB = ap.induced_drag(CL_BWB, AR_BWB, e_BWB)
print(CL_BWB, CD_induced_BWB, e_BWB)






