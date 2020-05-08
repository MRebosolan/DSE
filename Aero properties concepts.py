from math import cos, radians

def CLdes(q, W_start_cruise, W_end_cruise, S, sweep):
    CL = 1.1/1 * (0.5 * (W_start_cruise/S + W_end_cruise/S))
    Cl = CL/(cos(radians(sweep)))^2
    return (CL, cl)
