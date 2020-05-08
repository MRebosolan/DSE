from math import cos, radians


def CLdes(q, W_start_cruise, W_end_cruise, S, sweep): #finds design cruise CL for aircraft and airfoil
    CL = 1.1/1 * (0.5 * (W_start_cruise/S + W_end_cruise/S))
    Cl = CL/(cos(radians(sweep)))^2
    return (CL, Cl)


def taper (chord_r, chord_t):
    return chord_t/chord_r

def MAC (taper, chord_r):
    mac = (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
    return mac
