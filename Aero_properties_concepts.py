from math import cos, radians


#fuckthepolice

def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad*h
    P = 101300 * (T/288.15)**(-9.81/(T_grad*287))
    rho = P/(T*287)
    a = (1.4*287*T)**0.5
    return (T,P,rho,a)


def CLdes(q,W_start_cruise,W_end_cruise,S,sweep): #finds design cruise CL for aircraft and airfoil
    CL = 1.1/q * (0.5 * (W_start_cruise/S + W_end_cruise/S))
    cl = CL/(cos(radians(sweep)))**2
    return CL, cl

def MAC (chord_r,taper):

    return (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
   
    
def root_chord_tip_chord(S,wingspan, taper):
    chord_r = (2*S)/(wingspan*(1+taper))
    chord_t = chord_r * taper
    return chord_r,chord_t

def dynamic_pressure(rho, a, M):
    V = M * a
    return 0.5 * rho * V**2


sweep = 25
wingspan = 24
W_start_cruise = 31537
W_end_cruise = 28337
taper = 0.4
S = 70.6 #crj
q = 1
