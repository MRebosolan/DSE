from math import cos, radians
import math
import numpy as np

quarter_sweep_crj = 26.9
wingspan_crj = 23.2 #m
g = 9.81
W_start_cruise_crj = 31537 * g
W_end_cruise_crj = 28337 * g
taper_crj = 0.4
S_crj = 70.6 #crj

quarter_sweep = 25
wingspan = 24
g = 9.81
W_start_cruise = 32000 * g
W_end_cruise = 27000 * g
taper = 0.4
S = 70.6 #crj
wing_twist = -3 #degrees
AR = wingspan**2 /S

altitude = 11000
M = 0.8


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

def dynamic_pressure(M, h):
    T,P,rho,a = atmosphere_calculator(h)
    V = M * a
    return 0.5 * rho * V**2

def beta(M):
    return math.sqrt(1-M**2)

def half_sweep(quarter_sweep, AR,taper):
    return any_sweep(quarter_sweep, AR,taper, 0.25)

def LE_sweep(quarter_sweep,AR,taper):
    return any_sweep(quarter_sweep,AR,taper,0)

def any_sweep(quarter_sweep, AR,taper, chord_position):
    quarter_sweep = radians(quarter_sweep)
    any_sweep = np.tan(quarter_sweep) - (4/AR)*((0.5-chord_position)*(1-taper)/(1+taper))
    return any_sweep

def CL_alpha (AR, M, half_sweep, taper = taper):
    eff = 0.95
    b = beta(M)
    half_sweep = half_sweep(quarter_sweep, AR, taper)
    x = math.sqrt(4+((AR*b/eff)**2)*(1+ (math.tan(radians(half_sweep)))**2)/(b**2))
    return 2*np.pi*AR/(2+x)
   
    
def approxCLmax (clmax, quarter_sweep = quarter_sweep):
    return 0.9 * clmax *radians(quarter_sweep)

def c1_(taper):
    return  -4.2447*taper**4 +12.611*taper**3 - 12.8418 * taper**2 +4.50475*taper

def datcom( quarter_sweep = quarter_sweep, AR = AR, taper = taper):
    c1 = c1_(taper)
    LEsweep = any_sweep(quarter_sweep, AR,taper, 0.25)
    a = 4/((c1+1)*cos(LEsweep))
    if a > AR:
        print("low aspect ratio datcom method")
        datcom = "low"
    elif a < AR:
        print("high aspect ratio datcom method")
        datcom = 'high'
    return datcom, a, AR



# def CLclmax (quarter_sweep, AR, taper, delta_y):
    
          
q = dynamic_pressure(M, altitude)


CL,cl = CLdes(q,W_start_cruise,W_end_cruise,S,quarter_sweep)



crjCL, crjcl = CLdes(q,W_start_cruise_crj,W_end_cruise_crj,S_crj,quarter_sweep_crj)

CL_alpha = CL_alpha (AR, M, half_sweep, taper )
