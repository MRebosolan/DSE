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


W_start_blended_jet = 31537
W_end_blended_jet = 29146
S_blended_jet = 200

AR_blended_jet = wingspan**2/S_blended_jet
taper_blended_jet = 0.2
quarter_sweep_blended_jet = 25

def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad*h
    P = 101300 * (T/288.15)**(-9.81/(T_grad*287))
    rho = P/(T*287)
    a = (1.4*287*T)**0.5
    return (T,P,rho,a)


def CLdes(q,W_start_cruise,W_end_cruise,S,quarter_sweep): #finds design cruise CL for aircraft and airfoil
    CL = 1.1/q * (0.5 * (W_start_cruise/S + W_end_cruise/S))
    cl = CL/(cos(radians(quarter_sweep)))**2
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

def CL_alpha (AR, M, half_sweep, taper):
    eff = 0.95
    b = beta(M)
    half_sweep = half_sweep(quarter_sweep, AR, taper)
    x = math.sqrt(4+((AR*b/eff)**2)*(1+ (math.tan(radians(half_sweep)))**2)/(b**2))
    return 2*np.pi*AR/(2+x)
   
    
def approxCLmax (clmax, quarter_sweep):
    return 0.9 * clmax *radians(quarter_sweep)

def c1_(taper):
    return  -4.2447*taper**4 +12.611*taper**3 - 12.8418 * taper**2 +4.50475*taper

def datcom( quarter_sweep, AR, taper = taper):
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



def CLclmax (quarter_sweep, AR, taper):
    print("look at page 17 in adsee 2 summary, figure 2.3, if you want to fill in a new value")
    return 0.9

          
q = dynamic_pressure(M, altitude)


CL,cl = CLdes(q,W_start_cruise,W_end_cruise,S,quarter_sweep)



crjCL, crjcl = CLdes(q,W_start_cruise_crj,W_end_cruise_crj,S_crj,quarter_sweep_crj)

CL_alpha = CL_alpha (AR, M, half_sweep, taper )



#initial drag estimate:
def reynolds (rho, V, l, mu, k): #subsonic
    std = rho*V*l/mu
    adapted = 38.21*((l/k)**1.053)
    return min(std,adapted)

def laminar_friction(re):
    return 1.328/(math.sqrt(re))

def turbulent (re, M):
    lower = ((math.log(re))**2.58) *(1 + 0.144 * M*M)**0.65
    return 0.455/lower


def formfactor_wing (t_over_c, x_over_c, M, quarter_sweep, AR,taper):
	return (1+ 0.6*t_over_c/x_over_c+ 100*(t_over_c**4)) * (1.34*(M**0.18)*(cos(any_sweep(quarter_sweep, AR,taper, x_over_c))**0.28))
    
def formfactor_fuselage (length_fuselage, diameter_fuselage):
    f = length_fuselage/diameter_fuselage
    return (1+ 60/(f*f*f) + f/400)


def formfactor_nacelle(length_nacelle, diameter_nacelle):
    return 1 + 0.35*length_nacelle/diameter_nacelle
    
def fuselage_misc(upsweep, A_max, A_base,M):
    drag_over_mu_upsweep = 3.83*A_max*upsweep**2.5
    drag_over_mu_base = A_base *(0.139+0.419*(M-0.161)**2)
    return drag_over_mu_upsweep, drag_over_mu_base

def excrescence (CD0, factor): #typically 1.02-1.05 for transport
    return CD0*factor

def tc_streamwise(t_over_c, quarter_sweep):
    return t_over_c*cos(radians(quarter_sweep))

def Mdd (t_over_c, CL, quarter_sweep, kappa = 0.935):
    sweepcos = cos(radians(quarter_sweep))
    tc_stream = tc_streamwise(t_over_c, quarter_sweep)
    
    Mdd = kappa/sweepcos - tc_stream/(sweepcos*sweepcos) - CL / (10 *sweepcos**3)
    return Mdd

def wave(M, t_over_c, CL, quarter_sweep, kappa = 0.935):
    M_dd = Mdd(t_over_c, CL, quarter_sweep, kappa = 0.935)
    if M > M_dd:
        return (M-M_dd)*0.1
    elif M <= M_dd:
        return 0.002
    
# CD0_ = (1/Sref)*Cf_components*formfactor*interference*Swetted + CD_misc


