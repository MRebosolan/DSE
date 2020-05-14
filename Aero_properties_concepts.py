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
t_over_w = 0.4423
thrust = t_over_w*W_start_cruise

altitude = 10000
M = 0.8
<<<<<<< HEAD
<<<<<<< HEAD
Vstall = 63
=======
range = 2000

>>>>>>> 59993b55bd286507ba28154676e27aacb5dd3277
=======
Swet_Sw = 5.6
e = 0.85
cf = 0.003
ke = 0.5*(np.pi*e/cf)**0.5
Emax = ke*(AR/Swet_Sw)**0.5
Cd0 = np.pi*AR*e/(4*Emax*Emax)
=======
>>>>>>> 59993b55bd286507ba28154676e27aacb5dd3277
>>>>>>> 4af1c574d2446a0f73583648cdc495c13c4a16cd
>>>>>>> 8e10be0795901ef5f968ff3cb666b4a6340f22ea
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py


<<<<<<< HEAD
#fuckthepolice

=======
<<<<<<< HEAD

=======
>>>>>>> 8e10be0795901ef5f968ff3cb666b4a6340f22ea
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
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



<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
crjCL, crjcl = CLdes(q,W_start_cruise_crj,W_end_cruise_crj,S_crj,quarter_sweep_crj)

CL_alpha = CL_alpha (AR, M, half_sweep, taper )

#jetfuelcantmeltsteelbeams
Swet_Sw = 5.6
e = 0.85
cf = 0.003
ke = 0.5*(np.pi*e/cf)**0.5
Emax = ke*(AR/Swet_Sw)**0.5
Cd0 = np.pi*AR*e/(4*Emax*Emax)


<<<<<<< HEAD
=======
=======
>>>>>>> 8e10be0795901ef5f968ff3cb666b4a6340f22ea
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
T, P, rho, a = atmosphere_calculator(altitude)
q = 0.5 * rho * a**2 * M**2

def oswald_factor(LE_sweep, AR):
    if LE_sweep < 5:
        return 1.78 * (1 - 0.045 * AR**0.68) - 0.64
    else:
        return 4.61 * (1 - 0.045 * AR**0.68) * (cos(LE_sweep))**0.15 - 3.1
def induced_drag(CL, AR, e):
    return CL**2 / (math.pi * AR * e)
CL = CLdes(q, W_start_cruise, W_end_cruise, S, quarter_sweep)[0]
LEsweep = LE_sweep(quarter_sweep, AR, taper)
e = oswald_factor(LEsweep, AR)
CD_induced = induced_drag(CL, AR, e)
<<<<<<< HEAD

D_tot = q*S*(Cd0+(CL)**2/(np.pi*AR*e))
Cd1= Cd0+(CL)**2/(np.pi*AR*e)
=======
print(CL, AR, e, CD_induced)

T,P,rho,a = atmosphere_calculator(altitude)
V = M * a

def drag(V, rho, Cl, S = S, Cd0 = Cd0, AR = AR, e = e):
    
    CD = Cd0+ induced_drag(Cl, AR, e)
    D = 0.5*rho*V*V*S*CD
    return D, CD

<<<<<<< HEAD

# cruise_energy = D*range*1000 #Joules
# impulse = D*range*1000/V
# cruise_power = D*V
=======
D,CD = drag(V,rho,CL)
cruise_energy = D*range*1000 #Joules
impulse = D*range*1000/V
cruise_power = D*V
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py





<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
# for alt in range(altitude):
#     T,P,rho,a = atmosphere_calculator(altitude)
    
#     drag, CD = drag(V,rho)
#     CL = 2*W_start_cruise / (rho*S*V*V)
start = 50
speeds = np.arange(start, 250, 1)
drags = []
CDS = []
CLS = []
for V in speeds:
    T,P,rho,a = atmosphere_calculator(10000)
    CL = 2*W_start_cruise / (rho*S*V*V)
    D, CD = drag(V,rho, CL)
    drags.append(D)
    CDS.append(CD)
    CLS.append(CL)
    
D_min = min(drags)
V_opt = start + drags.index(min(drags))
CD_opt = CDS[drags.index(min(drags))]
CL_opt = CLS[drags.index(min(drags))]

climbrate = (thrust-D_min)*V_opt/W_start_cruise
climbrate_fpm = climbrate*60/0.3048
    
<<<<<<< HEAD
=======
    thrust_alt = thrust * (P/101325) * (288.15/T)**0.5
    thrust_reduced = thrust_alt *.9
    climbrate = (thrust_reduced-D_min)*V_opt/W_start_cruise
    climbrate_fpm = climbrate*60/0.3048
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
    
    




<<<<<<< HEAD

=======
# def CL_alpha (AR, M, half_sweep, taper = taper):
#     eff = 0.95
#     b = beta(M)
#     half_sweep = half_sweep(quarter_sweep, AR, taper)
#     x = math.sqrt(4+((AR*b/eff)**2)*(1+ (math.tan(radians(half_sweep)))**2)/(b**2))
#     return 2*np.pi*AR/(2+x)
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py


<<<<<<< HEAD
>>>>>>> 59993b55bd286507ba28154676e27aacb5dd3277
=======
<<<<<<< HEAD
# def approxCLmax (clmax, quarter_sweep = quarter_sweep):
#     return 0.9 * clmax *radians(quarter_sweep)

# def c1_(taper):
#     return  -4.2447*taper**4 +12.611*taper**3 - 12.8418 * taper**2 +4.50475*taper

# def MAC (chord_r,taper):

#     return (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
   
    
# def root_chord_tip_chord(S,wingspan, taper):
#     chord_r = (2*S)/(wingspan*(1+taper))
#     chord_t = chord_r * taper
#     return chord_r,chord_t
# def half_sweep(quarter_sweep, AR,taper):
#     return any_sweep(quarter_sweep, AR,taper, 0.25)
# quarter_sweep_crj = 26.9
# wingspan_crj = 23.2 #m
# g = 9.81
# W_start_cruise_crj = 31537 * g
# W_end_cruise_crj = 28337 * g
# taper_crj = 0.4
# S_crj = 70.6 #crj
=======
>>>>>>> 59993b55bd286507ba28154676e27aacb5dd3277
>>>>>>> 4af1c574d2446a0f73583648cdc495c13c4a16cd
>>>>>>> 8e10be0795901ef5f968ff3cb666b4a6340f22ea
>>>>>>> parent of b3f8c9c... Update Aero_properties_concepts.py
