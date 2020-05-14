from math import cos, radians
import math
import numpy as np
import matplotlib.pyplot as plt



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


range = 2000
Swet_Sw = 5.6


def aerodynamics(quarter_sweep, wingspan, W_start_cruise, taper,
                 S, t_over_w, AR, thrust, altitude, M, range, Swet_Sw):

    def atmosphere_calculator(h):
        T_grad = -0.0065
        T = 288.15 + T_grad*h
        P = 101325 * (T/288.15)**(-9.81/(T_grad*287))
        rho = P/(T*287)
        a = (1.4*287*T)**0.5
        return (T,P,rho,a)
    
    
    def CLdes(q,W_start_cruise,W_end_cruise,S,sweep): #finds design cruise CL for aircraft and airfoil
        CL = 1.1/q * (0.5 * (W_start_cruise/S + W_end_cruise/S))
        cl = CL/(cos(radians(sweep)))**2
        return CL, cl
    
    
    def LE_sweep(quarter_sweep,AR,taper):
        return any_sweep(quarter_sweep,AR,taper,0)
    
    def any_sweep(quarter_sweep, AR,taper, chord_position):
        quarter_sweep = radians(quarter_sweep)
        any_sweep = np.tan(quarter_sweep) - (4/AR)*((0.5-chord_position)*(1-taper)/(1+taper))
        return any_sweep
    
    
    
    
    def oswald_factor(LE_sweep, AR):
        if LE_sweep < 5:
            return 1.78 * (1 - 0.045 * AR**0.68) - 0.64
        else:
            return 4.61 * (1 - 0.045 * AR**0.68) * (cos(LE_sweep))**0.15 - 3.1
    
    def induced_drag(CL, AR = AR, taper = taper, quarter_sweep=quarter_sweep):
        LEsweep = LE_sweep(quarter_sweep, AR, taper)
        e = oswald_factor(LEsweep, AR)
        return CL**2 / (math.pi * AR * e)
    
    
    def CDZERO():
        e = oswald_factor(LE_sweep(quarter_sweep, AR, taper),AR)
        cf = 0.003
        ke = 0.5*(np.pi*e/cf)**0.5
        Emax = ke*(AR/Swet_Sw)**0.5
        Cd0 = np.pi*AR*e/(4*Emax*Emax)
        return Cd0
    
    def drag(V, rho, Cl, S = S, AR = AR):
        Cd0 = CDZERO()
        CD = Cd0 + induced_drag(Cl, AR)
        D = 0.5*rho*V*V*S*CD
        return D, CD
    
    

    
    def cruise():
        T, P, rho, a = atmosphere_calculator(altitude)
        V = M * a
        CL = CLdes(0.5*rho*V*V, W_start_cruise, W_end_cruise, S, quarter_sweep)[0]
        D,CD = drag(V,rho,CL)
        cruise_energy = D*range*1000 #Joules
        impulse = D*range*1000/V
        cruise_power = D*V
        print("Cruise energy: ", cruise_energy)
        print("Cruise impulse: ",impulse)
        return cruise_energy, impulse
    
    def climb(altitude = altitude):
        alt = 0
        t = 0
        impulse = 0
        altlist = []
        thrustlist = []
        v_list = []
        drag_list = []
        climb_list = []
        powerlist =[]
        energy = 0
        
        while alt < altitude and t < 1800:
            
            start = 80
            speeds = np.arange(start, 220, 1)
            drags = []
            CDS = []
            CLS = []
            for V in speeds:
                T,P,rho,a = atmosphere_calculator(alt)
                CL = 2*W_start_cruise / (rho*S*V*V)
                D, CD = drag(V,rho, CL)
                drags.append(D)
                CDS.append(CD)
                CLS.append(CL)
                
            D_min = min(drags)
            V_opt = start + drags.index(min(drags))
            CD_opt = CDS[drags.index(min(drags))]
            CL_opt = CLS[drags.index(min(drags))]
            
            thrust_alt = thrust * (P/101325) * (288.15/T)**0.5
            thrust_reduced = thrust_alt *.9
            climbrate = (thrust_reduced-D_min)*V_opt/W_start_cruise
            climbrate_fpm = climbrate*60/0.3048
            
            t = t + 1
            alt = alt + climbrate
            impulse = impulse + thrust_reduced
            
            power = thrust_reduced*V_opt
            powerlist.append(power)
            energy = energy + power
            altlist.append(alt)
            thrustlist.append(thrust_reduced)
            v_list.append(V_opt)
            drag_list.append(D_min)
            climb_list.append(climbrate)
        
        time = np.arange(0,len(altlist))
        
        #plt.plot(time, altlist)
        plt.plot(altlist,np.array(powerlist)/1000)
        plt.grid()
        print("Climb energy: ", energy)
        print("Climb impulse: ",impulse)
        
        return energy, impulse
    
    cr = cruise()
    cli = climb()    
    print(cr[1]/cli[1])
    
    total_energy = cr[0]+cli[0]
    total_impulse = cr=[1]+cli[1]
    
    return total_energy,total_impulse

aero = aerodynamics(quarter_sweep, wingspan, W_start_cruise, taper,
                 S, t_over_w, AR, thrust, altitude, M, range, Swet_Sw)

print(aero)






# def CL_alpha (AR, M, half_sweep, taper = taper):
#     eff = 0.95
#     b = beta(M)
#     half_sweep = half_sweep(quarter_sweep, AR, taper)
#     x = math.sqrt(4+((AR*b/eff)**2)*(1+ (math.tan(radians(half_sweep)))**2)/(b**2))
#     return 2*np.pi*AR/(2+x)

# def datcom( quarter_sweep = quarter_sweep, AR = AR, taper = taper):
#     c1 = c1_(taper)
#     LEsweep = any_sweep(quarter_sweep, AR,taper, 0.25)
#     a = 4/((c1+1)*cos(LEsweep))
#     if a > AR:
#         print("low aspect ratio datcom method")
#         datcom = "low"
#     elif a < AR:
#         print("high aspect ratio datcom method")
#         datcom = 'high'
#     return datcom, a, AR

    # def beta(M):
    #     return math.sqrt(1-M**2)
    


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

