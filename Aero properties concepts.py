from math import cos, radians

def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad*h
    P = 101300 * (T/288.15)**(-9.81/(T_grad*287))
    rho = P/(T*287)
    a = (1.4*287*T)**0.5
    return (T,P,rho,a)


class wing:

    def __init__(self, sweep, W_start_cruise, W_end_cruise, S, chord_t, chord_r, wingspan):
        self.sweep = sweep
        self.W_start_cruise = W_start_cruise
        self.W_end_cruise = W_end_cruise
        self.S = S
        self.chord_t = chord_t
        self.chord_r = chord_r
        self.wingspan = wingspan


    def CLdes(self, q, W_start_cruise, W_end_cruise, S, sweep): #finds design cruise CL for aircraft and airfoil
        CL = 1.1/q * (0.5 * (W_start_cruise/S + W_end_cruise/S))
        cl = CL/(cos(radians(sweep)))^2
        self.CL = CL
        self.cl = cl



    def taper (self, chord_r, chord_t):
        self.taper = chord_t/chord_r

    def MAC (self, taper, chord_r):
        MAC = (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
        self.MAC = MAC

    def root_chord_tip_chord(self, S, wingspan, taper):
        self.chord_r = (2*S)/(wingspan*(1+taper))
        self.chord_t = self.chord_r * taper
