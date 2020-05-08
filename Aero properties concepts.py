from math import cos, radians

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
        CL = 1.1/1 * (0.5 * (W_start_cruise/S + W_end_cruise/S))
        cl = CL/(cos(radians(sweep)))^2
        self.CL = CL
        self.cl = cl



    def taper (self, chord_r, chord_t):
        self.taper = chord_t/chord_r

    def MAC (self, taper, chord_r):
        MAC = (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
        self.MAC = MAC

