#file for calculating several wing parameters
# authors: Matteo & Jorn

wingspan = 23.94 #m
wingspan_2 = 27.0 #m, value if we use foldable wingtips
def cldesign(q, W_S_start, W_S_end):
    CL_design = 1.1*(1/q)*(0.5*(W_S_start+W_S_end))
    return CL_design

def taper (chord_r, chord_t):
    return chord_t/chord_r

def MAC (taper, chord_r):
    mac = (2/3)*chord_r*((1+taper+taper**2)/(1+taper))
    return mac
