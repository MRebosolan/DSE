#file for calculating several wing parameters
# authors: Matteo & Jorn

def cldesign(q, W_S_start, W_S_end):
    CL_design = 1.1*(1/q)*(0.5*(W_S_start+W_S_end))
    return CL_design