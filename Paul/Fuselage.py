from math import pi,sqrt,atan,degrees



m = 0.0254                  #multiply to convert to m
n = 75                      #number of passengers
seat = 18*m                 #seat width
aisle = 16*m                #aisle width
clearance = 0.05
height = 2.1                #height of aisle
skin_t = 4.7*m              #shell thickness

#horizontal diameter/2
horizontal = (5*seat + aisle + 7*clearance)/2


#max radius either height or width
r = max(horizontal,height) + skin_t

#diameter equation for noncylindrical fuselage
#d = sqrt((2*horizontal+skin_t)*(2*height+skin_t))

area= pi*r**2
circ = 2*pi*r

#min and max fuselage length for noncylindrical fuselage
#min_l_d = 5.6*d
#max_l_d = 10*d

#min and max fuselage length based on roskam
min_l_r = 5.6*r*2
max_l_r = 10*r*2
print(f"D: {round(2*r,1)}m")

print(f"minum fuselage: {round(min_l_r,1)}m")
print(f"maximum fuselage: {round(max_l_r,1)}m")

#baggage and cargo densities
baggage_d = 170
cargo_d = 160
#overhead densities
overhead = 0.065*n
k = 0.35            #percentage of fuselage length used for carge (regional suggest 0.35 or less)
cargo = 2400                #cargo in kilogram
V_cargo = 2400/cargo_d      #cargo in volume
l = V_cargo/(1.3*k)         #fuselage length determined by required cargo size (1.3 determined graphically)

print(f"min length for cargo: {round(l,1)}m")

l_fc_min = 2*2*r        #minimum tail length based on roskam
l_fc_max = 4*2*r        #maximum tail length based on roskam

theta_min = degrees(atan((2*r)/l_fc_max))    #minimum tail angle based on fuselage diameter and maximum tail length
theta_max = degrees(atan((2*r)/l_fc_min))    #maximum tail angle based on fuselage diameter and minimum tail length

print(f"minimum length tail: {round(l_fc_min,1)}m\nmaximum length tail: {round(l_fc_max,1)}m")
print(f"theta min: {round(theta_min,1)}°\ntheta max: {round(theta_max,1)}°")