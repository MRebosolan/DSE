# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:13:54 2020

@author: jarin malfliet
"""

n=75 #NUMBER OF PASSENGERS

crew=1#CABIN CREW

#INCHES TO METERS

#ROSKAM
seat_pitch=32*0.0254
seat_width=18*0.0254
arm_rest=2*0.0254
h_aisle=60*0.0254
w_aisle=15*0.0254
headroom=50*0.0254
clearance=1*0.0254
exit_aisle=91*0.01
cross_aisle=51*0.01

emergency_extra=cross_aisle-12*0.0254

#METERS
#GUESSES
l_cockpit_and_nose=3

width333=seat_width*9+w_aisle*2+clearance*2+12*arm_rest

width232=seat_width*7+w_aisle*2+clearance*2+10*arm_rest
width242=seat_width*8+w_aisle*2+clearance*2+11*arm_rest
width141=seat_width*6+w_aisle*2+clearance*2+9*arm_rest

l_tf=3
l_tb=1.1


u2=l_tb/l_tf
v1=width333/l_tf
v2=1.55

#2+4+2

l_rect=63/9*seat_pitch
u1=l_rect/l_tf
l_cabin=l_tf+l_rect+l_tb

wcp=2.7
w_rect=v1*l_tf
w_tb=v2*l_tf


S=(0.5*v1+u1*v1+0.5*u2*(v1+v2))*l_tf**2+0.5*wcp*l_tf


import matplotlib.pyplot as plt


aislex=[l_tf-2,l_tf+l_rect+l_tb]
aisley=[seat_width*1.5+2*arm_rest,seat_width*1.5+2*arm_rest]
plt.plot(aislex,aisley)

aislexu=[l_tf-2,l_tf+l_rect+l_tb]
aisleyu=[seat_width*1.5+2*arm_rest+w_aisle,seat_width*1.5+2*arm_rest+w_aisle]
plt.plot(aislexu,aisleyu)

xlist=[0,0, l_tf, l_tf+l_rect, l_cabin,l_cabin]
ylist=[0,wcp/2, w_rect/2, w_rect/2, w_tb/2,0]
nylist=[0,-wcp/2, -w_rect/2, -w_rect/2, -w_tb/2,0]
plt.plot(xlist,ylist)
plt.plot(xlist,nylist)

passengerlinex=[l_tf-seat_pitch,l_tf-seat_pitch]
passengerliney=[0,wcp/2+(w_rect-wcp)/2/l_tf*(l_tf-seat_pitch)]
plt.plot(passengerlinex,passengerliney)

exitx=[exit_aisle,exit_aisle]
exity=[0,wcp/2]
plt.plot(exitx,exity)


print(l_cabin, l_cockpit_and_nose, l_rect)
plt.axis('equal')
plt.ylabel('some numbers')
plt.show()