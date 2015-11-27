# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 10:33:37 2015

@author: Fermi
"""

# projectile1: input an initial speed, this code computes the maximum
# range, and outputs the ideal angle for this range.
# last update Nov.20, 2015 by Jason Harlow

# Press "CTRL-S" to save, then "F5" key to run the code.
# Press the "Esc" key to stop the code if it hangs.

#Import fixes to make this Python 2.7 run the code more like Python 3.x
from __future__ import print_function, division

# Import the visual library (vpython).
from visual import *

### Properties of Earth ###
# The acceleration due to gravity [m/s^2]
g = 9.8
# The density of air [kg/m^3]:
rho = 1.2

### Properties of the ball ###
# The mass of the baseball [kg]:
m = 0.145
# The radius of the baseball [m]:
r = 0.037

### degree-step and time-step ###
# The angle-step of the while-loop for choosing different angles [deg]:
dtheta = 0.5
# Set the minimum launch angle to consider [deg]:
thetamin = 20
# The time-step of the numerical integration [s]
dt = 0.01

# Input the initial speed of the ball: [m/s]
v0 = input('Please input the initial speed of the ball in m/s: ')

# Start a while loop that steps through launch angles:
theta0 = thetamin-dtheta
rangeold = -1
rangenew = 0
while rangenew > rangeold:
    # update the rangeold:
    rangeold=rangenew
    # update the launch angle
    theta0 = theta0 + dtheta

    # Set the initial x and y components of the velocity [m/s]
    vx = v0*cos(radians(theta0))
    vy = v0*sin(radians(theta0))
    # Start the ball at the origin:
    x = 0
    y = 0
    # Start an array in case we want to write this trajectory to a file:
    xarray = [x]
    yarray = [x]
    
    #Compute the new range: the x-distance for which the ball returns to y=0:
    #Use the index i to count the number of points in the trajectory:
    i=0
    while y >= 0:
        i=i+1
        
        # Using the velocity from the previous step, update the position:
        x = x + (vx*dt)
        y = y + (vy*dt)

        #Update the array where we are storing this trajectory:
        xarray.append(x)
        yarray.append(y)
        
        # Compute the force on this ball at this instant [N]: 
        fx = 0
        fy = -m*g

        # Compute the acceleration of the ball at this instant [m/s^2]
        ax = fx/m
        ay = fy/m

        # Using the acceleration just calculated, update the velocity:
        vx = vx + (ax*dt)
        vy = vy + (ay*dt)

    # Now that the loop is done, the value of x must be the new range:
    rangenew = x
    npoints=i

# Now to print out the results:
print('\n\n Okay.. all done.  Results of the numerical integration:\n')
print('Ball mass: ', m, ' kg.  Ball diameter: ', (r*200), ' cm.')
print('Time-step of numerical integration: ', (dt*1000), ' ms.')
print('Step-size of launch angle sweep: ', dtheta, ' degrees.')
print('Minimum launch angle considered: ', thetamin, ' degrees.')
print('Initial ball speed: ', v0, ' m/s.')
print('Maximum range computed: ', rangeold, ' m.')
print('The launch angle corresponding to this range: ', theta0-dtheta, ' degrees.')
if rangeold>100:
    print('You got a home run!!')
else:
    print('Sorry, not a home run.. you are OUT!')

# Also, you might want to graph the actual trajectory, so let's make
# a comma-separated text file you can read into excel
file = open('projectileA.csv', 'w')

i=0
while i<npoints:
    i=i+1
    file.write(str(xarray[i]))
    file.write(',')
    file.write(str(yarray[i]))
    file.write('\n')

file.close()
   
