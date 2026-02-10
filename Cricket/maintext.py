########################
#Author: Cricket Bergner
#Date: 02/10/2026
########################

#making this program again because I did the wrong thing on vscode earlier

#Step 1: Implement trapezoid rule approximation for the given integral. Start
# with one integration slice and work your way up from there. Calculate the 
# the integral to an accuracy of 10^-6. 

#########################################################################

# import libraries
import numpy as np
import math as m
from prettytable import PrettyTable as pt

# initialize variables
condition = N = itl = e = 0
t = pt(["Number of Slices", "Integral Estimate", "Estimate of Error"])

# trapezoid rule
def trap(f, a, b, N):
  dx = (b-a)/N
  s = 0 # sum  
  for i in range(N):
    xi = a + (i*dx)
    xip = a + ((i+1)*dx)
    s += ((f(xi**2 - xi) + f(xip**2 - xip))* (dx / 2))
  return s

# while loop to iterate until conditions of problem are met
while condition < 10:
  N += 2
  condition += 1
  itl = trap(lambda x: np.sin(m.sqrt(100*x))**2, 0, 2, N) # calculate integral
  e =  # calculate the error
  t.add_row([N, round(itl, 6), round(e, 6)]) # add to table
  
print(t) # table



#
#
#

