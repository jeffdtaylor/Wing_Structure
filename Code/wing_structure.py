#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import time

#User-Defined Module Containing functions used for
#-Reading Input File
#-Geometry Setup
#-Discretization
#-Moment Calculations
#-Wing Structure Calculations
#-Non-Structural Weight Calculations
#-Solver
import  wing_structure_m as ws

#This program calculates the structural weight distribution of a wing
#with a given lift distribution and non-structural weight distribution.
#The lift distriution should be input as a series of coefficients B_n
#If total non-structural weight is known, the non-structural weight
#distribution should be entered as a .json file. The wing geometry 
#should be also be input from a .json file format.
filename= sys.argv[-1]
#~ filename2 = sys.argv[-1]

#sets up the data structures   
start = time.clock() 
#case1=ws.plane_setup(filename)
case1=ws.Domain(filename,1)
#~ case2=ws.plane_setup(filename2)
elapsed = time.clock()-start
print('-----------------------------------------------------------------------------------------------------')
print(' ')
print( 'Setup Time: ',elapsed)
#Solver computes structural weight
start=time.clock()
case1.solver(1e-16)
#~ case=ws.solver(case2,1e-9)
elapsed = time.clock()-start
print( 'Solver Time: ',elapsed)
print(' ')
results="{0:<25}{1:<20}{2:<10}"

print(results.format('Induced Drag:',case1.D_i,'[N]'))
print(results.format('Structural Weight:',case1.W.s,'[N]'))
print(results.format('Non-Structural Weight:',case1.W.n,'[N]'))
print(results.format('Total Weight:',case1.W.tot,'[N]'))
print(results.format('S_b:',case1.spar.S_b[0],'[m^2]'))
print(results.format('b:',case1.wing.b,'[m]'))
print(results.format('W/S:',case1.W.tot/case1.wing.S,'[m^2]'))
print(results.format('R_n:',case1.W.r/case1.W.tot,' '))
print(results.format('S:',case1.wing.S,'[m^2]'))
#~ print('----------------------------------------------------------------')
#~ print(results.format('Induced Drag:',case2.D_i,'[N]'))
#~ print(results.format('Structural Weight:',case2.W.s,'[N]'))
#~ print(results.format('Non-Structural Weight:',case2.W.n,'[N]'))
#~ print(results.format('Total Weight:',case2.W.tot,'[N]'))
#~ print(results.format('S_b:',case2.spar.S_b[0],'[m^2]'))
#~ print(results.format('b:',case2.wing.b,'[m]'))
#~ print(results.format('W/S:',case2.W.tot/case.wing.S,'[m^2]'))
#~ print(results.format('R_n:',case2.W.r/case.W.tot,' '))
#~ print(results.format('S:',case2.wing.S,'[m^2]'))
#~ for i in range (0,case.wing.m+1):
	#~ print case.wing.z[i]/3.10896,',',case.L.nondim[i]

#~ plt.figure(1)
#~ plt.plot(case.wing.z,case.L.ratio*case.W.tot,'k')
#~ plt.ylabel('Lift Distribution [N/m]')
#~ plt.xlabel('Span [m]')
#~ ax=plt.gca()
#~ ax.set_aspect(.02)
#~ ax.set_ylim([0,55])
#~ ax.set_xlim([0.0,2.0])


plt.figure(1)
plt.plot(case1.wing.z,case1.L.ratio,'k')
plt.ylabel('Structural Weight Dist [N/m]')
plt.xlabel('Span [m]')
ax=plt.gca()
ax.set_aspect(.25)
#ax.set_ylim([0,55])
ax.set_xlim([0.0,2.0])
#~ plt.show()

plt.figure(2)
plt.plot(case1.wing.z,case1.W.ntilde,'k')
plt.ylabel('NS Weight Dist [N/m]')
plt.xlabel('Span [m]')
ax=plt.gca()
ax.set_aspect(.45)
#ax.set_ylim([0,55])
ax.set_xlim([0.0,2.0])
#~ plt.show()
