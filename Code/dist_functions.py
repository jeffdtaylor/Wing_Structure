#user_functions.py contains functions that define the following 
#distributions:
#	-Non-Structural Weight
#	-Thickness-to-chord ratio
#	-Chord Distribution

import numpy as np
import math as ma
import user_functions as uf

#######################################################################
#####                                                             #####
#####                       W_n Function                          #####
#####                                                             #####
#######################################################################
#W_n function sets up the non-sructural weight distribution	
def nonstruct_dist(plane) :
	
	#~ dist='even'
	#~ dist='equation_30'
	dist=plane.W.ns_dist_type
	
	if (dist=='custom') :
		uf.Wns(plane)
	elif (dist=='even'):
		for i in range (0,plane.wing.m+1) :
			#Eq. (44)
			if (plane.W.n_type=='variable'):
				plane.W.n=plane.W.tot-plane.W.s#8.8
			sumodd=0.0
			sumeven=0.0
			k=1
			while k < plane.wing.m :
				sumodd=sumodd+ma.sin(plane.wing.theta[k])
				k=k+2
			
			k=2
			while k < plane.wing.m :
				sumeven=sumeven+ma.sin(plane.wing.theta[k])
				k=k+2
			
			plane.W.ntilde[i]=-(3.0*np.real(plane.wing.m)*(plane.W.n-plane.W.r))/(plane.wing.b*(plane.wing.theta[plane.wing.m]-plane.wing.theta[0])*(ma.sin(plane.wing.theta[0])+4.0*sumodd+2.0*sumeven+ma.sin(plane.wing.theta[plane.wing.m])))


#######################################################################
#####                                                             #####
#####                       t_c Function                          #####
#####                                                             #####
#######################################################################
#t_c function sets up the t/c distribution	
def t_c(wing) :

	if (wing.t_c_type == 'root/tip') :
		for i in range (0,wing.m+1) :
			wing.t_c[i]=wing.root_tc+2.0*(wing.tip_tc-wing.root_tc)/wing.b*wing.z[i]
			
	elif (wing.t_c_type == 'custom'):
		uf.t_c(wing)


#######################################################################
#####                                                             #####
#####                       chord Function                        #####
#####                                                             #####
#######################################################################
#chord function sets up the chord distribution	
def chord(wing) :
	
	wing_shape=wing.c_type
	
	if (wing_shape=='taper') :
		tr=wing.taper_ratio
		Cr=(2.0*wing.S)/(wing.b*(1.0+tr))

		for i in range(0,wing.m+1) :
			wing.c[i]=Cr-((Cr-tr*Cr)/(wing.b/2.0))*wing.z[i]
			
	elif (wing_shape=='elliptic') :
		Ra=wing.b**2/wing.S
		for i in range(0,wing.m+1) :
			wing.c[i]=4.0*wing.b/(ma.pi*Ra)*ma.sqrt(1-(2.0*(wing.z[i])/wing.b)**2)
	
	elif (wing_shape=='rectangular'):
		for i in range(0,wing.m+1):
			wing.c[i]=wing.S/wing.b
			
	elif (wing_shape=='custom'):
		uf.chord(wing)


#######################################################################
#####                                                             #####
#####                      height Function                        #####
#####                                                             #####
#######################################################################
#chord function sets up the spar height distribution
def height(spar,wing) :
	
	if (spar.h_type == 'fill'):
		for i in range (0,wing.m+1):
			spar.h[i]=wing.t_max[i]*spar.fill_ratio
			
	elif (spar.h_type == 'min_fill'):
		minval=np.min(wing.t_max)
		for i in range (0,wing.m+1):
			spar.h[i]=minval
			
	elif (spar.h_type == 'custom'):
		uf.height(spar,wing)
	

