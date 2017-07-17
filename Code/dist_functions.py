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
def nonstruct_dist(case) :
	
	#~ dist='even'
	#~ dist='equation_30'
	dist=case.W.ns_dist_type
	
	case.W.ntilde=np.zeros((case.wing.m+1,), dtype=np.float64)
	
	if (dist=='function') :
		uf.Eq_30(case)
	elif (dist=='even'):
		for i in range (0,case.wing.m+1) :
			#Eq. (44)
			if (case.W.n_type=='variable'):
				case.W.n=case.W.tot-case.W.s#8.8
			sumodd=0.0
			sumeven=0.0
			k=1
			while k < case.wing.m :
				sumodd=sumodd+ma.sin(case.wing.theta[k])
				k=k+2
			
			k=2
			while k < case.wing.m :
				sumeven=sumeven+ma.sin(case.wing.theta[k])
				k=k+2
			
			case.W.ntilde[i]=-(3.0*np.real(case.wing.m)*(case.W.n-case.W.r))/(case.wing.b*(case.wing.theta[case.wing.m]-case.wing.theta[0])*(ma.sin(case.wing.theta[0])+4.0*sumodd+2.0*sumeven+ma.sin(case.wing.theta[case.wing.m])))
			
	return case.W.ntilde,case.W.n


#######################################################################
#####                                                             #####
#####                       t_c Function                          #####
#####                                                             #####
#######################################################################
#t_c function sets up the non-structural weight distribution	
def t_c(n) :
	
	t_c=np.zeros((n+1,), dtype=np.float64)
	for i in range(0,n+1) :
		t_c[i]=0.12
		
	return t_c 


#######################################################################
#####                                                             #####
#####                       chord Function                        #####
#####                                                             #####
#######################################################################
#chord function sets up the non-sructural weight distribution	
def chord(wing) :
	#~ wing_shape='taper'
	#~ wing_shape='elliptic'
	wing_shape='rectangular'
	chord=np.zeros((wing.m+1,), dtype=np.float64)
	if (wing_shape=='taper') :
		tr=0.7#(2.0*wing.S)/(wing.b*Cr)-1.0
	
		Cr=(2.0*wing.S)/(wing.b*(1.0+tr))
	
		
		for i in range(0,wing.m/2+1) :
			chord[i]=tr*Cr+(2.0*Cr*(1.0-tr))/(wing.b)*(wing.z[i]+wing.b/2.0)
		for i in range(0,wing.m/2+1) :
			chord[wing.m-i]=chord[i]
			
	elif (wing_shape=='elliptic') :
		Ra=wing.b**2/wing.S
		for i in range(0,wing.m/2+1) :
			chord[i]=4.0*wing.b/(ma.pi*Ra)*ma.sqrt(1-(2.0*(wing.z[i])/wing.b)**2)
		for i in range(0,wing.m/2+1) :
			chord[wing.m-i]=chord[i]
			
	elif (wing_shape=='rectangular'):
		for i in range(0,wing.m+1):
			chord[i]=wing.S/wing.b
			#chord[i]=0.22
	#~ for i in range (0,wing.m+1):
		#~ print wing.z[i],chord[i],','
	#~ print ' '
	return chord


#######################################################################
#####                                                             #####
#####                      height Function                        #####
#####                                                             #####
#######################################################################
#chord function sets up the spar height distribution
def height(wing) :
	h=np.zeros((wing.m+1,), dtype=np.float64)
	for i in range (0,wing.m+1):
		h[i]=wing.t_max[i]

	return h
