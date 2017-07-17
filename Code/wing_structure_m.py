#Wing_structure_m.py contains functions used for
#	-Reading Input File
#	-Geometry Setup
#	-Discretization
#	-Moment Calculations
#	-Wing Structure Calculations
#	-Non-Structural Weight Calculations
#	-Solver

import sys
import os
import numpy as np
import math as ma
import json
from collections import OrderedDict
import time
import integrators as ints
import dist_functions as df

class Domain(object):
	_slots_=['wing','spar','M','W','n','f','L','D_i']
	
	def __init__(self,filename):

		self.wing=self.Wing()
		self.spar=self.Spar()
		self.M=self.Moments()
		self.W=self.Weight()
		self.n=self.Limits()
		self.f=self.Flight()
		self.L=self.Lift_Distribution()
		
		#initialize data_structure
		self.init(filename)

		#discretize wing
		self.wing.discretize()	
		
		#determine remaining geometry
		self.geom()

		#Find Lift_Distribution
		self.lift_dist()
				
	def unpackplane(self):
		#######################################################################
		#####                                                             #####
		#####                   unpackplane Function                       #####
		#####                                                             #####
		#######################################################################
		#unpackplane function divides up the class plane into individual classes
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-wing=data structure containing wing geometry information
		#	-spar=data structure containing spar geometry information
		#	-M=data structure containing moment information
		#	-W=data structure containing weight information
		#	-n=data structure containing load factor information
		#	-f=data structure containing flight conditions
		#	-L=data structure containing lift distribution information
		#----------------------------------------------------------------------
		wing=self.wing
		spar=self.spar
		M=self.M
		W=self.W
		n=self.n
		f=self.f
		L=self.L
		
		return wing,spar,M,W,n,f,L	
			
	def packplane(self,wing,spar,M,W,n,f,L):
		#######################################################################
		#####                                                             #####
		#####                     packplane Function                       #####
		#####                                                             #####
		#######################################################################
		#packplane function combines all classes into plane class
		#----------------------------------------------------------------------
		#inputs:
		#	-wing=data structure containing wing geometry information
		#	-spar=data structure containing spar geometry information
		#	-M=data structure containing moment information
		#	-W=data structure containing weight information
		#	-n=data structure containing load factor information
		#	-f=data structure containing flight conditions
		#	-L=data structure containing lift distribution information
		#	-plane=data structure containing compiled problem scenario
		#outputs:
		#	-plane=data structure containing compiled problem scenario
		#----------------------------------------------------------------------
		self.wing=wing
		self.spar=spar
		self.M=M
		self.W=W
		self.n=n
		self.f=f
		self.L=L
					
	def init(self,filename):
		#######################################################################
		#####                                                             #####
		#####                       Init Function                         #####
		#####                                                             #####
		#######################################################################
		#init function initializes data structures with values from input file
		#----------------------------------------------------------------------
		#inputs:
		#	-filename=file containing aircraft geometry
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-
		#----------------------------------------------------------------------
		with open(filename) as file:
			data = json.load(file, object_pairs_hook = OrderedDict)
		
		#Initialize data from file
		self.file_values(data)
		
		#initialize distributions according to specification in input file
		self.wing.c,self.wing.c_type=distributions_init(filename,data["wing"]["chord"],self.wing.m)	
		self.wing.t_c,self.wing.t_c_type=distributions_init(filename,data["wing"]["thickness_chord"],self.wing.m)
		self.W.ntilde,self.W.ns_dist_type=distributions_init(filename,data["weight"]["nonstructural_distribution"],self.wing.m)
		self.spar.h,self.spar.h_type=distributions_init(filename,data["spar"]["height"],self.wing.m)
		
		#initialize weight
		self.W.init_weight(data)
		
		#initialize beam
		self.spar.init_beam(data)
		
		#set remaining values to zero
		self.init_zeros()

	def file_values(self,data):
		#######################################################################
		#####                                                             #####
		#####                  file_values Function                       #####
		#####                                                             #####
		#######################################################################
		#file_values function reads in non-zero values from the input file
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=struct containing plane setup information
		#	-data=input file data loaded from json.load
		#outputs:
		#	-
		#----------------------------------------------------------------------

		self.wing.m=data["wing"]["grid"]
		self.wing.S=data["wing"]["wing_area"]
		self.wing.b=data["wing"]["wing_span"]
		self.spar.sigma_max=data["spar"]["max_stress"]
		self.spar.gam=data["spar"]["specific_weight"]
		self.n.m=data["limits"]["maneuvering"]
		self.n.g=data["limits"]["hard_landing"]
		self.f.rho=data["flight"]["density"]
		self.f.V=data["flight"]["velocity"]
		self.L.B=data["lift_distribution"]["B"]
		self.spar.beam_type=data["spar"]["beam_type"]
		self.W.init_weight(data)
		if 'loading' in data["wing"]:
			self.wing.loading=data["wing"]["loading"]	
		else:
			self.wing.loading=self.W.tot/self.wing.S

	def init_zeros(self):
		#######################################################################
		#####                                                             #####
		#####                   init_zeros Function                       #####
		#####                                                             #####
		#######################################################################
		#init_zeros function initializes the all undefined values to zero
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=struct containing plane setup information

		#outputs:
		#	-
		#----------------------------------------------------------------------

		self.wing.t_max=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.wing.theta=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.wing.z=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.spar.C_sigma=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.spar.S_b=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.M.bnm=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.M.bng=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.M.b=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.W.stilde=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.L.tilde=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.L.ratio=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.L.nondim=np.zeros((self.wing.m+1,), dtype=np.float64)
		self.D_i=0.0	
			
	def geom(self):
		#######################################################################
		#####                                                             #####
		#####                       Geom Function                         #####
		#####                                                             #####
		#######################################################################	
		#Geom function calculates remaining geometry values
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-plane=data structure containing problem scenario
		#----------------------------------------------------------------------

		self.wing.set_wing_geom()
		
		self.spar.set_spar_geom(self.wing)
				
	def lift_dist(self) :
		#######################################################################
		#####                                                             #####
		#####                      lift_dist Function                     #####
		#####                                                             #####
		#######################################################################
		#lift_dist function returns the lift distribution
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-L=data structure containing lift distribution information
		#----------------------------------------------------------------------		
		wing,spar,M,W,n,f,L=self.unpackplane()
		
		for i in range (0,wing.m+1) :
			Bsinsum=0.0
			for k in range (0,len(L.B)) :
				Bsinsum=Bsinsum+L.B[k]*ma.sin((k+2)*wing.theta[i])
			
			#Nondimensional Lift Distribution Eq. (39)
			L.nondim[i]=(4.0/ma.pi)*(ma.sin(wing.theta[i])+Bsinsum)
			
			#Nondimensional Lift Distribution divided by Wingspan
			L.ratio[i]=L.nondim[i]/wing.b
			
		self.packplane(wing,spar,M,W,n,f,L)

	def bending_moment(self) :
		#######################################################################
		#####                                                             #####
		#####                  bending_moment Function                    #####
		#####                                                             #####
		#######################################################################
		#bending_moment function calculates the bending moments on the wing
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-M=data structure containing moment information
		#----------------------------------------------------------------------
		
		#if (plane.W.r>plane.W.tot*(plane.n.g-1.0)/(plane.n.g+plane.n.m)):
			#Maneuvering-Flight Bending moment distribution
		self.bending_moment_nm()
		#else:
			#Hard-Landing Bending moment distribution
		self.bending_moment_ng()	

	def bending_moment_nm(self) :
		#######################################################################
		#####                                                             #####
		#####                bending_moment_nm Function                   #####
		#####                                                             #####
		#######################################################################	
		#returns maneuvering-flight bending moments	
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-M=data structure containing moment information
		#----------------------------------------------------------------------
		wing,spar,M,W,n,f,L=self.unpackplane()
		
		ft=np.zeros((wing.m+1,), dtype=np.float64)
		
		#Maneuvering-Flight Bending moment Eq. (40)	
		for i in range (0,wing.m+1):
			#f(theta_k) Eq. (40)
			for k in range (i,wing.m+1) :
				ft[k]=(-wing.b*ma.sin(wing.theta[k])/2.0)*(W.tot*L.ratio[k]-W.ntilde[k]-W.stilde[k])*(ma.cos(wing.theta[k])-ma.cos(wing.theta[i]))
			#simpson's 1/3 rule Eq. (40)
			if (i != wing.m) :
				M.bnm[i]=n.m*wing.b/2.0*ints.simpson_third(i,wing.m,wing.theta,ft)
				#M.bnm[wing.m-i]=M.bnm[i]
			else :
				M.bnm[i]=0.0
				#M.bnm[wing.m-i]=M.bnm[i]
				
		self.packplane(wing,spar,M,W,n,f,L)	

	def bending_moment_ng(self) :
		#######################################################################
		#####                                                             #####
		#####                bending_moment_ng Function                   #####
		#####                                                             #####
		#######################################################################	
		#returns hard_landing bending moments	
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-M=data structure containing moment information
		#--------------------------------------------------------------------
		wing,spar,M,W,n,f,L=self.unpackplane()

		ft=np.zeros((wing.m+1,), dtype=np.float64)
		
		#Hard-Landing Bending moment Eq. (41)		
		for i in range (0,wing.m+1):
			#f(theta_k) Eq. (41)
			for k in range (i,wing.m+1) :
				ft[k]=(-wing.b*ma.sin(wing.theta[k])/2.0)*(W.ntilde[k]+W.stilde[k]-W.tot/n.g*L.ratio[k])*(ma.cos(wing.theta[k])-ma.cos(wing.theta[i]))
			#simpson's 1/3 rule Eq. (41)
			if (i != wing.m) :
				M.bng[i]=-n.g*wing.b/2.0*ints.simpson_third(i,wing.m,wing.theta,ft)
				M.bng[wing.m-i]=M.bng[i]
			else :
				M.bng[i]=0.0
				M.bng[wing.m-i]=M.bng[i]		
				
		self.packplane(wing,spar,M,W,n,f,L)
		
	def structure(self) :
		#######################################################################
		#####                                                             #####
		#####                    structure Function                       #####
		#####                                                             #####
		#######################################################################	
		#structure function finds the wing structure required to support the 
		#bending moments
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-W=data structure containing weight information
		#----------------------------------------------------------------------
		wing,spar,M,W,n,f,L=self.unpackplane()

		#Structural Weight Distribution Eq. (20)
		for i in range (0,wing.m+1) :
			if (wing.c[i]==0):
				W.stilde[i]=0.0
			else:
				W.stilde[i]=max(np.absolute(M.bnm[i])/spar.S_b[i],np.absolute(M.bng[i])/spar.S_b[i])
				#W.stilde[i]=np.absolute(M.bnm[i])/spar.S_b[i]
		#Total Structural Weight Eq. (42)
		ft=np.zeros((wing.m+1,), dtype=np.float64)
		for i in range (0,wing.m+1) :
			ft[i]=W.stilde[i]*ma.sin(wing.theta[i])
			
		W.s=-wing.b*ints.simpson_third(0,wing.m,wing.theta,ft)
		
		self.packplane(wing,spar,M,W,n,f,L)

	def nonstructure(self) :
		#######################################################################
		#####                                                             #####
		#####                   nonstructure Function                     #####
		#####                                                             #####
		#######################################################################
		#nonstructure function determines the nonstructural weight distribution
		#and total nonstructural weight
		#----------------------------------------------------------------------
		#inputs:
		#	-plane=data structure containing problem scenario
		#outputs:
		#	-W=data structure containing weight information
		#---------------------------------------------------------------------- 
		
		wing,spar,M,W,n,f,L=self.unpackplane()
		
		#Non-Structural Weight Distribution
		if ((W.ns_dist_type=='function') or (W.ns_dist_type=='even')):  
			W.ntilde,W.n=df.nonstruct_dist(self)
		else:
			W.ntilde=W.ntilde
		#Total Non-structural Weight Eq. (38)	
		ft=np.zeros((wing.m+1,), dtype=np.float64)
		for i in range (0,wing.m+1):
			ft[i]=W.ntilde[i]*ma.sin(wing.theta[i])
		
		W.n=-wing.b*ints.simpson_third(0,wing.m,wing.theta,ft)+W.r

		self.packplane(wing,spar,M,W,n,f,L)
		
		
	#Holds the wing geometry information
	class Wing(object):
		_slots_=['S','b','c','c_type','t_c','t_c_type','t_max','m','z','theta','loading']
		
		def set_t_c(self) :
			#######################################################################
			#####                                                             #####
			#####                       set_t_c Function                      #####
			#####                                                             #####
			#######################################################################	
			#set_chord function scales the chord to keep wing area constant
			#----------------------------------------------------------------------
			#inputs:
			#	-input_file_entry=.json entry containing distribution specification
			#					  (either 'function' or 'filename')
			#	-n=number of nodes
			#outputs:
			#	-t_c=array containing the thickness-to-chord ratio at each spanwise 
			#		 section
			#	-t_c_type=either 'function' or 'file'
			#----------------------------------------------------------------------
			if (self.t_c_type=='function'):
				self.t_c=df.t_c(self.m)	

		def wing_area(self) :
			#######################################################################
			#####                                                             #####
			#####                   wing_area Function                        #####
			#####                                                             #####
			#######################################################################
			#wing_area function finds the wing area from the chord distribution
			#----------------------------------------------------------------------
			#inputs:
			#	-wing=data_structure containing the wing information
			#outputs:
			#	-S=wing area determined from chord distribution
			#----------------------------------------------------------------------
			#Integrate chord in theta using simpson's rule
			ft=np.zeros((self.m+1,), dtype=np.float64)
			for i in range (0,self.m+1):
				ft[i]=self.c[i]*ma.sin(self.theta[i])
				
			self.S=-self.b*0.5*ints.simpson_third(0,self.m,self.theta,ft)

			return S
			
		def set_wing_geom(self):
			#######################################################################
			#####                                                             #####
			#####                  set_wing_geom Function                     #####
			#####                                                             #####
			#######################################################################	
			#set_wing_geom calculates wing geometry values
			#----------------------------------------------------------------------
			#inputs:
			#	-wing=struct containing wing geometry
			#outputs:
			#	-wing=struct containing wing geometryo
			#----------------------------------------------------------------------

			self.set_chord()
			
			self.set_t_c()

			#Set wing max thickness according to thickness-to-chord ratio
			#distribution and chord distribution
			for i in range (0,self.m+1):
				self.t_max[i]=self.t_c[i]*self.c[i]	
							
		def discretize(self) :
			#######################################################################
			#####                                                             #####
			#####                      discretize Function                    #####
			#####                                                             #####
			#######################################################################
			#discretize function discretizes the wing over the full span
			#----------------------------------------------------------------------
			#inputs:
			#	-wing=data structure containing wing information
			#outputs:
			#	-wing=data structure containing wing information
			#----------------------------------------------------------------------
			#Cosine Clustering
			for i in range (0,self.m+1) :
				self.theta[i]=ma.pi/2.0-ma.pi/2.0*(i)/float(self.m)
				self.z[i]=self.b/2.0*ma.cos(self.theta[i])

		def set_chord(self) :
			#######################################################################
			#####                                                             #####
			#####                      set_chord Function                     #####
			#####                                                             #####
			#######################################################################	
			#set_chord function scales the chord to keep wing area constant
			#----------------------------------------------------------------------
			#inputs:
			#	-wing=data structure containing wing information
			#	-b_old=wingspan before optimization
			#	-b_new=wingspan after optimization
			#	Note: if optimization is not performed, b_old=b_new-plane.wing.b
			#outputs:
			#	-c=array containing chord length at each spanwise section
			#----------------------------------------------------------------------
			#If chord is specified in function, use the function to find the 
			#new chord to maintain constant wing area
			
			if (self.c_type=='function'):
				self.c=df.chord(self)	
			
					
	#Holds the spar geometry information	
	class Spar(object):
		_slots_=['sigma_max','h','h_type','gam','C_sigma','S_b','beam_type','inner_width','inner_height','outer_width','flange_height','flange_width','web_width']

		def init_beam(self,data):
			#######################################################################
			#####                                                             #####
			#####                   init_beam Function                        #####
			#####                                                             #####
			#######################################################################
			#init_beam function initializes beam properties based on initial 
			#information
			#----------------------------------------------------------------------
			#inputs:
			#	-plane=data structure with the probelm scenario
			#	-data=input file data loaded from json.load

			#outputs:
			#	-
			#----------------------------------------------------------------------
			if (self.beam_type=='box'):
				self.inner_height=data["spar"]["inner_height"]
				self.inner_width=data["spar"]["inner_width"]
				self.outer_width=data["spar"]["outer_width"]
			elif (self.beam_type=='I'):
				self.flange_height=data["spar"]["flange_height"]
				self.web_width=data["spar"]["web_width"]
				self.flange_width=data["spar"]["flange_width"]	

		def set_spar_geom(self,wing):
			#######################################################################
			#####                                                             #####
			#####                  set_spar_geom Function                     #####
			#####                                                             #####
			#######################################################################	
			#set_spar_geom calculates spar geometry values
			#----------------------------------------------------------------------
			#inputs:
			#	-spar=struct containing spar geometry
			#	-wing=struct containing wing geometry
			#outputs:
			#	-spar=struct containing spar geometry
			#----------------------------------------------------------------------
			#set or scale spar height
			for i in range (0,wing.m+1):	
				self.set_height(wing)
				
			#Calculate shape factor according to beam type	
			if (self.beam_type == 'rectangular'):
				self.rectangular_beam(wing)
			elif (self.beam_type == "box"):
				self.box_beam(wing)
			elif (self.beam_type == 'I'):
				self.I_beam(wing)
				
			#Calculate proportionality Constant Eq. (15)
			for i in range (0,wing.m+1):
				if (wing.c[i]==0):
					self.S_b[i]=0.0
				else:
					self.S_b[i]=(self.C_sigma[i]*(wing.t_c[i])*wing.c[i]*self.sigma_max)/(self.gam)
					
		def rectangular_beam(self,wing) :
			#######################################################################
			#####                                                             #####
			#####                 rectangular_beam Function                   #####
			#####                                                             #####
			#######################################################################
			#calculates shape factor for rectangular beam
			#----------------------------------------------------------------------
			#inputs:
			#	-spar=data structure containing spar information
			#	-wing=data structure containing wing information
			#outputs:
			#	-C_sigma=array of shape factors at each spanwise section
			#----------------------------------------------------------------------
			self.C_sigma=np.zeros((wing.m+1,), dtype=np.float64)
			#Shape factor for rectangular beam Eq. (14)
			for i in range (0,wing.m+1):
				if (wing.c[i]==0):
					self.C_sigma[i]=0.0
				else:
					#self.C_sigma[i]=(self.h[i]/wing.t_max[i])/6.0
					self.C_sigma[i]=0.165
			
		def box_beam(self,wing) :
			#######################################################################
			#####                                                             #####
			#####                     box_beam Function                       #####
			#####                                                             #####
			#######################################################################
			#calculates shape factor for box beam
			#----------------------------------------------------------------------
			#inputs:
			#	-data=json information from input file containing beam information
			#	-spar=data structure containing spar information
			#	-wing=data structure containing wing information
			#outputs:
			#	-C_sigma=array of shape factors at each spanwise section
			#----------------------------------------------------------------------
			C_sigma=np.zeros((wing.m+1,), dtype=np.float64)
			h_i=np.zeros((wing.m+1,), dtype=np.float64)
			w_i=np.zeros((wing.m+1,), dtype=np.float64)
			w=np.zeros((wing.m+1,), dtype=np.float64)
			#Set beam dimensions to values from input file (SCALE??)
			h_i=self.inner_height#/self.h[wing.m/2])*self.h[i]
			w_i=self.inner_width#*h_i[i]/self.h[wing.m/2])
			w=self.outer_width#*h_i[i]/self.h[wing.m/2])
			h=self.h
			ratio=((1-w_i*h_i**3/(w*h[0]**3)))/(6.0*(1-w_i*h_i/(w*h[0])))
			
			#Shape factor for box beam Eq. (14)
			for i in range (0,wing.m+1):
				if (wing.c[i]==0):
					self.C_sigma[i]=0.0
				else:
					self.C_sigma[i]=ratio*(h[i]/wing.t_max[i])

		def I_beam(self,wing) :
			#######################################################################
			#####                                                             #####
			#####                      I_beam Function                        #####
			#####                                                             #####
			#######################################################################
			#calculates shape factor for box beam
			#----------------------------------------------------------------------
			#inputs:
			#	-data=json information from input file containing beam information
			#	-spar=data structure containing spar information
			#	-wing=data structure containing wing information
			#outputs:
			#	-C_sigma=array of shape factors at each spanwise section
			#----------------------------------------------------------------------
			C_sigma=np.zeros((wing.m+1,), dtype=np.float64)
			h_f=np.zeros((wing.m+1,), dtype=np.float64)
			w_w=np.zeros((wing.m+1,), dtype=np.float64)
			w=np.zeros((wing.m+1,), dtype=np.float64)
			ratio=np.zeros((wing.m+1,), dtype=np.float64)
			#Set beam dimensions to values from input file (SCALE??)
			h_f=self.flange_height#*(self.h[i]/self.h[wing.m/2])
			w_w=self.web_width#*(self.h[i]/self.h[wing.m/2])
			w=self.flange_width#*(self.h[i]/self.h[wing.m/2])
			h=self.h
			t_max=wing.t_max
			ratio=((2.0*(h_f/h[0])**3+6.0*(h_f/h[0])*(1-h_f/h[0])**2+(w_w/w)*(1-2.0*h_f/h[0])**3))/(6.0*(2.0*h_f/h[0]+(w_w/w)*(1-2.0*h_f/h[0])))
			
			#Shape factor for I beam Eq. (14)
			for i in range (0,wing.m+1):
				if (wing.c[i]==0):
					self.C_sigma[i]=0.0
				else:
					self.C_sigma[i]=ratio*(h[i]/t_max[i])	

		def set_height(self,wing) :
			#######################################################################
			#####                                                             #####
			#####                   set_height Function                       #####
			#####                                                             #####
			#######################################################################	
			#set_height function sets the spar height
			#----------------------------------------------------------------------
			#inputs:
			#	-wing=data structure containing wing information
			#	-spar=data structure containing spar information
			#	-t_max_old=t_max before chord change
			#	Note: if optimization is not performed, b_old=b_new-plane.wing.b
			#outputs:
			#	-h=array containing spar height at each spanwise section
			#----------------------------------------------------------------------
			#If height is specified in function, use the function to find the 
			#new height
			
			if (self.h_type=='function'):
				self.h=df.height(wing)
			
			#If height is specified in a file, scale the height
			else:
				for i in range (0,wing.m+1):
					self.h[i]=self.h[i]					

									
	#Holds the moment information	
	class Moments(object):
		_slots_=['bnm','bng','b']


	#Holds the weight information			
	class Weight(object):
		_slots_=['tot','r','n','n_type','s','R_n','r_type','ntilde','stilde','ns_dist_type']
		
		def init_weight(self,data):
			#######################################################################
			#####                                                             #####
			#####                  init_weight Function                       #####
			#####                                                             #####
			#######################################################################
			#init_weight function initializes weights based on initial information
			#----------------------------------------------------------------------
			#inputs:
			#	-W=data structure with the weight information
			#	-data=input file data loaded from json.load

			#outputs:
			#	-W
			#----------------------------------------------------------------------
			self.tot=0.0
			self.r=0.0
			self.s=0.0
			self.n=0.0

			#initialize Non-structural Weight
			if ('nonstructural_weight' in data["weight"]):
				self.n=data["weight"]["nonstructural_weight"]
				self.n_type='constant'
			else:
				self.n_type='variable'
				
			#initialize root weight	
			if ('root_total' in data["weight"]):
				self.R_n=data["weight"]["root_total"]
				self.r_type='ratio'
			else:
				self.r=data["weight"]["root_weight"]
				self.r_type='constant'

			#initialize total weight
			if ('total_weight' in data["weight"]):
				self.tot=data["weight"]["total_weight"]
			else:
				self.tot=self.n+self.r+self.s
			
			
	#Holds the load limits information	
	class Limits(object):
		_slots_=['m','g']


	#Holds the flight condition information	
	class Flight(object):
		_slots_=['rho','V']
		
		
	#Holds the lift distribution information
	class Lift_Distribution(object):
		_slots_=['B','tilde','ratio','nondim']

########################################################################
########################################################################
########################################################################
		
def setup():
#######################################################################
#####                                                             #####
#####                      Setup Function                         #####
#####                                                             #####
#######################################################################
#setup function sets up the test case.
#----------------------------------------------------------------------
#inputs:
#	-filename=file containing aircraft geometry
#outputs:
#	-plane=data structure containing initialized problem scenario
#----------------------------------------------------------------------
	#call function to setup data structure
	plane=plane_setup()
	

	
	return plane


def solver(plane,convergence) :
#######################################################################
#####                                                             #####
#####                      Solver Function   	                  #####
#####                                                             #####
#######################################################################
#solver function runs through the algorithm until a solution converges
#----------------------------------------------------------------------
#inputs:
#	-plane=data structure containing problem scenario
#	-convergence=convergence criterion. Iterations will continue until 
#	             the error in structural weight is below convergence 
#outputs:
#	-plane=data structure containing problem scenario
#----------------------------------------------------------------------
	W_S=plane.wing.loading
	error=1.0
	while error>convergence :
		
		if (plane.W.r_type=='ratio'):
			plane.W.r=plane.W.R_n*plane.W.tot
		plane.wing.S=plane.W.tot/W_S	
		plane.geom()
		
		
		prev=plane.W.s
			
		#Step (2) [3 new]
		plane.bending_moment()
		
		#Step (4)
		plane.structure()
		
		#Step (4b)
		plane.nonstructure()
		curr=plane.W.s
		
		#error=error+1
		if (prev==0.0):
			error=1.0
		else:
			error=np.absolute(((curr-prev)/prev)*100)
			
		#Step (1) [2 new]
		plane.W.tot=plane.W.s+plane.W.n

	#Total Weight Eq. (4)
	plane.W.tot=plane.W.s+plane.W.n

	#Induced Drag Eq. (3)
	Bsum=0.0
	for i in range (0,len(plane.L.B)) :
			Bsum=Bsum+(i+2)*plane.L.B[i]**2
	plane.D_i=(2.0*(plane.W.tot/plane.wing.b)**2)/(ma.pi*plane.f.rho*plane.f.V**2)*(1.0+Bsum)
	return plane	


def plane_setup(filename):
#######################################################################
#####                                                             #####
#####                   plane_setup Function                      #####
#####                                                             #####
#######################################################################
#plane_setup function creates an instance of the Domain class named plane
#----------------------------------------------------------------------
#inputs:
#	-
#outputs:
#	-plane=data structure containing empty problem scenario
#----------------------------------------------------------------------	
	#Create problem scenario instance	
	plane=Domain(filename)
	
	return plane
	
	
def distributions_init(input_file,input_file_entry,n) :
#######################################################################
#####                                                             #####
#####               distributions_init Function                   #####
#####                                                             #####
#######################################################################
#distributions_init function determines whether a distribution should be
#initialized using a function or file.
#----------------------------------------------------------------------
#inputs:
#	-input_file=.json file containing plane setup information
#	-input_file_entry=.json entry containing distribution specification
#					  (either 'function' or 'filename')
#	-n=number of nodes
#outputs:
#	-init_var=array containing initialized distribution values
#	-init_var_type=either 'function' or 'file'
#----------------------------------------------------------------------
	
	#Check to see if distribution is specified from function or file
	if (isinstance(input_file_entry,str)):
		init_var=input_file_entry
		
		#If from function, initialize distribution to zero and set type 
		#to 'function'
		if ((init_var=='function') or (init_var=='even')):
			init_var_type=init_var
			init_var=np.zeros((n+1,), dtype=np.float64)

		#Anything else is interpreted as a file name containing the 
		#distribution
		else :
			filename=init_var
			init_var,init_var_type=file_dist(filename)
			
	#If neither function nor file is specified, set distribution constant
	#value and print notification
	else:
		init_var,init_var_type=constant_dists(input_file,input_file_entry,n)

	return init_var,init_var_type


def constant_dists(input_file,input_file_entry,n):
#######################################################################
#####                                                             #####
#####                   constant_dists Function                   #####
#####                                                             #####
#######################################################################	
#constant_dists sets all distributions specified by a constant in the
#input file to that constant
#----------------------------------------------------------------------
#inputs:
#	-input_file=input distribution file in ".json" format.
#	-input_file_entry=constant value specified in the input distribution
#	 file to which the distribution is to be set.
#	-n=number of nodes along wing
#outputs:
#	-init_var=distribution now initialized to correct value
#	-init_var_type=if file format is correct, return 'file'
#----------------------------------------------------------------------

	with open(input_file) as which_dist:
		dist_data=json.load(which_dist,object_pairs_hook=OrderedDict)
		
		if (input_file_entry==dist_data["wing"]["chord"]):
			init_var,init_var_type=set_constant_dist(input_file_entry,n)
			print('Chord distribution initialized to spanwise-constant value: 			'+str(input_file_entry))
			
		elif (input_file_entry==dist_data["wing"]["thickness_chord"]):
			init_var,init_var_type=set_constant_dist(input_file_entry,n)
			print('Thickness-to-chord ratio distribution initialized to spanwise-constant value:	'+str(input_file_entry))
			
		elif (input_file_entry==dist_data["weight"]["nonstructural_distribution"]):
			init_var,init_var_type=set_constant_dist(input_file_entry,n)	
			print('Non-structural weight initialized to spanwise-constant value:			'+str(input_file_entry))
			
		elif (input_file_entry==dist_data["spar"]["height"]):
			init_var,init_var_type=set_constant_dist(input_file_entry,n)		
			print('Spar height initialized to spanwise-constant value:				'+str(input_file_entry))	

	return init_var,init_var_type
	
	
def set_constant_dist(input_file_entry,n):
#######################################################################
#####                                                             #####
#####                 set_constant_dist Function                  #####
#####                                                             #####
#######################################################################	
#set_constant_dist function sets a distribution to a constant value
#----------------------------------------------------------------------
#inputs:
#	-input_file_entry=constant value specified in the input file to
#	 which the distribution is to be set.
#	-n=number of nodes along wing
#outputs:
#	-init_var=distribution now initialized to correct value
#	-init_var_type=if file format is correct, return 'file'
#----------------------------------------------------------------------

	init_var=np.zeros((n+1,), dtype=np.float64)
	for i in range (0,n+1):
		init_var[i]=input_file_entry
	init_var_type='constant'
	
	return init_var,init_var_type


def file_dist(filename):
#######################################################################
#####                                                             #####
#####                     file_dist Function                      #####
#####                                                             #####
#######################################################################	
#file_dist checks for correct distribution file entry format
#----------------------------------------------------------------------
#inputs:
#	-filename=name of distribution input file
#outputs:
#	-init_var_type=if file format is correct, return 'file'
#----------------------------------------------------------------------
	#Check for correct file format
	if (filename.endswith('.json')):
		with open(filename) as file:
			data = json.load(file, object_pairs_hook = OrderedDict)
			#check for correct key in .json file
			if ("distribution" in data):
				#Set init_var to distribution in file
				init_var=data["distribution"]
			else:
				sys.exit('File read error! No "distribution" key in "'+filename+'"')
				exit
		init_var_type='file'
	else :
		sys.exit('File read error! "'+filename+'" is invalid distribution filename. Please specify distribution file in .json format')
		exit

	return init_var,init_var_type

