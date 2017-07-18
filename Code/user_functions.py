import numpy as np
import math as ma

#Equation 30
def Wns(case):
	for i in range (0,case.wing.m+1) :
	#Eq. (23)
		if (case.W.n_type=='constant'):
			case.W.ntilde[i]=(case.W.n+case.W.s-case.W.r)*(case.L.ratio[i])-case.W.stilde[i]
		else:
			case.W.ntilde[i]=(case.W.tot-case.W.r)*(case.L.ratio[i])-case.W.stilde[i]
