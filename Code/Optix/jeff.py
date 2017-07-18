import sys
import json
from collections import OrderedDict
sys.path.append('/home/aerolab/Dropbox/Transfer/Lift_Distribution/Wing_Structure/Code')
import wing_structure_m as ws
sys.path.append('.')

def script(args):
	design_vars = args[0]
	case_id = args[1]
	
	W=design_vars[0]
	#W_s=design_vars[0]
	b=abs(design_vars[1])
	

	# set up input files using design_vars
	with open('Phillips.json','r') as input_var:
		input_var=json.load(input_var,object_pairs_hook = OrderedDict)
	ng=input_var['limits']['hard_landing']
	nm=input_var['limits']['maneuvering']
	#input_var['weight']['structural_weight']=W_s
	input_var['weight']['total_weight']=W
	input_var['wing']['wing_span']=b
	input_var['weight']['root_weight']=(ng-1.0)/(nm+ng)*(W)
	input_var['wing']['wing_area']=(W)/15.0#((3000.0+W_s)/15.0)
	input_file='Phillips{}.json'.format(case_id)
	with open(input_file,'w') as file1:
		json.dump(input_var,file1,indent=4,sort_keys=False)

    # Execute your function
	#case=ws.plane_setup(input_file)
	case=ws.Domain(input_file,0)
	case.solver(1e-9)
	
	print('NS_weight:	',case.W.n)
	print('weight:	',case.W.tot)
	print('R_A:		',case.W.r/case.W.tot)
	print('wingspan:	',b)
	print('W/S		:',case.W.tot/case.wing.S)
	print('S		:',case.wing.S)
	print('R_n		:',case.W.r/case.W.tot)
    # Extract drag from results
	drag=case.D_i
	
	#~ target_Wn=3000.0
	#~ predicted_Wn=case.W.n+case.W.r

	#~ factor=0.01
	
	#~ penalty=((target_Wn-predicted_Wn)**2)*factor
	#~ print('penalty:	', penalty)
	#~ drag+=penalty

	return drag
