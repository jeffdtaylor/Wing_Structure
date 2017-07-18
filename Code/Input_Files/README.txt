Input files must be in .json format. There are many options for the input file.
Required entries are followed by [Req]. Items still under development are followed by [Dev].
The options are given below:

wing[Req] = variable containing wing data
	
	wing_area = variable containing wing area. If wing loading is specified, wing area is not required. If wing loading is not specified, wing_area is used to set the wing loading.
	wing_span[Req] = variable contianing wing span
	chord[Req] = variable containing chord length data
	loading = variable containing the design wing loading. If not specified, wing loading will be determined from total_weight and wing_area. If total_weight is not specified, loading must be given.
	
		function = chord length is specified by a function. Can have the following values:
			rectangular = chord length is defined such that the wing is rectangular with area specified by wing_area
			elliptic = chord is distributed using the elliptic distribution with wing area specified by wing_area
			custom = chord length is specified by a custom function defined in the chord("your_plane") function of the user_functions module
			
			or the following variable:
			taper = chord length is defined by the taper ratio given and scaled such that the wing has area specified by wing_area
	
		file[Dev] = chord is given in the file specified. File name must be in .json format
		constant = chord is set to the constant value specified.

	thickness_chord[Req] = variable containing t/c ratio data for the wing airfoil
		
		function = t/c is specified by a function Can have the following value:
			custom = t/c is specified by a custom function defined in the t_c("your_plane") function of the user_functions module

			or the following variable:
			root/tip = t/c ratio is linearly interpolated between the given root and tip values
				root[Req] = root t/c ratio
				tip[Req] = tip t/c ratio
			
		file[Dev] = t/c is given in the file specified. File name must be in .json format
		constant = t/c is set to the constant value specified.

	grid[Req] = number of nodes along the wing


spar[Req] = variable containing the spar data

	C_sigma = shape factor of the spar beam cross-section. If C_sigma is specified, beam_type and height need not be specified
	beam_type = defines the beam cross-section shape. Can have the following value:
		rectangular

		or the following variables:
		I = I beam
			flange_width[Req] = relative width of I beam flange
			flange_height[Req] = relative height of I beam flange
			web_width[Req] = relative width of I beam web
		box = box beam
			inner_height[Req] = relative inner height of box beam
			inner_width[Req] = relative inner width of box beam
			outer_width[Req] = relative outer width of box beam
		(note: the values of each of the beam dimensions above are relative values, meaning that only their values relative to one another matter)

	height = variable containing spar beam height distribution data.

		function = spar height is specified by a function. Can have the following values:
			fill = spar fills the max airfoil thickness at each section
			min_fill = spar height is set to the minimum airfoil section thickness at all sections
			custom = height is specified by a custom function defined in the height("spar","wing") function of the user_functions module
		
		file[Dev] = spar height is given in the file specified. File name must be in .json format.
		constant = spar height is set to the constant value specified.
	
	max_stress[Req] = maximum allowable stress in the beam. (Dependent on beam material)
	specific_weight[Req] = specific weight of beam material. (Dependent on beam material)


limits[Req] = design load factors for different flight phases
	
	maneuvering[Req] = maneuvering flight load limit
	hard_landing[Req] = hard-landing load limit


weight[Req] = variable containing weight data
	
	total_weight = variable that gives total aircraft weight. If total_weight is included, 
	root_weight = the weight carried at the wing root.
	nonstructural_weight = total nonstructural weight carried by the aircraft (includes both root weight and nonstructural weight distributed along wing)
	root_total = the ratio of root weight to total weight. If specified, root_weight is disregarded
	nonstructural_distribution[req] = variable containing the nonstructural weight distribution data.

		function = nonstructural weight distribution is specified by a function. Can have the following values:
			even = evenly distributes nonstructural weight across span
			custom = nonstructural weight is specified by a custom function defined in the Wns("your_plane") function of the user_functions module

		file[Dev] = nonstructural weight is given in the file specified. File name must be in .json format. If a file is used, the total_weight and nonstructural_weight will be disregarded.
		constant = nonstructural weight is set to the constant value specified at each section. If constant is used, the total_weight and nonstructural_weight will be disregarded.


flight[Req] = variable containing flight conditions
	
	density = air density
	velocity = freestream velocity


lift_distribution[Req] = variable containing lift distribution information

	B = array containing fourier coefficients beginning at n=2.

	
	
 
