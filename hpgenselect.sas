data;
	input count birth death medic;
	label
		birth = 'names on birth certificate'
		death = 'names on death certificate'
		medic = 'names on medical record';
	datalines;
	60 0 0 1
	49 0 1 0
	 4 0 1 1
   247 1 0 0
   112 1 0 1
   142 1 1 0
    12 1 1 1
run;
proc hpgenselect;
   class medic birth death;
   model count = medic|birth|death@2 / Distribution=Poisson;
   selection method=stepwise details=all HIERARCHY=singleclass;
run;

