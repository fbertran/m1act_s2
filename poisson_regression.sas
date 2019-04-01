title1 'Deaths from testicular cancer in Japan';
       data tcancer;
          input age pop1 c1 pop2 c2 pop3 c3 pop4 c4 pop5 c5;
          array p(5) pop1-pop5;       /* populations at each year group */
          array d(5) c1-c5;         /* cancer cases for each year group */
          label
             agegrp = 'age in 5yr interval'
             logpop = 'log-population';
            /* Produce a separate line for each year/age combination    */
           agegrp=age/5;                  /* age in five year intervals */
           do j=1 to 5;                   /* for each line read in . .  */
              deaths=d(j);                   /* number of cancer deaths */
              logpop=log(p(j));                    /* log of population */
              year=j-1;                  /* recode the years as 0,...,4 */
              yearc=year;                              /* year category */
              cohort=agegrp - year + 4; /* identify the diagonal cohort */
              output;     /* produce five output lines for each read in */
           end;
           drop j pop1-pop5 c1-c5;           /* omit unneeded variables */
           datalines;
 0 15501 17 26914 51 21027 65 20246 69 21596  74
 5 14236  . 25380  6 26613  7 20885  8 20051   7
10 13270  . 23492  3 25324  3 26540  7 20718  11
15 12658  2 21881  6 23211 15 24931 25 26182  39
20 10696  5 20402 27 21263 39 22228 56 24033  83
25  7563  5 17242 40 19994 58 20606 97 21805 125
30  7074  7 12609 18 17128 54 19864 77 20750 129
35  7038 10 11712 13 12476 36 17001 70 19890 101
40  6418  9 11478 26 11450 32 12275 29 16794  67
45  5981  7 10274 16 11157 26 11147 34 11962  37
50  4944  7  9325 16  9828 27 10705 27 10741  29
55  3994  7  7562 17  8718 19  9206 32 10086  39
60  3098  6  5902 13  6796 21  7869 21  8399  31
65  2317  4  4244 12  4911 26  5728 29  6715  34
70  1513  7  2845 17  3197 22  3737 25  4448  33
75   688  5  1587  9  1812 10  2061 25  2482  31
80   264  2   583  6   787  6   904 14  1068   9
85    73  2   179  2   246  3   335  3   419   3
 run;
       title2 'Model (5.1): Cohort, age, year; No population information';
       proc genmod;
          class yearc agegrp cohort;
          model deaths  = yearc agegrp cohort / type1 type3 dist=Poisson;
       run;
       proc genmod;
          class yearc agegrp cohort;
          model deaths  = agegrp cohort yearc/ type1 type3 dist=Poisson;
       run;
       title2 'Model (5.2): Cohort, age, year, and offset log(Pop)';
       proc genmod;
          class yearc agegrp cohort;
*not used with SAS Studio          
		  ods listing exclude obstats;  /* turn off the obstats listing */
          output out=fitted  reschi=reschi; /* create an output dataset */
          model deaths  = yearc agegrp cohort
              /  obstats type1 type3 offset=logpop dist=Poisson;
       run;
       proc sgplot data=fitted;              /* bubble plot of residuals */
          bubble y=year x=agegrp size=reschi;
       run;
       
       



title1 'Diversity of species';
      data species;
INFILE DATALINES DELIMITER=',' DSD;
         input name :$15. species area d2neigh d2sc adjsp;
		length name $20.;
		loga=log(area);
		area=area/100;
		adjsp=adjsp/100;
		d2sc=d2sc/100;
		d2neigh=d2neigh/100;
		isa=0; if name='Isabela' then isa=1;
		label 
			species = '# of species on island'
			area = 'in sq-km'
			loga = 'log area'
			d2neigh = 'distance to nearest neighbor'
			d2sc = 'distance to Santa Cruz'
			adjsp = '# of species on adjacent island';		
          cards;       
Baltra,58,25.09,.6,.6,44
Bartolome,31,1.24,.6,26.3,237
Caldwell,3,.21,2.8,58.7,5
Champion,25,.10,1.9,47.4,2
Coamafio,2,.05,1.9,1.9,444
Daphne Major,18,.34,8.0,8.0,44
Daphne Minor,24,.08,6.0,12.0,18
Darwin,10,2.33,34.1,290.2,21
Eden,8,.03,.4,.4,108
Enderby,2,.18,2.6,50.2,25
Espafiola,97,58.27,1.1,88.3,58
Fernandina,93,634.49,4.3,95.3,347
Gardner A,58,.57,1.1,93.1,97
Gardner B,5,.78,4.6,62.2,3
Genovesa,40,17.35,47.4,92.2,51
Isabela,347,4669.32,.7,28.1,91
Marchena,51,129.49,29.1,85.9,104
Onslow,2,.01,3.3,45.9,25
Pinta,104,59.56,29.1,119.6,51
Pinzon,108,17.95,10.7,10.7,8
Las Plazas,12,.23,.5,.6,58
Rabida,70,4.89,4.4,24.4,237
San Cristobal,280,551.62,45.2,66.6,58
San Salvador,237,572.33,.2,19.8,70
Santa Cruz,444,903.82,.6,.0,8
Santa Fe,62,24.08,16.5,16.5,444
Santa Maria,285,170.92,2.6,49.2,25
Seymour,44,1.84,.6,9.6,58
Tortuga,16,1.24,6.8,50.9,108
Wolf,21,2.85,34.1,254.7,10 
run;


proc print;
run;

title2 'Fit Poisson model with all pairwise interactions';
proc genmod data=species;
	ods output obstats=fit;
	model species = loga | d2neigh | adjsp | d2sc @2 isa 
	/ dist=Poisson obstats type1 type3;
run;

title2 'Fit QuasiPoisson model with all pairwise interactions';
proc genmod data=species;
	ods output obstats=fit;
	model species = loga | d2neigh | adjsp | d2sc @2 isa 
	/ dist=Poisson obstats type1 type3 scale=pearson;
run;

title2 'Fit Negative Binomial model with all pairwise interactions';
proc genmod data=species;
	ods output obstats=fit;
	model species = loga | d2neigh | adjsp | d2sc @2 isa 
	/ dist=Negbin obstats type1 type3;
run;



 


