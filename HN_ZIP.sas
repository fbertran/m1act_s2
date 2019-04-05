

************************************************
****  CODE DE LA MACRO ENVELOPPES SIMULEES  ****
****                   ZIP                  ****
************************************************
;

%macro halfnorm_zip(
	data = _last_,
	resp = ,
	coefg = ,
	coefb = ,
	g = ,
	b = ,
	gv = ,
	bv = ,
	out = _res_,
	seed = 0,
	nres = 19 );
	
options nonotes ;
data ZIPD;
set &data;
	ObsNumber+1;

%macro coefvar (Coef = , Vars = );
%local I cvarList stop ;
%if %length(&&&Vars)=0 %then %let cvarList=&&&coef;

%if %length(&&&Vars) ne 0 %then %do;
%let Stop = %eval( %length(&&&Coef)-%length(%sysfunc(compress(&&&Coef)))+1);
	%do I = 1 %to &stop;
		%if (&I = 1 ) %then
			%let cvarList = %scan(&&&Coef,&i,' ')+;
		%else %if (&I = &Stop) %then
			%let cvarList = &cvarList %scan(&&&Coef,&i,' ')*%scan(&&&Vars,&i-1,' ');
		%else %let cvarList = &cvarList %scan(&&&Coef,&i,' ')*%scan(&&&Vars,&i-1,' ')+;
	%end;
%end;
&cvarList
%mend coefvar;

%macro coefval(Coef = , Val = );
%local I cvalList stop;
%let Stop = %eval( %length(&&&Coef)-%length(%sysfunc(compress(&&&Coef)))+1);
	%do i = 1 %to &stop ;
		%let cvalList = &cvalList %scan(&&&Coef,&i,' ') = %scan(&&&Val,&i,' ') ;
	%end;
&cvalList
%mend coefval;

%macro nwords(text= );
	%eval(%sysfunc(count(%cmpres(&&&text) , %str( ) )) +1)
%mend nwords ;
%let abort=0;
%if ( %nwords( text = &coefg ) ^= %nwords( text = &gv ) ) %then %do;
	%put ERROR: The number of coefficients should equal to number of initial values;
	%left abort=1;
	%goto done;
	%end;
	
%if ( %nwords( text = &coefb ) ^= %nwords( text = &bv ) ) %then %do;
	%put ERROR: The number of coefficients should equal to number of initial values;
	%left abort=1;
	%goto done;
	%end;

%if ( %nwords( text = &coefg ) ^= %nwords( text = &g ) +1) %then %do;
	%put ERROR The number of coefficients should be one more than the number of variables;
	%left abort=1;
	%goto done;
	%end;

%if ( %nwords( text = &coefb ) ^= %nwords( text = &b ) +1) %then %do;
	%put ERROR: The number of coefficients should be one more than the number of variables;
	%left abort=1;
	%goto done;
	%end;

%if %length(&resp) = 0 %then %do;
	%put ERROR: No response has been specified.;
	%let abort=1;
	%goto done;
	%end;

	
proc nlmixed data=ZIPD;
	parms %coefval(Coef = coefg, val=gv ) %coefval(Coef = coefb, val=bv);
	gb = %coefvar(Coef = coefg , vars = g );
	p_0 = exp(gb)/(1 + exp(gb));
	bb = %coefvar(Coef = coefb , vars = b );
	p_1 = exp(bb);
		if &resp = 0 then prob = p_0 + (1-p_0)*exp(-p_1);
		if &resp ne 0 then prob= (1-p_0)*(((p_1**&resp)*exp(-p_1))/gamma(&resp+1));
		if &resp = 0 then ll = log(p_0 + (1-p_0)*exp(-p_1));
					 else ll = log((1-p_0)*(((p_1**&resp)*exp(-p_1))/gamma(&resp+1)));

model 	&resp ~ general(ll);
predict p_1 out  = lamb;
predict p_0 out  = psi;

proc sort data=lamb; by ObsNumber;

proc sort data=psi ; by ObsNumber;
run;

data final;
	merge lamb(rename=(Pred=p_lamb)) psi(rename=(pred=p_psi)); by obsnumber;
	mu    = (1-p_psi)*p_lamb;
	var   = mu+(p_psi/(1-p_psi))*mu**2;
	odres = abs((&resp-mu)/sqrt(var));

	
proc sort data=final out=final;
	by odres;
run;

data _obstat_(drop=p_lamb p_psi);
set final;
drop tValue Alpha DF Lower Probt StdErrPred Upper;
	array ybin{*} _x1 - _x&nres;
	array _y_{*} _z1 - _z&nres;
	array ynb{*} _w1 - _w&nres;
	drop i seed _x1 - _x&nres _z1 - _z&nres;
	
	retain seed &seed;
	do i=1 to dim(ynb);
					call ranbin(seed, 1, p_psi, ybin[i]);
					call ranpoi(seed,p_lamb,_y_[i]); /*Pas de paramËtre &trials. dans le cas de la poisson*/
					ynb[i]=0*ybin[i]+(1-ybin[i])*_y_[i];
	end;
run;
	
proc sort data=_obstat_;
	by ObsNumber;
run;

ods listing close;

%do i=1 %to &nres;

proc nlmixed data=_obstat_;
	parms %coefval(Coef = coefg, val=gv ) %coefval(Coef = coefb, val=bv);
	gb = %coefvar(Coef = coefg , vars = g );
	p_0 = exp(gb)/(1 + exp(gb));
	bb = %coefvar(Coef = coefb , vars = b );
	p_1 = exp(bb);
		if _w&i = 0 then prob = p_0 + (1-p_0)*exp(-p_1);
		if _w&i ne 0 then prob= (1-p_0)*(((p_1**_w&i)*exp(-p_1))/gamma(_w&i+1));
		if _w&i = 0 then ll = log(p_0 + (1-p_0)*exp(-p_1));
					 else ll = log((1-p_0)*(((p_1**_w&i)*exp(-p_1))/gamma(_w&i+1)));

model 	_w&i ~ general(ll);
predict p_1 out =lamb&i;
predict p_0 out =psi&i;

proc sort data=lamb&i;
	by ObsNumber;

proc sort data=psi&i;
	by ObsNumber;
run;

data final&i;
	merge lamb&i(rename=(Pred=p_lamb&i)) psi&i(rename=(pred=p_psi&i));
		by obsnumber;
	mu&i = (1-p_psi&i)*p_lamb&i;	
	var&i = mu&i*(p_psi&i/(1-p_psi&i))*mu&i**2;
	odres&i = abs((_w&i-mu&i)/sqrt(var&i));
run;

data final&i;
	set final&i;
	keep ObsNumber _w&i odres&i;
run;
	
proc sort data=final&i out=final&i;
	by ObsNumber;
run;

%end; /* end %do i */

data _obstat_;
	merge _obstat_
	%do i=1  %to &nres;
		final&i
		%end;
		;
run;

proc sort data=_obstat_;
	by odres;
run;


proc iml;
start sortcols(X);
	do i=1 to ncol(X);
		xi=x[,i];
		if any(xi=.) then do;
			mi=xi[loc(xi=.),];
			xi=xi[loc(xi^=.),];
			end;
		else free mi;
		t=xi;
		r=rank(xi);
		t[r]=xi;
		x[,i]=mi//t;
		end;
	finish;
start symput(macnm,scal);
	call execute('%let ',macnm,'=',char(scal),';');
finish;

	use _obstat_;
	read all var("odres1" : "odres&nres") into X;
	nc=0;
	do i=1 to ncol(X);
		if ^all(X[,i]=.) then do;
			Y=Y || X[,i];
			nc=nc+1;
			end;
		end;
	run symput('nc',nc);
	run sortcols(Y);
	create _sorted_ from Y;
	append from Y;
	quit;
	
data _obstat_;
	merge _obstat_ _sorted_;
	array _res_{*} col1-col&nc;
	drop odres1-odres&nres;
	
	resmin = min(of col1-col&nc);
	resmax = max(of col1-col&nc);
	resmean = mean(of col1-col&nc);
run;

data &out;
	set _obstat_ nobs=nobs end=eof;
	drop col1-col&nres;
	_z_ = probit((_n_ + nobs - .125)/(2*nobs+.5));
	label _z_='Expected value of half normal quantile';
	if oef then call symput('nobs',left(put(nobs,best8.)));
run;
ods listing;

/*
legend1 across=2 origin=(10,32) value=(tick=1 justify=center 'Observed'
			tick=2 justify=center 'Mean' tick=3 justify=center 'Minimum'
			tick=4 justify=center 'Maximum')
		position=(inside top)
		label=(justify=center 'Standardized Pearson Residual' position=top) mode=protect;
		
symbol1 i=join c=blue line=1 w=0.5 v=star height=0.5;
symbol2 i=join c=red line=3 w=2;
symbol3 i=join c=green line=5 w=2;
symbol4 i=join c=green line=7 w=2;

axis1 label=(a=90 "Standardized Pearson residuals") order=0 to 2 by 0.5;
axis2 label=('Half-normal scores') order=0 to 2.7 by 0.3;

proc greplay igout=gseg nofs;
	delete _all_;
quit;

options notes;
goptions ftext=swissb;
*/

ods listing gpath="/folders/myfolders";
proc sgplot data=&out ;/*gout=gseg*/
	series x=_z_ y=odres;
	series x=_z_ y=resmean;
	series x=_z_ y=resmin;
	series x=_z_ y=resmax;
/* vaxis=axis1 frame legend=legend1 haxis=axis2 name='ZIB';*/
run;


proc datasets nofs nolist nowarn library=work memtype=(data);
	delete _obstat_ ZIPD psi final lamb _sorted_
		%do i=1 %to &nres;
			final&i lamb&i psi&i
			%end;
		;
		run; quit;
%done:
%if &abort %then %put ERROR: The HALFNORM_ZIP macro ended abnormally.;
%mend halfnorm_zip;




















************************************************
****  TEST DE LA MACRO ENVELOPPES SIMULEES  ****
****                   ZIP                  ****
************************************************
;

data test;
	bp_0 = -0.5;
	bp_1 = 0.2;
	bll_0 = 0.3;
	bll_1 = 0.5;
	
	seed = 2329;
	
do x=0.12 to 3.2 by 0.02;
		call rannor(seed,m);
		m=abs(floor(10*m));
			if m>0 then do;
			eta_prob = bp_0 + bp_1*x;
			
p_0 = exp(eta_prob) / (1 + exp(eta_prob));
call ranbin(seed,1,p_0,zero_inflate);
		
eta_lambda = bll_0 + bll_1*x;
p_1 = exp(eta_lambda);
				call ranpoi(seed,p_1,ynb1);
				yzip=0*zero_inflate+(1-zero_inflate)*ynb1;
				output;
		end;
end;
keep yzip x;
run;

%halfnorm_zip(data= test, 
resp =yzip, 
coefg = bp_0 bp_1 , 
coefb = bll_0 bll_1, 
g=x, 
b=x, 
gv= 0 0, 
bv = 0 0, 
out = pp, 
seed = 2006, 
nres = 19);

