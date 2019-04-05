
************************************************
****  CODE DE LA MACRO ENVELOPPES SIMULEES  ****
****                  ZINB                  ****
************************************************
;

%macro halfnorm_zinb(
	data = _last_,
	resp = ,
	coefg = ,
	coefb = ,
	g = ,
	b = ,
	gv = ,
	bv = ,
	/*trials = , */
	out = _res_,
	seed = 0,
	nres = 19 );
	
options nonotes ;
data ZINBD;
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



proc nlmixed data=ZINBD;
	parms %coefval(Coef = coefg, val=gv ) %coefval(Coef = coefb, val=bv)  k=0.8 ; /* On ajoute ‡ la liste le coefficient de dispersion en attribuant une valeur intiale */
	gb = %coefvar(Coef = coefg , vars = g );
	p_0 = exp(gb)/(1 + exp(gb));
	bb = %coefvar(Coef = coefb , vars = b );
    lmdi=exp(bb);

    if &resp = 0 then prob = p_0 + (1-p_0)*((1+k*lmdi)**(1/k));
    if &resp ne 0 then prob=(1-p_0)*gamma(&resp+(1/k))*((k*lmdi)**&resp)/(gamma(&resp+1)*gamma(1/k)*((1+k*lmdi)**(&resp+(1/k)))); 
    if &resp = 0 then  ll = log(p_0 + (1-p_0)*((1+k*lmdi)**(1/k))) ;
	if &resp ne 0 then ll = log((1-p_0)) + lgamma(&resp+(1/k)) + &resp*log(k*lmdi) - lgamma(&resp+1)  - lgamma(1/k) - (&resp+(1/k))*log(1+k*lmdi);

    model 	&resp ~ general(ll);

predict lmdi out  = lamb;
predict p_0 out  = psi;
predict k out  = kdata;

proc sort data=lamb; by ObsNumber;

proc sort data=psi ; by ObsNumber;

proc sort data=kdata ; by ObsNumber;
run;

data final;
	merge lamb(rename=(Pred=pred_lmdi)) psi(rename=(pred=p_psi))kdata(rename=(pred=k)); by obsnumber; /* pred_lmdi remplace p_lamd*/
	mu    = (1-p_psi)*pred_lmdi;
	var   = mu*(1+pred_lmdi*(p_psi+k));
	odres = abs((&resp-mu)/sqrt(var));
	
proc sort data=final out=final;
	by odres;
run;


data _obstat_(drop=pred_lmdi p_psi k);
set final;
drop tValue Alpha DF Lower Probt StdErrPred Upper;
	array ynbin{*} _x1 - _x&nres; 
	array _y_{*} _z1 - _z&nres;
	array ynb{*} _w1 - _w&nres;
	drop i seed _x1 - _x&nres _z1 - _z&nres;
	
	retain seed &seed;
	do i=1 to dim(ynb);
		call ranbin(seed, 1, p_psi, ynbin[i]); 
		call streaminit(seed); 
		_y_[i]=RAND('POISSON', pred_lmdi*k*RAND('GAMMA', 1/k));
		ynb[i]=0*ynbin[i]+(1-ynbin[i])*_y_[i];
	end;
run;
	
proc sort data=_obstat_;
	by ObsNumber;
run;

ods listing close;

%do i=1 %to &nres;

proc nlmixed data=_obstat_  ; 
	parms %coefval(Coef = coefg, val=gv ) %coefval(Coef = coefb, val=bv) k=0.8 ;
	gb = %coefvar(Coef = coefg , vars = g );
	p_0 = exp(gb)/(1 + exp(gb));
	bb = %coefvar(Coef = coefb , vars = b );
    lmdi=exp(bb);
    
	if _w&i = 0 then prob = p_0 + (1-p_0)*((1+k*lmdi)**(1/k)); 
	if _w&i ne 0 then prob=(1-p_0)*gamma(_w&i+(1/k))*((k*lmdi)**_w&i)/(gamma(_w&i+1)*gamma(1/k)*((1+k*lmdi)**(_w&i+(1/k))));
	if _w&i = 0 then ll = log(p_0 + (1-p_0)*((1+k*lmdi)**(1/k)));
	if _w&i ne 0 then ll = log((1-p_0)) + lgamma(_w&i+(1/k)) + _w&i*log(k*lmdi) - lgamma(_w&i+1)  - lgamma(1/k) - (_w&i+(1/k))*log(1+k*lmdi); 

model 	_w&i ~ general(ll);

predict lmdi out =lamb&i;
predict p_0 out =psi&i;
predict k out  = kdata&i; 
%put "Attention: La table lamb&i, psi&i et  kdata&i ne seront pas créees si l'estimation correspondant à _w&i n'est pas convergente" ; 

proc sort data=lamb&i;
	by ObsNumber;

proc sort data=psi&i;
	by ObsNumber;

proc sort data=kdata&i;
	by ObsNumber;
run;

data final&i;
	merge lamb&i(rename=(Pred=pred_lmdi&i)) psi&i(rename=(pred=p_psi&i)) kdata&i(rename=(pred=k&i)); /* pred_lmdi remplace p_lamd*/
		by obsnumber;
	mu&i   = (1-p_psi&i)*pred_lmdi&i;
	var&i   = mu&i*(1+pred_lmdi&i*(p_psi&i+k&i));
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



ods listing gpath="/folders/myfolders";
proc sgplot data=&out ;/*gout=gseg*/
	series x=_z_ y=odres;
	series x=_z_ y=resmean;
	series x=_z_ y=resmin;
	series x=_z_ y=resmax;
run;


proc datasets nofs nolist nowarn library=work memtype=(data);
	delete _obstat_ ZINBD psi final lamb _sorted_
		%do i=1 %to &nres;
			final&i lamb&i psi&i
			%end;
		;
		run; quit;

%done:
%if &abort %then %put ERROR: The HALFNORM_ZINB macro ended abnormally.;
%mend halfnorm_zinb;

















************************************************
****  TEST DE LA MACRO ENVELOPPES SIMULEES  ****
****                  ZINB                  ****
************************************************
;

data test;
	bp_0 = -0.5;
	bp_1 = 0.2;
	bll_0 = 0.3;
	bll_1 = 0.5;
	k=0.5 ;

	seed = 2329;
	
do x=0.12 to 3.2 by 0.02;
		call rannor(seed,m);
		m=abs(floor(10*m));
			if m>0 then do;
			eta_prob = bp_0 + bp_1*x;
			
p_0 = exp(eta_prob) / (1 + exp(eta_prob)); 
call ranbin(seed,1,p_0,zero_inflate); 
		
eta_lambda = bll_0 + bll_1*x;
lmdi=exp(eta_lambda);
call streaminit(seed); /* fixation du seed */
ynb1=RAND('POISSON', lmdi*k*RAND('GAMMA', 1/k));
yzinb=0*zero_inflate+(1-zero_inflate)*ynb1;
		output;
		end;
end;
keep m yzinb x ;
run;

%halfnorm_zinb( data= test, 
resp =yzinb, 
coefg = bp_0 bp_1 , 
coefb = bll_0 bll_1, 
g=x, 
b=x, 
gv= 0 0 , 
bv = 0 0, 
/*trials = m, */
out = pp, 
seed = 2006, 
nres = 19);









