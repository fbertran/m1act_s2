%macro tpr;
 y=_RESP_;
    lambda=exp(_XBETA_);
 eml=exp(-lambda);

 /* Truncated Poisson deviance */
 dev=0;
 if(y>1) then do;
  %inv(y);
  dev=-lam+y*log(lam)-log(1-exp(-lam));
 end;
 dev=dev+lambda-y*log(lambda)+log(1-eml);
 deviance d = 2*dev;

/* GENMOD links */
 %inv(_MEAN_);
 fwdlink ey = log(lam);
 invlink link = lambda/(1-eml);
 variance var = (lambda*(lambda+1)*(1-eml)-lambda**2)/(1-eml)**2;
%mend tpr;

%macro inv(ev);
 if &ev LE 1 then lam=. ;
 else do;
  lamlo=&ev-1;
  lamhi=&ev;
  do until (abs(lamhi-lamlo)<1e-7);
   lam=(lamhi+lamlo)/2;
   mal=lam/(1-exp(-lam));
   if mal GE &ev then lamhi=lam;
   if mal LE &ev then lamlo=lam;
  end;
 end;
%mend inv;

data skunk;
INFILE DATALINES DELIMITER=' ';
 input count freq year sex $ @@;
 yr =year; sx=sex; /* make copies of YEAR and SEX */
 label
  count='count'
  freq='frequency'
  year='1977 or 1978'
  sex='F M';
 datalines;
1 1 77 F 2 2 77 F 3 4 77 F 4 2 77 F 5 1 77 F
1 3 77 M 2 0 77 M 3 3 77 M 4 3 77 M 5 2 77 M
6 2 77 M 1 7 78 F 2 7 78 F 3 3 78 F 4 1 78 F
5 2 78 F 1 4 78 M 2 3 78 M 3 1 78 M
run;


proc genmod data=skunk;
 %tpr;
 frequency freq;
 class sex year;
 model count = sex year year*sex / obstats type1 type3 scale=pearson;
* ods html exclude obstats;
 ods output obstats=fitted;
run;


data fit2;
 set fitted;
 drop sex year;
run;

data both;
 merge skunk fit2;
 lambda = exp(xbeta);
 p0 = exp(-lambda);
 lamup=exp(xbeta+1.96*std);
 lamlow=exp(xbeta-1.96*std);
 p0up=exp(-lamlow);
    p0low=exp(-lamup);
 se=std*lambda;
 sep0=std*lambda*p0;
run;

data both;
 set both;
 by yr sx;
 if not (first.yr|first.sx) then delete;
run;

proc print noobs;
 var yr sx lambda se lamlow lamup p0 sep0 p0low p0up;
run;












data rearrest;
INFILE DATALINES DELIMITER=' ';
 input count freq Expelled $ Type $ @@;
 ex =Expelled; ty=Type; /* make copies of YEAR and SEX */
 label
  count='count'
  freq='frequency'
  Expelled='Not Eff Other'
  Type='Obs Exp';
 datalines;
1 2755 Not Obs 1 2693 Not Exp 1 2431 Eff Obs 1 2423 Eff Exp 1 654 Other Obs
1 652 Other Exp 2 251 Not Obs 2 349 Not Exp 2 52 Eff Obs 2 65 Eff Exp
2 57 Other Obs 2 61 Other Exp 3 46 Not Obs 3 30 Not Exp 3 4 Eff Obs 3 1 Eff Exp
3 6 Other Obs 3 4 Other Exp 4 18 Not Obs 4 2 Not Exp 4 1 Eff Obs 5 2 Not Obs
5 1 Eff Obs 6 2 Not Obs
run;

proc genmod data=rearrest;
 %tpr;
 frequency freq;
 class Expelled;
 model count = Expelled / obstats type1 type3 noscale;
 ods output obstats=fitted;
run;


proc genmod data=rearrest;
 %tpr;
 frequency freq;
 class Expelled;
 model count = Expelled / obstats type1 type3 scale=pearson;
 ods output obstats=fitted;
run;

data fit2;
 set fitted;
 drop Expelled;
run;

data both;
 merge rearrest fit2;
 lambda = exp(xbeta);
 p0 = exp(-lambda);
 lamup=exp(xbeta+1.96*std);
 lamlow=exp(xbeta-1.96*std);
 p0up=exp(-lamlow);
    p0low=exp(-lamup);
 se=std*lambda;
 sep0=std*lambda*p0;
run;

proc sort data=rearrest;
	by ex;
run;


proc sort data=both;
	by ex;
run;

data both;
 set both;
 by ex;
 if not first.ex then delete;
run;

proc print noobs;
 var ex lambda se lamlow lamup p0 sep0 p0low p0up;
run;


