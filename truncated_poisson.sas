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
 model count = sex year year*sex / obstats type1 type3;
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
