data a; /* generate the data */
call streaminit(1234);
do kk=1 to 10000;
x1 = rannor(1234);
x2 = rannor(1234);
x3 = rannor(1234);
theta = 1;
mu = exp(1 + .3*x1 + .3*x2);
parm1 = 1/(1+mu/theta);
yneg = rand('NEGB',parm1,theta);
ypoi = ranpoi(1234,mu);
pzeroorig = cdf('LOGISTIC',x3*2);
if ranuni(1234)>pzeroorig then do;
ynegzim = yneg;
ypoizim = ypoi;
end;
else do;
ynegzim = 0;
ypoizim = 0;
end;
y=ynegzim;
output ;
end ;
run;














*compare Poisson simulated values with fit from P, NB, ZIP and ZINB.;
proc freq data=a;
   table ypoi / out=obs;
run;


proc genmod data=a;
model ypoi=x1 x2/dist=poisson;
output out=poisson predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ypoi=x1 x2/dist=negbin;
output out=negbin predicted=pred pzero=pzero;
ods output ParameterEstimates=nbparms;
run;

proc genmod data=a;
model ypoi=x1 x2/dist=zip;
zeromodel x3;
output out=zip predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ypoi=x1 x2/dist=zinb;
zeromodel x3;
output out=zinb predicted=pred pzero=pzero;
ods output ParameterEstimates=zinbparms;
run;





proc means data=a noprint;
   var ypoi;
   output out=maxcount max=max N=N;
run;

data _null_;
  set maxcount;
  call symput('N',N);
  call symput('max',max);
run;

%let max=%sysfunc(strip(&max));



data poisson(drop= i);
   set poisson;
   lambda=pred;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pdf('POISSON',i,lambda);
      else        ep{i}= pdf('POISSON',i,lambda);
      c{i}=ifn(ypoi=i,1,0);
   end;
run;


proc means data=poisson noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=poisson) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data poissonprob;
   merge ep p;
   poissondiff=p-poisson;
   ypoi=_N_ -1;
   label poisson='Poisson Probabilities'
         p='Relative Frequencies'
         poissondiff='Observed minus Predicted';
run;







data zip(drop= i);
   set zip;
   lambda=pred/(1-pzero);
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('POISSON',i,lambda);
      else        ep{i}=         (1-pzero)*pdf('POISSON',i,lambda);
      c{i}=ifn(ypoi=i,1,0);
   end;
run;


proc means data=zip noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zip) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zipprob;
   merge ep p;
   zipdiff=p-zip;
   ypoi=_N_ -1;
   label zip='ZIPoisson Probabilities'
         p='Relative Frequencies'
         zipdiff='Observed minus Predicted';
run;



data nbparms;
   set nbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('knb',estimate);
run;

data negbin(drop= i);
   set negbin;
   lambda=pred;
   knb=&knb;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb));
      else        ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb)); 
      c{i}=ifn(ypoi=i,1,0);  
   end;
run;


proc means data=negbin noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=nb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data negbinprob;
   merge ep p;
   nbdiff=p-nb;
   ypoi=_N_ -1;
   label nb='NB Probabilities'
         p='Relative Frequencies'
         nbdiff='Observed minus Predicted';
run;










data zinbparms;
   set zinbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('k',estimate);
run;

data zinb(drop= i);
   set zinb;
   lambda=pred/(1-pzero);
   k=&k;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k));
      else        ep{i}=         (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k)); 
      c{i}=ifn(ypoi=i,1,0);  
   end;
run;


proc means data=zinb noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zinb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zinbprob;
   merge ep p;
   zinbdiff=p-zinb;
   ypoi=_N_ -1;
   label zinb='ZINB Probabilities'
         p='Relative Frequencies'
         zinbdiff='Observed minus Predicted';
run;


data compare;
   merge poissonprob zipprob negbinprob zinbprob;
   by ypoi;
run;


proc sgplot;
   yaxis label='Probability';
   xaxis label='Simulated ypoi Value';
   series y=p  x=ypoi / name='obs' legendlabel='Observed'
      lineattrs=(color=black thickness=4px);
   series y=poisson x=ypoi / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nb x=ypoi/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zip x=ypoi/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinb x=ypoi/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;



proc sgplot;
   yaxis label='Probability differences';
   xaxis label='Simulated ypoi Value';
   series y=poissondiff x=ypoi / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nbdiff x=ypoi/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zipdiff x=ypoi/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinbdiff x=ypoi/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;
















*compare ZIPoisson simulated values with fit from P, NB, ZIP and ZINB.;
proc freq data=a;
   table ypoizim / out=obs;
run;


proc genmod data=a;
model ypoizim=x1 x2/dist=poisson;
output out=poisson predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ypoizim=x1 x2/dist=negbin;
output out=negbin predicted=pred pzero=pzero;
ods output ParameterEstimates=nbparms;
run;

proc genmod data=a;
model ypoizim=x1 x2/dist=zip;
zeromodel x3;
output out=zip predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ypoizim=x1 x2/dist=zinb;
zeromodel x3;
output out=zinb predicted=pred pzero=pzero;
ods output ParameterEstimates=zinbparms;
run;





proc means data=a noprint;
   var ypoizim;
   output out=maxcount max=max N=N;
run;

data _null_;
  set maxcount;
  call symput('N',N);
  call symput('max',max);
run;

%let max=%sysfunc(strip(&max));



data poisson(drop= i);
   set poisson;
   lambda=pred;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pdf('POISSON',i,lambda);
      else        ep{i}= pdf('POISSON',i,lambda);
      c{i}=ifn(ypoizim=i,1,0);
   end;
run;


proc means data=poisson noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=poisson) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data poissonprob;
   merge ep p;
   poissondiff=p-poisson;
   ypoizim=_N_ -1;
   label poisson='ZIPoisson Probabilities'
         p='Relative Frequencies'
         poissondiff='Observed minus Predicted';
run;







data zip(drop= i);
   set zip;
   lambda=pred/(1-pzero);
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('POISSON',i,lambda);
      else        ep{i}=         (1-pzero)*pdf('POISSON',i,lambda);
      c{i}=ifn(ypoizim=i,1,0);
   end;
run;


proc means data=zip noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zip) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zipprob;
   merge ep p;
   zipdiff=p-zip;
   ypoizim=_N_ -1;
   label zip='Poisson Probabilities'
         p='Relative Frequencies'
         zipdiff='Observed minus Predicted';
run;



data nbparms;
   set nbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('knb',estimate);
run;

data negbin(drop= i);
   set negbin;
   lambda=pred;
   knb=&knb;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb));
      else        ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb)); 
      c{i}=ifn(ypoizim=i,1,0);  
   end;
run;


proc means data=negbin noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=nb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data negbinprob;
   merge ep p;
   nbdiff=p-nb;
   ypoizim=_N_ -1;
   label nb='NB Probabilities'
         p='Relative Frequencies'
         nbdiff='Observed minus Predicted';
run;










data zinbparms;
   set zinbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('k',estimate);
run;

data zinb(drop= i);
   set zinb;
   lambda=pred/(1-pzero);
   k=&k;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k));
      else        ep{i}=         (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k)); 
      c{i}=ifn(ypoizim=i,1,0);  
   end;
run;


proc means data=zinb noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zinb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zinbprob;
   merge ep p;
   zinbdiff=p-zinb;
   ypoizim=_N_ -1;
   label zinb='ZINB Probabilities'
         p='Relative Frequencies'
         zinbdiff='Observed minus Predicted';
run;


data compare;
   merge poissonprob zipprob negbinprob zinbprob;
   by ypoizim;
run;


proc sgplot;
   yaxis label='Probability';
   xaxis label='Simulated ypoizim Value';
   series y=p  x=ypoizim / name='obs' legendlabel='Observed'
      lineattrs=(color=black thickness=4px);
   series y=poisson x=ypoizim / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nb x=ypoizim/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zip x=ypoizim/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinb x=ypoizim/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;



proc sgplot;
   yaxis label='Probability differences';
   xaxis label='Simulated ypoizim Value';
   series y=poissondiff x=ypoizim / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nbdiff x=ypoizim/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zipdiff x=ypoizim/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinbdiff x=ypoizim/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;






























*compare NegBin simulated values with fit from P, NB, ZIP and ZINB.;
proc freq data=a;
   table yneg / out=obs;
run;


proc genmod data=a;
model yneg=x1 x2/dist=poisson;
output out=poisson predicted=pred pzero=pzero;
run;

proc genmod data=a;
model yneg=x1 x2/dist=negbin;
output out=negbin predicted=pred pzero=pzero;
ods output ParameterEstimates=nbparms;
run;

proc genmod data=a;
model yneg=x1 x2/dist=zip;
zeromodel x3;
output out=zip predicted=pred pzero=pzero;
run;

proc genmod data=a;
model yneg=x1 x2/dist=zinb;
zeromodel x3;
output out=zinb predicted=pred pzero=pzero;
ods output ParameterEstimates=zinbparms;
run;





proc means data=a noprint;
   var yneg;
   output out=maxcount max=max N=N;
run;

data _null_;
  set maxcount;
  call symput('N',N);
  call symput('max',max);
run;

%let max=%sysfunc(strip(&max));



data poisson(drop= i);
   set poisson;
   lambda=pred;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pdf('POISSON',i,lambda);
      else        ep{i}= pdf('POISSON',i,lambda);
      c{i}=ifn(yneg=i,1,0);
   end;
run;


proc means data=poisson noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=poisson) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data poissonprob;
   merge ep p;
   poissondiff=p-poisson;
   yneg=_N_ -1;
   label poisson='Poisson Probabilities'
         p='Relative Frequencies'
         poissondiff='Observed minus Predicted';
run;







data zip(drop= i);
   set zip;
   lambda=pred/(1-pzero);
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('POISSON',i,lambda);
      else        ep{i}=         (1-pzero)*pdf('POISSON',i,lambda);
      c{i}=ifn(yneg=i,1,0);
   end;
run;


proc means data=zip noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zip) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zipprob;
   merge ep p;
   zipdiff=p-zip;
   yneg=_N_ -1;
   label zip='ZIPoisson Probabilities'
         p='Relative Frequencies'
         zipdiff='Observed minus Predicted';
run;






data nbparms;
   set nbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('knb',estimate);
run;

data negbin(drop= i);
   set negbin;
   lambda=pred;
   knb=&knb;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb));
      else        ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb)); 
      c{i}=ifn(yneg=i,1,0);  
   end;
run;


proc means data=negbin noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=nb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data negbinprob;
   merge ep p;
   nbdiff=p-nb;
   yneg=_N_ -1;
   label nb='NB Probabilities'
         p='Relative Frequencies'
         nbdiff='Observed minus Predicted';
run;










data zinbparms;
   set zinbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('k',estimate);
run;

data zinb(drop= i);
   set zinb;
   lambda=pred/(1-pzero);
   k=&k;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k));
      else        ep{i}=         (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k)); 
      c{i}=ifn(yneg=i,1,0);  
   end;
run;


proc means data=zinb noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zinb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zinbprob;
   merge ep p;
   zinbdiff=p-zinb;
   yneg=_N_ -1;
   label zinb='ZINB Probabilities'
         p='Relative Frequencies'
         zinbdiff='Observed minus Predicted';
run;


data compare;
   merge poissonprob zipprob negbinprob zinbprob;
   by yneg;
run;


proc sgplot;
   yaxis label='Probability';
   xaxis label='Simulated yneg Value';
   series y=p  x=yneg / name='obs' legendlabel='Observed'
      lineattrs=(color=black thickness=4px);
   series y=poisson x=yneg / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nb x=yneg/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zip x=yneg/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinb x=yneg/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;



proc sgplot;
   yaxis label='Probability differences';
   xaxis label='Simulated yneg Value';
   series y=poissondiff x=yneg / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nbdiff x=yneg/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zipdiff x=yneg/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinbdiff x=yneg/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;






























*compare NegBin simulated values with fit from P, NB, ZIP and ZINB.;
proc freq data=a;
   table ynegzim / out=obs;
run;


proc genmod data=a;
model ynegzim=x1 x2/dist=poisson;
output out=poisson predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ynegzim=x1 x2/dist=negbin;
output out=negbin predicted=pred pzero=pzero;
ods output ParameterEstimates=nbparms;
run;

proc genmod data=a;
model ynegzim=x1 x2/dist=zip;
zeromodel x3;
output out=zip predicted=pred pzero=pzero;
run;

proc genmod data=a;
model ynegzim=x1 x2/dist=zinb;
zeromodel x3;
output out=zinb predicted=pred pzero=pzero;
ods output ParameterEstimates=zinbparms;
run;





proc means data=a noprint;
   var ynegzim;
   output out=maxcount max=max N=N;
run;

data _null_;
  set maxcount;
  call symput('N',N);
  call symput('max',max);
run;

%let max=%sysfunc(strip(&max));



data poisson(drop= i);
   set poisson;
   lambda=pred;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pdf('POISSON',i,lambda);
      else        ep{i}= pdf('POISSON',i,lambda);
      c{i}=ifn(ynegzim=i,1,0);
   end;
run;


proc means data=poisson noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=poisson) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data poissonprob;
   merge ep p;
   poissondiff=p-poisson;
   ynegzim=_N_ -1;
   label poisson='Poisson Probabilities'
         p='Relative Frequencies'
         poissondiff='Observed minus Predicted';
run;







data zip(drop= i);
   set zip;
   lambda=pred/(1-pzero);
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('POISSON',i,lambda);
      else        ep{i}=         (1-pzero)*pdf('POISSON',i,lambda);
      c{i}=ifn(ynegzim=i,1,0);
   end;
run;


proc means data=zip noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zip) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zipprob;
   merge ep p;
   zipdiff=p-zip;
   ynegzim=_N_ -1;
   label zip='ZIPoisson Probabilities'
         p='Relative Frequencies'
         zipdiff='Observed minus Predicted';
run;






data nbparms;
   set nbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('knb',estimate);
run;

data negbin(drop= i);
   set negbin;
   lambda=pred;
   knb=&knb;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb));
      else        ep{i}=         pdf('NEGBINOMIAL',i,(1/(1+knb*lambda)),(1/knb)); 
      c{i}=ifn(ynegzim=i,1,0);  
   end;
run;


proc means data=negbin noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=nb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data negbinprob;
   merge ep p;
   nbdiff=p-nb;
   ynegzim=_N_ -1;
   label nb='NB Probabilities'
         p='Relative Frequencies'
         nbdiff='Observed minus Predicted';
run;










data zinbparms;
   set zinbparms(where=(Parameter="Dispersion"));
   keep estimate;
   call symput('k',estimate);
run;

data zinb(drop= i);
   set zinb;
   lambda=pred/(1-pzero);
   k=&k;
   array ep{0:&max} ep0-ep&max;
   array c{0:&max} c0-c&max;
   do i = 0 to &max;
      if i=0 then ep{i}= pzero + (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k));
      else        ep{i}=         (1-pzero)*pdf('NEGBINOMIAL',i,(1/(1+k*lambda)),(1/k)); 
      c{i}=ifn(ynegzim=i,1,0);  
   end;
run;


proc means data=zinb noprint;
   var ep0 - ep&max c0-c&max;
   output out=ep(drop=_TYPE_ _FREQ_) mean(ep0-ep&max)=ep0-ep&max;
   output out=p(drop=_TYPE_ _FREQ_) mean(c0-c&max)=p0-p&max;
run;


proc transpose data=ep out=ep(rename=(col1=zinb) drop=_NAME_);
run;

proc transpose data=p out=p(rename=(col1=p) drop=_NAME_);
run;


data zinbprob;
   merge ep p;
   zinbdiff=p-zinb;
   ynegzim=_N_ -1;
   label zinb='ZINB Probabilities'
         p='Relative Frequencies'
         zinbdiff='Observed minus Predicted';
run;


data compare;
   merge poissonprob zipprob negbinprob zinbprob;
   by ynegzim;
run;


proc sgplot;
   yaxis label='Probability';
   xaxis label='Simulated ynegzim Value';
   series y=p  x=ynegzim / name='obs' legendlabel='Observed'
      lineattrs=(color=black thickness=4px);
   series y=poisson x=ynegzim / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nb x=ynegzim/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zip x=ynegzim/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinb x=ynegzim/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;



proc sgplot;
   yaxis label='Probability differences';
   xaxis label='Simulated ynegzim Value';
   series y=poissondiff x=ynegzim / name='poi' legendlabel='Poisson'
      lineattrs=(color=blue);
   series y=nbdiff x=ynegzim/ name='nb' legendlabel='Negative Binomial'
      lineattrs=(color=red);
   series y=zipdiff x=ynegzim/ name='zip' legendlabel='ZIP'
      lineattrs=(color=blue pattern=2);
   series y=zinbdiff x=ynegzim/ name='zinb' legendlabel='ZINB'
      lineattrs=(color=red pattern=2);
   discretelegend 'poi' 'zip' 'nb' 'zinb' 'obs' / title='Models:'
      location=inside position=ne across=2 down=3;
run;


