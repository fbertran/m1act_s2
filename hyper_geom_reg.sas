data Data1;
	do row=1 to 2;
		do col=1 to 2;
			input count@@;
			output;
		end;
	end;
	datalines;
	4 5 12 74
run;

proc freq data=Data1;
	tables row * col/all expected;
	weight count;
	exact or chisq;
run;


%macro hyp(m,n,nn,x);
	lgamma(&m+1)-lgamma(&x+1)-lgamma(&m-&x+1)+lgamma(&nn-&m+1)-lgamma(&n-&x+1)-lgamma(&nn-&n-&m+&x+1)
%mend hyp;


%macro hypmean(m,n,nn,lor);
	den=0;
	hypmean=0;
	hypvar=0;
	do j=max(0,&m+&n-&nn) to min(&n,&m);
		dterm=exp(%hyp(&m,&n,&nn,j)+j*&lor);
		den=den+dterm;
		hypmean=hypmean+j*dterm;
		hypvar=hypvar+j**2*dterm;
	end;
	hypmean=hypmean/den;
	hypvar=hypvar/den-hypmean**2;
%mend hypmean;


%macro hypinv(m,n,nn,expv);
	lorlimit=15;
	lamlo=-lorlimit; lamhi=lorlimit;
	if &expv LE max(0,&m+&n-&nn) then lamhi=lamlo;
	if &expv GE min(&n,&m) then lamlo=lamhi;
	do until(abs(lamhi-lamlo)<1e-6);
		lamd=(lamhi+lamlo)/2;
		%hypmean(&m,&n,&nn,lamd);
		if hypmean GE &expv then lamhi=lamd;
		if hypmean LE &expv then lamlo=lamd;
	end;
	temp=hypmean;
%mend hypinv;


%macro fithyp(n,m,nn);
	%hypmean(&m,&n,&nn,_XBETA_);
	mean0=hypmean;

	devden=exp(%hyp(&m,&n,&nn,_RESP_)+_RESP_*_XBETA_)/den;

	%hypinv(&m,&n,&nn,_RESP_);
	var0=hypvar;

	devnum=exp(%hyp(&m,&n,&nn,_RESP_)+_RESP_*lamd)/den;

	devi=2*log(devnum/devden);

	invlink ilink = mean0;
	fwdlink link = log(_MEAN_);
	deviance dev=devi;
	variance vari = var0;

%mend fithyp;




data mle;
	%hypinv(9,16,95,4);	
	lamex=exp(lamd);
run;


proc print data=mle;
	var lamex;
run;


data avadex;
	input expt exposed tumors total strain $ sex $;
	datalines;
	4 16  9  95 X M
	2 14  5 101 X F
	4 18 14 108 Y M
	1 15  4  97 Y F
run;

proc genmod;
	%fithyp(exposed,tumors,total);
	class strain sex;
	model expt = strain sex /obstats maxiter=250 type1 type3;
run;

