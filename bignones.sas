FILENAME REFFILE '/folders/myfolders/bignones.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.bignones;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.bignones; RUN;

PROC REG DATA=bignones ;
	MODEL masse=masse_sec / SPEC HCCMETHOD=3;
	*Defaut d'homogeneite des variances. La methode HCCMETHOD=3 est
	recommandee si la taille de l'echantillon est inferieure 
	ou egale à 250.;
	OUTPUT OUT=WORK.RESIDUSbignones RESIDUAL=RES PREDICTED=PRED;
RUN;

PROC UnIvArIaTe DATA=work.residusbignones normaltest;
	Var RES;
RUN;