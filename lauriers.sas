FILENAME REFFILE '/folders/myfolders/lauriers.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.LAURIERS;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.LAURIERS; RUN;

PROC REG DATA=LAURIERS;
	MODEL taille=masse;
	OUTPUT OUT=WORK.RESIDUSLAURIERS RESIDUAL=RES PREDICTED=PRED;
RUN;

PROC UnIvArIaTe DATA=work.residuslauriers normaltest;
	Var RES;
RUN;

*Absence de normalité;

/*
Commentaires
sur
quatre
lignes
*/

PROC PRINT DATA=RES;
RUN;


