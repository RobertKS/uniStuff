clear all;
close all;

%Solverpath setzen

addpath(genpath('~/solver'));

%Daten importieren und foramtieren


%walker Daten
temp=importdata('../walker2d_Daten/lowerbound_it1.txt');
bl=temp.data;
temp=importdata('../walker2d_Daten/upperbound_it1.txt');
bu=temp.data;

M=importdata('../walker2d_Daten/constraintMatrix_it1.txt');

%lotka Daten

% temp=importdata('../lotkaTest_Daten/lowerbound_it1.txt');
% bl=temp.data;
% temp=importdata('../lotkaTest_Daten/upperbound_it1.txt');
% bu=temp.data;
% 
% M=importdata('../lotkaTest_Daten/constraintMatrix_it1.txt');





[rows, col]=size(M);

% Bounds aufsplitten

bl1=bl(1:col);
bl2=bl(col+1:end);

bu1=bu(1:col);
bu2=bu(col+1:end);


% Zeilen weglassen

N=0;

bl2=bl2([1:end-N]);
bu2=bu2([1:end-N]);

cm=M([1:end-N],:);

% Konfiguration und Optimierung

opt=sdpsettings('solver','lpsolve','verbose',10,'warning',1);

x = sdpvar(col,1);
d1= sdpvar(rows-N,1);
d2= sdpvar(rows-N,1);

Constraints = [ bl1<=x<=bu1 , ...
                cm*x <= bu2 + d1 , ...
                bl2-d2<=cm*x , ...
                d1 >=0 , ...
                d2 >= 0];

obj=sum(d1)+sum(d2);
solvesdp(Constraints,obj,opt)


