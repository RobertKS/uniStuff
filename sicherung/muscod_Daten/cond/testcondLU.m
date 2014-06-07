clear all;
close all;

%Solverpath setzen

addpath(genpath('~/solver'));

%Daten laden und formatieren


itData=load('RES/lotkaTest.mat');
%itData=load('../walker2d_Daten/walker2d.mat');

data=getData(itData,1);

% Condensing durchführen



resLU=condensingLU(data);

% Gleichungungen und Ungleichungen aufsplitten

cmeq=resLU.Eq(resLU.istate==1,:);
cmineq=resLU.Eq(resLU.istate==0,:);

veceq=resLU.vec(resLU.istate==1);
vecineq=resLU.vec(resLU.istate==0);

[~,dimx]=size(resLU.Eq);
[dimd1,~]=size(cmeq);
[dimd3,~]=size(cmineq);

% Konfiguration und Optimierung

opt=sdpsettings('solver','lpsolve','verbose',10,'warning',1);

x = sdpvar(dimx,1);
d1= sdpvar(dimd1,1);
d2= sdpvar(dimd1,1);
d3= sdpvar(dimd3,1);

Constraints = [ resLU.bl<=x<=resLU.bu , ...
                cmeq*x<=veceq+d1, ...
                d2+cmeq*x>=veceq, ...
                cmineq*x>=vecineq+d3, ...
                d1>=0, d2>=0, d3>=0];

obj=sum(d1)+sum(d2)+sum(d3);
solvesdp(Constraints,obj,opt)
