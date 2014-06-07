clear all;
close all;

%Solverpath setzen

addpath(genpath('~/solver'));

%Daten laden und formatieren
%itData=load('RES/lotkaTest.mat');
itData=load('../walker2d_Daten/walker2d.mat');

data=getData(itData,1);

% Condensing durchführen

resEnd=condensingEnd(data);

% Gleichungungen und Ungleichungen aufsplitten

cmeq=resEnd.Eq(resEnd.istate==1,:);
cmineq=resEnd.Eq(resEnd.istate==0,:);

veceq=resEnd.vec(resEnd.istate==1);
vecineq=resEnd.vec(resEnd.istate==0);

[~,dimx]=size(resEnd.Eq);
[dimd1,~]=size(cmeq);
[dimd3,~]=size(cmineq);

% Konfiguration und Optimierung

opt=sdpsettings('solver','lpsolve','verbose',10,'warning',1);

x = sdpvar(dimx,1);
d1= sdpvar(dimd1,1);
d2= sdpvar(dimd1,1);
d3= sdpvar(dimd3,1);

Constraints = [ resEnd.bl<=x<=resEnd.bu , ...
                cmeq*x<=veceq+d1, ...
                d2+cmeq*x>=veceq, ...
                cmineq*x>=vecineq+d3, ...
                d1>=0, d2>=0, d3>=0];

obj=sum(d1)+sum(d2)+sum(d3);
solvesdp(Constraints,obj,opt)
