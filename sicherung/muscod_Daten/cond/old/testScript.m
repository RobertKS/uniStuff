clc;
clear all;
close all;

[resStd,resLU,resEnd,data, itData]=main('../walker2d_Daten/walker2d.mat', 1);

temp=load('../walker2d_Daten/constraintMatrix.txt');

temp=temp([1:153],[1:298]);

div=abs(temp);
div(abs(temp) < 1e-10)=1;

figure(2)
temp=log10(abs((temp-resStd.Eq)./div));
temp(temp <=0.01)=0;
imagesc(temp)