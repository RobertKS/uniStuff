%% Script zum lösen von Muscod QP via qpDunes
%
%
clc
clear all;
close all;

% Daten einlesen
%path='RES/lotkaTest.mat';
path='../walker2d_Daten/walker2d.mat';
itData=load(path);
data=getData(itData, 1);

% Pfad setzen


addpath(genpath('~/solver'));


% Dimensionen

nI  = data.ndis-1;          % number of control intervals
nX  = data.nxd+data.np;     % number of states   %Parameter als Zustaende 
nU  = data.nu;              % number of controls
nZ  = nX+nU;
nD = zeros(nI+1,1);         % number of affine constraints



% H Matrizen zusammensetzen

H=zeros(nZ,nZ*nI);
g=zeros(nZ*nI+nX,1);


for i=1:nI
    
   H(:,[1+(i-1)*nZ:i*nZ])= ...
        [data.Bss{i}     data.Bps{i}  data.Bsq{i}';
         data.Bps{i}'    data.Bpp     data.Bpq{i}';
         data.Bsq{i}     data.Bpq{i}  data.Bqq{i}];
    
    g([1+(i-1)*nZ:i*nZ],1)=[data.fs{i}; data.fp; data.fq{i}];
     
end

P=[data.Bss{data.ndis} data.Bps{data.ndis};
    data.Bps{data.ndis}' data.Bpp];

g(nZ*nI+1:end)=[data.fs{data.ndis}; data.fp];

% Matching Condition zusammensetzen

C=zeros(nX,nZ*nI);
c=zeros(nX*nI,1);

for i=1:nI
    
    C(:,[1+(i-1)*nZ:i*nZ])= ...
        [data.Xs{i}                 data.Xp{i}     data.Xq{i};
         zeros(data.np,data.nxd)    eye(data.np)   zeros(data.np,data.nu)];
    
     c([1+(i-1)*nX:i*nX],1)=[data.Xc{i}; zeros(data.np,1)];
     
end


% Bounds setzen

zLow=zeros(nI*nZ+nX,1);
zUpp=zeros(nI*nZ+nX,1);

for i=1:nI
    
    zLow([1+(i-1)*nZ:i*nZ],1)=[data.xbl{i}' ; data.pbl' ; data.ubl{i}'];
    zUpp([1+(i-1)*nZ:i*nZ],1)=[data.xbu{i}' ; data.pbu' ; data.ubu{i}'];
end

zLow([1+nI*nZ:end],1)=[data.xbl{data.ndis}' ; data.pbl'];
zUpp([1+nI*nZ:end],1)=[data.xbu{data.ndis}' ; data.pbu'];

% SOLVE

qpOptions = qpDUNES_options( 'default', ...
                             'maxIter', 100, ...
                             'printLevel', 2, ...
                             'logLevel', 0, ...     % log all data
                             'lsType', 4, ...       % Accelerated gradient biscection LS
                             'stationarityTolerance', 1.e-6, ...
                             'regType', 2 ...       % regularize only singular directions; 1 is normalized Levenberg Marquardt
                             ...
                             );


% a) initialize data
qpDUNES( 'init', nI, ...
          H, P, g, ...
          C, -c, ...
          zLow, zUpp, ...
          [], [], [], ...        % affine constraints not yet supported
          qpOptions );

[zOpt, stat, lambda, mu, objFctnVal] = qpDUNES( 'solve' );



