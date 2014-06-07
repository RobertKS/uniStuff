clear all;
close all;

%itData=load('RES/lotkaTest.mat');
itData=load('../walker2d_Daten/walker2d.mat');

it=3;

data=getData(itData,it);

%% Als Referenz und zum Abgleich der Ergebnisse uncondensiertes Problem mit quadprog Lösen

H=itData.iterations{it,1}.qp.hessian;

f=itData.iterations{it,1}.qp.gradient;

A=itData.iterations{it,1}.qp.constraints.matrix;


iVektor=itData.iterations{it,1}.qp.istateStart;

iVektor=iVektor([data.nvar+1:end]);

Aeq=A(iVektor==3,:);

A=A(iVektor~=3,:);

b=itData.iterations{it,1}.qp.constraints.vector;

beq=b(iVektor==3);

b=b(iVektor~=3);

lb=itData.iterations{it,1}.qp.bounds.lower;

ub=itData.iterations{it,1}.qp.bounds.upper;

opts = optimoptions('quadprog','Algorithm','active-set','Display','none');

x0=itData.iterations{it,1}.qp.solution.value;

fprintf('**************************************************\n');
display('Loese uncondensiertes QP...');
fprintf('**************************************************\n');
tic;
x_std=quadprog(H,f,-A,-b,Aeq,beq,lb,ub,[],opts);
toc;
fprintf('\n\n');

%% Problem in condensierter Form loesen

fprintf('**************************************************\n');
display('Loese condensiertes QP...');
fprintf('**************************************************\n');


% Condensing durchfuehren
timecom=tic;
res=condensingQR(data);
timecon=toc(timecom);
display('Loese QP');

%QP Loesen
timesol=tic;
x_qr=quadprog( ...
    res.H, ...
    res.f, ...
    res.Con_ineq, ...
    res.c_ineq, ...
    res.Con_eq, ...
    res.c_eq, ...
    [], ...
    [], ...
    [], ...
    opts);
timesol=toc(timesol);


display('Transformiere Loesung zurueck');

% Abhaengigie Variablen berechnen
v=res.s_t-res.S_t*x_qr;
sol=[v;x_qr];
x_qr1=res.T'*sol;

x_qrt=x_qr1;

% Loesung umsortieren
for i=1:data.ndis
        
        x_qrt([1+data.np+(i-1)*(data.nu+data.nxd)+data.nu:data.np+i*(data.nu+data.nxd)])=...
            x_qr1([1+data.np+(i-1)*(data.nu+data.nxd):data.np+(i-1)*(data.nu+data.nxd)+data.nxd]);
        
        x_qrt([1+data.np+(i-1)*(data.nu+data.nxd):data.np+(i-1)*(data.nu+data.nxd)+data.nu])=...
            x_qr1([1+data.np+(i-1)*(data.nu+data.nxd)+data.nxd:data.np+i*(data.nu+data.nxd)]);
end
timecom=toc(timecom);

fprintf('**************** Zeitmessung **************\n Condensing : %f \n Solver : %f \n Gesamt : %f \n', ...
    timecon,timesol,timecom);

fprintf('*******************************************\n Differenz der Loesungen : %f \n', norm(x_std-x_qrt));


x0_t=x0;

for i=1:data.ndis
        
        x0_t([1+data.np+(i-1)*(data.nu+data.nxd):data.np+(i-1)*(data.nu+data.nxd)+data.nxd])=...
            x0([1+data.np+(i-1)*(data.nu+data.nxd)+data.nu:data.np+i*(data.nu+data.nxd)]);
        
        x0_t([1+data.np+(i-1)*(data.nu+data.nxd)+data.nxd:data.np+i*(data.nu+data.nxd)])=...
            x0([1+data.np+(i-1)*(data.nu+data.nxd):data.np+(i-1)*(data.nu+data.nxd)+data.nu]);
end