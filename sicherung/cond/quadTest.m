

itData=load('../walker2d_Daten/walker2d.mat');


H=itData.iterations{1,1}.qp.hessian;
f=itData.iterations{1,1}.qp.gradient;
Aeq=itData.iterations{1,1}.qp.constraints.matrix;
beq=itData.iterations{1,1}.qp.constraints.vector;
lb=itData.iterations{1,1}.qp.bounds.lower;
ub=itData.iterations{1,1}.qp.bounds.upper;



x_quad = quadprog(H,f,[],[],Aeq,beq,lb,ub);