%% Berechnung der Zulässigkeit von durch Condensing erzeugte Polyeder
%


function []=checkFeasible(path, it)

    %% Daten einlesen und Condensing durchführen
%     itData=load(path);
% 
%     data=getData(itData,it);
% 
%     resStd=condensing(data);
%     resLU=condensingLU(data);
%     resEnd=condensingEnd(data);
% 
% 
% 
%     % Allgemeine linprog Optionen
% 
%     options = optimoptions('linprog','MaxIter',10000,'Display','iter');
%


%% Berechnung der Zulässigkeit von resStd
%
    
    % Normierung der Zeilen von A
    
    
    
%     [rows,col]=size(resStd.Eq);
%     
%     A=diag(((resStd.Eq.^2)*ones(col,1)).^(-1/2))*resStd.Eq;
%     b=diag(((resStd.Eq.^2)*ones(col,1)).^(-1/2))*resStd.vec;
%     
%     
%     
%     % Berechnung der Zulässigkeit
%     
%     [x1,fval,exitflag,output,lambda] = ...
%         linprog([zeros(col,1); ones(rows+2*col,1)] , ...
%         [eye(col) zeros(col,rows) -eye(col) zeros(col); ...
%          -eye(col) zeros(col,rows) zeros(col) -eye(col)], ...
%         [resStd.bu; -resStd.bl], ...
%         [A eye(rows) zeros(rows,col) zeros(rows,col)], ... 
%         [b], ...
%         [-10e10*ones(col+rows,1) ; zeros(2*col,1)], ...
%         [10e10*ones(3*col+rows,1)], ...
%         [zeros(col,1) ; b ; abs(resStd.bu); abs(resStd.bl)], ...
%         options);
% 
%     
%     % Berechnung der Zulässigkeit von resLu
%     
%     [rows,col]=size(resLU.Eq);
%     
%     A=diag(((resLU.Eq.^2)*ones(col,1)).^(-1/2))*resLU.Eq;
%     b=diag(((resLU.Eq.^2)*ones(col,1)).^(-1/2))*resLU.vec;
%     
%     
%     [x2,fval,exitflag,output,lambda] = ...
%             linprog([zeros(col,1); ones(rows+2*col,1)] , ...
%             [eye(col) zeros(col,rows) -eye(col) zeros(col); ...
%              -eye(col) zeros(col,rows) zeros(col) -eye(col)], ...
%             [resLU.bu; -resLU.bl], ...
%             [A eye(rows) zeros(rows,col) zeros(rows,col)], ... 
%             [b], ...
%             [-10e10*ones(col+rows,1) ; zeros(2*col,1)], ...
%             [10e10*ones(3*col+rows,1)], ...
%             [zeros(col,1) ; b ; abs(resLU.bu); abs(resLU.bl)], ...
%             options);
% 
%     % Berechnung der Zulässigkeit von resEnd
%     
%     [rows,col]=size(resEnd.Eq);
%     
%     A=diag(((resEnd.Eq.^2)*ones(col,1)).^(-1/2))*resEnd.Eq;
%     b=diag(((resEnd.Eq.^2)*ones(col,1)).^(-1/2))*resEnd.vec;
%     
%     [x3,fval,exitflag,output,lambda] = ...
%         linprog([zeros(col,1); ones(rows+2*col,1)] , ...
%         [eye(col) zeros(col,rows) -eye(col) zeros(col); ...
%         -eye(col) zeros(col,rows) zeros(col) -eye(col)], ...
%         [resEnd.bu; -resEnd.bl], ...
%         [A eye(rows) zeros(rows,col) zeros(rows,col)], ... 
%         [b], ...
%         [-10e10*ones(col+rows,1) ; zeros(2*col,1)], ...
%         [10e10*ones(3*col+rows,1)], ...
%         [A eye(rows) zeros(rows,col) zeros(rows,col)], ...
%         options);
%        
%     
%     [A eye(rows) zeros(rows,col) zeros(rows,col)]*[A eye(rows) zeros(rows,col) zeros(rows,col)]
%   

%% Berechnung der Zulässigkeit des von Muscod condensierten Problems
    
    
    temp=importdata('../walker2d_Daten/lowerbound_it1.txt');
    bl=temp.data;
    temp=importdata('../walker2d_Daten/upperbound_it1.txt');
    bu=temp.data;

    cm=importdata('../walker2d_Daten/constraintMatrix_it1.txt');
    
    [rows, col]=size(cm);

    % Bounds aufsplitten

    bl1=bl(1:col);
    bl2=bl(col+1:end);

    bu1=bu(1:col);
    bu2=bu(col+1:end);

    % Optimierungsvariablen definieren und Normieren



    %A=diag(((cm.^2)*ones(col,1)).^(-1/2))*cm;
    %bl2=diag(((cm.^2)*ones(col,1)).^(-1/2))*bl2;
    %bu2=diag(((cm.^2)*ones(col,1)).^(-1/2))*bu2;


    Aineq=[ eye(col)  zeros(col,rows)   zeros(col, rows);
            -eye(col) zeros(col,rows)   zeros(col,rows);
             cm        -eye(rows)        zeros(rows);
            -cm        zeros(rows)       -eye(rows)];

    bineq=[bu1;-bl1;bu2;-bl2];

    f=[zeros(col,1) ; ones(2*rows,1)];

%     nonnegativ=[-10e20*ones(col,1); zeros(2*rows,1)];
    nonnegativ=[-10e20*ones(col,1); zeros(2*rows,1)];

    keyboard  
    [contain ok]=BuildMPS(Aineq,bineq,[],[],f,nonnegativ,[])
    
      

    options =optimoptions('linprog','Algorithm','simplex','MaxIter',10000,'Display','iter');



    [x,fval,exitflag,output,lambda] =linprog( ...
            f, ...
            Aineq, ...
            bineq, ...
            [], ...
            [], ...
            nonnegativ, ...
            [], ...
            [], ...
            options);
     keyboard   

end