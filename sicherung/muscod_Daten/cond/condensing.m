%% Standard Condensing 

% Die Reihenfolge der Variablen permutieren
% v=(s1,s2,....,sndis)
% w=(p,s0,q0,q1,...,qndis-1)

% Berechnung der Constraints. Matrix hat die Form:
%
% -I | Pm | S0m | Qm
% ------------------ 
%  0 | Pd | S0d | Qd
% ------------------
%  0 | Pc | S0c | Qc

% Output:   -------------------------------
%
%           Eq          Constraint-Matrix
%           vec         Rechte Seite
%           istate      Indexvektor
%                       -> 1 Gleichung
%                       -> 0 Ungleichung
%           bu          Obere Schranke 
%           bl          Untere Schranke
%           H           Hessematrix
%           f           Gradient

function [res] = condensing( data )
    %% Die Matrix Pm bauen :
    %
    % Variablen initialisieren
    temp=data.Xp{1};
    Pm=zeros(data.nxd*(data.ndis-1),data.np);
    
    Pm([1:data.nxd],[1:data.np])=temp;
    
    for m=2:data.ndis-1
        temp=data.Xp{m}+data.Xs{m}*temp;
        Pm([(m-1)*data.nxd+1:m*data.nxd],[1:data.np])=temp;
    end
    
    
    res.Pm=Pm;
    
    %% Die Matrix S0m bauen :
    %
    % Variablen initialisieren
    temp=data.Xs{1};
    S0m=zeros(data.nxd*(data.ndis-1),data.nxd);
    
    % Startelement setzen
    
    S0m([1:data.nxd],[1:data.nxd])=temp;
    
    for m=2:data.ndis-1
        
        temp=data.Xs{m}*temp;
        S0m([(m-1)*data.nxd+1:data.nxd*m],[1:data.nxd])=temp;
        
    end
    
    
    res.S0m=S0m;
    
    %% Die Matrix Qm bauen:
    %
    % Variablen initialisieren
    temp=zeros(data.nxd,data.nu);
    Qm=zeros(data.nxd*(data.ndis-1),data.nu*(data.ndis-1));
    
    for m=1:data.ndis-1
        % Diagonalelemente initialisieren
        temp=data.Xq{m};
        % Indizes initialisieren 
        row =(m-1)*data.nxd+1;
        col =(m-1)*data.nu+1;
        
        % Diagonalelement setzen
        Qm([row:row+data.nxd-1],[col:col+data.nu-1])= temp;
        
        % Füllen der Spalte
        for i= 1:data.ndis -1-m
            temp=data.Xs{m+i}*temp;
            Qm([row+i*data.nxd:row+data.nxd-1+i*data.nxd],[col:col+data.nu-1])=temp;
        end
        
        
    end
    
    res.Qm=Qm;
    
    %% Die Matrix Pd bauen
    %
    % Variabelen initialisieren
    
    row=1;
    
    % Ersten Eintrag setzen
    [dim,~]=size(data.Rp{1});
    Pd([row:row+dim-1],[1:data.np])=data.Rp{1};
    
    row=row+dim;
    
    for i=2:data.ndis
        
        [dim,~]=size(data.Rp{i});

        Pd([row:row+dim-1],[1:data.np])=data.Rp{i}+ ... 
            data.Rs{i}*Pm([(i-2)*data.nxd+1:(i-1)*data.nxd],:);
        
        row=row+dim;
        
    end
    
    res.Pd=Pd;
    
    %% Die Matrix S0d bauen
    %
    %
    
    row=1;
    
    % Ersten Eintrag setzen
    [dim,~]=size(data.Rs{1});
    S0d(row:row+dim-1,1:data.nxd)=data.Rs{1};
    
    row=row+dim;
    
    for i=2:data.ndis
        
        [dim,~]=size(data.Rs{i});
        
        S0d([row:row+dim-1],:)= ...
            data.Rs{i}*S0m([(i-2)*data.nxd+1:(i-1)*data.nxd],:);
        
        row=row+dim;
        
    end
    
    res.S0d=S0d;
    
    %% Die Matrix Qd bauen
    %
    %
    drow=1;
    col=1;
    
    for m=1:data.ndis-1
        
        
        % Diagonalelemente setzen
        [dim,~]=size(data.Rq{m});
        Qd([drow:drow+dim-1],[col:col+data.nu-1])=data.Rq{m};
        row=drow+dim;
        
        % Spalte füllen
        for i=m+1:data.ndis
            
            [dim,~]=size(data.Rs{i});
            Qd([row:row+dim-1],[col:col+data.nu-1])=data.Rs{i}* ...
                Qm([(i-2)*data.nxd+1:data.nxd*(i-1)],[col:col+data.nu-1]);
            row=row+dim;
            
        end
        % Indices setzen
        [dim,~]=size(data.Rq{m});
        drow=drow+dim;
        col=m*data.nu+1;
    end
    
    
    
    res.Qd=Qd;
    
    
    %% Die Matrix Pc bauen
    %
    %
    Pc=data.Cp;
    
    for i=2:data.ndis
        Pc=Pc+data.Cs{i}*Pm([(i-2)*data.nxd+1:(i-1)*data.nxd],:);
    end
    
    
    res.Pc=Pc;
    
    %% Die Matrix S0c bauen
    %
    %
    S0c=data.Cs{1};
    
    for i=2:data.ndis
        S0c=S0c+data.Cs{i}*S0m([(i-2)*data.nxd+1:(i-1)*data.nxd],:);
    end
    
    res.S0c=S0c;

    %% Die Matrix Qc Bauen
    %
    %
    for i=1:data.ndis-1
        
        Qc(:,[(i-1)*data.nu+1:i*data.nu])=data.Cq{i};
        
        for m=2:data.ndis
            Qc(:,[(i-1)*data.nu+1:i*data.nu])=Qc(:,[(i-1)*data.nu+1:i*data.nu])+ ...
                data.Cs{m}*Qm([(m-2)*data.nxd+1:(m-1)*data.nxd],[(i-1)*data.nu+1:i*data.nu]);
        end
    end
    
    res.Qc=Qc;
    
    %% Die Matrix zusammensetzen
    %
    %
    res.Eq=[Pd S0d Qd;
          Pc S0c Qc];
    %% Die Rechte Seite zusammenbauen
    %
    % Der Vektor c hat die Form
    %    |cm|
    %    |--|
    % c= |cd|
    %    |--|
    %    |cc|
    
    %% Den Vektor cm bauen
    cm=zeros((data.ndis-1)*data.nxd,1);
    cm([1:data.nxd])=data.Xc{1};
   
    for i=2:data.ndis-1
        cm([(i-1)*data.nxd+1:i*data.nxd])= cm([(i-1)*data.nxd+1:i*data.nxd]) + ...
            data.Xs{i}*cm([(i-2)*data.nxd+1:(i-1)*data.nxd]);
       
    end
    
    res.cm=cm;
    
    %% Den Vektor cd bauen
    %
    [dim,~]=size(data.Rs{1});
    cd([1:dim])=data.Rc{1};
    
    row=dim+1;
    
    for i=2:data.ndis
        
        [dim,~]=size(data.Rs{i});
        
        cd([row:row+dim-1])=data.Rc{i} + ...
            data.Rs{i}*cm([(i-2)*data.nxd+1:(i-1)*data.nxd]);
        
        row=row+dim;
    end
    
    
    
    res.cd=cd;
    
    %% Den Vektor cc bauen
    %
    %
    cc=data.Cc;
    
    for i=2:data.ndis
        
        cc=cc+data.Cs{i}*cm([(i-2)*data.nxd+1:(i-1)*data.nxd]);
        
    end
    
    res.cc=cc;
    
    %% Den Vektor zusammensetzen
    
    res.vec=[cd' ; cc];
    
    %% Indexvektor fuer Gleichungsnebenbedingungen
    
    res.istate=zeros(data.nrdc+data.nrcc,1);
    row=1;
    
    for i=1:data.ndis
        
        [dim,~]=size(data.Rs{i});
        size(data.iRc{i});
        res.istate([row:row+dim-1])=data.iRc{i}';
        
        
        row=row+dim;
    end
    
    res.istate([row:end])=data.iCc;
    
    %% Die Bounds zusammensetzen
    
    res.bl=zeros(data.np+data.nxd+data.nu*(data.ndis-1),1);
    res.bu=zeros(data.np+data.nxd+data.nu*(data.ndis-1),1);

    res.bl([1:data.np])=data.pbl;
    res.bu([1:data.np])=data.pbu;
    
    res.bl([1+data.np:data.np+data.nxd])=data.xbl{1};
    res.bu([1+data.np:data.np+data.nxd])=data.xbu{1};
    
    for i=1:data.ndis-1
            res.bl([1+data.nxd+data.np+(i-1)*data.nu:data.nxd+data.np+i*data.nu])=...
                data.ubl{i};
            res.bu([1+data.nxd+data.np+(i-1)*data.nu:data.nxd+data.np+i*data.nu])=...
                data.ubu{i};
    end
    
    %% Neue Hessematrix bauen
    %
    %
    
    % Variablen bauen
    Cm= [Pm S0m Qm];
    
    Bs=zeros(data.nxd*(data.ndis-1));
    for i=2:data.ndis
        Bs([(i-2)*data.nxd+1:(i-1)*data.nxd],[(i-2)*data.nxd+1:(i-1)*data.nxd]) ... 
            =data.Bss{i};
    end
    
    Bsq=zeros(data.np+data.nxd+(data.ndis-1)*data.nu,(data.ndis-1)*data.nxd);
    
    for i=2:data.ndis
        
        
        Bsq([1:data.np],[1+(i-2)*data.nxd:(i-1)*data.nxd])=data.Bps{i}';
        
    end
    
    for i=2:data.ndis-1
        
        Bsq([data.np+data.nxd+(i-1)*data.nu+1:data.np+data.nxd+i*data.nu], ...
            [(i-2)*data.nxd+1:(i-1)*data.nxd])=data.Bsq{i};
    end
    
    Bqq=zeros(data.np+data.nxd+(data.ndis-1)*data.nu);
    
    Bqq([1:data.np],[1:data.np])=data.Bpp;
    Bqq([data.np+1:data.np+data.nxd],[1:data.np])=data.Bps{1};
    Bqq([1:data.np],[data.np+1:data.np+data.nxd])=data.Bps{1}';
    
    Bqq([data.np+1:data.np+data.nxd],[data.np+1:data.np+data.nxd])=data.Bss{1};
    Bqq([data.np+data.nxd+1:data.np+data.nxd+data.nu],[data.np+1:data.np+data.nxd]) ...
        =data.Bsq{1};
    Bqq([data.np+1:data.np+data.nxd],[data.np+data.nxd+1:data.np+data.nxd+data.nu]) ...
        =data.Bsq{1}';
    
    for i=1:data.ndis-1
        Bqq([1+data.np+data.nxd+(i-1)*data.nu:data.np+data.nxd+i*data.nu], ... 
            [1+data.np+data.nxd+(i-1)*data.nu:data.np+data.nxd+i*data.nu])=data.Bqq{i};
    end
    
    
    %% Neuen Gradienten bauen
    %
    %
    bs=zeros(data.nxd*(data.ndis-1),1);
    
    for i=2:data.ndis
        bs([1+(i-2)*data.nxd:(i-1)*data.nxd])=data.fs{i};
    end
    
    bq=zeros(data.nxd+data.np+(data.ndis-1)*data.nu,1);
    
    bq([1:data.np])=data.fp;
    bq([1+data.np:data.np+data.nxd])=data.fs{1};
    res.bs=bs;
    for i=1:data.ndis-1
        bq([1+data.np+data.nxd+(i-1)*data.nu:data.np+data.nxd+i*data.nu])=data.fq{i};
    end
    res.bq=bq;
    
    %% Condensing ausfuehren
    res.f=Cm'*bs+bq-Cm'*Bs*cm -Bsq*cm;
    res.H=Cm'*Bs*Cm+Cm'*Bsq'+Bsq*Cm+Bqq;
end

