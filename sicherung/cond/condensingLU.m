%% Berechnung des Condensing über die LU Zerlegung
%
function [ res ] = condensingLU( data )

    % Die Reihenfolge der Variablen permutieren
    % v=(s1,s2,....,sndis)
    % w=(p,s0,q0,q1,...,qndis-1)
    
    % Die Matrizen entrsprechend dieser Variablen aufbauen:
    %
    %  |Sm Qm| |v| = |Cm|
    %  |Sc Qc| |w|   |Cc|
    
    %% Die Matrix Sm zusammensetzen
    %
    
    Sm=-1*eye((data.ndis-1)*data.nxd);
    
    for i=2:data.ndis-1
        Sm([(i-1)*data.nxd+1:i*data.nxd],[1+(i-2)*data.nxd:(i-1)*data.nxd])=data.Xs{i};
    end
    
    res.Sm=Sm;
    
    %% Die Matrix Qm zusammensetzen
    %
    
    Qm=zeros((data.ndis-1)*data.nxd,data.np+data.nxd+data.nu*(data.ndis-1));
    
    for i=1:data.ndis-1
        Qm([1+(i-1)*data.nxd:i*data.nxd],[1:data.np])=data.Xp{i};
    end
    
    Qm([1:data.nxd],[1+data.np:data.np+data.nxd])=data.Xs{1};
    
    for i=1:data.ndis-1
        Qm([1+(i-1)*data.nxd:i*data.nxd], ...
            [data.np+data.nxd+1+(i-1)*data.nu:data.np+data.nxd+i*data.nu])=data.Xq{i};
    end
    
    res.Qm=Qm;
    
    %% Die Matrix Sc zusammensetzen
    %
    
    
    [n,~]=size(data.Rs{1});
    
    irow=1+n;
    
    for i=2:data.ndis
        
        [n,~]=size(data.Rs{i});
        
        Sc([irow:irow+n-1],[1+(i-2)*data.nxd:(i-1)*data.nxd])=data.Rs{i};
        
        irow=irow+n;
        
    end
    
    [n,~]=size(data.Cs{1});
    
    for i=2:data.ndis
        Sc([irow:irow+n-1],[1+(i-2)*data.nxd:(i-1)*data.nxd])=data.Cs{i};
    end
    
    res.Sc=Sc;
    
    %% Die Matrix Qc zusammensetzen
    %
    
   
    
    irow=1;
    
    for i=1:data.ndis
        
        [n,~]=size(data.Rp{i});
        
        Qc([irow:irow+n-1],[1:data.np])=data.Rp{i};
        
        irow=irow+n;
       
    end
    
    [n,~]=size(data.Cp);
    
    Qc([irow:irow+n-1],[1:data.np])=data.Cp;
    
    [n,~]=size(data.Rs{1});
    
    Qc([1:n],[1+data.np:data.np+data.nxd])=data.Rs{1};
    
    [n,~]=size(data.Cs{1});
    
    Qc([data.nrdc+1:data.nrdc+n],[data.np+1:data.np+data.nxd])=data.Cs{1};
    
    irow=1;
    
    for i=1:data.ndis-1
        
        [n,~]=size(data.Rq{i});
        
        Qc([irow:irow+n-1],...
            [data.np+data.nxd+(i-1)*data.nu+1:data.np+data.nxd+i*data.nu])=data.Rq{i};
        
        irow=irow+n;
        
        [n,~]=size(data.Cq{1});
        
        Qc([data.nrdc+1:data.nrdc+n],...
            [data.np+data.nxd+(i-1)*data.nu+1:data.np+data.nxd+i*data.nu])=data.Cq{i};
    end
    
    res.Qc=Qc;
    
    %% Den Vektor Cm zusammensetzen
    %
    
    Cm=zeros((data.ndis-1)*data.nxd,1);
    
    for i=1:data.ndis-1
        Cm([1+(i-1)*data.nxd:i*data.nxd],1)=data.Xc{i};
    end
    
    %% Den Vektor Cc zusammenstzen
    %
    

    
    irow=1;
    
    for i=1:data.ndis
        
        [n,~]=size(data.Rc{i});
        
        Cc([irow:irow+n-1])=data.Rc{i};
        
        irow=irow+n;
    end
    
    [n, ~]=size(data.Cc);
    
    Cc([data.nrdc+1:data.nrdc+n])=data.Cc;
    
    res.Cc=Cc;
   
    
    res.Mat=[Sm Qm;
             Sc Qc];
     
    
    res.right=[Cm;Cc'];
    
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
    
        %% Indexvektor fuer Gleichungsnebenbedingungen
    
    res.istate=zeros(data.nrdc+data.nrcc,1);
    row=1;
    
    for i=1:data.ndis
        
        [dim,~]=size(data.Rs{i});
        
        res.istate([row:row+dim-1])=data.iRc{i}';
        
        
        row=row+dim;
    end
    
    res.istate([row:end])=data.iCc;
    
    
    %% Neue Hessematrix bauen
    %
    %
    
    % Variablen bauen
   
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
    
    %% Condensing durchfuehren
    %
    

    res.H=Qm'*(Sm'\Bs)*(Sm\Qm)-Qm'*(Sm'\Bsq')-Bsq*(Sm'\Qm)+Bqq;
    
    res.Eq=Qc-Sc*(Sm\Qm);
    res.vec=Cc'-Sc*(Sm\Cm);
    
    
    
end

