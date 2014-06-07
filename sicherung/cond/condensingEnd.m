%% Condensing des QPs über Abhängigkeit von Endwert
%
function [res]=condensingEnd(data)

    % Die Reihenfolge der Variablen permutieren
    % v=(sndis-1,sndis-2, ... , s1,s0)
    % w=(p,sndis,q0,q1,...,qndis-1)
    
    % Die Reihenfolge der Matching Conditions umdrehen:
    % Cm=(cndis,cndis-1,....,c1,c0)^T
    
    % Die Matrizen entrsprechend dieser Variablen aufbauen:
    %
    %  |Sm Qm| |v| = |Cm|
    %  |Sc Qc| |w|   |Cc|
    
    %% Die Matrix Sm zusammmensetzen
    %
    
    Sm=diag(-1*ones((data.ndis-2)*data.nxd,1),-data.nxd);
    
    for i=1:data.ndis-1
        Sm([1+(i-1)*data.nxd:i*data.nxd],[1+(i-1)*data.nxd:i*data.nxd]) ...
            =data.Xs{data.ndis-i};
    end
    res.Sm=Sm;
    
    %% Die Matrix Qm zusammensetzen
    %
    
    Qm=zeros((data.ndis-1)*data.nxd,data.np+data.nxd+data.nu*(data.ndis-1));
    
    for i=1:data.ndis-1
        Qm([1+(i-1)*data.nxd:i*data.nxd],[1:data.np])=data.Xp{data.ndis-i};
    end
    
    Qm([1:data.nxd],[1+data.np:data.np+data.nxd])=-1*eye(data.nxd);
    
    for i=1:data.ndis-1
        Qm([1+(i-1)*data.nxd:i*data.nxd], ...
           [1+data.nxd+data.np+(data.ndis-i-1)*data.nu: ... 
           data.nxd+data.np+(data.ndis-i)*data.nu]) = data.Xq{data.ndis-i};
    end
    
    
    res.Qm=Qm;
    
    %% Die Matrix Sc zusammensetzen
    %
    
    Sc=zeros(data.nrdc+data.nrcc,(data.ndis-1)*data.nxd);
    
    irow=1;
    
    for i=1:data.ndis -1
        
        [dim,~]=size(data.Rs{i});
        
        Sc([irow:irow+dim-1],...
           [1+(data.ndis-1-i)*data.nxd:(data.ndis-i)*data.nxd])...
           =data.Rs{i};
        
        irow=irow+dim;
    end
    
    for i=1:data.ndis-1
        Sc([data.nrdc+1:data.nrdc+data.nrcc],[1+(i-1)*data.nxd:i*data.nxd])...
            =data.Cs{data.ndis-i};
    end
    
    res.Sc=Sc;
    
    %% Die Matrix Qc zusammensetzen
    %
    
    Qc=zeros(data.nrdc+data.nrcc,data.np+data.nxd+(data.ndis-1)*data.nu);
    
    irow=1;
    
    for i=1:data.ndis
        
        [dim,~]=size(data.Rp{i});
        
        Qc([irow:irow+dim-1],[1:data.np])=data.Rp{i};
        
        irow=irow+dim;
    end
    
    Qc([data.nrdc+1:data.nrdc+data.nrcc],[1:data.np])=data.Cp;
    
    [dim,~]=size(data.Rs{data.ndis});
    
    Qc([data.nrdc-dim+1:data.nrdc],[data.np+1:data.np+data.nxd])=data.Rs{data.ndis};
    Qc([1+data.nrdc:data.nrdc+data.nrcc],[1+data.np:data.nxd+data.np])...
        =data.Cs{data.ndis};
    
    irow=1;
    
    for i=1:data.ndis-1
        
        [dim,~]=size(data.Rp{i});
        
        Qc([irow:irow+dim-1],...
            [1+data.np+data.nxd+(i-1)*data.nu:data.np+data.nxd+i*data.nu])...
            =data.Rq{i};
        
        irow=irow+dim;
        
    end
    
    for i=1:data.ndis-1
        Qc([data.nrdc+1:data.nrdc+data.nrcc], ... 
            [1+data.np+data.nxd+(i-1)*data.nu:data.np+data.nxd+i*data.nu]) ...
            =data.Cq{i};
    end
    
    res.Qc=Qc;
    
    %% Die rechte Seite Cm zusammensetzen
    % 
    Cm=zeros((data.ndis-1)*data.nxd,1);
    
    for i=1:data.ndis-1 
        Cm([1+(i-1)*data.nxd:i*data.nxd],1)=data.Xc{data.ndis-i};
    end
    
    res.Cm=Cm;
    
    %% Die rechte Seite Cc zusammensetzen
    % 
    
    Cc=zeros(data.nrdc+data.nrcc,1);
    
    irow =1;
    
    for i =1:data.ndis
        
        [dim,~]=size(data.Rc{i});
        
        Cc([irow:irow+dim-1],1)=data.Rc{i};
        
        irow=irow+dim;
        
    end
    
    Cc([data.nrdc+1:data.nrdc+data.nrcc])=data.Cc;
    

    %% Die Bounds zusammensetzen
    
    res.bl=zeros(data.np+data.nxd+data.nu*(data.ndis-1),1);
    res.bu=zeros(data.np+data.nxd+data.nu*(data.ndis-1),1);

    res.bl([1:data.np])=data.pbl;
    res.bu([1:data.np])=data.pbu;
    
    res.bl([1+data.np:data.np+data.nxd])=data.xbl{data.ndis};
    res.bu([1+data.np:data.np+data.nxd])=data.xbu{data.ndis};
    
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
    
    %% Condensing durchfuehren
    %   
    
    res.Eq=Qc-Sc*(Sm\Qm);
    res.vec=Cc-Sc*(Sm\Cm);
end    
