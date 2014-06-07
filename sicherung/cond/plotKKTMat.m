function [KKT] = plotKKTMat( res )

    %%KKT Matrix bauen
    
    
    [m,~]=size(res.Eq);
    KKT=[res.H res.Eq';
         res.Eq zeros(m) ];
    
    
    imagesc(log10(abs(KKT)));
    
    colorbar('location', 'southoutside');

end

