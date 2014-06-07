function [resStd,resLU,resEnd,data,itData]=main(path,it)
    %% Daten laden und Condensing durchführen

    itData=load(path);

    data=getData(itData,it);


    %% Standard Condensing ausführen
    %
    subplot(1,3,1);
    resStd=condensing(data);
    KKT=plotKKTMat(resStd);
    title('Standard Condensing');

    %% LU Condensing ausführen
    %
    subplot(1,3,2);
    resLU=condensingLU(data);
    KKT1=plotKKTMat(resLU);
    title('Condensing via LU Zerlegung');
    
    %% Condensing auf Endwert durchführen
    %
    
    resEnd=condensingEnd(data);

    %% Differenz plotten
    %

    subplot(1,3,3)
    M = KKT-KKT1;
    sca = abs(M);
    sca(sca < 1) = 1;
    imagesc(log10(abs(M) ./ sca));

    colorbar('location', 'southoutside');
    title('Differenz');

    

end
