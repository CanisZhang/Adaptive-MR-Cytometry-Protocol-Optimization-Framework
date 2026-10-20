function fitout = PostProcessIMPULSED(fitout)
%% post-process IMPULSED fitout to include rho (3D cell density), cellularity, and cellularity index (fin/d)

    %% post processing of MRI cytometry 
    % ------------------------rho (3D cell density: number of cells per unit vollume )------------------------
    fitout.mapLabels = [fitout.mapLabels, "rho"] ; 
    fitout.rho = single(fitout.data.mask) ; 
    fitout.rho(fitout.data.mask) = 6/pi * fitout.vin(fitout.data.mask) ./ (fitout.d(fitout.data.mask)).^3 ; 

    % ------------------------cellularity (2D: number of cells per unit area) ------------------------
    fitout.mapLabels = [fitout.mapLabels, "cellularity"] ; 
    fitout.cellularity = single(fitout.data.mask) ; 
    fitout.cellularity(fitout.data.mask) = 2* (3/2/pi * fitout.vin(fitout.data.mask)).^(2/3) ./ (fitout.d(fitout.data.mask)).^2 ;  

    % ------------------------cellularity index (vin/d)------------------------
    fitout.mapLabels = [fitout.mapLabels, "cellularity_index"] ; 
    fitout.cellularity_index = single(fitout.data.mask) ; 
    fitout.cellularity_index(fitout.data.mask) = fitout.vin(fitout.data.mask) ./ fitout.d(fitout.data.mask) ; 

end
