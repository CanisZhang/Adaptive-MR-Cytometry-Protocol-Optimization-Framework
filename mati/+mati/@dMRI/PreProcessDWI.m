function [out, fitout] = PreProcessDWI(fitpars, data)
%% mati.IMPULSED.PreProcessDWI to either normalize or remove IVIM effects in diffusion weighted MRI data
% If fitopts.flag.getADC == 'y', ADC(tdiff) will be calculated as well. 
% INPUTS:
%       fitpars: a mati.FitPars object
%       data: a mati.ImageData object
% OUTPUTS:
%       out: a mati.ImageData object
% -------------------------------------------------------------------------------------------------------------------------
% Author: Junzhong Xu
%
% -------------------------------------------------------------------------------------------------------------------------

    %% preliminary
    warning('off') ; 
    
    % check input 
    if ~isa(data,'mati.ImageData'), error('%s: the input data should be an ImageData object',mfilename) ; end
    
    % preparation 
    Img = data.img ; 
    pulse = fitpars.model.pulse ; 
    
    %% check tdiff
    tdiffs = unique(pulse.tdiff) ; tdiffs = sort(tdiffs, 'descend') ; Ntdiff = length(tdiffs) ; 
    % temporary Img reshaped for computing convenience
    tmpImgAll = reshape(Img, [data.Nx*data.Ny*data.Nz, data.Nt]) ; 
    % initialization of fitout of ADC 
    if fitpars.fitopts.flag.getADC == 'y'
        fitout.tdiff = tdiffs ; 
    end
    % deivim or normalization 
    for ntdiff = 1:Ntdiff
        ind_tdiff = (pulse.tdiff==tdiffs(ntdiff)) ; 
        ind_fit = ind_tdiff & (((pulse.b >= fitpars.fitopts.ADCfit_bmin) & (pulse.b <= fitpars.fitopts.ADCfit_bmax)) | pulse.n~=0) ; 
        ind_b0 = (pulse.b<fitpars.fitopts.S0_bmax) & ind_tdiff ;       % indices of all b=0 of the nth tdiff 

        if fitpars.fitopts.flag.deivim ~= 'y'
            % ----------------------------------- regular normalization -----------------------------
            if sum(ind_b0) == 0 % no b=0 for nth tdiff
                warning('%s: Cannot find a b0 image with tdiff = %f ms. Signal will be normalized using all b=0 images. It is highly recommended to acquire at last one b=0 image for each diffusion time (gradietn shape) for normalization purpose. ',mfilename, tdiffs(ntdiff)) 
                ind_b0 =  (pulse.b<fitpars.model.defaultFitopts.S0_bmax) ; 
                if sum(ind_b0) == 0 % no b=0 at all
                    error('%s: No b=0 (<%f) images can be found at all!',mfilename, fitpars.model.defaultFitopts.S0_bmax)
                end
            end
            % normalization
            S0 = mean(tmpImgAll(data.mask, ind_b0), 2) ; 
            tmpImgAll(data.mask, ind_tdiff) = tmpImgAll(data.mask, ind_tdiff) ./ S0 ; 
            tmpImgAll(~data.mask, ind_tdiff) = tmpImgAll(~data.mask, ind_tdiff) ./ mean(S0) ; 
         else 
            % -------------------------------------- deivim ------------------------------------------
            switch sum(ind_fit)
                case {0,1}      % not enough points to perform ADC fitting so use regular normalization instead
                    warning on ; warning('%s: Not sufficient points for ADC fitting to remove IVIM with tdiff=%.2f ms. The regular normalization is used instead. ',mfilename, tdiffs(ntdiff)) ; warning off ; 
                    if sum(ind_b0) == 0 % no b=0 for nth tdiff so use all b=0 images instead
                        warning('%s: Cannot find a b0 image with tdiff = %.2f ms. Signal will be normalized using all b=0 images. It is highly recommended to acquire at last one b=0 image for each diffusion time (gradietn shape) for normalization purpose. ',mfilename, tdiffs(ntdiff)) 
                        ind_b0 =  (pulse.b<fitpars.fitopts.S0_bmax) ; 
                        if sum(ind_b0) == 0 % no b=0 at all
                            error('%s: No b=0 (<%.2f) images with tdiff=.2f ms can be found at all!',mfilename, fitpars.fitopts.S0_bmax, tdiffs(ntdiff))
                        end
                    end
                    % normalization
                    S0 = mean(tmpImgAll(data.mask, ind_b0), 2) ; 
                    tmpImgAll(data.mask, ind_tdiff) = tmpImgAll(data.mask, ind_tdiff) ./ S0 ; 
                    tmpImgAll(~data.mask, ind_tdiff) = tmpImgAll(~data.mask, ind_tdiff) ./ mean(S0) ; 
                otherwise  % fit -log(signal) ~ -b*D to get S(b=0)
                    % fitting all ind_fit
                    tmpImg = tmpImgAll(data.mask, ind_fit) ; 
                    % mono-exponential fit
                    [tmpP] =arrayfun(@(n) polyfit(pulse.b(ind_fit), -log(tmpImg(n,:)), 1), [1:size(tmpImg,1)], 'UniformOutput', false) ; 
                    tmpS0 = cellfun(@(x) x(2), tmpP) ; tmpS0 = exp(-tmpS0) ; 

                    % scale all tmpImgAll voxels first
                    tmpImgAll(data.mask(:), ind_tdiff) = tmpImgAll(data.mask(:), ind_tdiff) ./ tmpS0(:) ; 
                    tmpImgAll(~data.mask(:), ind_tdiff) = tmpImgAll(~data.mask(:), ind_tdiff) ./ mean(tmpS0(:)) ;  
            end 
        end
    end % ntdiff

    Img = reshape(tmpImgAll, [data.Nx, data.Ny, data.Nz, data.Nt]) ; 
    % remove nan if any
    Img(isnan(Img)) = 0 ; 

    %% update image object
    out = data.UpdateImg(Img) ; 
    
end  % end of function
