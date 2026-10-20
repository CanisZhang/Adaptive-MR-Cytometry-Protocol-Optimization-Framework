function E = RestrictedDWISignal(parms, model)
%% Calculates analytical diffusion MRI signal attenuation of restricted diffusion 
% The geometry must be a member of {plane, cylinder, sphere, or spherical shell}.
% The separation of pulse sequence and microtructure parameters is based on "Janez Stepisnik, Physica tp, 183, 343-350 (1993)"
% The signal equations of sin or cos and microstructure expressions of spherical shell are based on " Xu J et al. Quantitative characterization of tissue microstructure with temporal diffusion spectroscopy. J Magn Reson 2009;200(2):189-197."
% The signal equations of tpgse or tcos shapes are based on "Xu J et al. Magnetic resonance imaging of mean cell size in human breast tumors. Magn Reson Med. 2020 Jun;83(6):2002-2014."
% 
% -------------------------------------------------------------------------------------------------------------------------
% INPUTS: 
%       model:
%           model.pulse: pulse sequence parameters. The elements such as f, T, delta, tau, w should be scalars or matrices with the same dimensions
%               model.pulse.shape: pgse, tpgse, sin, cos, tcos, tsin
%           model.structure: microstructural parameters: geometry can be plane, cylinder, sphere, or spherical shell
%       parms: 
%           for {'plane','cylinder','sphere'}, parms is a 2 x N marix [d; D]
%           for {'hollowSphere', 'sphericalShell'}, parms is a 3 x N matrix [din; dout; D]
% OUTPUTS:
%       E:   signal attenuation with the same dimension of pulse elements
% -------------------------------------------------------------------------------------------------------------------------
% NOTES:
%        1. For plane: d is the distance between planes
%        2. For cylinder or sphere: d is the diameter
%        3. All inputs are d (diameter)
%        4. The model.pulse can take any combinations of pulse sequences (pgse, tpgse, sin, cos, tcos, ...)
%        5. model.pulse.shape == cos is NOT apodised but an ideal cosine shape. It provides little bias for fitting size. 
%        6. xstep in the function get_hollow_sphere_roots(a,b) should be small enough to count for a small din (inner diameter)
% -------------------------------------------------------------------------------------------------------------------------
% Author: Junzhong Xu, Dec 08, 2020; Modified by Zhiyuan Zhang 
% -------------------------------------------------------------------------------------------------------------------------

    %% Check dimensions: d and D should be either scalars or row vectors of the same size
    structure = model.structure; 
    pulse = model.pulse; 

    switch structure.geometry
        case {'plane','cylinder','sphere'}
            d = parms(1,:); 
            D = parms(2,:); 
            % Check inputs
            if ~isvector(d) || ~isvector(D), error('%s: The input {d, D} should be vectors of the same size', mfilename); end
            if ~isscalar(d) && iscolumn(d), d = d'; end
            if ~isscalar(D) && iscolumn(D), D = D'; end
            if isscalar(d) && ~isscalar(D), d = repmat(d, [1 length(D)]); end
            if isscalar(D) && ~isscalar(d), D = repmat(D, [1 length(d)]); end
            if any(length(d)~=length(D)), error('%s: both d and D should be vectors of the same size', mfilename); end
            
            % Check values
            if any(d<=0), error('%s: The diameter d should be a positive number', mfilename); end
            if any(D<=0), error('%s: D should be a positive number and smaller than the diffusion coefficient of free water', mfilename); end
            
            % Convert to R (radius)
            R = d/2; 
            
            % Load Bk and ak
            [Bk, ak] = load_Bk_ak(R, structure); 
            
        case {'hollowSphere', 'sphericalShell'}
            din = parms(1,:); 
            dout = parms(2,:); 
            D = parms(3,:); 
            
            % Check inputs
            if any(din>=dout), error('din (inner diameter) should be smaller than dout (outer diameter)'); end
            
            % Convert d (diameter) to R (radius)
            Rin = din/2; 
            Rout = dout/2; 
            
            % Load Bk and ak
            [Bk, ak] = load_Bk_ak_hollow(Rin, Rout); 
            
        otherwise
            error('%s: The input model.structure.geometry should be {''plane'',''cylinder'',''sphere'',''hollow_sphere'',''spherical_shell''}', mfilename); 
    end
    
    %% Preliminary setup
    % Convert diameter to radius. All followings use radius
    for n=1:length(pulse.shape)
        if ~ismember(pulse.shape(n), ["pgse"; "tpgse"; "sin"; "cos"; "tcos"; "tsin"] ) 
            error('%s: The #%d shape (%s) is not supported. A shape should be one of [pgse, tpgse, sin, cos, tcos, tsin]', mfilename, n, pulse.shape(n)); 
        end
    end
    
    % Get indices of each gradient shape
    Ind.pgse    = (pulse.shape == 'pgse'); 
    Ind.tpgse   = (pulse.shape == 'tpgse'); 
    Ind.sin     = (pulse.shape == 'sin'); 
    Ind.cos     = (pulse.shape == 'cos'); 
    Ind.tcos1   = (pulse.shape == 'tcos' & pulse.n==1); 
    Ind.tcos2   = (pulse.shape == 'tcos' & pulse.n==2); 
    Ind.tcos1_5 = (pulse.shape == 'tcos' & pulse.n==1.5); 
    Ind.tcos2_5 = (pulse.shape == 'tcos' & pulse.n==2.5); 
    Ind.tcos3   = (pulse.shape == 'tcos' & pulse.n==3); 
    Ind.tsin1   = (pulse.shape == 'tsin' & pulse.n==1); 
    Ind.tsin2   = (pulse.shape == 'tsin' & pulse.n==2); 
    
    %% Calculate structural parameters
    [~, NR, Nk] = size(Bk); 
    Bk = repmat(Bk, [pulse.Nacq, 1, 1]); 
    ak = repmat(ak, [pulse.Nacq, 1, 1]); 
    D  = repmat(D,  [pulse.Nacq, 1, Nk]); 
    
    % Reshape pulse
    delta = pulse.delta; 
    Delta = pulse.Delta;             
    delta = repmat(delta,[1 NR, Nk]); 
    Delta = repmat(Delta,[1 NR, Nk]); 
    
    %% Calculate attenuation
    % First, pgse for all as a place holder
    Etmp = 1./ak.^2./D.^2.*(ak.*D.*delta - 1 + exp(-ak.*D.*delta) + exp(-ak.*D.*Delta)- 1/2*exp(-ak.*D.*(Delta+delta)) - 1/2*exp(-ak.*D.*(Delta-delta))); 
    
    % tpgse
    if any(Ind.tpgse)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tpgse = 1/4*(12*exp(-(ak.*D).*tp) - 24*exp(-(ak.*D).*trise) - 24*exp(-Delta.*(ak.*D)) - 24*(ak.*D).*trise + 12*exp(-(ak.*D).*(Delta-trise)) + 12*exp(-(ak.*D).*(Delta+trise)) - 6*exp(-(ak.*D).*(Delta-tp)) - 6*exp(-(ak.*D).*(Delta+tp)) - 24*exp(-(ak.*D).*(trise+tp)) + 12*exp(-(ak.*D).*(2*trise+tp)) + 8*(ak.*D).^3.*trise.^3 + 12*(ak.*D).^3.*trise.^2.*tp + 12*exp(-(ak.*D).*(Delta-trise-tp)) + 12*exp(-(ak.*D).*(Delta+trise+tp)) - 6*exp(-(ak.*D).*(Delta-2*trise-tp)) - 6*exp(-(ak.*D).*(Delta+2*trise+tp)) + 24)./(3*(ak.*D).^4.*trise.^2);
        Etmp(Ind.tpgse,:,:) = Etmp_tpgse(Ind.tpgse,:,:); 
    end
    
    % sin
    if any(Ind.sin)
        w = pulse.w; w = repmat(w,[1 NR, Nk]); 
        Etmp_sin = w.^2 ./(ak.^2.*D.^2 + w.^2).^2 .* (ak.*D.*(ak.^2.*D.^2 + w.^2).*delta./2./w.^2 +1 - exp(-ak.*D.*delta) - exp(-ak.*D.*Delta)  +0.5.*exp(-ak.*D.*(delta+Delta)) +0.5.*exp(-ak.*D.*(Delta-delta)));  
        Etmp(Ind.sin,:,:) = Etmp_sin(Ind.sin,:,:); 
    end
    
    % cos
    if any(Ind.cos)
        w = pulse.w; w = repmat(w,[1 NR, Nk]); 
        Etmp_cos = ak.^2 .*D.^2 ./(ak.^2.*D.^2 + w.^2).^2 .* ((ak.^2.*D.^2 + w.^2).*delta./2./ak./D -1 + exp(-ak.*D.*delta) + exp(-ak.*D.*Delta)   -0.5.*exp(-ak.*D.*(delta+Delta)) -0.5.*exp(-ak.*D.*(Delta-delta))); 
        Etmp(Ind.cos,:,:) = Etmp_cos(Ind.cos,:,:); 
    end
    
    % tcos and n=1
    if any(Ind.tcos1)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tcos1 = 1/4*(8*exp(-(ak.*D).*tp) - 8*exp(-(ak.*D).*trise) - 8.*exp(-2.*(ak.*D).*trise) - 16.*exp(-(ak.*D).*Delta) - 24.*(ak.*D).*trise + 4.*exp(-(ak.*D).*(Delta-trise)) + 4.*exp(-(ak.*D).*(Delta+trise)) + 4.*exp(-(ak.*D).*(Delta+2.*trise)) + 4.*exp(-(ak.*D).*(Delta-2.*trise)) - 4.*exp(-(ak.*D).*(Delta-tp)) - 4.*exp(-(ak.*D).*(Delta+tp)) - 8.*exp(-(ak.*D).*(trise+tp)) + 4.*exp(-(ak.*D).*(trise+2.*tp)) - 8.*exp(-(ak.*D).*(2.*trise+tp)) + 8.*exp(-(ak.*D).*(3.*trise+tp)) - 8.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 8.*exp(-(ak.*D).*(3.*trise+3.*tp)) + 8.*exp(-(ak.*D).*(4.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(5.*trise+2.*tp)) + 8.*exp(-(ak.*D).*(5.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(5.*trise+4.*tp)) - 8.*exp(-(ak.*D).*(6.*trise+3.*tp)) - 8.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 4.*exp(-(ak.*D).*(7.*trise+4.*tp)) + 12.*(ak.*D).^3.*trise.^3 + 16.*(ak.*D).^3.*trise.^2.*tp + 4.*exp(-(ak.*D).*(Delta-trise-tp)) - 2.*exp(-(ak.*D).*(Delta-trise-2.*tp)) + 4.*exp(-(ak.*D).*(Delta+trise+tp)) + 4.*exp(-(ak.*D).*(Delta-2.*trise-tp)) - 2.*exp(-(ak.*D).*(Delta+trise+2.*tp)) + 4.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 4.*exp(-(ak.*D).*(Delta-3.*trise-tp)) - 4.*exp(-(ak.*D).*(Delta+3.*trise+tp)) + 4.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) + 4.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 4.*exp(-(ak.*D).*(Delta+3.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(Delta-3.*trise-3.*tp)) - 4.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 4.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) - 2.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) - 2.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) - 4.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) - 4.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) - 2.*exp(-(ak.*D).*(Delta+5.*trise+4.*tp)) - 2.*exp(-(ak.*D).*(Delta-5.*trise-4.*tp)) + 4.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) + 4.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) + 4.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 2.*exp(-(ak.*D).*(Delta+7.*trise+4.*tp)) - 2.*exp(-(ak.*D).*(Delta-7.*trise-4.*tp)) + 16)./((ak.*D).^4.*trise.^2); 
        Etmp(Ind.tcos1,:,:) = Etmp_tcos1(Ind.tcos1,:,:);     
    end
    
    % tcos and n=2
    if any(Ind.tcos2)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tcos2 = 1/4*(24*exp(-(ak.*D).*tp) - 24*exp(-(ak.*D).*trise) - 48.*exp(-2.*(ak.*D).*trise) - 72.*exp(-(ak.*D).*Delta) - 120.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 24.*exp(-(ak.*D).*(Delta+2.*trise)) + 24.*exp(-(ak.*D).*(Delta-2.*trise)) - 12.*exp(-(ak.*D).*(Delta-tp)) - 12.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) + 36.*exp(-(ak.*D).*(trise+2.*tp)) - 24.*exp(-(ak.*D).*(2.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) - 72.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(3.*trise+3.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+3.*tp)) + 36.*exp(-(ak.*D).*(5.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(4.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(5.*trise+3.*tp)) - 24.*exp(-(ak.*D).*(6.*trise+3.*tp)) + 48.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(6.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(7.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(8.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(7.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(8.*trise+5.*tp)) + 24.*exp(-(ak.*D).*(9.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(9.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(9.*trise+7.*tp)) + 24.*exp(-(ak.*D).*(10.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(11.*trise+6.*tp)) + 24.*exp(-(ak.*D).*(11.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(11.*trise+8.*tp)) - 24.*exp(-(ak.*D).*(12.*trise+7.*tp)) - 24.*exp(-(ak.*D).*(12.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(13.*trise+8.*tp)) + 76.*(ak.*D).^3.*trise.^3 + 96.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) - 18.*exp(-(ak.*D).*(Delta-trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 12.*exp(-(ak.*D).*(Delta-2.*trise-tp)) - 18.*exp(-(ak.*D).*(Delta+trise+2.*tp)) + 12.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) + 36.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) + 36.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+3.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-3.*trise-3.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) - 18.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) - 18.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+4.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(Delta-4.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) + 12.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) - 24.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) - 24.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+6.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-6.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+7.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-7.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+8.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(Delta-8.*trise-4.*tp)) - 6.*exp(-(ak.*D).*(Delta+7.*trise+6.*tp)) - 6.*exp(-(ak.*D).*(Delta-7.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+8.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-8.*trise-5.*tp)) - 12.*exp(-(ak.*D).*(Delta+9.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-9.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+9.*trise+6.*tp)) + 12.*exp(-(ak.*D).*(Delta-9.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+9.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-9.*trise-7.*tp)) - 12.*exp(-(ak.*D).*(Delta+10.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-10.*trise-7.*tp)) - 6.*exp(-(ak.*D).*(Delta+11.*trise+6.*tp)) - 6.*exp(-(ak.*D).*(Delta-11.*trise-6.*tp)) - 12.*exp(-(ak.*D).*(Delta+11.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-11.*trise-7.*tp)) - 6.*exp(-(ak.*D).*(Delta+11.*trise+8.*tp)) - 6.*exp(-(ak.*D).*(Delta-11.*trise-8.*tp)) + 12.*exp(-(ak.*D).*(Delta+12.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-12.*trise-7.*tp)) + 12.*exp(-(ak.*D).*(Delta+12.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(Delta-12.*trise-8.*tp)) - 6.*exp(-(ak.*D).*(Delta+13.*trise+8.*tp)) - 6.*exp(-(ak.*D).*(Delta-13.*trise-8.*tp)) + 72)./(3.*(ak.*D).^4.*trise.^2); 
        Etmp(Ind.tcos2,:,:) = Etmp_tcos2(Ind.tcos2,:,:); 
    end
    
    % tcos and n=3
    if any(Ind.tcos3)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tcos3 = 1/4*(24.*exp(-(ak.*D).*tp) - 24.*exp(-(ak.*D).*trise) - 72.*exp(-2.*(ak.*D).*trise) - 96.*exp(-Delta.*(ak.*D)) - 168.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 36.*exp(-(ak.*D).*(Delta+2.*trise)) + 36.*exp(-(ak.*D).*(Delta-2.*trise)) - 12.*exp(-(ak.*D).*(Delta-tp)) - 12.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) + 60.*exp(-(ak.*D).*(trise+2.*tp)) - 24.*exp(-(ak.*D).*(2.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) - 120.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(3.*trise+3.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+3.*tp)) + 60.*exp(-(ak.*D).*(5.*trise+2.*tp)) - 48.*exp(-(ak.*D).*(4.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(5.*trise+3.*tp)) - 24.*exp(-(ak.*D).*(6.*trise+3.*tp)) + 96.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(6.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(7.*trise+5.*tp)) - 48.*exp(-(ak.*D).*(8.*trise+4.*tp)) + 36.*exp(-(ak.*D).*(7.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(8.*trise+5.*tp)) + 24.*exp(-(ak.*D).*(9.*trise+5.*tp)) - 72.*exp(-(ak.*D).*(9.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(9.*trise+7.*tp)) + 24.*exp(-(ak.*D).*(10.*trise+7.*tp)) + 36.*exp(-(ak.*D).*(11.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(10.*trise+8.*tp)) + 24.*exp(-(ak.*D).*(11.*trise+7.*tp)) - 24.*exp(-(ak.*D).*(12.*trise+7.*tp)) + 48.*exp(-(ak.*D).*(12.*trise+8.*tp)) + 24.*exp(-(ak.*D).*(12.*trise+9.*tp)) - 24.*exp(-(ak.*D).*(13.*trise+9.*tp)) - 24.*exp(-(ak.*D).*(14.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(13.*trise+10.*tp)) - 24.*exp(-(ak.*D).*(14.*trise+9.*tp)) + 24.*exp(-(ak.*D).*(15.*trise+9.*tp)) - 24.*exp(-(ak.*D).*(15.*trise+10.*tp)) - 24.*exp(-(ak.*D).*(15.*trise+11.*tp)) + 24.*exp(-(ak.*D).*(16.*trise+11.*tp)) + 12.*exp(-(ak.*D).*(17.*trise+10.*tp)) + 24.*exp(-(ak.*D).*(17.*trise+11.*tp)) + 12.*exp(-(ak.*D).*(17.*trise+12.*tp)) - 24.*exp(-(ak.*D).*(18.*trise+11.*tp)) - 24.*exp(-(ak.*D).*(18.*trise+12.*tp)) + 12.*exp(-(ak.*D).*(19.*trise+12.*tp)) + 116.*(ak.*D).^3.*trise.^3 + 144.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) - 30.*exp(-(ak.*D).*(Delta-trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 12.*exp(-(ak.*D).*(Delta-2.*trise-tp)) - 30.*exp(-(ak.*D).*(Delta+trise+2.*tp)) + 12.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) + 60.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) + 60.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+3.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-3.*trise-3.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) - 30.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) - 30.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) + 24.*exp(-(ak.*D).*(Delta+4.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(Delta-4.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) + 12.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) - 48.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) - 48.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+6.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-6.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+7.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-7.*trise-5.*tp)) + 24.*exp(-(ak.*D).*(Delta+8.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(Delta-8.*trise-4.*tp)) - 18.*exp(-(ak.*D).*(Delta+7.*trise+6.*tp)) - 18.*exp(-(ak.*D).*(Delta-7.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+8.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-8.*trise-5.*tp)) - 12.*exp(-(ak.*D).*(Delta+9.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-9.*trise-5.*tp)) + 36.*exp(-(ak.*D).*(Delta+9.*trise+6.*tp)) + 36.*exp(-(ak.*D).*(Delta-9.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+9.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-9.*trise-7.*tp)) - 12.*exp(-(ak.*D).*(Delta+10.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-10.*trise-7.*tp)) - 18.*exp(-(ak.*D).*(Delta+11.*trise+6.*tp)) - 18.*exp(-(ak.*D).*(Delta-11.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+10.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(Delta-10.*trise-8.*tp)) - 12.*exp(-(ak.*D).*(Delta+11.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-11.*trise-7.*tp)) + 12.*exp(-(ak.*D).*(Delta+12.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-12.*trise-7.*tp)) - 24.*exp(-(ak.*D).*(Delta+12.*trise+8.*tp)) - 24.*exp(-(ak.*D).*(Delta-12.*trise-8.*tp)) - 12.*exp(-(ak.*D).*(Delta+12.*trise+9.*tp)) - 12.*exp(-(ak.*D).*(Delta-12.*trise-9.*tp)) + 12.*exp(-(ak.*D).*(Delta+13.*trise+9.*tp)) + 12.*exp(-(ak.*D).*(Delta-13.*trise-9.*tp)) + 12.*exp(-(ak.*D).*(Delta+14.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(Delta-14.*trise-8.*tp)) - 6.*exp(-(ak.*D).*(Delta+13.*trise+10.*tp)) - 6.*exp(-(ak.*D).*(Delta-13.*trise-10.*tp)) + 12.*exp(-(ak.*D).*(Delta+14.*trise+9.*tp)) + 12.*exp(-(ak.*D).*(Delta-14.*trise-9.*tp)) - 12.*exp(-(ak.*D).*(Delta+15.*trise+9.*tp)) - 12.*exp(-(ak.*D).*(Delta-15.*trise-9.*tp)) + 12.*exp(-(ak.*D).*(Delta+15.*trise+10.*tp)) + 12.*exp(-(ak.*D).*(Delta-15.*trise-10.*tp)) + 12.*exp(-(ak.*D).*(Delta+15.*trise+11.*tp)) + 12.*exp(-(ak.*D).*(Delta-15.*trise-11.*tp)) - 12.*exp(-(ak.*D).*(Delta+16.*trise+11.*tp)) - 12.*exp(-(ak.*D).*(Delta-16.*trise-11.*tp)) - 6.*exp(-(ak.*D).*(Delta+17.*trise+10.*tp)) - 6.*exp(-(ak.*D).*(Delta-17.*trise-10.*tp)) - 12.*exp(-(ak.*D).*(Delta+17.*trise+11.*tp)) - 12.*exp(-(ak.*D).*(Delta-17.*trise-11.*tp)) - 6.*exp(-(ak.*D).*(Delta+17.*trise+12.*tp)) - 6.*exp(-(ak.*D).*(Delta-17.*trise-12.*tp)) + 12.*exp(-(ak.*D).*(Delta+18.*trise+11.*tp)) + 12.*exp(-(ak.*D).*(Delta-18.*trise-11.*tp)) + 12.*exp(-(ak.*D).*(Delta+18.*trise+12.*tp)) + 12.*exp(-(ak.*D).*(Delta-18.*trise-12.*tp)) - 6.*exp(-(ak.*D).*(Delta+19.*trise+12.*tp)) - 6.*exp(-(ak.*D).*(Delta-19.*trise-12.*tp)) + 96)./(3.*(ak.*D).^4.*trise.^2); 
        Etmp(Ind.tcos3,:,:) = Etmp_tcos3(Ind.tcos3,:,:); 
    end
    
    % tcos and n=1.5
    if any(Ind.tcos1_5) 
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tcos1_5_test = calculate_Etmp_1_5(D, ak, tp, trise, Delta);
        Etmp(Ind.tcos1_5,:,:) = Etmp_tcos1_5_test(Ind.tcos1_5,:,:);
    end
    
    % tcos and n=2.5
    if any(Ind.tcos2_5) 
        Etmp(Ind.tcos2_5,:,:) = integrate_tcos_2_5( ...
            ak(Ind.tcos2_5,:,:).*D(Ind.tcos2_5,:,:), ...
            pulse.trise(Ind.tcos2_5), pulse.tp(Ind.tcos2_5), ...
            pulse.delta(Ind.tcos2_5), pulse.Delta(Ind.tcos2_5));
    end
    
    % tsin and n=1
    if any(Ind.tsin1)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tsin1 = 1/4*(24.*exp(-(ak.*D).*tp) - 24.*exp(-(ak.*D).*trise) - 12.*exp(-2.*(ak.*D).*trise) - 36.*exp(-Delta.*(ak.*D)) - 48.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 6.*exp(-(ak.*D).*(Delta+2.*trise)) + 6.*exp(-(ak.*D).*(Delta-2.*trise)) - 12.*exp(-(ak.*D).*(Delta-tp)) - 12.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) - 24.*exp(-(ak.*D).*(2.*trise+tp)) - 12.*exp(-(ak.*D).*(2.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(4.*trise+2.*tp)) + 16.*(ak.*D).^3.*trise.^3 + 24.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 12.*exp(-(ak.*D).*(Delta-2.*trise-tp)) + 12.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) + 6.*exp(-(ak.*D).*(Delta+2.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-2.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 6.*exp(-(ak.*D).*(Delta+4.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-4.*trise-2.*tp)) + 36) ./ (3.*(ak.*D).^4.*trise.^2); 
        Etmp(Ind.tsin1,:,:) = Etmp_tsin1(Ind.tsin1,:,:); 
    end
    
    % tsin and n=2
    if any(Ind.tsin2)
        trise = pulse.trise; tp = pulse.tp; trise = repmat(trise,[1 NR, Nk]); tp = repmat(tp,[1 NR, Nk]); 
        Etmp_tsin2 = 1/4*(48.*exp(-(ak.*D).*tp) - 24.*exp(-(ak.*D).*trise) - 36.*exp(-2.*(ak.*D).*trise) - 60.*exp(-Delta.*(ak.*D)) - 96.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 18.*exp(-(ak.*D).*(Delta+2.*trise)) + 18.*exp(-(ak.*D).*(Delta-2.*trise)) - 24.*exp(-(ak.*D).*(Delta-tp)) - 24.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) - 72.*exp(-(ak.*D).*(2.*trise+tp)) - 36.*exp(-(ak.*D).*(2.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+tp)) + 48.*exp(-(ak.*D).*(4.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+3.*tp)) - 24.*exp(-(ak.*D).*(5.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(5.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(6.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(6.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(7.*trise+3.*tp)) + 24.*exp(-(ak.*D).*(7.*trise+4.*tp)) - 12.*exp(-(ak.*D).*(8.*trise+4.*tp)) + 32.*(ak.*D).^3.*trise.^3 + 48.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 36.*exp(-(ak.*D).*(Delta-2.*trise-tp)) + 36.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) + 18.*exp(-(ak.*D).*(Delta+2.*trise+2.*tp)) + 18.*exp(-(ak.*D).*(Delta-2.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+tp)) - 24.*exp(-(ak.*D).*(Delta+4.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(Delta-4.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) + 12.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) + 12.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) + 6.*exp(-(ak.*D).*(Delta+6.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-6.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) + 6.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) + 6.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+7.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-7.*trise-3.*tp)) - 12.*exp(-(ak.*D).*(Delta+7.*trise+4.*tp)) - 12.*exp(-(ak.*D).*(Delta-7.*trise-4.*tp)) + 6.*exp(-(ak.*D).*(Delta+8.*trise+4.*tp)) + 6.*exp(-(ak.*D).*(Delta-8.*trise-4.*tp)) + 60) ./ (3.*(ak.*D).^4.*trise.^2);  
        Etmp(Ind.tsin2,:,:) = Etmp_tsin2(Ind.tsin2,:,:); 
    end
        
    %% Final signal attenuation
    G = repmat(pulse.G, [1 NR]); 
    E = exp( - 2 .* (pulse.gamma*G).^2 .* sum(Bk .* Etmp,3)); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bk, ak] = load_Bk_ak(R, structure) 
%% Get Bk and ak 
    NR = length(R); 
    
    switch structure.geometry
        case 'plane'
            % R is the half distance between two planes
            krange = 1:30;             
            Nk = length(krange); 
            krange = reshape(krange,[1 1 Nk]); 
            krange = repmat(krange, [1 NR 1]); 
            R = repmat(R, [1 1 Nk]); 
            Bk = 8.*(2*R).^2./pi^4./(2.*krange-1).^4; 
            ak = pi^2.*(2.*krange-1).^2./(2*R).^2; 
        case 'cylinder'
            muk = [ 1.8412, 5.3315, 8.5364, 11.7061, 14.8636, 18.0156, 21.1644, 24.3114, 27.4571, 30.6020, 33.7462, 36.8900, 40.0335, 43.1767, 46.3196, 49.4624, 52.6051, 55.7476, 58.8901, 62.0324, 65.1747, 68.3169, 71.4590, 74.6011, 77.7432, 80.8852, 84.0272, 87.1692, 90.3111  ]'; 
            Nk = length(muk); 
            muk = reshape(muk, [1 1 Nk]); muk = repmat(muk,[1 NR 1]);             R = repmat(R,[1 1 Nk]); 
           Bk = 2.*(R./muk).^2./(muk.^2-1); 
           ak = (muk./R).^2; 
        case 'sphere'
            muk = [ 2.082,  5.941, 9.206, 12.405, 15.580, 18.743, 21.900, 25.053, 28.204, 31.353, 34.500, 37.646, 40.792, 43.937, 47.082, 50.226, 53.370, 56.514, 59.657, 62.801, 65.944, 69.087, 72.229, 75.372, 78.515, 81.657, 84.800, 87.942, 91.085  ]'; 
            Nk = length(muk); 
            muk = reshape(muk, [1 1 Nk]); muk = repmat(muk,[1 NR 1]);             R = repmat(R,[1 1 Nk]); 
            Bk = 2.*(R./muk).^2./(muk.^2-2); 
            ak = (muk./R).^2; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bk, ak] = load_Bk_ak_hollow(Rin, Rout) 
%% Get Bk and ak for sphericalShell or hollowSphere
    lambdas = get_hollow_sphere_roots(Rin, Rout); 
    NumLambda = size(lambdas,1); 
    ak = lambdas.^2; 

    Rin = repmat(Rin, [NumLambda, 1]); 
    Rout = repmat(Rout, [NumLambda, 1]); 
    ratiotemp = Rout.^3./Rin.^3 .* sphbesselj1_core_prim(lambdas.*Rin) ./ sphbesselj1_core_prim(lambdas.*Rout); 
    Bk = 2*Rin.^3.*Rout.^3./(Rout.^3-Rin.^3)./lambdas.^2 .* (ratiotemp - 1).^2 ...
            ./ (Rin.^3.*(lambdas.^2.*Rout.^2-2).*ratiotemp.^2 - Rout.^3.*(lambdas.^2.*Rin.^2-2)); 
    % Permute dimension
    ak = permute(ak, [3 2 1]); 
    Bk = permute(Bk, [3 2 1]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lambdan = get_hollow_sphere_roots(a,b)
%% Get lambda for hollow sphere structure
    xstep = 0.00001;     
    xmax = 100; 
    x = xstep:xstep:xmax; 
    x = x'; 
    numRoots = 20; 
    lambdan = zeros([numRoots, size(a,2)]); 
    
    y1 = sphbessely1_core_prim(x.*a).*sphbesselj1_core_prim(x.*b) - sphbesselj1_core_prim(x.*a).*sphbessely1_core_prim(x.*b); 
    
    x = x + xstep; 
    y2 = sphbessely1_core_prim(x.*a).*sphbesselj1_core_prim(x.*b) - sphbesselj1_core_prim(x.*a).*sphbessely1_core_prim(x.*b); 
    
    for n = 1:size(a,2)
        root_tmp = x(y1(:,n).*y2(:,n)<0); 
        lambdan(:,n) = root_tmp(1:numRoots); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spherical bessel functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sphbesselj1_core_prim(x)
    y = -2.*sin(x) + 2.*x.*cos(x) + x.^2.*sin(x); 
end
function y = sphbessely1_core_prim(x)
    y = 2.*cos(x) + 2.*x.*sin(x) - x.^2.*cos(x); 
end
function y = sphbesselj(n,x)
    y = besselj(n+0.5,x)./sqrt(x); 
end
function y = sphbessely(n,x)
    y = bessely(n+0.5,x)./sqrt(x); 
end
function y = sphbesselj_prim(n,x)
    y = n./x.*sphbesselj(n,x) - sphbesselj(n+1,x); 
end
function y = sphbessely_prim(n,x)
    y = n./x.*sphbessely(n,x) - sphbessely(n+1,x); 
end

function Etmp_tcos_test = calculate_Etmp_1_5(D, ak, tp, trise, Delta)
    % Initialize A and B
    A_sum = zeros(size(D));
    B_sum = zeros(size(D));
    
    % Common exponent base (for exponent adjustment)
    base_exp_A = D .* ak .* (6*tp + 10*trise);
    base_exp_B = D .* ak .* (Delta + 6*tp + 10*trise);
    
    % Define terms for A (27 terms)
    A_terms = {
        % Coefficient function                 Exponential function
        @(D,ak,tp,trise) 18*D.^3.*ak.^3.*tp.*trise.^2, @(tp,trise) 2*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) 14*D.^3.*ak.^3.*trise.^3,     @(tp,trise) 2*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -24*D.*ak.*trise,             @(tp,trise) 2*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -3,                           @(tp,trise) 2*D.*ak.*trise;
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*trise;
        @(D,ak,tp,trise) -3,                           @(tp,trise) 2*D.*ak.*(tp+trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(tp+trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) 5*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) 3*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) 2*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) -3,                           @(tp,trise) 2*D.*ak.*(tp+3*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(tp+3*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(tp+4*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) 3*D.*ak.*(2*tp+3*trise);
        @(D,ak,tp,trise) -9,                           @(tp,trise) 2*D.*ak.*(3*tp+4*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(3*tp+4*trise);
        @(D,ak,tp,trise) 15,                           @(tp,trise) 2*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(3*tp+7*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(4*tp+5*trise);
        @(D,ak,tp,trise) -12,                          @(tp,trise) D.*ak.*(4*tp+7*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(4*tp+9*trise);
        @(D,ak,tp,trise) 6,                            @(tp,trise) D.*ak.*(5*tp+7*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(5*tp+8*trise);
        @(D,ak,tp,trise) -6,                           @(tp,trise) D.*ak.*(5*tp+9*trise);
        @(D,ak,tp,trise) -3,                           @(tp,trise) 0;  % Constant term
    };
    % Calculate each term of A (adjust exponent)
    for i = 1:size(A_terms, 1)
        coef_fun = A_terms{i,1};
        exp_fun  = A_terms{i,2};
        
        coef_val = coef_fun(D, ak, tp, trise);
        exp_val  = exp_fun(tp, trise);
        
        % Exponent adjustment: subtract base exponent
        adjusted_exp = exp_val - base_exp_A;
        % Numerical stability processing
        adjusted_exp = min(adjusted_exp, 700);
        term_val = coef_val .* exp(adjusted_exp);
        
        A_sum = A_sum + term_val;
    end
    
    % Define terms for B (52 terms)
    B_terms = {
        % Coefficient function                 Exponential function
        @(D,ak,tp,trise) 1,   @(tp,trise) 2*D.*ak.*trise;
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*trise;
        @(D,ak,tp,trise) 1,   @(tp,trise) 2*D.*ak.*(tp+trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(tp+trise);
        @(D,ak,tp,trise) 3,   @(tp,trise) 6*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) 5*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) 3*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) 2*D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(tp+2*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 2*D.*ak.*(tp+3*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(tp+3*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(tp+4*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 6*D.*ak.*(2*tp+3*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) 3*D.*ak.*(2*tp+3*trise);
        @(D,ak,tp,trise) 3,   @(tp,trise) 2*D.*ak.*(3*tp+4*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(3*tp+4*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 4*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) 3*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -10, @(tp,trise) 2*D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(3*tp+5*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(3*tp+7*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(4*tp+5*trise);
        @(D,ak,tp,trise) 4,   @(tp,trise) D.*ak.*(4*tp+7*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(4*tp+9*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 2*D.*ak.*(5*tp+7*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(5*tp+7*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) 2*D.*ak.*(5*tp+8*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(5*tp+8*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 2*D.*ak.*(5*tp+9*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(5*tp+9*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(6*tp+11*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(7*tp+10*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(7*tp+11*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(7*tp+12*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(7*tp+13*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(8*tp+11*trise);
        @(D,ak,tp,trise) 4,   @(tp,trise) D.*ak.*(8*tp+13*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(8*tp+15*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(9*tp+13*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(9*tp+14*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(9*tp+16*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(11*tp+16*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(11*tp+17*trise);
        @(D,ak,tp,trise) 2,   @(tp,trise) D.*ak.*(11*tp+18*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(11*tp+19*trise);
        @(D,ak,tp,trise) -2,  @(tp,trise) D.*ak.*(12*tp+19*trise);
        @(D,ak,tp,trise) 1,   @(tp,trise) 0;  % Constant term
    };
    
    % Calculate each term of B (adjust exponent)
    for i = 1:size(B_terms, 1)
        coef_fun = B_terms{i,1};
        exp_fun  = B_terms{i,2};
        
        coef_val = coef_fun(D, ak, tp, trise);
        exp_val  = exp_fun(tp, trise);
        
        % Exponent adjustment: subtract base exponent
        adjusted_exp = exp_val - base_exp_B;
        % Numerical stability processing
        adjusted_exp = min(adjusted_exp, 700);
        term_val = coef_val .* exp(adjusted_exp);
        
        B_sum = B_sum + term_val;
    end
    
    % Calculate denominator
    denominator = 6 * (ak .* D).^4 .* trise.^2;
    
    % Final result
    Etmp_tcos_test = (2 * A_sum + 3 * B_sum) ./ denominator;
end

function Etmp = integrate_tcos_2_5(lambda, trise, tp, delta, Delta)
% Integrate g(t)g(s)exp(-lambda*(t-s)) exactly on each linear segment.
    history = zeros(size(lambda));
    integral = zeros(size(lambda));
    for lobe = [1, -1]
        [history, integral] = linear_segment(history, integral, lambda, trise, 0, lobe);
        [history, integral] = linear_segment(history, integral, lambda, tp, lobe, lobe);
        polarity = lobe;
        for reversal = 1:5
            [history, integral] = linear_segment(history, integral, lambda, ...
                2*trise, polarity, -polarity);
            polarity = -polarity;
            if reversal == 5
                plateau = tp;
            else
                plateau = trise + 2*tp;
            end
            [history, integral] = linear_segment(history, integral, lambda, ...
                plateau, polarity, polarity);
        end
        [history, integral] = linear_segment(history, integral, lambda, trise, polarity, 0);
        history = history .* exp(-lambda .* reshape(Delta-delta, [], 1, 1));
    end
    % The signal exponent uses 2*Etmp; the triangular double integral is integral.
    Etmp = integral/2;
end

function [history, integral] = linear_segment(history, integral, lambda, width, from, to)
    width = reshape(width, [], 1, 1);
    from = reshape(from, [], 1, 1);
    to = reshape(to, [], 1, 1);
    change = to-from;
    z = lambda.*width;
    [m0, m1, j0, j2] = kernel_moments(z);
    integral = integral + history.*width.*(from.*m0 + change.*m1) ...
        + width.^2.*((from.^2 + from.*change).*j0 + change.^2.*j2);
    history = exp(-z).*history + width.*(to.*m0 - change.*m1);
end

function [m0, m1, j0, j2] = kernel_moments(z)
    m0 = zeros(size(z)); m1 = m0; j0 = m0; j2 = m0;
    small = z < 0.2;
    zs = z(small);
    term = ones(size(zs));
    for order = 0:14
        if order > 0
            term = term.*(-zs)/order;
        end
        m0(small) = m0(small) + term/(order+1);
        m1(small) = m1(small) + term/(order+2);
        j0(small) = j0(small) + term/((order+1)*(order+2));
        j2(small) = j2(small) + term.* ...
            (2/(order+1)-3/(order+2)+1/(order+4))/6;
    end
    zl = z(~small);
    ez = exp(-zl);
    m0(~small) = -expm1(-zl)./zl;
    m1(~small) = (1-(1+zl).*ez)./zl.^2;
    j0(~small) = (zl+expm1(-zl))./zl.^2;
    m3 = 6*(1-ez.*(1+zl+zl.^2/2+zl.^3/6))./zl.^4;
    j2(~small) = (2*m0(~small)-3*m1(~small)+m3)/6;
end
