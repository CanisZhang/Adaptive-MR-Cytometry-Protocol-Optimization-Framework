classdef IMPULSED < mati.SignalModel
%% class definition of IMPULSED  
% IMPULSED can take different structure models
% Model choices: 
%           'impulsed_vin_d_Dex':   two-compartment with vin (intracellular proton fraction), d (intracelular size), and extracelluar diffusion coefficient Dex
%           'impulsed_vin_d_Dex_Din':    two-compartment with vin, d, Dex, Din
% -------------------------------------------------------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------------------------------------------------------

    %% properties
    properties (SetAccess = private)
        structure
        defaultStructure
    end
    
    %% static methods
    methods (Static)
        out = PostProcessIMPULSED(fitout) ; 
        out = Signal_IMPULSED_vin_d_Dex(parms, this) ; 
        out = Signal_IMPULSED_vin_d_Dex_Din(parms, this) ; 
    end

    %% other methods
    methods                
        % ------------------------------- Constructor -------------------------------------
        % IMPULSED should take two inputs: a struct (structure model) and an object (DiffusionPulseSequence)
        function this = IMPULSED(varargin)
            if nargin == 0 , return ; end      % return an empty object
            for n = 1:length(varargin)
                if isa(varargin{n},'mati.DiffusionPulseSequence')
                    this.pulse = varargin{n};
                elseif isstruct(varargin{n})
                    inputStructure = varargin{n};
                else
                    error('Inputs for IMPULSED should be StructureModel, DiffusionPulseSequence')
                end
            end    
            
            %------------------------ setting signal function handle and default fitopts --------------------------
            % !!! Note that the DefaultFitopts.x0, lb, ub should all be column vectors !!!
            switch inputStructure.modelName
                case 'IMPULSED_vin_d_Dex'                                                                                      % [vin, d, Dex]
                    this.parmsLabels = ["vin", "d", "Dex"] ; this.NumParms = 3 ; 
                    this.FcnSignal = @mati.IMPULSED.Signal_IMPULSED_vin_d_Dex; 
                    this.defaultFitopts.x0 = [0.5, 15, 1.5]' ;  this.defaultFitopts.lb = [0, 0, 0]' ; this.defaultFitopts.ub = [1, 50, 3]' ; 
                    this.defaultStructure.modelName = 'IMPULSED_vin_d_Dex'  ; this.defaultStructure.Din = 1.56 ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere';
                    this.defaultFitopts.dictionary.vins = [0.:0.02:1] ; this.defaultFitopts.dictionary.ds = [0.2:0.2:30] ; this.defaultFitopts.dictionary.Dexs = [0.2:0.2:3] ; 
                    this.defaultFitopts.dictionary.parms = {this.defaultFitopts.dictionary.vins, this.defaultFitopts.dictionary.ds, this.defaultFitopts.dictionary.Dexs} ; 
                    this.defaultFitopts.priorDist.distName = 'gaussian' ; 
                    this.defaultFitopts.priorDist.pd.vins = normpdf(this.defaultFitopts.dictionary.vins,0.7,10) ; this.defaultFitopts.priorDist.pd.ds = normpdf(this.defaultFitopts.dictionary.ds,14,200) ; this.defaultFitopts.priorDist.pd.Dexs = normpdf(this.defaultFitopts.dictionary.Dexs,1.8, 20) ; 
                case 'IMPULSED_vin_d_Dex_Din'                                                                               % [vin, d, Dex, Din]
                    this.parmsLabels = ["vin", "d", "Dex", "Din"] ; this.NumParms = 4 ; 
                    this.FcnSignal = @mati.IMPULSED.Signal_IMPULSED_vin_d_Dex_Din; 
                    this.defaultFitopts.x0 = [0.5, 10, 1.5, 1.5]' ; this.defaultFitopts.lb = [0, 0, 0, 0]' ; this.defaultFitopts.ub = [1, 30, 3, 3]' ; 
                    this.defaultStructure.modelName = 'IMPULSED_vin_d_Dex_Din'  ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere' ; 
                    this.defaultFitopts.dictionary.vins = [0.02:0.02:1] ; this.defaultFitopts.dictionary.ds = [0.2:0.2:30] ; this.defaultFitopts.dictionary.Dexs = [0.25:0.25:3] ; this.defaultFitopts.dictionary.Dins = [0.5:.25:3] ; 
                    this.defaultFitopts.dictionary.parms = {this.defaultFitopts.dictionary.vins, this.defaultFitopts.dictionary.ds, this.defaultFitopts.dictionary.Dexs, this.defaultFitopts.dictionary.Dins} ; 
                    this.defaultFitopts.priorDist.distName = 'gaussian' ; 
                    this.defaultFitopts.priorDist.pd.vins = normpdf(this.defaultFitopts.dictionary.vins,0.7,14) ; this.defaultFitopts.priorDist.pd.ds = normpdf(this.defaultFitopts.dictionary.ds,14,200) ; this.defaultFitopts.priorDist.pd.Dexs = normpdf(this.defaultFitopts.dictionary.Dexs,1.8,40) ; this.defaultFitopts.priorDist.pd.Dins = normpdf(this.defaultFitopts.dictionary.Dins,1.56,32) ; 
                otherwise
                    error('%s: This subset supports IMPULSED_vin_d_Dex and IMPULSED_vin_d_Dex_Din.',mfilename) ;
            end
            
            % ------------------ update structure -----------------------
            this.structure = this.defaultStructure ; 
            names = fieldnames(inputStructure) ; 
            for n=1:length(names)
                if ~isfield(this.structure, names{n}) , error('%s: The input structure has a field that is not a field of mig.IMPULSED.structure',mfilename) ; end
                this.structure.(names{n}) = inputStructure.(names{n}) ; 
            end

            % ------------------- set default postProcess function -------------------
            this.FcnPostProcess = @mati.IMPULSED.PostProcessIMPULSED ; 

            % -------------------- set other defaultFitopts fields --------------------------------
            this.defaultFitopts.fittingMethod = 'nllsq' ;
            this.defaultFitopts.solverName = 'lsqnonlin' ;
            this.defaultFitopts.noiseModel = 'none' ;       % {'none','simple','rician','meanrician'}
            this.defaultFitopts.options =optimset('Display','off') ; 
            this.defaultFitopts.S0_bmax = 0.020 ;                  % maximum b for S(b=0)
            this.defaultFitopts.ADCfit_bmin = 0.19 ;                % minimum b for ADC fitting
            this.defaultFitopts.ADCfit_bmax = 1.05 ;             % maximum b for ADC fitting
            this.defaultFitopts.NumStarts = 5 ; 
            % fmincon 
            this.defaultFitopts.Aeq = [] ; 
            this.defaultFitopts.beq = [] ; 
            this.defaultFitopts.Aineq = [] ; 
            this.defaultFitopts.bineq = [] ; 
            this.defaultFitopts.nonlcon = [] ; 
            % flags for many options
            this.defaultFitopts.flag.useGPU = 'n' ;
            this.defaultFitopts.flag.indGPU = 1 ; 
            this.defaultFitopts.flag.denoise = 'n' ;
            this.defaultFitopts.flag.degibbs = 'n' ; 
            this.defaultFitopts.flag.deivim = 'n' ;
            this.defaultFitopts.flag.getADC = 'n' ; 
            this.defaultFitopts.flag.postProcess = 'y' ; 
            this.defaultFitopts.flag.multistart = 'n' ; 
            this.defaultFitopts.flag.parfor = 'y' ; 
            
            if ~iscolumn(this.defaultFitopts.x0) && ~iscolumn(this.defaultFitopts.lb) && ~iscolumn(this.defaultFitopts.ub)
                error('%s: the DefaultFitopts.x0, lb, ub should all be column vectors',mfilename) ; 
            end
            
        end % End of Constructor
        
    end % End of Methods
    
end % End of class
