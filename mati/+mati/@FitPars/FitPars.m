classdef FitPars
%% class defination of FitPars to save fitting parameters
%
% -------------------------------------------------------------------------------------------------------------------------

    properties
        model
        fitopts
    end
    methods (Static)
        out = ObjectiveFcn(parms, this, ydata, sigma) ; 
    end
    
    
    methods        
        function this = set.fitopts(this,fitopts)
            this.fitopts = fitopts;
        end
        function this = set.model(this,model)
            this.model = model;
        end
        
        function this = FitPars(varargin)
            switch nargin
                case 0
                    return 
                case 1
                    if ~isa(varargin{1}, 'mati.SignalModel'), error('%s: FitPars should take a SignalModel object', mfilename) ; end
                    this.model = varargin{1} ; 
                    this.fitopts = this.model.defaultFitopts ; 
                case 2
                    for n = 1:nargin        % fitopts structure and SignalModel object can be in any order
                        if isa(varargin{n},'mati.SignalModel')
                            this.model = varargin{n} ; 
                        elseif isstruct(varargin{n})
                            inputFitopts = varargin{n} ;              % take input fitopts structure
                        else
                            error('FitPars can only take a fitopts structure and SignalModel object inputs')
                        end
                    end
                    % use input fitopts to update defaultFitopts fields
                    this.fitopts = this.model.defaultFitopts ; 
                    this.fitopts = mati.Tools.UpdateStruct(this.fitopts, inputFitopts) ; 
                    % update dictionary.parms if it is not a special fitting such as MRI-Cytometry 
                    if ~strcmp(this.fitopts.fittingMethod,'special')
                        for nparm=1:length(this.model.parmsLabels)
                            this.fitopts.dictionary.parms{nparm}=this.fitopts.dictionary.([char(this.model.parmsLabels(nparm)),'s']);
                        end
                    end
                otherwise
                    error('%s: Input should be either a SignalModel object,either with or without a fitopts structure', mfilename) ; 
            end % End of switch nargin
            
        end % End of constructor
        
    end % End of Methods
    
end % End of class
