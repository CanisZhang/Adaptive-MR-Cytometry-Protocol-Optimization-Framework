function [out,this] = Fit(this, data)
% Fit IMPULSED signals using bounded nonlinear least squares and MultiStart.
% exitflag is true for failed fits. solverExitflag retains the solver code.
if ~isa(this,'mati.FitPars') || ~isa(data,'mati.ImageData')
    error('%s: Inputs must be FitPars and ImageData objects.',mfilename);
end
if ~strcmp(this.fitopts.fittingMethod,'nllsq') || ...
        ~strcmp(this.fitopts.solverName,'lsqnonlin')
    error('%s: This subset supports nllsq with lsqnonlin.',mfilename);
end
for name = {'denoise','degibbs','deivim','getADC'}
    if any(strcmpi(this.fitopts.flag.(name{1}),{'y','yes'}))
        error('%s: %s is unavailable in this subset.',mfilename,name{1});
    end
end

data = mati.dMRI.PreProcessDWI(this,data);
out.data = data;
out.mapLabels = this.model.parmsLabels;
map0 = zeros([data.Nx,data.Ny,data.Nz]);
for n = 1:length(out.mapLabels)
    out.(out.mapLabels(n)) = map0;
end
out.resnorm = map0;
out.exitflag = false([data.Nx,data.Ny,data.Nz]);
out.solverExitflag = nan([data.Nx,data.Ny,data.Nz]);

if any(strcmpi(this.fitopts.flag.parfor,{'y','yes'}))
    poolobj = gcp;
    Ncore = poolobj.NumWorkers;
    parData = parallel.pool.Constant(data);
    parFitPars = parallel.pool.Constant(this);
else
    Ncore = 0;
    parData.Value = data;
    parFitPars.Value = this;
end
if any(strcmpi(this.fitopts.flag.multistart,{'y','yes'}))
    NumStarts = this.fitopts.NumStarts;
else
    NumStarts = 1;
end
[row,col,slice] = ind2sub(size(data.mask),find(data.mask>0));
Npixel = length(row);
parmap1d = cell(Npixel,1);
resnormmap1d = zeros(Npixel,1);
exitflag1d = zeros(Npixel,1);
solverExitflag1d = nan(Npixel,1);
parfor (npixel = 1:Npixel,Ncore)
    nx = row(npixel);
    ny = col(npixel);
    nz = slice(npixel);
    signal = squeeze(parData.Value.img(nx,ny,nz,:));
    if isempty(parData.Value.sigma) || isscalar(parData.Value.sigma)
        sigma = parData.Value.sigma;
    else
        sigma = parData.Value.sigma(nx,ny,nz);
    end
    try
        [x,resnorm,solverExitflag] = run_optimization( ...
            parFitPars.Value,signal,sigma,NumStarts);
        solverExitflag1d(npixel) = solverExitflag;
        exitflag1d(npixel) = solverExitflag<=0 || ...
            any(~isfinite(x(:))) || ~isfinite(resnorm);
    catch
        warning('Fitting in voxel %d (%d,%d,%d) failed.',npixel,nx,ny,nz);
        x = nan(parFitPars.Value.model.NumParms,1);
        resnorm = NaN;
        exitflag1d(npixel) = true;
    end
    parmap1d{npixel} = x;
    resnormmap1d(npixel) = resnorm;
end
for npixel = 1:Npixel
    for n = 1:length(out.mapLabels)
        out.(out.mapLabels(n))(row(npixel),col(npixel),slice(npixel)) = ...
            parmap1d{npixel}(n);
    end
    out.resnorm(row(npixel),col(npixel),slice(npixel)) = resnormmap1d(npixel);
    out.exitflag(row(npixel),col(npixel),slice(npixel)) = exitflag1d(npixel);
    out.solverExitflag(row(npixel),col(npixel),slice(npixel)) = ...
        solverExitflag1d(npixel);
end
if any(strcmpi(this.fitopts.flag.postProcess,{'y','yes'}))
    out = this.model.FcnPostProcess(out);
end
for n = 1:length(out.mapLabels)
    out.(out.mapLabels(n)) = single(out.(out.mapLabels(n)));
end
end

function [x,resnorm,solverExitflag] = run_optimization(this,signal,sigma,NumStarts)
ObjectiveFcn = @(parms)mati.FitPars.ObjectiveFcn(parms,this,signal,sigma);
problem = createOptimProblem('lsqnonlin','x0',this.fitopts.x0, ...
    'objective',ObjectiveFcn,'lb',this.fitopts.lb,'ub',this.fitopts.ub, ...
    'options',this.fitopts.options);
ms = MultiStart('UseParallel',false,'Display','off');
[x,resnorm,solverExitflag] = run(ms,problem,NumStarts);
if isempty(x)
    x = nan(this.model.NumParms,1);
    resnorm = NaN;
end
end
