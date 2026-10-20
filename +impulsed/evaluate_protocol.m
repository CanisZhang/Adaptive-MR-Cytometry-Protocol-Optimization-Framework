function out = evaluate_protocol(x, truth, n, seed, cfg, workers, model, protocol)
% Fixed seeds pair candidate noise and fit starts, with independent validation.
if nargin < 6, workers = 0; end
if nargin < 8
    [model, protocol] = impulsed.protocol_model(x,cfg);
end
if ~protocol.feasible
    error('optimization:InfeasibleProtocol','The requested protocol is infeasible.');
end
assert(size(truth,2) == 4 && (size(truth,3) == 1 || size(truth,3) == n), ...
    'Truth must be scenarios-by-4 or scenarios-by-4-by-realizations.');
fitpars = mati.FitPars(model,cfg.fitopts);
nt = size(truth,1);
np = model.NumParms;
estimates = nan(n,nt,np);
errors = nan(n,nt,3);
loss = nan(n,nt);
failed = false(n,nt);
residual = nan(n,nt);
exitflag = nan(n,nt);
state = rng;
warning_state = warning;
cleanup = onCleanup(@()restore_state(state,warning_state)); %#ok<NASGU>
parfor (r = 1:n,workers)
    target = truth(:,:,min(r,size(truth,3)));
    [estimates(r,:,:),errors(r,:,:),loss(r,:),failed(r,:), ...
        residual(r,:),exitflag(r,:)] = simulate_one(model,fitpars,target,protocol,cfg,seed+r);
end
out.protocol = protocol;
out.estimates = estimates;
out.relative_error = errors;
out.failed = failed;
out.normalized_residual = residual;
out.solver_exitflag = exitflag;
out.loss = loss;
out.objective = mean(loss(:));
out.error_mean = mean(loss(:)-cfg.failure_penalty*double(failed(:)));
out.error_std = std(loss(:)-cfg.failure_penalty*double(failed(:)));
out.failure_rate = mean(failed(:));
out.seed = seed;
out.repeats = n;
end

function [estimate,error,loss,failed,residual,exitflag] = simulate_one(model,fitpars,target,protocol,cfg,seed)
warning_state = warning;
cleanup = onCleanup(@()warning(warning_state)); %#ok<NASGU>
signal = model.FcnSignal(target(:,1:model.NumParms)',model);
S0 = exp(-(protocol.TE-cfg.TE0)/cfg.T2);
sigma = 1/(cfg.SNR0*sqrt(cfg.num_averages));
rng(seed,'twister');
noisy = mati.Physics.AddRicianNoise(S0*signal,sigma);
nt = size(target,1);
data = mati.ImageData(reshape(noisy',[nt 1 1 size(noisy,1)]),[]);
rng(seed+cfg.fit_seed_offset,'twister');
fit = fitpars.Fit(data);
estimate = [double(fit.vin(:)),double(fit.d(:)),double(fit.Dex(:))];
if model.NumParms == 4, estimate(:,4) = double(fit.Din(:)); end
residual = double(fit.resnorm(:));
exitflag = double(fit.solverExitflag(:));
[loss,failed,error] = impulsed.fit_error(estimate,target,exitflag,residual,cfg);
estimate = reshape(estimate,1,nt,model.NumParms);
error = reshape(error,1,nt,3);
loss = loss';
failed = failed';
residual = residual';
exitflag = exitflag';
end

function restore_state(state,warning_state)
rng(state);
warning(warning_state);
end
