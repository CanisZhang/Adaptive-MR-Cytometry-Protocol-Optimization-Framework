function truth = sample_tissues(cfg, n, seed)
% Truncated normal perturbations of the paired tissue scenarios.
state = rng;
cleanup = onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister');
truth = repmat(cfg.truth,1,1,n);
for j = 1:3
    center = cfg.truth(:,j);
    sigma = cfg.perturbation_fraction*center;
    bounds = cfg.perturbation_bounds(j,:);
    if any(sigma <= 0) || bounds(1) == bounds(2)
        continue
    end
    lower = .5*(1+erf((bounds(1)-center)./(sqrt(2)*sigma)));
    upper = .5*(1+erf((bounds(2)-center)./(sqrt(2)*sigma)));
    probability = lower + rand(numel(center),n).*(upper-lower);
    values = center + sqrt(2)*sigma.*erfinv(2*probability-1);
    truth(:,j,:) = reshape(values,numel(center),1,n);
end
end
