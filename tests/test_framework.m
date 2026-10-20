function report = test_framework
% Focused numerical checks for the submission framework.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'mati'), '-begin');
c = optimization_config('human_liver', 'full');
assert(c.trise == 1 && c.timeACQ == 11 && c.pgse_te_margin == 18);
assert(c.te_raster == 1);
assert(c.fitopts.NumStarts == 100);
assert(isequal(c.truth(:,1:3), [.20 20 1.8; .40 14 1.6; .55 12 1.1]));
assert(isequal(c.fitopts.lb, [0;5;0;0]));
assert(isequal(c.fitopts.ub, [1;30;3;3]));
assert(c.shortlist_size == 10);
assert(strcmp(c.model_name, 'IMPULSED_vin_d_Dex_Din'));
k = optimization_config('kidney', 'full');
assert(k.Din == 2 && numel(k.fitopts.lb) == 3);
b = optimization_config('brain', 'full');
assert(b.Din == 1 && numel(b.fitopts.lb) == 3);
p = optimization_config('prostate', 'full');
assert(p.Din == 1 && p.Gmax == 80e-5);

truth = [.4 16 1.6 1.6];
[loss, failed] = impulsed.fit_error(truth, truth, 1, 0, c);
assert(loss == 0 && ~failed);
changed = truth;
changed(4) = 2.8;
[loss, failed] = impulsed.fit_error(changed, truth, 1, 0, c);
assert(loss == 0 && ~failed, 'Din is fitted but is outside the three-term objective.');
changed = truth;
changed(1:3) = 1.1*changed(1:3);
[loss, failed] = impulsed.fit_error(changed, truth, 1, 0, c);
assert(abs(loss-.03) < 1e-12 && ~failed);
changed = truth;
changed(2) = 4.99;
[loss, failed] = impulsed.fit_error(changed, truth, 1, 0, c);
assert(failed && loss >= c.failure_penalty);
[~, failed] = impulsed.fit_error(truth, truth, 1, 10, c);
assert(~failed);
[~, failed] = impulsed.fit_error(truth, truth, 1, 10+eps(10), c);
assert(failed);
[loss, failed] = impulsed.fit_error(nan(1,4), truth, -1, NaN, c);
assert(failed && isfinite(loss) && loss >= c.failure_penalty);

x = [25 50 1000 300 74 1000];
[model, protocol] = impulsed.protocol_model(x, c);
assert(protocol.feasible && ~isempty(model));
assert(protocol.peak_gradient <= c.Gmax+1e-12);
assert(protocol.peak_gradient/c.trise <= c.slewMax+1e-12);
assert(all(model.pulse.delta+model.pulse.Delta <= protocol.TE+1e-10));
assert(max(abs(1000*model.pulse.b([4 8 12])'-[1000 1000 300])) < 1e-9);
assert(max(abs(1000*model.pulse.f([8 12])'-[25 50])) < 1e-9);
assert(all(model.FcnSignal(truth', model) > 0));
bad = x;
bad(6) = 10000;
[~, invalid] = impulsed.protocol_model(bad, c);
assert(~invalid.feasible);

cases = impulsed.s2_protocols(c);
calculated_TE = zeros(5,1);
for j = 1:numel(cases)
    [reference_model, reference] = impulsed.protocol_model(cases(j).parameters,cases(j).config);
    assert(reference.feasible);
    calculated_TE(j) = reference.TE;
    assert(reference.peak_gradient <= c.Gmax+1e-12);
    assert(reference.peak_gradient/c.trise <= c.slewMax+1e-12);
    expected_b = vertcat(cases(j).b_values{:});
    assert(max(abs(1000*reference_model.pulse.b-expected_b)) < 1e-9);
    assert(cases(j).config.timeACQ == c.timeACQ);
end
assert(isequal(calculated_TE,[142;104;99;82;87]));
assert(isequal(calculated_TE-[cases.reported_TE]',[1;0;0;0;0]));

samples = impulsed.sample_tissues(c, 25, 5100);
assert(isequal(size(samples), [3 4 25]));
for j = 1:3
    bounds = c.perturbation_bounds(j,:);
    column = samples(:,j,:);
    assert(all(column(:) >= bounds(1) & column(:) <= bounds(2)));
end
assert(all(reshape(samples(:,4,:), [], 1) == c.Din));
assert(isequal(samples, impulsed.sample_tissues(c, 25, 5100)));

quick = optimization_config('brain', 'smoke');
quick.SNR0 = Inf;
quick.fitopts.NumStarts = 5;
quick.truth = [.4 16 1.6 1];
out = impulsed.evaluate_protocol(x, quick.truth, 2, 6100, quick, 0);
assert(~any(out.failed(:)));
assert(max(abs(out.relative_error(:))) < .01);
assert(out.objective < 1e-4);
assert(isequal(size(out.estimates), [2 1 3]));
repeat = impulsed.evaluate_protocol(x, quick.truth, 2, 6100, quick, 0);
assert(isequal(out.estimates, repeat.estimates));
[cached_model,cached_protocol] = impulsed.protocol_model(x,quick);
cached = impulsed.evaluate_protocol(x,quick.truth,2,6100,quick,0,cached_model,cached_protocol);
assert(isequaln(out,cached));
report = struct('passed', true, 'noiseless_objective', out.objective, ...
    'protocol_TE_ms', protocol.TE, 'peak_gradient_mTm', protocol.peak_gradient*1e5, ...
    's2_TE_ms',calculated_TE);
fprintf('PASS: unified S2 timing, hardware, exact b grids, loss, perturbation and fitting.\n');
end
