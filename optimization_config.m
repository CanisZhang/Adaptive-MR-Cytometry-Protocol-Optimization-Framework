function cfg = optimization_config(tissue, mode)
if nargin < 1, tissue = 'human_liver'; end
if nargin < 2, mode = 'full'; end
cfg.tissue = char(tissue);
cfg.mode = char(mode);
cfg.example_settings = true;
cfg.Din = 1.6;
cfg.Gmax = 120e-5;
switch cfg.tissue
    case 'human_liver'
        scenarios = [.20 20 1.8; .40 14 1.6; .55 12 1.1];
        free_Din = true;
    case 'rat_liver'
        scenarios = [.70 22 2.5; .60 25 1.8; .55 18 1.4];
        free_Din = true;
    case 'kidney'
        scenarios = [.40 14 2.4; .20 11 1.6; .60 19 1.0];
        cfg.Din = 2;
        free_Din = false;
    case 'prostate'
        scenarios = [.30 18 2.3; .45 15 1.7; .65 12 1.2];
        cfg.Din = 1;
        cfg.Gmax = 80e-5;
        free_Din = false;
    case 'brain'
        scenarios = [.30 11 1.8; .40 14 1.3; .55 17 1.0];
        cfg.Din = 1;
        free_Din = false;
    otherwise
        error('optimization:UnknownTissue', 'Unknown tissue: %s', cfg.tissue);
end
cfg.truth = [scenarios, cfg.Din*ones(size(scenarios,1),1)];
cfg.perturbation_bounds = [min(scenarios,[],1); max(scenarios,[],1)]';
cfg.perturbation_fraction = .10;
cfg.model_name = 'IMPULSED_vin_d_Dex';
cfg.fitopts.x0 = [.5; 5; 1.5];
cfg.fitopts.lb = [0; 5; 0];
cfg.fitopts.ub = [1; 30; 3];
if free_Din
    cfg.model_name = 'IMPULSED_vin_d_Dex_Din';
    cfg.fitopts.x0(4) = 1.5;
    cfg.fitopts.lb(4) = 0;
    cfg.fitopts.ub(4) = 3;
end
cfg.fitopts.fittingMethod = 'nllsq';
cfg.fitopts.solverName = 'lsqnonlin';
cfg.fitopts.noiseModel = 'none';
cfg.fitopts.NumStarts = 100;
cfg.fitopts.flag = struct('denoise','n','degibbs','n','getADC','n', ...
    'postProcess','n','multistart','y','useGPU','n','parfor','n','deivim','n');
cfg.fitopts.options = optimoptions('lsqnonlin', 'Display','off', ...
    'FunctionTolerance',1e-8, 'OptimalityTolerance',1e-8, ...
    'StepTolerance',1e-8, 'MaxFunctionEvaluations',3000);
cfg.failure_penalty = 10;
cfg.nonfinite_error = 200;
cfg.residual_threshold = 10;
cfg.fit_seed_offset = 310000000;


cfg.slewMax = 200e-5;
cfg.trise = 1;
cfg.timeACQ = 11;
cfg.pgse_delta = 12;
cfg.pgse_te_margin = 18;
cfg.te_raster = 1;
cfg.ogse_separation = 6;
cfg.allowed_cycles = [1 1.5 2 2.5 3];
cfg.frequency_gap = 2;
cfg.b_fractions = repmat({linspace(0,1,4)'},1,3);
cfg.lb = [17 30 300 200 30 300];
cfg.ub = [40 65 1200 1200 90 1200];
cfg.T2 = 35;
cfg.TE0 = 135;
cfg.SNR0 = 15;
cfg.num_averages = 2;

cfg.search_seed = 20260128;
cfg.validation_seed = 13579;
cfg.tissue_seed = 24680;
cfg.search_repeats = 20;
cfg.validation_repeats = 200;
cfg.max_evaluations = 512;
cfg.initial_designs = 20;
cfg.shortlist_size = 10;
cfg.use_parallel = true;
cfg.workers = min(16, feature('numcores'));
cfg.display = 'iter';
switch cfg.mode
    case 'full'
    case 'smoke'
        cfg.fitopts.NumStarts = 2;
        cfg.search_repeats = 1;
        cfg.validation_repeats = 3;
        cfg.max_evaluations = 24;
        cfg.shortlist_size = 2;
        cfg.use_parallel = false;
        cfg.display = 'off';
    otherwise
        error('optimization:UnknownMode', 'Unknown mode: %s', cfg.mode);
end
cfg.output_dir = fullfile(fileparts(mfilename('fullpath')), 'results', ...
    [char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')) '_' cfg.tissue '_' cfg.mode]);
end
