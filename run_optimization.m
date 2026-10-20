function results = run_optimization(cfg)
% Surrogate search followed by independent evaluation of the best candidates.
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'mati'),'-begin');
if nargin == 0, cfg = optimization_config; end
state = rng;
cleanup_rng = onCleanup(@()rng(state)); %#ok<NASGU>
rng(cfg.search_seed,'twister');
initial = initial_points(cfg);

% Verify one complete simulation and fit before starting the parallel search.
pilot = impulsed.evaluate_protocol(initial(1,:),cfg.truth,1,cfg.search_seed,cfg,0);
assert(isfinite(pilot.objective) && any(~pilot.failed(:)), ...
    'The single-design pilot did not produce a valid fit.');
workers = 0;
pool = [];
owns_pool = false;
if cfg.use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('Processes',cfg.workers);
        owns_pool = true;
    elseif isa(pool,'parallel.ThreadPool')
        error('optimization:ProcessPoolRequired', ...
            'Use a process pool for surrogateopt or set cfg.use_parallel=false.');
    end
    workers = pool.NumWorkers;
    spmd
        addpath(root,fullfile(root,'mati'),'-begin');
    end
end
cleanup_pool = onCleanup(@()release_pool(pool,owns_pool)); %#ok<NASGU>
if isfolder(cfg.output_dir)
    error('optimization:ExistingOutput','Choose a new cfg.output_dir.');
end
mkdir(cfg.output_dir);
options = optimoptions('surrogateopt','Display',cfg.display, ...
    'UseParallel',cfg.use_parallel,'MaxFunctionEvaluations',cfg.max_evaluations, ...
    'MinSurrogatePoints',20,'InitialPoints',initial, ...
    'ConstraintTolerance',1e-8,'BatchUpdateInterval',max(workers,1));
[~,~,exitflag,output,trials] = surrogateopt(@(x)objective(x,cfg),cfg.lb,cfg.ub,options);
save(fullfile(cfg.output_dir,'search.mat'),'cfg','trials','exitflag','output','pilot');

[~,order] = sort(trials.Fval);
parameters = zeros(0,6);
search_score = zeros(0,1);
for j = order(:)'
    x = trials.X(j,:);
    [~,p] = impulsed.protocol_model(x,cfg);
    if p.feasible && isfinite(trials.Fval(j)) && ~ismember(x,parameters,'rows')
        parameters(end+1,:) = x; %#ok<AGROW>
        search_score(end+1,1) = trials.Fval(j); %#ok<AGROW>
        if size(parameters,1) == cfg.shortlist_size, break; end
    end
end
assert(~isempty(parameters),'The search returned no feasible candidates.');
perturbed = impulsed.sample_tissues(cfg,cfg.validation_repeats,cfg.tissue_seed);
runs = cell(size(parameters,1),1);
for j = 1:numel(runs)
    runs{j} = impulsed.evaluate_protocol(parameters(j,:),perturbed, ...
        cfg.validation_repeats,cfg.validation_seed,cfg,workers);
    fprintf('Validation %d/%d: objective %.6g, SD %.6g, failures %.2f%%\n', ...
        j,numel(runs),runs{j}.objective,runs{j}.error_std,100*runs{j}.failure_rate);
end
values = cellfun(@(r)[r.objective,r.error_std,r.failure_rate],runs,'UniformOutput',false);
[~,ranking] = sortrows(vertcat(values{:}),[1 2 3]);
results.cfg = cfg;
results.parameters = parameters;
results.search_score = search_score;
results.validation = runs;
results.validation_truth = perturbed;
results.ranking = ranking;
results.selected = parameters(ranking(1),:);
results.selection_rule = 'objective_then_error_SD_then_failure_rate';
results.example_settings = cfg.example_settings;
results.search_exitflag = exitflag;
results.search_output = output;
save(fullfile(cfg.output_dir,'results.mat'),'results');
export_table(results,cfg.output_dir);
plot_signal(results,cfg);
fprintf('Completed %s. Results: %s\n',cfg.mode,cfg.output_dir);
end

function out = objective(x,cfg)
[model,p] = impulsed.protocol_model(x,cfg);
out.Ineq = double(~p.feasible);
if p.feasible
    result = impulsed.evaluate_protocol(x,cfg.truth,cfg.search_repeats,cfg.search_seed,cfg,0,model,p);
    out.Fval = result.objective;
else
    out.Fval = 1e6;
end
end

function points = initial_points(cfg)
points = zeros(0,6);
for j = 1:5000
    x = cfg.lb+rand(1,6).*(cfg.ub-cfg.lb);
    [~,p] = impulsed.protocol_model(x,cfg);
    if p.feasible
        points(end+1,:) = x; %#ok<AGROW>
        if size(points,1) == cfg.initial_designs, return; end
    end
end
error('optimization:InitialDesigns','Insufficient feasible initial designs in the configured bounds.');
end

function export_table(results,output_dir)
n = size(results.parameters,1);
extra = zeros(n,7);
for j = 1:n
    r = results.validation{j};
    extra(j,:) = [r.protocol.TE,r.protocol.cycles,r.objective, ...
        r.error_mean,r.error_std,r.failure_rate];
end
rank = zeros(n,1);
rank(results.ranking) = (1:n)';
table_out = array2table([rank,results.parameters,results.search_score,extra], ...
    'VariableNames',{'Rank','LowFrequency_Hz','HighFrequency_Hz', ...
    'LowBmax_smm2','HighBmax_smm2','PGSEDelta_ms','PGSEBmax_smm2', ...
    'SearchObjective','TE_ms','LowCycles','HighCycles','ValidationObjective', ...
    'MeanRelativeSquaredError','ErrorSD','FailureRate'});
table_out = sortrows(table_out,'Rank');
writetable(table_out,fullfile(output_dir,'candidates.csv'));
end

function plot_signal(results,cfg)
[model,p] = impulsed.protocol_model(results.selected,cfg);
target = cfg.truth(ceil(size(cfg.truth,1)/2),:);
signal = model.FcnSignal(target(1:model.NumParms)',model);
fig = figure('Visible','off','Color','w');
cleanup = onCleanup(@()close(fig)); %#ok<NASGU>
hold on
colors = lines(3);
labels = {'PGSE','OGSE low','OGSE high'};
offset = [0,cumsum(p.block_sizes)];
for j = 1:3
    ix = offset(j)+(1:p.block_sizes(j));
    plot(1000*model.pulse.b(ix),signal(ix),'-o', ...
        'LineWidth',1.5,'Color',colors(j,:),'DisplayName',labels{j});
end
xlabel('b (s/mm^2)');
ylabel('Signal / S_0');
title(sprintf('%s, TE %.2f ms',strrep(cfg.tissue,'_',' '),p.TE));
legend('Location','northeast');
grid on
exportgraphics(fig,fullfile(cfg.output_dir,'selected_signal.png'),'Resolution',150);
end

function release_pool(pool,owns_pool)
if owns_pool && ~isempty(pool), delete(pool); end
end
