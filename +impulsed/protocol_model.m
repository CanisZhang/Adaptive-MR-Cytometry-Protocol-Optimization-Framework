function [model, p] = protocol_model(x, cfg)
% x = [f_low f_high b_low b_high Delta_PGSE b_PGSE], in Hz, s/mm^2 and ms.
x = x(:)';
model = [];
p = struct('feasible',false,'parameters',x,'TE',NaN,'cycles',[NaN NaN], ...
    'peak_gradient',NaN);
if numel(x) ~= 6 || any(~isfinite(x)) || any(x <= 0) || ...
        any(x < cfg.lb | x > cfg.ub) || x(2)-x(1) < cfg.frequency_gap
    return
end
cap = min(cfg.Gmax,cfg.slewMax*cfg.trise);
cycles = cfg.allowed_cycles(:);
duration = cfg.trise+1000*cycles./x(1:2);
maxb = zeros(numel(cycles),2);
for g = 1:2
    pulse = mati.DiffusionPulseSequence(numel(cycles), ...
        'TE',2*(duration(:,g)+cfg.timeACQ), ...
        'delta',duration(:,g),'Delta',duration(:,g)+cfg.ogse_separation, ...
        'G',cap,'n',cycles,'shape',"tcos",'trise',cfg.trise);
    maxb(:,g) = 1000*pulse.b;
end
pgse = mati.DiffusionPulseSequence(1,'TE',x(5)+cfg.pgse_te_margin, ...
    'delta',cfg.pgse_delta,'Delta',x(5),'G',cap, ...
    'shape',"tpgse",'trise',cfg.trise);
if 1000*pgse.b < x(6)-1e-8, return; end

choices = zeros(0,5);
for low = find(maxb(:,1) >= x(3)-1e-8)'
    for high = find(maxb(:,2) >= x(4)-1e-8)'
        delta = max(duration(low,1),duration(high,2));
        raw_TE = max([2*(delta+cfg.timeACQ),x(5)+cfg.pgse_te_margin, ...
            2*delta+cfg.ogse_separation,x(5)+cfg.pgse_delta]);
        TE = cfg.te_raster*ceil(raw_TE/cfg.te_raster-1e-12);
        choices(end+1,:) = [TE,-cycles(low),-cycles(high),low,high]; %#ok<AGROW>
    end
end
if isempty(choices), return; end
choices = sortrows(choices,[1 2 3]);
TE = choices(1,1);
index = choices(1,4:5);
N = [0;cycles(index(:))];
delta = [cfg.pgse_delta;duration(index(1),1);duration(index(2),2)];
separation = [x(5);delta(2:3)+cfg.ogse_separation];
bmax = x([6 3 4]);
shape = ["tpgse","tcos","tcos"];
pulses = cell(1,3);
for g = 1:3
    pulses{g} = mati.DiffusionPulseSequence(numel(cfg.b_fractions{g}), ...
        'TE',TE,'delta',delta(g),'Delta',separation(g), ...
        'b',bmax(g)*cfg.b_fractions{g}/1000,'n',N(g), ...
        'shape',shape(g),'trise',cfg.trise);
end
pulse = mati.PulseSequence.cat(pulses{:});
if max(pulse.G) > cap+1e-12 || any(pulse.delta+pulse.Delta > TE+1e-10)
    return
end
structure = struct('modelName',cfg.model_name,'betaex',0,'geometry','sphere');
if strcmp(cfg.model_name,'IMPULSED_vin_d_Dex'), structure.Din = cfg.Din; end
model = mati.IMPULSED(structure,pulse);
p.feasible = true;
p.TE = TE;
p.cycles = N(2:3)';
p.delta = delta(2:3)';
p.block_sizes = cellfun(@numel,cfg.b_fractions);
p.b_smm2 = 1000*pulse.b;
p.frequency_Hz = 1000*pulse.f;
p.peak_gradient = max(pulse.G);
p.effective_SNR = cfg.SNR0*sqrt(cfg.num_averages)*exp(-(TE-cfg.TE0)/cfg.T2);
end
