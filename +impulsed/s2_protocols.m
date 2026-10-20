function protocols = s2_protocols(cfg)
if nargin == 0, cfg = optimization_config; end
names = {'Gen.1','Gen.2','Cust.1','Cust.2','Cust.3'};
x = [17 33 1200 600 30 1200;25 50 1000 300 74 1000; ...
    27 47 1200 400 81 1200;35 53 960 400 48 860;32 50 1200 470 68 750];
reported_TE = [141 104 99 82 87];
pgse_delta = [10 12 12 12 12];
b = {{[0 400 800 1200]',[0 400 800 1200]',[0 300 600]'}, ...
    {[0 250 500 750 1000]',[0 250 500 750 1000]',[0 100 200 300]'}, ...
    {[0 400 800 1200]',[0 400 800 1200]',[0 135 265 400]'}, ...
    {[0 285 575 860]',[0 320 630 960]',[0 135 265 400]'}, ...
    {[0 250 500 750]',[0 400 800 1200]',[0 160 320 470]'}};
for j = 1:numel(names)
    setting = cfg;
    setting.pgse_delta = pgse_delta(j);
    setting.b_fractions = cellfun(@(v)v/v(end),b{j},'UniformOutput',false);
    protocols(j) = struct('name',names{j},'parameters',x(j,:), ... %#ok<AGROW>
        'reported_TE',reported_TE(j),'b_values',{b{j}},'config',setting);
end
end
