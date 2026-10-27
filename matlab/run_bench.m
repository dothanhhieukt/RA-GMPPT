function R = run_bench(group, nMax, retryFailed, onlyIds)
% =========================================================================
% RUN_BENCH  Batch simulation of the switched Simulink model (Charger_Buck.slx).
%   Results are saved after every chunk; an interrupted group CONTINUES where it stopped.
% =========================================================================
% Order:
%   run_bench('sweep')     % 1) quasi-static P-V sweep in the circuit -> reference P_reach (RUN FIRST)
%   run_bench('test')      % 2) quick test: 5 controllers x 72-C4-S60, 48-C1-S60 x 1 run
%   run_bench('static')    % 3) 14 scenarios x (P&O, DPSO, RA-DPSO x1; CS, RA-CS x10)
%   run_bench('ablRA')     %    ablation of the RA layer: R1 / R1+HT / no RA, 2 cores, 6 scenarios (CS x5)
%   run_bench('ablLim')    %    conventional limiter (RA kept), 3 CC scenarios
%   run_bench('dynamic')   %    4 irradiance steps, 1.1 s
%   run_bench('vrise')     %    rising battery voltage makes the peak unreachable, 1.12 s (with ablation variants)
%   run_bench('vrise2')    %    rising battery voltage makes an INTERIOR held peak unreachable (72-C7)
%   run_bench('cv')        %    CC -> CV, 1 scenario
% Additional groups:
%   run_bench('adaptRange')     %    adaptive-range baselines and SCAN: DPSO-AR, CS-AR, SCAN, SCAN-AR, RA-SCAN x 10 MPPT scenarios
%   run_bench('vriseAR')     %    rising battery voltage with SCAN, RA-SCAN, RA-SCAN-R1, DPSO-AR, CS-AR
%   run_bench('limiter')    %    limiter variants x 7 CC scenarios (3 heavy-shading scenarios with a low Ilim)
%   run_bench('dynLimiter')    %    irradiance steps + CC->CV with the fast left-flank limiter / PI min-select
%   run_bench('sens')   %    sensitivity of RA-DPSO (Dmax, Vbat offset, dVc, dVR, R2 filter) x 5 scenarios
%   run_bench('fs20k')     %    20 kHz sampling, DPSO and RA-DPSO x 12 scenarios
% Options:
%   run_bench('static', 8)            % run at most 8 runs in this call
%   ids = run_bench('static','ids')   % list the run ids of a group (no simulation)
%   run_bench('static', Inf, true)    % re-run the FAILED runs (smaller MaxStep)
%   run_bench('static', Inf, false, {'static|CS|72-C4-S60|r03'})   % re-run given runs and compare
% Output: <cfg.outDir>/<group>.mat (table R); waveforms: <cfg.outDir>/traces/<group>/<id>.mat
% Needs: Charger_Bench.m (MATLAB Function block), battery_35E.m, the model Charger_Buck.slx with the
%   workspace variables alg, abl, seed0, Gsig, Tsig, Ilim, Vcv, bat, Tc and the To Workspace blocks
%   Vpv, Ipv, Vbat, Ibat, Dlog, Vref, phase, mode, SOC.
% =========================================================================
if nargin >= 2 && ischar(nMax) && strcmp(nMax,'ids')
    jj = build_jobs(bench_config(), group);  R = {jj.id}';  return;
end
if nargin < 2, nMax = Inf; end
if nargin < 3, retryFailed = false; end
if nargin < 4, onlyIds = {}; end
cfg  = bench_config();
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
refs = [];
if ~strcmp(group,'sweep')
    fs = fullfile(cfg.outDir,'sweep.mat');
    if ~exist(fs,'file'), error('No reference: run run_bench(''sweep'') first.'); end
    Q = load(fs,'R');  refs = Q.R(Q.R.ok, :);
end
jobs = build_jobs(cfg, group);
if ~isempty(refs)
    need = unique(arrayfun(@(j) refkey(j.scen.pack, j.scen.caseB, j.scen.socRef), jobs, 'UniformOutput', false));
    miss = setdiff(need, refs.key);
    if ~isempty(miss)
        error('Missing P-V reference: %s -> run run_bench(''sweep'') first.', strjoin(miss, ', '));
    end
end
fres  = fullfile(cfg.outDir, [group '.mat']);
trDir = fullfile(cfg.outDir, 'traces', group);
if ~exist(trDir,'dir'), mkdir(trDir); end

R = table();
if exist(fres,'file'), S = load(fres,'R'); R = S.R; end
if ~isempty(R)
    bad = ~R.ok;
    if any(bad)
        fprintf('%d FAILED runs in the saved results:\n', nnz(bad));
        disp(R(bad, {'id','err'}));
        if retryFailed
            R = R(~bad, :);  save(fres, 'R');
            fprintf('-> Failed runs removed, they will be re-run.\n');
        else
            fprintf('-> Call run_bench(''%s'', Inf, true) to re-run these runs.\n', group);
        end
    end
end
oldRows = table();
if ~isempty(onlyIds) && ~isempty(R)
    sel = ismember(R.id, onlyIds);
    oldRows = R(sel, :);  R = R(~sel, :);
end
done = false(numel(jobs),1);
if ~isempty(R), done = ismember({jobs.id}', R.id); end
todo = find(~done);
if ~isempty(onlyIds), todo = todo(ismember({jobs(todo).id}', onlyIds)); end
todo = todo(1:min(numel(todo), nMax));
fprintf('Group "%s": %d runs, %d done, %d in this call.\n', group, numel(jobs), nnz(done), numel(todo));
if isempty(todo), return; end

ms = cfg.MaxStep;
if retryFailed
    ms = cfg.MaxStepRetry;  fprintf('Re-run with MaxStep = %s, adaptive zero-crossing\n', ms);
end
usePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
nChunk = 4;
if usePar
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool; end
    nChunk = cfg.nChunk;
    fprintf('Pool: %d workers, %d runs per chunk.\n', pool.NumWorkers, nChunk);
end
load_system(cfg.model);

% If the pool hangs or dies (e.g. "No workers are available for FevalQueue execution"), the run is NOT
%   recorded as failed but QUEUED AGAIN at the end, and the pool is recreated. The pool is also restarted
%   every POOL_RECYCLE chunks to release memory (run times grow when a pool runs for a long time).
POOL_RECYCLE = 3;  MAX_REQ = 2;
nReq = zeros(numel(jobs),1);  queue = todo(:)';  nDone = 0;  nTot = numel(todo);  kChunk = 0;
while ~isempty(queue)
    idx = queue(1:min(nChunk, numel(queue)));  queue(1:numel(idx)) = [];
    kChunk = kChunk + 1;
    if usePar
        pool = gcp('nocreate');
        if ~isempty(pool) && (~pool.Connected || (kChunk > 1 && mod(kChunk-1, POOL_RECYCLE) == 0))
            fprintf('  (restarting the pool)\n');  delete(pool);  pool = [];
        end
        if isempty(pool), pool = parpool; end %#ok<NASGU>
    end
    in  = Simulink.SimulationInput.empty(0,1);
    for i = 1:numel(idx)
        j = jobs(idx(i));
        bat = battery_35E(j.scen.pack, j.soc0, j.Qscale);
        s = Simulink.SimulationInput(cfg.model);
        s = s.setVariable('alg',   j.alg);
        s = s.setVariable('abl',   j.abl);
        s = s.setVariable('seed0', j.seed);
        s = s.setVariable('Gsig',  j.Gsig);
        s = s.setVariable('Tsig',  j.Tsig);
        IlimJ = cfg.Ilim;  if isfield(j.scen,'Ilim') && ~isempty(j.scen.Ilim), IlimJ = j.scen.Ilim; end
        s = s.setVariable('Ilim',  IlimJ);
        s = s.setVariable('Vcv',   bat.Vcv);
        s = s.setVariable('bat',   bat);
        s = s.setVariable('Tc',    cfg.Tc);
        s = s.setModelParameter('StopTime', num2str(j.Tend), 'MaxStep', ms);
        if retryFailed
            s = s.setModelParameter('ZeroCrossAlgorithm','Adaptive', 'IgnoredZcDiagnostic','none');
        end
        in(end+1,1) = s; %#ok<AGROW>
    end
    t0 = tic;
    if usePar
        out = parsim(in, 'ShowProgress','off', 'StopOnError','off');
    else
        out = Simulink.SimulationOutput.empty(0,1);
        for i = 1:numel(in)
            try
                out(i,1) = sim(in(i)); %#ok<AGROW>
            catch ME
                out(i,1) = Simulink.SimulationOutput; %#ok<AGROW>
                warning('Run %s failed: %s', jobs(idx(i)).id, ME.message);
            end
        end
    end
    infra = false(numel(idx),1);
    for i = 1:numel(idx)
        em = '';  try, em = char(out(i).ErrorMessage); catch, end
        infra(i) = ~isempty(regexpi(em, 'No workers|FevalQueue|parallel pool|worker.*(crash|abort|died|unavailable)|Out of memory', 'once'));
    end
    if any(infra)
        rq = idx(infra);  rq = rq(nReq(rq) < MAX_REQ);  nReq(rq) = nReq(rq) + 1;
        queue = [queue, rq]; %#ok<AGROW>
        fprintf('  %d runs failed because of the POOL (not the simulation) -> %d queued again, pool recreated.\n', nnz(infra), numel(rq));
        if usePar, delete(gcp('nocreate')); end
    end
    for i = 1:numel(idx)
        if infra(i) && nReq(idx(i)) <= MAX_REQ && ismember(idx(i), queue), continue; end
        [M, Tr] = bench_metrics(out(i), jobs(idx(i)), cfg, refs);
        R = addrow(R, struct2table(M, 'AsArray', true));
        if ~isempty(Tr)
            ftr = fullfile(trDir, [regexprep(jobs(idx(i)).id, '[|>&]', '_') '.mat']);
            save(ftr, '-struct', 'Tr');
        end
    end
    save(fres, 'R');
    tw = NaN(numel(idx),1);
    for i = 1:numel(idx)
        try, tw(i) = out(i).SimulationMetadata.TimingInfo.ExecutionElapsedWallTime; catch, end
    end
    nDone = nDone + nnz(~(infra & ismember(idx(:), queue)));
    fprintf('  done %d/%d  (%.1f min for this chunk; one simulation: median %.0f s, max %.0f s)\n', ...
            nDone, nTot, toc(t0)/60, median(tw,'omitnan'), max(tw));
end
if ~isempty(oldRows)
    fprintf('\nCOMPARISON of the old (before) and new (after) results:\n');
    new = R(ismember(R.id, oldRows.id), :);
    [~, ia, ib] = intersect(oldRows.id, new.id, 'stable');
    if strcmp(group,'sweep')
        C = table(oldRows.id(ia), oldRows.Preach(ia), new.Preach(ib), oldRows.Vreach(ia), new.Vreach(ib), ...
                  oldRows.Vbnd(ia), new.Vbnd(ib), 'VariableNames', ...
                  {'id','Preach_old','Preach_new','Vreach_old','Vreach_new','Vbnd_old','Vbnd_new'});
        disp(C);  return;
    end
    C = table(oldRows.id(ia), oldRows.Vss(ia), new.Vss(ib), oldRows.Pss(ia), new.Pss(ib), ...
              oldRows.etaE(ia), new.etaE(ib), oldRows.Ebat(ia), new.Ebat(ib), ...
              'VariableNames', {'id','Vss_old','Vss_new','Pss_old','Pss_new','etaE_old','etaE_new','Ebat_old','Ebat_new'});
    disp(C);
end
end

% =========================================================================
function jobs = build_jobs(cfg, group)
jobs = struct('id',{},'group',{},'algName',{},'alg',{},'abl',{},'scen',{}, ...
              'seed',{},'soc0',{},'Qscale',{},'Tend',{},'tsw',{},'Gsig',{},'Tsig',{});
A = cfg.algs;
switch group
    case 'sweep'
        allS = [cfg.static; cfg.dynamic; cfg.vrise; cfg.vrise2; cfg.cv];
        keys = {};
        for k = 1:numel(allS)
            s = allS(k);
            key = refkey(s.pack, s.caseB, s.socRef);
            if ~ismember(key, keys)
                keys{end+1} = key; %#ok<AGROW>
                ss = s;  ss.caseA = s.caseB;  ss.tsw = Inf;  ss.Tend = cfg.Tsweep;
                ss.soc0 = s.socRef;  ss.Qscale = 1;  ss.name = key;
                jobs(end+1) = mkjob(cfg, group, struct('name','SWEEP','alg',-2,'abl',0), ss, 1); %#ok<AGROW>
            end
        end
    case 'test'
        for a = 1:numel(A), for nm = {'72-C4-S60','48-C1-S60'}
            s = cfg.static(strcmp({cfg.static.name}, nm{1}));
            jobs(end+1) = mkjob(cfg, group, A(a), s, 1); %#ok<AGROW>
        end, end
    case 'static'
        for a = 1:numel(A), for k = 1:numel(cfg.static)
            for r = 1:nrun(cfg, A(a), cfg.nRun.static)
                jobs(end+1) = mkjob(cfg, group, A(a), cfg.static(k), r); %#ok<AGROW>
            end
        end, end
    case {'ablRA','ablLim'}
        VV = cfg.(group);  L = cfg.([group 'Scen']);
        for v = 1:numel(VV), for k = 1:numel(L)
            s = cfg.static(strcmp({cfg.static.name}, L{k}));
            for r = 1:nrun(cfg, VV(v), cfg.nRun.(group))
                jobs(end+1) = mkjob(cfg, group, VV(v), s, r); %#ok<AGROW>
            end
        end, end
    case {'dynamic','vrise','vrise2','cv'}
        SS = cfg.(group);
        switch group
            case 'dynamic', L = cfg.dynAlgs;
            case 'vrise',   L = cfg.vriseAlgs;
            case 'vrise2',  L = cfg.vrise2Algs;
            otherwise,      L = cfg.cvAlgs;
        end
        AL = [A; cfg.ablRA; cfg.ablNoHT];  anL = {AL.name};          % ablation variants allowed
        for a = 1:numel(L), for k = 1:numel(SS)
            Aa = AL(strcmp(anL, L{a}));
            for r = 1:nrun(cfg, Aa, cfg.nRun.(group))
                jobs(end+1) = mkjob(cfg, group, Aa, SS(k), r); %#ok<AGROW>
            end
        end, end
    case {'adaptRange','vriseAR','limiter','dynLimiter','sens','fs20k'}      % additional groups
        RV = cfg.(group);
        for q = 1:numel(RV)
            Aa = RV(q).alg;  SS = RV(q).scen;
            for a = 1:numel(Aa), for k = 1:numel(SS)
                for r = 1:nrun(cfg, Aa(a), RV(q).nRun)
                    jobs(end+1) = mkjob(cfg, group, Aa(a), SS(k), r); %#ok<AGROW>
                end
            end, end
        end
    otherwise
        error('Unknown group: %s', group);
end
end

function n = nrun(cfg, A, nDefault)
if ismember(A.name, cfg.detAlgs) || ismember(A.alg, [1 5 6 9 10]), n = 1; return; end   % deterministic controllers: one run
n = nDefault;
end

function j = mkjob(cfg, group, A, s, r)
ca = cfg.cases(strcmp({cfg.cases.name}, s.caseA));
cb = cfg.cases(strcmp({cfg.cases.name}, s.caseB));
if isinf(s.tsw)
    Gsig = [0 ca.G; s.Tend ca.G];   Tsig = [0 ca.T; s.Tend ca.T];
else
    Gsig = [0 ca.G; s.tsw ca.G; s.tsw+1e-6 cb.G; s.Tend cb.G];
    Tsig = [0 ca.T; s.tsw ca.T; s.tsw+1e-6 cb.T; s.Tend cb.T];
end
j = struct('id', sprintf('%s|%s|%s|r%02d', group, A.name, s.name, r), ...
           'group',group, 'algName',A.name, 'alg',A.alg, 'abl',A.abl, 'scen',s, ...
           'seed',r, 'soc0',s.soc0, 'Qscale',s.Qscale, 'Tend',s.Tend, 'tsw',s.tsw, ...
           'Gsig',Gsig, 'Tsig',Tsig);
end

function k = refkey(pack, caseName, soc)
k = sprintf('%d|%s|S%02d', pack, caseName, round(100*soc));
end

% =========================================================================
function R = addrow(R, r)
if isempty(R), R = r; return; end
for v = setdiff(r.Properties.VariableNames, R.Properties.VariableNames)
    R.(v{1}) = NaN(height(R), 1);
end
for v = setdiff(R.Properties.VariableNames, r.Properties.VariableNames)
    r.(v{1}) = NaN(height(r), 1);
end
R = [R; r(:, R.Properties.VariableNames)];
end
