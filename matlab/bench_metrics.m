function [M, Tr] = bench_metrics(o, j, cfg, refs)
% =========================================================================
% BENCH_METRICS  Indices of one switched-model run.
% =========================================================================
% Group 'sweep' (alg = -2): REFERENCE from the P-V curve swept in the circuit:
%   Preach, Vreach : largest reachable power and the voltage at which it occurs
%   Vlo, Vhi       : "basin" of the best peak (connected region around Vreach with P > P2)
%   P2             : second-best reachable peak (> 5 V away); NaN if there is only one peak
%   Vbnd, Pbnd     : physical bound (D = 0.95) and the power there
% Other groups (relative to the reference of the FINAL state):
%   isCC     : power-surplus scenario (Preach*ETA > Vbat*Ilim)
%   success  : MPPT: Pss >= effMin*Preach (inBasin: Vss in [Vlo, Vhi], for information only);  CC: |Ibat_ss - Ilim| <= tolI*Ilim
%   eff      : Pss / Preach (MPPT);  etaE: integral of Ppv / (Preach * duration)
%   Ebat     : energy into the battery (J) from tStart;  IbErr: Ibat_ss - Ilim (CC);  VbErr: Vbat_ss - Vcv
%   (group 'cv': success = |VbErr| <= 0.2 V and dVcv <= 0.5 V)
%   side     : CC: 1 = left flank (Vss < Vreach), -1 = right flank
%   tconv    : time until P (1 ms average) enters the +-2 % band around Pss and stays 20 ms
%   nEval    : number of Vref changes before convergence (meaningful in MPPT only)
%   tUnreach : time with Vref < Vbat/0.95 - 1.5 V (UNREACHABLE reference)
%   tDsat    : time with D >= 0.949 (saturation)
%   tOver    : time with Ibat (50 us average) > overI*Ilim;  Qover: integral of (Ibat - Ilim)+ (A.s)
%   (CC: a higher Ebat can come from EXCEEDING Ilim -> CC runs are compared by tOver/Qover, not by Ebat)
%   dVcv     : CV overshoot IN THE CV PHASE: max(Vbat, 1 ms average) - Vcv, counted from the moment CC is
%              established (|Ibat 1 ms average - Ilim| <= 2 % CONTINUOUSLY for >= 10 ms) and Vbat has reached Vcv - 0.02 V
%   Ipk      : peak charging current (Ibat 1 ms average, from tStart)
%   tSet     : time at which CC is established (see dVcv)
%   nLim     : number of entries into the limiter;  nRescan: number of global searches (hold -> first evaluation) after tStart
%   tHold    : first entry into the hold phase (from tStart);  fAudit: share of time outside the hold afterwards
% =========================================================================
isStep = ~isinf(j.tsw);
if isfield(j.scen,'Ilim') && ~isempty(j.scen.Ilim), cfg.Ilim = j.scen.Ilim; end   % Ilim of the scenario
tStart = 0;  if isStep, tStart = j.tsw; end
s = j.scen;
if strcmp(j.group,'sweep')
    M = struct('id',j.id,'key',refkey(s.pack,s.caseB,s.socRef),'pack',s.pack,'case',s.caseB, ...
               'soc',s.socRef,'ok',false,'Preach',NaN,'Vreach',NaN,'Vlo',NaN,'Vhi',NaN,'P2',NaN, ...
               'Vbnd',NaN,'Pbnd',NaN,'VbatSw',NaN,'Pgmpp',NaN,'err',"");
else
    M = struct('id',j.id,'group',j.group,'alg',j.algName,'abl',j.abl,'scen',s.name,'seed',j.seed, ...
               'ok',false,'isCC',false,'success',false,'inBasin',false,'Vss',NaN,'Pss',NaN,'eff',NaN,'etaE',NaN, ...
               'Ebat',NaN,'Ibss',NaN,'IbErr',NaN,'VbErr',NaN,'side',NaN,'tconv',NaN,'nEval',NaN, ...
               'tUnreach',NaN,'tDsat',NaN,'tOver',NaN,'Qover',NaN,'dVcv',NaN,'Ipk',NaN,'nLim',NaN,'nRescan',NaN,'tSet',NaN, ...
               'tHold',NaN,'fAudit',NaN,'Preach',NaN,'err',"");
end
Tr = [];
try
    if ~isempty(o.ErrorMessage), M.err = string(o.ErrorMessage); return; end
catch
end
try
    [t, V]  = getsig(o,'Vpv');
    Ip = sig_on(o,'Ipv',t);  Vb = sig_on(o,'Vbat',t);  Ib = sig_on(o,'Ibat',t);
    Dl = sig_on(o,'Dlog',t); Vr = sig_on(o,'Vref',t);
catch ME
    M.err = string(ME.message);  return;
end
ph = opt_on(o,'phase',t);  md = opt_on(o,'mode',t);  so = opt_on(o,'SOC',t);
if t(end) < j.Tend - 1e-4
    M.err = string(sprintf('Stopped early at %.4f s', t(end)));  return;
end
dt   = median(diff(t));
nMA1 = max(1, round(1e-3/dt));  nMA5 = max(1, round(50e-6/dt));
P    = V.*Ip;   Pb = Vb.*Ib;
P1   = movmean(P, nMA1);  V1 = movmean(V, nMA1);  Vb1 = movmean(Vb, nMA1);
Ib50 = movmean(Ib, nMA5);

% ======================= REFERENCE (P-V sweep) ============================
if strcmp(j.group,'sweep')
    ir = t >= 0.035;
    vv = V1(ir);  pp = P1(ir);
    Vg = (floor(min(vv)*4)/4 : 0.25 : ceil(max(vv)*4)/4)';
    bi = min(numel(Vg), max(1, round((vv - Vg(1))/0.25) + 1));
    Pg = accumarray(bi, pp, [numel(Vg) 1], @max, NaN);
    [Preach, ib] = max(Pg);
    Vbnd = min(vv);  Pbnd = median(pp(vv <= Vbnd + 0.2));
    % local peaks (+-2 V), and the bound if P rises towards it
    w = 8;  pkV = [];  pkP = [];
    for i = 1:numel(Vg)
        if isnan(Pg(i)) || Pg(i) < 0.3*Preach, continue; end
        a = max(1, i-w);  b = min(numel(Vg), i+w);
        if Pg(i) >= max(Pg(a:b)), pkV(end+1) = Vg(i); pkP(end+1) = Pg(i); end %#ok<AGROW>
    end
    iv = find(~isnan(Pg), 1);
    if iv < numel(Vg) && Pg(iv) >= 0.3*Preach && Pg(iv) > max(Pg(iv+1:min(numel(Vg), iv+w)))
        pkV(end+1) = Vg(iv);  pkP(end+1) = Pg(iv);
    end
    far = abs(pkV - Vg(ib)) > 5;
    P2 = NaN;  if any(far), P2 = max(pkP(far)); end
    thr = P2;  if isnan(thr), thr = 0.9*Preach; end
    lo = ib;  while lo > 1 && (isnan(Pg(lo-1)) || Pg(lo-1) > thr), lo = lo - 1; end
    hi = ib;  while hi < numel(Vg) && (isnan(Pg(hi+1)) || Pg(hi+1) > thr), hi = hi + 1; end
    cb = cfg.cases(strcmp({cfg.cases.name}, s.caseB));
    M.Preach = Preach;  M.Vreach = Vg(ib);  M.Vlo = Vg(lo);  M.Vhi = Vg(hi);  M.P2 = P2;
    M.Vbnd = Vbnd;  M.Pbnd = Pbnd;  M.VbatSw = median(Vb(ir));  M.Pgmpp = cb.Pg;
    M.ok = true;
    Tr = struct('id', j.id, 'Vg', single(Vg), 'Pg', single(Pg), 'pkV', pkV, 'pkP', pkP, ...
                't', single(dec(t,cfg.decim)), 'V', single(dec(V,cfg.decim)), 'P', single(dec(P,cfg.decim)), ...
                'Vb', single(dec(Vb,cfg.decim)), 'D', single(dec(Dl,cfg.decim)));
    return;
end

% ======================= CONTROLLER INDICES ===============================
key = refkey(s.pack, s.caseB, s.socRef);
rf  = refs(strcmp(refs.key, key), :);
if isempty(rf), M.err = "Thieu tham chieu " + key;  return; end
rf = rf(1,:);  M.Preach = rf.Preach;
bat = battery_35E(s.pack, s.soc0, s.Qscale);
iw  = t >= (j.Tend - cfg.winSS);   ie = t >= tStart;
M.Vss  = median(V(iw));  M.Pss = median(P1(iw));  M.Ibss = mean(Ib(iw));
Vbss   = mean(Vb(iw));
M.isCC = rf.Preach*cfg.ETA > 1.01*Vbss*cfg.Ilim;
M.Ebat = trapz(t(ie), Pb(ie));
Ib1  = movmean(Ib, nMA1);
M.Ipk  = max(Ib1(ie));
kCC  = cc_start(t, Ib1, ie, cfg);                                        % CC established (held >= 10 ms)
if ~isempty(kCC), M.tSet = t(kCC) - tStart; end                         
kCV  = [];
if ~isempty(kCC), kCV = find((1:numel(t))' >= kCC & Vb1 >= bat.Vcv - 0.02, 1); end
if isempty(kCV), M.dVcv = max(Vb1(ie)) - bat.Vcv; else, M.dVcv = max(Vb1(kCV:end)) - bat.Vcv; end
M.VbErr = Vbss - bat.Vcv;
if strcmp(j.group,'cv')                 % CC -> CV: success if Vbat tracks Vcv without exceeding it by > 0.5 V
    M.success = abs(M.VbErr) <= 0.2 && M.dVcv <= 0.5;
    M.IbErr   = M.Ibss - cfg.Ilim;
elseif M.isCC
    M.IbErr   = M.Ibss - cfg.Ilim;
    M.success = abs(M.IbErr) <= cfg.tolI*cfg.Ilim;
    M.side    = 1 - 2*double(M.Vss >= rf.Vreach);
else
    M.eff     = M.Pss/rf.Preach;
    M.etaE    = trapz(t(ie), P(ie)) / (rf.Preach*(t(end) - tStart));
    M.inBasin = (M.Vss >= rf.Vlo) && (M.Vss <= rf.Vhi);
    M.success = M.eff >= cfg.effMin;   % by POWER (e.g. 48-C7-S60: bound 351.95 W and 85 V peak 350.8 W are almost equal)
end
% --- convergence ---
inB = abs(P1 - M.Pss) <= cfg.bandConv*max(M.Pss, 1);
w   = max(1, round(cfg.holdConv/dt));
cs  = cumsum([0; inB(:)]);
k0  = find(ie, 1);  kc = NaN;
for k = k0:(numel(t)-w)
    if cs(k+w) - cs(k) == w, kc = k; break; end
end
if ~isnan(kc)
    M.tconv = t(kc) - tStart;
    M.nEval = nnz(abs(diff(Vr(k0:kc))) > 1e-3);
end
% --- time in unfavourable states ---
M.tUnreach = dt*nnz(ie & (Vr < Vb1/0.95 - 1.5));
M.tDsat    = dt*nnz(ie & (Dl >= 0.949));
M.tOver    = dt*nnz(ie & (Ib50 > cfg.overI*cfg.Ilim));
M.Qover    = trapz(t(ie), max(Ib50(ie) - cfg.Ilim, 0));   % charge above Ilim (A.s)
if any(~isnan(md))
    lim = md(ie) > 0;  M.nLim = nnz(diff([false; lim(:)]) == 1);
end
if any(~isnan(ph))
    % nRescan = number of global searches: transitions from the hold (6) to the first evaluation (1) after tStart
    pe = ph(ie);  pe = pe(:);
    M.nRescan = nnz(pe(2:end) == 1 & pe(1:end-1) == 6);
    if isStep                      % step: FIRST entry into the hold after tsw (ignores a hold left over from before the step)
        k6 = find(ie & ph == 6 & [false; ph(1:end-1) ~= 6], 1);
    else
        k6 = find(ie & ph == 6, 1);
    end
    if ~isempty(k6)
        M.tHold  = t(k6) - tStart;
        ia = t >= t(k6);  M.fAudit = nnz(ia & ph ~= 6)/nnz(ia);
    end
end
M.ok = true;
Tr = struct('id', j.id, 't', single(dec(t,cfg.decim)), 'V', single(dec(V,cfg.decim)), ...
            'P', single(dec(P,cfg.decim)), 'Vref', single(dec(Vr,cfg.decim)), ...
            'Vb', single(dec(Vb,cfg.decim)), 'Ib', single(dec(Ib,cfg.decim)), ...
            'D', single(dec(Dl,cfg.decim)), 'phase', single(decl(ph,cfg.decim)), ...
            'mode', single(decl(md,cfg.decim)), 'SOC', single(decl(so,cfg.decim)));
end

% =========================================================================
function [t, x] = getsig(o, name)
s = o.get(name);
if isa(s,'timeseries')
    t = s.Time(:);  x = squeeze(double(s.Data));  x = x(:);
elseif isstruct(s) && isfield(s,'signals')
    t = s.time(:);  x = double(s.signals.values(:));
else
    error('Signal "%s" not found in the simulation output.', name);
end
end

function x = sig_on(o, name, t)
[ts, xs] = getsig(o, name);
if numel(ts) == numel(t) && all(ts == t), x = xs;
else, x = interp1(ts, xs, t, 'previous', 'extrap'); end
end

function x = opt_on(o, name, t)
try, x = sig_on(o, name, t); catch, x = NaN(size(t)); end
end

function y = dec(x, d)          % mean of each block of d samples
n = floor(numel(x)/d)*d;  y = mean(reshape(x(1:n), d, []), 1)';
end

function y = decl(x, d)         % last sample of each block (discrete signals)
n = floor(numel(x)/d)*d;  y = x(d:d:n);
end

function k = refkey(pack, caseName, soc)
k = sprintf('%d|%s|S%02d', pack, caseName, round(100*soc));
end

% =========================================================================
function k = cc_start(t, Ib1, ie, cfg)
% index of the first sample from which |Ib1 - Ilim| <= tolI*Ilim holds continuously for >= 10 ms (ignores the start-up crossing)
k = [];
ok = ie(:) & abs(Ib1(:) - cfg.Ilim) <= cfg.tolI*cfg.Ilim;
w  = max(1, round(10e-3/median(diff(t))));
cs = cumsum([0; ok]);
i  = find(cs(1+w:end) - cs(1:end-w) == w, 1);
if ~isempty(i), k = i; end
end
