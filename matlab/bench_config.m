function cfg = bench_config()
% =========================================================================
% BENCH_CONFIG  Test configuration of the switched-model study.
%   The RA layer is applied to three published GMPPT cores:
%     DPSO (Ishaque & Salam, IEEE TIE 2013), CS (Ahmed & Salam, Applied Energy 2014)
%     and SCAN (Patel & Agarwal, IEEE TIE 2008); see Charger_Bench.m for alg and abl.
% =========================================================================
% Scenario = (battery pack, irradiance case, initial SOC). The reference P_reach is taken from a
% quasi-static P-V SWEEP in the same circuit (group 'sweep') -> run run_bench('sweep') FIRST.
% Case names vs. the paper: C1-C4 = Cases 1-4, C6 = Case 5, C7 = Case 6 (50 degC), C8 = Case 7 (65 degC);
% C5 is not used in the paper.
% =========================================================================
cfg.model   = 'Charger_Buck';
cfg.outDir  = 'results';
cfg.Tc      = 5e-6;       % controller sampling period (= ctrl_ts() in Charger_Bench.m)
% EVERY run uses the SAME solver settings. MaxStep 2e-6 s (1e-6 s is ~1.8 times slower).
%   A run that fails for numerical reasons is re-run with 1e-6 s, TOGETHER with all runs of the same
%   scenario and seed, so that paired comparisons use the same numerical settings.
cfg.MaxStep      = '2e-6';
cfg.MaxStepRetry = '1e-6';
cfg.nChunk  = 16;
cfg.Ilim    = 10;         % A
cfg.ETA     = 0.96;

% ------------------------- IRRADIANCE CASES (W/m2 per module, degC, GMPP power W) ----
C = [ ...
  mkc('C1', [1000 1000 1000 1000], 25, 919.2) ;
  mkc('C2', [1000 1000  600  600], 25, 601.6) ;
  mkc('C3', [1000  700  700  300], 25, 503.9) ;
  mkc('C4', [1000  850  500  300], 25, 400.0) ;
  mkc('C5', [1000  600  600  400], 25, 435.8) ;
  mkc('C6', [ 900  700  400  200], 25, 334.7) ;
  mkc('C7', [1000  850  500  300], 50, 358.1) ;
  mkc('C8', [1000  850  500  300], 65, 331.8) ];
cfg.cases = C;

% ------------------------- SCENARIOS --------------------------------------
cfg.static = [ ...
  sc(48,'C4','C4',0.60,1,0.60,Inf,0.8) ;
  sc(48,'C7','C7',0.20,1,0.20,Inf,0.8) ;
  sc(48,'C7','C7',0.60,1,0.60,Inf,0.8) ;
  sc(48,'C8','C8',0.60,1,0.60,Inf,0.8) ;
  sc(48,'C2','C2',0.60,1,0.60,Inf,0.8) ;
  sc(48,'C3','C3',0.20,1,0.20,Inf,0.8) ;   % CC
  sc(72,'C4','C4',0.60,1,0.60,Inf,0.8) ;
  sc(72,'C6','C6',0.60,1,0.60,Inf,0.8) ;
  sc(72,'C7','C7',0.60,1,0.60,Inf,0.8) ;
  sc(72,'C7','C7',0.95,1,0.95,Inf,0.8) ;
  sc(72,'C8','C8',0.20,1,0.20,Inf,0.8) ;
  sc(72,'C8','C8',0.60,1,0.60,Inf,0.8) ;
  sc(48,'C1','C1',0.60,1,0.60,Inf,0.8) ;   % CC
  sc(72,'C1','C1',0.60,1,0.60,Inf,0.8) ];  % CC
cfg.dynamic = [ ...
  sc(72,'C4','C1',0.60,1,0.60,0.5,1.1) ;
  sc(72,'C1','C4',0.60,1,0.60,0.5,1.1) ;
  sc(72,'C2','C4',0.60,1,0.60,0.5,1.1) ;
  sc(72,'C4','C2',0.60,1,0.60,0.5,1.1) ];
% Tend 1.12 s: SOC 0.60 -> 0.88 = SOC of the reference sweep
cfg.vrise = sc(72,'C8','C8',0.60,1.8e-4,0.88,Inf,1.12);
% vrise2: 72 V, C7, SOC 0.60 -> 0.90. The 85.1 V peak is REACHABLE at first and lies INSIDE the range (> 2 V from
%   the bound); the bound passes the peak at SOC ~0.80 -> an interior held point becomes unreachable (tests HT).
%   Qscale 1.9e-4: charging current ~4.2 A -> SOC rises by ~0.30 in ~1 s. Reference: sweep at SOC 0.90 (72|C7|S90).
cfg.vrise2 = sc(72,'C7','C7',0.60,1.9e-4,0.90,Inf,1.12);
cfg.cv = sc(48,'C1','C1',0.93,1e-3,0.99,Inf,0.8);

% ------------------------- CONTROLLERS ------------------------------------
cfg.algs = [ ...
  struct('name','P&O',     'alg',1,'abl',0) ;
  struct('name','DPSO',    'alg',5,'abl',0) ;   % conventional limiter (jumps towards Voc first)
  struct('name','CS',      'alg',7,'abl',0) ;
  struct('name','RA-DPSO', 'alg',6,'abl',0) ;   % proposed: RA layer + left-flank limiter + soft start
  struct('name','RA-CS',   'alg',8,'abl',0) ];
cfg.detAlgs = {'P&O','DPSO','RA-DPSO','RA-DPSO-R1','RA-DPSO-R1HT','RA-DPSO-noHT','RA-DPSO-noRA','RA-DPSO-convLim'};
cfg.excludeIds = {};

% ---- Ablation of the RA layer: components switched off, both cores ----
cfg.ablRA = [ ...
  struct('name','RA-DPSO-R1',   'alg',6,'abl',192) ;   % R1 only (R2 and HT off)
  struct('name','RA-DPSO-R1HT',  'alg',6,'abl',64)  ;   % R1 + HT (R2 off)
  struct('name','RA-DPSO-noRA', 'alg',6,'abl',16)  ;   % no RA (left-flank limiter kept)
  struct('name','RA-CS-R1',     'alg',8,'abl',192) ;
  struct('name','RA-CS-R1HT',    'alg',8,'abl',64)  ;
  struct('name','RA-CS-noRA',   'alg',8,'abl',16)  ];
cfg.ablRAScen = {'48-C7-S60','48-C8-S60','72-C4-S60','72-C7-S60','72-C7-S95','72-C8-S60'};
% ---- Ablation of the limiter: conventional limiter (RA kept) ----
cfg.ablLim = [ ...
  struct('name','RA-DPSO-convLim','alg',6,'abl',32) ;
  struct('name','RA-CS-convLim',  'alg',8,'abl',32) ];
cfg.ablLimScen = {'48-C1-S60','72-C1-S60','48-C3-S20'};

cfg.dynAlgs   = {'P&O','DPSO','CS','RA-DPSO','RA-CS'};
cfg.vriseAlgs = {'DPSO','CS','RA-DPSO','RA-CS','RA-DPSO-R1','RA-DPSO-R1HT','RA-CS-R1','RA-CS-R1HT'};
% Variant R1 + R2 (HT off, abl 128) - used for vrise2 only, to isolate the role of HT
cfg.ablNoHT = [ ...
  struct('name','RA-DPSO-noHT','alg',6,'abl',128) ;
  struct('name','RA-CS-noHT',  'alg',8,'abl',128) ];
cfg.vrise2Algs = [cfg.vriseAlgs, {'RA-DPSO-noHT','RA-CS-noHT'}];
cfg.cvAlgs    = {'P&O','DPSO','CS','RA-DPSO','RA-CS'};

% ------------------------- NUMBER OF RUNS ---------------------------------
cfg.nRun = struct('static',10,'ablRA',5,'ablLim',5,'dynamic',10,'vrise',10,'vrise2',10,'cv',10);
cfg.Tsweep = 0.75;

% ======================= ADDITIONAL GROUPS ========================================
% abl: 1 PI right-flank limiter | 2 fast left-flank limiter | 16 RA off | 32 conventional limiter |
%   64 R2 off | 128 HT off | 256 adaptive search range | 512 soft start of the right-flank limiter | 1024 20 kHz sampling |
%   4096*k: row k of the sensitivity table
mpptS = {'48-C4-S60','48-C7-S20','48-C7-S60','48-C8-S60','72-C4-S60','72-C6-S60','72-C7-S60','72-C7-S95','72-C8-S20','72-C8-S60'};
ccS   = {'48-C1-S60','72-C1-S60','48-C3-S20','48-C2-S60'};
heavy = [ sc(48,'C4','C4',0.60,1,0.60,Inf,0.8,6) ;      % CC under heavy shading: low Ilim
          sc(72,'C4','C4',0.60,1,0.60,Inf,0.8,4) ;
          sc(72,'C6','C6',0.60,1,0.60,Inf,0.8,3.5) ];
SM = pick(cfg, mpptS);  SC = [pick(cfg, ccS); heavy];
% adaptive search range (no special handling of the bound, no R2) and the third core SCAN
cfg.adaptRange = [ grp([al('DPSO-AR',5,256); al('SCAN',9,0); al('SCAN-AR',9,256); al('RA-SCAN',10,0)], SM, 1) ;
              grp(al('CS-AR',7,256), SM, 5) ];
cfg.vriseAR = [ grp([al('SCAN',9,0); al('RA-SCAN',10,0); al('RA-SCAN-R1',10,192); al('DPSO-AR',5,256)], cfg.vrise, 1) ;
              grp(al('CS-AR',7,256), cfg.vrise, 5) ];
% limiters - right-flank soft start, PI min-select, fast left flank; CC under heavy shading
cfg.limiter = [ grp([al('DPSO-SS',5,512); al('DPSO-PI',5,1); al('DPSO-PI-SS',5,513); ...
                    al('RA-DPSO-convSS',6,544); al('RA-DPSO-F',6,2)], SC, 1) ;
               grp([al('DPSO',5,0); al('RA-DPSO',6,0)], heavy, 1) ;
               grp([al('CS',7,0); al('RA-CS',8,0)], heavy, 3) ;
               grp(al('RA-CS-F',8,2), [pick(cfg, {'48-C3-S20','48-C2-S60'}); heavy], 3) ;
               grp(al('RA-CS-F',8,2), pick(cfg, {'48-C1-S60','72-C1-S60'}), 1) ];   % uniform irradiance: the stochastic core is not active
cfg.dynLimiter = [ grp([al('RA-DPSO-F',6,2); al('DPSO-PI-SS',5,513)], cfg.dynamic, 1) ;
               grp(al('RA-CS-F',8,2), cfg.dynamic, 5) ;
               grp([al('RA-DPSO-F',6,2); al('DPSO-PI-SS',5,513); al('DPSO-SS',5,512)], cfg.cv, 1) ;
               grp(al('RA-CS-F',8,2), cfg.cv, 1) ];
% sensitivity (RA-DPSO, deterministic)
% switched model: 5 main rows (all 10 rows were run with the averaged model)
sensN = {'Dmax0.93','Dmax0.97','Voff+1','Voff-1','dVc0','dVc1.5','dVR1','dVR4','tauR2-50ms','R2-OCV'};
sensK = [1 2 3 7 10];
AS = struct('name',{},'alg',{},'abl',{});
for k = sensK, AS(end+1,1) = al(['RA-DPSO-' sensN{k}], 6, 4096*k); end %#ok<AGROW>
cfg.sens = [ grp(AS, pick(cfg, {'72-C4-S60','72-C7-S95','72-C8-S60','48-C7-S60'}), 1) ;
                grp(AS, cfg.vrise, 1) ];
% 20 kHz sampling
cfg.fs20k = grp([al('DPSO-20k',5,1024); al('RA-DPSO-20k',6,1024)], [SM; pick(cfg, {'48-C1-S60','72-C1-S60'})], 1);

% ------------------------- EVALUATION CRITERIA ----------------------------
cfg.winSS    = 0.05;
cfg.bandConv = 0.02;
cfg.holdConv = 0.02;
cfg.effMin   = 0.99;
cfg.tolI     = 0.02;
cfg.overI    = 1.05;
cfg.decim    = 5;
end

% =========================================================================
function c = mkc(name, G, T, Pg)
c = struct('name',name,'G',G,'T',T*ones(1,4),'Pg',Pg);
end

function s = sc(pack, cA, cB, soc0, Qs, socRef, tsw, Tend, Ilim)
% Ilim (optional, default 10 A): CC scenario with a low current limit (heavy shading)
if nargin < 9, Ilim = 10; end
if strcmp(cA,cB), nm = sprintf('%d-%s-S%02d', pack, cA, round(100*soc0));
else,             nm = sprintf('%d-%s>%s-S%02d', pack, cA, cB, round(100*soc0)); end
if Qs < 1, nm = [nm sprintf('-Q%g', Qs)]; end
if Ilim ~= 10, nm = [nm sprintf('-I%g', Ilim)]; end
s = struct('name',nm,'pack',pack,'caseA',cA,'caseB',cB,'soc0',soc0,'Qscale',Qs, ...
           'socRef',socRef,'tsw',tsw,'Tend',Tend,'Ilim',Ilim);
end

function a = al(name, alg, abl)
a = struct('name',name,'alg',alg,'abl',abl);
end

function g = grp(alg, scen, nRun)
g = struct('alg',alg,'scen',scen,'nRun',nRun);
end

function S = pick(cfg, names)
S = cfg.static(ismember({cfg.static.name}, names));
end
