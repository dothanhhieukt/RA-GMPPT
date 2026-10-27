function [D, phaseOut, modeOut, VrefOut, VminOut, rstCount, duelCount] = Charger_Bench(Vpv,Ipv,Vbat,Ibat,seed0,alg,abl,Ilim,Vcv)
%#codegen
% =========================================================================
% CHARGER_BENCH  Controller of the buck-type PV-battery charger (MATLAB Function block)
%   Reachability-aware (RA) layer for global MPPT with CC/CV charging, as used in
%   "A Reachability-Aware Layer for Global MPPT in Buck-Type PV-Battery Chargers Under
%   Partial Shading: A Simulation and Processor-in-the-Loop Study".
%
% Inputs : Vpv, Ipv (PV voltage/current), Vbat, Ibat (battery voltage/current), seed0 (run seed:
%          the CS core seeds the random stream with rng(seed0*1000 + 7) at its first call),
%          alg, abl (see below), Ilim (CC limit, A), Vcv (CV voltage, V)
% Outputs: D (duty cycle), phaseOut (controller phase), modeOut (0 MPPT, 1 left-flank limiter,
%          2 right-flank limiter), VrefOut (PV-voltage reference), VminOut (reachability bound),
%          rstCount (number of global searches), duelCount (always 0, kept for the model ports)
%
% alg: -2 slow open-loop PV-voltage ramp 150 V -> 20 V (P-V sweep of the reference curve)
%      -1 / 0 open-loop duty-cycle / voltage steps (model checks)
%       1 P&O                    5 DPSO [Ishaque & Salam 2013]   6 RA-DPSO
%       7 CS [Ahmed & Salam 2014] 8 RA-CS                         9 SCAN [Patel & Agarwal 2008]
%      10 RA-SCAN
% abl (bitmask):
%       1 right-flank limiter as PI min-select (industrial variant; conventional limiter only)
%       2 fast left-flank limiter (extension F): enters the limiter 0.3 ms after the current
%         reaches 0.95*Ilim while the core is holding; no rate limit when reducing power;
%         early soft start (2 ms); CV loop 6 times faster
%       8 DPSO restart rule: first 10 % comparison after convergence skipped (modified rule)
%      16 RA layer off             32 conventional (right-flank) limiter
%      64 R2 off                  128 hold tracking (HT) off: a held point at the bound does not
%                                     follow the bound, and an interior held point that becomes
%                                     unreachable does not trigger a new search
%     256 adaptive search range for the core (alg 5/7/9): upper limit = Vbat/Dmax bound,
%         no special handling of the bound, no R2 (baselines DPSO-AR, CS-AR, SCAN-AR)
%     512 soft start of the right-flank limiter: warm-up at 150 V (~Voc), 19.5 ms into the limiter
%         from 150 V towards 119 V
%    1024 20 kHz sampling: the controller runs once per PWM period (50 us) on one sample taken
%         at mid on-time (t = D*T/2)
%   floor(abl/4096) = row of the sensitivity table [Dhat, Voff, dVc, dVR, tauR2, ocvR2]:
%      0 nominal; 1 Dhat 0.93; 2 Dhat 0.97; 3 Voff +1 V; 4 Voff -1 V; 5 dVc 0; 6 dVc 1.5 V;
%      7 dVR 1 V; 8 dVR 4 V; 9 R2 filter 50 ms; 10 R2 on the OCV estimate (Vbat - R*Ibat)
%   Dhat: Dmax used in the bound estimate (the true saturation of D stays at 0.95);
%   Voff: offset of the battery-voltage sensor (added to the measured Vbat).
% Coordinate: the cores work on x = 1 - Vref/200 (linear in the PV voltage).
% =========================================================================
% =========================================================================
persistent dInit sCnt lastP oD oPh oMd oVr oVm oRs oDu
ablU = uint32(abl);
dec20 = bitand(ablU, uint32(1024)) ~= 0;
rskip_flag(1, bitand(ablU, uint32(8)) ~= 0);
if isempty(dInit)
    dInit = true;  sCnt = int32(-1);  lastP = int32(-1);  oD = 0;  oPh = 0;  oMd = 0;  oVr = 0;  oVm = 0;  oRs = 0;  oDu = 0;
    if dec20, ctrl_ts(1, 50e-6); end
end
if dec20
    % one call per 50 us PWM period, sample taken at MID ON-TIME (t = D*T/2 from the start of the period),
    %        as with a PWM-triggered ADC: instantaneous value = average inductor current (CCM)
    sCnt = sCnt + 1;
    p  = idivide(sCnt, int32(10));            % index of the PWM period
    ph = sCnt - 10*p;                         % sample ph (0..9) within the period, 5 us step
    tg = int32(round(5*oD));                  % sample closest to D*T/2
    if p == lastP || ph < tg
        D = oD;  phaseOut = oPh;  modeOut = oMd;  VrefOut = oVr;  VminOut = oVm;  rstCount = oRs;  duelCount = oDu;
        return;
    end
    lastP = p;
end
[oD, oPh, oMd, oVr, oVm, oRs, oDu] = cb_main(Vpv, Ipv, Vbat, Ibat, seed0, alg, ablU, Ilim, Vcv);
D = oD;  phaseOut = oPh;  modeOut = oMd;  VrefOut = oVr;  VminOut = oVm;  rstCount = oRs;  duelCount = oDu;
end

% =========================================================================
function [D, phaseOut, modeOut, VrefOut, VminOut, rstCount, duelCount] = cb_main(Vpv,Ipv,Vbat,Ibat,seed0,alg,ablU,Ilim,Vcv)
VB    = 200;       % virtual bus (V), only for x <-> Vref: Vref = VB*(1 - x)
DMAXB = 0.95;      % maximum duty cycle of the buck (TRUE saturation)
HS    = 4;         % hold-interval factor: 650*HS us = 2.6 ms
TS = ctrl_ts(0, 0);  TAU = 5e-4;  aF = TS/TAU;
NPWM = max(1, round(50e-6/TS));  KP = 12;  KI = 1000;
LIM_PER = round(650e-6*HS/TS);
kL = 0.05;  kR = 0.01;  dxMax = 0.01;
ETA = 0.97;
SR  = 10e3;  dxSR = SR*NPWM*TS/VB;
kC  = 0.02;  NEX = 5;
T_ENT  = round(3e-3/TS);
T_FAST = max(1, round(0.3e-3/TS));          % extension F
KPL = 0.025;  KIL = 50;                     % PI right-flank limiter (x per unit error; 1/s)
% ---- sensitivity table (see header) ----
SENS = [0.95  0     0.5  2  5e-4  0 ;
        0.93  0     0.5  2  5e-4  0 ;
        0.97  0     0.5  2  5e-4  0 ;
        0.95  1     0.5  2  5e-4  0 ;
        0.95 -1     0.5  2  5e-4  0 ;
        0.95  0     0    2  5e-4  0 ;
        0.95  0     1.5  2  5e-4  0 ;
        0.95  0     0.5  1  5e-4  0 ;
        0.95  0     0.5  4  5e-4  0 ;
        0.95  0     0.5  2  50e-3 0 ;
        0.95  0     0.5  2  5e-4  1 ];
si = min(size(SENS,1), double(idivide(ablU, uint32(4096))) + 1);
Dhat = SENS(si,1);  Voff = SENS(si,2);  dVc = SENS(si,3);  dVR = SENS(si,4);  tauR2 = SENS(si,5);  ocvR2 = SENS(si,6) > 0;
Vbat = Vbat + Voff;                         % battery-voltage sensor offset

persistent init VbF IbF IpF lim xL xPk side cnt xLast olCnt poX poP poDir poCnt poW
persistent vAcc nAcc Dcmd integ restartCore lastPh lastRst lastDuel lastXh xTrim xCmd slewing xFFl
persistent Icv exCnt exhaust entCnt fromL Ip0 kT ssDone
persistent VbS fCnt piMode xI fInh
if isempty(init)
    init = true;  VbF = Vbat;  IbF = 0;  IpF = 0;  lim = false;  xL = 0.5;  xPk = 0.5;
    side = 1;  cnt = 0;  xLast = 0.4;  olCnt = 0;
    poX = 0.4;  poP = 0;  poDir = 1;  poCnt = 0;  poW = false;
    vAcc = 0;  nAcc = 0;  Dcmd = 0;  integ = 0;  restartCore = false;
    lastPh = 0;  lastRst = 0;  lastDuel = 0;  lastXh = 0.4;
    xTrim = 0;  xCmd = 0.4;  slewing = false;  xFFl = 0.4;
    Icv = 10;  exCnt = 0;  exhaust = false;  entCnt = 0;
    fromL = false;  Ip0 = 0;  kT = 0;  ssDone = false;
    VbS = Vbat;  fCnt = 0;  piMode = false;  xI = 0.4;  fInh = false;
end
VbF = VbF + aF*(Vbat - VbF);
IbF = IbF + aF*(Ibat - IbF);
IpF = IpF + aF*(Ipv - IpF);
VbS = VbS + (TS/tauR2)*(Vbat - VbS);
VpvMin  = max(VbF, 1)/Dhat + 1.0;           % +1 V margin (CC/CV limiter)
VpvMinC = max(VbF, 1)/Dhat - dVc;           % bound of the core (eq. 4)
xmaxC   = 1 - VpvMinC/VB;
xmax    = 1 - VpvMin/VB;
if Vcv > 70, Rhat = 0.15; else, Rhat = 0.0975; end
if ocvR2, VbR2 = VbF - Rhat*IbF; else, VbR2 = VbS; end
xR2 = 1 - (max(VbR2, 1)/Dhat - dVc)/VB;    % bound used by R2
propC  = (alg == 4) || (alg == 6) || (alg == 8) || (alg == 10);
rightLim = ~propC || bitand(ablU, uint32(32)) ~= 0;
raOn   = bitand(ablU, uint32(16)) == 0;
useR2  = raOn && bitand(ablU, uint32(64)) == 0;
useHT  = raOn && bitand(ablU, uint32(128)) == 0;
useAR  = bitand(ablU, uint32(256)) ~= 0;
piLim  = rightLim && bitand(ablU, uint32(1)) ~= 0;
fast   = ~rightLim && bitand(ablU, uint32(2)) ~= 0;
ssR    = rightLim && bitand(ablU, uint32(512)) ~= 0;
xWarm  = 0.4;  if ~rightLim && propC, xWarm = xmaxC; end
if ssR, xWarm = 1 - 150/VB; end
vAcc = vAcc + Vpv;  nAcc = nAcc + 1;

phaseOut = lastPh;  rstCount = lastRst;  duelCount = lastDuel;  modeOut = 0;
openCirc = false;

if alg == -1
    olCnt = olCnt + 1;
    seq = [0.60 0.65 0.55 0.70 0.50 0.80 0.45];
    kk = min(numel(seq), floor(olCnt*TS/20e-3) + 1);
    D = seq(kk);  VrefOut = 0;  VminOut = VpvMin;
    return;
end

if alg == -2
    olCnt = olCnt + 1;
    Vsw = 150 - 130*min(1, max(0, (olCnt*TS - 0.03)/0.7));
    x = 1 - Vsw/VB;
elseif alg == 0
    olCnt = olCnt + 1;
    seqV = [100 90 110 85 120 70 130];
    kk = min(numel(seqV), floor(olCnt*TS/20e-3) + 1);
    Vref = max(VpvMin + 2, seqV(kk));
    x = 1 - Vref/VB;
else
    IlimS = max(Ilim, 1e-3);
    e = max(IbF/IlimS - 1, (VbF - Vcv)/0.5);
    wasLim = lim;
    if e > 0.02, entCnt = entCnt + 1; else, entCnt = 0; end
    if lastPh == 1 || lastPh == 5, fInh = false; end            % core restarted its search -> fast entry allowed again
    if fast && ~fInh && lastPh == 6 && e > -0.05, fCnt = fCnt + 1; else, fCnt = 0; end
    kT = kT + 1;
    if fast, tSS = 2e-3; else, tSS = 19.5e-3; end
    if ~ssDone && kT >= round(tSS/TS) && ((~rightLim && propC) || ssR)
        ssDone = true;
        if ~lim
            lim = true;  cnt = 0;  xTrim = 0;  slewing = false;  exCnt = 0;  exhaust = false;  fromL = false;
            entCnt = 0;  Icv = IlimS;  xPk = 1 - 119/VB;
            if ssR                                   % right flank: from ~Voc towards 119 V
                side = -1;  xL = 1 - 150/VB;  xI = xL;  piMode = piLim;
            else
                side = 1;  xL = xmax;  xCmd = xmax;  xFFl = xmax;
            end
        end
    end
    if ~lim && (entCnt >= T_ENT || (fast && fCnt >= T_FAST))
        fastE = entCnt < T_ENT;
        entCnt = 0;  fCnt = 0;
        lim = true;  cnt = 0;  xTrim = 0;  slewing = false;  exCnt = 0;  exhaust = false;  fromL = false;
        if alg >= 2 && (lastPh == 3 || lastPh == 9 || lastPh == 6 || lastPh == 4)
            xPk = lastXh;
        else
            xPk = xLast;
        end
        xL = xLast;  xCmd = xLast;  xFFl = xPk;
        Icv = min(IlimS, max(IbF, 0.2));
        if fastE, Icv = IlimS;  xPk = min(xPk, 1 - 119/VB); end
        if rightLim
            if piLim                                 % PI min-select: no jump, starts at the held point
                side = -1;  piMode = true;  xI = xLast;  xPk = xLast;
            else
                side = -1;  xL = 1 - 150/VB;  fromL = true;  Ip0 = max(IpF, 1);
            end
        else
            side = 1;
        end
    end
    if lim
        Iset = max(0.2, min(IlimS, Icv));
        eI   = IbF/Iset - 1;
        cnt = cnt + 1;
        if side < 0 && piMode
            if nAcc >= NPWM                          % PI on the constraint error, every PWM period
                xI = min(xPk, max(0.15, xI - KIL*NPWM*TS*e));
                xL = min(xPk, max(0.15, xI - KPL*e));
                if xI >= xPk - 1e-9 && e < -0.02, lim = false;  piMode = false; end
            end
        elseif cnt >= LIM_PER
            cnt = 0;
            if side > 0
                if exhaust, exCnt = exCnt + 1; else, exCnt = 0; end
                if fast, kCe = 6*kC; else, kCe = kC; end
                Icv = min(IlimS, max(0, Icv - IlimS*kCe*(VbF - Vcv)/0.5));
                xUn = xFFl + xTrim;
                if ~slewing && ~(xUn >= xmax && eI > 0) && ~(xUn <= xPk && eI < 0)
                    xTrim = min(1, max(-1, xTrim + min(dxMax, max(-dxMax, kL*eI))));
                end
            else
                if fromL && IpF < 0.15*Ip0 && e < -0.1
                    xL = xL + 2*dxMax;
                else
                    xL = xL - min(dxMax, max(-dxMax, kR*e));
                end
                xL = min(xmax, max(xL, 0.15));
                if xL >= xPk && e < -0.02, lim = false; end
            end
        end
        if lim && side > 0 && nAcc >= NPWM
            if IpF > 0.5
                Vt  = VbF*Iset/(ETA*IpF);
                xFF = 1 - min(VB*(1 - xPk), max(VpvMin, Vt))/VB;
            else
                xFF = xPk;
            end
            xFFl = xFF;
            xDes = min(xmax, max(xPk, xFF + xTrim));
            if fast && xDes > xCmd                   % power-reducing direction: no rate limit
                slewing = false;  xCmd = xDes;
            else
                slewing = abs(xDes - xCmd) > dxSR;
                xCmd = xCmd + min(dxSR, max(-dxSR, xDes - xCmd));
            end
            exhaust = ~slewing && xCmd >= xmax - 1e-6 && e > 0.02;
            cvNoLeft = VbF >= Vcv - 0.1 && Icv < IlimS - 1e-6 && IpF > 0.5 && ...
                       VpvMin*IpF*ETA > 1.05*VbF*Iset;
            if fast && xCmd >= xmax - 1e-6 && IpF > 0.5 && Vpv >= VpvMin - 1.5 && Vpv*IpF*ETA > 1.05*VbF*Iset
                cvNoLeft = true;                     % no left-flank solution (CC included)
            end
            if (exhaust && exCnt >= NEX) || cvNoLeft
                side = -1;  xL = 1 - 150/VB;
                fromL = true;  Ip0 = max(IpF, 1);
            elseif xDes <= xPk + 1e-6 && xCmd <= xPk + 1e-6 && eI < -0.02 && Icv >= IlimS - 1e-6
                lim = false;
            else
                xL = xCmd;
            end
        end
        if wasLim && ~lim
            restartCore = true;  poCnt = 0;  poP = -1e9;  fInh = true;
        end
        x = xL;
        phaseOut = 20;
        if side > 0, modeOut = 1; else, modeOut = 2; end
    else
        if alg == 1
            [x, poX, poP, poDir, poCnt, poW] = po_step(Vpv*Ipv, poX, poP, poDir, poCnt, poW, xmax, HS);
            ph = 0;  rs = 0;  du = 0;  xh = x;
        elseif alg == 5 || alg == 6
            [x, ph, rs, xh] = dpso_core(Vpv, Ipv, alg == 6 && raOn, xmaxC, HS, restartCore, useR2, useHT, xWarm, useAR && alg == 5, dVR, xR2);
            du = 0;
        elseif alg == 7 || alg == 8
            [x, ph, rs, xh] = cs_core(Vpv, Ipv, seed0, alg == 8 && raOn, xmaxC, HS, restartCore, useR2, useHT, xWarm, useAR && alg == 7, dVR, xR2);
            du = 0;
        elseif alg == 9 || alg == 10
            [x, ph, rs, xh] = scan_core(Vpv, Ipv, alg == 10 && raOn, xmaxC, HS, restartCore, useR2, xWarm, useAR && alg == 9, dVR, xR2);
            du = 0;
        else
            x = xLast;  ph = 0;  rs = 0;  du = 0;  xh = x;   % unknown alg: hold the last point
        end
        restartCore = false;
        xLast = x;  lastPh = ph;  lastRst = rs;  lastDuel = du;  lastXh = xh;
        phaseOut = ph;  rstCount = rs;  duelCount = du;
        openCirc = (ph == 8);
    end
end

% ---------------- PV-voltage loop (updated every PWM period) ----------------
Vref = VB*(1 - x);
if nAcc >= NPWM
    vAvg = vAcc/nAcc;  vAcc = 0;  nAcc = 0;
    if openCirc
        Dcmd = 0;  integ = 0;
    else
        Dff  = VbF/max(Vref, 1);
        ev   = (vAvg - Vref)/max(Vref, 1);
        Dnew = Dff*(1 + KP*ev + integ);
        Dsat = min(DMAXB, max(0, Dnew));
        if (Dsat == Dnew) || (Dnew > DMAXB && ev < 0) || (Dnew < 0 && ev > 0)
            integ = integ + KI*NPWM*TS*ev;
        end
        Dcmd = Dsat;
    end
end
if openCirc, Dcmd = 0; end
D = Dcmd;
VrefOut = Vref;  VminOut = VpvMin;
end

% -------------------------------------------------------------------------
% classic P&O on Vref (step 0.005 in x ~ 1 V), updated every 650*HS samples
function [x, poX, poP, poDir, poCnt, poW] = po_step(P, poX, poP, poDir, poCnt, poW, xmax, HS)
poCnt = poCnt + 1;
if ~poW
    if poCnt >= round(20e-3/ctrl_ts(0, 0)), poW = true; poCnt = 0; poP = P; end   % WARMUP 20 ms
else
    if poCnt >= round(650e-6*HS/ctrl_ts(0, 0))
        poCnt = 0;
        if P < poP, poDir = -poDir; end
        poP = P;
        poX = min(xmax, max(0.15, poX + poDir*0.005));
    end
end
x = poX;
end

% ---------- SIGNED Levy step (Mantegna), used by CS ----------
function s = levy_signed()
beta = 1.5;
sigmau = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
uu = sigmau*randn(1);  vv = randn(1);
s = uu/(abs(vv)^(1/beta));
end

% =========================================================================
% DPSO - Deterministic PSO (K. Ishaque, Z. Salam, IEEE TIE 60(8), 2013, 3195-3206)
%   v = w*v + (Pbest - x) + (Gbest - x),  |v| <= VMAX,  NP = 3, w = 0.4, v0 = 0  (eq. 6)
%   Shading test (eqs. 11, 12) on the first 3 points: (I3 - I2)/I3 >= 0.1 and (V1 - V2)/V1 >= 0.2
%   (V1 > V2 > V3) -> shaded: global mode (DPSO); otherwise: local mode at once.
%   Stop: two particles closer than 0.45*Voc_module (paper: 30-60 %) and particle velocity < VSTOP,
%   or at most IT_MAX iterations. Local mode: variable-step hill climbing (paper: "variable step-size").
%   Restart when |P - P_previous window|/P > 10 % (window ~50 ms, as Pthr of the paper).
% Adapted to the buck + PV-voltage loop: particle = x (Vref = 200*(1 - x), linear in Vpv)
%   instead of the duty cycle of the buck-boost; VMAX set by the rule of the paper (voltage step
%   not larger than the peak spacing 0.8*Voc_module): 0.075 in x = 15 V = 0.5*0.8*Voc_module.
%   Hold interval per evaluation 2.6 ms (5.2 ms when |dx| > 0.05), P = mean of the last 200 us.
% RA (alg 6): upper limit of the search range and of every command = xmaxC (bound Vpv_min from Vbat);
%   the 3 initial particles are spread up to the BOUND (third particle at the bound = candidate standing for
%   the unreachable peaks); a held point in local mode follows the bound when Vbat changes. Deterministic.
% =========================================================================
function [x, phaseOut, rstCount, xHold] = dpso_core(Vpv, Ipv, useRA, xmaxDyn, HS, restart, useR2, useHT, xWarm, useAR, dVR, xR2)
NP = 3;  W = 0.4;  VMAX = 0.075;  VSTOP = 0.01;  DSPREAD = 0.45*37.2/200;  IT_MAX = 20;
XLO = 0.33;  XHI0 = 0.90;                       % Vref <= 134 V (0.9*Voc STC) ... >= 20 V
TSC = ctrl_ts(0, 0);
SETTLE = round(650e-6*HS/TSC);  LONG_SETTLE = round(1300e-6*HS/TSC);
WARMUP = round(20e-3/TSC);  NAVG = round(200e-6/TSC);
HC_S0 = 0.01;  HC_SMIN = 0.0025;  NWIN = 20;  PTHR = 0.10;
DXR = dVR/200;                                  % R2 bound-drift threshold (2 V nominal; sensitivity table)
persistent init st cnt hold xc X V Pb PbX Gb GbX k it Pm Vm Im pAcc vAcc iAcc nAcc
persistent hcS hcDir hcP Pw nw rst warmed xhiS atB lostI
xhi = XHI0;  if useRA || useAR, xhi = max(XLO + 0.05, min(XHI0, xmaxDyn)); end   % useAR: adaptive search range
xr = max(XLO + 0.05, min(XHI0, xR2));           % bound used by R2
if isempty(init)
    init = true;  st = 0;  cnt = 0;  hold = 0;  xc = 0.4;  warmed = false;
    X = zeros(NP,1);  V = zeros(NP,1);  Pb = -ones(NP,1);  PbX = zeros(NP,1);  Gb = -1;  GbX = 0.4;
    k = 1;  it = 0;  Pm = zeros(NP,1);  Vm = zeros(NP,1);  Im = zeros(NP,1);
    pAcc = 0;  vAcc = 0;  iAcc = 0;  nAcc = 0;
    hcS = HC_S0;  hcDir = 1;  hcP = 0;  Pw = -1;  nw = 0;  rst = 0;
    xhiS = xr;  atB = false;  lostI = false;
end
if restart, cnt = 0;  hold = min(hold, SETTLE);  pAcc = 0;  vAcc = 0;  iAcc = 0;  nAcc = 0;  end   % WARMUP not repeated
% ---------------- hold interval ----------------
if cnt < hold
    cnt = cnt + 1;
    if cnt > hold - NAVG, pAcc = pAcc + Vpv*Ipv;  vAcc = vAcc + Vpv;  iAcc = iAcc + Ipv;  nAcc = nAcc + 1; end
    if useRA && st == 3                                  % RA during the hold interval
        if atB && useHT, xc = xhi;                       % HT: a held point at the bound follows the bound
        else
            if ~atB && xc > xhi + 1e-6, lostI = true; end   % interior held point has just become unreachable
            xc = min(xc, xhi);                           % R1: projection
        end
    elseif useAR && st == 3
        xc = min(xc, xhi);                               % adaptive range: projection only
    end
    x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
    return;
end
P = Vpv*Ipv;  Vme = Vpv;  Ime = Ipv;
if nAcc > 0, P = pAcc/nAcc;  Vme = vAcc/nAcc;  Ime = iAcc/nAcc; end
pAcc = 0;  vAcc = 0;  iAcc = 0;  nAcc = 0;
xn = xc;  goGlobal = false;
switch st
    case 0                                        % WARMUP, then start of the search
        if ~warmed, warmed = true;  xc = xWarm;  cnt = 1;  hold = WARMUP;   % soft start
            x = xc;  phaseOut = 0;  rstCount = rst;  xHold = xc;  return;
        end
        goGlobal = true;
    case 1                                        % evaluate the first 3 points (initial Pbest)
        Pm(k) = P;  Vm(k) = Vme;  Im(k) = Ime;
        Pb(k) = P;  PbX(k) = X(k);
        if P > Gb, Gb = P;  GbX = X(k); end
        if k < NP
            k = k + 1;  xn = X(k);
        else
            % eqs. (11)-(12): X(1) has the highest voltage (smallest x)
            shaded = (Im(3) - Im(2))/max(Im(3), 1e-3) >= 0.1 && (Vm(1) - Vm(2))/max(Vm(1), 1) >= 0.2;
            if shaded
                st = 2;  it = 1;  k = 1;
                [X(1), V(1)] = dpso_step(X(1), V(1), PbX(1), GbX, Gb, it, XLO, xhi, W, VMAX);  xn = X(1);
            else
                st = 3;  xn = GbX;  hcS = HC_S0;  hcDir = 1;  hcP = Gb;  Pw = Gb;  if rskip_flag(0, false), Pw = -1; end;  nw = 0;  xhiS = xr;   % abl 8: first comparison after convergence skipped
            end
        end
    case 2                                        % global mode (DPSO)
        if P > Pb(k), Pb(k) = P;  PbX(k) = X(k); end
        if P > Gb, Gb = P;  GbX = X(k); end
        if k < NP
            k = k + 1;  [X(k), V(k)] = dpso_step(X(k), V(k), PbX(k), GbX, Gb, it, XLO, xhi, W, VMAX);  xn = X(k);
        else
            it = it + 1;
            doStop = it >= IT_MAX;
            for i = 1:NP
                for j = 1:NP
                    if i ~= j && abs(X(i) - X(j)) < DSPREAD && abs(V(i)) < VSTOP, doStop = true; end
                end
            end
            if doStop
                st = 3;  xn = GbX;  hcS = HC_S0;  hcDir = 1;  hcP = Gb;  Pw = Gb;  if rskip_flag(0, false), Pw = -1; end;  nw = 0;  xhiS = xr;   % abl 8: first comparison after convergence skipped
            else
                k = 1;  [X(1), V(1)] = dpso_step(X(1), V(1), PbX(1), GbX, Gb, it, XLO, xhi, W, VMAX);  xn = X(1);
            end
        end
    otherwise                                     % local mode: variable-step hill climbing
        if useRA && useR2 && atB && abs(xr - xhiS) > DXR
            rst = rst + 1;  goGlobal = true;      % R2: holding at the bound and the bound has drifted > DXR -> new search
        elseif useRA && useHT && lostI
            rst = rst + 1;  goGlobal = true;      % HT: interior held point became unreachable -> new search at once
        end
        lostI = false;
        nw = nw + 1;
        if nw >= NWIN
            if Pw > 0 && abs(P - Pw)/Pw > PTHR
                rst = rst + 1;  goGlobal = true;
            end
            Pw = P;  nw = 0;
        end
        if ~goGlobal
            if P < hcP, hcDir = -hcDir;  hcS = max(HC_SMIN, hcS/2); end
            hcP = P;
            xn = min(xhi, max(XLO, xc + hcDir*hcS));
        end
end
if goGlobal                                       % (re)start of the global search: 3 evenly spread particles
    st = 1;  k = 1;  Gb = -1;  it = 0;  xhiS = xr;
    for i = 1:NP
        if useRA                                  % RA: last particle AT THE BOUND (unreachable peaks -> the bound is a candidate)
            X(i) = XLO + i/NP*(xhi - XLO);
        else
            X(i) = XLO + (i - 0.5)/NP*(xhi - XLO);
        end
        V(i) = 0;  Pb(i) = -1;
    end
    xn = X(1);
end
if useRA || useAR, xn = min(xn, xhi); end
atB = useRA && st == 3 && xn >= xhi - 1e-6;      % held point at the bound
hold = SETTLE;  if abs(xn - xc) > 0.05, hold = LONG_SETTLE; end
xc = xn;  cnt = 1;
x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;

end

% =========================================================================
% CS - Cuckoo Search MPPT (J. Ahmed, Z. Salam, Applied Energy 119 (2014) 118-130)
%   Each sample = one value of x (Vref = 200*(1 - x)); x_best is held after convergence (no P&O).
%   Following Applied Energy 119 (2014): eqs. (11)-(14) step s = kappa*(u/|v|^(1/beta))*(V_best - V_i), kappa = 0.05;
%   a worse new sample is discarded (greedy); apart from the best sample, every sample is abandoned with
%   probability Pa and replaced by a random one; repeated until the samples converge (here: spread < DCONV ~ 2 V,
%   at most GMAX generations); random initialization in [0.1 Voc, 0.9 Voc] (eq. 17); re-initialized when |dP|/P > 0.1 (eq. 15).
% RA (useRA): upper limit of the range = bound (xmaxDyn); the bound is measured first as an ADDITIONAL
%   candidate that does not take one of the 3 nests (R1), so the core keeps its exploration; R2 and HT as in dpso_core.
% =========================================================================
function [x, phaseOut, rstCount, xHold] = cs_core(Vpv, Ipv, seed0, useRA, xmaxDyn, HS, restart, useR2, useHT, xWarm, useAR, dVR, xR2)
NN = 3;  KAPPA = 0.05;  PA = 0.25;  DCONV = 0.01;  GMAX = 25;   % Ahmed & Salam (2014), Section 7.1: n = 3, beta = 1.5, Pa = 0.25, kappa = 0.05
XLO = 0.33;  XHI0 = 0.90;                       % Vref <= 134 V ... >= 20 V (as DPSO)
TSC = ctrl_ts(0, 0);
SETTLE = round(650e-6*HS/TSC);  LONG_SETTLE = round(1300e-6*HS/TSC);
WARMUP = round(20e-3/TSC);  NAVG = round(200e-6/TSC);
NWIN = 20;  PTHR = 0.10;  DXR = dVR/200;
persistent init st cnt hold xc X P Xn k g Gb GbX pAcc nAcc Pw nw rst warmed xhiS atB lostI mode kA dropI bndM
xhi = XHI0;  if useRA || useAR, xhi = max(XLO + 0.05, min(XHI0, xmaxDyn)); end   % adaptive search range
xr = max(XLO + 0.05, min(XHI0, xR2));
if isempty(init)
    init = true;  rng(uint32(mod(seed0*1000 + 7, 4294967295)));
    st = 0;  cnt = 0;  hold = 0;  xc = 0.4;  warmed = false;
    X = zeros(NN,1);  P = -ones(NN,1);  Xn = zeros(NN,1);  k = 1;  g = 0;  Gb = -1;  GbX = 0.4;
    pAcc = 0;  nAcc = 0;  Pw = -1;  nw = 0;  rst = 0;  xhiS = xr;  atB = false;  lostI = false;
    mode = 1;  kA = 1;  dropI = false(NN,1);  bndM = false;
end
if restart, cnt = 0;  hold = min(hold, SETTLE);  pAcc = 0;  nAcc = 0;  end   % WARMUP not repeated
% ---------------- hold interval ----------------
if cnt < hold
    cnt = cnt + 1;
    if cnt > hold - NAVG, pAcc = pAcc + Vpv*Ipv;  nAcc = nAcc + 1; end
    if useRA && st == 3
        if atB && useHT, xc = xhi;
        else
            if ~atB && xc > xhi + 1e-6, lostI = true; end
            xc = min(xc, xhi);
        end
    elseif useAR && st == 3
        xc = min(xc, xhi);
    end
    x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
    return;
end
Pm = Vpv*Ipv;  if nAcc > 0, Pm = pAcc/nAcc; end
pAcc = 0;  nAcc = 0;
xn = xc;  goGlobal = false;
switch st
    case 0                                        % WARMUP
        if ~warmed, warmed = true;  xc = xWarm;  cnt = 1;  hold = WARMUP;
            x = xc;  phaseOut = 0;  rstCount = rst;  xHold = xc;  return;
        end
        goGlobal = true;
    case 1                                        % evaluate the initial nests
        if bndM                                   % R1: the bound is an ADDITIONAL candidate (does not take a nest)
            bndM = false;  Gb = Pm;  GbX = xc;  xn = X(1);  k = 1;
            hold = SETTLE;  if abs(xn - xc) > 0.05, hold = LONG_SETTLE; end
            xc = xn;  cnt = 1;  x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
            return;
        end
        P(k) = Pm;
        if Pm > Gb, Gb = Pm;  GbX = X(k); end
        if k < NN
            k = k + 1;  xn = X(k);
        else
            st = 2;  g = 0;  k = 1;  mode = 1;
            Xn = X;
            for i = 1:NN                          % Levy step for every nest
                Xn(i) = min(xhi, max(XLO, X(i) + KAPPA*(GbX - X(i))*levy_signed()));
            end
            xn = Xn(1);
        end
    case 2                                        % CS generation
        endGen = false;
        if mode == 1                              % new nests from the Levy step: greedy acceptance
            if Pm > P(k), P(k) = Pm;  X(k) = Xn(k); end
            if Pm > Gb, Gb = Pm;  GbX = Xn(k); end
            k = k + 1;
            if k <= NN
                xn = Xn(k);
            elseif (max(X) - min(X) < DCONV) || g + 1 >= GMAX   % samples converged -> hold x_best
                st = 3;  xn = GbX;  Pw = Gb;  nw = 0;  xhiS = xr;   % R2 reference = start of the hold
            else                                  % [2014] apart from the best sample, each sample is abandoned with probability PA
                nA = 0;
                for i = 1:NN
                    dropI(i) = abs(X(i) - GbX) > 1e-9 && rand(1) < PA;
                    if dropI(i), nA = nA + 1;  occB = false;                                % R1: another nest already at the bound?
                        for j = 1:NN
                            if j ~= i
                                xj = X(j);  if j < i && dropI(j), xj = Xn(j); end
                                if xj >= xhi - 1e-9, occB = true; end
                            end
                        end
                        if useAR, Xn(i) = XLO + rand(1)*(xhi - XLO);   % sample within the reduced range
                        else,     Xn(i) = ra_draw(XLO, XHI0, xhi, occB); end
                    end
                end
                if nA > 0
                    kA = 1;  while ~dropI(kA), kA = kA + 1; end
                    mode = 2;  xn = Xn(kA);
                else
                    endGen = true;
                end
            end
        else                                      % evaluate the replacement samples (random positions)
            X(kA) = Xn(kA);  P(kA) = Pm;  dropI(kA) = false;
            if Pm > Gb, Gb = Pm;  GbX = X(kA); end
            if any(dropI)
                kA = 1;  while ~dropI(kA), kA = kA + 1; end
                xn = Xn(kA);
            else
                endGen = true;
            end
        end
        if endGen                                 % new generation: Levy step for every sample
            g = g + 1;
            if g >= GMAX
                st = 3;  xn = GbX;  Pw = Gb;  nw = 0;  xhiS = xr;   % R2 reference = start of the hold
            else
                mode = 1;  k = 1;
                for i = 1:NN
                    Xn(i) = min(xhi, max(XLO, X(i) + KAPPA*(GbX - X(i))*levy_signed()));
                end
                xn = Xn(1);
            end
        end
    otherwise                                     % hold x_best, monitor the power
        if useRA && useR2 && atB && abs(xr - xhiS) > DXR
            rst = rst + 1;  goGlobal = true;
        elseif useRA && useHT && lostI
            rst = rst + 1;  goGlobal = true;
        end
        lostI = false;
        nw = nw + 1;
        if ~goGlobal && nw >= NWIN
            if Pw > 0 && abs(Pm - Pw)/Pw > PTHR
                rst = rst + 1;  goGlobal = true;
            end
            Pw = Pm;  nw = 0;
        end
        if ~goGlobal, xn = xc; end
end
if goGlobal                                       % initialization: random nests (RA: the bound is measured first)
    st = 1;  k = 1;  Gb = -1;  g = 0;  mode = 1;  xhiS = xr;
    for i = 1:NN                                  % [eq. 17] random within the range
        occB = false;
        for j = 1:i-1, if X(j) >= xhi - 1e-9, occB = true; end, end
        if useAR, X(i) = XLO + rand(1)*(xhi - XLO);   % adaptive search range
        else,     X(i) = ra_draw(XLO, XHI0, xhi, occB); end
        P(i) = -1;   % eq. 17
    end
    xn = X(1);
    if useRA, bndM = true;  xn = xhi; end         % R1: measure the bound first as an additional candidate
end
if useRA || useAR, xn = min(xn, xhi); end
atB = useRA && st == 3 && xn >= xhi - DXR;      % HT: held point within DXR (~2 V) of the bound = HOLD AT THE BOUND
if atB, xn = xhi; end                            %   (exact comparison let a 0.5 V bound ripple set lostI and restart the search)
hold = SETTLE;  if abs(xn - xc) > 0.05, hold = LONG_SETTLE; end
xc = xn;  cnt = 1;
x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
end

% =========================================================================
% SCAN - third core of a different kind (deterministic, no population): measures at fixed points
%   V_k = k*0.8*Voc,module (k = 4, 3, 2, 1 -> 119.0 / 89.3 / 59.5 / 29.8 V; Patel & Agarwal, IEEE TIE 55(4), 2008)
%   then P&O (1 V step every 2.6 ms) around the best point; new scan when |dP|/P > 10 % within 20 hold intervals.
% RA (alg 10): R1 - scan points below the bound are projected onto it and MERGED into ONE candidate (the bound
%   is always measured), every P&O step is projected; R2 - P&O point within DXR of the bound and the bound has
%   drifted > DXR since the start of the hold -> new scan. useAR (alg 9 + abl 256): scan points below the bound dropped.
% =========================================================================
function [x, phaseOut, rstCount, xHold] = scan_core(Vpv, Ipv, useRA, xmaxDyn, HS, restart, useR2, xWarm, useAR, dVR, xR2)
XLO = 0.33;  XHI0 = 0.90;
TSC = ctrl_ts(0, 0);
SETTLE = round(650e-6*HS/TSC);  LONG_SETTLE = round(1300e-6*HS/TSC);
WARMUP = round(20e-3/TSC);  NAVG = max(1, round(200e-6/TSC));
NWIN = 20;  PTHR = 0.10;  DXR = dVR/200;  DPO = 0.005;
persistent init st cnt hold xc warmed Xs ns k Gb GbX pAcc nAcc Pw nw rst xhiS atB poP poDir
xhi = XHI0;  if useRA || useAR, xhi = max(XLO + 0.05, min(XHI0, xmaxDyn)); end
xr = max(XLO + 0.05, min(XHI0, xR2));
if isempty(init)
    init = true;  st = 0;  cnt = 0;  hold = 0;  xc = 0.4;  warmed = false;
    Xs = zeros(5,1);  ns = 0;  k = 1;  Gb = -1;  GbX = 0.4;
    pAcc = 0;  nAcc = 0;  Pw = -1;  nw = 0;  rst = 0;  xhiS = xr;  atB = false;
    poP = 0;  poDir = 1;
end
if restart, cnt = 0;  hold = min(hold, SETTLE);  pAcc = 0;  nAcc = 0;  end
if cnt < hold
    cnt = cnt + 1;
    if cnt > hold - NAVG, pAcc = pAcc + Vpv*Ipv;  nAcc = nAcc + 1; end
    if (useRA || useAR) && st == 3, xc = min(xc, xhi); end
    x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
    return;
end
Pm = Vpv*Ipv;  if nAcc > 0, Pm = pAcc/nAcc; end
pAcc = 0;  nAcc = 0;
xn = xc;  goGlobal = false;
switch st
    case 0
        if ~warmed, warmed = true;  xc = xWarm;  cnt = 1;  hold = WARMUP;
            x = xc;  phaseOut = 0;  rstCount = rst;  xHold = xc;  return;
        end
        goGlobal = true;
    case 1                                        % measure the scan points
        if Pm > Gb, Gb = Pm;  GbX = Xs(k); end
        if k < ns
            k = k + 1;  xn = Xs(k);
        else
            st = 3;  xn = GbX;  poP = Gb;  poDir = 1;  Pw = Gb;  nw = 0;  xhiS = xr;
        end
    otherwise                                     % P&O around the best point
        if useRA && useR2 && atB && abs(xr - xhiS) > DXR
            rst = rst + 1;  goGlobal = true;
        end
        nw = nw + 1;
        if ~goGlobal && nw >= NWIN
            if Pw > 0 && abs(Pm - Pw)/Pw > PTHR
                rst = rst + 1;  goGlobal = true;
            end
            Pw = Pm;  nw = 0;
        end
        if ~goGlobal
            if Pm < poP, poDir = -poDir; end
            poP = Pm;
            xn = min(xhi, max(XLO, xc + poDir*DPO));
        end
end
if goGlobal
    st = 1;  k = 1;  Gb = -1;  xhiS = xr;  ns = 0;
    for kk = 4:-1:1
        p = 1 - kk*0.8*37.2/200;                  % 119.0 / 89.3 / 59.5 / 29.8 V
        if useRA
            if p < xhi - DXR, ns = ns + 1;  Xs(ns) = p; end
        elseif useAR
            if p <= xhi, ns = ns + 1;  Xs(ns) = p; end
        else
            ns = ns + 1;  Xs(ns) = p;
        end
    end
    if useRA || ns == 0, ns = ns + 1;  Xs(ns) = xhi; end   % R1: the bound is one candidate (once)
    xn = Xs(1);
end
if useRA || useAR, xn = min(xn, xhi); end
nowB = useRA && st == 3 && xn >= xhi - DXR;
if nowB && ~atB, xhiS = xr; end                  % R2 reference = start of the hold at the bound
atB = nowB;
hold = SETTLE;  if abs(xn - xc) > 0.05, hold = LONG_SETTLE; end
xc = xn;  cnt = 1;
x = xc;  phaseOut = phase_of(st);  rstCount = rst;  xHold = xc;
end

function xs = ra_draw(XLO, XHI0, xhi, occB)
% R1 detail (ii): random CS sample drawn in the ORIGINAL range [XLO, XHI0] and PROJECTED onto the bound xhi.
% If the projected sample lands on the bound while another nest is already there -> redrawn in [XLO, xhi]
% (at most ONE nest at the bound: keeps diversity, avoids false convergence of several nests at the bound).
% Without RA: xhi = XHI0 -> eq. (17) of Ahmed & Salam (2014) unchanged.
xs = min(xhi, XLO + rand(1)*(XHI0 - XLO));
if xs >= xhi - 1e-9 && xhi < XHI0 - 1e-9
    if occB
        xs = XLO + rand(1)*(xhi - XLO);
    end
end
end

function [Xi, Vi] = dpso_step(Xi, Vi, PbXi, GbX, Gb, it, XLO, xhi, W, VMAX)
% eq. (6) of the paper + VMAX limit
vi = W*Vi + (PbXi - Xi);
if Gb >= 0, vi = vi + (GbX - Xi); end
vi = min(VMAX, max(-VMAX, vi));
if abs(vi) < 1e-4                                 % particle at Gbest: small perturbation (paper, Fig. 5b)
    vi = 0.01*(2*mod(it,2) - 1);
end
Vi = vi;
Xi = min(xhi, max(XLO, Xi + vi));
end

function ph = phase_of(st)
% 0 WARMUP, 1 initial evaluations, 5 global search, 6 local mode / hold
ph = 0;
if st == 1, ph = 1; elseif st == 2, ph = 5; elseif st == 3, ph = 6; end
end

% -------------------------------------------------------------------------
% Sampling period of the controller (s). MUST equal the sample time of the input ZOH
% blocks and of the MATLAB Function block in the model (cfg.Tc in bench_config).
function ts = ctrl_ts(mode, val)
% mode = 1 sets ts = val (20 kHz mode, abl bit 1024); mode = 0 reads
persistent tsv
if isempty(tsv), tsv = 5e-6; end
if mode == 1, tsv = val; end
ts = tsv;
end

function f = rskip_flag(mode, val)
% abl bit 8 for dpso_core (mode = 1 sets, mode = 0 reads)
persistent fv
if isempty(fv), fv = false; end
if mode == 1, fv = val; end
f = fv;
end
