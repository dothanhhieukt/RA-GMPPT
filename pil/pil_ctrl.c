/* RA-GMPPT charger controller (see pil_ctrl.h). Line-by-line equivalent of matlab/Charger_Bench.m;
   constants derived from TS are computed once in ctrl_init. */
#include "pil_ctrl.h"
#include <math.h>
#include <string.h>

#define VB      200.0f
#define DMAXB   0.95f
#define HS      4
#define XLO     0.33f
#define XHI0    0.90f

static int32_t rnd_i(float x) { return (int32_t)(x + 0.5f); }
static float fminf_(float a, float b) { return a < b ? a : b; }
static float fmaxf_(float a, float b) { return a > b ? a : b; }
static float clampf(float x, float lo, float hi) { return x < lo ? lo : (x > hi ? hi : x); }
static int phase_of(int st) { return st == 1 ? 1 : (st == 2 ? 5 : (st == 3 ? 6 : 0)); }

/* ---------------- random numbers for CS (xorshift32, Box-Muller) ---------------- */
static float urand(core_t *s)
{
    uint32_t x = s->rng;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    s->rng = x;
    return ((float)(x >> 8) + 0.5f) * (1.0f / 16777216.0f);
}
static float nrand(core_t *s)
{
    float u1 = urand(s), u2 = urand(s);
    return sqrtf(-2.0f * logf(u1)) * cosf(6.2831853f * u2);
}
static float levy_signed(core_t *s)
{
    const float sig = 0.6965745f;      /* Mantegna sigma_u for beta = 1.5 */
    float u = nrand(s) * sig, v = nrand(s);
    return u / powf(fabsf(v), 1.0f / 1.5f);
}
static float ra_draw(core_t *s, float xhi, int occB)
{
    float xs = fminf_(xhi, XLO + urand(s) * (XHI0 - XLO));
    if (xs >= xhi - 1e-7f && xhi < XHI0 - 1e-7f && occB) xs = XLO + urand(s) * (xhi - XLO);
    return xs;
}

/* ---------------- DPSO (Ishaque & Salam 2013) with RA layer ---------------- */
static void dpso_step_p(core_t *s, int i, float xhi)
{
    const float W = 0.4f, VMAX = 0.075f;
    float vi = W * s->V[i] + (s->PbX[i] - s->X[i]);
    if (s->Gb >= 0.0f) vi += (s->GbX - s->X[i]);
    vi = clampf(vi, -VMAX, VMAX);
    if (fabsf(vi) < 1e-4f) vi = 0.01f * (float)(2 * (s->it % 2) - 1);
    s->V[i] = vi;
    s->X[i] = clampf(s->X[i] + vi, XLO, xhi);
}

static float dpso_core(ctrl_t *c, float Vpv, float Ipv, int useRA, float xmaxDyn, int restart, int useR2,
                       float xWarm, int *ph, int *rst)
{
    core_t *s = &c->core;
    const float VSTOP = 0.01f, DSPREAD = 0.45f * 37.2f / 200.0f, HC_S0 = 0.01f, HC_SMIN = 0.0025f;
    const float PTHR = 0.10f, DXR = 2.0f / 200.0f;
    const int IT_MAX = 20, NWIN = 20;
    const int32_t SETTLE = c->SETTLE, LONG_SETTLE = c->LONG_SETTLE, WARMUP = c->WARMUP, NAVG = c->NAVG;
    float xhi = XHI0, xr, P, Vme, Ime, xn; int goGlobal = 0, i, j;
    if (useRA) xhi = fmaxf_(XLO + 0.05f, fminf_(XHI0, xmaxDyn));
    xr = fmaxf_(XLO + 0.05f, fminf_(XHI0, xmaxDyn));
    if (restart) { s->cnt = 0; if (s->hold > SETTLE) s->hold = SETTLE; s->pAcc = s->vAcc = s->iAcc = 0; s->nAcc = 0; }
    if (s->cnt < s->hold) {
        s->cnt++;
        if (s->cnt > s->hold - NAVG) { s->pAcc += Vpv * Ipv; s->vAcc += Vpv; s->iAcc += Ipv; s->nAcc++; }
        if (useRA && s->st == 3) {                         /* RA: held point follows / is clipped */
            if (s->atB) s->xc = xhi;
            else { if (s->xc > xhi + 1e-6f) s->lostI = 1; s->xc = fminf_(s->xc, xhi); }
        }
        *ph = phase_of(s->st); *rst = s->rst; return s->xc;
    }
    P = Vpv * Ipv; Vme = Vpv; Ime = Ipv;
    if (s->nAcc > 0) { P = s->pAcc / s->nAcc; Vme = s->vAcc / s->nAcc; Ime = s->iAcc / s->nAcc; }
    s->pAcc = s->vAcc = s->iAcc = 0; s->nAcc = 0;
    xn = s->xc;
    switch (s->st) {
    case 0:
        if (!s->warmed) { s->warmed = 1; s->xc = xWarm; s->cnt = 1; s->hold = WARMUP; *ph = 0; *rst = s->rst; return s->xc; }
        goGlobal = 1; break;
    case 1:
        i = s->k;
        s->Pm[i] = P; s->Vm[i] = Vme; s->Im[i] = Ime; s->Pb[i] = P; s->PbX[i] = s->X[i];
        if (P > s->Gb) { s->Gb = P; s->GbX = s->X[i]; }
        if (i < NP - 1) { s->k++; xn = s->X[s->k]; }
        else {
            int shaded = (s->Im[2] - s->Im[1]) / fmaxf_(s->Im[2], 1e-3f) >= 0.1f && (s->Vm[0] - s->Vm[1]) / fmaxf_(s->Vm[0], 1.0f) >= 0.2f;
            if (shaded) { s->st = 2; s->it = 1; s->k = 0; dpso_step_p(s, 0, xhi); xn = s->X[0]; }
            else { s->st = 3; xn = s->GbX; s->hcS = HC_S0; s->hcDir = 1; s->hcP = s->Gb; s->Pw = (c->abl & 8u) ? -1.0f : s->Gb; s->nw = 0; s->xhiS = xr; }
        }
        break;
    case 2:
        i = s->k;
        if (P > s->Pb[i]) { s->Pb[i] = P; s->PbX[i] = s->X[i]; }
        if (P > s->Gb) { s->Gb = P; s->GbX = s->X[i]; }
        if (i < NP - 1) { s->k++; dpso_step_p(s, s->k, xhi); xn = s->X[s->k]; }
        else {
            int doStop;
            s->it++;
            doStop = s->it >= IT_MAX;
            for (i = 0; i < NP; i++) for (j = 0; j < NP; j++)
                if (i != j && fabsf(s->X[i] - s->X[j]) < DSPREAD && fabsf(s->V[i]) < VSTOP) doStop = 1;
            if (doStop) { s->st = 3; xn = s->GbX; s->hcS = HC_S0; s->hcDir = 1; s->hcP = s->Gb; s->Pw = (c->abl & 8u) ? -1.0f : s->Gb; s->nw = 0; s->xhiS = xr; }
            else { s->k = 0; dpso_step_p(s, 0, xhi); xn = s->X[0]; }
        }
        break;
    default:
        if (useRA && useR2 && s->atB && fabsf(xr - s->xhiS) > DXR) { s->rst++; s->r2++; goGlobal = 1; }   /* R2 */
        if (useRA && useR2 && (c->abl & 4u) && (s->atB || fabsf(s->xc - xhi) <= DXR)) { if (++s->tRe >= rnd_i(60.0f / (650e-6f * HS))) { s->tRe = 0; if (!goGlobal) { s->rst++; s->r2++; goGlobal = 1; } } } else s->tRe = 0;   /* time-based re-evaluation (abl 4) */
        s->lostI = 0;
        s->nw++;
        if (s->nw >= NWIN) {
            if (s->Pw > 0 && fabsf(P - s->Pw) / s->Pw > PTHR) { s->rst++; goGlobal = 1; }
            s->Pw = P; s->nw = 0;
        }
        if (!goGlobal) {
            if (P < s->hcP) { s->hcDir = -s->hcDir; s->hcS = fmaxf_(HC_SMIN, s->hcS * 0.5f); }
            s->hcP = P;
            xn = clampf(s->xc + s->hcDir * s->hcS, XLO, xhi);
        }
    }
    if (goGlobal) {
        s->st = 1; s->k = 0; s->Gb = -1; s->it = 0; s->xhiS = xr;
        for (i = 0; i < NP; i++) {
            s->X[i] = useRA ? XLO + (float)(i + 1) / NP * (xhi - XLO) : XLO + ((float)i + 0.5f) / NP * (xhi - XLO);
            s->V[i] = 0; s->Pb[i] = -1;
        }
        xn = s->X[0];
    }
    if (useRA) xn = fminf_(xn, xhi);                         /* R1: projection */
    s->atB = useRA && s->st == 3 && xn >= xhi - 1e-6f;
    s->hold = fabsf(xn - s->xc) > 0.05f ? LONG_SETTLE : SETTLE;
    s->xc = xn; s->cnt = 1;
    *ph = phase_of(s->st); *rst = s->rst; return s->xc;
}

/* ---------------- CS (Ahmed & Salam 2014) with RA layer ---------------- */
static float cs_core(ctrl_t *c, float Vpv, float Ipv, int useRA, float xmaxDyn, int restart, int useR2,
                     float xWarm, int *ph, int *rst)
{
    core_t *s = &c->core;
    const float KAPPA = 0.05f, PA = 0.25f, DCONV = 0.01f, PTHR = 0.10f, DXR = 2.0f / 200.0f;
    const int GMAX = 25, NWIN = 20;
    const int32_t SETTLE = c->SETTLE, LONG_SETTLE = c->LONG_SETTLE, WARMUP = c->WARMUP, NAVG = c->NAVG;
    float xhi = XHI0, xr, Pm, xn; int goGlobal = 0, i, j, endGen;
    if (useRA) xhi = fmaxf_(XLO + 0.05f, fminf_(XHI0, xmaxDyn));
    xr = fmaxf_(XLO + 0.05f, fminf_(XHI0, xmaxDyn));
    if (restart) { s->cnt = 0; if (s->hold > SETTLE) s->hold = SETTLE; s->pAcc = 0; s->nAcc = 0; }
    if (s->cnt < s->hold) {
        s->cnt++;
        if (s->cnt > s->hold - NAVG) { s->pAcc += Vpv * Ipv; s->nAcc++; }
        if (useRA && s->st == 3) {
            if (s->atB) s->xc = xhi;
            else { if (s->xc > xhi + 1e-6f) s->lostI = 1; s->xc = fminf_(s->xc, xhi); }
        }
        *ph = phase_of(s->st); *rst = s->rst; return s->xc;
    }
    Pm = Vpv * Ipv; if (s->nAcc > 0) Pm = s->pAcc / s->nAcc;
    s->pAcc = 0; s->nAcc = 0;
    xn = s->xc;
    switch (s->st) {
    case 0:
        if (!s->warmed) { s->warmed = 1; s->xc = xWarm; s->cnt = 1; s->hold = WARMUP; *ph = 0; *rst = s->rst; return s->xc; }
        goGlobal = 1; break;
    case 1:
        if (s->bndM) {                                       /* R1: boundary measured first, extra candidate */
            s->bndM = 0; s->Gb = Pm; s->GbX = s->xc; xn = s->X[0]; s->k = 0;
            s->hold = fabsf(xn - s->xc) > 0.05f ? LONG_SETTLE : SETTLE;
            s->xc = xn; s->cnt = 1; *ph = phase_of(s->st); *rst = s->rst; return s->xc;
        }
        i = s->k; s->P[i] = Pm;
        if (Pm > s->Gb) { s->Gb = Pm; s->GbX = s->X[i]; }
        if (i < NP - 1) { s->k++; xn = s->X[s->k]; }
        else {
            s->st = 2; s->g = 0; s->k = 0; s->mode = 1;
            for (i = 0; i < NP; i++) s->Xn[i] = clampf(s->X[i] + KAPPA * (s->GbX - s->X[i]) * levy_signed(s), XLO, xhi);
            xn = s->Xn[0];
        }
        break;
    case 2:
        endGen = 0;
        if (s->mode == 1) {
            i = s->k;
            if (Pm > s->P[i]) { s->P[i] = Pm; s->X[i] = s->Xn[i]; }
            if (Pm > s->Gb) { s->Gb = Pm; s->GbX = s->Xn[i]; }
            s->k++;
            if (s->k <= NP - 1) xn = s->Xn[s->k];
            else {
                float mx = s->X[0], mn = s->X[0];
                for (i = 1; i < NP; i++) { mx = fmaxf_(mx, s->X[i]); mn = fminf_(mn, s->X[i]); }
                if (mx - mn < DCONV || s->g + 1 >= GMAX) { s->st = 3; xn = s->GbX; s->Pw = s->Gb; s->nw = 0; s->xhiS = xr; }
                else {
                    int nA = 0;
                    for (i = 0; i < NP; i++) {
                        s->dropI[i] = fabsf(s->X[i] - s->GbX) > 1e-7f && urand(s) < PA;
                        if (s->dropI[i]) {
                            int occB = 0;
                            nA++;
                            for (j = 0; j < NP; j++) if (j != i) {
                                float xj = (j < i && s->dropI[j]) ? s->Xn[j] : s->X[j];
                                if (xj >= xhi - 1e-7f) occB = 1;
                            }
                            s->Xn[i] = ra_draw(s, xhi, occB);
                        }
                    }
                    if (nA > 0) { s->kA = 0; while (!s->dropI[s->kA]) s->kA++; s->mode = 2; xn = s->Xn[s->kA]; }
                    else endGen = 1;
                }
            }
        } else {
            i = s->kA;
            s->X[i] = s->Xn[i]; s->P[i] = Pm; s->dropI[i] = 0;
            if (Pm > s->Gb) { s->Gb = Pm; s->GbX = s->X[i]; }
            if (s->dropI[0] || s->dropI[1] || s->dropI[2]) { s->kA = 0; while (!s->dropI[s->kA]) s->kA++; xn = s->Xn[s->kA]; }
            else endGen = 1;
        }
        if (endGen) {
            s->g++;
            if (s->g >= GMAX) { s->st = 3; xn = s->GbX; s->Pw = s->Gb; s->nw = 0; s->xhiS = xr; }
            else {
                s->mode = 1; s->k = 0;
                for (i = 0; i < NP; i++) s->Xn[i] = clampf(s->X[i] + KAPPA * (s->GbX - s->X[i]) * levy_signed(s), XLO, xhi);
                xn = s->Xn[0];
            }
        }
        break;
    default:
        if (useRA && useR2 && s->atB && fabsf(xr - s->xhiS) > DXR) { s->rst++; s->r2++; goGlobal = 1; }   /* R2 */
        if (useRA && useR2 && (c->abl & 4u) && (s->atB || fabsf(s->xc - xhi) <= DXR)) { if (++s->tRe >= rnd_i(60.0f / (650e-6f * HS))) { s->tRe = 0; if (!goGlobal) { s->rst++; s->r2++; goGlobal = 1; } } } else s->tRe = 0;   /* time-based re-evaluation (abl 4) */
        s->lostI = 0;
        s->nw++;
        if (!goGlobal && s->nw >= NWIN) {
            if (s->Pw > 0 && fabsf(Pm - s->Pw) / s->Pw > PTHR) { s->rst++; goGlobal = 1; }
            s->Pw = Pm; s->nw = 0;
        }
        if (!goGlobal) xn = s->xc;
    }
    if (goGlobal) {
        s->st = 1; s->k = 0; s->Gb = -1; s->g = 0; s->mode = 1; s->xhiS = xr;
        for (i = 0; i < NP; i++) {
            int occB = 0;
            for (j = 0; j < i; j++) if (s->X[j] >= xhi - 1e-7f) occB = 1;
            s->X[i] = ra_draw(s, xhi, occB); s->P[i] = -1;
        }
        xn = s->X[0];
        if (useRA) { s->bndM = 1; xn = xhi; }
    }
    if (useRA) xn = fminf_(xn, xhi);
    s->atB = useRA && s->st == 3 && xn >= xhi - DXR;
    if (s->atB) xn = xhi;
    s->hold = fabsf(xn - s->xc) > 0.05f ? LONG_SETTLE : SETTLE;
    s->xc = xn; s->cnt = 1;
    *ph = phase_of(s->st); *rst = s->rst; return s->xc;
}

/* ---------------- top level: RA boundary, CC/CV limiter, PV-voltage loop ---------------- */
void ctrl_init(ctrl_t *c, int alg, uint32_t abl, uint32_t seed, float Ilim, float Vcv, float TS)
{
    memset(c, 0, sizeof(*c));
    c->alg = alg; c->abl = abl; c->Ilim = Ilim; c->Vcv = Vcv; c->TS = TS;
    /* TS-derived constants (same expressions as before, evaluated once) */
    c->NPWM = rnd_i(50e-6f / TS) < 1 ? 1 : rnd_i(50e-6f / TS);
    c->LIM_PER = rnd_i(650e-6f * HS / TS); c->T_ENT = rnd_i(3e-3f / TS);
    c->T_FAST = rnd_i(0.3e-3f / TS); if (c->T_FAST < 1) c->T_FAST = 1;
    c->SS_FAST = rnd_i(2e-3f / TS); c->SS_SLOW = rnd_i(19.5e-3f / TS);
    c->SETTLE = rnd_i(650e-6f * HS / TS); c->LONG_SETTLE = rnd_i(1300e-6f * HS / TS);
    c->WARMUP = rnd_i(20e-3f / TS); c->NAVG = rnd_i(200e-6f / TS); if (c->NAVG < 1) c->NAVG = 1;
    c->aF = TS / 5e-4f;
    c->dxSR = 10e3f * c->NPWM * TS / VB;
    c->KITS = 1000.0f * c->NPWM * TS;
    c->VbF = -1.0f;                       /* set at first call */
    c->xL = 0.5f; c->xPk = 0.5f; c->side = 1; c->xLast = 0.4f; c->lastXh = 0.4f;
    c->xCmd = 0.4f; c->xFFl = 0.4f; c->Icv = 10.0f;
    c->core.xc = 0.4f; c->core.GbX = 0.4f; c->core.Gb = -1; c->core.Pw = -1; c->core.hcS = 0.01f; c->core.hcDir = 1;
    c->core.rng = 2463534242u ^ (seed * 2654435761u); if (c->core.rng == 0) c->core.rng = 1;
    { int i; for (i = 0; i < NP; i++) { c->core.Pb[i] = -1; c->core.P[i] = -1; } }
}

void ctrl_step(ctrl_t *c, float Vpv, float Ipv, float Vbat, float Ibat, ctrl_out_t *o)
{
    const float aF = c->aF, KP = 12.0f;
    const int32_t NPWM = c->NPWM, LIM_PER = c->LIM_PER, T_ENT = c->T_ENT, T_FAST = c->T_FAST;
    const float kL = 0.05f, kR = 0.01f, dxMax = 0.01f, ETA = 0.97f, kC = 0.02f;
    const float dxSR = c->dxSR;
    const int NEX = 5, alg = c->alg;
    const uint32_t abl = c->abl;
    float VpvMin, VpvMinC, xmaxC, xmax, x, Vref, e, xWarm;
    int propC, rightLim, raOn, useR2, fast, wasLim, ph = c->lastPh, rs = c->lastRst;
    if (c->VbF < 0) c->VbF = Vbat;
    c->VbF += aF * (Vbat - c->VbF); c->IbF += aF * (Ibat - c->IbF); c->IpF += aF * (Ipv - c->IpF);
    /* ---- RA layer, eq. (4): reachability limit from the filtered battery voltage ---- */
    { float vd = fmaxf_(c->VbF, 1.0f) / DMAXB;          /* one division instead of two (identical result) */
      VpvMin = vd + 1.0f; VpvMinC = vd - 0.5f; }
    xmaxC = 1.0f - VpvMinC / VB; xmax = 1.0f - VpvMin / VB;
    propC = (alg == 6 || alg == 8);
    rightLim = !propC || (abl & 32u);
    raOn = (abl & 16u) == 0; useR2 = raOn && (abl & 64u) == 0;
    fast = !rightLim && (abl & 2u);
    xWarm = (!rightLim && propC) ? xmaxC : 0.4f;
    c->vAcc += Vpv; c->nAcc++;
    o->mode = 0;
    {
        float IlimS = fmaxf_(c->Ilim, 1e-3f);
        e = fmaxf_(c->IbF / IlimS - 1.0f, (c->VbF - c->Vcv) * 2.0f);   /* x/0.5 == 2x exactly */
        wasLim = c->lim;
        c->entCnt = e > 0.02f ? c->entCnt + 1 : 0;
        if (c->lastPh == 1 || c->lastPh == 5) c->fInh = 0;
        c->fCnt = (fast && !c->fInh && c->lastPh == 6 && e > -0.05f) ? c->fCnt + 1 : 0;
        c->kT++;
        if (!c->ssDone && c->kT >= (fast ? c->SS_FAST : c->SS_SLOW) && !rightLim && propC) {
            c->ssDone = 1;
            if (!c->lim) {
                c->lim = 1; c->cnt = 0; c->xTrim = 0; c->slewing = 0; c->exCnt = 0; c->exhaust = 0; c->fromL = 0;
                c->entCnt = 0; c->Icv = IlimS; c->xPk = 1.0f - 119.0f / VB;
                c->side = 1; c->xL = xmax; c->xCmd = xmax; c->xFFl = xmax;
            }
        }
        if (!c->lim && (c->entCnt >= T_ENT || (fast && c->fCnt >= T_FAST))) {
            int fastE = c->entCnt < T_ENT;
            c->entCnt = 0; c->fCnt = 0;
            c->lim = 1; c->cnt = 0; c->xTrim = 0; c->slewing = 0; c->exCnt = 0; c->exhaust = 0; c->fromL = 0;
            c->xPk = (c->lastPh == 6) ? c->lastXh : c->xLast;
            c->xL = c->xLast; c->xCmd = c->xLast; c->xFFl = c->xPk;
            c->Icv = fminf_(IlimS, fmaxf_(c->IbF, 0.2f));
            if (fastE) { c->Icv = IlimS; c->xPk = fminf_(c->xPk, 1.0f - 119.0f / VB); }
            if (rightLim) { c->side = -1; c->xL = 1.0f - 150.0f / VB; c->fromL = 1; c->Ip0 = fmaxf_(c->IpF, 1.0f); }
            else c->side = 1;
        }
        if (c->lim) {
            float Iset = fmaxf_(0.2f, fminf_(IlimS, c->Icv)), eI = c->IbF / Iset - 1.0f;
            c->cnt++;
            if (c->cnt >= LIM_PER) {
                c->cnt = 0;
                if (c->side > 0) {
                    float xUn;
                    c->exCnt = c->exhaust ? c->exCnt + 1 : 0;
                    c->Icv = fminf_(IlimS, fmaxf_(0.0f, c->Icv - IlimS * (fast ? 6.0f * kC : kC) * (c->VbF - c->Vcv) * 2.0f));
                    xUn = c->xFFl + c->xTrim;
                    if (!c->slewing && !(xUn >= xmax && eI > 0) && !(xUn <= c->xPk && eI < 0))
                        c->xTrim = clampf(c->xTrim + clampf(kL * eI, -dxMax, dxMax), -1.0f, 1.0f);
                } else {
                    if (c->fromL && c->IpF < 0.15f * c->Ip0 && e < -0.1f) c->xL += 2.0f * dxMax;
                    else c->xL -= clampf(kR * e, -dxMax, dxMax);
                    c->xL = clampf(c->xL, 0.15f, xmax);
                    if (c->xL >= c->xPk && e < -0.02f) c->lim = 0;
                }
            }
            if (c->lim && c->side > 0 && c->nAcc >= NPWM) {
                float xFF, xDes; int cvNoLeft;
                if (c->IpF > 0.5f) {
                    float Vt = c->VbF * Iset / (ETA * c->IpF);
                    xFF = 1.0f - fminf_(VB * (1.0f - c->xPk), fmaxf_(VpvMin, Vt)) / VB;
                } else xFF = c->xPk;
                c->xFFl = xFF;
                xDes = clampf(xFF + c->xTrim, c->xPk, xmax);
                if (fast && xDes > c->xCmd) { c->slewing = 0; c->xCmd = xDes; }
                else { c->slewing = fabsf(xDes - c->xCmd) > dxSR; c->xCmd += clampf(xDes - c->xCmd, -dxSR, dxSR); }
                c->exhaust = !c->slewing && c->xCmd >= xmax - 1e-6f && e > 0.02f;
                cvNoLeft = c->VbF >= c->Vcv - 0.1f && c->Icv < IlimS - 1e-6f && c->IpF > 0.5f && VpvMin * c->IpF * ETA > 1.05f * c->VbF * Iset;
                if (fast && c->xCmd >= xmax - 1e-6f && c->IpF > 0.5f && Vpv >= VpvMin - 1.5f && Vpv * c->IpF * ETA > 1.05f * c->VbF * Iset) cvNoLeft = 1;
                if ((c->exhaust && c->exCnt >= NEX) || cvNoLeft) { c->side = -1; c->xL = 1.0f - 150.0f / VB; c->fromL = 1; c->Ip0 = fmaxf_(c->IpF, 1.0f); }
                else if (xDes <= c->xPk + 1e-6f && c->xCmd <= c->xPk + 1e-6f && eI < -0.02f && c->Icv >= IlimS - 1e-6f) c->lim = 0;
                else c->xL = c->xCmd;
            }
            if (wasLim && !c->lim) { c->restartCore = 1; c->fInh = 1; }
            x = c->xL; ph = 20; o->mode = c->side > 0 ? 1 : 2;
        } else {
            int phc, rsc;
            if (alg == 5 || alg == 6) x = dpso_core(c, Vpv, Ipv, alg == 6 && raOn, xmaxC, c->restartCore, useR2, xWarm, &phc, &rsc);
            else x = cs_core(c, Vpv, Ipv, alg == 8 && raOn, xmaxC, c->restartCore, useR2, xWarm, &phc, &rsc);
            c->restartCore = 0;
            c->xLast = x; c->lastPh = phc; c->lastRst = rsc; c->lastXh = x; ph = phc; rs = rsc;
        }
    }
    /* ---- PV-voltage loop, eq. (2) ---- */
    Vref = VB * (1.0f - x);
    if (c->nAcc >= NPWM) {
        float vAvg = c->nAcc == 1 ? c->vAcc : c->vAcc / c->nAcc, Dff = c->VbF / fmaxf_(Vref, 1.0f), ev = (vAvg - Vref) / fmaxf_(Vref, 1.0f);
        float Dnew = Dff * (1.0f + KP * ev + c->integ), Dsat = clampf(Dnew, 0.0f, DMAXB);
        c->vAcc = 0; c->nAcc = 0;
        if (Dsat == Dnew || (Dnew > DMAXB && ev < 0) || (Dnew < 0 && ev > 0)) c->integ += c->KITS * ev;
        c->Dcmd = Dsat;
    }
    o->D = c->Dcmd; o->Vref = Vref; o->VpvMin = VpvMin; o->phase = ph; o->rst = rs;
}

/* RA layer alone (eq. 4 + R1 projection), used for the cycle micro-benchmark */
float ra_project(float VbF, float xcore)
{
    float VpvMinC = fmaxf_(VbF, 1.0f) / DMAXB - 0.5f;
    float xhi = fmaxf_(XLO + 0.05f, fminf_(XHI0, 1.0f - VpvMinC / VB));
    return fminf_(xcore, xhi);
}
