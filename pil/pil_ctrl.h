/* =========================================================================
 * RA-GMPPT charger controller - C implementation for processor-in-the-loop (PIL)
 * on TMS320F28335 (C28x + FPU32, single precision). Same algorithm as
 * matlab/Charger_Bench.m, controller part only (alg 5/6 DPSO/RA-DPSO, 7/8 CS/RA-CS),
 * abl bits 2 (fast left-flank limiter), 4 (time-based re-evaluation of a bound hold every 60 s), 8 (DPSO restart rule: first comparison after convergence skipped),
 * 16 (RA off), 32 (conventional limiter),
 * 64 (R2 off). The controller runs at 20 kHz (one sample per switching period).
 * ========================================================================= */
#ifndef PIL_CTRL_H
#define PIL_CTRL_H
#include <stdint.h>

#define NP 3

typedef struct {               /* DPSO / CS shared core state */
    int st, warmed;
    int32_t cnt, hold;
    float xc;
    float X[NP], V[NP], Pb[NP], PbX[NP], P[NP], Xn[NP];
    float Pm[NP], Vm[NP], Im[NP];
    float Gb, GbX;
    int k, it, g, mode, kA, bndM;
    int dropI[NP];
    float pAcc, vAcc, iAcc; int nAcc;
    float hcS, hcP, Pw; int hcDir, nw;
    int rst, atB, lostI, r2;
    float xhiS;
    int32_t tRe;               /* evaluations held at the bound (time-based re-evaluation) */
    uint32_t rng;              /* xorshift32 state (CS) */
} core_t;

typedef struct {
    int alg; uint32_t abl; float Ilim, Vcv, TS;
    /* constants derived from TS, computed once in ctrl_init */
    int32_t NPWM, LIM_PER, T_ENT, T_FAST, SS_FAST, SS_SLOW, SETTLE, LONG_SETTLE, WARMUP, NAVG;
    float aF, dxSR, KITS;
    /* filters, limiter */
    float VbF, IbF, IpF;
    int lim, side, fromL, exhaust, slewing, ssDone, fInh;
    int32_t cnt, exCnt, entCnt, fCnt, kT;
    float xL, xPk, xLast, xTrim, xCmd, xFFl, Icv, Ip0, lastXh;
    int lastPh, lastRst, restartCore;
    /* voltage loop */
    float vAcc; int nAcc; float Dcmd, integ;
    core_t core;
} ctrl_t;

typedef struct { float D, Vref, VpvMin; int phase, mode, rst; } ctrl_out_t;

void ctrl_init(ctrl_t *c, int alg, uint32_t abl, uint32_t seed, float Ilim, float Vcv, float TS);
void ctrl_step(ctrl_t *c, float Vpv, float Ipv, float Vbat, float Ibat, ctrl_out_t *o);
float ra_project(float VbF, float xcore);

#endif
