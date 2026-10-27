/* One-hour charging profile with the C controller (pil_ctrl.c, 20 kHz) and the averaged plant.
   usage: lt alg abl seed tend soc0 out.csv  -> prints summary JSON; writes per-second log */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pil_ctrl.h"
#define NTB 3201
#define DVT 0.05
static double Itab[NTB];
static const double q = 1.602176634e-19, kB = 1.380649e-23;
static const double IL0 = 8.388, I00 = 4.26e-10, NID = 1.019, RSH = 148.4, RS = 0.378, EG = 1.22, AISC = 0.0996e-2;
#define NVG 2000
static void module_curve(double G, double T, double *Vm, double *Im)   /* I(V) for V = -0.79..45, increasing I order out */
{
    double Tk = T + 273.15, a = NID * 60 * kB * Tk / q, Iph = IL0 * G / 1000 * (1 + AISC * (T - 25));
    double I0 = I00 * pow(Tk / 298.15, 3) * exp(q * EG / (NID * kB) * (1 / 298.15 - 1 / Tk)), rsh = RSH * 1000 / (G > 1 ? G : 1);
    int k, it;
    for (k = 0; k < NVG; k++) {
        double V = 45.0 - (45.79) * k / (NVG - 1), I = Iph;     /* V decreasing -> I increasing */
        for (it = 0; it < 60; it++) {
            double x = (V + I * RS) / a; if (x > 60) x = 60; if (x < -60) x = -60;
            double ex = exp(x), f = Iph - I0 * (ex - 1) - (V + I * RS) / rsh - I, df = -I0 * ex * RS / a - RS / rsh - 1;
            I -= f / df;
        }
        Vm[k] = V; Im[k] = I;
    }
}
static double vmod(double I, const double *Vm, const double *Im, double Iph0)
{
    double Vbp = -(0.8 + 1e-3 * (I > 0 ? I : 0)), V;
    int k;
    if (I <= Im[0]) V = Vm[0];
    else if (I >= Im[NVG - 1]) V = -10;
    else { for (k = 1; k < NVG && Im[k] < I; k++); V = Vm[k-1] + (I - Im[k-1]) * (Vm[k] - Vm[k-1]) / (Im[k] - Im[k-1]); }
    return V > Vbp ? V : Vbp;
}
static double Vm4[4][NVG], Im4[4][NVG];
static void build_table(const double *G, double T)
{
    int m, k, n = 3000; double Imax = 0, *Vs = malloc(n * sizeof(double)), *Is = malloc(n * sizeof(double));
    for (m = 0; m < 4; m++) { module_curve(G[m] > 1 ? G[m] : 1, T, Vm4[m], Im4[m]); double ip = IL0 * G[m] / 1000 * (1 + AISC * (T - 25)); if (ip > Imax) Imax = ip; }
    Imax = Imax * 1.02 + 0.05;
    for (k = 0; k < n; k++) { double I = Imax * k / (n - 1), V = 0; for (m = 0; m < 4; m++) V += vmod(I, Vm4[m], Im4[m], 0); Is[k] = I; Vs[k] = V; }
    /* Vs decreasing with k; build I(V) on uniform grid */
    for (k = 0; k < NTB; k++) {
        double V = k * DVT; int j;
        if (V >= Vs[0]) { Itab[k] = 0; continue; }
        if (V <= Vs[n - 1]) { Itab[k] = Is[n - 1]; continue; }
        for (j = 1; j < n && Vs[j] > V; j++);
        Itab[k] = Is[j-1] + (V - Vs[j-1]) * (Is[j] - Is[j-1]) / (Vs[j] - Vs[j-1]);
    }
    for (k = NTB - 2; k >= 0; k--) if (Itab[k] < Itab[k+1]) Itab[k] = Itab[k+1];
    free(Vs); free(Is);
}
static void ig(double v, double *i, double *g)
{
    if (v <= 0) { *i = Itab[0]; *g = 0; return; }
    double x = v / DVT; int k = (int)x;
    if (k >= NTB - 1) { *i = 0; *g = 0; return; }
    *g = (Itab[k+1] - Itab[k]) / DVT; *i = Itab[k] + (x - k) * (Itab[k+1] - Itab[k]);
}
static const double SOCt[17] = {0.0000,0.0413,0.0855,0.1150,0.1593,0.1888,0.2625,0.3363,0.4100,0.4838,0.5575,0.6313,0.7050,0.7788,0.8525,0.9263,1.0000};
static const double OCVt[17] = {2.807,3.007,3.207,3.307,3.427,3.447,3.517,3.587,3.657,3.727,3.807,3.867,3.917,3.997,4.067,4.107,4.167};
static double ocv(double s) { int i; if (s <= 0) return OCVt[0]; for (i = 1; i < 17; i++) if (s <= SOCt[i]) return OCVt[i-1] + (s - SOCt[i-1]) * (OCVt[i] - OCVt[i-1]) / (SOCt[i] - SOCt[i-1]); return OCVt[16]; }
static uint32_t nz;
static double nrnd(void) { double u1, u2; uint32_t x = nz; x ^= x << 13; x ^= x >> 17; x ^= x << 5; u1 = ((x >> 8) + 0.5) / 16777216.0; x ^= x << 13; x ^= x >> 17; x ^= x << 5; u2 = ((x >> 8) + 0.5) / 16777216.0; nz = x; return sqrt(-2 * log(u1)) * cos(6.283185307 * u2); }
static double gsun[20000]; static int ng = 0;
static void profile(void)
{
    int sl[10] = {10,20,50,100,10,30,5,100,50,20}, ho[10] = {60,60,30,30,120,60,120,20,60,60}, s, t, k; double G = 300;
    for (s = 0; s < 10; s++) for (t = 0; t < 2; t++) {
        double tg = t == 0 ? 1000 : 300; int n = (int)(fabs(tg - G) / sl[s]);
        for (k = 1; k <= n; k++) gsun[ng++] = G + (tg - G) * k / n;
        G = tg; for (k = 0; k < ho[s]; k++) gsun[ng++] = G;
    }
    while (ng < 3700) { int m = ng; for (k = 0; k < m && ng < 20000; k++) gsun[ng++] = gsun[k]; }
}
static const double PAT[12][4] = {{1,1,1,1},{1,1,.6,.6},{1,.85,.5,.3},{1,.7,.7,.3},{.9,.7,.4,.2},{1,.85,.5,.3},{1,1,1,.4},{1,.6,.6,.4},{1,.85,.5,.3},{1,.9,.5,.5},{1,1,.6,.6},{1,1,1,1}};
int main(int argc, char **argv)
{
    int alg = atoi(argv[1]); uint32_t abl = (uint32_t)atoi(argv[2]); int seed = atoi(argv[3]);
    double tend = atof(argv[4]), soc = atof(argv[5]); FILE *fo = fopen(argv[6], "w"); int steady = argc > 7 && !strcmp(argv[7], "steady");
    const double dt = 25e-6, L = 1e-3, CIN = 100e-6, RD = 3.3, CD = 150e-6, VF = 0.5, RON = 0.02;
    const int pack = 72; double Ns = 20, Rb = 20 * 0.045 / 6, Q = 6 * 3.39;
    ctrl_t C; ctrl_out_t o = {0}; long n, N = (long)(tend / dt + 0.5), nsec = (long)(1.0 / dt + 0.5);
    double v = 0, vd = 0, iL = 0, D = 0, oc, E = 0, Epv = 0, tUnr = 0, tB = 0, aP = 0, aPb = 0, aVb = 0, aV = 0;
    int lastPh = 0, nScan = 0, sec = 0; double Gs[4], T;
    double ib1[40] = {0}, ib1s = 0, tOv = 0; int ib1k = 0, nOv = 0, ovOn = 0; double ipk1 = 0;   /* 1 ms average of the charging current */
    double ib5[2] = {0}, ib5s = 0, tOv5 = 0; int ib5k = 0, nOv5 = 0, ov5On = 0;   /* 50 us average (2 samples of 25 us), as t_over */
    profile(); ctrl_init(&C, alg, abl, seed, 10.0f, 84.0f, 50e-6f); nz = 1000u + seed * 7919u;
    oc = Ns * ocv(soc);
    fprintf(fo, "sec,P,Pb,Vb,Vpv,G1,G2,G3,G4,T\n");
    for (n = 0; n < N; n++) {
        if (n % nsec == 0) {
            if (n > 0) fprintf(fo, "%d,%.3f,%.3f,%.4f,%.3f,%.1f,%.1f,%.1f,%.1f,%.2f\n", sec, aP / nsec, aPb / nsec, aVb / nsec, aV / nsec, Gs[0], Gs[1], Gs[2], Gs[3], T);
            sec = (int)(n / nsec); aP = aPb = aVb = aV = 0;
            int m; double gm = 0; const double *pt = PAT[(sec / 300) % 12];
            for (m = 0; m < 4; m++) { Gs[m] = gsun[sec] * pt[m]; gm += Gs[m] / 4; }
            T = 25 + 0.025 * gm;
            if (steady) { Gs[0] = 1000; Gs[1] = 850; Gs[2] = 500; Gs[3] = 300; T = 65; }
            if (!steady || sec == 0) build_table(Gs, T);
            oc = Ns * ocv(soc);
        }
        double ipv, g; ig(v, &ipv, &g);
        if (n % 2 == 1) {
            ctrl_step(&C, (float)(v + 0.05 * nrnd()), (float)(ipv + 0.005 * nrnd()), (float)(oc + Rb * iL), (float)iL, &o);
            D = o.D;
            if (o.phase == 1 && lastPh == 6) nScan++;
            lastPh = o.phase;
            if (o.Vref < (oc + Rb * iL) / 0.95 - 1.5) tUnr += 50e-6;
            if (fabs(v - o.VpvMin + 1.0) < 1.0 && o.phase == 6) tB += 50e-6;
        }
        double f = ipv - D * iL - (v - vd) / RD;
        v += dt * f / (CIN - dt * (g - 1 / RD));
        vd += dt * (v - vd) / (RD * CD);
        double vL = D * (v - VF - RON * iL) - (1 - D) * VF - (oc + Rb * iL);
        iL += dt * vL / L; if (iL < 0) iL = 0;
        soc += dt * iL / (3600 * Q);
        ib1s += iL - ib1[ib1k]; ib1[ib1k] = iL; ib1k = (ib1k + 1) % 40;
        if (n > 40 && ib1s / 40 > ipk1) ipk1 = ib1s / 40;
        if (ib1s / 40 > 1.05 * 10.0) { tOv += dt; if (!ovOn) { nOv++; ovOn = 1; } } else ovOn = 0;
        ib5s += iL - ib5[ib5k]; ib5[ib5k] = iL; ib5k = (ib5k + 1) % 2;
        if (n > 2 && ib5s / 2 > 1.05 * 10.0) { tOv5 += dt; if (!ov5On) { nOv5++; ov5On = 1; } } else ov5On = 0;
        double pb = (oc + Rb * iL) * iL, pp = v * ipv;
        E += pb * dt; Epv += pp * dt; aP += pp; aPb += pb; aVb += oc + Rb * iL; aV += v;
    }
    fclose(fo);
    printf("{\"alg\":%d,\"abl\":%u,\"seed\":%d,\"E\":%.2f,\"Epv\":%.2f,\"socEnd\":%.5f,\"tUnr\":%.3f,\"tBnd\":%.2f,\"nScan\":%d,\"nR2\":%d,\"nRst\":%d,\"nOver\":%d,\"tOver\":%.3f,\"Ipk\":%.2f,\"nOver50\":%d,\"tOver50\":%.3f}\n",
           alg, abl, seed, E, Epv, soc, tUnr, tB, nScan, C.core.r2, o.rst, nOv, tOv, ipk1, nOv5, tOv5);
    (void)pack; return 0;
}
