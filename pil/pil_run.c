/* Closed loop: controller (20 kHz) + averaged plant (5 us steps). Computes the indices on line. */
#include "pil_run.h"
#include "pil_ctrl.h"
#include "pil_plant.h"
#include <math.h>
static uint32_t ts_zero(void) { return 0; }
uint32_t (*pil_ts)(void) = ts_zero;
static uint32_t nz;
static float nrnd(void)
{
    float u1, u2; uint32_t x = nz;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5; u1 = ((float)(x >> 8) + 0.5f) / 16777216.0f;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5; u2 = ((float)(x >> 8) + 0.5f) / 16777216.0f; nz = x;
    return sqrtf(-2.0f * logf(u1)) * cosf(6.2831853f * u2);
}
static ctrl_t C;     /* static: large struct kept out of the stack */
static plant_t Pl;
void run_case(int alg, uint32_t abl, uint32_t seed, int pack, const float *Itab, int nTab, float dv,
              float soc0, float Tend, float Ilim, float noiseV, float noiseI, res_t *r, float *logP, int nLog)
{
    const float TSC = 50e-6f, DT = 5e-6f;
    const int32_t N = (int32_t)(Tend / TSC + 0.5f), NWIN = (int32_t)(0.05f / TSC + 0.5f);
    float pwin[50]; int nw = 0;             /* 1 ms averages in the last 50 ms */
    float acc1 = 0, accI = 0, accV = 0, ib1 = 0; int n1 = 0, i; int32_t k;
    uint64_t cycSum = 0, cycSumM[3] = {0, 0, 0}; uint32_t nM[3] = {0, 0, 0}; int md;
    ctrl_out_t o = {0};
    ctrl_init(&C, alg, abl, seed, Ilim, pack == 48 ? 54.6f : 84.0f, TSC);
    plant_init(&Pl, Itab, nTab, dv, pack, soc0);
    nz = 1000u + seed * 7919u;
    r->cycMaxM[0] = r->cycMaxM[1] = r->cycMaxM[2] = 0;
    r->Ebat = 0; r->Ipk = 0; r->tUnr = 0; r->nRst = 0; r->cycMax = 0; r->cycMin = 0xFFFFFFFFu; r->Ibss = 0; r->Vss = 0;
    for (k = 0; k < N; k++) {
        float ipv = plant_ipv(&Pl, Pl.v), vb = plant_vbat(&Pl);
        float vm = Pl.v + noiseV * nrnd(), im = ipv + noiseI * nrnd();
        uint32_t t0 = pil_ts(), t1, dc;
        ctrl_step(&C, vm, im, vb, Pl.iL, &o);
        t1 = pil_ts(); dc = t0 - t1;       /* down-counter */
        if (dc > r->cycMax) r->cycMax = dc;
        if (dc < r->cycMin) r->cycMin = dc;
        cycSum += dc;
        md = (o.phase == 1 || o.phase == 5) ? 0 : (o.phase == 20 ? 2 : 1);
        cycSumM[md] += dc; nM[md]++; if (dc > r->cycMaxM[md]) r->cycMaxM[md] = dc;
        for (i = 0; i < 10; i++) {          /* plant: 10 x 5 us per control period */
            plant_step(&Pl, o.D, DT, 2);
            r->Ebat += plant_vbat(&Pl) * Pl.iL * DT;
        }
        if (k % 200 == 199) { /* OCV update every 10 ms */
            extern float plant_ocv_update(plant_t *); plant_ocv_update(&Pl);
        }
        ipv = plant_ipv(&Pl, Pl.v);
        acc1 += Pl.v * ipv; ib1 += Pl.iL; n1++;
        if (o.Vref < plant_vbat(&Pl) / 0.95f - 1.5f) r->tUnr += TSC;
        if (n1 == 20) {                     /* 1 ms */
            float p1 = acc1 / n1, i1 = ib1 / n1;
            if (i1 > r->Ipk) r->Ipk = i1;
            if (k >= N - NWIN) { if (nw < 50) pwin[nw++] = p1; }
            if (logP && (k / 20) < nLog) logP[k / 20] = p1;
            acc1 = 0; ib1 = 0; n1 = 0;
        }
        if (k >= N - NWIN) { accV += Pl.v; accI += Pl.iL; }
    }
    /* median of the 1 ms powers in the last 50 ms */
    for (i = 1; i < nw; i++) { float t = pwin[i]; int j = i - 1; while (j >= 0 && pwin[j] > t) { pwin[j+1] = pwin[j]; j--; } pwin[j+1] = t; }
    r->Pss = nw ? (nw % 2 ? pwin[nw/2] : 0.5f * (pwin[nw/2-1] + pwin[nw/2])) : 0;
    r->Vss = accV / NWIN; r->Ibss = accI / NWIN; r->nRst = o.rst;
    r->cycMean = (uint32_t)(cycSum / (N ? N : 1));
    for (md = 0; md < 3; md++) r->cycMeanM[md] = nM[md] ? (uint32_t)(cycSumM[md] / nM[md]) : 0;
}
