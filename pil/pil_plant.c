#include "pil_plant.h"
#define L_   1e-3f
#define CIN  100e-6f
#define RD   3.3f
#define CD   150e-6f
#define VFDB 0.5f
#define VFDF 0.5f
#define RON  0.02f
static const float SOCt[17] = {0.0000f,0.0413f,0.0855f,0.1150f,0.1593f,0.1888f,0.2625f,0.3363f,0.4100f,0.4838f,0.5575f,0.6313f,0.7050f,0.7788f,0.8525f,0.9263f,1.0000f};
static const float OCVt[17] = {2.807f,3.007f,3.207f,3.307f,3.427f,3.447f,3.517f,3.587f,3.657f,3.727f,3.807f,3.867f,3.917f,3.997f,4.067f,4.107f,4.167f};
static float interp_ocv(float s)
{
    int i;
    if (s <= SOCt[0]) return OCVt[0];
    for (i = 1; i < 17; i++) if (s <= SOCt[i]) return OCVt[i-1] + (s - SOCt[i-1]) * (OCVt[i] - OCVt[i-1]) / (SOCt[i] - SOCt[i-1]);
    return OCVt[16];
}
void plant_init(plant_t *p, const float *Itab, int nTab, float dv, int pack, float soc0)
{
    p->Itab = Itab; p->nTab = nTab; p->dv = dv;
    p->v = 0; p->vd = 0; p->iL = 0; p->soc = soc0;
    p->Ns = pack == 48 ? 13.0f : 20.0f; p->Q = 6.0f * 3.39f; p->Rb = p->Ns * 0.045f / 6.0f;
    p->ocv = p->Ns * interp_ocv(soc0);
}
static void pv_ig(const plant_t *p, float v, float *i, float *g)
{
    float x; int k;
    if (v <= 0) { *i = p->Itab[0]; *g = 0; return; }
    x = v / p->dv; k = (int)x;
    if (k >= p->nTab - 1) { *i = 0; *g = 0; return; }
    *g = (p->Itab[k+1] - p->Itab[k]) / p->dv;
    *i = p->Itab[k] + (x - k) * (p->Itab[k+1] - p->Itab[k]);
}
float plant_ipv(const plant_t *p, float v) { float i, g; pv_ig(p, v, &i, &g); return i; }
float plant_vbat(const plant_t *p) { return p->ocv + p->Rb * p->iL; }
void plant_step(plant_t *p, float D, float dt, int nsub)
{
    float h = dt / nsub; int n;
    for (n = 0; n < nsub; n++) {
        float ipv, g, f, vL;
        pv_ig(p, p->v, &ipv, &g);
        f = ipv - D * p->iL - (p->v - p->vd) / RD;
        p->v += h * f / (CIN - h * (g - 1.0f / RD));
        p->vd += h * (p->v - p->vd) / (RD * CD);
        vL = D * (p->v - VFDB - RON * p->iL) - (1.0f - D) * VFDF - (p->ocv + p->Rb * p->iL);
        p->iL += h * vL / L_;
        if (p->iL < 0) p->iL = 0;
        p->soc += h * p->iL / (3600.0f * p->Q);
    }
}
float plant_ocv_update(plant_t *p) { p->ocv = p->Ns * interp_ocv(p->soc); return p->ocv; }
