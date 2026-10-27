/* Averaged plant (PV string, buck converter, battery) co-simulated on the target for PIL. */
#ifndef PIL_PLANT_H
#define PIL_PLANT_H
typedef struct {
    const float *Itab; int nTab; float dv;      /* PV current I(V) on a uniform grid */
    float v, vd, iL, soc, Q, Rb, Ns;
    const float *socT, *ocvT; int nOcv;
    float ocv;
} plant_t;
void plant_init(plant_t *p, const float *Itab, int nTab, float dv, int pack, float soc0);
void plant_step(plant_t *p, float D, float dt, int nsub);
float plant_ipv(const plant_t *p, float v);
float plant_vbat(const plant_t *p);
#endif
