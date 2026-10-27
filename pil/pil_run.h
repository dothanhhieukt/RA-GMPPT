#ifndef PIL_RUN_H
#define PIL_RUN_H
#include <stdint.h>
typedef struct {
    float Ebat, Pss, Ipk, tUnr, Vss, Ibss;
    int nRst;
    uint32_t cycMax, cycMean, cycMin;   /* controller cycles per call (target) */
    uint32_t cycMaxM[3], cycMeanM[3];   /* by mode: 0 search, 1 hold (incl. warm-up), 2 limiter */
} res_t;
/* timestamp hook: returns a free-running DOWN-counting cycle counter (target) or 0 (host) */
extern uint32_t (*pil_ts)(void);
void run_case(int alg, uint32_t abl, uint32_t seed, int pack, const float *Itab, int nTab, float dv,
              float soc0, float Tend, float Ilim, float noiseV, float noiseI, res_t *r, float *logP, int nLog);
#endif
