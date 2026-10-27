/* Host build of the PIL code (gcc): checks the C controller against the Python/MATLAB implementation. */
#include <stdio.h>
#include <stdlib.h>
#include "pil_run.h"
#include "pil_tables.h"
int main(int argc, char **argv)
{
    float nV = argc > 1 ? atof(argv[1]) : 0.05f, nI = argc > 2 ? atof(argv[2]) : 0.005f;
    struct { const char *nm; int alg; unsigned abl; int pack; const float *tab; float soc; float Ilim; } S[] = {
        {"DPSO 72-C4-S60",5,0,72,ITAB_C4,0.6f,10}, {"RA-DPSO 72-C4-S60",6,0,72,ITAB_C4,0.6f,10},
        {"DPSO 72-C8-S60",5,0,72,ITAB_C8,0.6f,10}, {"RA-DPSO 72-C8-S60",6,0,72,ITAB_C8,0.6f,10},
        {"DPSO 48-C1-S60",5,0,48,ITAB_C1,0.6f,10}, {"RA-DPSO 48-C1-S60",6,0,48,ITAB_C1,0.6f,10},
        {"RA-DPSO-F 48-C1-S60",6,2,48,ITAB_C1,0.6f,10},
        {"CS 72-C4-S60",7,0,72,ITAB_C4,0.6f,10}, {"RA-CS 72-C4-S60",8,0,72,ITAB_C4,0.6f,10},
        {"CS 72-C8-S60",7,0,72,ITAB_C8,0.6f,10}, {"RA-CS 72-C8-S60",8,0,72,ITAB_C8,0.6f,10}};
    int i, sd;
    for (i = 0; i < (int)(sizeof S / sizeof S[0]); i++) {
        int ns = S[i].alg >= 7 ? 5 : 1;
        for (sd = 1; sd <= ns; sd++) {
            res_t r;
            run_case(S[i].alg, S[i].abl, sd, S[i].pack, S[i].tab, NTAB, DVTAB, S[i].soc, 0.8f, S[i].Ilim, nV, nI, &r, 0, 0);
            printf("%-22s seed %d  Pss %7.2f W  Vss %6.2f V  Ebat %7.2f J  Ipk %6.2f A  tUnr %6.1f ms  nRst %d\n",
                   S[i].nm, sd, r.Pss, r.Vss, r.Ebat, r.Ipk, 1e3f * r.tUnr, r.nRst);
        }
    }
    return 0;
}
