/* =========================================================================
 * PIL on TMS320F28335: the controller (pil_ctrl.c) runs on the target at a 20 kHz control rate,
 * closed loop with the averaged plant model (pil_plant.c) co-simulated on the same processor.
 * CPU Timer 1 (SYSCLKOUT = 150 MHz, down-counting) measures the cycles of every controller call.
 * Results: CCS console (printf) and the array RES[] (Expressions window).
 * Build: C2000Ware DSP2833x headers + common (InitSysCtrl, InitPieCtrl...), --float_support=fpu32,
 *        -O2, linker cmd pil_28335_RAM.cmd + DSP2833x_Headers_nonBIOS.cmd; stack 0x380, heap 0x400.
 *        If .text does not fit in RAM: define NO_PRINTF (results then only in RES[], RA_CYC, TS_OVH).
 * ========================================================================= */
#include "DSP2833x_Device.h"
#include "DSP2833x_Examples.h"
#ifndef NO_PRINTF
#include <stdio.h>
#endif
#include "pil_run.h"
#include "pil_ctrl.h"
#include "pil_tables.h"

static uint32_t ts_timer1(void) { return CpuTimer1Regs.TIM.all; }

typedef struct { const char *nm; int alg; uint32_t abl; int pack; const float *tab; float soc; } case_t;
static const case_t CASES[] = {
    {"DPSO 72-C4-S60",5,0,72,ITAB_C4,0.6f}, {"RA-DPSO 72-C4-S60",6,0,72,ITAB_C4,0.6f},
    {"DPSO 72-C8-S60",5,0,72,ITAB_C8,0.6f}, {"RA-DPSO 72-C8-S60",6,0,72,ITAB_C8,0.6f},
    {"DPSO 48-C1-S60",5,0,48,ITAB_C1,0.6f}, {"RA-DPSO 48-C1-S60",6,0,48,ITAB_C1,0.6f},
    {"RA-DPSO-F 48-C1-S60",6,2,48,ITAB_C1,0.6f},
    {"CS 72-C4-S60",7,0,72,ITAB_C4,0.6f}, {"RA-CS 72-C4-S60",8,0,72,ITAB_C4,0.6f},
    {"CS 72-C8-S60",7,0,72,ITAB_C8,0.6f}, {"RA-CS 72-C8-S60",8,0,72,ITAB_C8,0.6f} };
#define NCASE (sizeof(CASES)/sizeof(CASES[0]))
volatile res_t RES[NCASE];
volatile uint32_t RA_CYC = 0, TS_OVH = 0, DIV_CYC = 0, MUL_CYC = 0;
volatile float src = 3.0f;
volatile int DONE = 0;
volatile float sink;

void main(void)
{
    int i, k; uint32_t t0, t1;
    InitSysCtrl();
    DINT; InitPieCtrl(); IER = 0x0000; IFR = 0x0000; InitPieVectTable();
    CpuTimer1Regs.TCR.bit.TSS = 1;            /* stop */
    CpuTimer1Regs.PRD.all = 0xFFFFFFFF;       /* free running, 1 count per SYSCLKOUT cycle */
    CpuTimer1Regs.TPR.all = 0; CpuTimer1Regs.TPRH.all = 0;
    CpuTimer1Regs.TCR.bit.TRB = 1;            /* reload */
    CpuTimer1Regs.TCR.bit.TSS = 0;            /* start */
    pil_ts = ts_timer1;
    /* timestamp overhead */
    t0 = ts_timer1(); t1 = ts_timer1(); TS_OVH = t0 - t1;
    /* micro-benchmark of the RA layer (eq. 4 + projection), 1000 calls */
    t0 = ts_timer1();
    for (k = 0; k < 1000; k++) sink = ra_project(78.0f + 0.001f * k, 0.55f);
    t1 = ts_timer1(); RA_CYC = (t0 - t1) / 1000;
    /* cost of one float division vs. one multiplication (loop of 1000, same loop body otherwise) */
    t0 = ts_timer1();
    for (k = 0; k < 1000; k++) sink = src / (1.0f + 0.001f * k);
    t1 = ts_timer1(); DIV_CYC = (t0 - t1) / 1000;
    t0 = ts_timer1();
    for (k = 0; k < 1000; k++) sink = src * (1.0f + 0.001f * k);
    t1 = ts_timer1(); MUL_CYC = (t0 - t1) / 1000;
    for (i = 0; i < (int)NCASE; i++) {
        res_t r;
        run_case(CASES[i].alg, CASES[i].abl, 1, CASES[i].pack, CASES[i].tab, NTAB, DVTAB, CASES[i].soc, 0.8f, 10.0f,
                 0.05f, 0.005f, &r, 0, 0);
        r.cycMax -= TS_OVH; r.cycMean -= TS_OVH; r.cycMin -= TS_OVH;
        for (k = 0; k < 3; k++) { if (r.cycMaxM[k]) r.cycMaxM[k] -= TS_OVH; if (r.cycMeanM[k]) r.cycMeanM[k] -= TS_OVH; }
        RES[i] = r;
#ifndef NO_PRINTF
        printf("%-22s Pss %7.2f W Ebat %7.2f J Ipk %6.2f A | cycles mean %lu max %lu (%.2f / %.2f us)\n",
               CASES[i].nm, r.Pss, r.Ebat, r.Ipk, (unsigned long)r.cycMean, (unsigned long)r.cycMax,
               r.cycMean / 150.0f, r.cycMax / 150.0f);
        printf("    by mode (mean/max cycles): search %lu/%lu  hold %lu/%lu  limiter %lu/%lu\n",
               (unsigned long)r.cycMeanM[0], (unsigned long)r.cycMaxM[0], (unsigned long)r.cycMeanM[1], (unsigned long)r.cycMaxM[1],
               (unsigned long)r.cycMeanM[2], (unsigned long)r.cycMaxM[2]);
#endif
    }
#ifndef NO_PRINTF
    printf("RA layer (eq. 4 + R1): %lu cycles = %.3f us; timestamp overhead %lu cycles\n",
           (unsigned long)RA_CYC, RA_CYC / 150.0f, (unsigned long)TS_OVH);
#endif
    DONE = 1;
    for (;;) { }
}
