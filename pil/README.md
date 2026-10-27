# Processor-in-the-loop (PIL) code for the TMS320F28335

The controller is written in C with single-precision arithmetic. It consists of:
- the RA layer;
- the DPSO and CS cores;
- the left-flank and conventional CC/CV limiters;
- the PV-voltage loop.

The controller runs on the F28335 in a closed loop with the averaged plant: PV string, buck converter and battery. The plant is co-simulated on the same processor between the controller calls. CPU timer 1 counts the clock cycles of every controller call. It counts them separately for the search, hold and limiter modes. The plant time is not included.

What the test does and does not cover:
- It verifies the computation time and the numerical equivalence with the PC build.
- It does not cover real-time I/O: ADC quantization, delays and real measurement noise are not part of the test.

## Files

| File | Content | In the CCS project |
|---|---|---|
| `pil_ctrl.c/.h` | controller (equivalent to `matlab/Charger_Bench.m`, alg 5–8) | yes |
| `pil_plant.c/.h` | averaged plant, 5 µs step | yes |
| `pil_run.c/.h` | closed loop, indices, cycle measurement | yes |
| `pil_tables.c/.h` | I(V) tables of the PV string (0.1 V grid) | yes |
| `target_main.c` | F28335 main program (11 runs, micro-benchmarks) | yes |
| `pil_28335_RAM.cmd` | linker file: code and data in RAM, PV tables in Flash | yes |
| `host_main.c` | the same runs on a PC | no |
| `lt_main.c`, `lt_ramp.c` | one-hour charging profiles on a PC (Section 5.9) | no |

## PC build

```bash
gcc -O2 -I. host_main.c pil_run.c pil_ctrl.c pil_plant.c pil_tables.c -lm -o host && ./host
gcc -O2 -I. lt_main.c pil_ctrl.c -lm -o lt     # ./lt alg abl seed tend soc0 out.csv [steady]
```

## F28335 build (Code Composer Studio 9.0, C2000 CGT 18.12.1.LTS)

1. Create an empty CCS project for the TMS320F28335 and add the six `.c` files, the headers and `pil_28335_RAM.cmd`.
2. Add from C2000Ware, `device_support/f2833x/`:
   - `headers/source/DSP2833x_GlobalVariableDefs.c`
   - `headers/cmd/DSP2833x_Headers_nonBIOS.cmd`
   - `common/source/`: `DSP2833x_SysCtrl.c`, `DSP2833x_PieCtrl.c`, `DSP2833x_PieVect.c`, `DSP2833x_DefaultIsr.c`, `DSP2833x_CodeStartBranch.asm`, `DSP2833x_usDelay.asm` and `DSP2833x_ADC_cal.asm`
3. Include paths: `headers/include` and `common/include`.
4. Compiler options: `--float_support=fpu32`, `-O2`, and the predefined symbol `NO_PRINTF`. The results are then read from memory; C I/O through the debug probe can stall the program.
5. Build, load with the debug probe, run, and wait until `DONE = 1`. Read the results in the Expressions view:
   - `RES[0..10]`: P_ss, E_bat, I_pk and the cycles per call (mean, max, min, and by mode);
   - `RA_CYC`, `TS_OVH`, `DIV_CYC`, `MUL_CYC`.
   Alternatively, save the memory block of `RES` with *Save Memory* in the *16-Bit Hex - TI Style* format. Use the address of `RES` and 512 words.

## Notes

- The F28335 has no hardware divider. A single-precision division costs about 245 cycles, including the loop overhead, against 11 cycles for a multiplication. For this reason, all constants that depend only on the sampling period are computed once in `ctrl_init`.
- The results on the F28335 agree with the PC build to the printed digits in all eleven runs.
- Measured cycles per controller call (mean / max):
  - DPSO-based controllers: about 2100–2850 / 3700;
  - CS-based controllers: about 2100–2150 / 10 400. The maximum occurs in the call that draws the Lévy flights of a new generation.
