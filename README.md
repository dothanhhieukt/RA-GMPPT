# RA-GMPPT: Reachability-Aware Global MPPT for Buck-Type PV–Battery Chargers

Simulation and processor-in-the-loop (PIL) code for the paper

> T. H. Do, "A Reachability-Aware Layer for Global MPPT in Buck-Type PV–Battery Chargers Under Partial Shading: A Simulation and Processor-in-the-Loop Study," *Electronics* (under review), 2026.

A buck converter can impose only PV voltages above about `V_bat/D_max`. Under partial shading the global maximum power point is often unreachable. The **reachability-aware (RA) layer** wraps a GMPPT core and has two parts:
- **R1**: projects every candidate onto the reachable range and makes the bound an explicit candidate.
- **R2**: restarts the search when the bound drifts while it is being held.

The repository covers:
- the switched Simulink model;
- an averaged model in Python;
- the C implementation that was run on a TMS320F28335;
- the scripts that produce the tables and figures.

## Repository layout

| Folder | Content |
|---|---|
| `matlab/` | Switched model (MATLAB/Simulink R2024a, Specialized Power Systems): controller `Charger_Bench.m` (MATLAB Function block), batch runner `run_bench.m`, test configuration `bench_config.m`, indices `bench_metrics.m`, battery model `battery_35E.m`, Simulink model `Charger_Buck.slx` |
| `averaged/` | Averaged model in Python: plant `plant.py`, controller `ctrl.py` (line-by-line equivalent of `Charger_Bench.m`), simulator `sim.py`, indices `metrics.py`, batch helper `batch.py`, plus the experiment and analysis scripts listed below |
| `pil/` | C implementation of the controller (`pil_ctrl.c`), averaged plant and test harness for the TMS320F28335 (`target_main.c`, linker file), PC builds (`host_main.c`, one-hour profiles `lt_main.c`, `lt_ramp.c`); see `pil/README.md` |
| `figures/` | `make_figures.py` (all figures of the paper) and the small data files it needs |

Raw results (switched-model result tables and waveforms, averaged-model CSV files, one-hour logs, PIL memory dumps) are archived on Zenodo: https://doi.org/10.5281/zenodo.23684746. Unpack the archive into `data/`, or set the environment variable `RA_DATA` to its location.

## Requirements

- MATLAB/Simulink R2024a with Simscape Electrical (Specialized Power Systems). The Parallel Computing Toolbox is optional; `run_bench` uses `parsim` when it is available.
- Python ≥ 3.10 with the packages in `requirements.txt`.
- A C compiler (gcc) for the PC builds. Code Composer Studio 9.0 with TI C2000 code-generation tools 18.12.1.LTS and C2000Ware for the F28335 build.

## Getting the code

```bash
git clone https://github.com/dothanhhieukt/RA-GMPPT.git
```

## Quick start

**Switched model (MATLAB):**
```matlab
cd matlab
run_bench('sweep')      % in-circuit P-V sweeps -> reference P_reach (run first)
run_bench('test')       % short check: 5 controllers x 2 scenarios
run_bench('static')     % static scenarios (see bench_config.m for all groups)
```

**Averaged model (Python):**
```bash
cd averaged
pip install -r ../requirements.txt
python -c "import sim; d = sim.run(6, 0, 1, 72, 'C4', Tend=0.8); print(d['P'][-2000:].mean())"
python camp_paired.py static      # 30 paired seeds, static MPPT scenarios (Tables 10, 11)
```

**PIL code on a PC:**
```bash
cd pil
gcc -O2 -I. host_main.c pil_run.c pil_ctrl.c pil_plant.c pil_tables.c -lm -o host && ./host
```

## Controller selection

`alg` selects the core:

| alg | Controller |
|---|---|
| 1 | P&O |
| 5 | DPSO |
| 6 | RA-DPSO |
| 7 | CS |
| 8 | RA-CS |
| 9 | SCAN |
| 10 | RA-SCAN |

`abl` is a bitmask of variants. Some examples:

| abl | Variant |
|---|---|
| 2 | fast left-flank limiter (extension F) |
| 8 | modified DPSO restart rule |
| 16 | RA layer off |
| 32 | conventional limiter |
| 64 | R2 off |
| 128 | hold tracking off |
| 256 | adaptive search range |
| 1024 | 20 kHz sampling |
| 4096·k | row k of the sensitivity table |

The header of `matlab/Charger_Bench.m` documents every bit.

**Irradiance cases:** C1–C4 correspond to Cases 1–4 of the paper. C6, C7 and C8 correspond to Cases 5, 6 and 7. C5 is not used in the paper.

**Scenario names:** a name such as `72-C4-S60` means the 72 V pack, case C4 and an initial SOC of 60 %.

## Paper → code

| Paper | Model | Script |
|---|---|---|
| Table 5 (reachable maximum) | switched | `run_bench('sweep')` |
| Table 8, 12, 14, A2, E1 | switched | `run_bench('static' / 'ablLim' / 'vrise' / 'dynamic' / 'cv')`, groups `adaptRange`, `vriseAR`, `limiter`, `dynLimiter`, `sens`, `fs20k` |
| Table 9 (adaptive range, SCAN) | averaged | `exp_adaptive_range.py` |
| Tables 10, 11 (CS, paired statistics) | averaged | `camp_paired.py static`, `an_static.py` |
| Table 13 (irradiance steps) | averaged | `camp_paired.py dyn`, `an_dyn.py`, `exp_irradiance_steps.py` |
| Section 5.6 (rising battery voltage) | averaged | `exp_vrise.py`, `camp_paired.py vrise` |
| Section 5.8, Table B1 (sensitivity) | averaged | `exp_sensitivity_det.py`, `camp_sens.py`, `an_sens.py`, `exp_sensor_offset.py` |
| Tables 15, 16, Figure 10 (one-hour profiles) | C controller + averaged plant | `pil/lt_main.c`, `pil/lt_ramp.c`, `averaged/an_lt.py` (Python reference: `longterm.py`) |
| Section 5.10, Table C1, Figure 11 (small-signal) | averaged | `bode.py`, `bode_tol.py`, `step_pm.py`, `exp_damping.py`, `exp_sampling_20k.py` |
| Table 17, Figure 12 (Monte Carlo) | static curve | `mc_static.py`, `mc_rskip.py`, `mc_dist.py` |
| Table 18 (PIL) | TMS320F28335 | `pil/target_main.c` |
| Table A1 (averaged vs. switched) | both | `validate_vs_switched.py` |
| Table D1 (limiter variants) | averaged | `exp_limiters.py`, `exp_fast_limiter.py` |
| Table D2 (ablation) | averaged | `exp_r1_only.py`, `camp_paired.py` |
| Restart rule, noise level, settling time | averaged | `exp_restart_rule.py`, `exp_noise.py`, `exp_tset.py` |
| Figures | – | `figures/make_figures.py` |

## Notes

- The stochastic core (CS) seeds the MATLAB random stream with `rng(seed0*1000 + 7)` at its first call (`seed0` = run number), so runs with the same seed are paired across controllers.
- `ctrl.py` and `pil_ctrl.c` reproduce `Charger_Bench.m` line by line for alg 1 and 5–10. The C build on the F28335 reproduces the PC build of the same sources.

## License

MIT (see `LICENSE`).

## Citation

See `CITATION.cff`.
