# Chromatin Clustering — Tethered Polymer Model with Stochastic Crosslinking

MATLAB implementation of a two-dimensional tethered bead–spring model of chromatin clustering, developed as a small, fully tracked member of the Bloom–Forest–Newhall family of chromatin polymer models (Vasquez 2016 → Hult 2017 → Walker 2019 → Coletti 2024). The current work asks how the timescale of crosslink turnover, relative to the polymer's own timescales, controls whether a chromatin segment settles into one cluster and keeps it, rearranges between clusters, or never clusters at all — and what sets the stability of a given cluster.

Author: Avery Zapata · Advisor: Dr. Katie Newhall · UNC Chapel Hill, Department of Mathematics

## Quick start

```matlab
% From the repository root in MATLAB (regime_band.m must be on the path):
section1_figs            % forces, timescales + measured bond lifetimes, collapse from a random start (~20 min)
model_schematic          % annotated N = 8 schematic of the model (instant)
fold_snapshots           % configuration snapshots per cluster size (~5 min)
state_occupancy_sweep    % main production sweep: 48 runs x 100 s at 14 lambda values
                         %   (~1.5 h from scratch on 4 workers; seconds if the .mat cache is present)
```

All drivers call `polymer_sim_2d.m`; the λ-axis figures also call `regime_band.m`. Figures are written to `figs/` as PNG; `.mat` results go to the working directory and are git-ignored. Requirements: base MATLAB; the Parallel Computing Toolbox is used by `state_occupancy_sweep` if present, otherwise it runs serially.

`state_occupancy_sweep` is cached: it loads `state_occupancy_N8.mat`, keeps every λ already there with the same protocol, simulates only λ values that are missing, and saves the merged set back (backing up the old file first). Adding a λ value to `lambda_list` therefore costs only that value. Runs are seeded by λ value, so existing points are never disturbed.

## Model (`polymer_sim_2d.m`)

### Setup

A chain of N beads (N = 8 in all production runs; the function default is 6) is confined to a disk of radius R_nuc = 175 nm with a hard wall. Beads 1 and N are pinned to opposite points on the boundary, (−R_nuc, 0) and (R_nuc, 0). Interior beads follow overdamped Langevin dynamics,

    ζ dx_i/dt = F_WLC + F_EV + F_cross + √(2 k_B T ζ) η_i(t),

integrated by Euler–Maruyama, with ζ = 2.5×10⁻³ pN·s/nm and k_BT = 4.1 pN·nm (D = k_BT/ζ ≈ 1640 nm²/s). Beads that step outside the disk are projected radially back onto it; the two end beads are re-pinned every step. The initial condition is a random self-avoiding walk between the pins (step 70 nm, minimum separation 15 nm).

### Forces

**Backbone (WLC).** Adjacent beads are joined by a worm-like-chain spring with zero rest length,

    F_WLC(r) = α ( −1 + 1/(1 − r/R₀)² + 4 r/R₀ ),    α = 0.2176 pN,  R₀ = N_k · 2 L_p = 1700 nm  (L_p = 50 nm, N_k = 17).

Within the disk r/R₀ < 0.2, so the spring is nearly linear with stiffness k_WLC = 6α/R₀ ≈ 7.7×10⁻⁴ pN/nm. The rms length of a backbone link is √(k_BT/k_WLC) ≈ 73 nm and the backbone relaxation time ζ/k_WLC ≈ 3.3 s.

**Excluded volume.** A soft Gaussian repulsion between every pair,

    F_EV(r) = c r e^{−a r²},    c = 8.305×10⁻⁵ pN/nm,  a = 3.268×10⁻⁵ nm⁻².

Its range 1/√a ≈ 175 nm and peak ≈ 6×10⁻³ pN; it is weaker than the WLC pull at every separation (6α/(R₀c) ≈ 9), so neighbour spacing is set thermally.

**Crosslinks.** A bonded pair is coupled by a linear spring of stiffness k_cross = 0.1 pN/nm and zero rest length. Bonded beads therefore sit ≈ √(k_BT/k_cross) ≈ 6 nm apart, and the spring pulls a pair together on the time τ_pull = ζ/k_cross = 25 ms.

### Crosslink kinetics

Each bead holds at most one crosslink (valency one). Every step, existing bonds are first released with probability k_off_eff·dt; then unbound, non-adjacent pairs (|i − j| ≥ 2) closer than r_elig = 90 nm bind with probability k_on_eff·dt·exp(−r²/σ²), σ = 60 nm, in pair order with no closest-pair priority. The 90 nm cutoff is the reach of one crosslinking complex (two ≈ 45 nm arms). Base rates are k_on = 1 s⁻¹ and k_off = 0.5 s⁻¹.

The control parameter λ_cross rescales both rates together:

    k_on_eff = k_on / λ_cross,    k_off_eff = k_off / λ_cross.

The binding affinity is unchanged; only the speed of the kinetics moves. Mean bond lifetime is λ_cross/k_off: 1 ms at λ_cross = 5×10⁻⁴, 0.16 s at 0.08. The measured mean lifetime matches λ_cross/k_off to 1–5 % at every sweep value. The simulator tightens dt so that max(k_on_eff, k_off_eff)·dt ≤ 0.05 (the drivers request 0.02) and rescales step counts to preserve total simulated time.

### Timescales

| quantity | value |
|---|---|
| time step (drivers) | min(10⁻³ s, 0.02 λ_cross) |
| bond lifetime λ_cross/k_off | 1 ms – 0.16 s over the sweep |
| crosslink pull-in ζ/k_cross | 25 ms |
| diffusion across r_thresh = 25 nm | ≈ 0.1 s |
| diffusion across σ = 60 nm | ≈ 0.5 s |
| backbone relaxation ζ/k_WLC | ≈ 3.3 s |
| collapse from a random chain | ≈ 50 s at λ = 10⁻³, ≈ 27 s at 0.015 |
| warm-up / production | 60 s / 100 s |

The bond lifetime crosses τ_pull at λ_cross = k_off·τ_pull = 0.0125. This is the only crossing between a kinetic and a mechanical timescale in the sweep, and it is where the rigid-to-flexible transition is observed.

## Analysis definitions

Vocabulary: **cluster** is the noun; **fold** is used only as a verb (the chain folds, refolds).

**Cluster.** Connected component of beads pairwise within r_thresh = 25 nm, of size ≥ 3. A configuration is labelled by the set of beads in the largest cluster; mirror images (bead i ↔ N+1−i) are pooled, so every label starting with "1" contains a pinned bead. Only a few clusters recur: the seven-bead cluster on one pin (1234567), the five-bead cluster on one pin with a three-bead loop at the other (12345), and two four-bead clusters, one at each pin (1234 + 5678). A six-bead cluster is a collapse in progress (234567 before bead 8 joins), not a settled configuration. *Pinned* sets contain a pinned bead; *interior* sets (2357, 23467, …) do not.

**State series.** Every 5 ms frame is labelled by its largest-cluster bead set, then smoothed by a 0.5 s sliding-window majority filter and a 0.5 s minimum-dwell rule so that threshold flicker on the bond-lifetime scale is not counted as switching. The raw per-frame labels are saved, so the smoothing can be changed without re-simulating. Because of the window, measured escape rates cannot exceed ≈ 1–2 s⁻¹; rates near that value are a floor of the analysis, not a measurement.

**Refold.** A change of the largest-cluster bead set to one that is neither a subset nor a superset of the previous set. This excludes a collapse completing, a pinned bead transiently joining, and (1234 + 5678) ↔ none flicker.

**Escape rate / dwell time.** Exits from a state divided by total time in that state, pooled over replicates; mean dwell = 1/rate. If no exit is observed, the 95 % upper bound 3/T is reported (a lower bound T/3 on the dwell). Poisson errors; Wilson intervals for fractions of runs.

## Entry points

| script | what it does | output |
|---|---|---|
| `polymer_sim_2d.m` | the simulator; returns trajectory, bond history, bond lifetimes, cluster snapshots | struct |
| `state_occupancy_sweep.m` | 48 runs × 100 s at 14 λ_cross values (5×10⁻⁴ … 0.08, dense in 5×10⁻³–0.03); occupancy by size class, escape rates, refold probability, dwell per bead set; cached | `state_occupancy_N8.mat`, `figs/fig26_time_occupancy`, `fig27_stability`, `fig28_refold_probability`, `fig30_escape_time_by_pattern` |
| `section1_figs.m` | forces and binding rate vs separation; timescales vs λ_cross with measured bond lifetimes (6 runs × 20 s per λ); crosslink count and largest cluster from a random start; prints the bound fraction per λ | `section1_figs.mat`, `figs/figI1_forces_vs_r`, `figI2_timescales`, `figI3_model_over_time` |
| `fold_snapshots.m` | 16 short runs at four λ_cross values; longest dwell in each size class drawn as chain + magnified cluster with the 25 nm proximity graph | `fold_snapshots_N8.mat`, `figs/snap_by_size` |
| `model_schematic.m` | annotated schematic of the N = 8 model | `figs/model_schematic` |
| `regime_band.m` | helper: draws the rigid / flexible / amorphic gradient behind any λ axis; boundaries (5×10⁻³, 2.5×10⁻²) and colours defined once here | — |
| `mixing_time_demo.m` | 10 × 1000 s runs at λ = 0.002 and 0.015 with no warm-up: fraction of runs that have refolded vs time, per-run occupancy vs time (frozen-basin demonstration; not yet run) | `mixing_time_N8.mat`, `figs/fig29_mixing_time` |
| `basin_census.m` | earlier whole-run census (superseded by `state_occupancy_sweep`) | `basin_census_N8.mat` |

`state_occupancy_sweep.m` is split into cells: CONFIG → CACHE → PREFLIGHT → RUNS → ANALYSIS → FIGURES. To re-analyse or redraw, `load('state_occupancy_N8.mat')` and run from the ANALYSIS or FIGURES cell.

All paper figures are drawn 6.5 in wide with 10 pt text (`paperfig` helper in each script) so that, included at `\textwidth`, figure text matches the caption size. Legends sit outside the axes; there are no figure titles (captions carry that information).

## Current findings (September 2026)

- **Rigid (λ_cross ≲ 5×10⁻³).** The cluster is decided during warm-up and never changes: zero refolds in 4 × 48 × 100 s of production. The seven-bead cluster's escape rate is < 1.2×10⁻³ s⁻¹ (no exit seen) at λ ≤ 2×10⁻³ and one exit in 3400 s at 5×10⁻³. Replicate spread is the mixture of clusters drawn at collapse, not sampling noise.
- **Transition.** The refold probability per 100 s rises from 0.04 at 6×10⁻³ to 0.81 at 0.02 and crosses ½ at λ_cross = 0.012, one grid point from the predicted k_off·τ_pull = 0.0125. The seven-bead dwell falls below the 100 s run at λ ≈ 8×10⁻³ (143 s → 47 s at 0.01 → 20 s at 0.012 → 5 s at 0.02).
- **Stability orders by size.** Escape rate decreases with cluster size at every λ: 7 > 5 > 4 beads. Over 5×10⁻³ → 0.02 the seven-bead rate rises ~650× while the four-bead rate rises ~3×: the transition is the collapse of the seven-bead cluster's lifetime, not a uniform speed-up. Clusters outlive bonds by > 5×10⁵ at λ = 10⁻³ and ~20 at 0.05.
- **Pinned vs interior.** Within a size class, a cluster that contains a pinned bead outlives interior clusters of the same size: 12345 vs 23467/23457 by 2–5×, 1234 vs 2357/2356 by 1.3–1.7×, at the same λ. Size and pinning together set stability; the bead count alone does not. Below λ ≈ 10⁻² only pinned sets accumulate enough time to be measured.
- **Occupancy lags the dynamics.** Time in the seven-bead cluster stays at 49–71 % from 5×10⁻⁴ to 0.015 while its lifetime falls a hundredfold, then drops (0.27 at 0.02, 0.01 at 0.05). Occupancy shows what forms, not how stable it is; regimes are located by escape rates and refold probability.
- **Bound fraction follows the cluster.** 0.61 ± 0.01 of beads are bound at every λ ≤ 0.01 (a twenty-fold range of bond lifetime), falling to 0.57 at 0.015–0.03 and 0.48–0.50 at 0.05–0.08 as the cluster is lost.
- **Amorphic (λ_cross ≳ 0.03).** No persistent cluster; unclustered 50 % of the time at 0.03, 70 % at 0.05, 81 % at 0.08; every escape rate sits at the ~1 s⁻¹ smoothing floor.

## Earlier scripts

The 1D build-up scripts (`One_Dimension_Tethered_Beads.m` → `Two_Dimensional_Tethered_Beads.m`), the binding-rate comparison scripts (`MatchingKons.m`, `ComparingKon*.m`), the sensitivity analyses (`AlphaSensitivityAnalysis.m`, `cSensivityAnalysis.m`, `Isolating_Parameter_Regime_WLC_and_EV.m`), the string-method utilities, and `TwoD_First_Pass.m` / `TimescaleFirstIntroduction.m` are retained for reference. They use earlier parameter values (k_cross = 0.01, σ = 20 nm, a wall-to-wall domain, and τ_cross for what is now λ_cross) and are not used for current results. The `1d/` directory shadows some root filenames; keep one or the other on the path.

## Workflow notes

- Generated `.mat`, `.fig`, `.avi`, `.mp4` and `figs/*.png` are git-ignored; `git add -f` to commit a canonical result. `state_occupancy_N8.mat` is the canonical sweep and is worth committing with `-f` — regenerating it from scratch takes hours.
- Tag a version before any breaking change. Checkpoints: `v0.2-tau-cross` (first verified timescale-separation implementation); [tag the 14-λ sweep before the paper submission].

## References

Primary modelling sources (UNC / Bloom–Forest group)

- Vasquez, P. A. et al. (2016). Entropy gives rise to topologically associating domains. *Nucleic Acids Research*, 44(12), 5540–5549.
- Hult, C. et al. (2017). Enrichment of dynamic chromosomal crosslinks drive phase separation of the nucleolus. *Nucleic Acids Research*, 45(19), 11159–11173.
- Walker, B. et al. (2019). Transient crosslinking kinetics optimize gene cluster interactions. *PLOS Computational Biology*, 15(8), e1007124.
- Walker, B. L. (2021). *Emergent Structure and Dynamics from Stochastic Pairwise Crosslinking in Chromosomal Polymer Models*. PhD thesis, UNC Chapel Hill.
- He, Y. et al. (2020). Statistical mechanics of chromosomes. *Nucleic Acids Research*, 48(20), 11284–11303.
- Kolbin, D. et al. (2023). Polymer modeling reveals interplay between physical properties of chromosomal DNA and the size and distribution of condensin-based chromatin loops. *Genes*, 14(12), 2193.
- Coletti, A., Newhall, K. A., Walker, B. L., & Bloom, K. (2024). Different relative scalings between transient forces and thermal fluctuations tune regimes of chromatin organization. arXiv:2401.06921.

Biology of SMC complexes

- Holmes, V. F., & Cozzarelli, N. R. (2000). Closing the ring: links between SMC proteins and chromosome partitioning, condensation, and supercoiling. *PNAS*, 97(4), 1322–1324.
- Eeftens, J. M. et al. (2016). Condensin Smc2–Smc4 dimers are flexible and dynamic. *Cell Reports*, 14(8), 1813–1818.

Other chromatin polymer models

- Bohn, M., & Heermann, D. W. (2010). Diffusion-driven looping provides a consistent framework for chromatin organization. *PLOS ONE*, 5(8), e12218.
- Barbieri, M. et al. (2012). Complexity of chromatin folding is captured by the strings and binders switch model. *PNAS*, 109(40), 16173–16178.
- Brackley, C. A. et al. (2013). Nonspecific bridging-induced attraction drives clustering of DNA-binding proteins and genome organization. *PNAS*, 110(38), E3605–E3611.
- Tjong, H. et al. (2012). Physical tethering and volume exclusion determine higher-order genome organization in budding yeast. *Genome Research*, 22, 1295–1305.

## License

See LICENSE.txt.
