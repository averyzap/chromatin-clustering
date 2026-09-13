# Chromatin Clustering — Tethered Polymer Model with Stochastic Crosslinking

MATLAB implementation of a two-dimensional tethered bead–spring model of chromatin clustering, developed as a simpler and more biologically grounded alternative to the quasi-string approach of Vasquez et al. (2016). The current work asks how the timescale of crosslink turnover, relative to polymer relaxation, controls whether a chromatin segment folds into a fixed arrangement or keeps rearranging.

Author: Avery Zapata · Advisor: Dr. Katie Newhall · UNC Chapel Hill, Department of Mathematics

## Quick start

```matlab
% From the repository root in MATLAB:
section1_figures         % analytic force/timescale figures + short runs (~5 min)
fold_snapshots           % configuration snapshots per cluster size / fold (~8 min)
state_occupancy_sweep    % main production sweep, 48 runs x 100 s at 10 lambda values (~1 h on 4 workers, ~3 h serial)
```

All three call `polymer_sim_2d.m`. Output figures go to `figs/`; `.mat` results are written to the working directory and are git-ignored by default. Requirements: base MATLAB; Parallel Computing Toolbox is used by `state_occupancy_sweep` if present, otherwise it runs serially.

## Model (`polymer_sim_2d.m`)

### Setup

A chain of N beads (N = 8 in all production runs; the function default is 6) is confined to a disk of radius R_nuc = 175 nm with a hard wall. Beads 1 and N are pinned to opposite points on the boundary, (−R_nuc, 0) and (R_nuc, 0). Interior beads follow overdamped Langevin dynamics,

    ζ dx_i/dt = F_WLC + F_EV + F_cross + √(2 k_B T ζ) η_i(t),

integrated by Euler–Maruyama, with ζ = 2.5×10⁻³ pN·s/nm and k_BT = 4.1 pN·nm (D = k_BT/ζ ≈ 1640 nm²/s). Beads that step outside the disk are projected radially back onto it; the two end beads are re-pinned every step.

### Forces

**Backbone (WLC).** Adjacent beads are joined by a worm-like-chain spring with zero rest length,

    F_WLC(r) = α ( −1 + 1/(1 − r/R₀)² + 4 r/R₀ ),    α = 0.2176 pN,  R₀ = N_k · 2 L_p = 1700 nm  (L_p = 50 nm, N_k = 17).

Within the disk r/R₀ < 0.2, so the spring is nearly linear with stiffness k_WLC = 6α/R₀ ≈ 7.7×10⁻⁴ pN/nm. The rms length of a backbone link is √(k_BT/k_WLC) ≈ 73 nm and the backbone relaxation time ζ/k_WLC ≈ 3.3 s.

**Excluded volume.** A soft Gaussian repulsion between every pair,

    F_EV(r) = c r e^{−a r²},    c = 8.305×10⁻⁵ pN/nm,  a = 3.268×10⁻⁵ nm⁻².

Its range 1/√a ≈ 175 nm and peak ≈ 6×10⁻³ pN; it is weaker than the WLC pull at every separation, so neighbour spacing is set thermally, not by a WLC/EV balance.

**Crosslinks.** A bonded pair is coupled by a linear spring of stiffness k_cross = 0.1 pN/nm and zero rest length. Bonded beads therefore sit ≈ √(k_BT/k_cross) ≈ 6 nm apart, and the spring pulls a pair together on the time τ_pull = ζ/k_cross = 25 ms.

### Crosslink kinetics

Each bead holds at most one crosslink (valency one). Every step, existing bonds are first released with probability k_off_eff·dt; then unbound, non-adjacent pairs (|i − j| ≥ 2) closer than r_elig = 90 nm bind with probability k_on_eff·dt·exp(−r²/σ²), σ = 60 nm, in pair order with no closest-pair priority. Base rates are k_on = 1 s⁻¹ and k_off = 0.5 s⁻¹.

The control parameter λ_cross rescales both rates together:

    k_on_eff = k_on / λ_cross,    k_off_eff = k_off / λ_cross.

The equilibrium bond fraction is unchanged; only the crosslink timescale moves. Mean bond lifetime is λ_cross/k_off: 1 ms at λ_cross = 5×10⁻⁴, 0.1 s at 0.05. The simulator tightens dt so that max(k_on_eff, k_off_eff)·dt ≤ 0.05 (the drivers request 0.02) and rescales step counts to preserve total simulated time.

### Timescales

| quantity | value |
|---|---|
| time step (drivers) | min(10⁻³ s, 0.02 λ_cross) |
| bond lifetime λ_cross/k_off | 1 ms – 0.1 s over the sweep |
| crosslink pull-in ζ/k_cross | 25 ms |
| diffusion across r_thresh = 25 nm | ≈ 0.1 s |
| diffusion across σ = 60 nm | ≈ 0.5 s |
| backbone relaxation ζ/k_WLC | ≈ 3.3 s |
| warm-up / production | 20 s / 100 s |

The bond lifetime crosses τ_pull at λ_cross = k_off·τ_pull = 0.0125, which is where the rigid-to-flexible transition is observed.

## Analysis definitions

**Cluster.** Connected component of beads pairwise within r_thresh = 25 nm, of size ≥ 3.

**Fold.** The chain's spatial arrangement, labelled by the set of beads in the largest cluster; mirror images (bead i ↔ N+1−i) are pooled. Only a few folds occur: the big fold (seven beads on one pin, 2345678/1234567), the five-bead fold (45678/12345), and the two-fold (1234 + 5678). A six-bead cluster is a collapse in progress.

**State series.** Every 5 ms frame is labelled by its largest-cluster bead set, then smoothed by a 0.5 s sliding-window majority filter and a 0.5 s minimum-dwell rule so that threshold flicker on the bond-lifetime scale is not counted as switching. The raw per-frame labels are saved, so the smoothing can be changed without re-simulating.

**Refold.** A change of the largest-cluster bead set to one that is neither a subset nor a superset of the previous set. This excludes a collapse completing, a pinned bead transiently joining, and two-fold ↔ none flicker.

**Escape rate.** Exits from a state divided by total time in that state, pooled over replicates. If no exit is observed, the 95% upper bound 3/T is reported. Poisson errors; Wilson intervals for fractions of runs.

## Entry points

| script | what it does | output |
|---|---|---|
| `polymer_sim_2d.m` | the simulator; returns trajectory, bond history, bond lifetimes, cluster snapshots | struct |
| `state_occupancy_sweep.m` | 48 runs × 100 s at 10 λ_cross values; time in each cluster-size class, escape rates, refold probability | `state_occupancy_N8.mat`, `figs/fig26–28` |
| `fold_snapshots.m` | short runs at four λ_cross values; longest dwell in each size class / named fold drawn as a snapshot | `fold_snapshots_N8.mat`, `figs/snap_by_*` |
| `section1_figures.m` | forces vs separation, timescales vs λ_cross, collapse from a random start, bond lifetime and bound fraction vs λ_cross | `section1_figs.mat`, `figs/figI1–I4` |
| `basin_census.m` | earlier whole-run fold census (superseded by `state_occupancy_sweep`) | `basin_census_N8.mat` |

`state_occupancy_sweep.m` is split into cells: CONFIG → PREFLIGHT → RUNS → ANALYSIS → FIGURES. To re-analyse or redraw, `load('state_occupancy_N8.mat')` and run from the ANALYSIS or FIGURES cell.

## Current findings (September 2026)

- Below λ_cross ≈ 5×10⁻³ the fold is decided during warm-up and never changes: 0–1 rearrangements in 4800 s of production per λ_cross; the seven-bead fold's escape rate is < 10⁻³ s⁻¹.
- Between 7×10⁻³ and 0.02 the refold probability per 100 s rises from 0.06 to 0.79 and the seven-bead escape rate rises ~300×; the crossover sits at λ_cross ≈ k_off τ_pull.
- At every λ_cross, escape rate decreases with cluster size: 7 > 5 > 4 beads. The transition is a collapse of the big fold's lifetime, not a uniform speed-up.
- Time-averaged occupancy of the seven-bead fold stays at 50–70% up to λ_cross = 0.015 and only then falls; occupancy lags the dynamics by a decade.
- Above λ_cross ≈ 0.03 no fold persists; the chain is unclustered > 50% of the time.
- Replicate spread at small λ_cross is the mixture of folds drawn, not sampling noise, and does not shrink with run length.

## Earlier scripts

The 1D build-up scripts (`One_Dimension_Tethered_Beads.m` → `Two_Dimensional_Tethered_Beads.m`), the kernel-comparison scripts (`MatchingKons.m`, `ComparingKon*.m`), the sensitivity analyses (`AlphaSensitivityAnalysis.m`, `cSensivityAnalysis.m`, `Isolating_Parameter_Regime_WLC_and_EV.m`), the string-method utilities, and `TwoD_First_Pass.m` / `TimescaleFirstIntroduction.m` are retained for reference. They use earlier parameter values (k_cross = 0.01, σ = 20 nm, a wall-to-wall domain, and τ_cross for what is now λ_cross) and are not used for current results. The `1d/` directory shadows some root filenames; keep one or the other on the path.

## Workflow notes

- Generated `.mat`, `.fig`, `.avi`, `.mp4` and `figs/*.png` are git-ignored; `git add -f` to commit a canonical result.
- Tag a version before any breaking change. Checkpoints: `v0.2-tau-cross` (first verified timescale-separation implementation); [tag the state_occupancy sweep before the paper submission].

## References

Primary modelling sources (UNC / Bloom–Forest group)

- Vasquez, P. A. et al. (2016). Entropy gives rise to topologically associating domains. *Nucleic Acids Research*, 44(12), 5540–5549.
- Hult, C. et al. (2017). Enrichment of dynamic chromosomal crosslinks drive phase separation of the nucleolus. *Nucleic Acids Research*, 45(19), 11159–11173.
- Walker, B. et al. (2019). Toward a dynamic model of the nucleolus: a stochastic model of transient gene–gene crosslinking. *PLOS Computational Biology*, 15(8), e1007124.
- He, Y. et al. (2020). Statistical mechanics of chromosomes. *Nucleic Acids Research*, 48(20), 11284–11303.
- Kolbin, D. et al. (2023). Polymer modeling reveals interplay between physical properties of chromosomal DNA and the size and distribution of condensin-based chromatin loops. *Genes*, 14(12), 2193.
- Coletti, A., Newhall, K. A., Walker, B. L., & Bloom, K. (2024). Different relative scalings between transient forces and thermal fluctuations tune regimes of chromatin organization. arXiv:2401.06921.
- Walker, B. L. *Emergent Structure and Dynamics from Stochastic Pairwise Crosslinking in Chromosomal Polymer Models* (thesis).

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
