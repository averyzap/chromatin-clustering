# Chromatin Clustering — Quasi-String Reproduction & Dynamic Model

MATLAB implementation of a tethered polymer model of chromatin clustering dynamics, developed as a simpler and more biologically grounded alternative to the quasi-string approach of Vasquez et al. (2016). The model is built up incrementally from 1D tethered beads through worm-like-chain elasticity, excluded volume, and dynamic crosslinking, and reproduces qualitative cluster-formation behavior observed in prior UNC/Bloom–Forest group models.

**Author:** Avery Zapata
**Advisor:** Dr. Katie Newhall
**Affiliation:** UNC Chapel Hill, Department of Mathematics

---

## Quick Start

```matlab
% From the repository root in MATLAB:
TwoD_First_Pass               % baseline 2D simulation with crosslinking
TimescaleFirstIntroduction    % 2D simulation with tau_cross timescale separation
```

Output figures and `.mat` files are written to the working directory (and are git-ignored by default).

**Requirements:** MATLAB <<< R20XXx or later >>>, no required toolboxes beyond base MATLAB.

---

## Model

### Setup

A linear chain of `N` beads is tethered between two fixed walls separated by length `L`. Beads evolve under overdamped Langevin dynamics:

```
ζ dx_i/dt = F_WLC + F_EV + F_crosslink + sqrt(2 k_B T ζ) η(t)
```

where `η(t)` is Gaussian white noise. The bond network between beads forms and breaks stochastically.

### Forces

**Worm-like chain (WLC) connectivity.** Each adjacent bead pair is connected by a WLC spring with persistence length `Lp = 50 nm` and Kuhn segment count `Nk = 17`, giving contour length `R0 = Nk · 2·Lp` per segment. The WLC force law is approximated using the interpolation factor `α = 0.2176`.

**Excluded volume (EV).** A short-range repulsion between all bead pairs with strength `cEV = 8.305e-5` and decay parameter `aEV = 3.268e-5`. This term sets the soft-core radius beyond which beads can pass through each other; values below the critical `c*` lead to collapse, while values above stabilize clustering.

**Dynamic crosslinking.** Bead pairs within range bind and unbind stochastically with a single-bond-per-bead constraint. Bound pairs contribute a harmonic spring force with stiffness `k_cross = 0.01`. Unbinding is a Poisson process with rate `koff = 0.5`. Binding rate `kon(r)` is defined below.

### Binding kernel kon(r)

Two functional forms are implemented and matched:

1. **Caitlin (Gaussian, reference):** `kon(r) = kon0 · exp(−r² / σ²)` with `σ = 20 nm`.
2. **Ben (logistic):** `kon(r) = scale · 2 / (1 + exp(steep · (r/r_half − 1)))`.

The logistic kernel parameters are matched to the Gaussian at a chosen radius `r_match` by adjusting `scale` so both kernels agree there. Current production values:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `r_match` | 31 nm | radius at which both kernels agree |
| `r_half`  | 21.5 nm | midpoint of the logistic drop-off |
| `steep`   | 5 | steepness of the logistic drop-off |
| `kon0`    | 1.0 | overall rate scale |

These were chosen so that the logistic kernel reaches ~0.91 at `r = 0` (vs. 0.37 for the Gaussian at the same point), making the Ben kernel substantially more binding-favorable at short range. See `MatchingKons.m` and `ComparingKon_With_Proper_BensKon.m`.

### Timescale separation (`tau_cross`)

To reproduce Anna's separation-of-timescales formulation, a single parameter `tau_cross` rescales both binding and unbinding rates:

```
kon_eff  = kon0 / tau_cross
koff_eff = koff / tau_cross
```

Because both rates scale identically, the equilibrium bound fraction `kon/(kon+koff)` is preserved while the *timescale* of crosslink dynamics relative to polymer relaxation changes:

- `tau_cross < 1` — crosslinks faster than polymer (quasi-static crosslink limit, Anna's regime)
- `tau_cross = 1` — baseline, matches earlier 1D runs
- `tau_cross > 1` — crosslinks slower than polymer

**Verified scaling** (from a `TimescaleFirstIntroduction.m` sweep at `dt = 1e-3`, 200k steps):

| `tau_cross` | Mean bond lifetime (s) | Expected `1/koff_eff` (s) | Mean bond fraction |
|---|---|---|---|
| 0.10 | 0.203 | 0.200 | 0.063 |
| 0.03 | 0.063 | 0.060 | 0.071 |
| 0.01 | 0.020 | 0.020 | 0.070 |

Lifetimes track the expected scaling to three significant figures; equilibrium bound fraction is invariant. Note: the safe lower bound on `tau_cross` is set by `max(kon_eff, koff_eff) · dt < 0.1` (otherwise the Poisson rate approximation breaks down); a runtime warning fires when this is violated.

### Numerical scheme

Euler–Maruyama integration with step `dt = 1e-3`, total steps `200000`. Drag coefficient `ζ = 2.5e-3`, thermal energy `kB·T = 4.1 pN·nm`.

---

## Repository Structure

### Entry points
- `TwoD_First_Pass.m` — **baseline 2D simulation** with Hult hard-cutoff crosslink kernel.
- `TimescaleFirstIntroduction.m` — **2D simulation with `tau_cross`** timescale-separation parameter. Use this to sweep across crosslink-vs-polymer dynamics regimes. Default `tau_cross = 1.0` reproduces `TwoD_First_Pass.m` exactly.

### Build-up scripts (1D → 2D, simple → full model)
- `One_Dimension_Tethered_Beads.m` — minimal 1D tethered chain, WLC only
- `OneD_Bead_Dynamics.m` — 1D dynamics study 
- `Tethered_OneD_With_Statistics.m` — adds cluster statistics collection
- `OneD_Tethered_with_Crosslinking.m` — adds dynamic bonds in 1D
- `Tethered_Crosslinking_OneBondPerBead_BenVCaitKon.m` — kernel comparison in 1D
- `Tethered_crosslinking_OneBondPerBead_SwitchableKon.m` — runtime-switchable kernel
- `Tethered_Crosslinking_and_Network_Analysis.m` — adds bond-network analysis
- `Two_Dimensional_Tethered_Beads.m` — 2D extension of the base model
- `Comprehensive_beads_first_pass.m` — Comprehensive 1D
- `Tethered_1D_FixedBeads_MoveWalls.m` — wall-motion variant for parameter sweeps
- `TetheredWithHultKernel.m` — Hult et al. (2017) kernel for comparison

### Analysis & sensitivity
- `AlphaSensitivityAnalysis.m` — sensitivity to WLC interpolation factor α
- `cSensivityAnalysis.m` — sensitivity to EV strength `cEV`
- `Isolating_Parameter_Regime_WLC_and_EV.m` — maps the WLC/EV phase boundary
- `BondFractionComparison.m` — fraction of bonded beads under different conditions
- `ComparingKonModels.m`, `ComparingKon_With_Proper_BensKon.m`, `MatchingKons.m` — kernel comparison & matching
- `ComparisonWithBondTrackingandLifetime.m`, `ComparisonWithLifetimeTracking.m` — bond-lifetime analysis
- `Exploring_Markov_Chains.m` — Markov state analysis of cluster configurations
- `stericcoreincluded.m` — variant with hard steric core
- `workoverbreak.m` — sensitivity analysis

### Helper functions
- `H.m`, `M.m`, `dM.m`, `dS.m`, `gradH.m` — Hamiltonian, mobility, derivatives
- `compute_W.m`, `compute_dW.m` — work functional and its gradient
- `compute_state_count.m`, `state_configurations.m` — combinatorial state enumeration
- `create_network.m`, `create_quasi_network_a2.m` — bond-network initialization
- `deterministic_U.m`, `deterministic_force.m` — deterministic potential / force
- `drift_mat.m` — drift matrix for the SDE system
- `escape_statistics.m`, `escape_stats_single.m` — escape-time statistics
- `find_minima.m` — locates potential minima
- `simulate_system.m`, `run_experiment.m` — generic simulation drivers
- `string_descend.m`, `StringMethodExample.m` — string method (transition path) utilities
- `switching.m` — kernel-switching helper
- `ind2sub4up.m` — index conversion for upper-triangular matrices
- `combine_mc_stats.m` — combines Monte Carlo runs
- `plot2dbeads.m`, `plotbeaddist.m` — visualization

### Subdirectories
- `1d/` — earlier 1D-only versions of `H.m`, `M.m`, switching code, and a well-escape analysis. **Note:** several filenames in `1d/` shadow files in the repository root; keep one or the other on the MATLAB path at a time.
- `notebooks/` — `Figure1.mlx` through `Figure5.mlx`, MATLAB Live Scripts that reproduce presentation figures.

### Figures (committed)
`InitialFlucturations.png`, `Initialtethered.png`, `Position1.png`, `fig2.png`–`fig6.png`, `first.png`, `forcebalance.png`, `stabledynamics.png` — canonical figures referenced in slides and write-ups.

---

## Workflow Notes

- Newly generated `.mat`, `.fig`, `.avi`, and `.mp4` files are git-ignored by default. If you want to commit a specific output (e.g., a canonical result `.mat`), use `git add -f <file>`.
- The `notebooks/` Live Scripts are intended to be the reproducible figure source — re-run them rather than copy-pasting from older scripts.
- Tag a version before any breaking change (`git tag v0.X-description`). Tagged checkpoints:
  - `v0.2-tau-cross` — first verified `tau_cross` timescale-separation implementation.

---

## References

**Primary modeling sources (UNC/Bloom–Forest group)**

- Vasquez, P. A. et al. (2016). Entropy gives rise to topologically associating domains. *Nucleic Acids Research*, 44(12), 5540–5549.
- Hult, C. et al. (2017). Enrichment of dynamic chromosomal crosslinks drive phase separation of the nucleolus. *Nucleic Acids Research*, 45(19), 11159–11173.
- Walker, B. et al. (2019). Toward a dynamic model of the nucleolus: a stochastic model of transient gene-gene crosslinking. *PLOS Computational Biology*, 15(8), e1007124.
- He, Y. et al. (2020). Statistical mechanics of chromosomes. *Nucleic Acids Research*, 48(20), 11284–11303.
- Kolbin, D. et al. (2023). Polymer modeling reveals interplay between physical properties of chromosomal DNA and the size and distribution of condensin-based chromatin loops. *Genes*, 14(12), 2193.
- Coletti, A., Newhall, K. A., Walker, B. L., & Bloom, K. (2024). Different relative scalings between transient forces and thermal fluctuations tune regimes of chromatin organization. *arXiv:2401.06921*.
- Walker, B. L. *Emergent Structure and Dynamics from Stochastic Pairwise Crosslinking in Chromosomal Polymer Models* (thesis).

**Biology of SMC complexes**

- Holmes, V. F., & Cozzarelli, N. R. (2000). Closing the ring: links between SMC proteins and chromosome partitioning, condensation, and supercoiling. *PNAS*, 97(4), 1322–1324.
- Eeftens, J. M. et al. (2016). Condensin Smc2-Smc4 dimers are flexible and dynamic. *Cell Reports*, 14(8), 1813–1818.

**Other chromatin polymer models**

- Bohn, M., & Heermann, D. W. (2010). Diffusion-driven looping provides a consistent framework for chromatin organization. *PLOS ONE*, 5(8), e12218.
- Barbieri, M. et al. (2012). Complexity of chromatin folding is captured by the strings and binders switch model. *PNAS*, 109(40), 16173–16178.
- Brackley, C. A. et al. (2013). Nonspecific bridging-induced attraction drives clustering of DNA-binding proteins and genome organization. *PNAS*, 110(38), E3605–E3611.
- Tjong, H. et al. (2012). Physical tethering and volume exclusion determine higher-order genome organization in budding yeast. *Genome Research*, 22, 1295–1305.

---

## License

See `LICENSE.txt`.
