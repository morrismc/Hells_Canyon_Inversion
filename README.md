# Hells Canyon Bayesian River Profile Inversion

Bayesian inversion of river profiles and cave burial ages to constrain the timing of drainage capture and canyon formation at Hells Canyon, North America's deepest river gorge.

Adapts the stream power inversion frameworks of [Gallen & Fernandez-Blanco (2021)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JF005651) and [Goren et al. (2014)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014JF003079) to incorporate cave sediment burial age constraints from [Morriss et al. (2025)](https://www.pnas.org/doi/10.1073/pnas.2413069122).

## Background

Hells Canyon was carved by rapid incision following a drainage capture event that established the modern Snake River route into the Columbia River system. Cave burial dating constrains a two-phase incision history:

| Phase | Period | Rate | Source |
|-------|--------|------|--------|
| Pre-capture | ~5.5 to ~2.1 Ma | ~0.01 mm/yr | Upper cave deposits |
| Post-capture | ~2.1 Ma to present | 0.09-0.16 mm/yr | Lower cave deposits |

The capture timing is constrained to **2.1 +/- 1.0 Ma** (Morriss et al. 2025).

This code uses these cave constraints as informative priors and likelihood data within a Bayesian MCMC framework to jointly estimate erodibility (K), slope exponent (n), concavity (m/n), and capture timing from tributary river profiles.

## Requirements

- **MATLAB** R2016b or later (for local functions in scripts)
- **[TopoToolbox](https://topotoolbox.wordpress.com/)** v2.4+
- **[Topographic Analysis Kit (TAK)](https://github.com/amforte/Topographic-Analysis-Kit)** (for basin processing)

## Input Data

### 1. Basin DEM (required for inversions)

Basin files in the `Basin_xx_Data.mat` format produced by TAK's `ProcessRiverBasins.m`. These contain:

| Variable | Type | Description |
|----------|------|-------------|
| `DEMoc` | GRIDobj | Original DEM cropped to basin |
| `DEMcc` | GRIDobj | Hydrologically conditioned DEM (preferred for stream analysis) |
| `FDc` | FLOWobj | Flow direction |
| `Ac` | GRIDobj | Drainage area in m^2 (flowacc * cellsize^2) |
| `Sc` | STREAMobj | Stream network |
| `Chic` | chiplot struct | Chi analysis from TopoToolbox |
| `theta_ref` | scalar | Reference concavity used by TAK |
| `RiverMouth` | [x, y, id] | Outlet coordinates and basin number |

To generate these files from a DEM:
```matlab
% In TAK:
ProcessRiverBasins(DEM, FD, river_mouths, 'theta_ref', 0.45);
```

### 2. ksn Table (required for K/n analysis)

An Excel file (`ksnTable.xlsx`) with at minimum these columns:

| Column | Description |
|--------|-------------|
| `Basin` | Basin identifier number |
| `Segment` | 1 = below knickpoint (adjusted), 2 = above knickpoint (relict) |
| `ksn` | Normalized channel steepness (m^0.9 for m/n=0.45) |

### 3. Cave Data (configured in scripts)

Cave burial ages and heights above the modern river, entered directly in the MATLAB scripts as arrays. Format:

```matlab
cave_data = [
    age_yr,  height_m,  age_err_yr,  height_err_m;
    5.5e6,   100,       0.5e6,       20;
    3.5e6,   120,       0.5e6,       20;
    1.5e6,   250,       0.3e6,       20;
];
```

## Repository Structure

```
Hells_Canyon_Inversion/
|-- README.md
|-- LICENSE
|-- HellsCanyon_BayesianInversion/     MCMC inversion code (this project)
|   |-- run_hc_analysis.m              Quick-start: runs all 3 analyses
|   |-- ksn_erosion_analysis.m         Analysis 1: estimate K and n from ksn
|   |-- hc_linear_inversion_with_caves.m  Analysis 2: linear inversion (n=1)
|   |-- main_hc_bayes_inversion.m      Analysis 3: full Bayesian MCMC
|   |-- hc_river_forward_model.m       Forward model (implicit FD, n!=1)
|   |-- cave_forward_model.m           Cave height predictions
|   |-- hc_loglikelihood.m             Combined stream + cave likelihood
|   |-- logprior_hc.m                  Informative priors from caves
|   |-- logproposal_hc.m              Gaussian random walk proposal
|   |-- prepare_hc_stream_data.m       Extract stream data from TAK Basin files
|   |-- plot_hc_results.m              MCMC diagnostics and result figures
|
|-- Block_Uplift_Linear_Inversion_Models/  Gallen (2018, 2020) codes
    |-- dimensional_block_uplift/      Linear inversion with known K
    |-- nondimensional_block_uplift/   Linear inversion in chi-space
    |-- block_uplift_ksn/              ksn inversion
```

## Workflow

### Quick Start

Edit paths in `run_hc_analysis.m` and run:

```matlab
% 1. Set your file paths
ksn_file = 'path/to/ksnTable.xlsx';
dem_file = 'path/to/Basin_80_Data.mat';

% 2. Update cave_data with your actual cave burial ages/heights

% 3. Run
run_hc_analysis
```

### Step-by-Step

#### Step 1: Estimate K and n from ksn data

```matlab
results = ksn_erosion_analysis('ksnTable.xlsx', ...
    'E_relict',     0.01, ...     % Pre-capture rate (mm/yr) from caves
    'E_adjusted',   0.125, ...    % Post-capture rate (mm/yr) from caves
    'n_samples',    100000);
```

This pairs ksn observations above and below knickpoints with cave-derived erosion rates via `E = K * ksn^n` to jointly estimate K and n. A key diagnostic is whether n = 1 is consistent with the data: if the ksn ratio (~2x) is much smaller than the erosion rate ratio (~9-16x), it implies either n >> 1 or spatially variable K.

**Outputs:** 6-panel figure (`ksn_erosion_analysis.png`), results struct with K, n estimates and credible intervals.

#### Step 2: Linear inversion with cave overlay (n=1)

```matlab
% Load Basin DEM from TAK
data = load('Basin_80_Data.mat');
DEM = data.DEMoc;

% Test multiple K values
K_range = [7e-8, 5.73e-7, 7e-7];

results = hc_linear_inversion_with_caves(DEM, K_range, ...
    'tau_inc',    1.5e6, ...
    'mn',         0.45, ...
    'flowOption', 'fill', ...
    't_capture',  2.1e6, ...
    't_capture_err', 1.0e6);
```

Wraps Gallen's `linear_inversion_block_uplift_erodibility.m` to run the n=1 inversion at multiple K values and overlay the cave-constrained capture window (shaded band) on the recovered uplift history.

**Outputs:** Sensitivity figure comparing K values, summary overlay figure, struct with peak timing for each K.

#### Step 3: Prepare stream data for MCMC

```matlab
sd = prepare_hc_stream_data('Basin_80_Data.mat', ...
    'mn',          0.45, ...
    'output_file', 'hc_stream_data.mat');
```

Reads a TAK `Basin_xx_Data.mat` file (or GeoTIFF), extracts stream network, computes chi coordinates, and saves a standardized struct for the MCMC. Automatically uses `DEMcc` (conditioned DEM) if available, and picks up `theta_ref` from the TAK file.

#### Step 4: Run full Bayesian MCMC (free n)

Edit `main_hc_bayes_inversion.m`:

```matlab
% Set paths
stream_data_file = 'hc_stream_data.mat';

% Enter your actual cave data
cave_data = [
    5.5e6, 100, 0.5e6, 20;
    3.5e6, 120, 0.5e6, 20;
    1.5e6, 250, 0.3e6, 20;
];

% For testing: fast settings
n_burnin   = 1e4;
n_postburn = 1e5;

% For publication: increase iterations
% n_burnin   = 3e5;
% n_postburn = 3e6;
```

Then run:

```matlab
main_hc_bayes_inversion
```

**Parameters estimated (6 total):**

| # | Parameter | Units | Prior Bounds | Cave Prior |
|---|-----------|-------|-------------|------------|
| 1 | U_pre | m/yr | [1e-6, 5e-4] | N(0.01e-3, 0.005e-3) |
| 2 | U_post | m/yr | [1e-5, 5e-3] | N(0.125e-3, 0.035e-3) |
| 3 | log10(K) | - | [-9, -4] | Uniform |
| 4 | n | - | [0.5, 10] | Uniform |
| 5 | m/n | - | [0.3, 0.7] | Uniform |
| 6 | t_capture | yr | [0.5e6, 5e6] | N(2.1e6, 1.0e6) |

Physical constraint: U_post > U_pre is enforced (capture must increase incision).

**Outputs:**
- `params_HC_capture.mat` -- Full MCMC chain and posterior samples
- `mMAP_HC_capture.mat` -- MAP model and predictions
- `diagnostics_HC_capture.png` -- Trace plots, acceptance rate, posterior histograms
- `posteriors_HC_capture.png` -- K-n trade-off, timing posterior vs cave prior
- `model_fit_HC_capture.png` -- River profile fit and cave data fit

## Model Description

### Forward Model

Solves the stream power incision equation:

```
dz/dt = U(t) - K * A^m * |dz/dx|^n
```

with a two-phase uplift/incision history:
- Phase 1 (before capture): uniform rate U_pre
- Phase 2 (after capture): uniform rate U_post

Uses the implicit finite-difference scheme of Braun & Willett (2013) with Newton-Raphson iteration for n != 1. Starts from steady state with U_pre and evolves forward through the capture event.

### Cave Forward Model

Predicts the height of a cave above the modern river as cumulative incision since the cave was at river level:

```
If cave_age > t_capture:
    height = U_post * t_capture + U_pre * (cave_age - t_capture)
If cave_age <= t_capture:
    height = U_post * cave_age
```

### Likelihood

Combined Gaussian log-likelihood for stream elevations and cave heights, with a weighting factor `W_s = n_cave / n_stream` applied to the stream component to prevent the more numerous stream observations from dominating.

### Priors

Cave constraints enter two ways:
1. **As informative Gaussian priors** on t_capture, U_pre, U_post (from `logprior_hc.m`)
2. **As data in the likelihood** via cave height predictions (from `cave_forward_model.m`)

Set `use_informative_priors = false` in `main_hc_bayes_inversion.m` to use only the likelihood pathway.

## Key Numerical Values

From previous analysis (Morriss et al. 2025 and ksn/cave rate analysis):

| Parameter | Value | Source |
|-----------|-------|--------|
| K (n=1, adjusted) | 5.73e-7 m^0.1/yr | ksn + cave E_post |
| K (n=1, relict) | 9.2e-8 m^0.1/yr | ksn + cave E_pre |
| ksn (adjusted) | ~159 m^0.9 | Below-KP mean |
| ksn (relict) | ~80 m^0.9 | Above-KP mean |
| Implied n (uniform K) | 3-4 | log(E ratio) / log(ksn ratio) |
| t_capture | 2.1 +/- 1.0 Ma | Cave burial ages |
| E_pre | ~0.01 mm/yr | Upper caves |
| E_post | 0.09-0.16 mm/yr | Lower caves |

## References

- Braun, J., & Willett, S.D. (2013). A very efficient O(n), implicit and parallel method to solve the stream power equation governing fluvial incision and landscape evolution. *Geomorphology*, 180-181, 170-179.
- Forte, A.M., & Whipple, K.X. (2019). Short communication: The Topographic Analysis Kit (TAK) for TopoToolbox. *Earth Surface Dynamics*, 7, 87-95. [GitHub](https://github.com/amforte/Topographic-Analysis-Kit)
- Gallen, S.F. (2018). Lithologic controls on landscape dynamics and aquatic species evolution in post-orogenic mountains. *Earth and Planetary Science Letters*, 493, 150-160.
- Gallen, S.F., & Fernandez-Blanco, D. (2021). A new data-driven Bayesian inversion of fluvial topography clarifies the tectonic history of the Corinth Rift and reveals a channel steepness threshold. *JGR Earth Surface*, 126, e2020JF005651. [GitHub](https://github.com/sfgallen/bayes_profiler)
- Goren, L., Fox, M., & Willett, S.D. (2014). Tectonics from fluvial topography using formal linear inversion: Theory and applications to the Inyo Mountains, California. *JGR Earth Surface*, 119, 1651-1681.
- Morriss, M.C., Mitchell, N.A., Yanites, B.J., Staisch, L.M., & Korup, O. (2025). Cave records reveal recent origin of North America's deepest canyon. *PNAS*, 122(21), e2413069122.
- Mudd, S.M., Attal, M., Milodowski, D.T., Grieve, S.W.D., & Valters, D.A. (2014). A statistical framework to quantify spatial variation in channel gradients using the integral method of channel profile analysis. *JGR Earth Surface*, 119, 138-152.
- Schwanghart, W., & Scherler, D. (2014). Short Communication: TopoToolbox 2 -- MATLAB-based software for topographic analysis and modeling in Earth surface sciences. *Earth Surface Dynamics*, 2, 1-7.

## License

GPL-3.0 (see [LICENSE](LICENSE))
