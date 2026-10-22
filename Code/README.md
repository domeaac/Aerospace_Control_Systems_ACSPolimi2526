# Code

MATLAB scripts for the ANT-X robust attitude control project. Run everything from this folder after calling `PathsInit` (which runs `addpath(genpath(pwd))`).

## Folder layout

```
Code/
├── ANT_X_data.m              # Model data: g, stability/control derivatives and their 1σ uncertainty (%)
├── PathsInit.m               # Adds Code/ and its subfolders to the MATLAB path
│
├── task1_lat_30.m            # Task 1 – nominal lateral (roll) model + hinfstruct design
├── task1_lon_30.m            # Task 1 – nominal longitudinal (pitch) model + hinfstruct design
├── task2_lat.m               # Task 2 – uncertain lateral model (ureal, ±3σ), usample
├── task2_lon.m               # Task 2 – uncertain longitudinal model
├── task3_latGIO.m            # Task 3 – lateral P–Δ form, N/M extraction, singular values
├── task3_lonGIO.m            # Task 3 – longitudinal P–Δ form, N/M extraction, singular values
├── task4GIO.m                # Task 4 – coupled lat+lon robustness: σ̄(M), D-scalings, μ (mussv)
├── task5.m                   # Task 5 – Gaussian uncertainty, Monte Carlo, CDFs, confidence intervals
├── uncertaintyRegions.m      # Uncertainty regions vs circular overbounding
│
├── refWorkspace_*.mat        # Saved workspaces so later tasks can run without re-running earlier ones
│
├── sweeping/                 # Task 1 sweep over the W_p weight parameters
│   ├── task1_lat_sweep.m     #   sweep driver (lateral)
│   ├── task1_lat_func.m      #   one design + metrics for one (M, ω_b, α) point
│   ├── task1_lon_sweep.m     #   sweep driver (longitudinal)
│   └── task1_lon_func.m
│
├── data/                     # Sweep results and their plots
│   ├── task1_{lat,lon}_sweep_results_{30,50}.mat
│   ├── plots_lat.m
│   └── plots_lon.m
│
├── src/
│   ├── plots/                # Plot chunks called from task1_{lat,lon}_30.m (when plot_nd = 1)
│   └── montecarlo/           # Plots of the Monte Carlo running means / standard deviations
│
└── auxiliary/                # Scratch and support scripts
    ├── xi_sweep.m            #   F and S shape vs damping ratio ξ (ω_n = 10 rad/s)
    ├── weights_S_comparison.m#   manual vs makeweight W_p compared with S (needs Task 1 workspace)
    └── test.m                #   quick robstab check (needs Task 3 workspace)
```

## Run order and dependencies

The tasks build on each other's workspace variables:

```
task1_lat_30 ─► task2_lat ─► task3_latGIO ─┐
                                           ├─► task4GIO ─► task5
task1_lon_30 ─► task2_lon ─► task3_lonGIO ─┘
```

To save time, each later script loads a saved workspace instead of re-running the earlier ones:

| Script | Loads | Produced by |
|--------|-------|-------------|
| `task3_latGIO.m` | `refWorkspace_t1t2_lat.mat` | `task1_lat_30` + `task2_lat` |
| `task3_lonGIO.m` | `refWorkspace_t1t2_lon.mat` | `task1_lon_30` + `task2_lon` |
| `task4GIO.m` | runs `task3_lonGIO` and `task3_latGIO`, then saves `refWorkspace_t1t2t3_latlon.mat` | — |
| `task5.m` | `refWorkspace_t4.mat` | `task4GIO` |

The `refWorkspace_t1t2t3_lat.mat` / `_lon.mat` files are intermediate snapshots of the same kind.

### Re-running from scratch

1. `task1_*_30.m` and `task4GIO.m` / `task5.m` start with `clc; clear; close all`. Before chaining the scripts, comment out that line in Tasks 1–3 (as the note at the top of `task4GIO.m` says), or the workspace is cleared halfway through.
2. In `task3_*GIO.m`, uncomment the `task1_*`, `task2_*` and `save(...)` lines to regenerate `refWorkspace_t1t2_*.mat`.
3. In `task5.m`, uncomment `task4` / `save('refWorkspace_t4.mat')` to regenerate the Task 4 workspace.

## Weight sweep (Task 1)

`sweeping/task1_*_sweep.m` loops over the performance-weight parameters (peak M, bandwidth ω_b, and α), calls `task1_*_func.m` for each point, and records settling time, peak control effort, achieved γ, overshoot and the peak of S. Results are stored in `data/task1_*_sweep_results_*.mat` (`_30` and `_50` refer to two sweep settings) and plotted with `data/plots_lat.m` / `data/plots_lon.m` against the requirements (T_set ≤ 0.47 s, γ ≤ 1, reference second-order model with ω_n = 10 rad/s, ξ = 0.9).

## Model data

All numbers come from `ANT_X_data.m`, which returns a struct `Aircraft` with:

- `Aircraft.Stability`: `Yv, Yp, Lp, Xu, Xq, Mq` and their `*_sigma_pct`
- `Aircraft.Control`: `Yd, Ld, Xd, Md` and their `*_sigma_pct`

Change the model or the uncertainty levels there; all tasks read from it.
