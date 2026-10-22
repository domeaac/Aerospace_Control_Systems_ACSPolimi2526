# Robust Attitude Control of the ANT-X Quadrotor

Robust H∞ control design for the ANT-X quadrotor (lateral and longitudinal dynamics), with uncertainty modelling, μ-analysis and Monte Carlo checks, in MATLAB.

*Politecnico di Milano, Aerospace Control Systems, A.Y. 2025/26.*

## Overview
The project designs and verifies cascaded attitude controllers for the roll (lateral) and pitch (longitudinal) axes of the ANT-X quadrotor, working from a linearised model near hover with uncertain stability and control derivatives.

| Task | Topic | What is done |
|------|-------|--------------|
| 1 | Nominal modelling and design | 3-state lateral/longitudinal models, structured H∞ tuning (`hinfstruct`) of a P outer loop + PID inner loop against a performance weight W_p and a control-effort weight W_Q; a sweep over the weight parameters |
| 2 | Uncertain modelling | Derivatives modelled as `ureal` parameters with ±3σ bounds, plant sampling with `usample` |
| 3 | Robust stability | P–Δ / LFT form of the closed loop (`lftdata`), the M matrix and its singular values |
| 4 | Robustness to coupling | Combined lat + lon analysis: unscaled σ̄(M), constant and frequency-dependent D-scalings, μ bounds (`mussv`) |
| 5 | Statistical robustness | Gaussian uncertainty, Monte Carlo campaign, convergence of the estimates, distributions/CDFs, confidence intervals, Bode envelopes compared with the nominal model |
| — | Uncertainty regions | Checks whether circular overbounding of the uncertainty sets is acceptable; it is too conservative, so the exact P–Δ approach is used |

## Repository contents

```
.
├── ACS_project_2526.pdf      # Project assignment
├── presentationACS26.pptx    # Final presentation
├── Code/                     # MATLAB code and saved workspaces (see Code/README.md)
└── README.md
```

## Requirements

- MATLAB (a recent release; `xscale` used in some plots needs R2023b or later)
- Control System Toolbox
- Robust Control Toolbox

## Getting started

```matlab
cd Code
PathsInit          % adds all subfolders to the path
task1_lat_30       % e.g. nominal lateral design
```

See [Code/README.md](Code/README.md) for the run order, how the scripts depend on each other and the saved workspaces.

## Authors

- Domenico Raffaele Acierno
- Giorgio Vedovi
- Miguel Maria Barros dos Santos
