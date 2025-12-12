# GOTM Local Lakes ISIMIP Simulations

**This repository contains 1D hydrodynamic simulations of local lakes using the General Ocean Turbulence Model (GOTM) as part of the ISIMIP framework.**
GOTM has been adapted for lakes to simulate water temperature and mixing processes under historical and future climate scenarios.([isimip.org][1])

---

## Overview

This project performs **ISIMIP2b local lake simulations** using GOTM.
The simulations include:

* calibration runs to fit model parameters to observed data,
* diagnostic plots to visualize the calibration performance,
* model output (daily vertical profiles of temperature and related variables).

## Repository Structure

```
.
├── run/                         # contains initial plots, calibration plots, and resulting .db files calibration
├── codes/                       # codes used during simulations
├── LakeCharacteristics.csv      # Basic information about the lakes used in this study
├── GOTM_MLE_performance.csv     # Final performance of the calibrated lake model for each lake using Maximum Likelihood Estimation
└── README.md                    # This file
```

---

## 🔧 Requirements

To run and analyse the simulations, you’ll need:

* **GOTM (General Ocean Turbulence Model)** — source code or compiled binary
  Download/build from: [https://gotm.net/](https://gotm.net/)
* **Python or R** environment (for plots & diagnostics)

  * Python: `xarray`, `matplotlib`, `netCDF4`
  * R: `gotmtools`, `ggplot2`, `ncdf4` (optional diagnostic plotting)

---

## Running Simulations

### 1. Prepare Input Data

Each lake must have:

* a *hypsograph* (depth vs area profile),
* meteorological forcing (daily climate data),
* initial temperature profile.

Place these in `run/`.

### 2. Configure GOTM calibration runs

Each run in `run` uses a configuration XML file:

```bash
Rscript codes/calibration.R
```

Adjust the paths and configs accordingly.

---

## Diagnostics

After running the calibration workflow, the folder `run/calibration_plot/` contains daily temperature profiles (surface to bottom) from calibrated simulations and observations:

### Example Plots
📄 **View the plot:**
[https://github.com/danielmerbet/GOTM_local/blob/main/run/calibration_plot/kinneret.pdf](https://github.com/danielmerbet/GOTM_local/blob/main/run/calibration_plot/kinneret.pdf)

### Example Metrics

* `calibration_metrics.csv` — key model performance statistics (e.g., RMSE, bias).

*Replace these placeholders with your actual figures and detailed captions.*

---

### Calibration Metrics

Calibration diagnostics included:

* RMSE (Root Mean Squared Error)
* Bias
* NSE
* r

The resulting metrics can be found in file GOTM_MLE_performance.csv
---

## 📎 References & Further Reading

* **ISIMIP Lake Sector** — a framework for impact modelling of lakes: ([gmd.copernicus.org][2])
* **GOTM Model Description** — turbulence model adapted for lake simulations.([gotm.net][3])

---

## 📄 License

This project is released under the **CC BY 4.0*.

---

[1]: https://www.isimip.org/impactmodels/details/242/?
[2]: https://gmd.copernicus.org/articles/15/4597/2022/gmd-15-4597-2022.html
[3]: https://gotm.net/
