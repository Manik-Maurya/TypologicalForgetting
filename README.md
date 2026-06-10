# Beyond the Single Forgetting Curve

**Simulation and figure-generation code for:**

Maurya, M. (2026). Beyond the single forgetting curve: A typologically-extended 
Wickelgren function, a simulation study, and a pre-registered test of 
error-type-specific memory decay. *PsyArXiv preprint.*  
DOI: [https://doi.org/10.31234/osf.io/INSERT](https://doi.org/10.31234/osf.io/INSERT)

**Pre-registration (OSF, frozen):**  
https://doi.org/10.17605/OSF.IO/YU6P5

**All materials (PDFs, priors table, OSF archive):**  
[https://osf.io/qm39j]

---

## Repository contents

| File | Description |
|------|-------------|
| `Maurya_2026_TypologicalForgetting_SimCode.R` | 7,200-fit operating-characteristics simulation (120 reps × 60 grid cells). Generates `outputs/power_grid.json`. |
| `Maurya_2026_TypologicalForgetting_FigureCode.R` | Produces Figure 1 (predicted curve divergence) and Figure 2 (power surface) from `power_grid.json`. |
| `outputs/power_grid.json` | Simulation output grid (deposited after full run). |

---

## Requirements

```r
# R >= 4.2
install.packages(c("brms", "loo", "dplyr", "purrr", "tibble", "tidyr", "ggplot2"))
# cmdstanr (see https://mc-stan.org/cmdstanr/)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Runtime note

The full 500-replication loop takes **12–25 hours** on a 4-core laptop.  
For a fast test run, set `N_REPS <- 5` at line [see script].

---

## Citation

If you use this code, please cite the preprint above.  
ORCID: [0009-0005-3554-693X](https://orcid.org/0009-0005-3554-693X)

## Licence

Code: CC BY 4.0 — Manik Maurya, 2026.
