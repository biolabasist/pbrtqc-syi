# pbrtqc-syi

R code for the analyses reported in

> Şara B, Pınar A. *No one rule fits all: Youden index heatmaps for analyte-specific patient-based real-time quality control parameter selection in 23 routine chemistry tests using biological variation.*

The code selects exponentially weighted moving average (EWMA) and moving average (MA) settings for patient-based real-time quality control (PBRTQC) from a laboratory's own chronologically ordered patient results. Truncation limits are derived from biological variation, systematic errors are introduced into the result stream before truncation, and each configuration is scored with the per-result Youden index summed over errors within the allowable total error (Summed Youden Index, SYI).

## Contents

```
config/analytes.csv           biological variation, CLIA 2024 criteria, units
R/pbrtqc.R                    truncation, transformation, charts, error simulation
R/analysis.R                  stream analysis, Youden index, SYI with bootstrap, detection delay
R/setup.R                     paths and allowable total error rules
R/plots.R                     plotting helpers
scripts/01_prepare_data.R     builds de-identified analysis files from LIS exports
scripts/02_iqc.R              imprecision, sigma metrics, results between IQC events
scripts/03_run_scenarios.R    simulations (primary and sensitivity analyses)
scripts/03b_direction.R       detection counted only in the direction of the error; chance signals
scripts/sim_direction.R       error simulation recording upward and downward signals
scripts/03c_recalibration.R   temporal validation with monthly re-estimation of mu0 and sigma0
scripts/03d_exclusion_monitor.R  moving count of results outside the truncation limits
scripts/03e_grid_extension.R  configurations beyond the edges of the grid
scripts/04_summaries.R        summary tables
scripts/04b_direction_summary.R  summary of the direction-specific analysis
scripts/04c_additional_summaries.R  alarm-limit, re-estimation and grid-extension summaries
scripts/05_real_world.R       alarm burden and agreement with the IQC record
scripts/06_workflow_figure.R  Figure 1
scripts/07_figures.R          Figures 3-5 and supplementary figures
scripts/08_control_charts.R   Figure 2 and data for the interactive charts
scripts/09_interactive_html.py  interactive full-year charts (self-contained HTML)
```

## Requirements

R (≥ 4.4) with `openxlsx`, `ggplot2`, `patchwork`, `scales` and `viridisLite`. The interactive charts are written with Python (≥ 3.10) and `plotly` (≥ 6) and `pandas`; they show the control statistic only, not individual patient results.

## Data

Patient-level data are not included. The scripts expect two files in the data directory:

`stream_2023.rds`, one row per patient result:

| column | content |
|---|---|
| `analyte` | analyte code as in `config/analytes.csv` |
| `time` | result time (POSIXct) |
| `value` | result in the units of `config/analytes.csv` |
| `analyzer` | analyzer identifier (integer) |

`iqc_2023.rds`, one row per internal quality control result, with columns `iqc_name`, `analyzer`, `module`, `level`, `time`, `value`, `westgard` and `lot`.

`scripts/01_prepare_data.R` shows how these files were built from the laboratory information system exports; other laboratories only need to produce the two files above.

## Running

From the repository root:

```sh
export PBRTQC_DATA=/path/to/data      # default ../data
export PBRTQC_OUT=/path/to/output     # default ../output

Rscript scripts/02_iqc.R
for s in primary stream_a1 stream_a2 derivation validation cva_measured median logsym \
         tea_greater tea_eflm permuted; do
  Rscript scripts/03_run_scenarios.R $s
done
Rscript scripts/03_run_scenarios.R analyzer1   # uses the onsets of the primary analysis
Rscript scripts/03b_direction.R
Rscript scripts/03c_recalibration.R
Rscript scripts/04_summaries.R
Rscript scripts/04b_direction_summary.R
Rscript scripts/03d_exclusion_monitor.R          # needs recommended.csv from 04_summaries.R
Rscript scripts/03e_grid_extension.R
Rscript scripts/04c_additional_summaries.R
Rscript scripts/05_real_world.R
Rscript scripts/06_workflow_figure.R
Rscript scripts/07_figures.R
Rscript scripts/08_control_charts.R
python scripts/09_interactive_html.py "$PBRTQC_OUT"
```

Each scenario takes roughly 20–80 minutes for 23 analytes on a desktop computer; scenarios are independent and can be run in parallel.

## Method in brief

- Truncation: RCVg = √2 × 1.96 × √(CVg² + CVi² + CVa²), CVa = 0.5·CVi, applied around the annual mean; natural logarithm of the retained results (ln(x + 1) when values below 1 occur).
- Charts: EWMA with λ ∈ {0.02, 0.05, 0.10, 0.20} and L ∈ {2.5, 2.7, 3.0, 3.5, 4.0, 4.5, 5.0}; MA with w ∈ {10, 25, 50, 75, 100} and k ∈ {3.0, 3.5, 4.0, 4.5, 5.0}. Results outside the truncation limits do not update the statistic.
- In-control performance: each chart is run over the whole year; specificity is one minus the proportion of results at which the chart is out of control; alarm episodes per day are counted.
- Error simulation: 100 onset positions spread over the year, the same for all configurations; proportional errors of ±0.25 to ±3 × TEa applied to raw results before truncation; the chart continues from its in-control state; follow-up of B results.
- Scoring: Se = (B − d + 1)/B for a first signal (either direction) at result d; J = Se + Sp − 1; SYI = sum of J over |bias| ≤ TEa, with bootstrap confidence intervals over onsets.
- Selection: highest SYI among configurations with at most one alarm episode per day.
- Large errors: a moving count of results outside the truncation limits among the last 50 results, with a limit allowing at most 0.1 alarm episodes per day, complements the chart.

## License

MIT, see `LICENSE`.
