# Loblaw Bio: Miraclib Clinical Trial Analysis

## 1. Project Overview

Bob Loblaw needs to understand how the drug **miraclib** affects immune cell
populations in melanoma patients enrolled in Loblaw Bio's clinical trials.
This project turns the raw `cell-count.csv` export into a reproducible
analysis pipeline and interactive dashboard that:

- loads the trial data into a normalized SQLite database,
- computes immune-cell population frequencies for every sample,
- statistically compares miraclib responders vs non-responders among
  melanoma PBMC samples,
- summarizes the baseline (day 0) miraclib melanoma PBMC cohort by project,
  response, and sex, and
- presents all of the above in an interactive Streamlit dashboard.

## 2. Dataset

Each row of `cell-count.csv` is one **sample**: a blood draw from one subject
at one point in time. Subjects are nested under projects, and each subject
contributes up to three longitudinal samples (day 0 / day 7 / day 14).

```text
Project
   └── Subject (age, sex, condition)
          ├── Sample @ day 0  (treatment, response, sample_type, cell counts)
          ├── Sample @ day 7
          └── Sample @ day 14
```

Clinical metadata per sample: `project`, `subject`, `condition`, `age`,
`sex`, `treatment`, `response`, `sample`, `sample_type`,
`time_from_treatment_start`.

Five immune-cell populations are measured (raw cell counts) per sample:

- `b_cell`
- `cd8_t_cell`
- `cd4_t_cell`
- `nk_cell`
- `monocyte`

The dataset contains 3 projects, 3,500 subjects, and 10,500 samples (exactly
3 per subject). It was inspected for missing values, duplicate sample/subject
IDs, negative cell counts, zero-total samples, and unexpected
treatment/response values before any code was written — the data was clean;
the only "missing" values are `response = NULL` for untreated
(`treatment = 'none'`) subjects, which is expected.

## 3. Setup

```bash
make setup
```

Installs `pandas`, `numpy`, `scipy`, `statsmodels`, `plotly`, and
`streamlit` from `requirements.txt`.

## 4. Run the Pipeline

```bash
make pipeline
```

Runs `load_data.py` then `analysis.py`. This regenerates, from scratch,
every generated artifact:

- `clinical_trial.db` — the SQLite database
- `outputs/summary_table.csv` — per-sample population frequencies (Part 2)
- `outputs/statistical_results.csv` — Mann-Whitney U + BH-adjusted results (Part 3)
- `outputs/plots/responder_boxplots.html` — responder vs non-responder boxplots (Part 3)
- `outputs/baseline_samples.csv` — baseline cohort summaries (Part 4)

No manual steps or pre-generated files are required; the pipeline is fully
reproducible from `cell-count.csv` alone.

## 5. Run the Dashboard

```bash
make dashboard
```

Starts Streamlit at **http://localhost:8501**. The dashboard reads only from
`clinical_trial.db` and `outputs/`, so `make pipeline` must be run first.

## 6. Database Schema

Defined in [schema.sql](schema.sql), four normalized tables:

- **`projects`** (`project_id` PK) — the set of trial projects.
- **`subjects`** (`subject_id` PK, FK → `projects`) — one row per patient:
  `age`, `sex`, `condition`. These attributes belong to the patient, not the
  sample, so they are stored once per subject instead of being repeated
  (and risking drift) on every one of that subject's samples.
- **`samples`** (`sample_id` PK, FK → `subjects`, `projects`) — one row per
  blood draw: `treatment`, `response`, `sample_type`,
  `time_from_treatment_start`.
- **`cell_counts`** (`sample_id`, `population` composite PK, FK → `samples`)
  — one row per (sample, population) pair, storing the raw `count`. This is
  a **long/tidy format** rather than five fixed columns
  (`b_cell`, `cd8_t_cell`, ...) on `samples`. The advantage: adding a new
  population in the future (e.g. `treg`, `dendritic_cell`, `macrophage`)
  only requires inserting new rows — no `ALTER TABLE`, no schema migration,
  and no code that has to change shape when the assay panel changes.

Indexes exist on `samples(subject_id)`, `samples(project_id)`,
`samples(treatment)`, `samples(sample_type)`, `samples(response)`,
`samples(time_from_treatment_start)`, and `cell_counts(population)`. At
10,500 samples these aren't strictly necessary for query speed, but they
demonstrate how the schema would scale to a much larger trial where the
`samples`/`cell_counts` tables are filtered and joined on these columns
constantly.

A SQL view, `population_frequencies`, computes `total_count` and
`percentage` for every (sample, population) pair on the fly, so the
frequency calculation lives in the database rather than being duplicated in
every downstream script.

## 7. Analysis Methodology

**Population frequency (Part 2).** For each sample,
`total_count = b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte`, and
each population's `percentage = count / total_count * 100`. Computed by the
`population_frequencies` SQL view and exported to
`outputs/summary_table.csv`.

**Cohort filtering (Part 3).** The responder analysis is restricted to
samples where `condition = 'melanoma'`, `treatment = 'miraclib'`, and
`sample_type = 'PBMC'`, then split by `response` (`yes` / `no`). Baseline is
**not** applied here — all matching timepoints are included, per Bob's
request in Part 3.

**Mann-Whitney U test.** For each of the five populations, the relative
frequencies of responders vs non-responders are compared with a two-sided
Mann-Whitney U test (a non-parametric test, appropriate since frequency
distributions aren't assumed to be normal and sample sizes per group differ
slightly).

**Benjamini-Hochberg correction.** Because five populations are tested
simultaneously, raw p-values are adjusted with the Benjamini-Hochberg FDR
procedure (`statsmodels.stats.multitest.multipletests`, `method='fdr_bh'`).
A population is called significant if its adjusted p-value is below
`alpha = 0.05`.

**Baseline filtering (Part 4).** Restricted to
`condition = 'melanoma'`, `sample_type = 'PBMC'`, `treatment = 'miraclib'`,
`time_from_treatment_start = 0`. Samples are counted with `COUNT(*)`
(the question asks about samples), while subjects are counted with
`COUNT(DISTINCT subject_id)` (the question asks about subjects, and each
subject appears once at baseline in this filtered set, but the distinct
count is used deliberately since a subject can otherwise contribute
multiple samples).

## 8. Code Structure

- **`load_data.py`** — locates `cell-count.csv`, validates it (required
  columns present, sample IDs present and unique, cell counts numeric and
  non-negative, no unexpected treatment/response values), recreates
  `clinical_trial.db` from `schema.sql`, and inserts projects, subjects,
  samples, and long-format cell counts inside a single transaction. Runs
  standalone with `python load_data.py` — no arguments, no manual DB setup.
- **`analysis.py`** — reads `clinical_trial.db` and writes every file in
  `outputs/`: the Part 2 summary table, the Part 3 filtered
  responder/non-responder dataset, statistics, and boxplots, and the Part 4
  baseline summaries.
- **`schema.sql`** — table definitions, foreign keys, indexes, and the
  `population_frequencies` view.
- **`dashboard.py`** — Streamlit app with three tabs (Data Overview,
  Treatment Response, Baseline Analysis) that reads exclusively from the
  database and `outputs/` — it performs no analysis of its own.

## 9. Outputs

| File | Description |
|---|---|
| `clinical_trial.db` | SQLite database: `projects`, `subjects`, `samples`, `cell_counts`, plus the `population_frequencies` view |
| `outputs/summary_table.csv` | `sample, total_count, population, count, percentage` — 5 rows per sample |
| `outputs/statistical_results.csv` | `population, responder_n, non_responder_n, responder_median, non_responder_median, u_statistic, p_value, adjusted_p_value, significant` |
| `outputs/baseline_samples.csv` | Long-format: `breakdown, group_value, value` for samples-by-project, subjects-by-response, and subjects-by-sex |
| `outputs/plots/responder_boxplots.html` | Interactive Plotly boxplots, one per population, responders vs non-responders |

## 10. Statistical Limitations

The Part 3 Mann-Whitney comparison is run at the **sample** level, but each
subject contributes up to three longitudinal samples (day 0, 7, 14). Samples
from the same subject are correlated with each other, so treating all
samples as independent observations understates the true variance and can
make results look more significant than they are. Because of this, the
comparison should be read as an **exploratory analysis** — a starting point
for hypothesis generation, not confirmation that a population's frequency
predicts treatment response. A confirmatory analysis would need to account
for repeated measures (e.g. a mixed-effects model with subject as a random
effect) or otherwise ensure each subject contributes only one independent
observation to the test.
