"""Run all analytical steps against clinical_trial.db and write outputs/.

Usage:
    python analysis.py

Requires clinical_trial.db to already exist (run load_data.py first).
Produces:
    outputs/summary_table.csv        -- Part 2: per-sample population frequencies
    outputs/statistical_results.csv  -- Part 3: Mann-Whitney U + BH correction
    outputs/plots/responder_boxplots.html -- Part 3: boxplots per population
    outputs/baseline_samples.csv     -- Part 4: baseline sample/subject summaries
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = REPO_ROOT / "data" / "clinical_trial.db"
OUTPUTS_DIR = REPO_ROOT / "outputs"
PLOTS_DIR = OUTPUTS_DIR / "plots"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
ALPHA = 0.05


def get_connection() -> sqlite3.Connection:
    if not DB_PATH.exists():
        sys.exit(f"ERROR: {DB_PATH} not found. Run load_data.py first.")
    return sqlite3.connect(DB_PATH)


# ---------------------------------------------------------------------------
# Part 2: population frequencies (reads the population_frequencies SQL view)
# ---------------------------------------------------------------------------

def part2_summary_table(conn: sqlite3.Connection) -> pd.DataFrame:
    df = pd.read_sql_query(
        "SELECT sample, total_count, population, count, percentage "
        "FROM population_frequencies ORDER BY sample, population",
        conn,
    )
    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUTS_DIR / "summary_table.csv", index=False)
    print(f"  wrote outputs/summary_table.csv ({len(df)} rows)")
    return df


# ---------------------------------------------------------------------------
# Part 3: melanoma + miraclib + PBMC responder vs non-responder comparison
# ---------------------------------------------------------------------------

def part3_analysis_dataset(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
        SELECT s.sample_id AS sample, s.response, pf.population, pf.percentage
        FROM samples s
        JOIN population_frequencies pf ON pf.sample = s.sample_id
        JOIN subjects sub ON sub.subject_id = s.subject_id
        WHERE sub.condition = 'melanoma'
          AND s.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
    """
    return pd.read_sql_query(query, conn)


def part3_statistics(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pop in POPULATIONS:
        pop_df = df[df["population"] == pop]
        responders = pop_df.loc[pop_df["response"] == "yes", "percentage"]
        non_responders = pop_df.loc[pop_df["response"] == "no", "percentage"]

        u_stat, p_value = mannwhitneyu(responders, non_responders, alternative="two-sided")

        rows.append({
            "population": pop,
            "responder_n": len(responders),
            "non_responder_n": len(non_responders),
            "responder_median": responders.median(),
            "non_responder_median": non_responders.median(),
            "u_statistic": u_stat,
            "p_value": p_value,
        })

    result = pd.DataFrame(rows)
    _, adj_p, _, _ = multipletests(result["p_value"], alpha=ALPHA, method="fdr_bh")
    result["adjusted_p_value"] = adj_p
    result["significant"] = result["adjusted_p_value"] < ALPHA

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUTS_DIR / "statistical_results.csv", index=False)
    print(f"  wrote outputs/statistical_results.csv ({len(result)} rows)")
    return result


def part3_boxplots(df: pd.DataFrame, stats: pd.DataFrame) -> None:
    labels = {
        "b_cell": "B Cell", "cd8_t_cell": "CD8 T Cell", "cd4_t_cell": "CD4 T Cell",
        "nk_cell": "NK Cell", "monocyte": "Monocyte",
    }
    fig = make_subplots(rows=1, cols=len(POPULATIONS), subplot_titles=[labels[p] for p in POPULATIONS])

    for i, pop in enumerate(POPULATIONS, start=1):
        pop_df = df[df["population"] == pop]
        stat_row = stats.loc[stats["population"] == pop].iloc[0]
        for response, color in [("yes", "#2E86AB"), ("no", "#E07A5F")]:
            group = pop_df[pop_df["response"] == response]
            fig.add_trace(
                go.Box(
                    y=group["percentage"],
                    name="Responder" if response == "yes" else "Non-responder",
                    marker_color=color,
                    boxpoints="all",
                    jitter=0.4,
                    pointpos=0,
                    text=group["sample"],
                    hovertemplate="Sample: %{text}<br>Frequency: %{y:.2f}%<extra></extra>",
                ),
                row=1, col=i,
            )
        fig.update_yaxes(title_text="Relative frequency (%)" if i == 1 else None, row=1, col=i)
        sig = "significant" if stat_row["significant"] else "not significant"
        fig.layout.annotations[i - 1].update(
            text=f"{labels[pop]}<br><span style='font-size:11px'>adj. p={stat_row['adjusted_p_value']:.3g} ({sig})</span>"
        )

    fig.update_layout(
        title=dict(
            text="Immune Cell Population Frequencies: Responders vs Non-Responders<br>"
                 "<sub>Melanoma / Miraclib / PBMC samples</sub>",
        ),
        height=600,
        width=1500,
        boxmode="group",
        showlegend=False,
        margin=dict(t=130, b=60),
    )

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOTS_DIR / "responder_boxplots.html"
    fig.write_html(out_path)
    print(f"  wrote {out_path.relative_to(REPO_ROOT)}")


# ---------------------------------------------------------------------------
# Part 4: baseline melanoma PBMC miraclib samples
# ---------------------------------------------------------------------------

BASELINE_FILTER = """
    FROM samples s
    JOIN subjects sub ON sub.subject_id = s.subject_id
    WHERE sub.condition = 'melanoma'
      AND s.sample_type = 'PBMC'
      AND s.treatment = 'miraclib'
      AND s.time_from_treatment_start = 0
"""


def part4_baseline(conn: sqlite3.Connection) -> dict:
    by_project = pd.read_sql_query(
        f"SELECT s.project_id AS project, COUNT(*) AS sample_count {BASELINE_FILTER} GROUP BY s.project_id ORDER BY s.project_id",
        conn,
    )
    by_response = pd.read_sql_query(
        f"SELECT s.response AS response, COUNT(DISTINCT s.subject_id) AS subject_count {BASELINE_FILTER} GROUP BY s.response ORDER BY s.response",
        conn,
    )
    by_sex = pd.read_sql_query(
        f"SELECT sub.sex AS sex, COUNT(DISTINCT s.subject_id) AS subject_count {BASELINE_FILTER} GROUP BY sub.sex ORDER BY sub.sex",
        conn,
    )

    by_project.insert(0, "breakdown", "samples_by_project")
    by_response.insert(0, "breakdown", "subjects_by_response")
    by_sex.insert(0, "breakdown", "subjects_by_sex")

    by_project = by_project.rename(columns={"project": "group_value", "sample_count": "value"})
    by_response = by_response.rename(columns={"response": "group_value", "subject_count": "value"})
    by_sex = by_sex.rename(columns={"sex": "group_value", "subject_count": "value"})

    combined = pd.concat([by_project, by_response, by_sex], ignore_index=True)

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUTPUTS_DIR / "baseline_samples.csv", index=False)
    print(f"  wrote outputs/baseline_samples.csv ({len(combined)} rows)")

    return {"by_project": by_project, "by_response": by_response, "by_sex": by_sex}


def main() -> None:
    conn = get_connection()
    try:
        print("Part 2: population frequencies ...")
        part2_summary_table(conn)

        print("Part 3: responder vs non-responder analysis ...")
        analysis_df = part3_analysis_dataset(conn)
        stats = part3_statistics(analysis_df)
        part3_boxplots(analysis_df, stats)

        print("Part 4: baseline melanoma/PBMC/miraclib summaries ...")
        part4_baseline(conn)
    finally:
        conn.close()

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
