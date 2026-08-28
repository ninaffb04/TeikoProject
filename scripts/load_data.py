"""Load cell-count.csv into a SQLite relational database.

Usage:
    python load_data.py

Recreates clinical_trial.db from scratch every run, so the pipeline is
always reproducible from the CSV alone. Paths are resolved relative to
this script's location, not the caller's working directory.
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
CSV_PATH = REPO_ROOT / "cell-count.csv"
SCHEMA_PATH = REPO_ROOT / "schema.sql"
DB_PATH = REPO_ROOT / "clinical_trial.db"

REQUIRED_COLUMNS = [
    "project", "subject", "condition", "age", "sex", "treatment",
    "response", "sample", "sample_type", "time_from_treatment_start",
    "b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte",
]
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
VALID_RESPONSES = {"yes", "no"}


def load_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        sys.exit(f"ERROR: could not find input CSV at {path}")
    return pd.read_csv(path)


def validate(df: pd.DataFrame) -> None:
    """Fail fast on structural or data-quality problems before writing anything."""
    errors = []

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        sys.exit(f"ERROR: missing required columns: {missing_cols}")

    if df["sample"].isnull().any():
        errors.append("some rows have a missing sample id")
    if df["sample"].duplicated().any():
        dupes = df.loc[df["sample"].duplicated(), "sample"].tolist()
        errors.append(f"duplicate sample ids found: {dupes}")

    for col in POPULATIONS:
        if not pd.api.types.is_numeric_dtype(df[col]):
            errors.append(f"population column '{col}' is not numeric")
        elif df[col].isnull().any():
            errors.append(f"population column '{col}' has missing values")
        elif (df[col] < 0).any():
            errors.append(f"population column '{col}' has negative counts")

    bad_treatment = set(df["treatment"].dropna().unique()) - {"miraclib", "phauximab", "none"}
    if bad_treatment:
        errors.append(f"unexpected treatment values: {bad_treatment}")

    bad_response = set(df["response"].dropna().unique()) - VALID_RESPONSES
    if bad_response:
        errors.append(f"unexpected response values: {bad_response}")

    # response should be null exactly when there is no treatment, and set otherwise
    treated = df[df["treatment"] != "none"]
    if treated["response"].isnull().any():
        errors.append("some treated samples (treatment != 'none') are missing a response")

    if errors:
        sys.exit("ERROR: validation failed:\n  - " + "\n  - ".join(errors))


def build_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA_PATH.read_text())


def insert_data(conn: sqlite3.Connection, df: pd.DataFrame) -> dict:
    cur = conn.cursor()

    projects = sorted(df["project"].unique())
    cur.executemany(
        "INSERT INTO projects (project_id) VALUES (?)",
        [(p,) for p in projects],
    )

    subjects = df.drop_duplicates(subset=["subject"])[
        ["subject", "project", "age", "sex", "condition"]
    ]
    cur.executemany(
        "INSERT INTO subjects (subject_id, project_id, age, sex, condition) "
        "VALUES (?, ?, ?, ?, ?)",
        list(subjects.itertuples(index=False, name=None)),
    )

    samples = df[[
        "sample", "subject", "project", "treatment", "response",
        "sample_type", "time_from_treatment_start",
    ]]
    cur.executemany(
        "INSERT INTO samples (sample_id, subject_id, project_id, treatment, "
        "response, sample_type, time_from_treatment_start) VALUES (?, ?, ?, ?, ?, ?, ?)",
        list(samples.itertuples(index=False, name=None)),
    )

    long_df = df.melt(
        id_vars=["sample"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    cur.executemany(
        "INSERT INTO cell_counts (sample_id, population, count) VALUES (?, ?, ?)",
        list(long_df.itertuples(index=False, name=None)),
    )

    conn.commit()

    return {
        "projects": len(projects),
        "subjects": len(subjects),
        "samples": len(samples),
        "cell_counts": len(long_df),
    }


def main() -> None:
    print(f"Loading {CSV_PATH} ...")
    df = load_csv(CSV_PATH)
    print(f"  {len(df)} rows read")

    print("Validating data ...")
    validate(df)
    print("  validation passed")

    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        conn.execute("PRAGMA foreign_keys = ON")
        print(f"Creating schema at {DB_PATH} ...")
        build_schema(conn)

        print("Inserting data ...")
        counts = insert_data(conn, df)
    finally:
        conn.close()

    print("\nLoad summary:")
    for table, n in counts.items():
        print(f"  {table:12s} {n:>6d} rows")
    print(f"\nDatabase written to {DB_PATH}")


if __name__ == "__main__":
    main()
