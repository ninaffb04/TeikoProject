PRAGMA foreign_keys = ON;

DROP VIEW IF EXISTS population_frequencies;
DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    age        INTEGER NOT NULL,
    sex        TEXT NOT NULL,
    condition  TEXT NOT NULL,
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL,
    project_id                 TEXT NOT NULL,
    treatment                  TEXT NOT NULL,
    response                   TEXT,
    sample_type                TEXT NOT NULL,
    time_from_treatment_start  INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id),
    FOREIGN KEY (project_id) REFERENCES projects(project_id)
);

CREATE TABLE cell_counts (
    sample_id   TEXT NOT NULL,
    population  TEXT NOT NULL,
    count       INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id)
);

CREATE INDEX idx_subjects_project_id ON subjects(project_id);
CREATE INDEX idx_samples_subject_id ON samples(subject_id);
CREATE INDEX idx_samples_project_id ON samples(project_id);
CREATE INDEX idx_samples_treatment ON samples(treatment);
CREATE INDEX idx_samples_sample_type ON samples(sample_type);
CREATE INDEX idx_samples_response ON samples(response);
CREATE INDEX idx_samples_time_from_treatment_start ON samples(time_from_treatment_start);
CREATE INDEX idx_cell_counts_population ON cell_counts(population);

-- View: per-sample, per-population relative frequency.
-- total_count is the sum of all five population counts for that sample.
CREATE VIEW population_frequencies AS
SELECT
    cc.sample_id AS sample,
    tot.total_count AS total_count,
    cc.population AS population,
    cc.count AS count,
    ROUND(100.0 * cc.count / tot.total_count, 4) AS percentage
FROM cell_counts cc
JOIN (
    SELECT sample_id, SUM(count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
) tot ON tot.sample_id = cc.sample_id;
