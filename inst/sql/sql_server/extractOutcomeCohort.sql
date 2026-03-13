-- Extract outcome cohort for Mendelian Randomization analysis
-- Parameters:
--   @cdm_database_schema: Schema containing OMOP CDM tables
--   @cohort_database_schema: Schema containing cohort table
--   @cohort_table: Name of the cohort table
--   @outcome_cohort_id: Cohort definition ID for the outcome
--   @washout_days: Required days of prior observation (default 365)
--   @exclude_prior_outcome: Whether to exclude persons with prior outcome (1=yes, 0=no)

SELECT
  c.subject_id AS person_id,
  c.cohort_start_date AS index_date,
  1 AS outcome,
  YEAR(c.cohort_start_date) - p.year_of_birth AS age_at_index,
  p.gender_concept_id,
  op.observation_period_start_date,
  op.observation_period_end_date
INTO #outcome_persons
FROM @cohort_database_schema.@cohort_table c
INNER JOIN @cdm_database_schema.person p
  ON c.subject_id = p.person_id
INNER JOIN @cdm_database_schema.observation_period op
  ON c.subject_id = op.person_id
  AND c.cohort_start_date >= op.observation_period_start_date
  AND c.cohort_start_date <= op.observation_period_end_date
WHERE c.cohort_definition_id = @outcome_cohort_id
  AND DATEDIFF(DAY, op.observation_period_start_date, c.cohort_start_date) >= @washout_days
{@exclude_prior_outcome == 1} ? {
  AND NOT EXISTS (
    SELECT 1
    FROM @cohort_database_schema.@cohort_table c2
    WHERE c2.subject_id = c.subject_id
      AND c2.cohort_definition_id = @outcome_cohort_id
      AND c2.cohort_start_date < c.cohort_start_date
  )
}
;

-- Also get non-outcome controls from the observation period table.
-- Controls are persons with sufficient observation who are NOT in the outcome cohort.
-- Assign observation_period_end_date as the index date so downstream covariate
-- extraction (FeatureExtraction lookback windows) has a valid anchor date.
-- Age is computed at this same reference date for consistency with cases.
SELECT
  p.person_id,
  op.observation_period_end_date AS index_date,
  0 AS outcome,
  YEAR(op.observation_period_end_date) - p.year_of_birth AS age_at_index,
  p.gender_concept_id,
  op.observation_period_start_date,
  op.observation_period_end_date
INTO #control_persons
FROM @cdm_database_schema.person p
INNER JOIN @cdm_database_schema.observation_period op
  ON p.person_id = op.person_id
WHERE DATEDIFF(DAY, op.observation_period_start_date, op.observation_period_end_date) >= @washout_days
  AND p.person_id NOT IN (
    SELECT person_id FROM #outcome_persons
  )
;

-- Combine cases and controls
SELECT * INTO #mr_cohort FROM (
  SELECT * FROM #outcome_persons
  UNION ALL
  SELECT * FROM #control_persons
) all_persons
;

TRUNCATE TABLE #outcome_persons;
DROP TABLE #outcome_persons;

TRUNCATE TABLE #control_persons;
DROP TABLE #control_persons;
