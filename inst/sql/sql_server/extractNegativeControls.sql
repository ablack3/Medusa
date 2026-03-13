-- Extract negative control outcomes for instrument validation
-- Returns long-format table of person x negative control outcome presence
-- Parameters:
--   @cdm_database_schema: Schema containing OMOP CDM tables
--   @cohort_database_schema: Schema containing the cohort table
--   @cohort_table: Name of the cohort table
--   @outcome_cohort_id: Cohort definition ID for the main outcome
--   @negative_control_cohort_ids: Comma-separated list of negative control cohort IDs

-- Use a subquery to collapse to one flag per person per negative control
-- outcome, avoiding duplicate rows when a person has multiple cohort entries
-- for the same negative control definition.
SELECT
  mr.person_id,
  nc_flag.outcome_cohort_id,
  nc_flag.has_outcome
FROM #mr_cohort mr
CROSS JOIN (
  SELECT DISTINCT cohort_definition_id AS outcome_cohort_id
  FROM @cohort_database_schema.@cohort_table
  WHERE cohort_definition_id IN (@negative_control_cohort_ids)
) nc_ids
LEFT JOIN (
  SELECT DISTINCT
    nc.subject_id AS person_id,
    nc.cohort_definition_id AS outcome_cohort_id,
    1 AS has_outcome
  FROM @cohort_database_schema.@cohort_table nc
  INNER JOIN #mr_cohort mr2
    ON nc.subject_id = mr2.person_id
  WHERE nc.cohort_definition_id IN (@negative_control_cohort_ids)
{@use_index_date == 1} ? {
    AND nc.cohort_start_date >= mr2.index_date
}
) nc_flag
  ON mr.person_id = nc_flag.person_id
  AND nc_ids.outcome_cohort_id = nc_flag.outcome_cohort_id
;
