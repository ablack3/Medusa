-- Extract negative control outcomes for instrument validation
-- Returns long-format table of person x negative control outcome presence
-- Parameters:
--   @cdm_database_schema: Schema containing OMOP CDM tables
--   @cohort_database_schema: Schema containing the cohort table
--   @cohort_table: Name of the cohort table
--   @outcome_cohort_id: Cohort definition ID for the main outcome
--   @negative_control_cohort_ids: Comma-separated list of negative control cohort IDs

SELECT
  mr.person_id,
  nc.cohort_definition_id AS outcome_cohort_id,
  CASE WHEN nc.subject_id IS NOT NULL THEN 1 ELSE 0 END AS has_outcome
FROM #mr_cohort mr
LEFT JOIN @cohort_database_schema.@cohort_table nc
  ON mr.person_id = nc.subject_id
  AND nc.cohort_definition_id IN (@negative_control_cohort_ids)
{@use_index_date == 1} ? {
  AND nc.cohort_start_date >= mr.index_date
}
;
