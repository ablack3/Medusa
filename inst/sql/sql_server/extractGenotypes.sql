-- Extract genotypes for MR instrument SNPs from genomic linkage table
-- Parameters:
--   @genomic_linkage_schema: Schema containing the genomic linkage table
--   @genomic_linkage_table: Name of the genomic linkage table
--   @genomic_person_id_column: Column name for person identifier (default: person_id)
--   @cohort_database_schema: Schema containing the cohort table
--   @cohort_table: Name of the cohort table
--   @outcome_cohort_id: Cohort definition ID for the outcome
--   @snp_ids: Comma-separated list of SNP IDs to extract

SELECT
  g.@genomic_person_id_column AS person_id,
  g.snp_id,
  g.genotype
FROM @genomic_linkage_schema.@genomic_linkage_table g
INNER JOIN (
  SELECT DISTINCT subject_id
  FROM @cohort_database_schema.@cohort_table
  WHERE cohort_definition_id = @outcome_cohort_id

  UNION

  SELECT DISTINCT person_id AS subject_id
  FROM #mr_cohort
  WHERE outcome = 0
) cohort_persons
  ON g.@genomic_person_id_column = cohort_persons.subject_id
WHERE g.snp_id IN (@snp_ids)
;
