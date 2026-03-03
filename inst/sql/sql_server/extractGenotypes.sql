-- Extract genotypes for MR instrument SNPs from OMOP Genomic Extension
-- OMOP Genomic CDM: https://github.com/OHDSI/Genomic-CDM
--
-- Uses the VARIANT_OCCURRENCE table from the OMOP CDM Genomic Extension.
-- Minimal required columns from VARIANT_OCCURRENCE: person_id, rs_id, genotype.
-- Also extracts reference_allele and alternate_allele for allele harmonization.
--
-- The genotype column in VARIANT_OCCURRENCE is VARCHAR and may contain VCF-style
-- values (e.g., "0/0", "0/1", "1/1") or plain integers ("0", "1", "2").
-- Conversion to integer allele dosage is performed in R after extraction.
--
-- Parameters:
--   @genomic_schema: Schema containing the VARIANT_OCCURRENCE table
--   @cohort_database_schema: Schema containing the cohort table
--   @cohort_table: Name of the cohort table
--   @outcome_cohort_id: Cohort definition ID for the outcome
--   @snp_ids: Comma-separated list of rs IDs to extract

SELECT
  vo.person_id,
  vo.rs_id AS snp_id,
  vo.genotype AS genotype_raw,
  vo.reference_allele,
  vo.alternate_allele
FROM @genomic_schema.VARIANT_OCCURRENCE vo
INNER JOIN (
  SELECT DISTINCT subject_id
  FROM @cohort_database_schema.@cohort_table
  WHERE cohort_definition_id = @outcome_cohort_id

  UNION

  SELECT DISTINCT person_id AS subject_id
  FROM #mr_cohort
  WHERE outcome = 0
) cohort_persons
  ON vo.person_id = cohort_persons.subject_id
WHERE vo.rs_id IN (@snp_ids)
;
