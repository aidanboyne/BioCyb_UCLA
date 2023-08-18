WITH dcd AS (
	--Select all diagnosis events where patient is diagnosed with Hashimoto's Disease
	SELECT DiagnosisID
	FROM rpt.DIM_DiagnosisCodeDim dcd
	WHERE DiagnosisCode = 'E06.3'
),

dpf AS (
	SELECT dcd.*,
	dpf.PatientID,
	dpf.StartDate
	FROM rpt.CDM_DiagnosisProblemFact dpf
	INNER JOIN dcd ON dcd.DiagnosisID = dpf.DiagnosisID
	WHERE dpf.PatientEncounterCSNID IS NOT NULL
	AND dpf.StartDate IS NOT NULL
),

minimum_dates AS (
	SELECT DISTINCT PatientID, MIN(StartDate) AS StartDate
	FROM dpf
	GROUP BY PatientID
),

-- At this point, we have the ID number for all patients
-- diagnosed with Hashimoto's disease, and the date at which
-- this diagnosis began

pd AS (
	SELECT minimum_dates.*,
	pd.Sex FROM rpt.DIM_PatientdIM pd
	INNER JOIN minimum_dates ON minimum_dates.PatientID = pd.PatientID
	WHERE pd.IsValid = 1
),

lof AS (
	SELECT pd.*, lof.PatientEncounterCSNID, lof.OrderID, lof.OrderDescription, lof.ProcedureID, lof.ProcedureStartDateTime
	FROM rpt.CDM_LabOrderFact lof
	INNER JOIN pd ON pd.PatientID = lof.PatientID
	WHERE lof.ProcedureID IN (304294, 906, 107749, 107737, 304294, 107753, --TSH
							  902, 900, 107887, 64850, 108099, 		      --FT4
							  922, 108081,64740, 108093,			      --FT3
							  920,		     --T3Total
							  108103,	     --T4Total
							  17791, 108131, 103233, 109925, 246978, 109927, --TGAb
							  64758, 108139, 109935,	     --TPOAb
							  246074,	     --PlasmaCell
							  244183,104329,	     --B&TCells
							  56290,	     --TRAb
							  1774,		     --IL6
							  55102)	     --PlasmaCytokines)
	 AND OrderStatusCode IN (3, 5) --Order Was completed
),

lorf AS (
	-- Pull the value and units for the labs ordered
	SELECT lof.*, lorf.OrderResultValue, lorf.ReferenceUnit
	FROM rpt.CDM_LabOrderResultFact lorf
	INNER JOIN lof ON lof.OrderID = lorf.OrderID
	WHERE lorf.OrderResultValue IS NOT NULL
	AND lorf.ReferenceUnit != 'ng/mL' AND lorf.ReferenceUnit != '% pos'
),

ef AS (
	-- Get date of this encounter, patient birthday and BMI
	SELECT lorf.*,
	ef.EncounterContactDate,
	ef.PatientBirthDate,
	ef.PatientBMI
	FROM rpt.CDM_EncounterFact ef
	INNER JOIN lorf ON lorf.PatientEncounterCSNID = ef.PatientEncounterCSNID
	WHERE ef.EncounterContactDate IS NOT NULL
),

allResult AS (
-- Data for all patients
SELECT DISTINCT 
	ef.PatientID, ef.PatientEncounterCSNID, DateDiff(year, ef.PatientBirthDate, ef.EncounterContactDate) AS Age,
	ef.sex, ef.PatientBMI, ef.StartDate AS DiagnosisDate, ef.EncounterContactDate AS LabDate,
	ef.OrderDescription, ef.ProcedureID, ef.OrderResultValue, ef.ReferenceUnit AS Units
FROM ef
WHERE ef.StartDate > '2012-11-02 00:00:00'
),

qualifiedPatient AS (
-- Get all qualified patientIDs that have all three types of orders we need
SELECT DISTINCT a1.PatientEncounterCSNID
FROM allResult a1
JOIN allResult a2 ON a1.PatientEncounterCSNID = a2.PatientEncounterCSNID
WHERE a1.OrderDescription = 'THYROID PEROXIDASE ANTIBODY'
  -- AND a2.OrderDescription = 'THYROGLOBULIN ANTIBODY'
  AND EXISTS (SELECT 1 FROM allResult WHERE PatientEncounterCSNID = a1.PatientEncounterCSNID AND OrderDescription LIKE 'TSH%')
  AND EXISTS (SELECT 1 FROM allResult WHERE PatientEncounterCSNID = a1.PatientEncounterCSNID AND OrderDescription LIKE 'Free T4%')
)
/*
We now get all patient that have all three (TPO, TSH, Free T4) lab results
For TSH%, it includes 'TSH', 'TSH WITH REFLEX FT4, FT3', and any other order starts with TSH
For T4%, it includes 'FREE T4', 'FREE T4 INDEX PROFILE', etc.
Let me know if any of these lab name should not be included 
(I remembered that FREE T4 INDEX PROFILE should not be included? If yes, I'll need to work on it)
*/

SELECT allResult.*
FROM allResult
INNER JOIN qualifiedPatient ON allResult.PatientEncounterCSNID = qualifiedPatient.PatientEncounterCSNID
	-- Only data for qualified patients left
WHERE allResult.PatientEncounterCSNID IN (SELECT PatientEncounterCSNID FROM allResult WHERE PatientBMI IS NOT NULL)
	-- includes all patients that have at least one BMI recorded
  AND allResult.LabDate > '2018-01-01' -- Lab performed after 2018
ORDER BY allResult.PatientEncounterCSNID
	-- 104960 rows left when using PatientID
	-- 11772 rows left when using PatientEncounterCSNID
	-- Number of rows: PatientID > PatientEncounterCSNID. more restrictive..?