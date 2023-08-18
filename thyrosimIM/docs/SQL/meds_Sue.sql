WITH dcd AS (
	SELECT DiagnosisID
	FROM rpt.DIM_DiagnosisCodeDim dcd
	WHERE DiagnosisCode IN ('E06.3', 'E03.9', 'E03', 'E03.3', 'E03.1') 
	--All non-cancer/surgery related hypothyroid patients
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

mof AS (
	SELECT minimum_dates.*
	  ,mof.PatientEncounterCSNID
	  ,mof.MedicationID
      ,mof.MedicationOrderID
      ,mof.OrderDescription AS DrugDescription
      ,mof.OrderDateTime
      ,mof.SigInstructions
      ,mof.PrescriptionQuantity
      ,mof.MedciationRouteCode
      ,mof.MedicationDoseUnitCode
      ,mof.Refills
      ,mof.OrderStatusCode AS DrugStatusCode
      ,mof.OrderingModeCode
      ,mof.DiscreteFrequencyID
      ,mof.DiscreteFrequencyName
      ,mof.OrderStartDate
      ,mof.OrderEndDate
      ,mof.IsCurrentActiveMedicationFlag
      ,mof.OrderingDate
      ,mof.PrescriptionWrittenDate
      ,mof.OrderClassCode
	  ,DATEDIFF(day, mof.OrderStartDate, mof.OrderEndDate) AS OrderDuration
  FROM rpt.CDM_MedicationOrderFact mof
  INNER JOIN minimum_dates on minimum_dates.PatientID = mof.PatientID
  WHERE
  mof.DiscreteFrequencyID NOT IN (200926, 200905, 200554, 200960, 200533, 40801002, 200902, 200913, 200537, 200914, 10010)
  OR mof.DiscreteFrequencyID IS NULL
  --mof.DiscreteFrequencyID IS NOT NULL
  --AND mof.PrescriptionQuantity IS NOT NULL
  --AND mof.Refills IS NOT NULL
  --AND mof.OrderStartDate IS NOT NULL
  --AND mof.OrderEndDate IS NOT NULL
),

md AS (
	SELECT mof.*
		  ,md.MedicationKey
		  ,md.MedicationName
		  ,md.GenericName
		  ,md.SimpleGenericName
		  ,md.GenericProductId
		  ,md.Strength
		  ,md.Form
		  ,md.Route
	  FROM rpt.DIM_MedicationDim md
	  INNER JOIN mof on mof.MedicationID = md.MedicationID
	  WHERE (md.MedicationName LIKE '%Thyroxine%'
	  OR md.MedicationName LIKE '%T4%'
	  OR md.SimpleGenericName LIKE '%Thyroid%'
	  OR md.SimpleGenericName LIKE '%Levothyroxine%')
	  AND md.SimpleGenericName NOT IN ('Parathyroid Hormone (Recomb)', 'Thyroid (Pork)')
	  AND md.Strength IS NOT NULL
	  AND OrderStartDate IS NOT NULL
	  AND StartDate > '2018-01-01'
),

multiDose AS (
	SELECT m1.PatientID
	FROM md m1
	JOIN md m2 ON m1.PatientID = m2.PatientID
	WHERE m1.OrderStartDate < m2.OrderEndDate
	  AND m1.OrderEndDate > m2.OrderStartDate
)

SELECT * FROM md
WHERE NOT EXISTS (
	SELECT 1
	FROM multiDose
	WHERE multiDose.PatientID = md.PatientID
	)
