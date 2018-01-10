# Parameters For Sample Files

## Description
	EBprotWithRatioConversion uses EBprotWithConv_Sample.txt as dataset
	EBprotV2 uses EBprotV2_Sample.txt as dataset

## EBprotWithRatioConversion Parameters
	Working Directory: <preferred directory where results should be saved>
	Experimental Design: Independent
	Data Input Form: log_2
	[x] Outlier Filtering:
		Min K: 5
		Bayesian False Discovery Rate Cutoff: 0.05
	Lower Bound: 0.20
	Upper Bound: 0.80
	Min # of Peptides: 1
	Min # of Samples: 4
	Peptide Identifier Column: Sequence
	Protein Indentifier Column: Proteins
	# of Groups: 3
		Group 1 Label: ERPR
		Group 1 Data: 
			BC1
			BC14
			... (from BC1 to BC5 in the given order)
			BC5
		Group 2 Label: HER2
		Group 2 Data:
			BC10
			BC13
			... (from BC10 to BC9 in the given order)
			BC9
		Group 3 Label: TN
		Group 3 Data:
			BC12
			BC15
			... (from BC12 to BC46 in the given order)
			BC46
		Group 1 Contrast: None
		Group 2 Contrast: Group 3/Group 2
		Group 3 Contrast: Group 3/Group 1 

Note that choosing `Group 1 Contrast: Group 3/Group 1` and `Group 3 Contrast: None` gives the same result since the same contrast is made.

## EBprotV2 Parameters
	Working Directory: <preferred directory where results should be saved>
	Experimental Design: Independent
	Data Input Form: log_2
	[x] Outlier Filtering:
		Min K: 5
		Bayesian False Discovery Rate Cutoff: 0.05
	Lower Bound: 0.20
	Upper Bound: 0.80
	Min # of Peptides: 1
	Peptide Identifier Column: Sequence
	Protein Indentifier Column: Proteins
	# of Labels: 2
		Label 1 Name: TNvsERPR
		Label 1 Data: 
			TN/ERPR
		Label 2 Name: TNvsHER2
		Label 2 Data:
			TN/HER2

		
