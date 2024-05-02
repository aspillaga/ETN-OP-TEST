# ETN-OP-TEST
#### *Data and scripts from OP interoperability tests conducted in European waters*

From the manuscript: ***Open Protocols, the new standard for acoustic tracking: results from interoperability and performance tests in European waters***

Authors: 
Eneko Aspillaga, Stijn Bruneel, Josep Alós, Pieterjan Verhelst, David Abecasis, Kim Aarestrup, Kim Birnie-Gauvin, Pedro Afonso, Miquel Palmer and Jan Reubens

---

### Contents:

- ***./data/***: directory with the range test data. Within "./data/range_test/etn_data/", a script to download the data from the ETN database is included ("download_etn_data.R"). The outputs of all the scripts are also stored in this directory.
- ***./plots/***: directory where the plots of all the scripts are stored.
- ***"01_range_tests_prepare_data.R"***: Script to prepare the acoustic range data for the Bayesian analysis.
- ***02_range_test_bayesian_model.R***: Script to apply the Bayesian range model to the detection data.
- ***03_range_test_false_detection_analysis.R***: Script to analyze the occurrence of false detections.
- ***04_smolt_migration_analysis_OPi_vs_R64K.R***: Script to analyze the OPi vs R64K protocol comparison in migrating smolts.

