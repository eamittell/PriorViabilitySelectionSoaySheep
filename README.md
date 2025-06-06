# Quantitative genetic analysis of early life viability resolves a case of the paradox of stasis in body size traits of Soay sheep

The code in this repository is to run the models presented in the manuscript above. In this study we assess quantitative genetic parameters to gain an evolutionary explanation as to why adult size traits are not evolving in the direction of current estimates of the response to selection within this wild population of Soay sheep.

Elizabeth A. Mittell1,2,*, Josephine M. Pemberton1, Loeske E. B. Kruuk1 & Michael B. Morrissey2

1 - Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, UK; 2 - School of Biology, University of St. Andrews, UK

*Corresponding author: Elizabeth A. Mittell, e.mittell@ed.ac.uk, e.mittell@gmail.com; Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, EH9 2LD, UK

J.M.P, M.B.M, L.E.B.K, and E.A.M collected the data. E.A.M and M.B.M analysed the data. E.A.M is responsible for the code in this repository.

## Data and Scripts
In the file names ABM, AMTL, MHL, MCS and FYS refer to August body mass, August metatarsal length, male horn length, scrotal circumference (males) and first year survival. SheepPedigreeRandomIds_May2025.csv contains the pedigree for the Soay sheep population. The IDs within all these files have been randomised. This means that they should link to each other, but will not link to any other data published from the Soay sheep Project.

All the dataframes have columns in common, which are: id/ID -- randomised ID; animal -- randomised ID that links to pedigree and is in the pedigree; dam -- randomised maternal ID in the pedigree that links to Mum ID in the data; sire -- randomised paternal ID in the pedigree; time/CapYear -- year of measurement for adult traits; y -- trait value; Twin -- whether the individual was a twin or singleton; Sex -- 1 indicates female, 2 indicates male; BirthYear -- birth year; jdateMC -- measurement Julian date mean centred; MumAgeMC -- maternal age mean centred; PopDensMC -- population density mean centred; MumID -- randomised maternal ID that links to dam in the pedigree; CapAgeYears -- age at time of measurement; trait -- r indicates adult trait, s indicates first year survival; measure -- S1 indicates first year survival, others indicate trait name; family -- distribution of the variable specified for the covu models; Weight -- August body mass (kg); HindLeg -- August metatarsal length (mm); HornLen -- Male normal horn length (mm); BolCirc -- scrotal circumference (mm); LambAugWeight -- Lamb August body mass (kg); Sur1 -- 1 is lamb survived to 1st May in year following year of birth, 0 is died before; birthjdateMC -- Julian birth date mean centred; LambHindLeg -- Lamb August metatarsal length (mm); LambHornLen -- male lamb normal horn length (mm); LambBolCirc -- lambe scrotal circumference (mm); PopDensBirth -- population density in year of birth; Mum Age -- maternal age; CapJDate -- Julian date of measurement.

### Section one -- genetic signatures of prior viability selection
_Prior viability selection_

CovuDataABMFYS_120525_RandomIDs.csv, CovuDataAMTLFYS_120525_RandomIDs.csv, CovuDataMHLFYS_120525_RandomIDs.csv, and CovuDataMSCFYS_120525_RandomIDs.csv, contain the phenotypic data for the bivariate animal models in the first section of the analyses.

_Response to selection using only information on adults_

DataUniVarABM_h2_RandomIDs_May2025.csv, DataUniVarAMTL_h2_RandomIDs_May2025.csv, DataUniVarMHL_h2_RandomIDs_May2025.csv, and DataUniVarMSC_h2_RandomIDs_May2025.csv, contain the phenotypic data for the univariate animal models used to estimate heritability of adult size traits. DataAdultPhenoSel_ABM_RandomIDs_May2025.csv, DataAdultPhenoSel_AMTL_RandomIDs_May2025.csv, DataAdultPhenoSel_MHL_RandomIDs_May2025.csv, and DataAdultPhenoSel_MSC_RandomIDs_May2025.csv, contain the phenotypic information used to estimate selection on adult size traits using information on adults alone. 

Section1_GeneticSignaturesPriorViabilitySelection.R contains the models that were used in section 1 of the manuscript. The packages used are shown within the script. These models were run in various versions of R. All run and are installable in R version 4.3.2 as of 13th May 2025 on macOS Monterey version 12.7.1.

### Section two -- phenotypic episodes of selection
_Viability selection lamb size traits_

DataPhenoSel_LABM_RandomIds_May2025.csv, DataPhenoSel_LAMTL_RandomIds_May2025.csv, DataPhenoSel_LMHL_RandomIds_May2025.csv, and DataPhenoSel_LMSC_RandomIds_May2025.csv, contain these phenotypic data for the bivariate models.

_Phenotypic covariance between homologous lamb and adult size traits_
DataPhenoCovu_ABM_LABM_RandomIds_May2025.csv, DataPhenoCovu_AMTL_LMTL_RandomIds_May2025.csv, DataPhenoCovu_AMHL_LMHL_RandomIds_May2025.csv, and DataPhenoCovu_AMSC_LMSC_RandomIds_May2025.csv, contain the phenotypic data for these bivariate models.

_Phenotypic covariance among lamb size traits_
DataPhenoCovarLambTraits_RandomIds_May2025.csv contains the phenotypic data for the multiresponse model.

Section2_PhenotypicEpisodesofSelection.R contains the models that were used in section 2 of the manuscript. The packages used are shown within the script. These models were run in various versions of R. All run and are installable in R version 4.3.2 as of 13th May 2025 on macOS Monterey version 12.7.1.

### SOAY SHEEP PROJECT DATA REUSE STATEMENT:

The attached files contain data derived from the long term field project monitoring individual Soay sheep on St Kilda and their environment. This is a request to please let us know if you use them. Several people have spent the best part of their careers collecting the data. If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact the authors or Professor Dan Nussey (dan.nussey@ed.ac.uk) before doing so.

[NB. If you are interested in analysing the detailed project data in any depth you may find it helpful to have our full relational database rather than the file(s) available here. If so, then we have a simple process for bringing you onto the project as a collaborator.]

1) The data can be subject to change due to updates in the pedigree, merging of records, occasional errors and so on.

2) The data are complex and workers who do not know the study system may benefit from advice when interpreting it.

3) At any one time a number of people within the existing project collaboration are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable to prevent duplication of effort. Frequently these projects will be the only research of an early-career researcher.

4) In order to maintain funding for the project, every few years we have to write proposals for original analyses to funding agencies. It is therefore very helpful for those running the project to know what data analyses are in progress.

5) Individual identifiers may vary relative to other data archives from papers using the individual-level data. Therefore, contacting us is advised.

6) Data are here for the explicit purpose of reproducible research.
