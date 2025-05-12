# Quantitative genetic analysis of early life viability resolves a case of the paradox of stasis in body size traits of Soay sheep

The code in this repository is to run the models presented in the manuscript above. In this study we assess quantitative genetic parameters to gain an evolutionary explanation as to why adult size traits are not evolving in the direction of current estimates of selection within this wild population of Soay sheep.

Elizabeth A. Mittell1,2,*, Josephine M. Pemberton1, Loeske E. B. Kruuk1 & Michael B. Morrissey2

1 - Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, UK; 2 - School of Biology, University of St. Andrews, UK

*Corresponding author: Elizabeth A. Mittell, e.mittell@ed.ac.uk, e.mittell@gmail.com; Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, EH9 2LD, UK

J.M.P, M.B.M, L.E.B.K, and E.A.M collected the data. E.A.M and M.B.M analysed the data. E.A.M is responsible for the code in this repository.

## Data and Scripts

CovuDataABMFYS_120525_RandomIDs.csv, CovuDataAMTLFYS_120525_RandomIDs, CovuDataMHLFYS_120525_RandomIDs, and CovuDataMSCFYS_120525_RandomIDs, contain the phenotypic data for the bivariate animal models in the first section of the analyses. The names ABM, AMTL, MHL, MCS and FYS refer to August body mass, August metatarsal length, male horn length, scrotal circumference (males) and first year survival. SheepPedigreeRandomIds_May2025.csv contains the pedigree for the Soay sheep population. The IDs within all these files have been randomised. This means that they should link to each other, but will not link to any other data published from the Soay sheep Project.

### Section one -- genetic signatures of prior viability selection
CovuDataABMFYS_120525_RandomIDs.csv columns are (same across these Covu files): id -- randomised ID; animal -- randomised ID that links to pedigree; time -- year of measurement for adult traits; y -- trait value; Twin -- whether the individual was a twin or singleton; Sex -- 1 indicates female, 2 indicates male; Horn -- 1 indicates scurred, 3 indicates normal horns; BirthYear -- birth year; PopDensBirth -- population size in the year of birth; MumID -- randomised maternal ID; jdate -- measurement Julian date for adult traits; PopDensCap -- population size in the year of measurement; CapAgeYears -- age of sheep; horngeno -- horn genotype; trait -- r indicates adult trait, s indicates first year survival; measure -- S1 indicates first year survival, others indicate adult trait name; family -- distribution of the variable specified for the covu models; MumAge -- age of the dam; AgeMC -- age mean centred within adults; jdateMC -- measurement Julian date mean centred within adults; AgeMC -- age mean centred within adults; MumAgeMC -- maternal age mean centred within juveniles; PopDensMS -- population density mean centred within juvenile and adult traits.

.R contains the models that were used in the manuscript. The packages used are shown within the script. These models were run in various versions of R. All run and are installable in R version 4.3.2 as of 5th April 2024 on macOS Monterey version 12.7.1.

### Section two -- phenotypic episodes of selection


### SOAY SHEEP PROJECT DATA REUSE STATEMENT:

The attached files contain data derived from the long term field project monitoring individual Soay sheep on St Kilda and their environment. This is a request to please let us know if you use them. Several people have spent the best part of their careers collecting the data. If you plan to analyse the data, there are a number of reasons why it would be very helpful if you could contact the authors or Professor Dan Nussey (dan.nussey@ed.ac.uk) before doing so.

[NB. If you are interested in analysing the detailed project data in any depth you may find it helpful to have our full relational database rather than the file(s) available here. If so, then we have a simple process for bringing you onto the project as a collaborator.]

1) The data can be subject to change due to updates in the pedigree, merging of records, occasional errors and so on.

2) The data are complex and workers who do not know the study system may benefit from advice when interpreting it.

3) At any one time a number of people within the existing project collaboration are analysing data from this project. Someone else may already be conducting the analysis you have in mind and it is desirable to prevent duplication of effort. Frequently these projects will be the only research of an early-career researcher.

4) In order to maintain funding for the project, every few years we have to write proposals for original analyses to funding agencies. It is therefore very helpful for those running the project to know what data analyses are in progress.

5) Individual identifiers may vary relative to other data archives from papers using the individual-level data. Therefore, contacting us is advised.

6) Data are here for the explicit purpose of reproducible research.
