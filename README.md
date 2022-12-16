# A comparison of machine learning and statistical species distribution models applied to freshwater macroinvertebrates in Swiss rivers

This data package provides the data, scripts and plots needed to produce the results of the manuscript *"A comparison of machine learning and statistical species distribution models: when overfitting hurts interpretation"*. In this study, we apply eight statistical and machine learning models with different complexity to predict the probability of occurrence of freshwater macroinvertabrates in Swiss rivers using nine environmental factors as explanatory variables. We compare the models in terms of predictive performance, overfitting degree and inferred response shape, during cross-validation and for generalization out of the calibration domain ("extrapolation"). 

*Authors:* Emma Chollet Ramampiandra (ECR), Andreas Scheidegger (AS), Jonas Wydler (JW), Nele Schuwirth (NS)
*Correspondence:* emma.chollet@eawag.ch

-----------------------------------------------------------------------------------------------------------------------------------

## Data

The input data is already pre-processed.
- The macroinvertebrate dataset comes from the MIDAT database, run by info fauna CSCF (Centre Suisse de Cartographie de la Faune) & karch. It consists of 2729 observations of the presence (1) and absence (0) of 60 taxa taken at 1802 different sites covering the whole of Switzerland.
- The environmental dataset consists of nine environmental factors selected based on expert knowledge and previous studies (see Manuscript for more details). They are derived from the Swiss Federal Office for the Environment (water quality monitoring data, hydrological data, hydromorphological data, land use data), the Swiss Federal Office of Topography (topographical data), and the Swiss Federal Statistical Office (population statistics).

-----------------------------------------------------------------------------------------------------------------------------------

## Models application and analysis

**Main script** (Directory: *R_scripts*)
- *main.r* : Is the entry point to the code. It automatically instals all R packages in the right version using `checkpont`.

**Utilities** (Directory: *R_scripts*)
- *stat_model_functions.r* : Functions to run hierachical statistical models (hGLM and chGLM) using `rstan` package.
- *ml_model_functions.r* : Functions to train statistical and machine learning models (iGLM, GAM, SVM, BCT, RF) using `caret` package (and specific algorithms packages).
- *ann_model_functions.r* : Functions to train multilayer perceptron (ANN) using `tensorflow` and `keras` packages.
- *plot_functions.r* : Functions to produce all data and models analysis.
- *utilities.r* : Utilities functions used in all scripts.

**Input data** (Directory: *Input_data*)
- */Swiss.map.gdb/* * : GDB files for plotting rivers and lake on Swiss map
- *All_2729samples_9envfact_lme.area.elev_ModelInputs.csv* : Models input data (2729 observations of presence [1] and absence [0] of 60 taxa and nine environmental factors)
- *All_2729samples_9envfact_lme.area.elev_PrevalenceTaxa.csv* : Information on prevalence of the 60 taxa