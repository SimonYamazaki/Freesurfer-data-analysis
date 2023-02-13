# Freesurfer-data-analysis

Statistical models used to investigate group differences in freesurfer structural brain data from VIA11 project.

If nothing else is specified the models are run with group, sex and site as factors.
If nothing else is specified the models are run with age and eulernumber as covariates.
All scripts run separate models with and without a global covariate, e.g. ICV. 


## Cortical brain regions 

All scripts are run on cortical brain regions if nothing else is specified.

```/behav/``` is statistical models run for with behavoiural measures. Both models with the behavoiural measure as covariate and as response variable is included. 

```/bilateral/``` models run on bilateral brain regions, that is the average of the two hemispheres for each brain region. There are 34 brain regions on each hemisphere.

```/global_measures/``` models run on global brain measures.

```/lateral/``` models run on lateral brain regions, that is separately run models for each brain region on each hemisphere.

```/normative_modelling/``` investigation of normative modelling results on the VIA11 brain data. Very explorative and not well documented.  

```/other/``` Other.


## Subcortical brain regions 

```lateral/freesurfer_analysis_lateral_subcortical.R``` is a script to run models on lateral (separately on left and right hemisphere) brain regions.


## Note on residuals

The residuals of each model can be extracted by extracting the residuals attribute from the model object as so: ```model_object["residuals"]```




