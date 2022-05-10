# Freesurfer-data-analysis

If nothing else is specified the models are run with group, sex and site as factors.
If nothing else is specified the models are run with age and eulernumber as covariates.
All scripts run separate models with and without a global covariate, e.g. ICV. 


## Cortical brain regions 

```freesurfer_analysis.R``` is an archive script with many different ways of running models and ways to do model diagnostics.
 
```freesurfer_analysis_bilateral.R``` is a script to run models with sex as factor on bilateral brain regions. 
  
```freesurfer_analysis_bilateral_sex_divided.R``` is a script to run models separately on each sex for bilateral brain regions. 

```freesurfer_analysis_global_behavioural_analysis.R``` is a script to run models on global measures with behavioural measures as covariates. Also includes interactions between groups and behavioural measures. The script also runs models on behavioural measures to find groups differences. 

```freesurfer_analysis_global_measures.R``` is a script to run models on global measures both with sex as factor and separately for each sex. Also runs multivariate anova models on global measures. Also estimates effect sizes for anova tests and post-hoc contrast tests. 

```freesurfer_analysis_global_measures_axis1.R``` is a script to run models on global measures with axis1 diagnosis as a factor.

```freesurfer_analysis_global_measures_euler.R``` is a script to run models on global measures with euler number as a covariate.

```freesurfer_analysis_lateral.R``` is a script to run models on lateral (separately on left and right hemisphere) brain regions. Runs models on sex as factor and also run models separately on each sex.


## Subcortical brain regions 

```freesurfer_analysis_lateral_subcortical.R``` is a script to run models on lateral (separately on left and right hemisphere) brain regions.
