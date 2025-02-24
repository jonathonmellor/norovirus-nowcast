# Norovirus Nowcast England 2023/24

This is the repository for the code and data associated with the paper entitled: ["An Application of Nowcasting Methods: Cases of Norovirus during the Winter 2023/2024 in England"](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012849) now published in PLOS Computational Biology, https://doi.org/10.1371/journal.pcbi.1012849. 

## Authors:

Jonathon Mellor, Maria L. Tang, Emilie Finch, Rachel Christie, Oliver Polhill, Christopher E. Overton, Ann Hoban, Amy Douglas, Sarah R. Deeny, Thomas Ward.

There are two key points to consider when exploring this repository:

1. the data shown are partially synthetic, based on the true data with added statistical noise. This is done to help protect individual data, and is a requirement for us to open source UKHSA data appropriately. Data processing for individual line list data have not been included in this repository.
2. some nowcasting models are highly computationally expensive. This research was developed in a high compute environment, therefore decisions around parallelisation, may not be appropriate for your computing system. Please run models with caution and once you are comfortable with the relevant packages involved.



## scripts

- `depends.R` - all package dependencies

### data_processing

Files in this director have been removed as a condition for open sourcing.

### descriptive

- `figure1.R` - generates descriptive time series plots
- `figure2.R` - generates plots on reporting delays
- `111_online.R` - generates time series of 111 data.

### run_models

- `run_gam.R` - runs the GAM nowcast and saves locally
- `norovirus_nowcast_config.yaml` - config and model parameters (currently used only for GAM model)
- `run_bsts_model.R` - runs the BSTS model and saves locally
- `run_bsts_111_online_model.R` - runs the BSTS model using NHS 111 online data and saves locally
- `run_epinowcast.R` - runs the epinowcast model and saves locally

#### exploration

- `gam_variations.R` - run and score different variations on the GAM norovirus nowcast
- `111_online_collinearity_analysis.R` - check the correlation of 111 online variables
- `tune_*.R` - run the tuning scripts for each model

#### functions

- `gam.R` - function for running GAM model
- `model_running_functions.R` - assorted functions needed for running models
- `plotting.R` - plotting functions
- `scoring.R` - scoring functions
- `epinowcast.R` function for running epinowcast model

### compare_models

- `plot_and_score_all_models.R` - loads in prediction dataframes from all models and scores and plots them together
- `plot_tuning_scores.R` - produce scoring tables


## outputs

### data

Contains the 111 and cases data used in this research.

### tuning

Outputs relating to model tuning are stored here.

### plots

Visualisations for the paper are stored here.

### scoring

Outputs relating to model scoring are stored here.
