# SIR-Ebola-DA-MCMC-Code

MATLAB code supporting the manuscript:

**"Bayesian and Frequentist Inference for Partially Observed Stochastic SIR Epidemic Models with Application to Ebola in Sierra Leone"**

## Description

This repository contains the MATLAB code and data supporting the simulation study and the real-data application presented in the manuscript.

The implemented methodology focuses on statistical inference for partially observed stochastic SIR epidemic models. In particular, the repository includes:

- Simulation of stochastic SIR epidemic trajectories using Gillespie's stochastic simulation algorithm.
- Bayesian inference using a data-augmentation Markov chain Monte Carlo (DA-MCMC) algorithm.
- Reconstruction of latent infectious trajectories under partial observation.
- Estimation of the basic reproduction number \(R_0\).
- Frequentist estimation for comparison with the Bayesian approach.
- Application to data from the 2014–2016 Ebola outbreak in Sierra Leone.

## Repository Structure

SIR-Ebola-DA-MCMC-Code/
│
├── README.md
├── LICENSE
│
├── simulation/
│   ├── simulation_SIR_DA_MCMC.m
│   └── gillespie_SIR.m
│
└── real_data/
    ├── DA_MCMC_Ebola_Sierra_Leone.m
    ├── ebola_daily_data_full.xlsx
    └── sierra_leone_real_values.xlsx

## File Description

### Simulation

- **`simulation/simulation_SIR_DA_MCMC.m`**  
  Main MATLAB script for the simulation study of the partially observed stochastic SIR epidemic model. It implements the DA-MCMC inference procedure and estimates the basic reproduction number \(R_0\).

- **`simulation/gillespie_SIR.m`**  
  MATLAB implementation of Gillespie's stochastic simulation algorithm for generating sample paths of the stochastic SIR epidemic model.

### Real-Data Application

- **`real_data/DA_MCMC_Ebola_Sierra_Leone.m`**  
  Main MATLAB script for the statistical analysis of the 2014–2016 Ebola outbreak in Sierra Leone using the proposed inference methodology.

- **`real_data/ebola_daily_data_full.xlsx`**  
  Data file containing the Ebola epidemic observations used in the real-data analysis.

- **`real_data/sierra_leone_real_values.xlsx`**  
  Processed data file used as input for the statistical analysis of the Sierra Leone Ebola outbreak.

## Software Requirements

The code was developed in MATLAB.

Required software:

- MATLAB

Additional MATLAB toolboxes may be required depending on the functions used in the scripts.

## Running the Code

### Simulation Study

Open MATLAB, navigate to the `simulation` directory, and run:

matlab - simulation_SIR_DA_MCMC


The script uses the function:

```matlab
gillespie_SIR


to generate stochastic SIR epidemic trajectories.

### Ebola Real-Data Analysis

Open MATLAB, navigate to the `real_data` directory, and run:

matlab    DA_MCMC_Ebola_Sierra_Leone


Make sure that the following data files remain in the same `real_data` directory:


ebola_daily_data_full.xlsx
sierra_leone_real_values.xlsx


## Reproducibility

The numerical results may exhibit small differences between runs because the simulation and MCMC procedures rely on random-number generation. For exact reproducibility, a fixed random seed may be specified in MATLAB before running the analysis, for example:

```matlab
rng(1);
```

## Data Availability

The data files used for the Ebola application are provided in the `real_data` directory of this repository. Details regarding the original data sources and preprocessing procedures are described in the associated manuscript.

## Code Availability

The MATLAB source code supporting the simulation study, the DA-MCMC algorithm, and the real-data analysis is publicly available in this repository. An archived version with a persistent DOI will be made available through Zenodo.

## Citation

If you use this code, please cite the associated manuscript:

> H. El Maroufy, A. Merbouha, and A. Lagzini,  
> *Bayesian and Frequentist Inference for Partially Observed Stochastic SIR Epidemic Models with Application to Ebola in Sierra Leone.*

Publication details and DOI will be added upon publication.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions regarding the code or the associated manuscript, please contact:

**Abdelati Lagzini**  
Sultan Moulay Slimane University  
Email: abdelati.lagzini@usms.ma
