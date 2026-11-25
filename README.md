Supplementary material for the article:

#### Global and Local Sensitivity Analysis of a Minimal T-Cell Activation Model to Inform Experimental Design for Cytokine Release
Authors: Darian Elia Tennessee Schirr<sup>1,2</sup>, Jan-Erik Busse<sup>1</sup>, Andreas Reichel<sup>1</sup>

Affiliations: 

<sup>1</sup> Preclinical Modeling & Simulation, Preclinical Development, Bayer AG, Berlin, Germany

<sup>2</sup> Institute of Pathology, Charité-Universitätsmedizin Berlin, Berlin, Germany

## Included files

#### (i) Complete Morris analysis code for the T-cell activation model 

#### (ii) Complete Sobol analysis code for the T-cell activation model 

#### (iii) Shiny app code for the T-cell activation
Explore parameter interdependencies through dose-response curves by simultaneous adjustment of multiple parameters.
To display pictures, download www folder with both png files and place in the same directory as the Shiny app code

#### (iv) local sensitivity analysis

## How to reproduce
Run each script top to bottom in a fresh R session. Code parts within a script depend on earlier parts. Different scripts do not depend on each other.

#### (i) Part 1: Morris indices with 95% confidence intervals (Fig. 2, Fig. S3, Table S3). Part 2: µ* over trajectories plot (Fig. S2).

#### (ii) Part 1: Si, STi and Dij from 25 Owen-scrambled replicates, with bootstrap-t confidence intervals (Fig. 3). Part 2: STi over sample size (Fig. S4). Part 3: Robustness to sampling ranges and input distributions. Especially part 2 is computationally demanding, it can be skipped to execute part 3.

#### (iii) Open in RStudio and click Run App (Figs. S5–S7).

#### (iv) Single and pairwise perturbations (Fig. 4, Fig. 5, Tables S4–S5).

## License 
This project is licensed under the terms of GPLv3.

