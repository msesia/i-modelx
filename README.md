# Local model-X inference

This repository contains R and Python code accompanying the following paper:


[Searching for local associations while controlling the false discovery rate](https://arxiv.org/abs/2412.02182)


## Paper abstract

We introduce  local conditional hypotheses that express how the  relation between explanatory variables and outcomes changes across different contexts, described by covariates. By expanding upon the model-X knockoff filter, we show how to adaptively discover these local  associations, all while controlling the false discovery rate. Our enhanced inferences can help explain sample heterogeneity and uncover interactions, making better use of the capabilities offered by modern machine learning models. Specifically, our method is able to leverage any model for the identification of data-driven hypotheses pertaining to different contexts. Then, it rigorously test these hypotheses without succumbing to selection bias. Importantly, our approach is efficient and does not require sample splitting. We demonstrate the effectiveness of our method through  numerical experiments and by studying the genetic architecture of Waist-Hip-Ratio across different sexes in the UKBiobank.

## Reproducing the results 

This repository contains the following folders: 

- **genetic_applications**: This folder contains instructions and all scripts required to reproduce the results from *section 3.2, 3.3, appendix A6 and A7*. The folder does not contain any data files. We have applied for the UK Biobank data and interested parties can [apply to the UK Biobank](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access) for data access. 

- **experiments**: This folder contains the required scripts to reproduce the synthetic experiment from *section 3.1* and the transfer experiment from *appendix A8.2* (subfolder *synthetic*)

- **tutorials**: This folder contains two tutorials: One first tutorial demonstrates the basic usage of the adaptive Local Knockoff Filter. The second tutorial demonstrates its inner workings. 

- **LKF**: Contains the function to run the aLKF on synthetic data and as demonstrated in the tutorials. 
