## Scripts for Imputation project

In this repository you will find all the scripts I used to write my thesis.

There are four main steps in the analysis:

1. Creating the NATS imputation reference panel (panel_creation folder)
2. Running imputation and evaluating performance (imp_performance folder)
3. Demographic simulations (demographic_sims folder)
4. Figure creation (figure_creation folder)

Now I will detail the contents of each folder

### Panel creation

In this folder you will find the sprits used to clean and merge the different data sources.

+ clean[data source name].sh: These scripts removes indels, non-biallelic snps and missing positions in each data source. It also saves the removed positions so they can be also removed from the final panel.
+ merge_hgdp_sgmx.sh: This script merges HGDP, SGDP and INMEGEN data into a single panel using bcftools.
+ merge_mx50_hgdp.sh: This script merges HGDP-SGDP-INMEGEN with the 50 Mexican genomes using bcftools.

### Running imputation and evaluating performance

In this folder you will find the necessary scripts for running and evaluating imputation performance using different reference panels.

+ create_hapref.py: Script used to create the leave-one-out 1KGP reference panels in hapref format.
+ impuatation_pipe_alt.py: This is the master script that controls the imputation run. It creates and launches other relevant scripts for the pipeline. It uses a config file to specify the mode of running (one or two reference panels), the paths to relevant data and scripts and where to save the results. An example of the config file can be found under the name config.example.txt
+ bind_and_annotate_ancestry.R: This script binds the imputation chunks and annotates the local ancestry of each variant.
+ make_summaries.R: This script calculates the r2 between imputed and real allele dosages by allele frequency bin and local ancestry.
+ infomaf_sum.R: This script summarizes the info scores and maf of each imputed variant in each individual.
+ plotnfilt.R: This script plots the results from infomaf_sum.R.
+ compute_concordance: This script calculates concordance between the most likely imputed allele and the real allele.

### Demographic simulations

In this folder you will find the scripts used to generate the demographic simulation data and the subsequent analysis.

+ make_data.sh: Runs the command from stdpopsim to create the simulated data.
+ make_natrefpanel.sh: Creates the different imputation reference panels with varying number of native individuals (0, 100, ..., 3000)
+ subsamp60k.R: Simulate a SNP array based on the MEGA array frequencies in the 1KGP European population.
+ get_varfreqs: Script to calculate the variant frequencies of the simulated data.
+ imp[number].sh: Scripts to run imputation with the different references.

### Figure creation

In this folder you will find the scripts used to generate the figures.

+ plot_natpancomp.R: Plot the imputation r2 results with the two references.
+ plot_supp.R. Supplementary plots
+ simulation_plot: Plot the results of the demographic simulations section
+ create_fig1.R: Includes NATS admixture plot, map of the samples and SNP category in 1KGP + NATS
