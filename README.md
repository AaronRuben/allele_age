# Estimating the age of the _HOXB13_ X285K variant in African ancestry men
The variant is almost exclusively found in West Africa at a frequency of ~0.5% (check out the publication [here](https://doi.org/10.1016/j.eururo.2021.12.023)).
We want to estimate its age.

- The observed allele frequency comes with some uncertainty depending on the sample size n. We can use a Binomial proportion confidence interval to determine the 95% CI.

- Slatkin & Rannala (2000) provide the theory to estimate the age of an allele in terms of 2N from its frequency.

- To estimate the distribution of the effective population size, we can use the heterozygosity observed in the Yoruba and Esan population.
The heterozygosity is defined as H = theta / (theta + 1).Thus, it allows us to compute theta. Making assumptions about the mutation rate u we can solve theta = 4Nu for N.
So that we need to extract all Yoruba and Esan samples from the 1000 Genome project by running the following script:

<code>./extract_yoruba_esan_1000G_phase_3.sh <path_to_1000G_phase3_vcf_files></code>

- Using these samples, we can determine the heterozygosity and the effective population size with the following script:

<code>./estimate_heterozygosity_plink.sh</code>

- The age of the allele also depends on the generation time. Here, we assume a generation time of 25 years.

- Sampling from the different distributions allows us to estimate the allele ages.

To estimate the age of an allele with frequency p, which is derived from n samples.
The estimate is refined by the probability of not observing a copy in the Ugandan cohort of size u if the variant arose before the Bantu migration b years ago.
The allele age is scaled by the generation time g:

<code>./estimate_allele_age.py -p p -n n -u u -b b -g g</code>

In the case of the YRI and ESN population in Ibadan p=0.0037, n=1071. The Bantu migration occurred b=3000 years ago, the Ugandan cohort size is u=677,
and g is assumed to be 25 years.
