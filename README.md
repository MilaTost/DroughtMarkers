# DroughtMarkers
This repository contains scripts to analyze natural forestry populations. The European beech populations are located at five locations in the South-Eastern Romanian Carpathians along a steep environmental gradient associated with precipitation and temperature. A Genome-wide association studies (GWAS) was conducted based on this data set to identify SNP markers associated with different drought stress traits. The polygenic test [Ghat](https://academic.oup.com/genetics/article/209/1/321/5931021?login=true) was also tested in the same population. Additionally, the same data set was tested in an environmental association analysis (EAA) to test for associations with 53 environmental variables, which may have posed as selective pressure. These environmental variables comprised frost frequency change, temperature and precipitation. <br />
<img width="400" alt="Figure_1" src="https://github.com/user-attachments/assets/5e5e0dd0-ee9d-4bd9-a38b-aadbb622ea16"> <br /> <br />

## Genome-wide association studies (GWAS)
Genome-wide association studies (GWAS) aim to identify loci associated with drought stress, tree physiology, growth or wood quality traits, which could be prioritized in breeding programs. Genotypic and phenotypic data were collected from approximately 100 adult beech trees per stand in five locations in the South-Eastern Romanian Carpathians along an altitudinal gradient associated with precipitation and temperature. We performed GWAS using PLINK to identify SNP markers associated with traits related to drought stress. Additionally, permutation testing was conducted to determine significance thresholds. This procedure is described in the script `GWAS_with_PLINK.R` script. <br />
The plotting of the GWAS results (Manhattan plots) is described in the script `Plot_GWAS_PLINK_results.R`. <br /> <br />

## Environmental association analysis (EAA)
Environmental association analysis (EAA) focuses on the identification of environmental variables which may have acted as selective pressure. For this EAA, the [latent factor mixed model (LFMM)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12382) implemented in LEA (Landscape and Ecological Association Studies) R package was used. LFMM contains a mixed effect model, which includes correction methods for population structure. This procedure is described in the script `EAA_with_LFMM2_LEA.R`. <br />

To test for inflation of false positive results, the environemental data was randomized and tested for remaining associations. This procedure is available in the `Randomization_of_EAA_analysis.R` script. <br />

The creation of MAF plots is included in the scipt `Create_MAF_plots.R`.
<br /> <br />

## Polygenic test Ghat
The polygenic test Ghat is a novel approach to study initially the impact of selection in real life breeding programs. The theoretical and statistical details of this polygenic test were previously published by [Beissinger et al. (2018)](https://academic.oup.com/genetics/article/209/1/321/5931021?login=true). The introduction as an R package can be found in [Mahmoud et al. (2023)](https://academic.oup.com/g3journal/article/13/2/jkac319/6858947). Ghat considers all SNP markers together to retrieve a signal of polygenic selection as described by Fisher’s ‘infinitesimal model’ (1918). The application of Ghat in the forestry population can be found in the script `Calculate_polygenic_test_Ghat.R`. <br /> <br />
