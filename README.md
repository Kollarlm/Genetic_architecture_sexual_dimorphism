# Genetic_architecture_sexual_dimorphism
All code corresponds to https://royalsocietypublishing.org/doi/full/10.1098/rspb.2020.2908

## Citation for dryad dataset
Kollar, Leslie (2021), Moss growth, development, morphology, and physiology dataset and code, Dryad, Dataset, https://doi.org/10.5061/dryad.59zw3r266

## Abstract
A central problem in evolutionary biology is to identify the forces that maintain genetic variation for fitness in natural populations. Sexual antagonism, in which selection favors different variants in males and females, can slow the transit of a polymorphism through a population or can actively maintain fitness variation. The amount of sexually antagonistic variation to be expected depends in part on the genetic architecture of sexual dimorphism, about which we know relatively little. Here, we used a multivariate quantitative genetic approach to examine the genetic architecture of sexual dimorphism in a scent-based fertilization syndrome of the moss  Ceratodon purpureus. We found sexual dimorphism in numerous traits, consistent with a history of sexually antagonistic selection. The cross-sex genetic correlations (rmf) were generally heterogeneous with many values indistinguishable from zero, which typically suggests that genetic constraints do not limit the response to sexually antagonistic selection. However, we detected no differentiation between the female- and male-specific trait (co)variance matrices (Gf and Gm, respectively), meaning the evolution of sexual dimorphism may be constrained. The cross-sex cross-trait covariance matrix B contained both symmetric and asymmetric elements, indicating that the response to sexually antagonistic or sexually concordant selection, and the constraint to sexual dimorphism, is highly dependent on the traits experiencing selection. The patterns of genetic variances and covariances among these fitness components is consistent with partly sex-specific genetic architectures having evolved in order to partially resolve multivariate genetic constraints (i.e. sexual conflict), enabling the sexes to evolve toward their sex-specific multivariate trait optima.

## Methods
This MossQuantGenreadme.txt file was generated on 2021-02-05 by Leslie M. Kollar

#GENERAL INFORMATION
1. Moss growth, development, morphology, and physiology dataset and code:

2. Author Information

               A. Principal Investigator Contact Information
                               Name: Leslie M. Kollar
                               Institution: Michigan State University
                               Email: lesliemkollar@gmail.com

               B. Associate or Co-investigator Contact Information
                               Name: Karl Grieshop
                               Institution: Stockholm University
                               Email: karlgrieshop@gmail.com

3. Date of data collection (single date, range, approximate date): Spring 2016

4. Geographic location of data collection: Portland, OR (Portland State University)

5. Information about funding sources that supported the collection of the data: This work was supported by the National Science Foundation Doctoral Dissertation Improvement Grant (NSF DEB 1701915) to LMK and SFM, NSF grants to SFM (DEB 1541005 and 1542609); EDEN: Evo-Devo-Eco Network Training Grant to LMK, MicroMorph Cross-Disciplinary Training Grant to LMK, the University of Florida’s Biology Department grants to LMK, and by the Swedish Research Council (2018-06775 to KG).

## DATA & FILE OVERVIEW

File List:
LH.traits.data.NONA.csv
Data.homemade.traits.csv
Multivariate_analyses_Kollar.L.R
Univariate_analyses_Kollar.L.R
scaled.threeLHtraits.parexp.obj
Final.VOC.4.repro.obj
Relationship between files, if important:
 
LH.traits.data.NONA.csv contains traits categorized as “growth and development” and Data.homemade.traits.csv contains traits categorized and “morphology and physiology”. We fit models for each category of traits and included both of the models. Scaled.threeLHtraits.parexp.obj is the model fit for the growth and development traits while Final.VOC.4.repro.obj is the model for the morphology and physiology traits. The majority of analyses can be found in Multivariate_analyses_Kollar.L.R while testing for significant genetic variation can be found in the file Univariate_analyses_Kollar.L.R.
 
## Methods for processing the data:

Data collection and processing is discussed in detail in the manuscript and the supplemental methods. We include some brief points here to clarify.

PTR Data: We identified 75 different masses using the PTR-TOF-MS in mature sex expressing gametophytes. We represented the masses as total volatile output and as number of compounds. The raw PTR data files were analyzed using the PTR viewer and background/blank cuvette air was subtracted from the sample readings to account for noise in the signal.

Leaf traits: After leaves were mounted flat on a slide and images were taken, leaf traits including (length, area, and perimeter) were calculated using a custom script in ImageJ.

Reproductive output: For female reproduction we counted the number of archegonia (eggs) in each female sex structure. To account for possible differences in placement of the sex structure along the stem, we sampled sex structures from three different locations on the stem (top, middle, and bottom) per sample. For male reproduction, we counted the total number of male sex structures per 10 stems per sample. To make the male and female reproductive units comparable we combined male and female reproductive units into a single column and mean centered and variance standardized reproductive units. Standardization methods are in the R scripts.

Growth and developmental data: We collected many traits encompassing growth and development in juvenile moss tissue. These samples were grown in a growth chamber in 12 well plates. Images were taken every 7 days and analyzed using ImageJ software. Some measurements such as total number of gametophores were measured simply by counting the presence of mature gametophytes.

## Instrument- or software-specific information needed to interpret the data:

Proton Transfer Reaction Time of Flight Mass Spectrometer (PTR-TOF-MS 1000, Ionicon)
PTR-MS Viewer 3.1 (Ionicon)
R (version 4.0.2; R Development Core Team 2020)
MCMCglmm’ (v. 2.29)
