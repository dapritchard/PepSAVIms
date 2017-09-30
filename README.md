# Software for the prioritization of putative bioactive compounds

The PepSAVIms R package provides a collection of software tools used to
facilitate the prioritization of putative bioactive compounds from a complex
biological matrix. The package was constructed to provide an implementation of
the statistical portion of the laboratory and statistical procedure proposed in
[_The PepSAVI-MS pipeline for natural product bioactive peptide discovery_](http://dx.doi.org/10.1021/acs.analchem.6b03625), by
Kirkpatrick et al.


## Data analysis pipeline

The software in this package aims to perform the following steps, described in
more detail below.

1.  Consolidation
2.  Filtering
3.  Compound ranking


#### Consolidation

The mass spectrometry abundance data can optionally undergo two preprocessing
steps. The first step is a consolidation step: the goal is to to consolidate
mass spectrometry observations in the data that are believed to belong to the
same underlying compound. In other words, the instrumentation may have obtained
multiple reads of mass spectrometry abundances that in actuality belong to the
same compound - in which case we wish to attribute all of those observations to
a single compound.


#### Filtering

The second optional preprocessing step for the mass spectrometry abundance data
is a filtering step. The goal of the filtering step is to further reduce the
data set to focus on only those compounds that could plausibly be contributing
to the bioactivity area of interest. Furthermore, these criteria aim to filter
out some of the noise detected in the dataset. By filtering the candidate set
prior to statistical analysis, the ability of the analysis to effectively
differentiate such compounds is greatly increased.


#### Data analysis

Once the mass spectrometry abundance data has optionally undergone any
preprocessing steps, a statistical procedure to search for putative bioactive
peptides is performed. The procedure works by specifying the level of the L2
penalty parameter in the elastic net penalty, and tracking the inclusion of the
coefficients corresponding to compounds into the nonzero set along the elastic
net path. An ordered list of candidate compounds is obtained by providing the
order in which the coefficients corresponding to compounds entered the nonzero
set.


## Further information

Please see the R function documentation or the package vignettes for more
information regarding the use of this package
at [CRAN](https://CRAN.R-project.org/package=PepSAVIms).
