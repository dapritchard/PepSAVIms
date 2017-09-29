
#' Bioactivity data
#'
#' The relative relative growth inhibition of bioactivity levels for the
#'     bacteria and virus strains studies in Kirkpatrick et al. (2016).
#'
#' @docType data
#'
#' @usage data(bioact)
#'
#' @format A \code{list} containing relative growth inhibition of bioactivity
#'     levels for the bacteria and virus strains listed below.  Each of the
#'     following elements in the \code{list} is a \code{data.frame} with 3 rows
#'     and 44 columns (with the exception of \emph{F. graminearum} which has 2
#'     rows).  The rows in each \code{data.frame} correspond to replications of
#'     the data collection process, while the columns correspond to relative
#'     growth inhibition bioactivity levels when subject to peptide libraries
#'     across fractions 1-43 and fraction 47.
#'
#' \describe{
#'
#'     \item{ec}{\emph{E. coli}}
#'
#'     \item{bc}{\emph{S. aureus}}
#'
#'     \item{pc}{\emph{K. pneumoniae}}
#'
#'     \item{oc}{\emph{A. baumannii}}
#'
#'     \item{ab}{\emph{A. baumannii}}
#'
#'     \item{pa}{\emph{P. aeruginosa}}
#'
#'     \item{fg}{\emph{F. graminearum}}
#'
#' }
"bioact"




#' Mass spectrometry data
#'
#' The mass spectrometry data collected for and described in Kirkpatrick et
#' al. (2016).  See paper for a full description of the data collection process,
#' or the package vignette for an abridged description.
#'
#' @docType data
#'
#' @usage data(mass_spec)
#'
#' @format A \code{data.frame} with 30,799 mass spectrometry levels and 38
#'     variables:
#'
#' \describe{
#'
#'     \item{m/z}{mass-to-charge ratio}
#'
#'     \item{Retention time (min)}{The time in minutes at which the peak
#'         retention time was achieved}
#'
#'     \item{Mass}{mass in daltons}
#'
#'     \item{Charge}{charge state}
#'
#'     \item{20150207_CLK_BAP_VO_11}{intensity of each MS feature in fraction 11}
#'
#'     \item{20150207_CLK_BAP_VO_12}{intensity of each MS feature in fraction 12}
#'
#'     \item{20150207_CLK_BAP_VO_13}{intensity of each MS feature in fraction 13}
#'
#'     \item{20150207_CLK_BAP_VO_14}{intensity of each MS feature in fraction 14}
#'
#'     \item{20150207_CLK_BAP_VO_15}{intensity of each MS feature in fraction 15}
#'
#'     \item{20150207_CLK_BAP_VO_16}{intensity of each MS feature in fraction 16}
#'
#'     \item{20150207_CLK_BAP_VO_17}{intensity of each MS feature in fraction 17}
#'
#'     \item{20150207_CLK_BAP_VO_18}{intensity of each MS feature in fraction 18}
#'
#'     \item{20150207_CLK_BAP_VO_19}{intensity of each MS feature in fraction 19}
#'
#'     \item{20150207_CLK_BAP_VO_20}{intensity of each MS feature in fraction 20}
#'
#'     \item{20150207_CLK_BAP_VO_21}{intensity of each MS feature in fraction 21}
#'
#'     \item{20150207_CLK_BAP_VO_22}{intensity of each MS feature in fraction 22}
#'
#'     \item{20150207_CLK_BAP_VO_23}{intensity of each MS feature in fraction 23}
#'
#'     \item{20150207_CLK_BAP_VO_24}{intensity of each MS feature in fraction 24}
#'
#'     \item{20150207_CLK_BAP_VO_25}{intensity of each MS feature in fraction 25}
#'
#'     \item{20150207_CLK_BAP_VO_26}{intensity of each MS feature in fraction 26}
#'
#'     \item{20150207_CLK_BAP_VO_27}{intensity of each MS feature in fraction 27}
#'
#'     \item{20150207_CLK_BAP_VO_28}{intensity of each MS feature in fraction 28}
#'
#'     \item{20150207_CLK_BAP_VO_29}{intensity of each MS feature in fraction 29}
#'
#'     \item{20150207_CLK_BAP_VO_30}{intensity of each MS feature in fraction 30}
#'
#'     \item{20150207_CLK_BAP_VO_31}{intensity of each MS feature in fraction 31}
#'
#'     \item{20150207_CLK_BAP_VO_32}{intensity of each MS feature in fraction 32}
#'
#'     \item{20150207_CLK_BAP_VO_33}{intensity of each MS feature in fraction 33}
#'
#'     \item{20150207_CLK_BAP_VO_34}{intensity of each MS feature in fraction 34}
#'
#'     \item{20150207_CLK_BAP_VO_35}{intensity of each MS feature in fraction 35}
#'
#'     \item{20150207_CLK_BAP_VO_36}{intensity of each MS feature in fraction 36}
#'
#'     \item{20150207_CLK_BAP_VO_37}{intensity of each MS feature in fraction 37}
#'
#'     \item{20150207_CLK_BAP_VO_38}{intensity of each MS feature in fraction 38}
#'
#'     \item{20150207_CLK_BAP_VO_39}{intensity of each MS feature in fraction 39}
#'
#'     \item{20150207_CLK_BAP_VO_40}{intensity of each MS feature in fraction 40}
#'
#'     \item{20150207_CLK_BAP_VO_41}{intensity of each MS feature in fraction 41}
#'
#'     \item{20150207_CLK_BAP_VO_42}{intensity of each MS feature in fraction 42}
#'
#'     \item{20150207_CLK_BAP_VO_43}{intensity of each MS feature in fraction 43}
#'
#'     \item{20150207_CLK_BAP_VO_47}{intensity of each MS feature in fraction 47}
#' }
"mass_spec"
