# FastLORS
An R package for joint modeling to perform eQTL mapping

This R package allows users to perform joint modeling of SNPs and the expression of genes. The package provides an implementation of LORS from Yang et al. (2013) as well as Fast-LORS from Rhyne et al. (2018). Fast-LORS uses the proximal gradient method to solve the LORS problem and improves upon the computational time required for joint modeling.

The package also includes two screening methods for reducing the number of SNPs included in the modeling. LORS-Screening from Yang et al. (2013) and HC-Screening from Rhyne et al. (2018) are the two screening methods included in the package. If the data is small enough so that modeling is feasible without screening, the screening option can simply be set to "None."
