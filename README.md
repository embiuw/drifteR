# drifterR
This package can be used to detect drift dives in dive data obtained using Satellite Relay Data Loggers (SRDLs) manufactured by the Sea Mammal Research Unit.

The package includes the function \code{rbs} for determining the order in which the broken stick algorithm has extracted the most salient points in the time depth profile, by applying a 'Reverse Broken Stick' algorithm. This is based on methods described in Photopoulou et al. (2015).

Using output from the \code{rbs} function, the \code{pDrift} function then calculates the putative drift rate, along with a series of weighting variables for the probability of a dive being a true drfift dive.

The \code{fitDrift} function then fits a state-space model to the weighted data to obtain an estimate of the drift rate change trajectory through time.

# Installation
devtools::install_github('enbiuw/drifteR')

# References
Biuw, M., McConnell, B., Bradshaw, C. J. A., Burton, H., & Fedak, M. (2003). Blubber and buoyancy: monitoring the body condition of free-ranging seals using simple dive characteristics. Journal of Experimental Biology, 206(19), 3405–3423. https://doi.org/10.1242/jeb.00583

Photopoulou, T., Lovell, P., Fedak, M. A., Thomas, L., & Matthiopoulos, J. (2015). Efficient abstracting of dive profiles using a broken-stick model. Methods in Ecology and Evolution, 6(3), 278–288. https://doi.org/10.1111/2041-210X.12328
