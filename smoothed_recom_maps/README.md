## Smoothed recombination maps

Recombination data were derived from four independent F2 QTL crosses as outlined in Whiting et al 2021 (https://www.biorxiv.org/content/10.1101/2021.03.18.435980v1.full)


To extract information on recombination rate, maps were smoothed using the smooth.spline function in R.
Smoothing was necessary to estimate recombination over regions where individual markers were not syntenic with the published reference genome (micro-rearrangements).
Larger rearrangements (such as assembly errors and inversions), where multiple markers were out of order, were manually re-oriented resulting in smoothed maps.
Smoothed maps have an average of 242 markers per chromosome.

Note that due to the paritcularly low recombination rate on LG12, only two points are present on this chromosome. 
