# QC plots
These QC plots pertain to the ribosome profiling data collected in Wu et al., and also Hsu et al..

## Directory contents
Read length distribution plot: *RPF_longest_CDS_read_dist.pdf*
Periodicity by region: *RPF_period_region.pdf*
Periodicity heatmap: *RPF_period_region_length.pdf*

## How to interpret the plots

### Read length distribution

Each histogram corresponds to one sample, and shows the distribution of RPF read lengths. 

A sharp peak in the 28-30 nt region generally indicates canonical ribosome protected fragments. 

The secondary peak around 20-21 nt corresponds to ribosomes in an alternative conformational ribosome state. (Reference: Lareau, L.F., Hite, D.H., Hogan, G.J. and Brown, P.O., 2014. Distinct stages of the translation elongation cycle revealed by sequencing ribosome-protected mRNA fragments. *elife*, 3, p.e01257. [https://pmc.ncbi.nlm.nih.gov/articles/PMC4052883/](https://pmc.ncbi.nlm.nih.gov/articles/PMC4052883))

### Periodicity by region

These bar p[lots show the percentage of P-sites in each reading frame, split by transcript region (5'UTR, CDS, 3'UTR), aggregrated across all read lengths. It is generally expected that frame 0 should be the dominant signal in the CDS, but no in the 5'UTR or 3'UTR.

### Periodicity heatmap

Similar to the "Periodicity by region" plots, this heatmap shows percentage of P-sites in the reading frame. It is expected that the majority of these should be in-frame at the CDS and generally centered around the expected size for ribosome footprints as well as potentially at the 20-21 size (as it corresponds to an alternative ribosome conformation).
