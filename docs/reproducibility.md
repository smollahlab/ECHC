# Reproducibility Guide

## System
MacOS, Windows, Linux

## Major Packages
gen3DNet
plsdof, MASS, NbClust, cli, progress

## R
Depends R (>= 3.1.0)

## NMF
NMF > 0.23.0

## gen3DNet with specific parameters
  nmf_nrun = 100, (number of NMF iterations)
  p_val_threshold = 0.01, (p-value threshold)
  k_picker = max_ward_kl (Ward algorithm)
  
## Random Seeds
seed = nndsvd (a robust initialization strategy using the seeding algorithm, that is based on a non-negative double singular value decomposition)

## Runtime

Approx. 5 min
