# FYP: My Final-Year Project

## What is This?

This is a repo which will store code I have done/am doing for my final-year project.  The aim of the project as a whole is to investigate
B --> X mu<sup>+</sup>mu<sup>-</sup> decays at the LHCb, which seem to be (from [previous searches](http://lhcb-public.web.cern.ch/lhcb-public/Welcome.html#P5p))
a potential avenue towards new physics.  So the idea is to conduct an inclusive search (previous forays have primarily been exclusive),
which provides a better environment for the theorists to make predictions.  However, on the experimental/data-analysis side, inclusive
searches are a lot more difficult than exclusive ones.  The main issue is the background, which mostly stems from cascade decays e.g. B0 --> D (K mu nu) mu, which can masquerade as a signal, because in the inclusive search all we look for are the two muons.

Currently we have identified that there are two distinct regions:

* High Q<sup>2</sup>, in which backgrounds are small and stem from combinatorial sources, and we can apply control data sets and extrapolation to estimate the backgrounds (we hope)
* Low Q<sup>2</sup>, in which background comes from the more-thorny cascade decays.

## Current Goals

* Comparison of high- and low-Q<sup>2</sup> regions.
* Generation of plots overlaying simulated signal and background in both regions for analysis and comparison.
* Creation of a Boosted Decision Tree or Neural Network to conduct a Multivariate Analysis to separate the signal from background.
* Train MVA on Monte Carlo-simulated signal and background
* Apply MVA to real data.

## Useful Commands

* In ROOT:
    - To Open TMVAGui: `.x /cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_69/ROOT/5.34.21/x86_64.../tmva/test/TMVAGui.C(<MVA output file name>)` 
