# Analysis of a collection of samples from The Gambia obtained between 2014 and 2017 (gambiarcodes)

The scripts used here were used to analyze data and produce figures for the corresponding article.

The details of this pipeline are described in the methods section. More information are available as comments in the scripts directly.

The main script `gambiarcodes` is the main script controlling each step of the pipeline which can run independently from each other. However, it is worth mentioning that outputs produced at a step might be used as input for other steps.

The step of the pipeline leading to the generation of a network of pairwise IBD values includes the connection to a distant server which will run this rather long task thus avoiding having to dedicate a substantial part of the resources of the user device. This behavior can be changed to run the step locally by a few modifications.

A few scripts use files that are only accessible to the original owner of this repository. This scripts are not needed by the pipeline; their purpose was to make consistent raw data tables from different source files.

# Contact

Marc-Antoine Guery, PhD student at Antoine Claessens lab, LPHI, University of Montpellier.