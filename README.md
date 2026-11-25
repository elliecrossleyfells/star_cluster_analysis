# Star Cluster Analysis
Observational astronomy project using the PIRATE robotic telescope to investigate stellar populations, evolution, ages and distances of galactic star clusters (M39, M52 and M15) through aperture photometry, colour–magnitude diagrams, main-sequence fitting and Python-based FITS data analysis.

## Photometric Analysis of Galactic Star Clusters Images
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24534135.svg)](https://doi.org/10.5281/zenodo.24534135)


### Overview
An undergraduate astrophysics project involving acquiring observations from the PIRATE telescope to study the stellar populations, evolution, relative ages and distances of two open clusters (M39 and M52) and one globular cluster (M15).

### Methodology
- Acquired astronomical observations using B- and V-band filters at multiple exposure times.
- Extracted FITS files from ZIP archives using Python and performed aperture photometry using AstroImageJ.
- Processed relative stellar flux measurements using Python to calculate apparent magnitudes and B−V colour indices.
- Constructed colour–magnitude diagrams (CMDs) to investigate stellar populations and relative cluster ages.
- Estimated open cluster distances using main-sequence fitting, accounting for interstellar reddening and extinction, and propagated magnitude uncertainties into distance estimates.
  
### Findings
The CMDs showed main-sequence trends for M39 and M52 and a distinct giant branch for M15. The inferred order of increasing age was M52, M39 and M15. Limitations included potential field-star contamination, photometric uncertainties and elevated sky background.

### Repository contents
- notebooks/file-extraction-code/file_extraction.py – Batch extraction of FITS observations from ZIP archives.
- notebooks/star_cluster_analysis.ipynb – Photometric data processing, CMD generation, distance estimation and uncertainty analysis.
- requirements.txt – Python dependencies.

### Data and reproducibility
The observational FITS dataset is archived on Zenodo: [(https://doi.org/10.5281/zenodo.24534135)]. Aperture photometry was performed separately in AstroImageJ to produce the csv files located in aperture-photometry-data; the notebook uses exported photometric measurements rather than raw FITS images.
