# Mass spectrometry data visualization

This folder contains a script used for visualization of mass spectrometry–based proteomics data generated from RLRb immunoprecipitation experiments in Nematostella vectensis.

MS/MS raw files and MaxQuant output files were deposited to the ProteomeXchange Consortium via the PRIDE partner repository under the dataset identifier **PXD064383**.

The visualization script uses the processed input file ***MaxQuant-Perseus-stats-Rubi-20210201.xlsx***, which was generated using MaxQuant (v1.5.3.12) followed by statistical analysis in Perseus (n = 3 technical replicates per RLRb-IP). Protein identification and quantification were performed using label-free quantification (LFQ) with standard filtering criteria, including a 1% false discovery rate (FDR) at the peptide and protein levels and a minimum of two unique or razor peptides per protein.

This code is intended for downstream data exploration and figure generation and does not reproduce the full MaxQuant or Perseus analysis workflow.
