# Code Function

There are 2 main parts to this folder:
### A) Calculate equivalent tRNA abundance of each anticodon
v1: Manually Excluded initiator Met and Selenocysteine tRNAs manually from tRNA Excel files.
v2: Manually include initiator Met and SeC tRNAs (need to differentiate from normal other tRNAs). Further statistical analysis and plotting carried out.

### B) Calculate tACI of genes' coding sequence from

tAI_obj.py: define a tAI object to calculate the tACI of gene sequences given a tRNA pool. 
calculatetACIallgenes.py: calculate tACI for all the gene sequences contained in Table EV4 from Eraslan et al (2019). Repeat calculation of tACI for gene CDS with the reading frame shifted by 1 or 2 nucleotides for altering codon usage. Carry out hypothesis testing on tACi values obtained to test for bias in method.

# Requirements
Python 3 or above.
Python libraries:
numpy, pandas, scipy, matplotlib

Files:
#### A) tRNA abundance
- Normal_Tissues_Only.xlsx: file containing tRNA probe measurements and tRNA gene correponding to each probe for normal tissue samples from the brain, bladder, prostate, B cells and colon.

#### B) tACI
- Table_EV4_excel.xlsx : supplementary data from Eraslan et al (2019), saved to an Excel format for easier manipulation. Contains the gene name, Ensemble IDs, genomic and protein sequences for 11,575 genes.

# Warnings

Beware broken links for file locations when importing data into Python as all Excel files have been moved to a single folder "Excel_Files".
Some intermediate files used during computation may not have been included in this Git for concision, however they should be recreated easily by running the original generating code.

# References

Dana, A. & Tuller, T. (2014) The effect of tRNA levels on decoding times of mRNA codons. Nucleic Acids Research. 42 (14), 9171-9181. Available from: https://www.ncbi.nlm.nih.gov/pubmed/25056313. Available from: doi: 10.1093/nar/gku646. 

dos Reis, M., Savva, R. & Wernisch, L. (2004) Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Research. 32 (17), 5036-5044. Available from: https://doi.org/10.1093/nar/gkh834. Available from: doi: 10.1093/nar/gkh834. [Accessed 11/16/2019]. 

Gingold, H., Tehler, D., Christoffersen, N., Nielsen, M., Asmar, F., Kooistra, S., Christophersen, N., Christensen, L. L., Borre, M., Sørensen, K., Andersen, L., Andersen, C., Hulleman, E., Wurdinger, T., Ralfkiær, E., Helin, K., Grønbæk, K., Ørntoft, T., Waszak, S., Dahan, O., Pedersen, J., Lund, A. & Pilpel, Y. (2014) A Dual Program for Translation Regulation in Cellular Proliferation and Differentiation. Available from: http://www.sciencedirect.com/science/article/pii/S0092867414010423. 

Eraslan, B., Wang, D., Gusic, M., Prokisch, H., Hallström, B. M., Uhlen, M., Asplund, A., Pontén, F., Wieland, T., Hopf, T., Hahne, H., Kuster, B. & Gagneur, J. (2019) Quantification and discovery of sequence determinants of protein-per-mRNA amount in 29 human tissues. Molecular Systems Biology. 15 (2), Available from: https://www.embopress.org/doi/abs/10.15252/msb.20188513.

Sabi, R. & Tuller, T. (2014) Modelling the Efficiency of Codon–tRNA Interactions Based on Codon Usage Bias. DNA Research. 21 (5), 511-526. Available from: https://doi.org/10.1093/dnares/dsu017. Available from: doi: 10.1093/dnares/dsu017. 

Python script tAI_obj.py (tAI calculation object) adapted from Shyam Saladi (2016) _tAI_ [code], available from: https://github.com/smsaladi/tAI [accessed 24/03/2020]
Copyright (C) 2016 Shyam Saladi - Caltech

# Authors 
Cerys Barclay, Marianne Buat
