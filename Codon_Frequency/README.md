# Code Function
Extract coding sequence from gene ID and create excel file containing codon table, RSCU and CAI values.

# Requirements
Pandas
```
  pip install pandas
```
Requests
```
  pip install requests
```
Biopython
```
  pip install biopython
```
Openpyxl 
```
  pip install openpyxl
```
Cai
```
  pip install CAI
```

# Warnings
loops through whole table so only use once on college computer!
Takes time to run! (Run once overnight)

# RSCU
Codons with values > 1 = abundant
Codons with values < 1 = less-abundant

# CAI
higher CAI = higher gene expression

# References
- CAI code and package: (https://github.com/Benjamin-Lee/CodonAdaptationIndex)
- Biopython package: (https://biopython.org/wiki/Download)
- Openpyxl: (https://openpyxl.readthedocs.io/en/stable/)
- Reference sequence used: **_ecol.heg.fasta_** (from https://github.com/Benjamin-Lee/CodonAdaptationIndex/tree/master/example_seqs)

# Authors 
Marianne Buat, Ruben Weitzman, Maria Arranz
