# Code Function
This code includes 2 python files and one excel workbook:

**GOscatter.py:**  
- reads in excel file **GOtACIdata.xlsx**, and generates 

**GOtACI_all_tissues.py:** 
- Code to get tACI values for all 4 tissues (colon, brain, bladder, prostate) with tRNA data.
- reads in **Table_EV1.xlsx** to get gene information (from Wang, et.al).
- reads in **tACI_values.xlsx** to get tACI values for all tissues.
- generates excel file: **alltissuestRNAdisplay3.xlsx** with all tissues displayed

## Background Knowledge and Equations
Gene ontology (GO) annotations are curated from gene product experiments and provide a statement linking a gene product to a molecular function, biological process or cellular component type. Furthermore, GO slims are a subset of GO terms that give a broad, species specific overview of the ontology of large datasets1. To compare tACI values between process-related gene sets and to explore the relation between the ontologies of tissue enriched gene sets, biological process GO slim terms were assigned to the genes enriched in each tissue for which tRNA data was available. Enriched gene sets were derived as previously described and Generic GO Term Mapper was used to map each gene to one of 71 biological process homosapien GOslim terms.2 For each term, the mean tACI for all genes assigned to it was calculated and values for the each tissue enriched gene set compared.

# Requirements
Pandas
```
  pip install pandas
```
Numpy
```
  pip install numpy
```
Openpyxl (for working with excel workbooks)
```
  pip install openpyxl
```
Matplotlib (for figures)
```
  python -m pip install matplotlib
```
# References
Wang, D., Eraslan, B., Wieland, T., Hallström, B., Hopf, T., Zolg, D. P., Zecha, J., Asplund, A., Li, L., Meng, C., Frejno, M., Schmidt, T., Schnatbaum, K., Wilhelm, M., Ponten, F., Uhlen, M., Gagneur, J., Hahne, H. & Kuster, B. (2019) A deep proteome and transcriptome abundance atlas of 29 healthy human tissues. Molecular Systems Biology. 15 (2), Available from: https://www.embopress.org/doi/abs/10.15252/msb.20188503. Available from: doi: 10.15252/msb.20188503. 

# Authors 
Cerys Barclay
