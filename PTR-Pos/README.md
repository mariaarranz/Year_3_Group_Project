# Code Function
This code includes two python files:
1. [**PTR-AI values**:](PTR.py)
  - This code will iterate through the file *'Table_EV4.xlsx' (from Eraslan, B, et.al)* and generate a new excel workbook called *'PTRdata.xlsx'* with the codon frequencies for all genes. 
  - Then, this will create a new excel sheet in the *'PTRdata.xlsx'* file and calculate the PTR-AI (protein to mRNA ratio adaption index) using the equation from the article *Quantification and discovery of sequence determinants of protein‐per‐mRNA amount in 29 human tissues (Eraslan, B, et.al)* - as shown below - and setting it to the log of the values.
  - Using these values, the code will then generate a box-plot of [**PTR-AI values vs Codons**](Figures/Box_plot.png)
2. [**PTR-AI values for specific codon vs. Position of specific Codon**:](CPosition.py)
  - This code will use the original data file with the sequences: *'Table_EV4.xlsx' (from Eraslan, B, et.al)*, and divide the sequences into codons, in order to generate a list with all the codons of a sequences in order (**starting from 0**).
  - The code will then loop through the codon list and find the secific codon (currently the selected codon is: ***AAA***), generating a new list containing the position of said codon on the list - returning '00' if the codon does not appear on the sequence.
  - Then the code takes in the PTR-AI values for said codon from the file generated from the preivou code: *'PTRdata.xlsx'* and generates a new excel file with two sheets: *'Pos.xlsx'* (sheet 1 = selected codon position; sheet 2 = selected codon PTR-AI values)
  - Using these lists, the code will generate a scatter plot of [**PTR-AI values of Codon vs Position of Codon**](Figures/Scatter_plot.png)
  
## Background Knowledge and Equations
From reading the article: *Quantification and discovery of sequence determinants of protein‐per‐mRNA amount in 29 human tissues (Eraslan, B, et.al)*

PTR = protein-ot-mRNA ratio = *log10(normalised protein levels/normalised mRNA levels)* - [Click here to see PTR values for liver](PTR_liver.xlsx)

PTR-AI = PTR ratio fold-change associated with doubling the frequency of a codon in a gene.
  - PTR-AI = 10^log2(frequency of codon)

## Understanding of Results
As expected, the results from the [**PTR-AI values vs Codons Box PLot**](Figures/Box_plot.png) will show a similar analysis as that of the article: *Quantification and discovery of sequence determinants of protein‐per‐mRNA amount in 29 human tissues (Eraslan, B, et.al)*. 

This states that: 
> PTR-AI of individual codons showed consistent amplitudes and directions across tissues which contests the hypothesis of widespread tissue specific post transcriptional regulation due to a varying tRNA pool among different tissues
This means that there was a consistent spread accross all tisues, meaning that codons are important in the effect of translation.

Furthermore, when looking at the [**PTR-AI values of Codon vs Position of Codon Scatter Plot**](Figures/Scatter_plot.png), we can conclude that there are higher PTR-AI values that correlate to the codon at earlier positions, meaning that codons located earlier in the sequence have more of an effect than those later in the sequence.

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
Biopython
```
  pip install biopython
```
Matplotlib (for figures)
```
python -m pip install matplotlib
```

# Warnings
I suggest running the [**PTR.py**:](PTR.py) code first in order to generate the *'PTRdata.xlsx'* file and be able to use it in the second code, [**CPosition.py**:](CPosition.py).

The codes take a bit to generate the excel workbooks, however once they are generated once, you may comment the code to only produce the figures. However, the code as a whole does not take more than 3 min to completely run through so it is not necessary to comment out the code.

# References
Eraslan, B., Wang, D., Gusic, M., Prokisch, H., Hallström, B. M., Uhlen, M., Asplund, A., Pontén, F., Wieland, T., Hopf, T., Hahne, H., Kuster, B. & Gagneur, J. (2019) Quantification and discovery of sequence determinants of protein-per-mRNA amount in 29 human tissues. Molecular Systems Biology. 15 (2), Available from: https://www.embopress.org/doi/abs/10.15252/msb.20188513.

# Author 
Maria Arranz Fombellida
