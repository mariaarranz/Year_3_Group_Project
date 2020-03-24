# Code Function
This code is composed of 4 different Codes:

**BetaCorrelation-CODE.py:** 
- Uses only the log2(frequencies) of the codons as the 'xi' of the beta coefficient in relation to PTR of Brain Tissue.

**PearsonCorrelation-CODE.py:** 
- Generates the **Pearson's Correlation Coefficient** relating the PTR values of Brain Tissue and the Codon Frequencies
- Generates a Heat Map that demonstrates these correlations
- Creates a Scatter Plot that demonstrates the correlation coefficient of all codons as points (*importance is shown below*).

**CodonPosition-CODE.py:**
- Generates a scatter plot relating the position of the first appearance of codon **AAG** vs the PTR values for each gene sequence of brain tissue.

**GroupedBoxplots-CODE.py:**
- Finds the top and bottom 5% PTR values of Brain Tissue Genes
- Plots the Codon Frequencies of these genes with the top and bottom 5% PTR values as multiple boxplots.

## Background Knowledge and Equations
From reading the article: *Quantification and discovery of sequence determinants of protein‐per‐mRNA amount in 29 human tissues (Eraslan, B, et.al)*

**PTR** = protein-ot-mRNA ratio
  - PTR = *log10(normalised protein levels/normalised mRNA levels)*

**Multivaritive linear model:** *Yi = beta0 + beta*Xi + e*
  - Yi = tissue specific PTR ration (Table_EV3)
  - beta0 = intercept varied by tissue
  - beta = coefficient of sequence features (equal accross all tissues)
  - Xi = matrix of sequence predictors
  - e = errors

**Beta correlation**:
  - PTR-AI = *10^beta*
    - beta = *[(Xi - Xbar)(Yi - Ybar)] / (Xi - Xbar)^2
        
**Pearson's Correlation Coefficient:** can be used to summarize the strength of the linear relationship between two data samples.
  - coefficient returns a value between -1 and 1 that represents the limits of correlation from a full negative correlation to a full positive correlation.
        - If value = 0, there is no correlation.
        - If -0.5 > value > 0.5, there is a notable correlation.
        - If it is below these values, there is a less notable correlation.
  - Pearson's Corr Coeff = *[Covariance(X,Y)] / [std(X)*std(Y)]*

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
Seaborn
```
  pip install seaborn
```

# References
Eraslan, B., Wang, D., Gusic, M., Prokisch, H., Hallström, B. M., Uhlen, M., Asplund, A., Pontén, F., Wieland, T., Hopf, T., Hahne, H., Kuster, B. & Gagneur, J. (2019) Quantification and discovery of sequence determinants of protein-per-mRNA amount in 29 human tissues. Molecular Systems Biology. 15 (2), Available from: https://www.embopress.org/doi/abs/10.15252/msb.20188513.
How to Calculate Correlation Between Variables in Python. (n.d.). Retrieved from https://machinelearningmastery.com/how-to-use-correlation-to-understand-the-relationship-between-variables/ 

Machine Learning glossary of terms. (n.d.). Retrieved from https://deepai.org/machine-learning-glossary-and-terms/feature-extraction 

The genetic code & codon table (article) | Khan Academy. (n.d.). Retrieved from https://www.khanacademy.org/science/biology/gene-expression-central-dogma/central-dogma-transcription/a/the-genetic-code-discovery-and-properties 

# Author 
Maria Arranz Fombellida
