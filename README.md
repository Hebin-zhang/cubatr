# cubatr <img src='docs/mstile-150x150.png' align="right" height="139" />

Codon Usage Bias Analysis Toolkit (CUBAT) is a powerful package for codon usage bias (CUB) analysis.
And cubatr is its R version.
### Installation
You can install `cubatr`  from from GitHub: 
```{r}
devtools::install_github("Hebin-zhang/cubatr")
```

### Main features
cubatr is still under development, more features will be updated later. 
Here is a summary of indices supported by CUBAT and cubatr in comparison with other popular software.

|  Index                                 | CUBAT |cubatr| codonW | DAMBE | EMBOSS |
| -------------------------------------  | ----- |-----| ------ | ----- | ------ |
| RSCU (relative synonymous codon usage) | √     | √     | √      | √     |        |
| Nc (effective number of codons)        |       |       | √      |       | √      |
| Nc (effective number of codons,SYX13)  | √     |       |        | √     |        |
| CAI (codon adaptation index)           | √     | √     | √      | √     | √      |
| CAI2 (Xuhua Xia,2007)                  | √     | √     |        | √     |        |
| CBI (codon bias index)                 | √     |       | √      | √     |        |
| Fop (frequency of optimal codons)      | √     | √     | √      | √     |        |
| TAI (tRNA adaptation index)            | √     |       |        |       |        |
| CSC (codon stabilization coefficient)  | √     |       |        |       |        |
| Scaled χ<sup>2<sup>                    | √     |       |        |       |        |
| **Other Features**                     |
| Amino acid usage                       | √     | √     | √      | √     |        |
| Custom genetic code                    | √     | √     |        | √     |        |
| Cross-platform                         | √     | √     |        | √     | √      |
