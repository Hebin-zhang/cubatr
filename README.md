## <p style="color:orange;">Attention: This Project Is No Longer Maintained</p>

I am pleased to announce that my professor has developed a more comprehensive and powerful R package, which covers the same functionality as this project. Therefore, I have decided to discontinue updates and maintenance for this project.

I strongly encourage you to switch to using my professor's R package for improved functionality and support. You can find their package at the following link: [cubar](https://github.com/mt1022/cubar)

Thank you for your past support and usage of this project! If you have any questions or need further assistance, you can still contact me, and I will do my best to provide support.

Once again, thank you!


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
