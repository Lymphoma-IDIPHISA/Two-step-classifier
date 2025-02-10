# Two-step DLBCL genetic classification
Simplified genetic classification based on those proposed by Schmitz-Wright and Lacy and their research teams. Alterations of specific genes were used to classify the samples as **N1<sup>2-S</sup>**, **BN2<sup>2-S</sup>**, **EZB<sup>2-S</sup>**, **MCD<sup>2-S</sup>**, or **ST2<sup>2-S</sup>**. The A53 LymphGen subtype was not included due to the lack of copy number variation data in our series and most published series. The genes from the Schmitz-Wright and Lacy studies were selected by calculating their power to classify determined by the Fisher's exact test, as featured in these studies. Several tests were then carried out to select the best combination of genes, taking into account the sensitivity and specificity in the PdH and HMRN cohorts (Lacy et. al., 2020) according to the LymphGen and Lacy classification (AIC cluster). In the first step, at least one of the top genes should be mutated: NOTCH1, for N1<sup>2-S</sup> (NOTCH1 mutated samples were classified as N1<sup>2-S</sup>, regardless the presence of other alterations); MYD88, CD79B, and PIM1 for MCD<sup>2-S</sup>; BCL6 translocation, NOTCH2, BCL10, and TNFAIP3 mutations for BN2<sup>2-S</sup>; BCL2, EZH2, and CREBBP for EZB<sup>2-S</sup>; and SGK1, TET2, and SOCS1 for ST2<sup>2-S</sup>. Samples with the same score for two or more subtypes, or samples with no mutations, were classified in the second step, in which we added the following genes for each subtype: PRMD1, BTG1, PIM2, and CD58 for MCD<sup>2-S</sup>; UBE2A, CD70, CCND3, and DTX1 for BN2<sup>2-S</sup>; TNFRSF14, KMT2D, IRF8, and EP300 for EZB<sup>2-S</sup>, and STAT3 for ST2<sup>2-S</sup>. In this second step, at least two genes from any of the genes that define each subtype should be mutated to assign the sample to the corresponding subtype.

For NOTCH1 and NOTCH2 we only consider truncation mutations that affect the C-terminal PEST domain.

## Usage
This script is made for running in R with a mutational table as the example.txt file. 

The columns and rows represent genes and samples, respectively.
1 indicates mutated or translocated gene in the sample and 0 indicates non-mutated or non-translocated gene.

## References
- Pedrosa, L. et al. Proposal and validation of a method to classify genetic subtypes of diffuse large B cell lymphoma. Sci. Rep. 11, 1886 (2021)
- Schmitz, R. et al. Genetics and Pathogenesis of Diffuse Large B-Cell Lymphoma. N. Engl. J. Med. 378, 1396–1407 (2018).
- Wright, G. W. et al. A Probabilistic Classification Tool for Genetic Subtypes of Diffuse Large B Cell Lymphoma with Therapeutic Implications. *Cancer Cell* 37, 551-568.e14 (2020).
- Chapuy, B. et al. Molecular subtypes of diffuse large B cell lymphoma are associated with distinct pathogenic mechanisms and outcomes. *Nat. Med*. 24, 679–690 (2018).
- Lacy, S. E. et al. Targeted sequencing in DLBCL, molecular subtypes, and outcomes: a Haematological Malignancy Research Network report. *Blood* 135, 1759–1771 (2020).
