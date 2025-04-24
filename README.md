# evopatchr

**evopatchr** is an R package for structure-aware population genetics, enabling patch-level evolutionary analysis of protein surfaces. It integrates codon-level selection metrics with 3D structural data to identify spatially clustered signals of diversity, selection, and epitope targeting.

Designed with applications in immunology, virology, and comparative genomics, `evopatchr` supports multi-chain proteins, multi-model structures, and antibody–antigen complexes.

---

## Key Features

- Computes selection metrics (π, haplotype diversity, Tajima’s D) at surface-defined patches
- Detects and integrates antibody epitopes from PDB structures
- Supports multi-chain and multi-model PDBs for robust structural inference
- Maps structure-defined residues to codon-aligned MSA windows
- Generates B-factor-encoded PDBs and per-residue selection tables
- Fully self-contained: no reliance on DSSP, MAFFT, or external alignment tools

---

## Installation

To install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yourusername/evopatchr")

# install.packages("BiocManager")
BiocManager::install("msa")
```

## Quick Example

```r
library(evopatchr)

msa_path <- system.file("extdata", "example.fasta", package = "evopatchr")
pdb_path <- system.file("extdata", "example.pdb", package = "evopatchr")

result <- run_patchr_single(msa_path = msa_path,
                             pdb_path = pdb_path,
                             chain = 'A')

write_stat_to_bfactor(result$selection_df,
                      result$pdb_info$pdb,
                      stat_name = "tajima",
                      outfile = "example_tajima.pdb")
```

## Package Structure

R/ – R functions for patch detection, alignment, epitope analysis, and visualization

src/ – C++ code for solvent accessibility (DSSP-style)

inst/extdata/ – Example FASTA and PDB files for reproducible workflows

## License

This package is released under the MIT License. 
Portions of the structural accessibility code are adapted from open-source DSSP
 implementations

## Conatct

Brad Broyles
PhD Candidate, Computational Biology
Purdue University
bbroyle@purdue.edu