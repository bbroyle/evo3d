#' Run PatchR Workflow on a Single MSA–PDB Pair
#'
#' Wrapper for full pipeline: loads alignment and structure, identifies patches, maps codons, and calculates selection statistics.
#'
#' @param msa A character matrix or list representing the multiple sequence alignment.
#' @param pdb A \code{bio3d} PDB object.
#' @param chain Optional chain ID to restrict structure analysis.
#' @param msa_path Optional path to a FASTA file if \code{msa} is not provided directly.
#' @param pdb_path Optional path to a PDB file if \code{pdb} is not provided directly.
#' @param patch.dist.cutoff Maximum distance (Å) for residues to be considered part of a patch.
#' @param patch.rsa.cutoff Minimum relative solvent accessibility for a residue to serve as a patch center.
#'
#' @return A list containing:
#' \item{selection_df}{Residue-level data frame with nucleotide diversity, Tajima’s D, and haplotype diversity.}
#' \item{msa_info}{Reference and peptide info from \code{WRAPPER_msa_to_ref()}.}
#' \item{pdb_info}{Structure and patch data from \code{WRAPPER_pdb_to_patch()}.}
#' \item{nuc_info}{Codon mapping and patch-level MSAs from \code{WRAPPER_align_msa_pdb()}.}
#' @export
# RUN WHOLE PIPELINE # ----
run_patchr_single = function(msa, pdb, chain = NULL, msa_path = NULL, pdb_path = NULL,
                             patch.dist.cutoff = 15, patch.rsa.cutoff = 0.1){
  # read in MSA #
  msa_info = WRAPPER_msa_to_ref(msa, msa_path = msa_path)
  
  # read in pdb #
  pdb_info = WRAPPER_pdb_to_patch(pdb, 
                                  pdb_path = pdb_path,
                                  chain = chain,
                                  patch.dist.cutoff = patch.dist.cutoff, 
                                  patch.rsa.cutoff = patch.rsa.cutoff)
  
  # get nuc positions (have to say chain -- should be able to avoid?)
  nuc_info = WRAPPER_align_msa_pdb(msa_info = msa_info,
                                   pdb_info = pdb_info, 
                                   chain = chain)
  
  # calculate selection #
  selection_df = suppressWarnings(run_pegas_three(nuc_info$msa_subsets, pdb_info$residue_df))
  
  return(list(
    selection_df = selection_df,
    msa_info = msa_info,
    pdb_info = pdb_info,
    nuc_info = nuc_info
  ))
  
}