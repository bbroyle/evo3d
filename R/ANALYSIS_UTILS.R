# OUTPUT UTILS ----
#' Write Patch FASTAs to Disk
#'
#' Writes MSA subsets (by patch) to individual FASTA files for downstream analysis.
#'
#' @param msa_subset A named list of character matrices (one per patch).
#' @param residue_df Data frame with residue patch info (from \code{map_patches_to_nucleotides()}).
#' @param fasta_dir Output directory to write FASTA files.
#' @param write_fastas Logical. If FALSE, returns the FASTA sequences without writing files.
#'
#' @return NULL (files written to disk) or list of FASTA sequences.
#' @export
write_patch_fastas = function(msa_subset, residue_df, fasta_dir = 'temp_fastas', write_fastas = F){
  # nuc_patches is from map_patches_to_nucleotides()
  # fasta dir for saving files
  # save.file to save files or if F return list of fastas
  
  # create dir if it doesn't exist
  if(!dir.exists(fasta_dir)){
    dir.create(fasta_dir, recursive = T)
  }
  
  for(i in 1:length(fasta_list)){
    fa = fasta_list[[i]]
    fa = apply(fa, 1, paste, collapse = '')
    fname = paste0(fasta_dir, '/codon_', names(fasta_list)[[i]], '.fa')
    
    # fix fasta seqs for writing
    fa = paste0('>', names(fa), '\n', fa, '\n')
    cat(fa, sep = '', file = fname)
  }
  
  
}


#' Generate PyMOL Patch Commands
#'
#' Converts residue patch definitions into PyMOL selection and coloring commands.
#'
#' @param patches Data frame with \code{patch} and \code{residue_id} columns.
#'
#' @return A list of command strings: \code{patches}, \code{colors}, and \code{centroid} selection.
#' @export
write_patch_pymol = function(patches){
  # expecting a data.frame with patch column #
  # patches will be converted to pymol selection commands
  # patch numbered to order in patches data.frame
  
  # for each row generate pymol command for selecting patch
  cmd_list = c()
  for(i in 1:nrow(patches)){
    patch = patches$patch[i]
    patch = data.frame(id = unlist(strsplit(patch, '\\+')))
    patch$resno = gsub('_.*', '', patch$id)
    patch$chain = gsub('^[^_]+_', '', patch$id)
    patch$ins = gsub('^[^_]+_', '', patch$chain)
    patch$chain = gsub('_.*', '', patch$chain)
    
    # each resi set should be paired with chain and patch
    patch_cmd = c()
    for(ch in unique(patch$chain)){
      resi = paste0(patch$resno[patch$chain == ch], patch$ins[patch$chain == ch], collapse = '+')
      cmd = paste0('resi ', resi, ' and chain ', ch)
      patch_cmd[length(patch_cmd)+1] = cmd
    }
    
    # combine patch_cmd into final output 
    cmd = paste0('select patch_', i, ', ', paste0('(', patch_cmd, ')', collapse = ' or '))
    
    cmd_list[i] = cmd
  }
  
  # copy and paste this into pymol patches #
  #cat(cmd_list, sep = '\n')
  
  # get color per patch #
  colors <- grDevices::hcl.colors(nrow(patches), palette = 'Dynamic')
  colors = gsub('#', '', colors)
  colors = colors[sample(1:length(colors))]
  color_cmd = paste0('color 0x', colors, ', patch_', 1:nrow(patches))
  
  # copy and paste this into pymol for colors #
  #cat(color_cmd, sep = '\n')
  
  # last generate color black for all centroids #
  centroid = data.frame(id = patches$residue_id)
  centroid$resno = gsub('_.*', '', centroid$id)
  centroid$chain = gsub('^[^_]+_', '', centroid$id)
  centroid$ins = gsub('^[^_]+_', '', centroid$chain)
  centroid$chain = gsub('_.*', '', centroid$chain)
  
  # could exist on different chains #
  centroid_cmd = list()
  for(ch in unique(centroid$chain)){
    resi = paste0(centroid$resno[centroid$chain == ch], centroid$ins[centroid$chain == ch], collapse = '+')
    cmd = paste0('resi ', resi, ' and chain ', ch)
    centroid_cmd[length(centroid_cmd)+1] = cmd
  }
  
  # combine patch_cmd into final output 
  centroid_cmd = paste0('select centroid, ', paste0('(', centroid_cmd, ')', collapse = ' or '))
  centroid_cmd[2] = 'color black, centroid'
  
  # copy and paste this into pymol patches #
  #cat(centroid_cmd, sep = '\n')
  
  # return these cmds for pymol #
  return(list(patches = cmd_list, 
              colors = color_cmd, 
              centroid = centroid_cmd))
}


#' Write Selection Stats to PDB B-Factors
#'
#' Embeds residue-level statistics into the B-factor column of a PDB file for visualization.
#'
#' @param patch_df Data frame with residue-level statistics and \code{residue_id}.
#' @param pdb A \code{bio3d} PDB object.
#' @param stat_name Name of the column in \code{patch_df} to write to the B-factor.
#' @param outfile Output path for the modified PDB file.
#'
#' @return NULL. Writes a PDB file to disk with updated B-factors.
#' @export
write_stat_to_bfactor = function(patch_df, pdb, stat_name = 'tajima', outfile = 'test.pdb'){
  # expecting patches as dataframe #
  # expecting pdb as pdb object #
  
  # if stat is NA convert to 0
  patch_df[is.na(patch_df[stat_name]), stat_name] = 0
  
  # set pdb insert to '' if NA
  pdb$atom$ins[is.na(pdb$atom$ins)] = ''
  pdb$atom$residue_id = paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$ins, sep = '_')
  
  # use residue id (resno_chain_ins) to match to pdb #
  pdb$atom$b = patch_df[match(pdb$atom$residue_id, patch_df$residue_id), stat_name]
  pdb$atom$b = round(pdb$atom$b, 2)
  pdb$atom$b = ifelse(is.na(pdb$atom$b), 0, pdb$atom$b)
  
  # drop residue ID, and write
  pdb$atom = pdb$atom[,!names(pdb$atom) %in% c('residue_id')]
  bio3d::write.pdb(pdb = pdb, b = pdb$atom$b, file = outfile)
  
}


# SELECTION UTILS ----
#' Calculate Diversity and Neutrality Stats per Patch
#'
#' Runs nucleotide diversity, Tajima’s D, and haplotype diversity on a single or list of MSAs.
#'
#' @param msa A single alignment matrix or a list of patch-specific alignments.
#' @param residue_df Optional residue data frame to append results to.
#'
#' @return A data frame with one row per patch, including \code{pi}, \code{tajima}, and \code{hap}.
#' @export
run_pegas_three = function(msa, residue_df = NULL){
  # check if input is a list of msa or just single msa matrix #
  if(is.list(msa)){
    seqs = ape::as.DNAbin.list(msa)
    pi = lapply(seqs, pegas::nuc.div)
    taj = lapply(seqs, pegas::tajima.test)
    hap = lapply(seqs, pegas::hap.div) # gives warnings if gaps present (okay?)
    
    # if residue_df is provided we can add results to residue_df // else make a dataframe #
    if(is.null(residue_df)){
      residue_df = data.frame(residue_id = names(msa), stringsAsFactors = F)
    }
    
    residue_df$pi = NA
    residue_df$tajima = NA
    residue_df$tajima_pnormal = NA
    residue_df$tajima_pbeta = NA
    residue_df$hap = NA
    
    residue_df$pi[match(names(pi), residue_df$residue_id)] = unlist(pi)
    residue_df$tajima[match(names(taj), residue_df$residue_id)] = unlist(lapply(taj, function(x) x[1]))
    residue_df$tajima_pnormal[match(names(taj), residue_df$residue_id)] = unlist(lapply(taj, function(x) x[2]))
    residue_df$tajima_pbeta[match(names(taj), residue_df$residue_id)] = unlist(lapply(taj, function(x) x[3]))
    residue_df$hap[match(names(hap), residue_df$residue_id)] = unlist(hap)
    
    return(residue_df)
    
  } else {
    # individual msa #
    seqs = ape::as.DNAbin(msa)
    pi = pegas::nuc.div(seqs)
    taj = pegas::tajima.test(seqs)
    hap = pegas::hap.div(seqs) # gives warnings if gaps present (okay?)
    
    # just was single sequence -- return results #
    return(data.frame(
      residue_id = names(msa),
      pi = pi,
      tajima = taj$D,
      tajima_pnormal = taj$Pval.normal,
      tajima_pbeta = taj$Pval.beta,
      hap = hap
    ))
  }
  
}