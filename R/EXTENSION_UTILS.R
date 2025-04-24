# MULTICHAIN / MUTLIMODEL UTILS ---- 

#' Extend Multi-Chain MSA
#'
#' Combines two MSAs across chains (or independently aligned segments) either by matching sample names or by row order.
#'
#' @param msa1 First alignment matrix or named list of matrices.
#' @param msa2 Second alignment matrix or named list of matrices.
#' @param use_sample_names Logical. If TRUE, attempts to match sample rows by name before combining.
#'
#' @return Combined alignment matrix or list of matrices, with each sample’s sequence extended across chains.
#' @export
extend_msa = function(msa1, msa2, use_sample_names = T){
  # if msa1 and msa2 are both lists of MSAs use names to connect them #
  if(is.list(msa1) && is.list(msa2)){
    # should be same order but use names instead #
    combined_msas = list()
    
    x = names(msa1)
    y = names(msa2)
    z = intersect(x, y)
    
    for(i in 1:length(z)){
      id = z[i]
      
      if(use_sample_names){
        x = msa1[[id]]
        y = msa2[[id]]
        
        # order y on x
        y = y[match(rownames(x), rownames(y)),]
      } # if not using sample names we assume order #
      
      # cbind and save #
      combined_msas[[i]] = cbind(x, y)
      names(combined_msas)[i] = id
    }
    
    # add back unique to x
    x = setdiff(names(msa1), z)
    for(i in 1:length(x)){
      id = x[i]
      combined_msas[[id]] = msa1[[id]]
      names(combined_msas)[length(combined_msas)] = id
    }
    
    # add back unique to y
    y = setdiff(names(msa2), z)
    for(i in 1:length(y)){
      id = y[i]
      combined_msas[[id]] = msa2[[id]]
      names(combined_msas)[length(combined_msas)] = id
    }
    
  } else {
    # just single msas #
    if(use_sample_names){
      msa1 = msa1[match(rownames(msa2), rownames(msa1)),]
    }
    
    combined_msas = cbind(msa1, msa2)
  }
  
  return(combined_msas)
}



#' Merge Nucleotide Windows Across Models
#'
#' Combines nucleotide patch windows from different structural models (e.g. chain A + chain B, or 2 structures).
#'
#' @param nuc_patch1 A data frame of nucleotide patches with \code{codon} and \code{nuc} columns.
#' @param nuc_patch2 A second patch data frame to merge.
#'
#' @return A merged data frame with codon identifiers and unioned nucleotide windows.
#' @export
extend_nuc_windows = function(nuc_patch1, nuc_patch2){
  # -- codon is same between both nuc_patches -- #
  codons = unique(c(nuc_patch1$codon, nuc_patch2$codon))
  
  # for now drop '' and '-' ~ maybe later '-' wont be here #
  codons = codons[codons != '' & codons != '-']
  
  merged_patch = data.frame(
    residue_id = NA,
    nuc = NA,
    codon = codons
  )
  
  for(i in 1:nrow(merged_patch)){
    w1 = nuc_patch1$nuc[nuc_patch1$codon == merged_patch$codon[i]]
    w2 = nuc_patch2$nuc[nuc_patch2$codon == merged_patch$codon[i]]
    
    # split by and join unique #
    w1 = unlist(strsplit(w1, '\\+'))
    w2 = unlist(strsplit(w2, '\\+'))
    
    w = unique(c(w1, w2))
    
    # arrange on number before :
    w = w[order(as.numeric(gsub('^([0-9]+):.*', '\\1', w)))]
    
    w = paste(w, collapse = '+')
    merged_patch$nuc[i] = w
  }
  
  # update residue id (its now meaningless given we are not tied to a structure)
  merged_patch$residue_id = paste0('codon_', merged_patch$codon)
  
  return(merged_patch)
  
}
# AB/EPI UTILS ----
#' Identify Antibody–Antigen Contacts
#'
#' Finds epitope and paratope residues from a PDB structure using chain-based atom distances.
#'
#' @param pdb A \code{bio3d} PDB object.
#' @param ag_chain Antigen chain ID.
#' @param h_chain Heavy chain ID.
#' @param l_chain Light chain ID.
#' @param dist_cutoff Maximum Å distance for contact (default = 5).
#'
#' @return A list with:
#' \item{epitope}{Concatenated residue IDs (e.g., "35_A_ 42_A_") of antigen residues near antibody.}
#' \item{paratope_h}{Heavy chain contact residues.}
#' \item{paratope_l}{Light chain contact residues.}
#' \item{contacts}{Matrix of all atom-level contacts.}
#' @export
identify_epitopes = function(pdb, ag_chain = NULL, h_chain = NULL, l_chain = NULL, dist_cutoff = 5){
  
  # remove H and HETATM
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')
  
  # set insert to '' if NA
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''
  
  # Grab chains of interest #
  ag_c = ag_chain
  h_c = h_chain
  l_c = l_chain
  
  pdb = bio3d::trim.pdb(pdb, chain = c(ag_c, h_c, l_c))
  
  pdb$atom$extra = paste(pdb$atom$resno, pdb$atom$chain, pdb$atom$insert, sep = '_')
  
  dist_mat = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom[,'extra'])
  
  colnames(dist_mat) = rownames(dist_mat) = unique(pdb$atom$extra)
  
  # grab close contacts
  close = which(dist_mat <= dist_cutoff, arr.ind = T)
  # replace pos with names
  close[, 1] = rownames(dist_mat)[close[, 1]]
  close[, 2] = colnames(dist_mat)[as.numeric(close[, 2])]
  
  # filter out heavy and light from column 1
  # filter out epi from column 2
  # ** probably break if l_c or h_c missing
  close = close[-grep(paste0('_', h_c, '_'), close[, 1]),]
  close = close[-grep(paste0('_', l_c, '_'), close[, 1]),]
  close = close[-grep(paste0('_', ag_c, '_'), close[, 2]),]
  rownames(close) = NULL
  
  # grab epi and paratope
  epi = paste0(unique(close[,1]), collapse = '+')
  para_h = paste0(close[grepl(paste0('_', h_c, '_'), close[, 2]), 2], collapse = '+')
  para_l = paste0(close[grepl(paste0('_', l_c, '_'), close[, 2]), 2], collapse = '+')
  
  return(list(
    epitope = epi, 
    paratope_h = para_h, 
    paratope_l = para_l,
    contacts = close)
  )
}