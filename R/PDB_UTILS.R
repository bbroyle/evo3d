# DOWNLOAD PDB UTILS ----
#' Find Matching PDB Structures
#'
#' Uses BLAST to find and optionally download PDB structures matching a query peptide sequence.
#'
#' @param pep A protein sequence (character string).
#' @param identity_cutoff Minimum identity percentage (default 80).
#' @param max_hits Maximum number of structures to return.
#' @param generate_plot Logical, whether to return a ggplot summary.
#' @param download_pdbs Logical, whether to download matching PDBs.
#' @param output_dir Directory to save downloaded structures.
#'
#' @return A list with BLAST results, plots, and optionally downloaded structure paths.
#' @export
find_matching_structures = function(pep, identity_cutoff = 80, max_hits = 5, generate_plot = T,
                                    download_pdbs = F, output_dir = 'retrieved_pdbs'){
  # If no PDB is provided to program -- this function will blast RCSB pdb for matching structures #
  # We want to return multiple structures for more robust window and solvent accessibility calculations #
  # essentially a wrapper for bio3d::blast.pdb() & bio3d::get.pdb() with helpful plot and table output

  message("\n⚠️ This function makes live queries to external databases ⚠️\n",
          "Please avoid running in parallel or sending requests too quickly\n",
          "You may be rate-limited or blocked\n")

  # blast pep (expecting single protein sequence) against PDB database
  bl = bio3d::blast.pdb(pep)
  hits = bl$hit.tbl

  # -- if identity_cutoff <= 1, we assume it is a percentage and convert to 0-100 -- #
  # -- no one wants 1% identity hits -- #
  if(identity_cutoff <= 1){
    identity_cutoff = identity_cutoff * 100
  }

  # sort by mlog.evalue and filter by identity_cutoff
  hits = hits[order(hits$mlog.evalue, decreasing = T),]
  filtered_hits = hits[hits$identity >= identity_cutoff,] # can be empty -- proceed anyway

  # get top hits or all if nrow(filtered_hits) is less than max_hits
  top_hits = filtered_hits[seq_len(min(max_hits, nrow(filtered_hits))),]

  # add status to hits
  hits$pdb_status = 'raw_hit'
  hits$pdb_status[hits$subjectids %in% filtered_hits$subjectids] = 'filtered_hit'
  hits$pdb_status[hits$subjectids %in% top_hits$subjectids] = 'top_hit'

  # generate plot if desired
  if(generate_plot){
    hit_plot = ggplot2::ggplot(hits, ggplot2::aes(alignmentlength, identity, color = pdb_status))+
      ggplot2::geom_jitter(height = 0.1, width = 0.5)+
      ggrepel::geom_text_repel(ggplot2::aes(label = subjectids), show.legend = F)+
      ggplot2::geom_hline(yintercept = identity_cutoff, linetype = 'dashed')+
      ggplot2::theme_bw()+
      ggplot2::scale_color_manual(values = c('raw_hit' = '#d95f02',
                                             'filtered_hit' = '#7570b3',
                                             'top_hit' = '#1b9e77'))+
      ggplot2::labs(title = 'PDB Blast Results',
                    subtitle = paste0('max_hits = ',
                                      max_hits, ', identity_cutoff = ',
                                      identity_cutoff, '%'),
                    x = 'Reported Alignment Length',
                    y = 'Identity (%)')
  } else {
    hit_plot = NULL
  }

  # ---------------------------------------------------------------------------- #
  # check if any hits passed filters -- if not stop here and return data at hand #
  if(nrow(top_hits) == 0){
    message('\n🛑 No hits found >= identity cutoff: ', identity_cutoff, '% 🛑',
            '\nInspect `all_hits` and `hit_plot` for more information\n',
            'If you want to continue manually provide `all_hits` to:\n',
            'download_structures(output$all_hits)\n')

    return(list(
      all_hits = hits,
      hit_plot = hit_plot)
    )
  }

  # ---------------------------------------------------------------------------- #
  # if download is desired, call download_structures() #
  # return hit_tables and hit_plot #
  if(download_pdbs){
    hit_table = suppressMessages(download_structures(top_hits, output_dir))

    return(list(
      top_hits = top_hits,
      all_hits = hits,
      hit_plot = hit_plot,
      pdb_table = hit_table)
    )
  } else {
    # return top_hits, all_hits, and hit_plot #
    return(list(
      top_hits = top_hits,
      all_hits = hits,
      hit_plot = hit_plot)
    )
  }

}


#' Download PDB Structures
#'
#' Downloads structures from a table of BLAST hits or PDB IDs.
#'
#' @param hit_table Data frame of BLAST hits or character vector of PDB IDs.
#' @param output_dir Output directory to save downloaded PDB files.
#'
#' @return A data frame with PDB IDs, chain IDs, and local file paths.
#' @export
download_structures = function(hit_table, output_dir = 'retrieved_pdbs'){
  # hit_table is one of the two table outputs of find_matching_structures() #
  # hit_list$all_hits or hit_list$top_hits #
  # or a user provided table with at least subjectids column #
  # or can be a character vector of subjectids #

  # subjectids are four letter pdb codes -- anything longer is trimmed in bio3d::get.pdb() #
  # user may provide unique subjectids that are duplicated in the first 4 characters #
  # thats fine I will add downloaded file path to both #

  # !!! output_dir cannot be empty !!! # -- need to check #
  if(output_dir == ''){
    output_dir = '.'
  } else if (grepl('/$', output_dir)){
    output_dir = substr(output_dir, 1, nchar(output_dir)-1)
  }

  message("\n⚠️ This function makes live queries to external databases ⚠️\n",
          "Please avoid running in parallel or sending requests too quickly\n",
          "You may be rate-limited or blocked\n")

  # check if hit_table is a character vector #
  if(is.character(hit_table)){
    hit_table = data.frame(subjectids = hit_table)
  }

  # check if subjectids column is available #
  if(!'subjectids' %in% colnames(hit_table)){
    message('\n🛑 hit_table must have a subjectids column or be a character vector 🛑\n')
    return(NULL)
  }

  # make output directory (also handled by bio3d::get.pdb())
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)
  }

  # output is variable depending on content of output directory and success of download #
  # if all successful (or previously downloaded) output is unnamed vector of file paths #
  # if some failed "1" and newly download get "0", while previous download get path #
  file_paths = suppressWarnings(bio3d::get.pdb(hit_table$subjectids, path = output_dir))

  # if all download then unnamed vector -- if some missing then named vector output #
  if(is.null(names(file_paths))){
    names(file_paths) = file_paths
  }

  # use names as paths -- drop paths that failed
  file_paths = file_paths[file_paths != '1']

  # add file_paths to hit_table (force ordering)
  ins = toupper(substr(hit_table$subjectids, 1, 4))
  colord = paste0(output_dir, '/', ins, '.pdb')

  hit_table$file_paths = names(file_paths[colord])

  # return a table of subject ids and file paths
  results = data.frame(
    pdb_id = gsub('_.+', '', hit_table$subjectids),
    chain_id = gsub('^[^_]+|_', '', hit_table$subjectids),
    file_path = hit_table$file_paths,
    stringsAsFactors = F
  )

  return(results)

}

# PDB UTILS ----
#' Standardize PDB Input
#'
#' Reads a PDB or CIF file and returns a trimmed structure for specified chains.
#'
#' @param pdb A PDB object (optional).
#' @param chain Character vector of chain IDs, or \code{"all"}.
#' @param pdb_path Path to a PDB or CIF file.
#'
#' @return A trimmed \code{bio3d} PDB object containing only selected chains.
#' @export
.standardize_pdb_input = function(pdb, chain = NULL, pdb_path = NULL){
  # if pdb_file is provided read in pdb
  if(!is.null(pdb_path)){
    # if file ends with .cif use read.cif()
    if(grepl('cif$', pdb_path)){
      pdb = bio3d::read.cif(pdb_path)
    } else {
      pdb = bio3d::read.pdb(pdb_path)
    }
  }

  # if chain is NULL - default to first chain / can also us "all"
  if(is.null(chain)){
    chain = unique(pdb$atom$chain)[1]
  } else if('all' %in% chain){
    chain = unique(pdb$atom$chain)
  }

  # trim pdb to chains #
  pdb = bio3d::trim.pdb(pdb, chain = chain)

  # could remove hydrogens and hetatm (but maybe hetatm will be useful) #

  # pass along pdb #
  return(pdb)

}



#' Extract Sequence from PDB
#'
#' Retrieves the amino acid sequence from a PDB file or object.
#'
#' @param pdb A PDB object from \code{bio3d}.
#' @param chain Optional chain ID(s).
#' @param pdb_path Path to a PDB file.
#'
#' @return A named character vector of sequences, one per chain.
#' @export
get_pdb_sequence = function(pdb, chain = NULL, pdb_path = NULL){
  pdb = .standardize_pdb_input(pdb, chain = chain, pdb_path = pdb_path)

  # -- THIS IS HANDLED BY .standardize_pdb_input() -- #
  # if is NULL or all # -- should should find a better option here #
  # if chain is NULL - default to first chain / can also us "all"
  #if(is.null(chain) | 'all' %in% chain){
  #  # chain options
  #  chain = unique(pdb$atom$chain)
  #
  #  # if chain is null default to first chain
  #  if(is.null(chain)){
  #    chain = unique(pdb$atom$chain)[1]
  #  }
  #}

  # get pdb sequences with names as chains
  aa_seq = sapply(chain, function(x){
    paste0(
      bio3d::pdbseq(bio3d::trim.pdb(pdb, chain = x)),
      collapse = '')
  }, USE.NAMES = T)

  return(aa_seq)
}



#' Residue-wise Distance Matrix
#'
#' Computes a pairwise distance matrix between residues based on 3D coordinates.
#'
#' @param pdb PDB object or NULL.
#' @param chain Chain ID(s).
#' @param pdb_path Path to a PDB file.
#' @param distance_method One of \code{'all'}, \code{'ca'}, \code{'backbone'}, or \code{'sidechain'}.
#'
#' @return A square numeric matrix of distances.
#' @export
calculate_residue_distance = function(pdb, chain = NULL, pdb_path = NULL, distance_method = 'all'){

  # distance_method can be 'backbone', 'sidechain', 'ca', 'all'
  # hydrogens always removed

  # format pdb
  pdb = .standardize_pdb_input(pdb, chain = chain, pdb_path = pdb_path)

  # pdb should be protein only and no H (what about glycan distance?)
  pdb = bio3d::trim.pdb(pdb, 'protein')
  pdb = bio3d::trim.pdb(pdb, 'noh')

  # METHOD IS TOO MUCH #
  # check method, maybe we are only doing backbone or sidechain
  if(distance_method == 'backbone'){
    pdb = bio3d::trim.pdb(pdb, 'backbone')

  } else if(distance_method == 'sidechain'){
    pdb = bio3d::trim.pdb(pdb, 'sidechain')

  } else if(distance_method == 'ca'){
    pdb = bio3d::trim.pdb(pdb, 'calpha')
  }

  # create residue / atom identifiers
  pdb$atom$insert[is.na(pdb$atom$insert)] = ''
  pdb$atom$residue_id = paste0(pdb$atom$resno, '_', pdb$atom$chain, '_', pdb$atom$insert)

  # get atom distance matrix
  res_dist = bio3d::dm.xyz(pdb$xyz, grpby = pdb$atom$residue_id)

  # set row and column names
  colnames(res_dist) = rownames(res_dist) = unique(pdb$atom$residue_id)

  # set diag to 0 and make symmetrical
  res_dist <- pmin(res_dist, t(res_dist), na.rm = TRUE)
  diag(res_dist) <- 0

  # return residue wise distance matrix
  return(res_dist)

}



#' Calculate Solvent Accessibility
#'
#' Estimates residue-wise solvent accessibility using a DSSP rewrite.
#'
#' @param pdb A bio3d PDB object.
#' @param chain Chain ID(s).
#' @param pdb_path Path to the PDB file.
#' @param method One of \code{'rose'}, \code{'miller'}, \code{'theoretical_tien'}, or \code{'empirical_tien'}.
#' @param drop_incomplete Logical, drop residues missing backbone atoms.
#'
#' @return A data frame with residue indices, exposure values, and metadata.
#' @export
calculate_accessibility = function(pdb, chain = NULL, pdb_path = NULL,
                                   method = 'rose', drop_incomplete = F){
  # return residue wise solvent accessibility

  # format pdb
  pdb = .standardize_pdb_input(pdb, chain = chain, pdb_path = pdb_path)

  # check method is valid #
  if(!method %in% c('rose', 'miller', 'theoretical_tien', 'empirical_tien')){
    message('method not recognized - setting to default: rose')
  }


  # MAX ACC VALUE TABLE ------------------------------------------- #
  max_acc = data.frame(
    residue = c('A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'),
    theoretical_tien = c(129.0, 274.0, 195.0, 193.0, 167.0, 223.0, 225.0, 104.0, 224.0, 197.0,
                         201.0, 236.0, 224.0, 240.0, 159.0, 155.0, 172.0, 285.0, 263.0, 174.0),
    empirical_tien = c(121.0, 265.0, 187.0, 187.0, 148.0, 214.0, 214.0, 97.0, 216.0, 195.0,
                       191.0, 230.0, 203.0, 228.0, 154.0, 143.0, 163.0, 264.0, 255.0, 165.0),
    miller = c(113.0, 241.0, 158.0, 151.0, 140.0, 183.0, 189.0, 85.0, 194.0, 182.0,
               180.0, 211.0, 204.0, 218.0, 143.0, 122.0, 146.0, 259.0, 229.0, 160.0),
    rose = c(118.1, 256.0, 165.5, 158.7, 146.1, 186.2, 193.2, 88.1, 202.5, 181.0,
             193.1, 225.8, 203.4, 222.8, 146.8, 129.8, 152.5, 266.3, 236.8, 164.5),
    stringsAsFactors = FALSE
  )

  # PREP PDB ------------------------------------------------------ #
  # Extract atom data
  atoms = pdb$atom

  # Drop hydrogen atoms and non-ATOM records
  atoms = atoms[atoms$elesy != 'H' & atoms$type == 'ATOM', ]

  # Handle residue numbering (each aa needs unique)
  atoms$orig_resno = atoms$resno
  atoms$resno = 1
  atoms$insert[is.na(atoms$insert)] = ''
  atoms$chain = paste0(atoms$chain, '_', atoms$insert)

  # renumber sequentially (keep atoms from same residue together)
  prev = atoms[1,]
  for(i in 2:nrow(atoms)){
    if(atoms$orig_resno[i] == prev$orig_resno &&
       atoms$chain[i] == prev$chain &&
       atoms$insert[i] == prev$insert){
      atoms$resno[i] = prev$resno
    } else {
      next_res = prev$resno + 1
      atoms$resno[i] = next_res
      prev = atoms[i,]
    }
  }

  # drop incomplete residues (same as original DSSP)
  if(drop_incomplete){
    # function to check completeness
    complete_check = function(resno){
      types = atoms$elety[atoms$resno == resno]
      all(c('N', 'CA', 'C', 'O') %in% types)
    }

    # check for complete residues
    resno_list = unique(atoms$resno)
    incomplete = sapply(resno_list, function(x) !complete_check(x))

    if(any(incomplete)){
      atoms = atoms[!atoms$resno %in% resno_list[incomplete],]
    }
  }

  # SETTING DATASETS FOR DSSP - rewrite  --------------------------- #
  resno_list = unique(atoms$resno)
  residue_df = data.frame(
    residue_index = resno_list,
    chain_id = sapply(resno_list, function(i) atoms$chain[atoms$resno == i][1]),
    residue_name = sapply(resno_list, function(i) bio3d::aa321(atoms$resid[atoms$resno == i][1])),
    orig_resno = sapply(resno_list, function(i) atoms$orig_resno[atoms$resno == i][1]),
    stringsAsFactors = FALSE
  )

  atom_df = data.frame(
    residue_index = atoms$resno,
    atom_type = atoms$elety,
    x = atoms$x,
    y = atoms$y,
    z = atoms$z,
    stringsAsFactors = FALSE
  )

  # RUNNING DSSP REWRITE ------------------------------------------- #
  #Rcpp::sourceCpp('src/dssp_sasa_rewrite.cpp')  # should move outside
  accessibility = calculateDSSPAccessibility(atom_df, residue_df)

  # Add results to the residue dataframe (round to 2 decimals) ** doesnt save mem **
  residue_df$sasa = round(accessibility, 2)

  # add rsa
  residue_df$max_acc = max_acc[match(residue_df$residue_name, max_acc$residue), method]
  residue_df$rsa = residue_df$sasa / residue_df$max_acc
  residue_df$rsa[residue_df$rsa > 1] = 1 # if rsa above 1 (set at 1)
  residue_df$rsa = round(residue_df$rsa, 2)

  # return with residue id
  residue_df$residue_id = paste0(residue_df$orig_resno, '_', residue_df$chain_id)

  # clean up a little bit #
  residue_df$aa = residue_df$residue_name
  residue_df$orig_chain = gsub('_.*', '', residue_df$chain_id)
  residue_df$orig_insert = gsub('.+_', '', residue_df$chain_id)

  # simplify output !! might want to have orig_resno, chain, insert seperate #
  residue_df = residue_df[, c('residue_index', 'aa', 'orig_resno', 'orig_chain', 'orig_insert', 'sasa', 'rsa', 'residue_id')]

  return(residue_df)
}

#' Identify Surface Patches
#'
#' Defines residue patches around surface-exposed centroids based on RSA, SASA, and spatial distance.
#'
#' @param dist_mat A residue-residue distance matrix.
#' @param accessibility_df Data frame of residue accessibility values (e.g., from \code{calculate_accessibility}).
#' @param rsa_cutoff Minimum RSA for centroid selection (default = 0.1).
#' @param sasa_cutoff Optional minimum SASA for centroid selection.
#' @param dist_cutoff Maximum distance (in Å) for neighbors.
#' @param max_patch Optional maximum number of neighbors per patch.
#' @param only_exposed_in_patch If TRUE, restricts patch members to surface residues only.
#'
#' @return Updated \code{accessibility_df} with a \code{patch} column listing neighbors.
#' @export
identify_patches = function(dist_mat, accessibility_df,
                            rsa_cutoff = 0.1, sasa_cutoff = NULL,
                            dist_cutoff = 15, max_patch = NULL,
                            only_exposed_in_patch = F){

  # returns patch around each residue #
  # dist_mat is residue_dist from calculate_distance_matrix
  # number of ways to define a window #
  # patch centroids defined by
  #     rsa_cutoff, or sasa_cutoff, or both NULL all residue
  # neighbors defined by
  #    dist_cutoff, or min_patch, or max_patch, or NULL all residues
  #    further filtered by only_exposed_in_patch

  # !! THERE ARE OTHER METHODS TO ADD -- COULD SEPERATE THESE PATCHERS
  # INTO DIFFERENT FUNCTIONS. OR Supply centroid_method = rsa, sasa, none
  # acc_cutoff. Window method = dist, min_patch, max_patch, none #


  # expect dist_mat and accessibility_df to have same amount of residues
  # could be different if residues dropped in calculate_accessibility()
  # dont really need this res_id checking ***
  res_ids = intersect(rownames(dist_mat), accessibility_df$residue_id)

  # get centroids
  if(!is.null(rsa_cutoff)){
    centroids = accessibility_df[accessibility_df$rsa >= rsa_cutoff &
                                   accessibility_df$residue_id %in% res_ids,]
  } else if(!is.null(sasa_cutoff)){
    centroids = accessibility_df[accessibility_df$sasa >= sasa_cutoff &
                                   accessibility_df$residue_id %in% res_ids,]
  } else {
    centroids = accessibility_df[accessibility_df$residue_id %in% res_ids,]
  }

  # shrink dist mat to only exposed residues if desired
  if(only_exposed_in_patch){
    in_set = centroids$residue_id
    dist_mat = dist_mat[in_set, in_set]
  }

  # get neighbors (in dist_cut and also check max_patch)
  centroids$patch = NA

  for(i in 1:nrow(centroids)){
    center = centroids$residue_id[i]

    # get neighbors
    neighbors = dist_mat[center,
                         names(which(dist_mat[center,] <= dist_cutoff))]

    # check max_patch
    if(!is.null(max_patch) && length(neighbors) > max_patch){
      # sort neighbors and take the first max_patch
      neighbors = sort(neighbors)
      neighbors = neighbors[1:max_patch]
    }

    # add to centroids
    centroids$patch[i] = paste0(names(neighbors), collapse = '+')
  }

  # add back to accessibility_df
  accessibility_df$patch = NA

  # join back
  accessibility_df$patch[match(centroids$residue_id, accessibility_df$residue_id)] = centroids$patch

  #accessibility_df$patch_method = paste0('-rsa ', rsa_cutoff, ' -sasa ', sasa_cutoff,
  #                                       ' -dist ', dist_cutoff, ' -max ', max_patch,
  #                                       ' -only_exposed ', only_exposed_in_patch)
  #

  return(accessibility_df)
}

#' Cluster Overlapping Patches
#'
#' Groups overlapping residue patches and selects one representative per cluster.
#'
#' @param patch_df A data frame with \code{patch} strings and \code{residue_id}.
#' @param overlap_cutoff Numeric threshold (0–1) for clustering based on neighbor overlap.
#'
#' @return A list with a residue data frame and the overlap matrix.
#' @export
cluster_patches = function(patch_df, overlap_cutoff = 0.5){
  # reduce number of windows by clustering on shared aa #
  # the most surface exposed centroid will be representative for cluster

  # remove positions with no patch
  patch_df = patch_df[!is.na(patch_df$patch),]
  neighbors = unique(unlist(strsplit(patch_df$patch, '\\+')))

  # set up clust_mat (one hot encode possible neigbors)
  clust_mat = matrix(0, nrow = nrow(patch_df), ncol = length(neighbors))
  rownames(clust_mat) = patch_df$residue_id
  colnames(clust_mat) = neighbors

  # fill with presence abscence (a better way might be fill with distance)
  for(i in 1:nrow(patch_df)){
    n = unlist(strsplit(patch_df$patch[i], '\\+'))
    clust_mat[i,n] = 1
  }

  # use my own dist() function #
  # 1 - shared / total ~ not symmetrical
  # 1 - shared / max(total1, total2) ~ symmetrical ** this one for now **
  # 1 - shared / min(total1, total2) ~ symmetrical

  overlap = matrix(0, nrow = nrow(clust_mat), ncol = nrow(clust_mat))
  rownames(overlap) = rownames(clust_mat)
  colnames(overlap) = rownames(clust_mat)

  for(i in 1:nrow(clust_mat)){
    for(j in 1:nrow(clust_mat)){
      # get shared neighbors
      shared = sum(clust_mat[i,] * clust_mat[j,])

      # get total neighbors
      #total = sum(clust_mat[i,])

      # min / max totals
      total = min(sum(clust_mat[i,]), sum(clust_mat[j,]))
      #total = max(sum(clust_mat[i,]), sum(clust_mat[j,]))

      # calculate distance
      overlap[i,j] = (shared / total)
    }
  }

  dist_mat = as.dist(1 - overlap)

  # clustering (many methods to chose ward.D is working fine)
  hc = hclust(dist_mat, method = "ward.D")

  # split clusters at overlap_cutoff
  clusters = cutree(hc, h = overlap_cutoff)

  patch_df$cluster = clusters[as.character(patch_df$residue_id)]

  # select representative for each cluster -- will use rsa
  unique_clusters = unique(clusters)
  patch_df$is_representative = FALSE

  # probably not the best -- update how it is selected
  # could try largest out of cluster distance (seems good) #
  for(c in unique(clusters)){
    ro = which(patch_df$cluster == c)
    # what is out of cluster distance for each of these rows #
    # overlap is % overlap (so lets minimize)
    out_overlap = overlap[ro, -ro]
    rmeans = rowMeans(out_overlap)
    rep = names(which.min(rmeans))

    patch_df$is_representative[patch_df$residue_id == rep] = TRUE
  }

  return(list(
    residue_df = patch_df,
    overlap_mat = overlap)
  )
}


#' Extract Surface Patches from PDB
#'
#' High-level wrapper for computing sequences, distances, solvent accessibility, and patches from a PDB structure.
#'
#' @param pdb A \code{bio3d} PDB object.
#' @param chain Chain ID(s) or \code{"all"}.
#' @param pdb_path Path to a PDB or CIF file.
#' @param distance_method One of \code{"all"}, \code{"ca"}, \code{"backbone"}, \code{"sidechain"}.
#' @param drop_incomplete_residue Logical. Drop residues with incomplete backbone atoms.
#' @param rsa.method Method for RSA normalization (e.g., \code{"rose"}).
#' @param patch.dist.cutoff Max Å distance between residues in a patch.
#' @param patch.rsa.cutoff RSA cutoff for defining patch centroids.
#' @param patch.sasa.cutoff Optional SASA cutoff for patch centroids.
#' @param patch.only.exposed Logical, restrict patches to exposed residues.
#' @param max.patch Optional max number of neighbors per patch.
#'
#' @return A list with \code{pdb}, \code{seq_set}, \code{residue_dist}, and \code{residue_df} (with patches).
#' @export
WRAPPER_pdb_to_patch = function(pdb, chain = NULL, pdb_path = NULL,
                                distance_method = 'all',
                                drop_incomplete_residue = F, rsa.method = 'rose',
                                patch.dist.cutoff = 15, patch.rsa.cutoff = 0.1,
                                patch.sasa.cutoff = NULL, patch.only.exposed = T,
                                max.patch = NULL){

  # step 0: validate pdb
  #print('Step 0: Validating PDB input')
  pdb = .standardize_pdb_input(pdb, chain = chain, pdb_path = pdb_path)

  # step 1: retrieve sequences for chains of interest
  #print('Step 1: Retrieving sequences')
  seq_set = get_pdb_sequence(pdb, chain = chain)

  # step 2: calculate residue-wise distance matrix
  #print('Step 2: Calculating residue-wise distance matrix')
  residue_dist = calculate_residue_distance(pdb, chain = chain,
                                            distance_method = distance_method)

  # step 3: calculate residue-wsie accessibility
  #print('Step 3: Calculating residue-wise accessibility')
  residue_df = calculate_accessibility(pdb, chain = chain,
                                       drop_incomplete = drop_incomplete_residue,
                                       method = rsa.method)

  # step 4: identify surface patches (expands residue_df)
  #print('Step 4: Identifying surface patches')
  residue_df = identify_patches(residue_dist,
                                residue_df, only_exposed_in_patch = patch.only.exposed,
                                dist_cutoff = patch.dist.cutoff,
                                rsa_cutoff = patch.rsa.cutoff,
                                sasa_cutoff = patch.sasa.cutoff)

  # return list object
  return(list(
    pdb = pdb,
    chain = chain,
    seq_set = seq_set,
    residue_dist = residue_dist,
    residue_df = residue_df
  ))

}
