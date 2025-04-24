

#' Standardize MSA Input
#'
#' Ensures that the input multiple sequence alignment (MSA) is formatted as a matrix.
#'
#' @param msa_mat A matrix or list representing the MSA.
#' @param msa_path Optional file path to a FASTA file, used if \code{msa_mat} is not provided.
#'
#' @return A character matrix with rows as sequences and columns as alignment positions.
.standardize_msa_input = function(msa_mat, msa_path = NULL){
  # check if msa is a path or a matrix #
  if(!is.null(msa_path)){
    msa = bio3d::read.fasta(msa_path)
    msa_mat = msa$ali
  }
  
  # if msa object from bio3d::read.fasta() is given -- extract msa_mat
  if(inherits(msa_mat, 'fasta')){
    msa_mat = msa_mat$ali
  }
  
  return(msa_mat)
}



#' Get Reference Sequence from MSA
#'
#' Extracts a reference sequence from a multiple sequence alignment based on completeness, consensus, or a specified row.
#'
#' @param msa_mat Character matrix of sequences (rows) by positions (columns).
#' @param method Character or numeric. One of \code{"most_complete"}, \code{"consensus"}, or an integer indicating the row number.
#' @param msa_path Optional path to FASTA file (used if \code{msa_mat} not provided).
#'
#' @return A named character string representing the reference sequence.
#' @export
get_reference_sequence = function(msa_mat, method = 1, msa_path = NULL){
  # msa_mat is a matrix with rows as sequences and columns as positions #
  # method can be "most_complete", "consensus", or a numeric value (row number)
  
  # format msa_mat
  msa_mat = .standardize_msa_input(msa_mat, msa_path)
  
  # get reference sequence
  if (method == "most_complete") {
    
    # count gaps in each sequence or 'N' **3/25 -- this is problematic for proteins (N is valid) -- but DNA might have N **
    gaps = apply(msa_mat, 1, function(x) sum(!x %in% c('-', 'N', 'X')))
    
    # return sequence with least_gaps
    min_gap_idx = which.min(gaps)
    ref = paste0(msa_mat[min_gap_idx, ], collapse = '')
    names(ref) = paste0('ref.',rownames(msa_mat)[min_gap_idx])
    return(ref)
    
  } else if (method == "consensus") {
    
    # calculate most frequent non-gap character at each position  
    x = apply(msa_mat, 2, function(x){
      y = as.data.frame(table(x), stringsAsFactors = F)
      
      # if any row is the gap row just remove from here # **3/25 -- this is problematic for proteins (N is valid) -- but DNA might have N **
      ro = which(y$x %in% c('-', 'N', 'X'))
      if(length(ro) > 0){
        y = y[-ro,]
      }
      
      # sort y by freq and return the first row
      y = y[order(y$Freq, decreasing = T),]
      
      return(y$x[1])
    })
    
    # paste consensus sequence
    ref = paste(x, collapse = '')
    names(ref) = 'ref.consensus'
    return(ref)
    
  } else if (is.numeric(method)) {
    
    # check that row exists
    if (method > nrow(msa_mat)) stop("Position exceeds number of sequences in MSA")
    
    # return sequence by number
    ref = paste0(msa_mat[method,], collapse = '')
    names(ref) = paste0('ref.', rownames(msa_mat)[method])
    return(ref)
    
  } else {
    stop('method must be one of "most_complete", "consensus", or a numeric value (row number)')
  }
}



#' Translate DNA to Amino Acids
#'
#' Converts a DNA reference sequence to a peptide sequence using the standard codon table.
#'
#' @param ref A named character string representing a DNA sequence.
#'
#' @return A named character string representing the translated amino acid sequence.
#' @export
translate_dna = function(ref){
  # translate a DNA sequence into an amino acid sequence
  # assumes that the sequence is in the correct frame
  
  # remove gaps if present (might not be good) ('A-T' goes to X)
  #ref = gsub('-', '', ref)
  
  # translate (NNN and --- are treated the same 'X')
  pep = seqinr::translate(strsplit(ref, '')[[1]])
  
  # check for internal stops (can proceed)
  if(any(pep == '*')){
    message('Internal stop codon(s) found in reference sequence\ntry different ref seq or try different frame')
  }
  
  # return as character vector
  pep = paste0(pep, collapse = '')
  names(pep) = names(ref)
  
  return(pep)
} 



#' Detect Sequence Type
#'
#' Determines whether a given sequence is nucleotide or protein.
#'
#' @param ref A character string representing a sequence.
#'
#' @return A string: either \code{"nucleotide"} or \code{"protein"}.
#' @export
detect_sequence_type = function(ref){
  # look for nuc or gap characters (if anything else we assume protein) #
  # similar to MUSCLE approach for nuc vs protein detection #
  if(grepl('[^ACGTUN-]', toupper(ref))){
    return('protein')
  } else {
    return('nucleotide')
  }
}


#' Extract Reference and Peptide Sequence from MSA
#'
#' Wrapper function that standardizes input, extracts a reference sequence, detects or applies a sequence type, and translates to peptide if needed.
#'
#' @param msa_mat Character matrix of the MSA, or NULL if \code{msa_path} is used.
#' @param msa_path Optional file path to FASTA file.
#' @param ref_method Reference extraction method: \code{"most_complete"}, \code{"consensus"}, or numeric row index.
#' @param seqtype Optional sequence type: \code{"protein"}, \code{"nucleotide"}, or NULL for auto-detection.
#'
#' @return A list with \code{msa_mat}, \code{ref} (reference sequence), and \code{pep} (amino acid sequence).
#' @export
WRAPPER_msa_to_ref = function(msa_mat, msa_path = NULL, ref_method = 1, seqtype = NULL){
  
  # format msa_mat
  msa_mat = .standardize_msa_input(msa_mat, msa_path)
  
  ref = get_reference_sequence(msa_mat, ref_method)
  
  # if seqtype is given, use that to determine if we need to translate #
  # make sure seqtype is protein or nucleotide #
  if(!is.null(seqtype)){
    if(seqtype == 'protein'){
      pep = ref
    } else if (seqtype == 'nucleotide'){
      pep = translate_dna(ref)
    } else {
      stop('seqtype must be "protein", "nucleotide", or NULL for auto-detection')
    }
  } else {
    # auto-detect sequence type #
    seqtype = detect_sequence_type(ref)
    
    # translate if nucleotide #
    if(seqtype == 'nucleotide'){
      pep = translate_dna(ref)
    } else {
      pep = ref
    }
  }
  
  return(list(msa_mat = msa_mat, ref = ref, pep = pep))
}