#'Similar to veganifyOTU from phyloseq.
#'@keywords internal
#'@export
reverseASV= function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(physeq)
}
