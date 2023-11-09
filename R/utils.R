#' Run the pacman 00LOCK file destroyer
#'
#' @description Stop 00LOCK files getting in the way of package builds
#'
#' @export
unlock <- function() {
  for(.lib in .libPaths()) {
    ff <- grep("00LOCK", list.files(.lib), value=TRUE)
    print(ff)
    if(length(ff)>0) unlink(paste0(.lib,"/",ff), recursive=TRUE, force=TRUE)
  }
  #pacman::p_unlock()
  invisible()
}
