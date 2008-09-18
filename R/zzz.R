.onLoad <- function(libname,pkgname){
  require("tools",quietly=TRUE)
  if(interactive()) pkg2WinMenu(pkgname)
}

.onAttach <- function(libname,pkgname){
  library.dynam(pkgname,pkgname,lib.loc=libname)
  options(warn = -1)
  cat("\n")
  cat("---------------------------------------------------------------------------------\n")
  pkg.info <- drop(read.dcf(file=file.path(.find.package(package=pkgname,lib.loc=libname),"DESCRIPTION")))
  cat(paste("Package",sQuote(pkg.info["Package"]),"versione",pkg.info["Version"],"loaded and ready to use.",sep=" "),"\n\n")
  cat(pkg.info["Description"],"\n")
  cat(paste("\n",sQuote(pkg.info["Package"])," by ",pkg.info["Author"],".",sep=""),"\n")
  cat("=================================================================================\n")
  #if(packageHasNamespace(package=pkgname,package.lib=libname))
  #cat("NAMESPACE:","si","\n")
  #else
  #cat("NAMESPACE:","no","\n")
  #if(any(names(getLoadedDLLs())==pkgname))
  #cat("C/Fortran:","si","\n")
  #else
  #cat("C/Fortran:","no","\n")
  #cat("=================================================================================\n")
  filename <- file.path(.find.package(package=pkgname,lib.loc=libname),"doc")
  result <- list_files_with_type(filename,"vignette")
  if(length(result)==0){
  cat(paste(sQuote(pkg.info["Package"])," not contains Vignettes.",sep=""))
  cat("\n")
  }
  if(length(result)>0){
  cat(paste(sQuote(pkg.info["Package"])," contains Vignettes.",sep=""))
  cat("\nType",sQuote(paste("vignette(package = \"",pkg.info["Package"],"\")",sep="")),
  "to list the Vignettes.")
  cat("\nChoose",sQuote(pkgname),"from menu bar to recall a Vignette.\n")
  }
  cat("=================================================================================\n")
  cat("Type",sQuote(paste("help(package = \"",pkg.info["Package"],"\")",sep="")),
  "or",sQuote(paste("package?",pkg.info["Package"],sep="")),"to have an overview.\n")
  #cat("Type",sQuote(paste("data(package = \"",pkg.info["Package"],"\")",sep="")),
  #"to view a list of data frames.\n")
  cat("Type",sQuote(paste("ls(\"package:",pkg.info["Package"],"\")",sep="")),
  "to view a list of functions.\n")
  #cat("Digitare",sQuote(paste("citation(package=\"",pkg.info["Package"],"\")",sep="")),
  #"per vedere la bibliografia.\n")
  cat("Type",sQuote(paste("detach(package:",pkg.info["Package"],")",sep="")),
  "to remove it from the search() path.\n")
  cat("---------------------------------------------------------------------------------\n")
  cat("\n")
  options(warn=0)
  return(invisible(0))  
}

.Last.lib <- function(libname){
  stringa <- unlist(strsplit(x=libname,split="/"))
  pkgname <- stringa[length(stringa)]
  library.dynam.unload(pkgname,libpath=libname)
  txt <- paste("Thank you to use",sQuote(pkgname),"package. See you.")
  writeLines(txt)
}
