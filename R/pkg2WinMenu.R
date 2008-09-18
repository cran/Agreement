  pkg2WinMenu <- function(pkgname){
  winMenuAdd(pkgname)
  vigFile <- system.file("Meta", "vignette.rds", package = pkgname)
  if (file.exists(vigFile)){
  vigMtrx <- .readRDS(vigFile)
  vigs <- file.path(.find.package(pkgname), "doc", vigMtrx[,"PDF"])
  names(vigs) <- vigMtrx[, "Title"]
  mykey <- unique(unlist(vigMtrx[,"Keywords"]))
  nkey <- length(mykey)
  nkey 
  mylist <- list()
  for(i in 1:nkey)
  mylist[[i]] <- vigs[vigMtrx[,"Keywords"] == mykey[i]]
  for(j in 1:nkey){
  pkgMenu <- paste(pkgname, mykey[j], sep = "/")
  winMenuAdd(pkgMenu)
  count <- 0
  for (i in mylist[[j]]) {
  count <- count + 1
  partI <- names(mylist[[j]][count])
  partII <- sub(".pdf","",basename(i))
  item <- paste(partI,"   ","(",partII,".pdf)",sep="")
  winMenuAddItem(pkgMenu, item, paste("shell.exec(\"", 
  as.character(i), "\")", sep = ""))
  }
  }
  }
  winMenuAddItem(pkgname, "-", "")
  mypath <- file.path(.libPaths(),pkgname,"Manual",paste(pkgname,".pdf",sep=""))
  winMenuAddItem(pkgname,"Manual",paste("shell.exec(\"",as.character(mypath), "\")", sep = ""))
  }



