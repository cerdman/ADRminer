.onAttach=function(libname, pkgname){ 
  packageStartupMessage("Loaded ADRminer ", as.character(packageDescription("ADRminer")[["Version"]]),"\n")
}