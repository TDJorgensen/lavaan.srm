.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###################################################################")
    packageStartupMessage("This is ", paste(pkgname, version), ", a package for fitting")
    packageStartupMessage("structural equation models to round-robin data.")
    packageStartupMessage("See the NEWS file for latest features, bug reports.")
    packageStartupMessage(pkgname, " was first developed at the University of Amsterdam")
    packageStartupMessage("partially funded by NWO grants 016.Veni.195.457 and 406.XS.01.078")
    packageStartupMessage("###################################################################")
}


