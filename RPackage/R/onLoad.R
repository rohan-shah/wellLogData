.onLoad <- function(libname, pkgname)
{
	library.dynam(package="wellLogData", chname="wellLogData", lib.loc = .libPaths())
}
