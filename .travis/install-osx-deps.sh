# Install required software
brew update
brew cask uninstall --force oclint
brew install git
brew upgrade cmake
brew tap homebrew/science
#brew unlink gcc
brew install gcc@7
#brew install open-mpi
#brew install netcdf

# Make sure the weird gfortran library links are in place.
#ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
#ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

