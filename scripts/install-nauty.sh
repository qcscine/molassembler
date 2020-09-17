set -euo pipefail

if [[ ! -f nauty27r1.tar.gz ]]; then
  wget http://pallini.di.uniroma1.it/nauty27r1.tar.gz
  echo "Pulled nauty. Stopping."
  exit 0
fi

if [[ $# -ne 1 ]]; then
  echo "You need to supply an installation prefix as a command line argument!"
  exit 1
fi

tar -xf nauty27r1.tar.gz
mv nauty27r1 nauty
cp nauty.cmake nauty/CMakeLists.txt
mkdir build-nauty

cd build-nauty
cmake -DCMAKE_INSTALL_PREFIX=$1 ../nauty
make -j 4 install
cd ..
