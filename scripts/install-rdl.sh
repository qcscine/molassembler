set -euo pipefail

if [[ ! -d RingDecomposerLib ]]; then
  git clone --depth 1 https://github.com/rareylab/RingDecomposerLib.git
  echo "Pulled RingDecomposerLib. Stopping."
  exit 0
fi

if [[ $# -ne 1 ]]; then
  echo "You need to supply an installation prefix as a command line argument!"
  exit 1
fi

(cd RingDecomposerLib && git apply ../rdl.patch)

mkdir build-rdl

cd build-rdl
cmake -DCMAKE_INSTALL_PREFIX=$1 ../RingDecomposerLib
make -j 4 install
cd ..
