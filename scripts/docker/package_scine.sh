git clone git@gitlab.chab.ethz.ch:scine/scine.git
cd scine
git fetch --tags
git checkout v0.2.0
rm -rf .git
cd ..
tar -cf scine.tar scine
rm -rf scine

cp scine.tar clang6/
mv scine.tar gcc7/
