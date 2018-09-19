#!/bin/bash
# adapted from https://stackoverflow.com/questions/151677/tool-for-adding-license-headers-to-source-files

# First create a temporary file from which to copy the header
cat > temp_copyright_header.txt << EOF
// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

EOF

# insert the header if it's not already here
for i in $(find src -not \( -path src/extern -prune \) -name '*.cpp' -or -name '*.h' -or -name '*.hxx')
do
  if ! grep -q "Copyright ETH Zurich" $i
  then
    cat temp_copyright_header.txt $i >$i.new && mv $i.new $i
  fi
done

# remove the temporary file
rm temp_copyright_header.txt
