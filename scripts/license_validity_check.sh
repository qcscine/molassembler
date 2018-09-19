#!/bin/bash

# Look for all source files, and then check if the copyright / license header is present
for i in $(find src -not \( -path src/extern -prune \) -name '*.cpp' -or -name '*.h' -or -name '*.hxx')
do
  if ! grep -q "Copyright ETH Zurich" $i
  then
    echo "ETH license is not present in all files!"
    echo "Please run the scripts/add_license.sh script."
    exit 1
  fi
done
