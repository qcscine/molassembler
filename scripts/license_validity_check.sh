#!/bin/bash

# Look for all source files, and then check if the copyright / license header is present
let retvalue=0
for i in $(find src -not \( -path src/extern -prune \) -name '*.cpp' -or -name '*.h' -or -name '*.hxx')
do
  if ! grep -q -i "Copyright ETH Zurich" $i
  then
    echo $i
    let retvalue=1
  fi
done

if [ $retvalue -eq 1 ]; then
  exit 1
fi

exit 0
