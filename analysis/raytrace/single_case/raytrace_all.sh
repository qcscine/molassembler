#!/bin/bash

for filename in ./*.pov; do
  `povray -I$filename -D0 -H768 -W1366`
done
