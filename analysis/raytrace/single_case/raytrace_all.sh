#!/bin/bash

for filename in ./*.pov; do
  `povray -I$filename -D0 -H1440 -W2560`
done
