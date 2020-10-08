#!/bin/bash

for filename in ./*progress*.csv; do
  `Rscript plot_progress.R $filename`
done
