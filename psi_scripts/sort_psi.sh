#!/usr/bin/bash

dir=Cells_Cultured_Fibroblasts

sort -t_ -k1.4,1n -k2,2n $dir/all.A.psi.filtered.tsv > $dir/all.A.psi.filtered.sorted.tsv
