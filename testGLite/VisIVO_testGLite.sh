#!/bin/sh
export LFC_HOST=lfc-01.ct.trigrid.it
hostname
lfc-ls -l /grid/cometa/ube

VisIVOImporter --fformat ascii --VO cometa  --lfnout lfn://grid/cometa/ube/mrvbt16.bin  lfn://grid/cometa/ube/mrvbt16.ascii
VisIVOFilters --op pointproperty --resolution 32 32 32 --points X Y Z --outcol density --append --VO cometa --file lfn://grid/cometa/ube/mrvbt16.bin
VisIVOViewer --x X --y Y --z Z --cycle cycle.par --showbox --showaxes --color --colorscalar density --VO cometa --lfnout lfn://grid/cometa/ube/VVImage  lfn://grid/cometa/ube/mrvbt16.bin


