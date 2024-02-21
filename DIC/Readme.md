# DIC scripts from Thomas

The scripts here are using a sub-program of the [Micmac](https://github.com/micmacIGN/micmac), called `mm3d MM2DPosSism` to perform DIC.

## `parMicmac.sh`
This small bash script performs DIC on a sequence of images, sequentially, all with respect to one frame (but this can easily be changed).
The line performing DIC actually is (with more generic arguments):

`mm3d MM2DPosSism ref.tif moving.tif DirMEC=../outFolder/$i/ SzW=$ss`
where:
  - `ref.tif` is the reference image (the first one usually)
  - `moving.tif` is the moved image
  - `../outFolder/$i/` is the folder will be put. It has to be looped over otherwise subsequent results are erased
  - `SzW` is the window size with 2 meaning a 5x5 (1+2+2) window, 3 a 7x7 (1+3+3)
  
This produces a number of files with absurd names, the only ones that matters are:
  - `Px1_Num6_DeZoom1_LeChantier.tif`: the moving field in the X direction
  - `Px2_Num6_DeZoom1_LeChantier.tif`: the moving field in the Y direction
  
All the other files are erased by the `parMicmac` script.

Parameters of this script (`ss`, `prename`, the folder after `cd` and the loop boundaries need to be adapted for every use case.

## `gatherData.sh`
This is a small script that copies and renames the results of `parMicmac` into two folders for the X and Y fields.
