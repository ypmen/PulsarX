## Requirements:
- To generate figures, python,matplotlib,numpy or PlotX should be installed.
- To calculate GL and GB, sofa library should be installed.
- To calculate YMW16 DM and distance, environment variable YMW16_DIR should be set.

## Installation:
1) ./bootstrap
2) ./configure --prefix=/install\_path LDFLAGS="-L/path_to_sofa" CPPFLAGS="-I/path_to_sofa"
3) make and make install

## Usage:
- use -h to see help
- set environment variable YMW16_DIR to $top_src/src/ymw16
- --template /path_to_fits_template ($top_src/include/template contains some examples).
- --candfile /path_to candfile, eg  
```
   #id   dm acc  F0 F1 S/N  
   1  100   0  1000  0  10  
   2  200   0  1  0  100  
```
## Docker:
```
docker pull ypmen/pulsarx
```
## Acknowledgement
Thanks for very helpful suggestions and bug reports from Ewan Barr, Colin Clark, Emma Carli, Shalini Sengupta, Vivek Venkatraman Krishnan,  Miquel Colom i Bernadich  and TRAPUM project.
