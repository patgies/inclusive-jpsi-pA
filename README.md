# Inclusive J/psi in pp and pA

# Inclusive D0 photoproduction

This code provides the cross section for inclusive J/psi production in proton-proton and proton-lead collisions. 

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [Phys.Rev.D 110 (2024) 9, 094035](https://inspirehep.net/literature/2824777) 

Dipoles amplitudes from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole) 
***

## Building
Requires
- Cmake
- GSL

How to compile:
```
mkdir build
cd build
cmake ..
make
```
This generates a library build/lib/libamplitude.a that you can link in your own program.

### Questions and comments
Please send an email to patricia.p.gimenoestivill@jyu.fi