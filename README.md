# Bound-States

Using this code, you can calculate the energy levels for any bound state problem. This code is written in FORTRAN 90. To use this code you need to install `gfortran` in your system.

If you use Ubuntu, you need `gfortran` , `liblapack-dev` , `libblas-dev` and `libopenblas-dev` packages as prerequisites. To install them use this command 

```
sudo apt install liblapack-dev libblas-dev libopenblas-dev -y
```

Next, to run this code, 

```
gfortran -o w1D well1D.f90 -lopenblas -fopenmp
```
You will get the energy levels for a one-dimensional potential well. For any other potential define V(i) accordingly.
