# Manifold-valued Image Restoration Toolbox

----

:information_source:

This package is no longer maintained :zzz:.
If you want to do optimisation on manifolds, depending on the programming language you use, please consider [Manopt.jl](https://manoptjl.org), which is the successor of this project written in Julia, [Manopt](https://www.manopt.org) in Matlab or [pymanopt](https://www.pymanopt.org/) in Python.

----

![](docs/MVIRT_banner.png)A subset of the [Camino Toolkit](http://camino.cs.ucl.ac.uk) Human Head DT-MRI data set using the [viridis](http://bids.github.io/colormap/) colormap.


maintained by
Ronny Bergmann

written by
Ronny Bergmann
Johannes Persch

In many application measured data appears nonlinear, i.e. restricted in a certain range or equipped with a different distance measure.

This toolbox provides an easy access to image processing tasks for such data,
where the values of each measurement are points on a manifold.
You can get started by downloading the
[source code from Github](https://github.com/kellertuer/MVIRT/archive/master.zip)
or by cloning the [git-repository](https://github.com/kellertuer/MVIRT) using

````
    git clone git@github.com:kellertuer/MVIRT.git
````

Examples are
[InSAR imaging](http://de.wikipedia.org/wiki/Interferometric_Synthetic_Aperture_Radar)
or when working with phase data, i.e. on the circle 𝕊¹ or [Diffusion Tensor Imaging](https://en.wikipedia.org/wiki/Diffusion_MRI) (DTI), where
every data items are an n ⨉ n symmetric positive definite matrices, 𝒫(n).
