# Manifold-valued Image Restoration Toolbox

![](docs/MVIRT_banner.png)A subset of the [Camino Toolkit](http://camino.cs.ucl.ac.uk) Human Head DT-MRI data set using the [viridis](http://bids.github.io/colormap/) colormap.

maintained by
[Ronny Bergmann](imagepro/members/bergmann/)

written by
[Ronny Bergmann](http://www.mathematik.uni-kl.de/imagepro/members/bergmann/)
[Johannes Persch](http://www.mathematik.uni-kl.de/imagepro/members/persch/)

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
or when working with phase data, i.e. on [the circle 𝕊¹](manifolds/S1.md), or [Diffusion Tensor Imaging](https://en.wikipedia.org/wiki/Diffusion_MRI) (DTI), where
every data items are an n ⨉ n
[symmetric positive definite matrices, 𝒫(n)](manifolds/SymPosDef).

see [ronnybergmann.net/mvirt/](http://ronnybergmann.net/mvirt) for the complete
documentation.
