# ManImRes – MANifold-valued IMage REStoration

In many applications data is constrained to certain manifolds. The simplest
example is the circle \(\mathbb S^1\), i.e. phase valued data like in InSAR
imaging. Every pixel of an image is here given as an angle or a point on the
circle. Combining these with real valued data we obtain data like in the HSV
color space, where the first value is phase/angle valued and the other two are
real valued. We obtain the product manifold of
\(\mathbb S^1\times \mathbb R^2\)-valued data. 

This `MatLab` package provides algorithms for restoring, i.e. denoising and inpainting images


## Available Manifolds
* The Circle `S1` for phase valued data
* The vector-valued (product) manifold `S1mRn`

## Installation, Initialization

For installation just place the folder ManImRes somewhere in your preferred  
directory. On Startup be sure to have `ManImRes` as your base directory or  
call `initMVIRT();` from anywhere. By default it adds all necessary folders  
to the MATLAB path and initializes the debug helping functions.

## Getting Started
To get started with the remaining features have a look at the `examples/` folder
or take a look at the `algorithms/` and `manifolds/` functions and classes.

## Authors
The package was initialized by Ronny Bergmann <bergmann@mathematik.uni-kl.de>
and some algorithms were implemented by Johannes Persch <persch@mathematik.uni-kl.de>.

## References
The examples provided in this toolbox are published in the following articles

### Journal Papers
1. M. Bačák, R. Bergmann, G. Steidl, A. Weinmann *A Second Order Non-Smooth Variational Model for Restoring Manifold-Valued Images.*, accepted for publication to SIAM Journal on Scientific Computing, 2016.

3. R. Bergmann, F. Laus, G. Steidl, A. Weinmann: *Second Order Differences of S1-valued Data and Applications in Variational Image Denoising.*, J. SIAM Imaging Sci. 7(4), 2916–2953, 2015. [DOI](http://dx.doi.org/10.1137/140969993), [arXiv](http://arxiv.org/pdf/1405.5349v1.pdf).

### Conference Proceedings
1. R. Bergmann, A. Weinmann: *Inpainting of Cyclic Data Using First and Second Order Differences Second Order Differences, in: EMMCVPR 2015.*, 2015. [DOI](http://dx.doi.org/10.1007/978-3-319-14612-6_12), [arXiv](http://arxiv.org/pdf/1410.1998v1.pdf)

### Data Sources
Both the images “sailboat on lake” and “Pappers” are taken from the SIPI image database (4.2.06 and 4.2.07) available at http://sipi.usc.edu/DATABASE/database.php?volume=misc.

The Mount Vesuvius Data is taken from https://earth.esa.int/workshops/ers97/program-details/speeches/rocca-et-al/

and the Camino examples use a part of the Camino data set available from http://camino.cs.ucl.ac.uk