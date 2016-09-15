/*
 *
 * Manifold functions of Symmetric positive definite matrices, elementary,
 *not working with vectors/arrays of matrices, in Eigen-Library notation for convenience
 * ---
 * ManImRes ~R. Bergmann, 2015-11-06
 */
#ifndef MANIFOLDHN_H
#define MANIFOLDHN_H

#include "math.h"
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

VectorXd mHnExp(VectorXd X, VectorXd V);
VectorXd mHnExp(VectorXd X, VectorXd V,double t);
/* Y = mHnExp(X,V) Compute one exponential map for
 * INPUT: X: Hn Vector, V a vector in TXM
 * OUTPUT: Y resulting Hn Vector Y. V may be shorten by t before, if not given, t=1.0 is used.
 */
VectorXd mHnLog(VectorXd X, VectorXd Y);
/* V = mHnLog(X,Y)
 * Compute one logarithmic map (inverts Exp) for two points X,Y on Hn,
 * Computes the vector in TXM pointing to Y
 * INPUT: X,Y points (unit vector in Rn+1).
 * OUTPUT: Direction V
 */
double mHnDot(VectorXd V, VectorXd W);
/* d = mHnDist(X,Y)
 * Compute the dot product of the two points X,Y on Hn.
 */
double mHnDist(VectorXd X, VectorXd Y);
/* d = mHnDist(X,Y)
 * Compute the distance of the two points X,Y on Hn.
 */
/* VectorXd mHnGradX(VectorXd X, VectorXd Y, VectorXd Z); */
/* compute the (sub) gradient w.r.t. to X 
 * INPUT: X,Y,Z points (unit vectors in Rn+1)
 * OUTPUT: Direction of the gradient of the second order difference at X.
 */
VectorXd mHnMean(VectorXd *F, double *W, double E, double I,size_t L);
// Caluculates the mean of the j points in F with weights W stops after I iteration or if change less than E
#endif /* MANIFOLDHn_H */