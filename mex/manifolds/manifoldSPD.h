/*
 *
 * Manifold functions of Symmetric positive definite matrices, elementary,
 *not working with vectors/arrays of matrices, in Eigen-Library notation for convenience
 * ---
 * Manifold Valued Image Restoration 1.0
 * R. Bergmann, 2015-04-12
 */
#ifndef MANIFOLDSPD_H
#define MANIFOLDSPD_H

#include "math.h"
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

/* Y = mSPDExp(X,V) Compute one exponential map for
 * INPUT: X: SPD Matrix, V a matrix in TXM
 * OUTPUT: Y resulting SPD Y. V may be shorten by t before, if not given, 1.0 is used.
 */
MatrixXd mSPDExp(MatrixXd X, MatrixXd V);
MatrixXd mSPDExp(MatrixXd X, MatrixXd V,double t);
/* V = mSPDLog(X,Y)
 * Compute one logarithmic maps (inverts Exp) for two points X,Y SPD,
 * Compute the vector in TXM pointing to Y
 * INPUT: X,Y points (SPDs). 
 * OUTPUT: Direction V
 */
MatrixXd mSPDLog(MatrixXd X, MatrixXd Y);
/* d = mSPDDist(X,Y)
 * Compute the distance of SPDs X,Y
 */
double mSPDDist(MatrixXd X, MatrixXd Y);
/* d = mSPDDot(X,V,W)
 * compute dot product in the TXM of the elements V and W.
 */
double mSPDDot(MatrixXd X, MatrixXd V, MatrixXd W);
MatrixXd mSPDGradX(MatrixXd X, MatrixXd Y, MatrixXd Z);
/* compute the (sub) gradient w.r.t. to X 
 * INPUT: X,Y,Z points (SPDs)
 * OUTPUT: Direction of the gradient of the second order difference at X.
 */
MatrixXd mSPDParallelTransport(MatrixXd X, MatrixXd Y, MatrixXd V);
MatrixXd mSPDParallelTransport(MatrixXd X, MatrixXd Y, MatrixXd V, double t);

MatrixXd mSPDMean(MatrixXd *F, double *W, double E, double I,size_t L);
// Caluculates the mean of the j points in F with weights W stops after I iteration or if change less than E
#endif /* MANIFOLDSPD_H */