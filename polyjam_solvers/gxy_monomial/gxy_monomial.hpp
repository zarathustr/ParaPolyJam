#ifndef POLYJAM_GXY_MONOMIAL_PUBLIC_HPP_
#define POLYJAM_GXY_MONOMIAL_PUBLIC_HPP_

#include <Eigen/Eigen>
#include <vector>

namespace polyjam
{
namespace gxy_monomial
{

// Solves the six-equation MATLAB monomial system
//   eqs(1:3) = M1 * mons1(gX,gY) = 0,
//   eqs(4:6) = M2 * mons4(gX,gY) = 0,
// with unknown order [gX1,gX2,gX3,gY1,gY2,gY3]^T.
// M1 is 3x40 and M2 is 3x56, using the exact monomial ordering stated in
// the corresponding generator and README.
void solve(
    const Eigen::Matrix<double,3,40> & M1,
    const Eigen::Matrix<double,3,56> & M2,
    std::vector<Eigen::Matrix<double,6,1> > & solutions);

}
}

#endif /* POLYJAM_GXY_MONOMIAL_PUBLIC_HPP_ */
