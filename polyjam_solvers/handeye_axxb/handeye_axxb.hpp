#ifndef POLYJAM_HANDEYE_AXXB_PUBLIC_HPP_
#define POLYJAM_HANDEYE_AXXB_PUBLIC_HPP_

#include <Eigen/Eigen>
#include <vector>

namespace polyjam
{
namespace handeye_axxb
{

// Public hand-eye solver wrapper.
// The generated Groebner solver first finds quaternion critical points of the
// determinant-cleared reduced objective.  This wrapper then reconstructs the
// least-squares translation from the translation normal equations.
void solve(
    std::vector<Eigen::Matrix4d> & A,
    std::vector<Eigen::Matrix4d> & B,
    std::vector<Eigen::Matrix<double,7,1> > & solutions);

}
}

#endif /* POLYJAM_HANDEYE_AXXB_PUBLIC_HPP_ */
