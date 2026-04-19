#ifndef POLYJAM_WU_Q_LAMBDA_PUBLIC_HPP_
#define POLYJAM_WU_Q_LAMBDA_PUBLIC_HPP_

#include <Eigen/Eigen>
#include <vector>

namespace polyjam
{
namespace wu_q_lambda
{

// Solves W*u(q) = (Q - lambda*I_4)*q with q^T q = 1.
// Here u(q) is the 64-vector q kron q kron q in row-major lexicographic order:
//   u[(a*4+b)*4+c] = q[a]*q[b]*q[c], q=[qw,qx,qy,qz]^T.
// Each returned vector is [lambda,qw,qx,qy,qz]^T.
void solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    std::vector<Eigen::Matrix<double,5,1> > & solutions);

}
}

#endif /* POLYJAM_WU_Q_LAMBDA_PUBLIC_HPP_ */
