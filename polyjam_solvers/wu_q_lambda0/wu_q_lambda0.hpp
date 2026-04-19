#ifndef POLYJAM_WU_Q_LAMBDA0_PUBLIC_HPP_
#define POLYJAM_WU_Q_LAMBDA0_PUBLIC_HPP_

#include <Eigen/Eigen>
#include <array>
#include <vector>

namespace polyjam
{
namespace wu_q_lambda0
{

// Solves the lambda-free subsystem obtained from
//   W*u(q) = Q*q,  q^T q = 1,
// after selecting three of the four residual rows.  Here u(q) is the
// 64-vector q kron q kron q in row-major lexicographic order:
//   u[(a*4+b)*4+c] = q[a]*q[b]*q[c], q=[qw,qx,qy,qz]^T.
// Each returned vector is [qw,qx,qy,qz]^T and is normalized.
void solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    std::vector<Eigen::Matrix<double,4,1> > & solutions);

// Same solver, but with explicit residual-row selection.  The generated core
// is row-agnostic: this wrapper copies the requested rows into the three-row
// coefficient vector expected by the core solver.
void solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const std::array<int,3> & selected_rows,
    std::vector<Eigen::Matrix<double,4,1> > & solutions);

}
}

#endif /* POLYJAM_WU_Q_LAMBDA0_PUBLIC_HPP_ */
