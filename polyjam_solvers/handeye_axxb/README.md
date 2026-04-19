# Hand-eye AX=XB solver

This directory contains the public wrapper and the synthetic test harness for the generated Polyjam hand-eye solver.

The top-level generator creates the quaternion-only Groebner core in the sibling directory:

```text
polyjam_solvers/handeye_axxb_q/handeye_axxb_q.cpp
polyjam_solvers/handeye_axxb_q/handeye_axxb_q.hpp
```

This directory contains the hand-written public wrapper:

```text
polyjam_solvers/handeye_axxb/handeye_axxb.cpp
polyjam_solvers/handeye_axxb/handeye_axxb.hpp
```

The wrapper first calls the generated quaternion core and then reconstructs translation from the normal equations of the original Euclidean least-squares objective.

Build the generated core with:

```bash
cmake --build build --target generate_handeye_axxb_solver
```

Build and run the synthetic test with:

```bash
cmake --build build --target test_handeye_axxb
cmake --build build --target run_handeye_axxb_test
```

The public solver function accepts `std::vector<Eigen::Matrix4d> & A` and `std::vector<Eigen::Matrix4d> & B`, using `POLYJAM_HAND_EYE_PAIRS` pairs compiled into the generator. Each returned `Eigen::Matrix<double,7,1>` is `[q_w, q_x, q_y, q_z, t_x, t_y, t_z]^T`, where the rotation is represented by a unit quaternion.
