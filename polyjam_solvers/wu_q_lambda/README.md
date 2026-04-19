# W*u = (Q - lambda I) q solver

This directory contains the public wrapper and synthetic test for the Polyjam-generated solver of

```text
W * u(q) = (Q - lambda * I_4) * q
q^T q = 1
```

where `W` is `4x64`, `Q` is `4x4`, `lambda` is a scalar unknown, and
`q=[qw,qx,qy,qz]^T` is a unit quaternion.  The 64-vector is assumed to be the cubic Kronecker monomial vector

```text
u[(a*4+b)*4+c] = q[a] * q[b] * q[c],  a,b,c in {0,1,2,3}.
```

The public wrapper returns candidate vectors in the order

```text
[lambda, qw, qx, qy, qz]^T
```

Generate and test from the repository root with:

```bash
cmake -S . -B build -DPOLYJAM_MACAULAY_COMMAND=M2
cmake --build build --target generate_wu_q_lambda_solver
cmake --build build --target test_wu_q_lambda
cmake --build build --target run_wu_q_lambda_test
```
