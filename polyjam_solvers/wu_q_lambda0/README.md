# Lambda-free W*u = Q*q quaternion solver

This directory contains the public wrapper and synthetic test for the Polyjam-generated solver obtained by setting `lambda = 0` in

```text
W * u(q) = (Q - lambda * I_4) * q,
q^T q = 1.
```

After zeroing `lambda`, the four scalar residual rows become

```text
W_i * u(q) - Q_i * q = 0,  i = 0..3,
```

and the fifth equation is the unit-quaternion constraint.  The generated Groebner core solves a square four-equation system made from three selected residual rows plus

```text
q_w^2 + q_x^2 + q_y^2 + q_z^2 = 1.
```

The 64-vector is assumed to be the cubic Kronecker monomial vector

```text
u[(a*4+b)*4+c] = q[a] * q[b] * q[c],  a,b,c in {0,1,2,3}.
```

The public wrapper returns normalized candidate quaternions in the order

```text
[q_w, q_x, q_y, q_z]^T.
```

The default wrapper uses rows `{0,1,2}`.  An overload accepts any three distinct rows from `{0,1,2,3}` and packs the selected rows into the coefficient vector used by the generated core.

Generate and test from the repository root with:

```bash
cmake -S . -B build -DPOLYJAM_MACAULAY_COMMAND=M2
cmake --build build --target generate_wu_q_lambda0_solver
cmake --build build --target test_wu_q_lambda0
cmake --build build --target run_wu_q_lambda0_test
```
