# gxy_monomial solver

This folder contains the public Eigen wrapper and synthetic test for the MATLAB-style six-equation monomial system

```text
M1(1:3,:) * mons1(gX,gY)^T = 0
M2(1:3,:) * mons4(gX,gY)^T = 0
```

with unknown vector

```text
[gX1,gX2,gX3,gY1,gY2,gY3]^T.
```

`M1` is `3x40` and `M2` is `3x56`.  The monomial ordering is exactly the MATLAB ordering implemented in `polyjam_workspace/gxy_monomial/gxy_monomial.cpp`.

From the repository root:

```bash
cmake -S . -B build -DPOLYJAM_MACAULAY_COMMAND=M2
cmake --build build --target polyjam_generate_gxy_monomial
cmake --build build --target generate_gxy_monomial_solver
cmake --build build --target test_gxy_monomial
cmake --build build --target run_gxy_monomial_test
```

The test loops over RNG seeds, constructs consistent synthetic `M1` and `M2`, and reports min/max/mean/median statistics for the best solution error, polynomial residual, number of returned solutions, and solve time in milliseconds.  You can pass an optional seed count and first seed:

```bash
./test_gxy_monomial 100 0
```
