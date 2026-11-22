# Changelog

All notable changes to this package. Versions follow semver, independent of the Nabi toolchain.

## 0.1.1 (2026-11-22)

- The crate docs open with a tutorial: building a matrix from triplets
  and solving a large system with `cg`, `ic0` and `cholesky`. It moved
  here from the Nabi docs site; each step is a doc test.

## 0.1.0 (2026-11-21)

First release.

- `Csr`, `Csc` and `Coo` matrices of f64, with conversions to each other
  and to dense.
- Building: `from_triplets` (repeated positions add up), `from_dense`,
  `eye`, `diags`, `kron`.
- Products: `matvec`, `rmatvec`, sparse x sparse `matmul`, `matmul_dense`,
  `add`, `sub`, `scale`.
- Krylov solvers `cg`, `gmres` and `bicgstab` on any `Operator`, with
  `Jacobi`, `ilu0` and `ic0` preconditioners.
- Direct `cholesky` and `lu` after a fill-reducing ordering (`rcm`,
  `min_degree`).
- The `ops` feature (default): linalg's `@`, `\` and `†` on sparse matrices.
- `SparseError` for bad input. Pure Nabi, every target including wasm.
  Needs Nabi 0.1 or later; depends on `linalg`.
