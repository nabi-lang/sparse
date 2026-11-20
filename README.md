# sparse

[![CI](https://github.com/nabi-lang/sparse/actions/workflows/ci.yml/badge.svg)](https://github.com/nabi-lang/sparse/actions/workflows/ci.yml) [![packages.nabi-lang.org](https://img.shields.io/github/v/tag/nabi-lang/sparse?label=packages.nabi-lang.org&sort=semver)](https://packages.nabi-lang.org/packages/sparse) [![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

Sparse matrices for Nabi: compressed rows (`Csr`), compressed columns
(`Csc`) and triplets (`Coo`), with conversions, products, the Krylov
solvers CG, GMRES and BiCGSTAB with Jacobi, ILU(0) and IC(0)
preconditioners, and sparse Cholesky and LU after a fill-reducing
ordering. Values are f64. Pure Nabi, so it runs on every target, wasm
included.

## Add it

```sh
nabi add sparse
```

Needs Nabi 0.1 or later.

## Usage

```nabi
use core::array::Vector;

// The 2-d Laplacian on a 100 x 100 grid: 10 000 unknowns, 49 600 entries.
let t = try sparse::diags([[2.0], [-1.0], [-1.0]], [0, 1, -1], [100, 100]);
let i = sparse::eye(100);
let a = sparse::kron(i, t) + sparse::kron(t, i);
let b = Vector::<f64>::ones([10000]);

let s = try sparse::cg(a, b, pre: try sparse::ic0(a));  // iterative
println(f"{s.iterations} steps, residual {s.residual:e}");

let f = try sparse::cholesky(a);                       // direct, minimum degree order
let x = try f.solve(b);
let y = try (try sparse::lu(a)).solve(b);              // any square matrix
```

Assemble with triplets (repeated positions add up, as in finite-element
assembly), then convert once:

```nabi
let c = sparse::Coo::new([3, 3]);
for k in 0..3 { try c.push(k, k, 2.0); }
try c.push(0, 2, 1.0);
let a = c.to_csr();
let y = try a.matvec([1, 2, 3]);                       // lists, ranges, Vectors
```

Every function that takes a vector takes `linalg::AsVector<f64>`: a
`Vector<f64>`, a `List` or fixed array of `f64`, `f32`, `i64`, `i32` or
`usize`, or an `i64` range.

The solvers take any `Operator`: `Csr`, `Csc`, a dense `Matrix<f64>`, or
a function wrapped by `sparse::operator(n, f)` for a matrix you never
form. A preconditioner is any `Preconditioner`; `Cholesky` and `Lu` of a
nearby matrix work as one too.

Errors are one `SparseError` (`ShapeMismatch`, `OutOfBounds`,
`BadStructure`, `NotSquare`, `NotSymmetric`, `Singular`,
`NotPositiveDefinite`, `NoConvergence`, `Breakdown`, `NotFinite`); each
message says what to call instead.

## Operators

`+`, `-`, unary `-` and `*` by a scalar work on `Csr` and `Csc`. With
the `ops` feature (on by default) linalg's operator tokens work too,
after `use linalg::ops::*`:

| Operator         | Meaning                                         | Named form              |
|------------------|-------------------------------------------------|-------------------------|
| `a @ x`, `a ⊡ x` | product with a vector, dense or sparse matrix   | `matvec`, `matmul`, `matmul_dense` |
| `a \ b`, `a ⧵ b` | solve `a x = b` by sparse LU (a `Result`)       | `lu(a).solve(b)`        |
| `a†`             | transpose, sharing the arrays                   | `t`                     |

Operators cannot throw, so `@` and `+` panic on a shape mismatch, as
core's broadcast operators do; the named forms throw.

## What is in it

| Area            | Items                                                               |
|-----------------|---------------------------------------------------------------------|
| Types           | `Csr`, `Csc`, `Coo`, `Sparse` (any of the three)                    |
| Building        | `from_triplets`, `from_dense`, `new` (compressed arrays), `eye`, `diags`, `kron` |
| Conversions     | `to_csr`, `to_csc`, `to_coo`, `to_dense`, `t`, `transpose`          |
| Products        | `matvec`, `rmatvec`, `matmul` (sparse x sparse), `matmul_dense`, `add`, `sub`, `scale` |
| Access          | `get`, `rows`, `cols`, `diagonal`, `indptr`, `indices`, `data`, `nnz` |
| Iterative       | `cg`, `gmres`, `bicgstab`, `Solution`, `Operator`, `operator`       |
| Preconditioners | `Identity`, `Jacobi`, `ilu0`, `ic0`, `Preconditioner`               |
| Direct          | `cholesky`, `lu`, `Cholesky`, `Lu`, `Order`, `rcm`, `min_degree`    |

Methods: triplets to compressed form by two counting sorts; Gustavson's
row-by-row sparse product; CG (Hestenes-Stiefel), restarted GMRES with
modified Gram-Schmidt and Givens rotations and BiCGSTAB (van der Vorst),
both right preconditioned so the residual tested is the true one;
ILU(0) and IC(0) without fill (Saad); up-looking sparse Cholesky through
the elimination tree and left-looking Gilbert-Peierls LU with threshold
partial pivoting (Davis, Direct Methods for Sparse Linear Systems);
reverse Cuthill-McKee from a pseudo-peripheral vertex and exact minimum
degree on the elimination graph.

## Speed

One Apple M-series core, `nabi bench` (-O2), on the 2-d five-point
Laplacian:

| Bench                                   | Time     |
|-----------------------------------------|----------|
| SpMV, 300 x 300 grid (448 800 entries)  | 0.70 ms (CSR), 0.55 ms (CSC); about 650-800 M entries/s |
| `a.matmul(a)`, 100 x 100 grid           | 1.6 ms   |
| `cg`, 100 x 100 grid, rtol 1e-8         | 20 ms (19 ms with `ic0`) |
| `gmres` with `ilu0`, 100 x 100 grid     | 38 ms    |
| `min_degree`, 100 x 100 grid            | 51 ms    |
| `cholesky` (minimum degree, order included) | 50 ms |
| `cholesky` with `Order::Rcm`            | 26 ms    |
| `lu` (minimum degree, order included)   | 68 ms    |

Ballpark reference, not measured side by side: scipy.sparse's CSR
matvec runs at a similar rate per entry; CHOLMOD and SuperLU factor this
matrix in a few milliseconds, most of the gap being the ordering, which
here is exact minimum degree without AMD's supervariables and
approximate degrees. On big 3-d problems pass `order: Order::Rcm` or use
an iterative solver.

## Not in 0.1.0

Complex or f32 values; AMD/COLAMD and nested dissection orderings;
supernodal or multifrontal factorizations; eigenvalue solvers (Lanczos,
Arnoldi); MINRES and LSQR; block or parallel kernels.

## Docs

https://packages.nabi-lang.org/packages/sparse/latest/docs

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at
your option.
