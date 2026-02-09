# Linear Algebra Compendium

<div style="text-align: center; margin: 2em 0;">
<strong style="font-size: 1.3em;">From freshman to PhD — a comprehensive, systematic linear algebra knowledge base</strong>
</div>

---

## About This Site

This site is a **comprehensive, systematic, and self-contained** linear algebra knowledge base, covering all core linear algebra content from undergraduate introductory courses through doctoral-level research. Whether you are a beginner or a researcher, you can find the knowledge you need here.

The content is organized into eight parts with 90 chapters, arranged by increasing difficulty and logical coherence. Each chapter contains complete definitions, theorems, proofs, and examples, striving for rigor and clarity.

---

## Content Guide

### Part I: Foundations of Linear Algebra <span class="difficulty-tag beginner">Undergraduate Basics</span>

Suitable for freshmen and sophomores, covering all core content of an introductory linear algebra course.

| Chapter | Overview |
|---------|----------|
| [Chapter 0  Polynomial Algebra](part1/ch00-polynomials.md) | Polynomial rings, divisibility, GCD, irreducible factorization |
| [Chapter 1  Systems of Linear Equations](part1/ch01-linear-equations.md) | Solving linear systems, Gaussian elimination, solution structure |
| [Chapter 2  Matrices and Matrix Operations](part1/ch02-matrices.md) | Matrix operations, inverse matrices, block matrices, elementary matrices, rank |
| [Chapter 3  Determinants](part1/ch03-determinants.md) | Definition and properties of determinants, cofactor expansion, Cramer's rule |
| [Chapter 4  Vector Spaces](part1/ch04-vector-spaces.md) | Vector space axioms, subspaces, bases and dimension, rank-nullity theorem |
| [Chapter 5  Linear Transformations](part1/ch05-linear-transformations.md) | Linear maps, kernel and image, matrix representation, change of basis |
| [Chapter 6  Eigenvalues and Eigenvectors](part1/ch06-eigenvalues.md) | Characteristic polynomial, diagonalization, Cayley–Hamilton theorem |
| [Chapter 7  Orthogonality and Least Squares](part1/ch07-orthogonality.md) | Orthogonal sets, Gram–Schmidt, orthogonal projection, least squares |

### Part II: Intermediate Linear Algebra <span class="difficulty-tag intermediate">Advanced Undergraduate</span>

Suitable for sophomores through seniors and early graduate students, delving deeper into core linear algebra theory.

| Chapter | Overview |
|---------|----------|
| [Chapter 8  Inner Product Spaces](part2/ch08-inner-product-spaces.md) | General inner product spaces, orthogonal complements, adjoint operators, spectral theorem |
| [Chapter 9  Quadratic and Bilinear Forms](part2/ch09-quadratic-forms.md) | Quadratic forms, bilinear forms, symplectic spaces, Hermitian forms |
| [Chapter 10  Matrix Decompositions](part2/ch10-matrix-decompositions.md) | LU, Cholesky, QR, Schur decompositions |
| [Chapter 11  Singular Value Decomposition](part2/ch11-svd.md) | SVD theory and applications, low-rank approximation, pseudoinverse |
| [Chapter 12  Jordan Normal Form](part2/ch12-jordan-form.md) | Generalized eigenvectors, Jordan blocks, minimal polynomial |
| [Chapter 13  Matrix Functions](part2/ch13-matrix-functions.md) | Matrix exponential, matrix logarithm, matrix power series |
| [Chapter 13A  Quotient Spaces and Dual Spaces](part2/ch13a-quotient-dual-spaces.md) | Quotient spaces, dual spaces, annihilators, transpose maps, canonical isomorphism |
| [Chapter 13B  Lambda-Matrices and Rational Canonical Form](part2/ch13b-lambda-rational-form.md) | Lambda-matrices, Smith normal form, invariant factors, rational canonical form |

### Part III: Advanced Linear Algebra <span class="difficulty-tag advanced">Graduate</span>

Suitable for master's and doctoral students, covering matrix analysis and advanced theory.

| Chapter | Overview |
|---------|----------|
| [Chapter 14  Matrix Analysis](part3/ch14-matrix-analysis.md) | Matrix sequences and series, spectral radius, Gershgorin's theorem |
| [Chapter 15  Norms and Perturbation Theory](part3/ch15-norms-perturbation.md) | Matrix norms, condition numbers, eigenvalue perturbation |
| [Chapter 16  Positive Definite Matrices](part3/ch16-positive-definite.md) | Equivalent conditions for positive definiteness, Schur complement, Löwner partial order |
| [Chapter 17  Nonnegative Matrices and Perron–Frobenius Theory](part3/ch17-nonnegative-matrices.md) | Perron–Frobenius theorem, irreducible matrices, stochastic matrices |
| [Chapter 18  Matrix Inequalities](part3/ch18-matrix-inequalities.md) | Eigenvalue inequalities, trace inequalities, determinantal inequalities, majorization |
| [Chapter 19  Kronecker Product and Vec Operator](part3/ch19-kronecker.md) | Kronecker product, Vec operator, and applications to matrix equations |
| [Chapter 20  Matrix Equations](part3/ch20-matrix-equations.md) | Sylvester equation, Lyapunov equation, Riccati equation |

### Part IV: Algebraic Structures & Generalizations <span class="difficulty-tag research">Doctoral/Research</span>

Algebraic foundations and structural generalizations of linear algebra.

| Chapter | Overview |
|---------|----------|
| [Chapter 21  Multilinear Algebra and Tensors](part4/ch21-multilinear-algebra.md) | Dual spaces, tensor products, exterior algebra, tensor decomposition |
| [Chapter 48  Modules over PIDs](part7/ch48-modules-pid.md) | Structure theorem for modules, invariant factors, unified Jordan/rational form |
| [Chapter 49  Exterior Algebra & Grassmannian](part7/ch49-exterior-algebra.md) | Wedge product, exterior powers, compound matrices, Plücker coordinates |
| [Chapter 50  Clifford Algebras & Geometric Algebra](part7/ch50-clifford-algebra.md) | Geometric product, spinors, Spin groups, Bott periodicity |
| [Chapter 51  Quaternion Matrices](part7/ch51-quaternion-matrices.md) | Quaternion algebra, left/right eigenvalues, quaternion SVD |
| [Chapter 52  Linear Algebra over Finite Fields](part7/ch52-finite-fields.md) | Vector spaces over GF(q), Gaussian binomials, GL(n,q), linear codes |
| [Chapter 53  Symplectic & Hamiltonian Matrices](part7/ch53-symplectic-matrices.md) | Symplectic group, Hamiltonian matrices, Williamson's theorem |
| [Chapter 54  Quiver Representations](part7/ch54-quiver-representations.md) | Quiver representations, Gabriel's theorem, Dynkin diagrams |
| [Chapter 55  Matrix Groups & Classical Lie Groups](part7/ch55-matrix-groups.md) | GL, SL, O, SO, U, SU, Sp, Lie algebras, exponential map, BCH formula |
| [Chapter 56  Pfaffian](part7/ch56-pfaffian.md) | Pfaffian properties, pf(A)²=det(A), FKT algorithm, perfect matchings |

### Part V–VIII: Advanced Topics

Parts V through VIII cover matrix structures, operator theory, modern frontiers, and applications across 66 chapters. These parts are currently available in Chinese; English translations are in progress.

---

## Notation Conventions

This site uses the following standard mathematical notation:

| Symbol | Meaning |
|--------|---------|
| $\mathbb{R}, \mathbb{C}, \mathbb{F}$ | Real numbers, complex numbers, general field |
| $\mathbb{R}^n, \mathbb{R}^{m \times n}$ | $n$-dimensional real vector space, $m \times n$ real matrix space |
| $A, B, C$ | Matrices (uppercase letters) |
| $\mathbf{v}, \mathbf{u}, \mathbf{w}$ | Vectors (boldface lowercase letters) |
| $a, b, \lambda, \alpha$ | Scalars (lowercase or Greek letters) |
| $V, W, U$ | Vector spaces |
| $T, S$ | Linear transformations |
| $A^T, A^H$ | Transpose, conjugate transpose (Hermitian transpose) |
| $A^{-1}$ | Inverse matrix |
| $\det(A)$ | Determinant |
| $\operatorname{tr}(A)$ | Trace |
| $\operatorname{rank}(A)$ | Rank |
| $\dim(V)$ | Dimension |
| $\ker(T), \operatorname{im}(T)$ | Kernel (null space), image (range) |
| $\langle \mathbf{u}, \mathbf{v} \rangle$ | Inner product |
| $\|\mathbf{v}\|$ | Norm |
| $\sigma_i(A)$ | The $i$-th singular value |
| $\lambda_i(A)$ | The $i$-th eigenvalue |
| $I_n$ or $I$ | $n \times n$ identity matrix |
| $O$ | Zero matrix |
| $\mathbf{0}$ | Zero vector |
| $\oplus$ | Direct sum |
| $\otimes$ | Tensor product / Kronecker product |

---

## How to Use This Site

- **Search**: Use the search bar at the top of the page to quickly find any concept, theorem, or keyword. Supports both Chinese and English search.
- **Navigation**: Browse by chapter using the left sidebar, or use the right-hand table of contents to jump to a specific section on the current page.
- **Light/Dark Mode**: Click the brightness icon at the top to switch between light and dark modes.
- **Math Formulas**: All formulas on this site are rendered with MathJax. You can copy the LaTeX source by right-clicking a formula.

---

## References

The content of this site draws from the following classic textbooks:

1. **Strang, G.** *Introduction to Linear Algebra*. Wellesley-Cambridge Press.
2. **Lay, D.C.** *Linear Algebra and Its Applications*. Pearson.
3. **Axler, S.** *Linear Algebra Done Right*. Springer.
4. **Hoffman, K. & Kunze, R.** *Linear Algebra*. Prentice Hall.
5. **Horn, R.A. & Johnson, C.R.** *Matrix Analysis*. Cambridge University Press.
6. **Horn, R.A. & Johnson, C.R.** *Topics in Matrix Analysis*. Cambridge University Press.
7. **Meyer, C.D.** *Matrix Analysis and Applied Linear Algebra*. SIAM.
8. **Golub, G.H. & Van Loan, C.F.** *Matrix Computations*. Johns Hopkins University Press.
9. **Halmos, P.R.** *Finite-Dimensional Vector Spaces*. Springer.
10. **Lax, P.D.** *Linear Algebra and Its Applications*. Wiley.
