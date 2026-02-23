# Chapter 39: Totally Positive Matrices

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Positive Definite Matrices (Ch16) · Matrix Analysis (Ch14)

**Chapter Outline**: Definition of Totally Non-negative (TN) and Totally Positive (TP) Matrices → Key Examples (Vandermonde, Cauchy, Pascal) → Bidiagonal Factorization → Spectral Properties: Positivity, Reality, and Distinctness of Eigenvalues → Variation-Diminishing Property (VDP) → Oscillatory Matrices → Applications: Splines, Statistics (Log-concave Densities), and Combinatorics

**Extension**: Total positivity is a far more restrictive condition than positive definiteness; it requires every "detail" (minor) of the matrix to be positive, revealing the operator's talent for preserving data monotonicity and smoothness.

</div>

If positive definite matrices require the quadratic form to be positive, then **Totally Positive Matrices** require much more: every single minor (not just the principal ones) must be non-negative. This extreme structural constraint endows TP matrices with a set of astonishing properties, such as having strictly positive and distinct eigenvalues, and the ability to reduce the number of sign changes in a vector. This chapter dives into this field where algebraic beauty meets combinatorial depth.

---

## 39.1 Definitions and Examples

!!! definition "Definition 39.1 (TN and TP Matrices)"
    A matrix $A \in M_{m \times n}(\mathbb{R})$ is:
    1.  **Totally Non-negative (TN)**: if all its minors (determinants of submatrices of any order) are $\ge 0$.
    2.  **Totally Positive (TP)**: if all its minors are $> 0$.

!!! example "Example 39.1 (Classic Instances)"
    1.  **Vandermonde Matrix**: With nodes $0 < x_1 < x_2 < \cdots < x_n$, total positivity arises from the monotonicity of polynomials.
    2.  **Cauchy Matrix**: $a_{ij} = 1/(x_i + y_j)$ is TP under appropriate parameter constraints.
    3.  **Pascal Matrix**: The matrix formed by binomial coefficients $\binom{i+j}{i}$ is TN.

---

## 39.2 Bidiagonal Factorization

!!! theorem "Theorem 39.1 (Bidiagonal Factorization)"
    Every non-singular totally non-negative matrix $A$ can be uniquely factored into a product of elementary non-negative lower triangular bidiagonal matrices and upper triangular bidiagonal matrices.
    **Significance**: This structure reveals that total positivity can be generated through extremely simple local positive operations (linear combinations of adjacent rows/columns).

---

## 39.3 Spectral Properties and Sign Reduction

!!! theorem "Theorem 39.2 (Spectrum of TP Matrices)"
    If $A$ is an $n \times n$ totally positive matrix, its eigenvalues $\lambda_1, \lambda_2, \ldots, \lambda_n$ satisfy:
    $$\lambda_1 > \lambda_2 > \cdots > \lambda_n > 0$$
    The eigenvalues are all **positive, real, and distinct**.

!!! theorem "Theorem 39.3 (Variation-Diminishing Property - VDP)"
    Let $A$ be a TN matrix. If $\mathbf{y} = A\mathbf{x}$, the number of sign changes in the components of $\mathbf{y}$ does not exceed the number of sign changes in the components of $\mathbf{x}$.
    **Application**: This is the fundamental mathematical reason why spline approximations preserve the smoothness of functions.

---

## Exercises

1.  **[Criteria] Determine if $\begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$ is totally positive.**
    ??? success "Solution"
        Minors: $1, 1, 1, 2$ (all positive); $\det = 1 > 0$. Since all minors are positive, it is TP.

2.  **[Counter-example] Is $\begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$ TP?**
    ??? success "Solution"
        Yes. All entries are positive and the determinant is 1. All minors are positive.

3.  **[Minors] How many minors does a $3 \times 3$ matrix have in total?**
    ??? success "Solution"
        9 minors of size $1 \times 1$, 9 minors of size $2 \times 2$, and 1 minor of size $3 \times 3$. Total = 19.

4.  **[Property] Prove: If $A$ is TP, then $A^T$ is also TP.**
    ??? success "Solution"
        The determinant of a submatrix is invariant under transposition, so the positivity of all minors is preserved.

5.  **[Eigenvalues] If $A$ is a $4 \times 4$ TP matrix, can its eigenvalues be $5, 5, 2, 1$?**
    ??? success "Solution"
        No. The eigenvalues of a TP matrix must be strictly distinct.

6.  **[Pascal] Write the $2 \times 2$ Pascal matrix and verify it is TN.**
    ??? success "Solution"
        $\begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$. All minors $\ge 0$.

7.  **[VDP] What is the number of sign changes in the vector $(1, -1, 1)$?**
    ??? success "Solution"
        2 changes (from 1 to -1, then from -1 to 1).

8.  **[Oscillatory] What is an oscillatory matrix?**
    ??? success "Solution"
        A TN matrix $A$ such that some power $A^k$ is TP.

9.  **[Product] Is the product of two TN matrices always TN?**
    ??? success "Solution"
        Yes. By the Cauchy-Binet formula, the minors of the product are weighted sums of the minors of the factors, so positivity is maintained.

10. **[Statistics] Why is total positivity related to stochastic dominance?**

   ??? success "Solution"
        TP ensures that distribution functions maintain monotonicity and extrema stability under transformations, serving as a core tool in multivariate statistical ordering theory.

## Chapter Summary

Totally positive matrices represent the ultimate "order" among linear operators:

1.  **Global Positivity**: Every local determinant detail is positive, leading to extreme stability in the global spectral structure (distinct, positive real eigenvalues).
2.  **Guardian of Smoothness**: The variation-diminishing property proves that TN matrices are essentially "anti-oscillatory," suppressing noise and fluctuations during signal propagation—the cornerstone of spline theory.
3.  **Combinatorial Aesthetics**: Through bidiagonal factorization, total positivity simplifies complex determinant relations into adjacent element interactions, bridging continuous analysis and discrete combinatorics.
