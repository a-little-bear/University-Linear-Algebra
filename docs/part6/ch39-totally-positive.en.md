# Chapter 39: Totally Positive Matrices

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Positive Definite Matrices (Ch16) · Matrix Analysis (Ch14)

**Chapter Outline**: Definition of Totally Non-negative (TN) and Totally Positive (TP) Matrices (Positivity of all Minors) → Classic Examples (Vandermonde, Cauchy, Pascal) → Core Structure: Bidiagonal Factorization → Spectral Properties: Reality, Positivity, and Distinctness of Eigenvalues → Variation-Diminishing Property (VDP) → Oscillatory Matrices → Applications: Spline Approximation, Statistics (Log-concave Densities), and Combinatorics

**Extension**: Total positivity is a far more stringent constraint than positive definiteness; it requires every "detail" (minor) of the matrix to possess positivity, revealing the operator's immense talent for preserving data monotonicity and smoothness—a bridge between continuous analysis and discrete combinatorics.

</div>

If positive definite matrices require the quadratic form to be positive, then **Totally Positive Matrices** require much more: every single minor (not just the principal ones) must be non-negative. This extreme structural constraint endows TP matrices with a set of astonishing properties, such as having strictly positive and distinct eigenvalues, and the ability to reduce the number of sign changes in a vector. This chapter explores this field where algebraic beauty meets combinatorial depth.

---

## 39.1 Definitions and Classic Examples

!!! definition "Definition 39.1 (TN and TP Matrices)"
    A matrix $A \in M_{m \times n}(\mathbb{R})$ is:
    1.  **Totally Non-negative (TN)**: if all its minors (determinants of submatrices of any order) are $\ge 0$.
    2.  **Totally Positive (TP)**: if all its minors are $> 0$.

!!! example "Example 39.1 (Typical TP Matrices)"
    1.  **Vandermonde Matrix**: With nodes $0 < x_1 < x_2 < \cdots < x_n$, its total positivity stems from the monotonicity of polynomials.
    2.  **Cauchy Matrix**: Entries $a_{ij} = 1/(x_i + y_j)$ form a TP matrix under specific node constraints.
    3.  **Pascal Matrix**: Matrices formed by binomial coefficients are totally non-negative.

---

## 39.2 Core Property: Variation-Diminishing (VDP)

!!! theorem "Theorem 39.1 (Variation-Diminishing Property)"
    Let $A$ be a TN matrix. If $\mathbf{y} = A\mathbf{x}$, then the number of sign changes in the components of $\mathbf{y}$ does not exceed the number of sign changes in $\mathbf{x}$:
    $$V(A\mathbf{x}) \le V(\mathbf{x})$$
    **Physical Meaning**: TN matrices act as powerful "smoothing" operators, suppressing oscillations in a signal.

---

## 39.3 Spectral Properties and Bidiagonal Factorization

!!! theorem "Theorem 39.2 (Spectrum of TP Matrices)"
    If $A$ is an $n \times n$ totally positive matrix, its eigenvalues $\lambda_1, \lambda_2, \ldots, \lambda_n$ satisfy:
    $$\lambda_1 > \lambda_2 > \cdots > \lambda_n > 0$$
    **Conclusion**: Eigenvalues are all **positive, real, and strictly distinct**.

!!! technique "Technique: Bidiagonal Factorization"
    Every non-singular TN matrix can be factored into a product of simple non-negative bidiagonal matrices. This proves that total positivity can be generated from local positive combinations of adjacent elements.

---

## Exercises

**1. [Basics] Determine if $\begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$ is totally positive.**

??? success "Solution"
    **Check all minors:**
    1. $1 \times 1$ minors: $1, 1, 1, 2$. All positive.
    2. $2 \times 2$ minor: $\det = 1 \cdot 2 - 1 \cdot 1 = 1 > 0$.
    **Conclusion**: Since all minors are positive, it is a Totally Positive (TP) matrix.

**2. [Counter-example] If a matrix is symmetric positive definite, is it necessarily totally positive?**

??? success "Solution"
    **Conclusion: Not necessarily.**
    **Counter-example**: $A = \begin{pmatrix} 2 & 1.5 \\ 1.5 & 2 \end{pmatrix}$.
    While PD, if we increase dimensions and include negative elements, principal minors may stay positive (maintaining PD-ness) while $1 \times 1$ non-principal minors become negative, violating TP-ness. TP-ness requires **all** minors to be positive.

**3. [VDP] What is the number of sign changes in the vector $\mathbf{x} = (1, -1, 1)^T$?**

??? success "Solution"
    **Counting:**
    1. From 1 to -1: 1st change.
    2. From -1 to 1: 2nd change.
    **Conclusion**: Number of sign changes $V(\mathbf{x}) = 2$.

**4. [Property] Prove: If $A$ is TP, then its transpose $A^T$ is also TP.**

??? success "Solution"
    **Proof:**
    1. The minors of $A^T$ are the determinants of the transposes of submatrices of $A$.
    2. Transposition does not change the determinant.
    3. Since all minors of $A$ are positive, all minors of $A^T$ must be positive as well.

**5. [Eigenvalues] If $A$ is a $3 \times 3$ TP matrix, can its eigenvalues be $\{4, 4, 1\}$?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: According to the spectral properties of TP matrices, the eigenvalues must be strictly distinct simple roots. The existence of a repeated root would imply some form of degeneracy inconsistent with the strict constraints of total positivity.

**6. [Pascal] Write the $2 \times 2$ lower triangular Pascal matrix and verify it is TN.**

??? success "Solution"
    **Construction:**
    $L = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$.
    **Minors check**: $1, 0, 1, 1 \ge 0$ and $\det=1 \ge 0$.
    **Conclusion**: It is a totally non-negative (TN) matrix.

**7. [Product] Is the product of two TN matrices always TN?**

??? success "Solution"
    **Conclusion: Yes.**
    **Proof Key**: Use the **Binet-Cauchy formula**. A $k$-th order minor of the product is a weighted combination of $k$-th order minors of the factors. Since all factor minors are non-negative, the weighted sum remains non-negative.

**8. [Oscillatory] What is an oscillatory matrix and how does it relate to TN?**

??? success "Solution"
    **Definition**: A TN matrix $A$ is oscillatory if there exists an integer $k$ such that $A^k$ is totally positive (TP).
    **Criterion**: An irreducible TN matrix with strictly positive diagonal entries is generally oscillatory.

**9. [Application] Why are TP matrices used in Splines?**

??? success "Solution"
    **Core Reason:**
    The Variation-Diminishing Property (VDP) of TP matrices ensures that the interpolating curve does not produce unnecessary oscillations between data points. It guarantees that local monotonicity of the function is perfectly preserved by the algebraic operator.

**10. [Limits] As the order $n \to \infty$, how does the condition number of a TP matrix typically behave?**

??? success "Solution"
    **Numerical Trait:**
    TP matrices are typically **extremely ill-conditioned**.
    Due to the strict decay chain of eigenvalues $\lambda_1 > \lambda_2 > \cdots$, the minimum eigenvalue tends toward zero rapidly for many classic instances (like Hilbert matrices), causing the condition number to explode.

## Chapter Summary

Totally positive matrices represent the ultimate sense of "order" among linear operators:

1.  **Global Positivity**: It requires every local determinant detail to point in the positive direction, leading to extreme stability in the global spectral structure.
2.  **Guardians of Smoothness**: The variation-diminishing property proves that TP matrices are essentially "anti-oscillatory," serving as the mathematical bedrock for spline theory and approximation.
3.  **Combinatorial Aesthetics**: Through bidiagonal factorization, total positivity simplifies complex determinant relations into adjacent element interactions, bridging the gap between continuous analysis and discrete combinatorics.
