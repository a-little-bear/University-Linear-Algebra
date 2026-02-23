# Chapter 18: Matrix Inequalities

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Norms (Ch15) · SVD (Ch11) · Eigenvalues (Ch06)

**Chapter Outline**: From Scalar to Operator Inequalities → Eigenvalue Inequalities (Weyl’s Inequalities, Cauchy Interlacing Theorem) → Determinant Inequalities (Hadamard, Fischer Inequalities) → Trace Inequalities (von Neumann, Golden-Thompson) → Singular Value Inequalities (Ky Fan Norms) → Majorization Theory ($\prec$) and its Matrix Connections → Applications: Information Theory (Concavity of Entropy), Uncertainty Relations in Physics, and Bound Estimation in Optimization

**Extension**: Matrix inequalities are the language for describing the "boundaries" of a system; they transform exact equalities into restricted inclusions, proving that even under perturbation, the core energy distribution (spectrum) of an operator maintains strong geometric inertia. They are the mathematical pillars of quantum information and modern statistics.

</div>

Matrix inequalities are among the most elegant branches of matrix analysis. They study the constraints and trade-offs between matrix properties (such as eigenvalues, traces, and determinants) rather than specific numerical values. Just as inequalities on the real line characterize relative magnitudes, matrix inequalities characterize the relative distribution of energy and information within operators. This chapter establishes a series of classical criteria for describing the "magnitude" of global matrix properties.

---

## 18.1 Eigenvalue Inequalities

!!! theorem "Theorem 18.1 (Weyl’s Inequality)"
    Let $A, B$ be $n \times n$ Hermitian matrices, and let $C = A + B$. Arrange their eigenvalues in descending order $\lambda_1 \ge \lambda_2 \ge \cdots$. Then for any $1 \le j, k \le n$ such that $j+k-1 \le n$:
    $$\lambda_{j+k-1}(A+B) \le \lambda_j(A) + \lambda_k(B)$$
    **Physical Meaning**: This indicates that the impact of a perturbation on a system's energy levels (eigenvalues) is dually constrained by the strength of the perturbation and the spectral structure of the original system.

!!! theorem "Theorem 18.2 (Cauchy Interlacing Theorem)"
    Let $B$ be an $(n-1) \times (n-1)$ principal submatrix of an $n \times n$ Hermitian matrix $A$. Then:
    $$\lambda_1(A) \ge \lambda_1(B) \ge \lambda_2(A) \ge \lambda_2(B) \ge \cdots \ge \lambda_{n-1}(B) \ge \lambda_n(A)$$
    **Geometric Insight**: The compression of space leads to a "contraction" of the spectrum; the eigenvalues of the submatrix are strictly bracketed by those of the original matrix.

---

## 18.2 Determinant and Trace Inequalities

!!! theorem "Theorem 18.3 (Hadamard’s Inequality)"
    For any positive definite matrix $A \in M_n(\mathbb{C})$:
    $$\det(A) \le \prod_{i=1}^n a_{ii}$$
    Equality holds iff $A$ is a diagonal matrix.
    **Geometric Meaning**: The volume of a parallelotope is less than or equal to the product of its side lengths (maximum volume is attained only when the sides are orthogonal).

---

## 18.3 Majorization Theory ($\prec$)

!!! definition "Definition 18.1 (Majorization)"
    Let $x, y \in \mathbb{R}^n$, with components sorted in non-increasing order. We say $y$ **majorizes** $x$ (written $x \prec y$) if:
    1.  For $k=1, \ldots, n-1$, $\sum_{i=1}^k x_{[i]} \le \sum_{i=1}^k y_{[i]}$.
    2.  The total sums are equal: $\sum x_i = \sum y_i$.

!!! theorem "Theorem 18.4 (Schur-Horn Theorem)"
    For any Hermitian matrix $A$, the vector of its diagonal entries $\mathbf{d}$ is majorized by the vector of its eigenvalues $\boldsymbol{\lambda}$:
    $$\mathbf{d} \prec \boldsymbol{\lambda}$$
    This reveals a deep containment relationship between internal entry distribution and external spectral performance.

---

## Exercises

**1. [Weyl] If $\|E\|_2 = 0.1$ and the maximum eigenvalue of $A$ is 5, find the range of the maximum eigenvalue of $A+E$.**

??? success "Solution"
    **Calculation:**
    According to a corollary of Weyl’s inequality: $|\lambda_i(A+E) - \lambda_i(A)| \le \|E\|_2$.
    1. $\lambda_{\max}(A) = 5$.
    2. The maximum shift is 0.1.
    **Conclusion**: $\lambda_{\max}(A+E)$ lies in the interval $[4.9, 5.1]$. This reflects the Lipschitz continuity of eigenvalues under perturbation.

**2. [Hadamard] Calculate the determinant of $\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ and verify Hadamard's inequality.**

??? success "Solution"
    **Steps:**
    1. Determinant: $\det = 2\cdot 2 - 1\cdot 1 = 3$.
    2. Product of diagonals: $a_{11}a_{22} = 2\cdot 2 = 4$.
    3. Check inequality: $3 \le 4$.
    **Conclusion**: Verification successful. Since off-diagonal entries are non-zero (basis vectors are not orthogonal), the volume is reduced.

**3. [Interlacing] A $3 \times 3$ Hermitian matrix has eigenvalues $10, 5, 1$. Can the maximum eigenvalue of its $2 \times 2$ principal submatrix be 12? Can it be 0.5?**

??? success "Solution"
    **Determination:**
    1. By the interlacing theorem, the submatrix's largest eigenvalue $\mu_1$ must satisfy $\lambda_2(A) \le \mu_1 \le \lambda_1(A)$.
    2. Here, $5 \le \mu_1 \le 10$.
    **Conclusion**: It cannot be 12 (exceeds the upper bound) or 0.5 (falls below the lower bound). A subspace cannot jump outside the spectral span of the original space.

**4. [Trace] Prove that for positive definite matrices $A, B$, $\operatorname{tr}(AB) \ge 0$.**

??? success "Solution"
    **Proof:**
    1. Use the cyclic property: $\operatorname{tr}(AB) = \operatorname{tr}(A^{1/2} A^{1/2} B) = \operatorname{tr}(A^{1/2} B A^{1/2})$.
    2. Let $C = A^{1/2} B A^{1/2}$. Since $B \succ 0$, $C$ must also be positive definite (congruence transformation).
    3. The diagonal entries of a PD matrix are positive, so its trace must be positive.
    **Conclusion**: $\operatorname{tr}(AB) \ge 0$. This shows the product of positive operators preserves positivity in an average sense.

**5. [Majorization] Determine the majorization relationship between $(1, 1, 1)$ and $(3, 0, 0)$.**

??? success "Solution"
    **Analysis:**
    1. Total sum: $1+1+1=3$ and $3+0+0=3$. Equal.
    2. Prefix sums:
       - $k=1$: $1 < 3$. (Satisfied)
       - $k=2$: $1+1 < 3+0$. (Satisfied)
    **Conclusion**: $(1, 1, 1) \prec (3, 0, 0)$. Intuitively, uniform distributions are majorized by extreme ones.

**6. [Fischer] What is Fischer's Inequality and how does it generalize Hadamard's?**

??? success "Solution"
    **Definition:**
    For a partitioned positive definite matrix $M = \begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$, we have $\det(M) \le \det(A)\det(C)$.
    **Generalization**: Hadamard's inequality is the special case where each block is $1 \times 1$. It states that partitioning a matrix into independent blocks increases its generalized variance (volume).

**7. [Ky Fan] What is the Ky Fan $k$-norm and its relation to majorization?**

??? success "Solution"
    **Definition**: $\|A\|_{(k)} = \sum_{i=1}^k \sigma_i(A)$ (sum of the $k$ largest singular values).
    **Relation**: For singular value vectors, $x \prec y$ is equivalent to $\|x\|_{(k)} \le \|y\|_{(k)}$ for all $k$. This establishes majorization as the core for defining operator norms.

**8. [Arithmetic-Geometric] Prove $\det(A)^{1/n} \le \frac{1}{n} \operatorname{tr}(A)$ for $A \succ 0$.**

??? success "Solution"
    **Proof:**
    1. Let $\lambda_1, \ldots, \lambda_n$ be the eigenvalues of $A$.
    2. $\det(A) = \prod \lambda_i$.
    3. $\operatorname{tr}(A) = \sum \lambda_i$.
    4. Apply the classic **Arithmetic-Geometric Mean Inequality**: $(\prod \lambda_i)^{1/n} \le \frac{1}{n} \sum \lambda_i$.
    **Conclusion**: The volume of a matrix is strictly limited by its average energy.

**9. [Concavity] Is the mapping $A \mapsto \log \det A$ convex or concave on the set of positive definite matrices?**

??? success "Solution"
    **Conclusion**: It is a **concave function**.
    **Significance**: This corresponds to the entropy property in information theory. The mixture (convex combination) of two systems always contains at least as much information (logarithm of volume) as the weighted sum of their individual information contents.

**10. [Application] How do matrix inequalities manifest the Uncertainty Principle in quantum mechanics?**

??? success "Solution"
    **Explanation:**
    The Heisenberg Uncertainty Principle is algebraically expressed as a lower bound on the product of variances of two observable operators $A, B$:
    $\operatorname{Var}(A)\operatorname{Var}(B) \ge \frac{1}{4} |\langle [A, B] \rangle|^2$.
    This is a direct consequence of the commutator inequality in inner product spaces, proving that non-commutativity leads to mutually exclusive measurement boundaries.

## Chapter Summary

Matrix inequalities define the "energy boundaries" of linear systems:

1.  **Spectral Stability**: Weyl’s inequality and interlacing theorems prove that matrix eigenvalues possess extreme geometric inertia; small structural changes can only cause controlled spectral drifts.
2.  **Informational Extrema**: Hadamard and trace inequalities reveal the patterns of information loss in non-diagonalized (coupled) states, providing algebraic upper bounds for entropy estimation in information theory.
3.  **Quantification of Distribution**: Majorization theory provides a powerful tool for comparing the "dispersion" of vectors, revealing the deep containment relationship between a matrix's diagonal elements and its spectrum—the cornerstone of modern compressed sensing and statistical physics.
