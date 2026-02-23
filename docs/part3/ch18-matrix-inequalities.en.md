# Chapter 18: Matrix Inequalities

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Norms (Ch15) · Singular Value Decomposition (Ch11)

**Chapter Outline**: From Scalar to Operator Inequalities → Eigenvalue Inequalities (Weyl's Inequalities, Interlacing Theorems) → Determinant Inequalities (Hadamard, Fischer Inequalities) → Trace Inequalities (von Neumann, Golden-Thompson) → Singular Value Inequalities (Ky Fan Norms) → Majorization Theory ($\prec$) → Introduction to Operator Monotone and Convex Functions (Ch46)

**Extension**: Matrix inequalities are the mathematical pillars of Information Theory (concavity of entropy), Compressed Sensing, and the Uncertainty Principle in Quantum Mechanics; they transform exact "equalities" into restricted "inclusions" or "bounds."

</div>

Matrix inequalities are among the most sophisticated branches of matrix analysis. They study the constraints and trade-offs between matrix properties (such as eigenvalues, singular values, and traces) rather than exact values. Just as inequalities on the real line characterize relative magnitudes, matrix inequalities characterize the relative distribution of energy and information within operators.

---

## 18.1 Eigenvalue Inequalities

!!! theorem "Theorem 18.1 (Weyl's Inequalities)"
    Let $A, B$ be Hermitian matrices, and let $C = A + B$. Arrange their eigenvalues in descending order. For all $j+k-1 \le n$:
    $$\lambda_{j+k-1}(A+B) \le \lambda_j(A) + \lambda_k(B)$$
    **Physical Meaning**: The impact of a perturbation on a system's eigenvalues (energy levels) is strictly limited by the scale of the perturbation matrix's spectrum.

!!! theorem "Theorem 18.2 (Cauchy Interlacing Theorem)"
    Let $B$ be an $(n-1) \times (n-1)$ principal submatrix of an $n \times n$ Hermitian matrix $A$. Then:
    $$\lambda_1(A) \ge \lambda_1(B) \ge \lambda_2(A) \ge \lambda_2(B) \ge \cdots \ge \lambda_{n-1}(B) \ge \lambda_n(A)$$

---

## 18.2 Determinant and Trace Inequalities

!!! theorem "Theorem 18.3 (Hadamard's Inequality)"
    For any positive definite matrix $A \succ 0$:
    $$\det(A) \le \prod_{i=1}^n a_{ii}$$
    Equality holds if and only if $A$ is diagonal.
    **Geometric Interpretation**: The volume of a parallelotope is less than or equal to the product of its side lengths (attaining equality only when sides are orthogonal).

!!! theorem "Theorem 18.4 (Golden-Thompson Inequality)"
    For Hermitian matrices $A, B$:
    $$\operatorname{tr}(e^{A+B}) \le \operatorname{tr}(e^A e^B)$$
    This is a central inequality in quantum statistical mechanics.

---

## 18.3 Majorization Theory ($\prec$)

!!! definition "Definition 18.1 (Majorization)"
    Let $x, y \in \mathbb{R}^n$. $y$ **majorizes** $x$ (written $x \prec y$) if $\sum_{i=1}^k x_{(i)} \le \sum_{i=1}^k y_{(i)}$ for $k=1,\ldots,n-1$ and the totals are equal.
    **Schur-Horn Theorem**: The vector of diagonal entries of a Hermitian matrix is majorized by the vector of its eigenvalues: $\operatorname{diag}(A) \prec \lambda(A)$.

---

## Exercises

1. **[Weyl] Given $\|E\|_2 = 0.1$, if an eigenvalue of $A$ is 5, in what interval must the corresponding eigenvalue of $A+E$ lie?**
   ??? success "Solution"
       By Weyl's inequality, $|\lambda_i(A+E) - \lambda_i(A)| \le \|E\|_2$. Thus, it lies in $[4.9, 5.1]$.

2. **[Hadamard] Calculate the determinant of $\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ and verify Hadamard's inequality.**
   ??? success "Solution"
       $\det = 3$. Product of diagonals $2 \cdot 2 = 4$. $3 \le 4$ is verified.

3. **[Interlacing] If a $3 \times 3$ matrix has eigenvalues 10, 5, 1, can its $2 \times 2$ principal submatrix have a maximum eigenvalue of 12?**
   ??? success "Solution"
       No. By the interlacing theorem, $\lambda_1(B) \le \lambda_1(A) = 10$.

4. **[Trace] Prove for positive definite $A, B$ that $\operatorname{tr}(AB) \le \operatorname{tr}(A)\operatorname{tr}(B)$.**
   ??? success "Solution"
       $\operatorname{tr}(AB) = \sum \lambda_i(AB) \le \sum \sigma_i(A)\sigma_i(B) \le (\sum \sigma_i(A))(\sum \sigma_i(B)) = \operatorname{tr}(A)\operatorname{tr}(B)$.

5. **[Majorization] Determine the majorization relationship between $(1, 1)$ and $(2, 0)$.**
   ??? success "Solution"
       $(1, 1) \prec (2, 0)$ because $1 < 2$ and $1+1 = 2+0$.

6. **[Fischer] State Fischer's Inequality.**
   ??? success "Solution"
       For a partitioned positive definite matrix $\begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$, $\det \begin{pmatrix} A & B \\ B^T & C \end{pmatrix} \le \det(A)\det(C)$.

7. **[Ky Fan] What is the Ky Fan $k$-norm?**
   ??? success "Solution"
       The sum of the $k$ largest singular values: $\|A\|_{(k)} = \sum_{i=1}^k \sigma_i(A)$.

8. **[Arithmetic-Geometric] Prove $\det(A)^{1/n} \le \frac{1}{n} \operatorname{tr}(A)$ for $A \succ 0$.**
   ??? success "Solution"
       This is the arithmetic-geometric mean inequality applied to eigenvalues: $(\prod \lambda_i)^{1/n} \le \frac{1}{n} \sum \lambda_i$.

9. **[Concavity] Prove the mapping $A \mapsto \log \det A$ is concave on the positive definite cone.**
   ??? success "Solution"
       This is equivalent to verifying $\det(\lambda A + (1-\lambda)B) \ge (\det A)^\lambda (\det B)^{1-\lambda}$, which is the matrix version of the Brunn-Minkowski inequality.

10. **[Application] How are matrix inequalities used in quantum information?**
    ??? success "Solution"
        They are used to prove the Strong Subadditivity of quantum entropy, establishing upper limits on quantum communication capacity.

## Chapter Summary

Matrix inequalities define the "energy boundaries" of linear systems:

1.  **Spectral Stability**: Weyl's inequalities and interlacing theorems prove that matrix eigenvalues possess significant geometric inertia; small structural changes can only cause controlled spectral shifts.
2.  **Informational Extrema**: Hadamard and trace inequalities reveal the patterns of information loss in non-diagonalized (coupled) states, providing algebraic upper bounds for entropy estimation in information theory.
3.  **Quantification of Distribution**: Majorization theory provides a powerful tool for comparing the "dispersion" of vectors, revealing the deep containment relationship between a matrix's diagonal entries and its spectrum.
