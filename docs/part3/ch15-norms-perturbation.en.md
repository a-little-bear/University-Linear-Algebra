# Chapter 15: Norms and Perturbation Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Matrix Decompositions (Ch10) · Singular Value Decomposition (Ch11)

**Chapter Outline**: Motivation for Norms (Measuring Magnitude) → Vector Norms ($L_1, L_2, L_\infty, L_p$) → Induced Matrix Norms (Operator Norms) → Frobenius and Schatten Norms → Norm Equivalence Theorem → Algebraic and Geometric Meanings of Condition Number $\kappa(A)$ → Perturbation Analysis of Linear Systems → Eigenvalue Perturbation and the Bauer-Fike Theorem → Classic Ill-conditioned Matrices (The Hilbert Matrix) → Applications: Numerical Stability Assessment and Regularization Methods

**Extension**: Norms are the "rulers" for measuring the scale of mathematical objects, while perturbation theory is the science of investigating how fluctuations in input (noise) propagate to output; it is the bottom line of Numerical Linear Algebra (Ch22).

</div>

In pure mathematics, we often speak of exact solutions. In numerical computation, however, every input carries rounding errors or observation noise. **Norms** quantify the magnitude of these errors, while the **Condition Number** reveals the magnification factor of the system on those errors. This chapter establishes a rigorous framework for assessing the reliability of numerical computations, explaining why some problems, though solvable in theory, are impossible to solve on a computer.

---

## 15.1 Vector and Matrix Norms

!!! definition "Definition 15.1 (Common Vector Norms)"
    For $\mathbf{x} \in \mathbb{C}^n$:
    1.  **1-norm**: $\|\mathbf{x}\|_1 = \sum_{i=1}^n |x_i|$
    2.  **2-norm** (Euclidean norm): $\|\mathbf{x}\|_2 = \sqrt{\sum_{i=1}^n |x_i|^2}$
    3.  **$\infty$-norm**: $\|\mathbf{x}\|_\infty = \max_{i} |x_i|$

!!! definition "Definition 15.2 (Induced Matrix Norms)"
    The matrix norm defined by a vector norm is the induced norm: $\|A\| = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|}$.
    - **Spectral Norm**: $\|A\|_2 = \sigma_{\max}(A)$ (Maximum singular value).
    - **1-norm**: $\|A\|_1 = \max_{j} \sum_{i} |a_{ij}|$ (Maximum absolute column sum).
    - **$\infty$-norm**: $\|A\|_\infty = \max_{i} \sum_{j} |a_{ij}|$ (Maximum absolute row sum).

---

## 15.2 Condition Numbers and Stability

!!! definition "Definition 15.3 (Condition Number $\kappa(A)$)"
    The condition number of a square matrix $A$ is:
    $$\kappa(A) = \|A\| \|A^{-1}\|$$
    $\kappa(A) \ge 1$ always holds. The larger the condition number, the more **ill-conditioned** the system, meaning tiny input changes result in massive output fluctuations.

!!! theorem "Theorem 15.1 (Perturbation Bound for Linear Systems)"
    Consider $(A+\Delta A)(x+\Delta x) = b+\Delta b$. The relative error of the solution satisfies:
    $$\frac{\|\Delta x\|}{\|x\|} \le \frac{\kappa(A)}{1 - \kappa(A) \frac{\|\Delta A\|}{\|A\|}} \left( \frac{\|\Delta A\|}{\|A\|} + \frac{\|\Delta b\|}{\|b\|} \right)$$
    This indicates that the condition number is the **amplification factor** for error propagation.

---

## 15.3 Eigenvalue Perturbation

!!! theorem "Theorem 15.2 (The Bauer-Fike Theorem)"
    Let $A = V \Lambda V^{-1}$ be diagonalizable. If $\mu$ is an eigenvalue of $A+E$, there exists an eigenvalue $\lambda$ of $A$ such that:
    $$| \mu - \lambda | \le \kappa_p(V) \|E\|_p$$
    **Significance**: If the condition number of the diagonalizing matrix $V$ is large (i.e., the eigenvector basis is near linear dependence), the eigenvalues are extremely sensitive to perturbations. For normal matrices, $V$ is unitary and $\kappa_2(V)=1$, making eigenvalues very stable.

---

## Exercises

**1. [Calculation] Compute the $L_1, L_2, L_\infty$ norms of $\mathbf{x} = (3, -4)^T$.**

??? success "Solution"
    **Steps:**
    1. $\|\mathbf{x}\|_1 = |3| + |-4| = 7$.
    2. $\|\mathbf{x}\|_2 = \sqrt{3^2 + (-4)^2} = 5$.
    3. $\|\mathbf{x}\|_\infty = \max(3, 4) = 4$.
    **Geometric Insight**: Different norms define differently shaped "unit balls" (squares, circles, diamonds).

**2. [Matrix Norm] Find the $\infty$-norm of $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. By definition, the $\infty$-norm is the maximum absolute row sum.
    2. Row 1 sum: $|1| + |2| = 3$.
    3. Row 2 sum: $|0| + |3| = 3$.
    **Conclusion**: $\|A\|_\infty = 3$.

**3. [Conditioning] If $A = \operatorname{diag}(10, 0.1)$, compute its condition number relative to the 2-norm.**

??? success "Solution"
    **Steps:**
    1. $\|A\|_2 = \sigma_{\max} = 10$.
    2. $A^{-1} = \operatorname{diag}(0.1, 10) \implies \|A^{-1}\|_2 = 10$.
    3. $\kappa_2(A) = 10 \cdot 10 = 100$.
    **Assessment**: A condition number of 100 means input errors can be amplified 100 times. For 64-bit double precision (~16 digits), we lose about 2 digits of accuracy.

**4. [Property] Prove the induced matrix norm satisfies sub-multiplicativity: $\|AB\| \le \|A\| \|B\|$.**

??? success "Solution"
    **Proof:**
    1. By definition, $\|AB\| = \max_{\|x\|=1} \|ABx\|$.
    2. Use the vector norm inequality: $\|ABx\| \le \|A\| \|Bx\|$.
    3. Use the definition again: $\|Bx\| \le \|B\| \|x\|$.
    4. Combine: $\|ABx\| \le \|A\|\|B\|\|x\|$.
    5. Since $\|x\|=1$, the inequality holds after taking the supremum.

**5. [Error Analysis] If $\kappa(A)=10^4$ and the relative error in input $b$ is $10^{-6}$, what is the upper bound for the relative error in $x$?**

??? success "Solution"
    **Calculation:**
    From the perturbation theorem: $\frac{\|\Delta x\|}{\|x\|} \le \kappa(A) \frac{\|\Delta b\|}{\|b\|}$.
    $= 10^4 \cdot 10^{-6} = 10^{-2}$.
    **Significance**: Even if the input error is only one-in-a-million, the solution can deviate by 1%.

**6. [Unitary] Prove that the 2-condition number of an orthogonal matrix is always 1.**

??? success "Solution"
    **Proof:**
    1. An orthogonal matrix $Q$ satisfies $Q^T Q = I$.
    2. All its singular values are exactly 1. Thus $\|Q\|_2 = 1$.
    3. $Q^{-1} = Q^T$ is also orthogonal, so $\|Q^{-1}\|_2 = 1$.
    4. $\kappa_2(Q) = 1 \cdot 1 = 1$.
    **Numerical Fact**: Orthogonal transformations are the safest operations in numerical computing as they do not amplify errors at all.

**7. [Bauer-Fike] Why are eigenvalues of symmetric matrices more stable than those of highly non-normal matrices?**

??? success "Solution"
    **Analysis:**
    1. Symmetric matrices are unitarily diagonalizable: $A = Q\Lambda Q^T$.
    2. In the Bauer-Fike inequality, $V=Q$, so $\kappa_2(V) = 1$.
    3. The error bound is simply $\|E\|_2$.
    4. For non-symmetric (especially near-defective) matrices, the columns of $V$ are near linear dependence, leading to a massive $\kappa(V)$, which amplifies the spectral shift.

**8. [Frobenius] Calculate the Frobenius norm of $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**

??? success "Solution"
    **Formula:**
    $\|A\|_F = \sqrt{\sum a_{ij}^2} = \sqrt{\operatorname{tr}(A^* A)}$.
    **Calculation:**
    $\sqrt{1^2 + 1^2 + 1^2 + 1^2} = \sqrt{4} = 2$.
    Note: The spectral norm $\|A\|_2$ is also 2 here (eigenvalues are 2 and 0).

**9. [Ill-conditioning] Name a famous ill-conditioned matrix and explain why it is so.**

??? success "Solution"
    **Example: The Hilbert matrix** $H_{ij} = \frac{1}{i+j-1}$.
    **Reason**: The row vectors of this matrix are extremely "close" (almost parallel). As the order $n$ increases, the minimum singular value decreases exponentially, causing the condition number to explode (e.g., for $n=10$, $\kappa \approx 10^{13}$).

**10. [Application] Why use norm regularization (like Ridge regression) in optimization?**

??? success "Solution"
    **Algebraic Reason:**
    1. In least squares, if $A^T A$ is near singular (ill-conditioned), the solution $x$ can become massive and sensitive to noise.
    2. Adding the term $\lambda \|x\|_2^2$ is equivalent to adding a small positive $\lambda$ to the diagonal of $A^T A$.
    3. This effectively lowers the condition number of the matrix, smoothing the output and preventing overfitting.

## Chapter Summary

Norms and perturbation theory are the "lifeline" of computational mathematics:

1.  **Relativity of Magnitude**: Norms transform abstract operators into comparable numbers, establishing a measure for error assessment.
2.  **Stability Barometer**: The condition number is the unique indicator of algorithmic reliability, defining the boundaries of simulation—not every mathematically correct problem is solvable on a computer.
3.  **Cost of Diagonalization**: The Bauer-Fike theorem reveals that spectral stability depends on the orthogonality of the eigenvector basis, emphasizing the central role of orthonormal bases (unitary transforms) in numerical linear algebra.
