# Chapter 15: Matrix Norms and Perturbation Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Matrix Decompositions (Ch10) · Singular Value Decomposition (Ch11)

**Chapter Outline**: Motivation for Norms (Measuring Magnitude) → Vector Norms ($L_1, L_2, L_\infty$) → Induced Matrix Norms (Operator Norms) → Frobenius Norm → Equivalence of Norms → Condition Number $\kappa(A)$ → Perturbation Analysis of Linear Systems → Eigenvalue Perturbation & Bauer-Fike Theorem → Classic Ill-conditioned Matrices (Hilbert Matrix)

**Extension**: Norms are the "rulers" used to measure the scale of mathematical objects, while perturbation theory is the science of how computational results fluctuate when inputs are disturbed by real-world noise; it is the red line of Numerical Linear Algebra (Ch22).

</div>

In pure mathematics, we speak of exact solutions, but in numerical computation, every input carries rounding errors or observation noise. **Norms** quantify the size of these errors, while **Condition Numbers** reveal the magnification factor the system applies to those errors. This chapter establishes a rigorous framework for assessing the reliability of numerical computations.

---

## 15.1 Vector and Matrix Norms

!!! definition "Definition 15.1 (Common Vector Norms)"
    For $\mathbf{x} \in \mathbb{C}^n$:
    1.  **1-norm**: $\|\mathbf{x}\|_1 = \sum |x_i|$
    2.  **2-norm** (Euclidean norm): $\|\mathbf{x}\|_2 = \sqrt{\sum |x_i|^2}$
    3.  **$\infty$-norm**: $\|\mathbf{x}\|_\infty = \max |x_i|$

!!! definition "Definition 15.2 (Induced Matrix Norms)"
    The matrix norm defined by a vector norm is called an induced norm: $\|A\| = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|}{\|\mathbf{x}\|}$.
    - **Spectral Norm**: $\|A\|_2 = \sigma_{\max}(A)$ (the largest singular value).
    - **1-norm**: Maximum absolute column sum.
    - **$\infty$-norm**: Maximum absolute row sum.

---

## 15.2 Condition Numbers and Stability

!!! definition "Definition 15.3 (Condition Number $\kappa(A)$)"
    The condition number of a square matrix $A$ is defined as:
    $$\kappa(A) = \|A\| \|A^{-1}\|$$
    $\kappa(A) \ge 1$. The larger the condition number, the more **ill-conditioned** the system, meaning small changes in input lead to large fluctuations in output.

!!! theorem "Theorem 15.1 (Perturbation Bound for Linear Systems)"
    Consider $(A+\Delta A)(x+\Delta x) = b+\Delta b$. The relative error satisfies:
    $$\frac{\|\Delta x\|}{\|x\|} \le \frac{\kappa(A)}{1 - \kappa(A) \frac{\|\Delta A\|}{\|A\|}} \left( \frac{\|\Delta A\|}{\|A\|} + \frac{\|\Delta b\|}{\|b\|} \right)$$
    This shows the condition number is the **magnification factor** for error propagation.

---

## 15.3 Eigenvalue Perturbation

!!! theorem "Theorem 15.2 (Bauer-Fike Theorem)"
    Let $A = V \Lambda V^{-1}$ be diagonalizable. If $\mu$ is an eigenvalue of $A+E$, then there exists an eigenvalue $\lambda$ of $A$ such that:
    $$|\mu - \lambda| \le \kappa_p(V) \|E\|_p$$
    **Significance**: If the diagonalizing matrix $V$ is ill-conditioned (near a non-normal matrix), eigenvalues are extremely sensitive to perturbations. For normal matrices (e.g., symmetric), $\kappa_2(V)=1$, and eigenvalues are very stable.

---

## Exercises

1. **[Calculation] Compute the $L_1, L_2, L_\infty$ norms of $\mathbf{x} = (3, -4)^T$.**
   ??? success "Solution"
       $\|\mathbf{x}\|_1 = 7, \|\mathbf{x}\|_2 = 5, \|\mathbf{x}\|_\infty = 4$.

2. **[Matrix Norm] Find the $\infty$-norm of $A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}$.**
   ??? success "Solution"
       Maximum row sum: $\max(1+2, 0+3) = 3$.

3. **[Condition] If $A = \operatorname{diag}(10, 0.1)$, compute its condition number relative to the 2-norm.**
   ??? success "Solution"
       $\|A\|_2 = 10, \|A^{-1}\|_2 = 10 \implies \kappa_2(A) = 100$.

4. **[Properties] Prove the induced matrix norm satisfies sub-multiplicativity: $\|AB\| \le \|A\| \|B\|$.**
   ??? success "Solution"
       $\|AB\mathbf{x}\| \le \|A\|\|B\mathbf{x}\| \le \|A\|\|B\|\|\mathbf{x}\|$. The definition follows.

5. **[Error] If $\kappa(A)=10^4$ and input error is $10^{-6}$, what is the approximate upper bound for the relative error of the solution?**
   ??? success "Solution"
       Approximately $10^4 \cdot 10^{-6} = 10^{-2}$.

6. **[Unitary] Prove the 2-condition number of an orthogonal matrix is always 1.**
   ??? success "Solution"
       Since orthogonal matrices preserve length, $\|Q\|_2 = 1$ and $\|Q^{-1}\|_2 = \|Q^T\|_2 = 1$.

7. **[Bauer-Fike] Why are eigenvalues of symmetric matrices more stable than those of non-normal matrices?**
   ??? success "Solution"
       Symmetric matrices are unitarily diagonalizable, meaning $V$ is orthogonal and its condition number is 1, so the error is not magnified.

8. **[Frobenius] Calculate the Frobenius norm of $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $\sqrt{1^2+1^2+1^2+1^2} = 2$.

9. **[Ill-conditioned] Name a well-known ill-conditioned matrix.**
   ??? success "Solution"
       The Hilbert matrix $H_{ij} = \frac{1}{i+j-1}$.

10. **[Application] Why are unitary transformations (like Householder) favored in numerical computation?**
    ??? success "Solution"
        Because their condition number is always 1, they do not magnify rounding errors, preserving the stability of the algorithm.

## Chapter Summary

Norms and perturbation theory are the "lifeline" of computational mathematics:

1.  **Relative Magnitude**: Norms transform abstract matrices into comparable numerical values, establishing a metric for error assessment.
2.  **Stability Indicator**: The condition number is the ultimate barometer of algorithmic reliability, marking the boundaries of numerical simulation—not every mathematically correct problem is solvable on a computer.
3.  **Cost of Diagonalization**: The Bauer-Fike theorem reveals that the stability of the diagonalization process itself is limited by the degree of orthogonality of the basis vectors, emphasizing the central role of orthonormal bases in numerical linear algebra.
