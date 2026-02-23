# Chapter 12: Jordan Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues and Multiplicity (Ch6) · Similarity (Ch5) · Minimal Polynomials (Ch13b)

**Chapter Outline**: Limitations of Diagonalization → Generalized Eigenvectors → Jordan Blocks $J_k(\lambda)$ → Jordan Decomposition Theorem → Calculation Steps (Chain Construction) → Uniqueness (Elementary Divisors) → Applications (Matrix Powers, Stability of DE Solutions)

**Extension**: Jordan Canonical Form is the most refined structure under similarity; it reveals the "break-and-repair" mechanism (the 1s on the super-diagonal) when a matrix cannot be diagonalized.

</div>

When geometric multiplicity is less than algebraic multiplicity, diagonalization fails. The **Jordan Canonical Form** is the optimal algebraic compensation for this failure. By introducing 1s above the diagonal (coupling generalized eigenvectors), it establishes a universal endpoint for similarity transformations.

---

## 12.1 Core Structure and Theorem

!!! definition "Definition 12.1 (Jordan Block)"
    A Jordan block $J_k(\lambda)$ of order $k$ is a square matrix with $\lambda$ on the diagonal and 1s on the super-diagonal:
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & & \\ & \lambda & \ddots & \\ & & \ddots & 1 \\ & & & \lambda \end{pmatrix}$$

!!! theorem "Theorem 12.1 (Jordan Canonical Form Theorem)"
    Every complex square matrix $A$ is similar to a Jordan form $J$, which is unique up to the ordering of the blocks.

---

## Exercises

1. **[Fundamentals] Write the explicit form of $J_2(3)$.**
   ??? success "Solution"
       $J_2(3) = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$.

2. **[Diagonalizability] What are the eigenvalues of $\begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$? Can it be diagonalized?**
   ??? success "Solution"
       Eigenvalue is 2 (multiplicity 2). It is a Jordan block itself. Since there is only one independent eigenvector (geometric multiplicity 1 < algebraic multiplicity 2), it cannot be diagonalized.

3. **[Block Count] What determines the total number of Jordan blocks for eigenvalue $\lambda$ in the Jordan form?**
   ??? success "Solution"
       The number of blocks equals the **geometric multiplicity** of $\lambda$, which is the dimension of the nullspace $\ker(\lambda I - A)$. Each block corresponds to an independent chain of generalized eigenvectors.

4. **[Calculation] Find the Jordan form of $A = \begin{pmatrix} 5 & 4 \\ -1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Characteristic poly: $(\lambda-5)(\lambda-1)+4 = \lambda^2-6\lambda+9 = (\lambda-3)^2$.
       Eigenvalue $\lambda=3$ (multiplicity 2).
       Geometric multiplicity: $3I-A = \begin{pmatrix} -2 & -4 \\ 1 & 2 \end{pmatrix}$, rank is 1, so nullity is 1.
       The Jordan form is $J_2(3) = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$.

5. **[Generalized Eigenvectors] Define a generalized eigenvector of order 2.**
   ??? success "Solution"
       A vector $v$ such that $(A-\lambda I)^2 v = 0$ but $(A-\lambda I) v \neq 0$. It is the end of the chain $[(A-\lambda I)v, v]$.

6. **[Matrix Power] Calculate the $n$-th power of $J = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$.**
   ??? success "Solution"
       $J^n = \begin{pmatrix} \lambda^n & n\lambda^{n-1} \\ 0 & \lambda^n \end{pmatrix}$. The upper-right term is obtained using the binomial expansion or derivative properties.

7. **[Uniqueness] Why is the Jordan form unique under similarity?**
   ??? success "Solution"
       The size and number of blocks are uniquely determined by the sequence of dimensions of kernels $\dim \ker(A-\lambda I)^k$. These dimensions are similarity invariants.

8. **[Trace and Det] If $A$ has Jordan form composed of $J_2(1)$ and $J_1(2)$, find $\det(A)$.**
   ??? success "Solution"
       Eigenvalues are $\{1, 1, 2\}$. $\det(A) = 1 \cdot 1 \cdot 2 = 2$.

9. **[Minimal Polynomial] If the largest block for $\lambda$ has size $k$, what is the power of $(\lambda-x)$ in the minimal polynomial?**
   ??? success "Solution"
       The power is exactly $k$. The minimal polynomial captures the minimum power needed to annihilate each Jordan block.

10. **[Application] How does Jordan form explain the $t e^{\lambda t}$ term in the solution to $\dot{x} = Ax$?**
    ??? success "Solution"
        When $A$ is not diagonalizable, the matrix exponential $e^{At} = P e^{Jt} P^{-1}$. The exponential of a Jordan block $J_k(\lambda)$ contains $t^m e^{\lambda t}$ terms ($m < k$), corresponding to resonance or critical damping behavior.

## Chapter Summary

Jordan form is the finale of similarity theory:

1. **Structural Completeness**: Solved the classification problem for all matrices under similarity.
2. **Chain Logic**: Smoothed out degenerate operators via chains of generalized eigenvectors.
3. **Calculus Basis**: Provided the simplest model for studying matrix functions (especially exponentials) and long-term dynamics.
