# Chapter 12: Jordan Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Decompositions (Ch10) · Polynomial Algebra (Ch00)

**Chapter Outline**: Limitations of Diagonalization (Defective Matrices) → Definition of Jordan Blocks → Generalized Eigenvectors and Jordan Chains → Existence and Uniqueness of Jordan Canonical Form (JCF) → Relationship between Minimal and Characteristic Polynomials → Steps to Determine JCF (Rank Methods, Weyr Characteristic) → Jordan Analysis of Matrix Powers and Series → Numerical Instability

**Extension**: The Jordan Canonical Form is the ultimate representation under similarity transformations; it is indispensable for the theoretical analysis of linear differential equations (Ch26) and matrix functions (Ch13).

</div>

Not all square matrices can be diagonalized. When the geometric multiplicity of an eigenvalue is less than its algebraic multiplicity, the matrix is called "defective." The Jordan Canonical Form (JCF) provides the structure closest to diagonal for such matrices. It not only reveals the deep structure of linear operators but also serves as the most important theoretical tool in matrix analysis.

---

## 12.1 Jordan Blocks and Generalized Eigenvectors

!!! definition "Definition 12.1 (Jordan Block)"
    A **Jordan block** $J_k(\lambda)$ of order $k$ is a square matrix with $\lambda$ on the main diagonal, 1s on the super-diagonal, and 0s elsewhere:
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & 0 & \lambda & 1 \\ 0 & 0 & 0 & 0 & \lambda \end{pmatrix}$$

!!! definition "Definition 12.2 (Generalized Eigenvector)"
    A vector $\mathbf{v}$ is a **generalized eigenvector of rank $k$** corresponding to $\lambda$ if $(A - \lambda I)^k \mathbf{v} = \mathbf{0}$ but $(A - \lambda I)^{k-1} \mathbf{v} \neq \mathbf{0}$.

---

## 12.2 The Jordan Canonical Form Theorem

!!! theorem "Theorem 12.1 (Jordan Canonical Form Theorem)"
    Every complex square matrix $A$ is similar to a Jordan Canonical Form $J$, and $J$ is unique up to the permutation of its blocks:
    $$P^{-1} A P = J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_m}(\lambda_m))$$
    - Each Jordan block corresponds to one linearly independent eigenvector (geometric multiplicity).
    - The sum of the orders of all Jordan blocks for a given eigenvalue equals its algebraic multiplicity.

---

## 12.3 The Minimal Polynomial

!!! definition "Definition 12.3 (Minimal Polynomial)"
    The **minimal polynomial** $m(\lambda)$ is the monic polynomial of lowest degree such that $m(A) = O$.
    **Properties**:
    1.  $m(\lambda)$ divides the characteristic polynomial $p(\lambda)$.
    2.  $A$ is diagonalizable $\iff$ $m(\lambda)$ has no repeated roots.
    3.  The multiplicity of $\lambda_i$ in $m(\lambda)$ is equal to the size of the largest Jordan block corresponding to that eigenvalue.

---

## Exercises

1. **[Jordan Block] Write the square of $J_2(5)$.**
   ??? success "Solution"
       $\begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix}^2 = \begin{pmatrix} 25 & 10 \\ 0 & 25 \end{pmatrix}$.

2. **[Diagonalization] Determine if $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ is diagonalizable.**
   ??? success "Solution"
       No. The eigenvalue 2 has algebraic multiplicity 2, but the dimension of the eigenspace (geometric multiplicity) is 1.

3. **[JCF Search] If a $3 \times 3$ matrix $A$ has all eigenvalues equal to 0 and $\operatorname{rank}(A)=1$, what is its JCF?**
   ??? success "Solution"
       $\operatorname{rank}(A)=1 \implies \dim \ker(A) = 2$, meaning there are 2 Jordan blocks. Since the total dimension is 3, it must be $\operatorname{diag}(J_2(0), J_1(0))$.

4. **[Minimal] Given $J = \operatorname{diag}(J_3(2), J_2(2))$, find its minimal polynomial.**
   ??? success "Solution"
       The size of the largest block for eigenvalue 2 is 3, so $m(\lambda) = (\lambda-2)^3$.

5. **[Possibility] If the characteristic polynomial is $(\lambda-1)^4$ and the minimal polynomial is $(\lambda-1)^2$, find possible JCFs.**
   ??? success "Solution"
       The largest block size is 2. Possible forms are $\operatorname{diag}(J_2(1), J_2(1))$ or $\operatorname{diag}(J_2(1), J_1(1), J_1(1))$.

6. **[Nilpotency] Describe $J_k(0)^n$ for $n \ge k$.**
   ??? success "Solution"
       It results in the zero matrix (Nilpotency).

7. **[Rank Method] How do you determine the number of blocks of size $\ge 2$ for eigenvalue $\lambda$?**
   ??? success "Solution"
       It is given by $\operatorname{rank}(A-\lambda I) - \operatorname{rank}(A-\lambda I)^2$.

8. **[Spaces] What is the difference between an eigenspace and a generalized eigenspace?**
   ??? success "Solution"
       The eigenspace is $\ker(A-\lambda I)$, while the generalized eigenspace is $\ker(A-\lambda I)^k$ (where $k$ is the algebraic multiplicity).

9. **[Uniqueness] If $A$ and $B$ have the same JCF, are they similar?**
   ??? success "Solution"
       Yes. The JCF is a complete invariant under similarity.

10. **[Numerical] Why is JCF rarely computed directly in numerical software?**
    ??? success "Solution"
        The JCF is extremely sensitive to small perturbations in matrix entries (discontinuous), making it unstable under floating-point arithmetic.

## Chapter Summary

The Jordan Canonical Form is the final verdict on the structure of square matrices:

1.  **Completing the Defective**: It perfectly fills the gap in eigenvectors for non-diagonalizable matrices by introducing the "1" step structure (Jordan block).
2.  **Polynomial Depth**: The correspondence between the minimal polynomial and JCF block sizes reveals the geometric depth of a matrix as a root of a polynomial.
3.  **Structural Uniqueness**: JCF establishes the classification standard for matrix similarity classes, providing the most precise theoretical framework for analyzing matrix functions and power series convergence.
