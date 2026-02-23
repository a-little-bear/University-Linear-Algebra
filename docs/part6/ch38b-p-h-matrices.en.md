# Chapter 38B: P-Matrices, H-Matrices and Related Classes

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · M-Matrices (Ch38A) · Eigenvalues (Ch6) · LCP Basics

**Chapter Outline**: P-Matrices & Principal Minors → LCP Unique Solvability → P0-Matrices → Comparison Matrices & H-Matrices → Diagonal Dominance Hierarchy → N-Matrices → Semipositive Matrices → Ostrowski-Reich Theorem → Iterative Stability

**Extension**: P-matrices are central to mathematical programming (LCP theory); H-matrices generalize M-matrix theory to complex matrices and are vital for the convergence analysis of iterative methods.

</div>

This chapter explores broader matrix classes than Z-matrices. **P-matrices** require all principal minors to be positive, with no restriction on sign patterns. **H-matrices** utilize the comparison matrix to extend M-matrix properties to general complex matrices. These classes provide the analytical foundation for Linear Complementarity Problems (LCP) and numerical stability.

---

## 38B.1 P-Matrices

!!! definition "Definition 38B.1 (P-Matrix)"
    A matrix $A \in M_n(\mathbb{R})$ is a **P-matrix** if all its principal minors are positive:
    $$\det(A[\alpha, \alpha]) > 0, \quad \forall\, \alpha \subseteq \{1, \dots, n\}, \alpha 
e \emptyset.$$

!!! theorem "Theorem 38B.3 (LCP Solvability)"
    $A$ is a P-matrix if and only if the Linear Complementarity Problem $\operatorname{LCP}(q, A)$ has a unique solution for every $q \in \mathbb{R}^n$.

---

## Exercises

1. **[Fundamentals] Determine if $A = \begin{pmatrix} 1 & -2 \ 1 & 1 \end{pmatrix}$ is a P-matrix.**

   ??? success "Solution"
       Principal minors: $a_{11}=1 > 0$, $a_{22}=1 > 0$, and $\det A = 1 - (-2) = 3 > 0$. Since all are positive, $A$ is a P-matrix. Note that it is not a Z-matrix.

2. **[Transpose] Prove: If $A$ is a P-matrix, then $A^T$ is also a P-matrix.**

   ??? success "Solution"
       The principal minors of $A^T$ are the determinants of $(A[\alpha, \alpha])^T$, which equal the determinants of $A[\alpha, \alpha]$. Since the latter are positive, so are the former.

3. **[Real Eigenvalues] Show that all real eigenvalues of a P-matrix must be positive.**

   ??? success "Solution"
       If $\lambda$ is a real eigenvalue with eigenvector $v$, then $Av = \lambda v$. By property P2, there exists $i$ such that $v_i(Av)_i = \lambda v_i^2 > 0$. This forces $\lambda > 0$.

4. **[M-Matrices] What is the intersection of the class of P-matrices and Z-matrices?**

   ??? success "Solution"
       The intersection is exactly the class of non-singular M-matrices.

5. **[Comparison Matrix] Compute the comparison matrix $\mathcal{M}(A)$ for $A = \begin{pmatrix} 4 & 1-i \ -1+2i & 5 \end{pmatrix}$.**

   ??? success "Solution"
       $\mathcal{M}(A)_{ii} = |a_{ii}|$, $\mathcal{M}(A)_{ij} = -|a_{ij}|$.
       $\mathcal{M}(A) = \begin{pmatrix} 4 & -\sqrt{2} \ -\sqrt{5} & 5 \end{pmatrix}$.

6. **[H-Matrix Criteria] Determine if $A$ from the previous exercise is an H-matrix.**

   ??? success "Solution"
       $A$ is an H-matrix if $\mathcal{M}(A)$ is an M-matrix. $\det \mathcal{M}(A) = 20 - \sqrt{10} > 0$. Since the diagonal is positive and the determinant is positive, it is a non-singular M-matrix. Thus, $A$ is an H-matrix.

7. **[Stability] State the Ostrowski-Reich Theorem regarding SOR iteration.**

   ??? success "Solution"
       For a symmetric matrix, the SOR iteration converges for $\omega \in (0, 2)$ if and only if $A$ is positive definite. For H-matrices, convergence is guaranteed for $\omega \in (0, 1]$.

8. **[N-Matrices] Define an N-matrix.**

   ??? success "Solution"
       A matrix whose all principal minors are negative. An example is $A = \begin{pmatrix} -1 & 2 \ 3 & -2 \end{pmatrix}$.

9. **[Semipositive] What does it mean for a matrix to be semipositive?**

   ??? success "Solution"
       There exists a vector $x \ge 0$ ($x 
eq 0$) such that $Ax > 0$. Every non-singular M-matrix is semipositive.

10. **[LCP Example] Solve $\operatorname{LCP}(q, A)$ for $A = \begin{pmatrix} 2 & -1 \ -1 & 2 \end{pmatrix}$ and $q = \begin{pmatrix} -1 \ -1 \end{pmatrix}$.**

   ??? success "Solution"
        Trying $z > 0$, we solve $Az = -q \implies z = A^{-1} \begin{pmatrix} 1 \ 1 \end{pmatrix} = \begin{pmatrix} 1 & 1/2 \ 1/2 & 1 \end{pmatrix} \begin{pmatrix} 1 \ 1 \end{pmatrix} \frac{2}{3} = \begin{pmatrix} 1 \ 1 \end{pmatrix}$. The solution is $z = (1, 1)^T, w = (0, 0)^T$.

## Chapter Summary

This chapter classifies matrices based on the positivity of their sub-structures:

1. **Sign Agnostic Positivity**: Defined P-matrices through principal minors, establishing their role in LCP solvability.
2. **Comparison Theory**: Developed H-matrices as the complex generalization of M-matrices via magnitude-based dominance.
3. **Hierarchical Convergence**: Linked these classes to the numerical stability of iterative solvers and error bounds in interval arithmetic.
