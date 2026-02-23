# Chapter 10: Matrix Factorizations

<div class="context-flow" markdown>

**Prerequisites**: Matrix Multiplication (Ch2) · Gaussian Elimination (Ch1) · Positive Definiteness (Ch16) · QR Decomposition (Ch7)

**Chapter Outline**: LU Factorization ($A = LU$) → Partial Pivoting ($PA = LU$) → Cholesky Factorization ($A = LL^T$) → LDL Factorization ($A = LDL^T$) → Spectral Decomposition ($A = Q \Lambda Q^T$) → Singular Value Decomposition (SVD) Overview → Computational Efficiency → Uniqueness of Factorizations

**Extension**: Matrix factorizations are the "atomic decompositions" of operators; they break down a complex transformation into a sequence of simpler, more manageable steps.

</div>

Matrix factorizations (or decompositions) express a matrix as a product of simpler matrices with specific structures (triangular, orthogonal, diagonal). These representations are the foundation of numerical linear algebra. The **LU factorization** automates Gaussian elimination for repeated solving, while the **Cholesky factorization** provides a highly efficient method for symmetric positive definite matrices. This chapter details the derivation and properties of these core factorizations and discusses their role in optimizing computational resources.

---

## 10.1 Triangular Factorizations

!!! definition "Definition 10.1 (LU Factorization)"
    The LU factorization of $A$ is $A = LU$, where $L$ is a lower triangular matrix with 1s on the diagonal and $U$ is an upper triangular matrix. It corresponds to the row operations of Gaussian elimination.

!!! theorem "Theorem 10.1 (Cholesky Existence)"
    If $A$ is a symmetric positive definite matrix, there exists a unique lower triangular matrix $L$ with positive diagonal entries such that $A = LL^T$.

---

## Exercises

1. **[Fundamentals] Find the LU factorization of $A = \begin{pmatrix} 2 & 1 \\ 4 & 4 \end{pmatrix}$.**
   ??? success "Solution"
       Step 1: $R_2 \to R_2 - 2R_1$. The upper triangular result is $U = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$. The multiplier was 2, so $L = \begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$.

2. **[Cholesky] Compute the Cholesky factorization of $A = \begin{pmatrix} 4 & 2 \\ 2 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       Let $L = \begin{pmatrix} l_{11} & 0 \\ l_{21} & l_{22} \end{pmatrix}$. $l_{11} = \sqrt{4} = 2$. $l_{21} = 2/2 = 1$. $l_{22} = \sqrt{2 - 1^2} = 1$. Thus $L = \begin{pmatrix} 2 & 0 \\ 1 & 1 \end{pmatrix}$.

3. **[Pivoting] Why is $PA = LU$ used instead of $A = LU$?**
   ??? success "Solution"
       For numerical stability. Partial pivoting ($P$) ensures that the largest available entry is used as a pivot, avoiding division by zero or very small numbers that amplify rounding errors.

4. **[Complexity] Compare the FLOP count of Cholesky versus LU.**
   ??? success "Solution"
       Cholesky requires approximately $n^3/3$ operations, whereas LU requires $2n^3/3$. By exploiting symmetry, Cholesky is twice as fast.

5. **[Solving] Solve $Ax = b$ using $A = LU$.**
   ??? success "Solution"
       1. Solve $Ly = b$ via forward substitution. 2. Solve $Ux = y$ via back substitution. This approach is efficient when solving for multiple $b$ vectors.

6. **[LDL] What is the LDL factorization?**
   ??? success "Solution"
       $A = LDL^T$ where $L$ is unit lower triangular and $D$ is diagonal. It generalizes Cholesky to symmetric matrices that are not necessarily positive definite (it avoids square roots).

7. **[Uniqueness] Under what condition is the LU factorization unique?**
   ??? success "Solution"
       When all leading principal minors of $A$ are non-zero. This ensures that Gaussian elimination can proceed without row swaps.

8. **[Determinant] Use LU to compute $\det A$.**
   ??? success "Solution"
       $\det A = \det L \cdot \det U = (1) \cdot \prod u_{ii}$. The determinant is simply the product of the pivots.

9. **[Symmetry] Show that if $A = LU$ and $A$ is symmetric, then $U = DL^T$.**
   ??? success "Solution"
       $A = A^T \implies LU = U^T L^T$. Since $L$ is unit lower triangular and $U^T$ is unit lower triangular (if scaled), the unique $LDU$ factorization implies the relation.

10. **[Memory] How can $L$ and $U$ be stored in the same memory space as $A$?**
    ??? success "Solution"
        Since $L$ has 1s on the diagonal, only the off-diagonal entries of $L$ need to be stored. These can occupy the lower triangular part of the original matrix $A$, while $U$ occupies the upper triangular part (including the diagonal).

## Chapter Summary

This chapter establishes the core structural decompositions of matrices:

1. **Procedural Encoding**: Formulated LU factorization as the persistent record of Gaussian elimination steps.
2. **Symmetry Optimization**: Developed Cholesky and LDL factorizations to halve computational costs for symmetric systems.
3. **Numerical Robustness**: Introduced partial pivoting as the essential mechanism for stable matrix computation.
4. **Hierarchical Solving**: Linked factorizations to forward and back substitution, providing the standard template for large-scale linear solvers.
