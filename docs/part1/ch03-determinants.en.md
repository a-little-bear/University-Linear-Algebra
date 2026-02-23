# Chapter 03: Determinants

<div class="context-flow" markdown>

**Prerequisites**: Matrix Basics (Ch2)

**Chapter Outline**: Definition of Determinant (Permutations and Inversions) → 2nd and 3rd Order Expansion → Properties of Determinants (Effect of Row/Column Ops) → Cofactors and Laplace Expansion → Cramer's Rule → Adjoint Matrix $A^*$ → Determinant Formula for Inverses

**Extension**: Geometrically, the determinant represents the "volume scaling factor" caused by a linear mapping.

</div>

The determinant is a function that maps a square matrix to a scalar. Although complex to compute, it plays an irreplaceable role in determining matrix invertibility, solving linear systems, and studying eigenvalues.

---

## 03.1 Properties of Determinants

!!! theorem "Theorem 03.1 (Key Properties)"
    1. $\det(AB) = \det(A) \det(B)$
    2. $\det(A^T) = \det(A)$
    3. If a row is multiplied by $k$, the determinant is multiplied by $k$.
    4. Swapping two rows flips the sign.
    5. Adding a multiple of one row to another **leaves the determinant unchanged**.

---

## Exercises

1. **[Calculation] Compute the $2 \times 2$ determinant $\begin{vmatrix} a & b \\ c & d \end{vmatrix}$.**
   ??? success "Solution"
       $\det = ad - bc$$.

2. **[Property App] If $\det(A) = 5$ and $A$ is a $3 \times 3$ matrix, calculate $\det(2A)$.**
   ??? success "Solution"
       According to the properties, each row pulls out a factor of 2. For an $n \times n$ matrix, $\det(kA) = k^n \det(A)$.
       Thus $\det(2A) = 2^3 \cdot 5 = 8 \cdot 5 = 40$.

3. **[Triangular] Compute $\begin{vmatrix} 1 & 2 & 3 \\ 0 & 4 & 5 \\ 0 & 0 & 6 \end{vmatrix}$.**
   ??? success "Solution"
       For upper/lower triangular matrices, the determinant is the product of the diagonal entries.
       $\det = 1 \cdot 4 \cdot 6 = 24$.

4. **[Laplace Expansion] Expand along the first column to calculate $\begin{vmatrix} 1 & 0 & 2 \\ 0 & 3 & 0 \\ 4 & 0 & 5 \end{vmatrix}$.**
   ??? success "Solution"
       $\det = 1 \cdot \begin{vmatrix} 3 & 0 \\ 0 & 5 \end{vmatrix} - 0 + 4 \cdot \begin{vmatrix} 0 & 2 \\ 3 & 0 \end{vmatrix} = 1(15) + 4(-6) = 15 - 24 = -9$.

5. **[Invertibility] If $\det(A) = 0$, is $A$ invertible?**
   ??? success "Solution"
       No (singular matrix). A matrix is invertible if and only if $\det(A) \neq 0$.

6. **[Cramer's Rule] Use Cramer's Rule to determine if $x+y=1, x+y=2$ has a unique solution.**
   ??? success "Solution"
       The coefficient determinant $D = \begin{vmatrix} 1 & 1 \\ 1 & 1 \end{vmatrix} = 0$. Since $D=0$, the system has no unique solution.

7. **[Adjoint Matrix] Known that $AA^* = \det(A)I$. If $\det(A) = 2$, find $AA^*$.**
   ??? success "Solution"
       $AA^* = 2I = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$ (assuming $2 \times 2$).

8. **[Transpose] Prove: If $A$ is skew-symmetric ($A^T = -A$) and has odd order, then $\det(A) = 0$.**
   ??? success "Solution"
       $\det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)$.
       When $n$ is odd, $\det(A) = -\det(A) \implies 2\det(A) = 0 \implies \det(A) = 0$.

9. **[Product] Prove $\det(A^{-1}) = 1/\det(A)$.**
   ??? success "Solution"
       Since $AA^{-1} = I$, taking the determinant gives $\det(A)\det(A^{-1}) = \det(I) = 1$.
       Thus $\det(A^{-1}) = 1/\det(A)$.

10. **[Permutations] What is the number of inversions in $(3, 1, 2)$? What sign does its term carry in the expansion?**
    ??? success "Solution"
        Inversion pairs: (3,1), (3,2). Number of inversions is 2.
        Since 2 is even, the term carries a positive $(+)$ sign.

## Chapter Summary

The determinant is the numerical essence of a square matrix:

1. **Volumetric Meaning**: It captures how much a linear transformation changes the size of a space.
2. **Existence Indicator**: A non-zero determinant is the mark of a "well-behaved" system (invertible, unique solution).
3. **Computational Core**: Laplace expansion and elementary transformation properties form the two wings of determinant calculation.
