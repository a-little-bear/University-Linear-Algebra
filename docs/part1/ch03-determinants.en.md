# Chapter 03: Determinants

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02)

**Chapter Outline**: Definition of Determinant (Permutations & Inversions) → 2x2 and 3x3 Intuition → Key Properties (Effect of Row/Column Operations) → Minors and Cofactors → Laplace Expansion Theorem → Determinants of Products and Inverses → Adjoint Matrix $A^*$ & Inverse Formula → Cramer's Rule → Geometric Interpretation (Volume Scaling Factor)

**Extension**: The determinant is the unique scalar indicator of matrix invertibility and serves as the polynomial engine for studying Eigenvalues (Ch06).

</div>

The determinant is a function that maps a square matrix to a scalar value. Although solving systems via determinants is computationally inefficient, they hold an irreplaceable position in theoretical analysis, area/volume calculation, and eigenvalue determination. This chapter begins with the algebraic definition using permutations and ultimately reveals its geometric essence.

---

## 03.1 Definition of the Determinant

!!! definition "Definition 03.1 (Permutations and Inversions)"
    The number of **inversions** $\tau(\sigma)$ in a permutation of $n$ numbers is the count of pairs where a larger number precedes a smaller one. The sign of a permutation is defined as $\operatorname{sgn}(\sigma) = (-1)^{\tau(\sigma)}$.

!!! definition "Definition 03.2 (Leibniz Formula)"
    The **determinant** of an $n \times n$ matrix $A$, denoted $\det(A)$ or $|A|$, is defined as:
    $$\det(A) = \sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) a_{1,\sigma(1)} a_{2,\sigma(2)} \cdots a_{n,\sigma(n)}$$

---

## 03.2 Properties of Determinants

!!! theorem "Theorem 03.1 (Core Properties)"
    1.  **Transpose Invariance**: $\det(A^T) = \det(A)$.
    2.  **Multiplicative Property**: $\det(AB) = \det(A) \det(B)$.
    3.  **Row Operation Effects**:
        - Swapping two rows flips the sign.
        - Scaling a row by $k$ scales the determinant by $k$.
        - Adding a multiple of one row to another **does not change** the determinant.
    4.  **Invertibility**: $A$ is invertible $\iff \det(A) \neq 0$.

---

## 03.3 Minors and Laplace Expansion

!!! definition "Definition 03.3 (Minors and Cofactors)"
    The determinant of the submatrix formed by deleting the $i$-th row and $j$-th column of $A$ is the **minor** $M_{ij}$.
    The **cofactor** is defined as $C_{ij} = (-1)^{i+j} M_{ij}$.

!!! theorem "Theorem 03.2 (Laplace Expansion)"
    The determinant equals the sum of the products of elements of any row (or column) and their corresponding cofactors:
    $$\det(A) = \sum_{j=1}^n a_{ij} C_{ij} \quad (\text{for a fixed } i)$$

---

## 03.4 The Adjoint Matrix and Cramer's Rule

!!! definition "Definition 03.4 (Adjoint Matrix $A^*$)"
    The transpose of the matrix of cofactors is called the **adjoint** (or adjugate) matrix of $A$:
    $$(A^*)_{ij} = C_{ji}$$
    It satisfies the identity: $AA^* = A^*A = \det(A)I$. If $\det(A) \neq 0$, then $A^{-1} = \frac{1}{\det(A)} A^*$.

!!! technique "Application: Cramer's Rule"
    For a system $Ax = b$, if $\det(A) \neq 0$, then $x_j = \frac{\det(A_j)}{\det(A)}$, where $A_j$ is obtained by replacing the $j$-th column of $A$ with $b$.

---

## Exercises

1. **[Calculation] Compute the $2 \times 2$ determinant $\begin{vmatrix} a & b \\ c & d \end{vmatrix}$.**

   ??? success "Solution"
       $\det = ad - bc$.

2. **[Scaling] If $\det(A) = 5$ and $A$ is $3 \times 3$, calculate $\det(2A)$.**

   ??? success "Solution"
       For an $n \times n$ matrix, $\det(kA) = k^n \det(A)$. Thus $\det(2A) = 2^3 \cdot 5 = 40$.

3. **[Special] Compute the determinant of $\begin{pmatrix} 1 & 2 & 3 \\ 0 & 4 & 5 \\ 0 & 0 & 6 \end{pmatrix}$.**

   ??? success "Solution"
       The determinant of a triangular matrix is the product of its diagonal entries: $1 \cdot 4 \cdot 6 = 24$.

4. **[Laplace] Use expansion along the first column to compute $\begin{vmatrix} 1 & 0 & 2 \\ 0 & 3 & 0 \\ 4 & 0 & 5 \end{vmatrix}$.**

   ??? success "Solution"
       $\det = 1 \cdot \begin{vmatrix} 3 & 0 \\ 0 & 5 \end{vmatrix} + 4 \cdot \begin{vmatrix} 0 & 2 \\ 3 & 0 \end{vmatrix} = 15 + 4(-6) = -9$.

5. **[Powers] If $\det(A^2) = 9$, find possible values for $\det(A)$.**

   ??? success "Solution"
       $\det(A)^2 = 9 \implies \det(A) = \pm 3$.

6. **[Adjoint] Prove $\det(A^*) = (\det A)^{n-1}$.**

   ??? success "Solution"
       From $AA^* = (\det A)I$, taking the determinant gives $(\det A)\det(A^*) = (\det A)^n$. If $\det A \neq 0$, the result follows. By continuity, it holds for singular matrices too.

7. **[Inverse] Prove $\det(A^{-1}) = 1/\det(A)$.**

   ??? success "Solution"
       $\det(A A^{-1}) = \det(I) = 1$. Using the multiplicative property, $\det(A)\det(A^{-1}) = 1$.

8. **[Inversion] What is the number of inversions in the permutation $(3, 1, 2)$?**

   ??? success "Solution"
       The pairs are (3,1) and (3,2). The count is 2.

9. **[Cramer] Use Cramer's rule to evaluate $x+y=1, x+y=2$.**

   ??? success "Solution"
       The coefficient determinant $D = \begin{vmatrix} 1 & 1 \\ 1 & 1 \end{vmatrix} = 0$. Cramer's rule is not applicable (no unique solution).

10. **[Geometry] What does $\det \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}$ represent geometrically?**

   ??? success "Solution"
        It represents the area of a rectangle with sides $a$ and $b$.

## Chapter Summary

The determinant is the "volume tag" of a square matrix:

1.  **Classification**: Whether the determinant is zero is the ultimate divide between invertible and singular matrices, and between systems with or without unique solutions.
2.  **Structural Properties**: Laplace expansion recursion connects high-order determinants to lower-order submatrices, while the multiplicative property demonstrates the accumulation of volume after the composition of linear operators.
3.  **Adjoint Tools**: The adjoint matrix not only provides an explicit formula for inversion but also deeply reveals the algebraic relationship between a matrix and its internal sub-blocks.
