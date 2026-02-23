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

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ## Chapter Summary

The determinant is the "volume tag" of a square matrix:


****: Whether the determinant is zero is the ultimate divide between invertible and singular matrices, and between systems with or without unique solutions.

****: Laplace expansion recursion connects high-order determinants to lower-order submatrices, while the multiplicative property demonstrates the accumulation of volume after the composition of linear operators.

****: The adjoint matrix not only provides an explicit formula for inversion but also deeply reveals the algebraic relationship between a matrix and its internal sub-blocks.
