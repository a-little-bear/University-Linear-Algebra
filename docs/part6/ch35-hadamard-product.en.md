# Chapter 35: Hadamard Product

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Positive Definite Matrices (Ch16) · Matrix Inequalities (Ch18)

**Chapter Outline**: Definition of the Hadamard Product (Element-wise Product) → Basic Algebraic Properties → The Schur Product Theorem (Preservation of Positivity) → Hadamard Product Inequalities (Oppenheim, Hadamard Inequalities) → Spectral Properties and Singular Value Bounds → Applications in Statistics (Correlation Matrices) → Windowing in Signal Processing → Relation to the Kronecker Product

**Extension**: The Hadamard product introduces pointwise scalar multiplication into matrix space; it is key to understanding non-linear combinations in Kernel Methods (Ch29) and sparsification operations in modern compressed sensing algorithms.

</div>

Unlike standard matrix multiplication (which reflects operator composition), the **Hadamard Product** (also known as the Schur product) is performed element-wise. Although it appears simpler algebraically, its preservation of positive semi-definiteness (the Schur Product Theorem) and the rich inequalities derived from it give it immense theoretical value in statistical modeling, image processing, and numerical preconditioning.

---

## 35.1 Definition and Basic Properties

!!! definition "Definition 35.1 (Hadamard Product)"
    Let $A$ and $B$ be matrices of the same dimension. Their **Hadamard Product** $A \circ B$ is a matrix of the same dimension with entries:
    $$(A \circ B)_{ij} = a_{ij} b_{ij}$$

!!! note "Algebraic Properties"
    1.  **Commutativity**: $A \circ B = B \circ A$.
    2.  **Distributivity**: $A \circ (B + C) = A \circ B + A \circ C$.
    3.  **Identity**: The all-ones matrix $J$ is the identity element for the Hadamard product.

---

## 35.2 The Schur Product Theorem

!!! theorem "Theorem 35.1 (Schur Product Theorem)"
    If $A \succeq 0$ and $B \succeq 0$ are positive semi-definite matrices, then their Hadamard product is also positive semi-definite:
    $$A \circ B \succeq 0$$
    **Significance**: This property ensures that in kernel methods, the pointwise product of two valid kernels is itself a valid kernel.

---

## 35.3 Hadamard Product Inequalities

!!! theorem "Theorem 35.2 (Oppenheim's Inequality)"
    For positive definite matrices $A, B \succ 0$:
    $$\det(A \circ B) \ge \left( \prod_{i=1}^n a_{ii} \right) \det(B) \ge \det(A) \det(B)$$

!!! theorem "Theorem 35.3 (Hadamard's Inequality)"
    As a special case of the Schur product, for $A \succ 0$:
    $$\det(A) \le \prod_{i=1}^n a_{ii}$$
    This can be viewed as a generalized Hadamard interaction between $A$ and the identity matrix $I$.

---

## Exercises

1.  **[Basic] Calculate $\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix} \circ \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
    ??? success "Solution"
        $\begin{pmatrix} 0 & 2 \\ 3 & 0 \end{pmatrix}$.

2.  **[Schur] If $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$, find $A \circ A$.**
    ??? success "Solution"
        It is still $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.

3.  **[Property] Prove $\operatorname{tr}(A \circ B) = \operatorname{tr}(A^T B)$.**
    ??? success "Solution"
        $\sum a_{ii} b_{ii} = \sum a_{ij} b_{ij}$ (only for diagonal summation). More generally, $\sum a_{ij} b_{ij} = \operatorname{tr}(A B^T)$. For symmetric matrices, this reduces to the stated trace property.

4.  **[Positivity] Determine if $\begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix} \circ \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$ is positive definite.**
    ??? success "Solution"
        Yes. By the Schur Product Theorem, the Hadamard product of two PD matrices is PD.

5.  **[Bound] Does the inequality $\|A \circ B\|_2 \le \|A\|_2 \|B\|_2$ always hold?**
    ??? success "Solution"
        Not necessarily for the spectral norm. It holds for the Frobenius norm ($\|A \circ B\|_F \le \|A\|_F \|B\|_F$), but spectral norm bounds usually involve diagonal entries.

6.  **[Application] Why is windowing in signal processing a Hadamard product?**
    ??? success "Solution"
        Windowing involves scaling each sample of a signal by a corresponding weight from a window function, which is exactly element-wise multiplication.

7.  **[Kronecker] Is the Hadamard product $A \circ B$ a submatrix of $A \otimes B$?**
    ??? success "Solution"
        Yes. It is precisely the principal submatrix of the Kronecker product corresponding to specific row and column indices.

8.  **[Determinant] Verify Oppenheim's inequality for diagonal matrices.**
    ??? success "Solution"
        For diagonal matrices, $A \circ B = AB$, so $\det(AB) = \det(A)\det(B)$. The inequality becomes an equality.

9.  **[Rank] Prove $\operatorname{rank}(A \circ B) \le \operatorname{rank}(A) \operatorname{rank}(B)$.**
    ??? success "Solution"
        Since $A \circ B$ is a submatrix of $A \otimes B$, its rank cannot exceed $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \operatorname{rank}(B)$.

10. **[Correlation] Is the Hadamard product of two correlation matrices still a correlation matrix?**

   ??? success "Solution"
        Yes. It preserves symmetry, the property of having ones on the diagonal, and positive semi-definiteness.

## Chapter Summary

The Hadamard product is a refined operator in matrix analysis:

1.  **Pointwise Harmony**: It preserves the simplicity of scalar multiplication in matrix space while revealing profound non-trivial preservation laws at the level of positive semi-definiteness (Schur Product Theorem).
2.  **Information Compression**: Inequalities like Oppenheim’s demonstrate how element-wise coupling affects the global integrity of a matrix (e.g., determinant), providing algebraic descriptions for feature interactions in information theory.
3.  **Computational Bridge**: As a projected microcosm of the Kronecker product, the Hadamard product plays an irreplaceable role in dimension reduction when handling high-dimensional sparse data and structured correlation models.
