# Chapter 35: The Hadamard Product

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Positive Definiteness (Ch16) · Inner Product (Ch8)

**Chapter Outline**: Definition of Hadamard (Schur) Product → Properties → Schur Product Theorem (PSD Preservation) → Oppenheim's Inequality → Hadamard Product and Rank → Spectral Norm Bounds → Application in Kernel Methods and Statistics

**Extension**: The Hadamard product is the basis for element-wise operations in neural networks and the construction of positive definite kernels in machine learning.

</div>

The **Hadamard product** (also known as the Schur product) is the element-wise multiplication of two matrices of the same dimension. While standard matrix multiplication represents the composition of linear maps, the Hadamard product arises in situations requiring entry-specific scaling, such as mask operations in image processing or kernel functions in statistics. The most profound result in this area is the **Schur Product Theorem**, which states that the Hadamard product of two positive semi-definite matrices remains positive semi-definite.

---

## 35.1 Definitions and the Schur Product Theorem

!!! definition "Definition 35.1 (Hadamard Product)"
    The Hadamard product of $A = (a_{ij})$ and $B = (b_{ij})$, denoted $A \circ B$, is the matrix:
    $$(A \circ B)_{ij} = a_{ij} b_{ij}$$

!!! theorem "Theorem 35.1 (Schur Product Theorem)"
    If $A \succeq 0$ and $B \succeq 0$ are $n \times n$ matrices, then $A \circ B \succeq 0$.

---

## Exercises

1. **[Fundamentals] Compute $A \circ B$ for $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ and $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $A \circ B = \begin{pmatrix} 1 \cdot 0 & 2 \cdot 1 \\ 3 \cdot 1 & 4 \cdot 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ 3 & 0 \end{pmatrix}$.

2. **[Schur Theorem] Prove the Schur Product Theorem using the property that $A \succeq 0$ can be written as a sum of rank-1 matrices $v_k v_k^*$.**
   ??? success "Solution"
       Let $A = \sum v_k v_k^*$ and $B = \sum w_l w_l^*$. Then $A \circ B = \sum_{k,l} (v_k v_k^*) \circ (w_l w_l^*) = \sum_{k,l} (v_k \circ w_l) (v_k \circ w_l)^*$. Each term in the sum is a rank-1 PSD matrix, so the sum is PSD.

3. **[Oppenheim] State Oppenheim's Inequality for $A, B \succeq 0$.**
   ??? success "Solution"
       $\det(A \circ B) \ge \det A \cdot \prod b_{ii}$. This provides a lower bound on the determinant of the Hadamard product.

4. **[Eigenvalues] How do the eigenvalues of $A \circ B$ relate to those of $A$ and $B$?**
   ??? success "Solution"
       The spectral radius satisfies $\rho(A \circ B) \le \rho(A) \rho(B)$ is false. However, it is true that $\lambda_{\max}(A \circ B) \le \lambda_{\max}(A) \cdot \max b_{ii}$.

5. **[Spectral Norm] Prove that $\|A \circ B\| \le \|A\| \cdot \|B\|$.**
   ??? success "Solution"
       This follows from the fact that $A \circ B$ is a principal submatrix of the Kronecker product $A \otimes B$. Since the spectral norm of $A \otimes B$ is $\|A\| \cdot \|B\|$, and taking submatrices does not increase the spectral norm, the inequality holds.

6. **[Rank] What is the maximum possible rank of $A \circ B$?**
   ??? success "Solution"
       $\operatorname{rank}(A \circ B) \le \operatorname{rank}(A) \cdot \operatorname{rank}(B)$.

7. **[Kernel Methods] Why is the Hadamard product important for Gaussian kernels?**
   ??? success "Solution"
       The Gaussian kernel $K(x, y) = \exp(-\gamma \|x-y\|^2)$ can be viewed as a Hadamard product of simpler kernels. The Schur Product Theorem ensures that the product of valid kernels remains a valid positive definite kernel.

8. **[Diagonal] Express the diagonal of the product $AB$ using the Hadamard product.**
   ??? success "Solution"
       $\operatorname{diag}(AB) = (A \circ B^T) \mathbf{1}$.

9. **[Identity] What is the identity element for the Hadamard product?**
   ??? success "Solution"
       The all-ones matrix $J$, where $J_{ij} = 1$ for all $i, j$.

10. **[Trace] Prove the Schur-Hadamard trace identity $\operatorname{tr}((A \circ B)C) = \operatorname{tr}(A(B \circ C))$ for symmetric $B$.**
    ??? success "Solution"
        $\sum_{i,j} a_{ij} b_{ij} c_{ji} = \sum_{i,j} a_{ij} (b_{ij} c_{ij})$ assuming symmetry of $B$ and $C$. Both sides expand to the same sum of element-wise products.

## Chapter Summary

This chapter explores the element-wise calculus of matrices:

1. **Analytical Power**: Established the Schur Product Theorem as the defining positivity property of the Hadamard product.
2. **Determinantal Bounds**: Utilized Oppenheim's inequality to constrain the volume of Hadamard-transformed operators.
3. **Subspace Dynamics**: Analyzed rank and norm relations, positioning the Hadamard product as a contraction of the Kronecker product.
4. **Statistical Relevance**: Highlighted its role in covariance structure modeling and the construction of positive definite kernels.
