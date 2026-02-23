# Chapter 11: Singular Value Decomposition (SVD)

<div class="context-flow" markdown>

**Prerequisites**: Matrix Decompositions (Ch10) · Eigenvalues (Ch06) · Orthogonality (Ch07)

**Chapter Outline**: Existence Theorem of SVD → Geometric Interpretation of Singular Values and Vectors → Compact vs. Truncated SVD → Relation to Eigenvalues of $A^T A$ → Best Low-rank Approximation (Eckart-Young Theorem) → Moore-Penrose Pseudoinverse ($A^+$) → Applications: PCA, Denoising, and Image Compression

**Extension**: SVD is arguably the most powerful decomposition tool in linear algebra, applicable to any matrix (square or rectangular, full-rank or deficient); it is the mathematical bedrock of modern Data Science and Machine Learning.

</div>

The Singular Value Decomposition (SVD) is often called the "pinnacle" of linear algebra. While eigenvalue decomposition deconstructs square matrices, SVD provides the ultimate analysis for *any* matrix. It not only reveals the rank structure but also defines how a matrix, as a linear mapping, stretches and rotates space.

---

## 11.1 Definition and Existence

!!! theorem "Theorem 11.1 (Existence of SVD)"
    For any $m \times n$ matrix $A$, there exists a factorization:
    $$A = U \Sigma V^T$$
    where:
    - $U$ is an $m \times m$ orthogonal matrix, whose columns are **left singular vectors**.
    - $V$ is an $n \times n$ orthogonal matrix, whose columns are **right singular vectors**.
    - $\Sigma$ is an $m \times n$ rectangular diagonal matrix with diagonal entries $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$, known as **singular values**.

---

## 11.2 Geometric Interpretation

!!! technique "Geometry: Transformation of the Hyper-ellipsoid"
    SVD shows that any linear mapping can be decomposed into three steps:
    1.  **Rotation**: $V^T$ rotates the input space to align with the principal axes.
    2.  **Scaling**: $\Sigma$ stretches or compresses space along these axes.
    3.  **Final Rotation**: $U$ rotates the result to its final position in the output space.
    Thus, $A$ maps a unit sphere to a hyper-ellipsoid, where the singular values are the lengths of the semi-axes.

---

## 11.3 Best Low-rank Approximation

!!! theorem "Theorem 11.2 (Eckart-Young Theorem)"
    Among all matrices of rank $k$ ($k < r$), the matrix $A_k$ that minimizes $\|A - A_k\|_F$ is the truncated SVD:
    $$A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$$
    **Application**: This is the theoretical core of data compression, noise reduction, and Principal Component Analysis (PCA).

---

## 11.4 Moore-Penrose Pseudoinverse

!!! definition "Definition 11.1 (Pseudoinverse $A^+$)"
    SVD allows for the definition of a **generalized inverse** for non-square or singular matrices:
    $$A^+ = V \Sigma^+ U^T$$
    where $\Sigma^+$ is obtained by taking the reciprocal of non-zero diagonal entries of $\Sigma$ and transposing the result.
    **Property**: $A^+ \mathbf{b}$ is the minimum-norm solution to the linear least squares problem.

---

## Exercises

1. **[Calculation] Find the singular values of $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$.**
   ??? success "Solution"
       Singular values are magnitudes (positive), so $\sigma_1 = 3, \sigma_2 = 2$.

2. **[Relation] Prove that the squares of the singular values of $A$ are the eigenvalues of $A^T A$.**
   ??? success "Solution"
       $A^T A = (V \Sigma^T U^T)(U \Sigma V^T) = V (\Sigma^T \Sigma) V^T$. This is an eigenvalue decomposition with eigenvalues $\sigma_i^2$.

3. **[Rank] If the singular values of $A$ are $5, 3, 0, 0$, what is $\operatorname{rank}(A)$?**
   ??? success "Solution"
       Rank equals the number of non-zero singular values, which is 2.

4. **[Norm] Prove $\|A\|_2 = \sigma_1$.**
   ??? success "Solution"
       The spectral norm is defined as the maximum stretch factor. From SVD, the maximum stretch is clearly the largest singular value $\sigma_1$.

5. **[Low-rank] For a matrix with singular values $10, 8, 1$, what is the Frobenius norm error of the best rank-2 approximation?**
   ??? success "Solution"
       The error is the square root of the sum of squares of discarded singular values: $\sqrt{1^2} = 1$.

6. **[Pseudoinverse] If $A = \operatorname{diag}(2, 0)$, find $A^+$.**
   ??? success "Solution"
       $A^+ = \operatorname{diag}(1/2, 0)$.

7. **[Uniqueness] Are $U$ and $V$ in the SVD unique?**
   ??? success "Solution"
       No (e.g., signs of columns can be flipped), but the singular values in $\Sigma$ are unique.

8. **[Compact] What is the Compact SVD?**
   ??? success "Solution"
       It retains only the columns/rows corresponding to non-zero singular values: $A = U_r \Sigma_r V_r^T$, where $\Sigma_r$ is an $r \times r$ positive definite diagonal matrix.

9. **[Compression] Why is SVD used for image compression?**
   ??? success "Solution"
       Image matrices are often low-rank or approximately so. By keeping only the top $k$ singular values, one can reconstruct the image's main features using much less data.

10. **[Condition] How is the condition number $\kappa(A)$ expressed in terms of singular values?**
    ??? success "Solution"
        $\kappa(A) = \sigma_{\max} / \sigma_{\min}$ (for square matrices).

## Chapter Summary

SVD is the ultimate weapon for solving real-world linear algebra problems:

1.  **Universality**: It bridges the gap between square and rectangular, full-rank and deficient matrices, providing a unified analytic perspective for all linear mappings.
2.  **Energy Concentration**: By ordering singular values, SVD identifies the "most important" components of data, providing mathematical criteria for compression and dimensionality reduction.
3.  **Numerical Robustness**: As the basis for pseudoinverses and least squares solutions, SVD exhibits high stability against ill-conditioned matrices, serving as the final line of defense in engineering computations.
