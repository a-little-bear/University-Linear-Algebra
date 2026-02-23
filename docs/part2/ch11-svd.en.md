# Chapter 11: Singular Value Decomposition (SVD)

<div class="context-flow" markdown>

**Prerequisites**: Matrix Decompositions (Ch10) · Eigenvalues (Ch06) · Orthogonality (Ch07) · Matrix Norms (Ch15)

**Chapter Outline**: Existence and Uniqueness Theorems of SVD → Geometric Interpretation of Singular Values and Vectors (The Hyper-ellipsoid Transform) → Compact vs. Truncated SVD → Deep Relationship with Eigen-decomposition ($A^* A$ and $AA^*$) → Best Low-rank Approximation (Eckart-Young Theorem) → Moore-Penrose Pseudoinverse ($A^+$) → Statistical Applications: The Algebraic Core of PCA → Signal Processing: Image Compression and Denoising

**Extension**: The SVD is widely considered the "pinnacle" of linear algebra; it removes the divide between square/rectangular and full-rank/deficient matrices, providing the ultimate analytic tool for describing any linear operator. It is the bedrock of modern Data Science and AI.

</div>

Singular Value Decomposition (SVD) is the ultimate deconstruction of any matrix. While eigenvalue decomposition is an "internal deconstruction" for square matrices, SVD provides a "global assessment" for any linear mapping. It not only reveals the rank structure of a matrix but also specifies how the matrix, as an operator, stretches space anisotropically. This chapter demonstrates how SVD reduces a matrix to a concise sum of its energy components.

---

## 11.1 Definition and Geometric Essence

!!! theorem "Theorem 11.1 (Existence of SVD)"
    For any $m \times n$ complex matrix $A$, there exists a factorization:
    $$A = U \Sigma V^*$$
    where:
    - $U \in M_m$ is a unitary matrix, whose columns are **left singular vectors**.
    - $V \in M_n$ is a unitary matrix, whose columns are **right singular vectors**.
    - $\Sigma \in M_{m \times n}$ is a rectangular diagonal matrix with entries $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$, known as **singular values**.

!!! technique "Geometry: Hyper-ellipsoid Transformation"
    SVD shows that any linear mapping can be broken into three steps:
    1.  **Rotation** ($V^*$): Aligning input space axes with principal component directions.
    2.  **Scaling** ($\Sigma$): Stretching or compressing along these axes (singular values are the scaling factors).
    3.  **Final Rotation** ($U$): Rotating the result to its final pose in the output space.
    Thus, $A$ maps a unit sphere to a **hyper-ellipsoid**, where singular values correspond to the semi-axis lengths.

---

## 11.2 Best Low-Rank Approximation

!!! theorem "Theorem 11.2 (Eckart-Young Theorem)"
    Among all matrices of rank $k$ ($k < r$), the matrix $A_k$ that minimizes the Frobenius norm distance $\|A - A_k\|_F$ is the truncated SVD:
    $$A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^*$$
    **Application**: This is the mathematical criterion for data compression (retaining only large singular values) and Principal Component Analysis (PCA).

---

## 11.3 Moore-Penrose Pseudoinverse

!!! definition "Definition 11.1 (Pseudoinverse $A^+$)"
    SVD allows for the definition of a **generalized inverse** for any matrix:
    $$A^+ = V \Sigma^+ U^*$$
    where $\Sigma^+$ is obtained by taking the reciprocals of the non-zero entries of $\Sigma$ and transposing the result.
    **Property**: For any $b$, $\hat{x} = A^+ b$ is the minimum-norm solution to the linear least squares problem.

---

## Exercises

**1. [Calculation] Find the singular values and SVD of $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. Singular values must be non-negative.
    2. Diagonals are 3 and -2.
    3. The singular values are $\sigma_1 = 3, \sigma_2 = 2$.
    4. To keep $\Sigma = \operatorname{diag}(3, 2)$ positive, $U$ must flip the second axis.
    **SVD**:
    $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix} \begin{pmatrix} 3 & 0 \\ 0 & 2 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.
    Here $U = \operatorname{diag}(1, -1)$ and $V = I$.

**2. [Relation] Prove that the squares of the singular values of $A$ are the eigenvalues of $A^* A$.**

??? success "Solution"
    **Proof:**
    1. Substitute SVD: $A^* A = (V \Sigma^T U^*)(U \Sigma V^*) = V (\Sigma^T \Sigma) V^*$.
    2. Since $V$ is unitary and $V^*$ its inverse, this is exactly the **Eigen-decomposition** form of $A^* A$.
    3. The eigenvalue matrix is $\Sigma^T \Sigma = \operatorname{diag}(\sigma_1^2, \sigma_2^2, \ldots)$.
    **Conclusion**: Singular values are the arithmetic square roots of the eigenvalues of $A^* A$.

**3. [Rank] If the singular values of $A$ are $5, 3, 0, 0$, what is $\operatorname{rank}(A)$ and its nullity?**

??? success "Solution"
    **Determination:**
    1. Rank equals the number of non-zero singular values: $\operatorname{rank}(A) = 2$.
    2. Nullity equals the number of columns minus the rank. If $A$ is $4 \times 4$, then $4 - 2 = 2$.
    **Conclusion**: Non-zero singular values reveal the "effective dimensions" through which information propagates.

**4. [Norm] Prove $\|A\|_2 = \sigma_1$.**

??? success "Solution"
    **Proof:**
    1. Operator norm $\|A\|_2 = \max \frac{\|Ax\|}{\|x\|}$.
    2. Substitute SVD: $\|U \Sigma V^* x\| = \|\Sigma (V^* x)\|$ (since $U$ preserves norms).
    3. Let $y = V^* x$. Since $V^*$ is unitary, $\|y\| = \|x\|$.
    4. Problem becomes maximizing $\|\Sigma y\|$ such that $\|y\|=1$.
    5. Since $\Sigma$ is diagonal and $\sigma_1$ is largest, the maximum occurs at $y = e_1$, yielding $\sigma_1$.

**5. [Low-rank] Given singular values $10, 8, 1$. What is the Frobenius norm error of the best rank-2 approximation?**

??? success "Solution"
    **Analysis:**
    1. The best rank-$k$ approximation is obtained by discarding smaller singular values.
    2. The error $\|A - A_k\|_F$ is the square root of the sum of squares of the discarded values.
    3. The discarded value is $\sigma_3 = 1$.
    **Conclusion**: Error is $\sqrt{1^2} = 1$. Most of the "energy" is captured by the rank-2 matrix.

**6. [Pseudoinverse] If $A = \operatorname{diag}(2, 0)$, find $A^+$ and verify $AA^+ A = A$.**

??? success "Solution"
    **Calculation:**
    1. Take the reciprocal of non-zeros, keep zeros at 0.
    2. $A^+ = \operatorname{diag}(0.5, 0)$.
    **Verification:**
    $AA^+ = \operatorname{diag}(1, 0)$.
    $A A^+ A = \operatorname{diag}(1, 0) \operatorname{diag}(2, 0) = \operatorname{diag}(2, 0) = A$.
    Satisfies Penrose Condition 1.

**7. [Uniqueness] Are $U$ and $V$ in the SVD unique?**

??? success "Solution"
    **Conclusion:**
    **No**.
    **Reasoning**:
    1. Signs of columns (or phases in $\mathbb{C}$) can be flipped without changing the product.
    2. If there are repeated singular values, the corresponding vectors can be rotated within their eigenspan.
    However, the sequence of singular values in $\Sigma$ is **unique**.

**8. [Compact SVD] What is the Compact SVD and how does it differ from the Full SVD?**

??? success "Solution"
    **Difference:**
    - **Full SVD**: $U$ is $m \times m$, $V$ is $n \times n$.
    - **Compact SVD**: Retains only columns corresponding to non-zero singular values. If rank is $r$, $U_r$ is $m \times r$ and $V_r$ is $n \times r$.
    **Significance**: Compact SVD discards basis vectors that do not contribute to the result, making it more computationally efficient.

**9. [Application] Why is SVD used for image compression?**

??? success "Solution"
    **Principle:**
    1. An image is a matrix $A$ of pixel values.
    2. Natural images have high local correlation, meaning singular values decay rapidly.
    3. By keeping only the top $k$ singular values (truncated SVD), we reconstruct a visually near-lossless image with very little data.
    4. This reduces storage from $mn$ to $(m+n)k$.

**10. [Condition] Express the condition number $\kappa(A)$ using singular values.**

??? success "Solution"
    **Conclusion:**
    For a non-singular square matrix, $\kappa(A) = \sigma_{\max} / \sigma_{\min}$.
    **Significance**: This shows that if a matrix has a singular value very close to 0, it is "ill-conditioned," making its inversion extremely sensitive to noise.

## Chapter Summary

SVD is the ultimate weapon for solving real-world linear algebra problems:

1.  **Universality**: It bridges the gap between square and rectangular, full-rank and deficient matrices, providing a unified analytic perspective for all linear mappings.
2.  **Energy Concentration**: By ordering singular values, SVD identifies the "most important" components of data, providing mathematical criteria for compression and dimensionality reduction.
3.  **Numerical Robustness**: As the basis for pseudoinverses and least squares solutions, SVD exhibits high stability against ill-conditioned matrices, serving as the final line of defense in engineering computations.
