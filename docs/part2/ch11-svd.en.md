# Chapter 11  Singular Value Decomposition

<div class="context-flow" markdown>

**Prerequisites**: Ch8 Orthogonality/Spectral theorem · Ch10 Schur/Spectral decomposition · **Chapter arc**: $\sigma_i = \sqrt{\lambda_i(A^TA)}$ → $A = U\Sigma V^T$ → Geometry (rotation × scaling × rotation) → **Eckart-Young** (optimal low-rank approximation) → Pseudoinverse/PCA/Condition number
Essence: SVD is the generalization of the spectral theorem to **arbitrary matrices** — simultaneously diagonalizing $A$ with two orthonormal bases

</div>

The Singular Value Decomposition (SVD) is one of the most important matrix decompositions in linear algebra. Unlike eigenvalue decomposition, SVD applies to **any** matrix — it requires neither a square matrix nor diagonalizability. SVD factors an $m \times n$ matrix into a product of two orthogonal matrices and a diagonal matrix, clearly revealing the geometric essence of a linear transformation: rotation, scaling, then rotation again. SVD has extremely broad applications in data science, signal processing, numerical computation, and other fields, making it a core tool of modern applied mathematics.

---

## 11.1 Definition of Singular Values

<div class="context-flow" markdown>

$A^TA$ positive semidefinite → eigenvalues $\lambda_i \ge 0$ → $\sigma_i = \sqrt{\lambda_i}$ → Singular values unify rank, spectral norm $\|A\|_2 = \sigma_1$, Frobenius norm

</div>

For any $m \times n$ real matrix $A$, the matrix $A^T A$ is an $n \times n$ symmetric positive semidefinite matrix, so all its eigenvalues are nonnegative real numbers. Singular values are extracted from these eigenvalues.

!!! definition "Definition 11.1 (Singular Value)"
    Let $A$ be an $m \times n$ real matrix, and let $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n \ge 0$ be the eigenvalues of $A^T A$. The **singular values** of $A$ are defined as

    $$
    \sigma_i = \sqrt{\lambda_i}, \quad i = 1, 2, \ldots, n.
    $$

    They are arranged in descending order $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_n \ge 0$.

!!! note "Note"
    Proof that $A^T A$ is positive semidefinite: for any $\mathbf{x} \in \mathbb{R}^n$,

    $$
    \mathbf{x}^T (A^T A) \mathbf{x} = (A\mathbf{x})^T (A\mathbf{x}) = \|A\mathbf{x}\|^2 \ge 0.
    $$

    Therefore all eigenvalues of $A^T A$ are nonnegative, and the definition of singular values is well-posed.

!!! theorem "Theorem 11.1 (Relationship between singular values and matrix norms)"
    Let $A$ be an $m \times n$ real matrix with singular values $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_n \ge 0$. Then:

    1. $\|A\|_2 = \sigma_1$ (the spectral norm equals the largest singular value);
    2. $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_n^2}$ (the Frobenius norm equals the square root of the sum of squared singular values);
    3. $\operatorname{rank}(A) = $ the number of nonzero singular values.

??? proof "Proof"
    **(1)** The spectral norm is defined as

    $$
    \|A\|_2 = \max_{\|\mathbf{x}\|=1} \|A\mathbf{x}\|.
    $$

    Since $\|A\mathbf{x}\|^2 = \mathbf{x}^T A^T A \mathbf{x}$, and $A^T A$ is a symmetric positive semidefinite matrix with eigenvalues $\lambda_1 \ge \cdots \ge \lambda_n \ge 0$ and corresponding orthonormal eigenvectors $\mathbf{v}_1, \ldots, \mathbf{v}_n$. For any unit vector $\mathbf{x} = \sum c_i \mathbf{v}_i$ ($\sum c_i^2 = 1$),

    $$
    \|A\mathbf{x}\|^2 = \sum_{i=1}^n \lambda_i c_i^2 \le \lambda_1 \sum c_i^2 = \lambda_1.
    $$

    Equality holds when $\mathbf{x} = \mathbf{v}_1$. Therefore $\|A\|_2 = \sqrt{\lambda_1} = \sigma_1$.

    **(2)** $\|A\|_F^2 = \operatorname{tr}(A^T A) = \sum_{i=1}^n \lambda_i = \sum_{i=1}^n \sigma_i^2$.

    **(3)** $\operatorname{rank}(A) = \operatorname{rank}(A^T A)$ (since $A\mathbf{x} = \mathbf{0}$ if and only if $A^T A \mathbf{x} = \mathbf{0}$), and the rank of a symmetric matrix equals the number of nonzero eigenvalues.

!!! definition "Definition 11.2 (Left and right singular vectors)"
    Let $A$ be an $m \times n$ real matrix.

    - The orthonormal eigenvectors $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^n$ of $A^T A$ are called the **right singular vectors** of $A$;
    - The orthonormal eigenvectors $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m \in \mathbb{R}^m$ of $A A^T$ are called the **left singular vectors** of $A$.

!!! theorem "Theorem 11.2 ($A^TA$ and $AA^T$ have the same nonzero eigenvalues)"
    Let $A$ be an $m \times n$ real matrix. Then $A^T A$ and $A A^T$ have the same nonzero eigenvalues (counting multiplicities).

??? proof "Proof"
    Let $\lambda \neq 0$ be an eigenvalue of $A^T A$ with eigenvector $\mathbf{v}$, so $A^T A \mathbf{v} = \lambda \mathbf{v}$. Let $\mathbf{u} = A\mathbf{v}$. Then $\mathbf{u} \neq \mathbf{0}$ (since if $A\mathbf{v} = \mathbf{0}$, then $\lambda \mathbf{v} = A^T A \mathbf{v} = \mathbf{0}$, contradicting $\lambda \neq 0$ and $\mathbf{v} \neq \mathbf{0}$), and

    $$
    A A^T \mathbf{u} = A A^T (A \mathbf{v}) = A (A^T A \mathbf{v}) = A (\lambda \mathbf{v}) = \lambda (A\mathbf{v}) = \lambda \mathbf{u}.
    $$

    Therefore $\lambda$ is also an eigenvalue of $A A^T$. The reverse direction is proved analogously.

!!! example "Example 11.1"
    Find the singular values of the matrix $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}$.

    **Solution:** Compute

    $$
    A^T A = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}.
    $$

    The characteristic polynomial is

    $$
    \det(A^T A - \lambda I) = (2 - \lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 3)(\lambda - 1).
    $$

    The eigenvalues are $\lambda_1 = 3$ and $\lambda_2 = 1$.

    Therefore the singular values are $\sigma_1 = \sqrt{3}$ and $\sigma_2 = 1$.

---

## 11.2 The SVD Theorem

<div class="context-flow" markdown>

Construction path: Apply spectral theorem to $A^TA$ to get $V$ → $\mathbf{u}_i = \frac{1}{\sigma_i}A\mathbf{v}_i$ gives $U$ → Verify $AV = U\Sigma$ → Outer product expansion $A = \sum \sigma_i \mathbf{u}_i\mathbf{v}_i^T$ (rank-one decomposition)

</div>

The SVD theorem is the central result of this chapter, guaranteeing that any matrix can be decomposed into a standard form.

<div class="context-flow" markdown>

**Insight**: The essence of SVD existence — apply the spectral theorem to $A^TA$ (a square, positive semidefinite matrix), then "map" the right singular vectors to left singular vectors through $A$

</div>

!!! theorem "Theorem 11.3 (SVD Theorem)"
    Let $A$ be an $m \times n$ real matrix with $\operatorname{rank}(A) = r$. Then there exist an $m \times m$ orthogonal matrix $U$, an $n \times n$ orthogonal matrix $V$, and an $m \times n$ matrix $\Sigma$ such that

    $$
    A = U \Sigma V^T,
    $$

    where $\Sigma$ has the form

    $$
    \Sigma = \begin{pmatrix} \sigma_1 & & & 0 & \cdots & 0 \\ & \sigma_2 & & 0 & \cdots & 0 \\ & & \ddots & \vdots & & \vdots \\ & & & \sigma_r & \cdots & 0 \\ 0 & \cdots & 0 & 0 & \cdots & 0 \\ \vdots & & \vdots & \vdots & & \vdots \\ 0 & \cdots & 0 & 0 & \cdots & 0 \end{pmatrix}_{m \times n},
    $$

    with $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$.

??? proof "Proof"
    **Constructive proof:**

    **Step 1:** Since $A^T A$ is an $n \times n$ symmetric positive semidefinite matrix, by the spectral theorem, there exists an orthogonal matrix $V = [\mathbf{v}_1, \ldots, \mathbf{v}_n]$ such that

    $$
    A^T A = V \operatorname{diag}(\lambda_1, \ldots, \lambda_n) V^T,
    $$

    where $\lambda_1 \ge \cdots \ge \lambda_r > 0 = \lambda_{r+1} = \cdots = \lambda_n$. Let $\sigma_i = \sqrt{\lambda_i}$ ($i = 1, \ldots, r$).

    **Step 2:** For $i = 1, \ldots, r$, define

    $$
    \mathbf{u}_i = \frac{1}{\sigma_i} A \mathbf{v}_i.
    $$

    Verify that $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ is an orthonormal set:

    $$
    \mathbf{u}_i^T \mathbf{u}_j = \frac{1}{\sigma_i \sigma_j} \mathbf{v}_i^T A^T A \mathbf{v}_j = \frac{1}{\sigma_i \sigma_j} \mathbf{v}_i^T (\lambda_j \mathbf{v}_j) = \frac{\lambda_j}{\sigma_i \sigma_j} \delta_{ij} = \delta_{ij}.
    $$

    **Step 3:** Extend $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ to an orthonormal basis $\{\mathbf{u}_1, \ldots, \mathbf{u}_m\}$ of $\mathbb{R}^m$, and let $U = [\mathbf{u}_1, \ldots, \mathbf{u}_m]$.

    **Step 4:** Verify $A = U \Sigma V^T$. By construction, $A\mathbf{v}_i = \sigma_i \mathbf{u}_i$ ($i \le r$), and for $i > r$, $A\mathbf{v}_i = \mathbf{0}$ (since $\|A\mathbf{v}_i\|^2 = \mathbf{v}_i^T A^T A \mathbf{v}_i = \lambda_i = 0$). Therefore

    $$
    AV = U\Sigma, \quad \text{i.e.,} \quad A = U\Sigma V^T. \qquad \blacksquare
    $$

!!! definition "Definition 11.3 (Components of the SVD)"
    In the SVD $A = U\Sigma V^T$:

    - $U = [\mathbf{u}_1, \ldots, \mathbf{u}_m]$ is an $m \times m$ orthogonal matrix whose columns are the **left singular vectors**;
    - $V = [\mathbf{v}_1, \ldots, \mathbf{v}_n]$ is an $n \times n$ orthogonal matrix whose columns are the **right singular vectors**;
    - $\Sigma$ is an $m \times n$ diagonal matrix with diagonal entries being the singular values $\sigma_1 \ge \cdots \ge \sigma_r > 0$.

!!! theorem "Theorem 11.4 (Outer product expansion of SVD)"
    Let $A = U\Sigma V^T$ be the SVD of $A$ with $\operatorname{rank}(A) = r$. Then

    $$
    A = \sum_{i=1}^{r} \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$

    That is, $A$ can be expressed as a weighted sum of $r$ rank-one matrices.

??? proof "Proof"
    From $A = U\Sigma V^T$, expanding gives

    $$
    A = \sum_{i=1}^{\min(m,n)} \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$

    Since $\sigma_i = 0$ for $i > r$, the sum effectively runs only up to $r$. Each $\mathbf{u}_i \mathbf{v}_i^T$ is an $m \times n$ matrix with $\operatorname{rank}(\mathbf{u}_i \mathbf{v}_i^T) = 1$.

!!! example "Example 11.2"
    Find the SVD of the matrix $A = \begin{pmatrix} 3 & 0 \\ 0 & 2 \end{pmatrix}$.

    **Solution:**

    $$
    A^T A = \begin{pmatrix} 9 & 0 \\ 0 & 4 \end{pmatrix},
    $$

    eigenvalues $\lambda_1 = 9$, $\lambda_2 = 4$; eigenvectors $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$.

    Singular values $\sigma_1 = 3$, $\sigma_2 = 2$.

    Left singular vectors: $\mathbf{u}_1 = \frac{1}{3}A\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$, $\mathbf{u}_2 = \frac{1}{2}A\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$.

    Therefore

    $$
    A = \begin{pmatrix}1&0\\0&1\end{pmatrix}\begin{pmatrix}3&0\\0&2\end{pmatrix}\begin{pmatrix}1&0\\0&1\end{pmatrix} = I \cdot \Sigma \cdot I.
    $$

    A diagonal matrix is its own SVD ($U = V = I$).

---

## 11.3 Geometric Meaning of SVD

<div class="context-flow" markdown>

$V^T$ (rotation) → $\Sigma$ (axis-aligned scaling) → $U$ (rotation): unit sphere → hyperellipsoid → semi-axis lengths = $\sigma_i$, directions = $\mathbf{u}_i$

</div>

The geometric interpretation of SVD is extremely intuitive: any linear transformation can be decomposed into **rotation (or reflection)**, **scaling**, and another **rotation (or reflection)**.

!!! definition "Definition 11.4 (Geometric meaning of orthogonal transformations)"
    The linear transformation represented by an orthogonal matrix $Q$ ($Q^TQ = I$) is isometric, i.e., $\|Q\mathbf{x}\| = \|\mathbf{x}\|$. When $\det Q = 1$ it is a rotation; when $\det Q = -1$ it is a rotation composed with a reflection.

!!! theorem "Theorem 11.5 (Geometric decomposition of a linear transformation)"
    Let $A = U\Sigma V^T$ be the SVD of an $m \times n$ matrix $A$. Then the linear transformation $\mathbf{x} \mapsto A\mathbf{x}$ can be decomposed into three steps:

    1. $V^T$: Rotate in $\mathbb{R}^n$ (transforming the standard basis to the right singular vector directions);
    2. $\Sigma$: Scale along coordinate axes by $\sigma_i$ (from $\mathbb{R}^n$ to $\mathbb{R}^m$);
    3. $U$: Rotate in $\mathbb{R}^m$ (transforming coordinate axis directions to the left singular vector directions).

??? proof "Proof"
    For any $\mathbf{x} \in \mathbb{R}^n$, let $\mathbf{y} = V^T \mathbf{x}$. Since $V$ is orthogonal, $\|\mathbf{y}\| = \|\mathbf{x}\|$, so the first step is an isometry (rotation or reflection).

    In the second step, $\mathbf{z} = \Sigma \mathbf{y}$ multiplies the $i$-th component of $\mathbf{y}$ by $\sigma_i$, which is axis-aligned scaling.

    In the third step, $A\mathbf{x} = U\mathbf{z}$, and $U$ is orthogonal, so this is also an isometry.

    Combining the three steps, the unit sphere $\{\mathbf{x} : \|\mathbf{x}\| = 1\}$ remains a unit sphere under $V^T$, becomes a hyperellipsoid with semi-axis lengths $\sigma_i$ under $\Sigma$, and is rotated to its final position under $U$.

!!! note "Note"
    In the two-dimensional case, the unit circle $\|\mathbf{x}\| = 1$ is transformed by $A$ into an ellipse, whose semi-major and semi-minor axis lengths are exactly $\sigma_1$ and $\sigma_2$, with directions $\mathbf{u}_1$ and $\mathbf{u}_2$ respectively.

!!! example "Example 11.3"
    Let $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$. Describe the action of $A$ on the unit circle.

    **Solution:** $A$ is already a diagonal matrix, and its SVD is $A = I \cdot \begin{pmatrix}2&0\\0&1\end{pmatrix} \cdot I$.

    The unit circle $x^2 + y^2 = 1$ is transformed by $A$ into the ellipse $\frac{x^2}{4} + y^2 = 1$, with semi-major axis $\sigma_1 = 2$ (along the $x$-axis) and semi-minor axis $\sigma_2 = 1$ (along the $y$-axis).

---

## 11.4 Computing the SVD

<div class="context-flow" markdown>

Standard procedure: $A^TA$ → eigenvalues/eigenvectors → singular values → construct left singular vectors via $\mathbf{u}_i = \frac{1}{\sigma_i}A\mathbf{v}_i$ → assemble

</div>

This section demonstrates the steps of computing the SVD through detailed examples.

!!! definition "Definition 11.5 (Standard procedure for computing the SVD)"
    Given an $m \times n$ matrix $A$, the SVD is computed as follows:

    1. Compute $A^T A$;
    2. Find the eigenvalues $\lambda_1 \ge \cdots \ge \lambda_n \ge 0$ and corresponding orthonormal eigenvectors $\mathbf{v}_1, \ldots, \mathbf{v}_n$ of $A^T A$;
    3. Singular values $\sigma_i = \sqrt{\lambda_i}$;
    4. For each $\sigma_i > 0$, compute $\mathbf{u}_i = \frac{1}{\sigma_i} A \mathbf{v}_i$;
    5. Extend $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ to an orthonormal basis of $\mathbb{R}^m$;
    6. Assemble $U$, $\Sigma$, $V$.

!!! example "Example 11.4"
    Find the SVD of the matrix $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \\ 0 & 0 \end{pmatrix}$.

    **Solution:**

    **Step 1:** Compute $A^T A$:

    $$
    A^T A = \begin{pmatrix} 1&1&0 \\ 1&-1&0 \end{pmatrix}\begin{pmatrix} 1&1 \\ 1&-1 \\ 0&0 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}.
    $$

    **Step 2:** Eigenvalues $\lambda_1 = \lambda_2 = 2$; take orthogonal eigenvectors $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$, $\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$.

    **Step 3:** Singular values $\sigma_1 = \sigma_2 = \sqrt{2}$.

    **Step 4:** Compute left singular vectors:

    $$
    \mathbf{u}_1 = \frac{1}{\sqrt{2}} A \mathbf{v}_1 = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\1\\0\end{pmatrix}, \quad
    \mathbf{u}_2 = \frac{1}{\sqrt{2}} A \mathbf{v}_2 = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-1\\0\end{pmatrix}.
    $$

    **Step 5:** Extend $\{\mathbf{u}_1, \mathbf{u}_2\}$; take $\mathbf{u}_3 = \begin{pmatrix}0\\0\\1\end{pmatrix}$.

    **Step 6:** Therefore

    $$
    A = \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 0 \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 0 \\ 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} \sqrt{2} & 0 \\ 0 & \sqrt{2} \\ 0 & 0 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}^T.
    $$

!!! example "Example 11.5"
    Find the SVD of the matrix $A = \begin{pmatrix} 4 & 0 \\ 3 & -5 \end{pmatrix}$.

    **Solution:**

    **Step 1:**

    $$
    A^T A = \begin{pmatrix} 4&3 \\ 0&-5 \end{pmatrix}\begin{pmatrix} 4&0 \\ 3&-5 \end{pmatrix} = \begin{pmatrix} 25 & -15 \\ -15 & 25 \end{pmatrix}.
    $$

    **Step 2:** Characteristic polynomial:

    $$
    \det(A^T A - \lambda I) = (25-\lambda)^2 - 225 = \lambda^2 - 50\lambda + 400 = (\lambda - 40)(\lambda - 10).
    $$

    Eigenvalues $\lambda_1 = 40$, $\lambda_2 = 10$.

    For $\lambda_1 = 40$: $(A^T A - 40I)\mathbf{v} = \mathbf{0}$ gives $\begin{pmatrix}-15&-15\\-15&-15\end{pmatrix}\mathbf{v} = \mathbf{0}$, yielding $\mathbf{v}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\-1\end{pmatrix}$.

    For $\lambda_2 = 10$: $(A^T A - 10I)\mathbf{v} = \mathbf{0}$ gives $\begin{pmatrix}15&-15\\-15&15\end{pmatrix}\mathbf{v} = \mathbf{0}$, yielding $\mathbf{v}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$.

    **Step 3:** $\sigma_1 = \sqrt{40} = 2\sqrt{10}$, $\sigma_2 = \sqrt{10}$.

    **Step 4:**

    $$
    \mathbf{u}_1 = \frac{1}{2\sqrt{10}} A \mathbf{v}_1 = \frac{1}{2\sqrt{10}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}4\\8\end{pmatrix} = \frac{1}{\sqrt{5}}\begin{pmatrix}1\\2\end{pmatrix},
    $$

    $$
    \mathbf{u}_2 = \frac{1}{\sqrt{10}} A \mathbf{v}_2 = \frac{1}{\sqrt{10}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}4\\-2\end{pmatrix} = \frac{1}{\sqrt{5}}\begin{pmatrix}2\\-1\end{pmatrix}.
    $$

    **Step 5:** Assemble:

    $$
    A = \frac{1}{\sqrt{5}}\begin{pmatrix}1&2\\2&-1\end{pmatrix} \begin{pmatrix}2\sqrt{10}&0\\0&\sqrt{10}\end{pmatrix} \frac{1}{\sqrt{2}}\begin{pmatrix}1&1\\-1&1\end{pmatrix}^T.
    $$

---

## 11.5 Compact SVD and Truncated SVD

<div class="context-flow" markdown>

Compact SVD: remove zero singular values → $A = U_r\Sigma_r V_r^T$ · Truncated SVD: keep only the top $k$ → $A_k$ is a rank-$k$ approximation, error $\|A - A_k\|_2 = \sigma_{k+1}$

</div>

When the rank of a matrix is much smaller than its dimensions, the full SVD contains a large amount of redundant information. The compact SVD and truncated SVD provide more economical representations.

!!! definition "Definition 11.6 (Compact SVD)"
    Let $A = U\Sigma V^T$ be the full SVD of an $m \times n$ matrix $A$ with $\operatorname{rank}(A) = r$. Let $U_r = [\mathbf{u}_1, \ldots, \mathbf{u}_r]$ ($m \times r$), $\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$ ($r \times r$), $V_r = [\mathbf{v}_1, \ldots, \mathbf{v}_r]$ ($n \times r$). Then

    $$
    A = U_r \Sigma_r V_r^T
    $$

    is called the **compact SVD** (or **economy SVD**) of $A$.

!!! definition "Definition 11.7 (Truncated SVD)"
    Let $\operatorname{rank}(A) = r$. For $k < r$, let

    $$
    A_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T = U_k \Sigma_k V_k^T,
    $$

    where $U_k$, $\Sigma_k$, $V_k$ consist of the first $k$ columns/rows respectively. $A_k$ is called the **rank-$k$ truncated SVD** of $A$.

!!! theorem "Theorem 11.6 (Approximation error of the truncated SVD)"
    Let the singular values of $A$ be $\sigma_1 \ge \cdots \ge \sigma_r > 0$, and let $A_k$ be the rank-$k$ truncated approximation. Then

    $$
    \|A - A_k\|_2 = \sigma_{k+1}, \qquad \|A - A_k\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}.
    $$

??? proof "Proof"
    From $A = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ and $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$, we get

    $$
    A - A_k = \sum_{i=k+1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$

    This is itself the SVD of $A - A_k$ (with singular values $\sigma_{k+1}, \ldots, \sigma_r$), so

    $$
    \|A - A_k\|_2 = \sigma_{k+1}, \qquad \|A - A_k\|_F = \sqrt{\sum_{i=k+1}^r \sigma_i^2}. \qquad \blacksquare
    $$

!!! example "Example 11.6"
    Suppose the singular values of $A$ are $10, 5, 2, 0.1$. What are the spectral norm error and the Frobenius norm error when using the rank-2 truncated approximation? What is the relative Frobenius error?

    **Solution:**

    $$
    \|A - A_2\|_2 = \sigma_3 = 2, \qquad \|A - A_2\|_F = \sqrt{2^2 + 0.1^2} = \sqrt{4.01} \approx 2.0025.
    $$

    $$
    \|A\|_F = \sqrt{100 + 25 + 4 + 0.01} = \sqrt{129.01} \approx 11.359.
    $$

    The relative error is $\frac{\|A - A_2\|_F}{\|A\|_F} \approx \frac{2.0025}{11.359} \approx 17.6\%$.

---

## 11.6 Low-Rank Approximation

<div class="context-flow" markdown>

**Eckart-Young-Mirsky**: the truncated SVD is the **best approximation** among all matrices of rank $\le k$ (optimal under both spectral and Frobenius norms)

</div>

SVD provides the optimal low-rank approximation, a result precisely stated by the Eckart-Young-Mirsky theorem.

<div class="context-flow" markdown>

**Insight**: The dimension argument in the proof $\dim W + \dim\ker(B) > n$ is analogous to the inertia theorem in Ch9 — "two large subspaces must intersect"

</div>

!!! theorem "Theorem 11.7 (Eckart-Young-Mirsky theorem)"
    Let $A$ be an $m \times n$ real matrix with singular values $\sigma_1 \ge \cdots \ge \sigma_r > 0$. For any matrix $B$ of rank at most $k$ ($k < r$),

    $$
    \|A - A_k\|_2 \le \|A - B\|_2, \qquad \|A - A_k\|_F \le \|A - B\|_F,
    $$

    where $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ is the rank-$k$ truncated SVD of $A$. That is, $A_k$ is the **best rank-$k$ approximation** of $A$ under both spectral and Frobenius norms.

??? proof "Proof"
    **Proof for the spectral norm case:**

    Let $B$ be any matrix of rank at most $k$. The dimension of $\ker(B)$ is at least $n - k$. Consider the subspace $W = \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{k+1}\}$ with $\dim W = k + 1$.

    By the dimension formula, $W \cap \ker(B) \neq \{\mathbf{0}\}$ (since $(k+1) + (n-k) = n + 1 > n$). Take a unit vector $\mathbf{w} \in W \cap \ker(B)$, so $B\mathbf{w} = \mathbf{0}$, and therefore

    $$
    \|A - B\|_2^2 \ge \|(A-B)\mathbf{w}\|^2 = \|A\mathbf{w}\|^2.
    $$

    Writing $\mathbf{w} = \sum_{i=1}^{k+1} c_i \mathbf{v}_i$ ($\sum c_i^2 = 1$), we get

    $$
    \|A\mathbf{w}\|^2 = \left\|\sum_{i=1}^{k+1} c_i \sigma_i \mathbf{u}_i\right\|^2 = \sum_{i=1}^{k+1} c_i^2 \sigma_i^2 \ge \sigma_{k+1}^2 \sum c_i^2 = \sigma_{k+1}^2.
    $$

    Therefore $\|A - B\|_2 \ge \sigma_{k+1} = \|A - A_k\|_2$. $\blacksquare$

!!! example "Example 11.7"
    Let $A = \begin{pmatrix} 3 & 2 & 2 \\ 2 & 3 & -2 \end{pmatrix}$. Find the best rank-one approximation of $A$.

    **Solution:** Compute

    $$
    A A^T = \begin{pmatrix} 17 & 8 \\ 8 & 17 \end{pmatrix}.
    $$

    Eigenvalues $\lambda_1 = 25$, $\lambda_2 = 9$, so $\sigma_1 = 5$, $\sigma_2 = 3$.

    For $\lambda_1 = 25$: $\mathbf{u}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$.

    The corresponding right singular vector:

    $$
    \mathbf{v}_1 = \frac{1}{\sigma_1} A^T \mathbf{u}_1 = \frac{1}{5} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}5\\5\\0\end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\\0\end{pmatrix}.
    $$

    The best rank-one approximation is

    $$
    A_1 = \sigma_1 \mathbf{u}_1 \mathbf{v}_1^T = 5 \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1&1&0\end{pmatrix} = \begin{pmatrix} \frac{5}{2} & \frac{5}{2} & 0 \\ \frac{5}{2} & \frac{5}{2} & 0 \end{pmatrix}.
    $$

    Approximation error $\|A - A_1\|_2 = \sigma_2 = 3$.

---

## 11.7 The Moore-Penrose Pseudoinverse

<div class="context-flow" markdown>

$A^+ = V\Sigma^+ U^T$: invert the nonzero singular values in the SVD → gives the **minimum-norm least-squares solution** → unifies Ch8 orthogonal projection and least-squares theory

</div>

For matrices that are noninvertible or nonsquare, the Moore-Penrose pseudoinverse provides the concept of a "best inverse," which can be concisely constructed via SVD.

!!! definition "Definition 11.8 (Moore-Penrose pseudoinverse)"
    Let $A$ be an $m \times n$ real matrix. The **Moore-Penrose pseudoinverse** $A^+$ is the unique $n \times m$ matrix satisfying the following four conditions:

    1. $A A^+ A = A$;
    2. $A^+ A A^+ = A^+$;
    3. $(A A^+)^T = A A^+$ (i.e., $A A^+$ is symmetric);
    4. $(A^+ A)^T = A^+ A$ (i.e., $A^+ A$ is symmetric).

!!! theorem "Theorem 11.8 (Existence and uniqueness of the Moore-Penrose pseudoinverse)"
    For any $m \times n$ real matrix $A$, the matrix $A^+$ satisfying the four conditions of Definition 11.8 exists and is unique.

??? proof "Proof"
    **Existence:** Let $A = U\Sigma V^T$ be the SVD of $A$ with $\operatorname{rank}(A) = r$. Define

    $$
    \Sigma^+ = \begin{pmatrix} \sigma_1^{-1} & & \\ & \ddots & \\ & & \sigma_r^{-1} \\ & \mathbf{0} & \end{pmatrix}_{n \times m},
    $$

    and let $A^+ = V \Sigma^+ U^T$. We verify the four conditions.

    (1) $A A^+ A = (U\Sigma V^T)(V\Sigma^+ U^T)(U\Sigma V^T) = U\Sigma\Sigma^+\Sigma V^T = U\Sigma V^T = A$.

    Here $\Sigma\Sigma^+\Sigma = \Sigma$ holds because the diagonal entries satisfy $\sigma_i \cdot \sigma_i^{-1} \cdot \sigma_i = \sigma_i$ ($i \le r$), and the rest are zero.

    (2)-(4) are verified similarly.

    **Uniqueness:** Suppose $B_1, B_2$ both satisfy the four conditions. Using conditions (1)(3) one can show $AB_1 = AB_2$, then using (1)(4) one shows $B_1A = B_2A$, and finally by (2), $B_1 = B_1AB_1 = B_2AB_2 = B_2$. $\blacksquare$

!!! definition "Definition 11.9 (SVD formula for the pseudoinverse)"
    Let $A = U\Sigma V^T$ be the SVD of $A$ with $\operatorname{rank}(A) = r$. Then

    $$
    A^+ = V \Sigma^+ U^T = \sum_{i=1}^{r} \frac{1}{\sigma_i} \mathbf{v}_i \mathbf{u}_i^T.
    $$

!!! theorem "Theorem 11.9 (Pseudoinverse and least squares)"
    Let $A$ be an $m \times n$ matrix and $\mathbf{b} \in \mathbb{R}^m$. Then $\mathbf{x}^* = A^+ \mathbf{b}$ is the **minimum-norm least-squares solution** of the linear system $A\mathbf{x} = \mathbf{b}$, i.e., among all $\mathbf{x}$ that minimize $\|A\mathbf{x} - \mathbf{b}\|$, $\mathbf{x}^*$ has the smallest norm $\|\mathbf{x}^*\|$.

??? proof "Proof"
    Let $A = U\Sigma V^T$ with $\operatorname{rank}(A) = r$. Let $\mathbf{c} = U^T \mathbf{b}$ and $\mathbf{y} = V^T \mathbf{x}$. Since $U, V$ are orthogonal,

    $$
    \|A\mathbf{x} - \mathbf{b}\|^2 = \|\Sigma\mathbf{y} - \mathbf{c}\|^2 = \sum_{i=1}^r (\sigma_i y_i - c_i)^2 + \sum_{i=r+1}^m c_i^2.
    $$

    The second term is independent of $\mathbf{x}$. Minimizing the first term gives $y_i = c_i / \sigma_i$ ($i = 1, \ldots, r$).

    Under this constraint, $\|\mathbf{x}\|^2 = \|\mathbf{y}\|^2 = \sum_{i=1}^r y_i^2 + \sum_{i=r+1}^n y_i^2$. Minimizing $\|\mathbf{x}\|$ requires $y_i = 0$ ($i = r+1, \ldots, n$).

    Therefore the minimum-norm least-squares solution is $\mathbf{y}^* = (\frac{c_1}{\sigma_1}, \ldots, \frac{c_r}{\sigma_r}, 0, \ldots, 0)^T$, i.e., $\mathbf{x}^* = V\mathbf{y}^* = V\Sigma^+ U^T \mathbf{b} = A^+ \mathbf{b}$. $\blacksquare$

!!! example "Example 11.8"
    Let $A = \begin{pmatrix}1\\2\end{pmatrix}$. Find $A^+$ and solve $A\mathbf{x} = \begin{pmatrix}3\\4\end{pmatrix}$ in the least-squares sense.

    **Solution:** $A^T A = (5)$, $\sigma_1 = \sqrt{5}$. $\mathbf{v}_1 = (1)$, $\mathbf{u}_1 = \frac{1}{\sqrt{5}}\begin{pmatrix}1\\2\end{pmatrix}$.

    $$
    A^+ = \frac{1}{\sigma_1} \mathbf{v}_1 \mathbf{u}_1^T = \frac{1}{\sqrt{5}} \cdot (1) \cdot \frac{1}{\sqrt{5}}\begin{pmatrix}1&2\end{pmatrix} = \frac{1}{5}\begin{pmatrix}1&2\end{pmatrix} = \begin{pmatrix}\frac{1}{5}&\frac{2}{5}\end{pmatrix}.
    $$

    Least-squares solution:

    $$
    x^* = A^+ \mathbf{b} = \frac{1}{5}(1 \cdot 3 + 2 \cdot 4) = \frac{11}{5} = 2.2.
    $$

    Verification: $A x^* = \begin{pmatrix}2.2\\4.4\end{pmatrix}$, residual $\mathbf{b} - Ax^* = \begin{pmatrix}0.8\\-0.4\end{pmatrix}$, which is orthogonal to the column space of $A$: $A^T(\mathbf{b} - Ax^*) = (1)(0.8) + (2)(-0.4) = 0$.

---

## 11.8 Applications of SVD

<div class="context-flow" markdown>

Data compression (low-rank approximation) · PCA (columns of $V$ = principal component directions) · Condition number $\kappa = \sigma_1/\sigma_r$ (measures ill-conditioning) → SVD is the Swiss army knife of applied mathematics

</div>

SVD has important applications in numerous fields. This section introduces several typical scenarios.

### 11.8.1 Data Compression

!!! definition "Definition 11.10 (SVD image compression)"
    Viewing an $m \times n$ grayscale image as a matrix $A$, its rank-$k$ truncated SVD approximation $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ can be used for image compression. Storing $A_k$ requires $k(m + n + 1)$ numbers, while the original matrix requires $mn$ numbers. When $k \ll \min(m,n)$, the compression ratio is

    $$
    \rho = \frac{k(m + n + 1)}{mn}.
    $$

!!! example "Example 11.9"
    For a $1000 \times 800$ grayscale image, what is the compression ratio when using a rank-$k = 50$ truncated SVD approximation?

    **Solution:**

    $$
    \rho = \frac{50 \times (1000 + 800 + 1)}{1000 \times 800} = \frac{50 \times 1801}{800000} = \frac{90050}{800000} \approx 11.3\%.
    $$

    That is, only about $11.3\%$ of the original data needs to be stored.

### 11.8.2 Introduction to Principal Component Analysis

!!! theorem "Theorem 11.10 (Relationship between SVD and PCA)"
    Let the data matrix $X$ ($n \times p$, centered) have SVD $X = U\Sigma V^T$. Then:

    1. The covariance matrix of $X$ is $S = \frac{1}{n-1}X^TX = \frac{1}{n-1}V\Sigma^2 V^T$;
    2. The columns of $V$ are the principal component directions (principal axes);
    3. The variance of the $i$-th principal component is $\frac{\sigma_i^2}{n-1}$;
    4. The principal component score matrix is $XV = U\Sigma$.

??? proof "Proof"
    From $X = U\Sigma V^T$, we get $X^T X = V\Sigma^T U^T U \Sigma V^T = V\Sigma^T\Sigma V^T = V\Sigma^2 V^T$ (where $\Sigma^2$ denotes $\Sigma^T\Sigma$, a $p \times p$ diagonal matrix). Therefore

    $$
    S = \frac{1}{n-1}X^T X = V \left(\frac{\Sigma^2}{n-1}\right) V^T.
    $$

    This is precisely the spectral decomposition of the covariance matrix. The columns of $V$ are the eigenvectors of $S$ (principal component directions), with corresponding eigenvalues $\frac{\sigma_i^2}{n-1}$ (the variance of the $i$-th principal component).

    Principal component scores $Z = XV = U\Sigma V^T V = U\Sigma$. $\blacksquare$

### 11.8.3 Condition Number

!!! definition "Definition 11.11 (Condition Number)"
    Let $A$ be an $m \times n$ matrix with $\operatorname{rank}(A) = r$. The **condition number** of $A$ is defined as

    $$
    \kappa(A) = \frac{\sigma_1}{\sigma_r} = \frac{\sigma_{\max}}{\sigma_{\min}},
    $$

    where $\sigma_1$ and $\sigma_r$ are the largest and smallest nonzero singular values of $A$, respectively.

!!! theorem "Theorem 11.11 (Condition number and numerical stability)"
    For the linear system $A\mathbf{x} = \mathbf{b}$ ($A$ invertible), if $\mathbf{b}$ has a relative perturbation $\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}$, then the relative perturbation of the solution satisfies

    $$
    \frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \kappa(A) \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.
    $$

    The larger the condition number, the more **ill-conditioned** the problem; when the condition number is close to 1, the problem is **well-conditioned**.

??? proof "Proof"
    Let $A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$, so $A\delta\mathbf{x} = \delta\mathbf{b}$ and $\delta\mathbf{x} = A^{-1}\delta\mathbf{b}$.

    $$
    \|\delta\mathbf{x}\| = \|A^{-1}\delta\mathbf{b}\| \le \|A^{-1}\|_2 \|\delta\mathbf{b}\|.
    $$

    Also $\|\mathbf{b}\| = \|A\mathbf{x}\| \le \|A\|_2 \|\mathbf{x}\|$, so $\frac{1}{\|\mathbf{x}\|} \le \frac{\|A\|_2}{\|\mathbf{b}\|}$. Therefore

    $$
    \frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \|A\|_2 \|A^{-1}\|_2 \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = \frac{\sigma_1}{\sigma_n} \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = \kappa(A) \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}. \qquad \blacksquare
    $$

!!! example "Example 11.10"
    Let $A = \begin{pmatrix}1&1\\1&1.0001\end{pmatrix}$. Find its condition number and analyze the numerical stability.

    **Solution:**

    $$
    A^T A = \begin{pmatrix}2&2.0001\\2.0001&2.00020001\end{pmatrix}.
    $$

    The eigenvalues are approximately $\lambda_1 \approx 4.00020001$ and $\lambda_2 \approx 0.000000005$ (exact values can be computed).

    Rough estimates: $\sigma_1 \approx 2.00005$, $\sigma_2 \approx 0.00005$.

    $$
    \kappa(A) = \frac{\sigma_1}{\sigma_2} \approx \frac{2.00005}{0.00005} \approx 40000.
    $$

    The condition number is approximately $40000$, indicating that this matrix is ill-conditioned. A small perturbation in the right-hand side $\mathbf{b}$ can cause a large change in the solution.

---

## Chapter Summary

This chapter systematically introduced the Singular Value Decomposition (SVD), including:

1. **Singular values** come from the square roots of the eigenvalues of $A^T A$;
2. The **SVD theorem** guarantees that any matrix satisfies $A = U\Sigma V^T$;
3. **Geometric meaning**: a linear transformation = rotation + scaling + rotation;
4. The **Eckart-Young-Mirsky theorem**: the truncated SVD gives the optimal low-rank approximation;
5. The **Moore-Penrose pseudoinverse** $A^+ = V\Sigma^+ U^T$ gives the minimum-norm least-squares solution;
6. **Applications** span data compression, PCA, condition number analysis, and other fields.

SVD perfectly unifies the theoretical elegance of linear algebra with the effectiveness of practical computation, and it is the cornerstone for subsequent study of numerical linear algebra, machine learning, and other courses.
