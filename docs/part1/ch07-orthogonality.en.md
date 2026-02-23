# Chapter 07: Orthogonality and Least Squares

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05)

**Chapter Outline**: Inner Products and Norms → Definition of Orthogonality → Orthogonal Sets and Orthonormal Bases → Orthogonal Projection Matrices → The Gram-Schmidt Process → QR Decomposition → Method of Least Squares → Normal Equations → Properties of Orthogonal Matrices → Best Approximation Theorem

**Extension**: Orthogonality introduces metrics (length and angle) into abstract vector spaces; it is the cornerstone of signal processing (Wavelets, Fourier), numerical stability, and the SVD (Ch11).

</div>

If linear independence is the "skeleton" of a space, then orthogonality is its "ruler." Orthogonality not only simplifies coordinate representations but also solves "unsolvable" systems resulting from data noise through the mechanism of projection. This chapter demonstrates how to transform any basis into a perfect computational basis through orthogonalization.

---

## 07.1 Inner Product, Norm, and Orthogonality

!!! definition "Definition 07.1 (Dot Product and Norm)"
    For vectors $\mathbf{u}, \mathbf{v}$ in $\mathbb{R}^n$, the **dot product** (inner product) is $\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T \mathbf{v}$. The **length** (norm) of a vector is $\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}$.

!!! definition "Definition 07.2 (Orthogonality)"
    Two vectors $\mathbf{u}, \mathbf{v}$ are **orthogonal** if $\mathbf{u} \cdot \mathbf{v} = 0$.
    **Pythagorean Theorem**: $\mathbf{u}, \mathbf{v}$ are orthogonal $\iff$ $\|\mathbf{u} + \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$.

---

## 07.2 Orthogonal Projections and Gram-Schmidt

!!! theorem "Theorem 07.1 (Best Approximation Theorem)"
    Let $W$ be a subspace of $V$. For any vector $\mathbf{y}$ in $V$, the orthogonal projection $\hat{\mathbf{y}} = \operatorname{proj}_W \mathbf{y}$ is the vector in $W$ closest to $\mathbf{y}$.

!!! algorithm "Algorithm 07.1 (The Gram-Schmidt Process)"
    Given a basis $\{\mathbf{x}_1, \ldots, \mathbf{x}_n\}$, construct an orthonormal basis $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$:
    1.  $\mathbf{v}_1 = \mathbf{x}_1$
    2.  $\mathbf{v}_2 = \mathbf{x}_2 - \frac{\mathbf{x}_2 \cdot \mathbf{v}_1}{\mathbf{v}_1 \cdot \mathbf{v}_1} \mathbf{v}_1$
    3.  Continue similarly, then normalize: $\mathbf{q}_i = \mathbf{v}_i / \|\mathbf{v}_i\|$.

---

## 07.3 QR Decomposition

!!! theorem "Theorem 07.2 (QR Decomposition)"
    If the columns of an $m \times n$ matrix $A$ are linearly independent, then $A$ can be factored as $A = QR$.
    - $Q$ is an $m \times n$ matrix whose columns form an orthonormal basis for $C(A)$.
    - $R$ is an $n \times n$ upper triangular invertible matrix.

---

## 07.4 Least Squares Problems

!!! definition "Definition 07.3 (Least Squares Solution)"
    For the system $Ax = \mathbf{b}$ that has no solution, we seek $\hat{x}$ that minimizes $\|\mathbf{b} - A\hat{x}\|$.
    $\hat{x}$ satisfies the **Normal Equations**:
    $$A^T A \hat{x} = A^T \mathbf{b}$$

---

## Exercises

1. **[Dot Product] Determine if $(1, 2, 3)$ and $(1, 0, -1)$ are orthogonal.**

   ??? success "Solution"
       $1(1) + 2(0) + 3(-1) = -2 \neq 0$. Thus, they are not orthogonal.

2. **[Unit] Normalize the vector $(3, 4)$.**

   ??? success "Solution"
       Length is $\sqrt{3^2+4^2}=5$. The normalized vector is $(0.6, 0.8)$.

3. **[Projection] Find the projection of $(1, 2)$ onto the direction of $(1, 0)$.**

   ??? success "Solution"
       $\hat{y} = \frac{(1,2)\cdot(1,0)}{(1,0)\cdot(1,0)}(1,0) = (1, 0)$.

4. **[Orthogonal Matrix] Prove that an orthogonal matrix $Q$ preserves length: $\|Q\mathbf{x}\| = \|\mathbf{x}\|$.**

   ??? success "Solution"
       $\|Q\mathbf{x}\|^2 = (Q\mathbf{x})^T(Q\mathbf{x}) = \mathbf{x}^T Q^T Q \mathbf{x} = \mathbf{x}^T I \mathbf{x} = \|\mathbf{x}\|^2$.

5. **[Gram-Schmidt] Orthogonalize $(1, 1)$ and $(0, 1)$.**

   ??? success "Solution"
       $v_1 = (1, 1)$. $v_2 = (0, 1) - \frac{1}{2}(1, 1) = (-0.5, 0.5)$. Normalizing gives $(\frac{1}{\sqrt{2}}, \frac{1}{\sqrt{2}})$ and $(-\frac{1}{\sqrt{2}}, \frac{1}{\sqrt{2}})$.

6. **[Least Squares] Why is $A^T A$ invertible when $A$ has full column rank?**

   ??? success "Solution"
       If $A^T A x = 0 \implies x^T A^T A x = 0 \implies \|Ax\|^2 = 0 \implies Ax=0$. Since $A$ has full column rank, $x=0$. Thus $A^T A$ is non-singular.

7. **[QR] If $A=QR$, what is $A^T A$?**

   ??? success "Solution"
       $(QR)^T(QR) = R^T Q^T Q R = R^T R$.

8. **[Determinant] Prove the determinant of an orthogonal matrix must be $\pm 1$.**

   ??? success "Solution"
       $\det(Q^T Q) = \det(I) = 1 \implies \det(Q)^2 = 1 \implies \det(Q) = \pm 1$.

9. **[Geometry] What is the relationship between the residual $\mathbf{b} - A\hat{x}$ and $C(A)$?**

   ??? success "Solution"
       The residual vector is orthogonal to the column space of $A$ (it lies in the left nullspace).

10. **[Numerical] Why is QR decomposition more stable than the normal equations for least squares?**

   ??? success "Solution"
        The condition number of $A^T A$ is the square of the condition number of $A$, which amplifies errors; whereas $Q$ is an isometry, maintaining numerical stability.

## Chapter Summary

Orthogonality establishes the metric geometry of linear spaces:

1.  **Coordinate Purification**: The Gram-Schmidt process shows how to extract independent, standardized "pure" dimensions from a cluttered basis.
2.  **Optimal Error Handling**: Through orthogonal projection, the method of least squares finds the most reasonable "compromise" for inconsistent observation data caused by noise.
3.  **Stability Assurance**: Orthogonal matrices $Q$, due to their energy-preserving (norm-preserving) nature, form the underlying core of modern numerical linear algebra algorithms (e.g., Householder transforms, QR algorithm).
