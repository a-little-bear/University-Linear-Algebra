# Chapter 24: Matrix Manifolds

<div class="context-flow" markdown>

**Prerequisites**: Matrix Decompositions (Ch10) · Orthogonality (Ch07) · Matrix Analysis (Ch14) · Differential Geometry Basics

**Chapter Outline**: Introduction to Matrix Manifolds → The Stiefel Manifold → The Grassmann Manifold → Orthogonal Groups $O(n)$ and $SO(n)$ → Positive Definite Manifold $\mathcal{P}_n$ → Riemannian Geometry on Matrices (Tangent Spaces, Metrics) → Projections and Retractions → Optimization on Manifolds (Riemannian Gradient Descent) → Geodesics & Exponential Mapping → Applications: Computer Vision, Signal Processing, and Low-rank Matrix Recovery

</div>

In most of linear algebra, we work in flat vector spaces. However, many important sets of matrices—such as orthogonal matrices ($Q^T Q = I$) or fixed-rank matrices—do not form vector spaces because they are not closed under addition. Instead, they form curved surfaces known as **Matrix Manifolds**. Optimization on these manifolds allows us to solve constrained problems (like finding an optimal rotation or a best-fit subspace) using the tools of calculus and Riemannian geometry.

---

## 24.1 Important Matrix Manifolds

<div class="context-flow" markdown>

**Motivation**: How do we mathematically describe the "space of all rotations" or "the space of all $k$-dimensional subspaces"?

</div>

!!! definition "Definition 24.1 (The Stiefel Manifold)"
    The **Stiefel Manifold** $St(n, k)$ is the set of all $n \times k$ matrices with orthonormal columns:
    $$St(n, k) = \{X \in \mathbb{R}^{n \times k} : X^T X = I_k\}$$
    - If $k = n$, this is the **Orthogonal Group** $O(n)$.
    - If $k = 1$, this is the **Unit Sphere** $S^{n-1}$.

!!! definition "Definition 24.2 (The Grassmann Manifold)"
    The **Grassmann Manifold** $Gr(n, k)$ is the set of all $k$-dimensional subspaces of $\mathbb{R}^n$. Each point in the Grassmannian is a subspace, not a single matrix. We often represent a point by a matrix $X \in St(n, k)$ whose columns span the subspace.

!!! definition "Definition 24.3 (The Positive Definite Manifold)"
    The set of all $n \times n$ symmetric positive definite matrices $\mathcal{P}_n$ forms a manifold. Unlike the Stiefel manifold, $\mathcal{P}_n$ is an open set in the space of symmetric matrices, but its natural distance is curved (Riemannian).

---

## 24.2 Riemannian Geometry on Matrices

!!! definition "Definition 24.4 (Tangent Space)"
    The **tangent space** $T_X \mathcal{M}$ at a point $X$ on a manifold $\mathcal{M}$ is the set of all possible "velocity vectors" of paths passing through $X$.
    - For $St(n, k)$, the tangent space consists of matrices $Z$ such that $Z^T X + X^T Z = 0$.

!!! definition "Definition 24.5 (Riemannian Metric)"
    A **metric** $\langle \cdot, \cdot \rangle_X$ is an inner product on the tangent space. The most common metric is the Frobenius inner product: $\langle \eta, \zeta \rangle = \operatorname{tr}(\eta^T \zeta)$.

---

## 24.3 Optimization on Manifolds

<div class="context-flow" markdown>

**Problem**: How to update $X$ to minimize $f(X)$ while staying on the manifold? A standard step $X - \alpha \nabla f$ will fall off the surface.

</div>

!!! technique "Retraction: Staying on the Surface"
    A **retraction** $R_X(\xi)$ is a mapping from the tangent space back to the manifold. For the Stiefel manifold, common retractions include:
    1.  **QR Projection**: $R_X(\xi) = \operatorname{qf}(X + \xi)$ (the $Q$ part of the QR decomposition).
    2.  **Polar Projection**: $R_X(\xi) = (X + \xi)((X + \xi)^T (X + \xi))^{-1/2}$.

!!! algorithm "Algorithm 24.1 (Riemannian Gradient Descent)"
    1.  Compute the Euclidean gradient $\nabla f(X)$.
    2.  Project $\nabla f(X)$ onto the tangent space to get the **Riemannian gradient** $\operatorname{grad} f(X)$.
    3.  Take a step in the tangent space: $\xi = -\alpha \operatorname{grad} f(X)$.
    4.  Apply retraction: $X_{\text{new}} = R_X(\xi)$.

---

## 24.4 Geodesics and Distances

!!! theorem "Theorem 24.1 (Geodesics on $SO(n)$)"
    The shortest path (geodesic) between two rotation matrices $R_1$ and $R_2$ is a constant velocity rotation in the plane formed by their difference. This is used in computer graphics for smooth camera transitions (SLERP).

---

## Exercises

1.  **[Stiefel] What is the dimension of the Stiefel manifold $St(n, k)$?**
    ??? success "Solution"
        The dimension is $nk - \frac{1}{2}k(k+1)$. This accounts for the $nk$ total entries minus the $\frac{1}{2}k(k+1)$ independent constraints imposed by $X^T X = I_k$.

2.  **[Tangent] Show that for the orthogonal group $O(n)$, the tangent space at the identity $I$ is the set of skew-symmetric matrices.**
    ??? success "Solution"
        Let $Q(t)$ be a path in $O(n)$ with $Q(0)=I$. $Q(t)^T Q(t) = I$.
        Differentiating at $t=0$: $\dot{Q}(0)^T Q(0) + Q(0)^T \dot{Q}(0) = 0 \implies \dot{Q}^T + \dot{Q} = 0$.
        Thus, the tangent vectors $\dot{Q}$ must satisfy $A^T = -A$.

3.  **[Grassmann] Why is $Gr(n, k)$ considered a "quotient space"?**
    ??? success "Solution"
        Because many different orthonormal matrices in $St(n, k)$ span the same subspace. Specifically, $Gr(n, k) \cong St(n, k) / O(k)$, where we treat two matrices as equivalent if they differ only by a rotation within the subspace.

4.  **[Projection] How do you project a general matrix $Z$ onto the tangent space of $St(n, k)$ at $X$?**
    ??? success "Solution"
        $P_X(Z) = Z - X \operatorname{sym}(X^T Z)$, where $\operatorname{sym}(M) = (M + M^T)/2$.

5.  **[Metric] What is the natural distance between two positive definite matrices $A$ and $B$?**
    ??? success "Solution"
        It is the Riemannian distance: $d(A, B) = \|\log(A^{-1/2} B A^{-1/2})\|_F$. This distance is invariant under inversion and congruent transformations.

6.  **[Applications] Give an example of manifold optimization in computer vision.**
    ??? success "Solution"
        **Camera Pose Estimation**: Finding the rotation $R \in SO(3)$ and translation $t$ that minimizes the projection error of 3D points onto a 2D image.

7.  **[Exponential] Define the exponential map on a manifold.**
    ??? success "Solution"
        The exponential map $\operatorname{Exp}_X(\xi)$ takes a tangent vector $\xi$ and follows the geodesic (the "straightest" possible path on the curved surface) for unit time to reach a new point on the manifold.

8.  **[Low-rank] Is the set of matrices with rank *exactly* $k$ a manifold?**
    ??? success "Solution"
        Yes, it is a smooth manifold of dimension $(m+n-k)k$. However, the set of matrices with rank *at most* $k$ is not a manifold because it has singularities at points where the rank drops.

9.  **[Hessian] What is the Riemannian Hessian?**
    ??? success "Solution"
        It is the second-order derivative of a function on a manifold, accounting for the curvature of the surface. It is essential for Riemannian Newton methods.

****

??? success "Solution"
    

## Chapter Summary

Matrix manifolds transform constrained linear algebra into unconstrained geometric optimization:

1.  **Non-linear Structure**: Established that important sets of matrices (rotations, subspaces) are not vector spaces but curved manifolds.
2.  **Geometric Calculus**: Developed the concepts of tangent spaces and Riemannian gradients to allow calculus on curved surfaces.
3.  **Optimization Framework**: Introduced retractions and Riemannian gradient descent as the standard tools for solving optimization problems with matrix constraints.
4.  **Physical Relevance**: Demonstrated how manifolds provide the natural language for rigid body dynamics, camera geometry, and stable neural network architectures.
