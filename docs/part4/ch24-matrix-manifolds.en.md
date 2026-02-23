# Chapter 24: Matrix Manifolds

<div class="context-flow" markdown>

**Prerequisites**: Matrix Groups (Ch55A) · Matrix Calculus (Ch47A) · Basics of Differential Geometry

**Chapter Outline**: From Euclidean Space to Surface Constraints → Definition of Matrix Manifolds → Typical Manifolds: Stiefel Manifolds (Orthonormal Bases), Grassmann Manifolds (Subspaces), and Positive Definite Manifolds → Tangent and Normal Spaces → Riemannian Metrics → Exponential and Logarithmic Maps (Geodesics) → Optimization on Manifolds: Riemannian Gradient Descent (RGD) → Applications: Low-rank Matrix Recovery, Pose Estimation in Computer Vision, and Geometric Generalization of PCA

**Extension**: Matrix manifolds provide the ultimate geometric framework for studying "constrained" matrix evolution; they transform complex matrix constraints (such as $Q^T Q = I$) into intrinsic properties of a manifold, enabling "smooth" gradient optimization in curved operator spaces. This is at the heart of modern non-linear signal processing.

</div>

In many practical problems, we seek optimal solutions within sets of matrices that satisfy specific constraints, such as orthogonality. **Matrix Manifolds** treat these sets as smooth surfaces in high-dimensional space. By introducing tools from **Riemannian Geometry**, we can perform gradient descent or Newton iterations directly on these surfaces without worrying about violating the constraints. This chapter introduces this geometric perspective that underpins modern high-level optimization algorithms.

---

## 24.1 Typical Matrix Manifolds

!!! definition "Definition 24.1 (Stiefel Manifold $St(k, n)$)"
    The set of all $n \times k$ matrices with orthonormal columns:
    $$St(k, n) = \{ X \in \mathbb{R}^{n \times k} : X^T X = I_k \}$$
    When $k=n$, this reduces to the Orthogonal Group $O(n)$.

!!! definition "Definition 24.2 (Positive Definite Manifold $\mathcal{P}_n$)"
    The manifold of all $n \times n$ symmetric positive definite matrices. Its natural Riemannian metric leads to the matrix geometric mean (see Ch46B) as its geodesic midpoint.

---

## 24.2 Tangent Spaces and Projections

!!! technique "Technique: Gradients on Manifolds"
    Optimization on a manifold involves three key steps:
    1.  **Compute Euclidean Gradient** $\nabla f$.
    2.  **Orthogonal Projection**: Project $\nabla f$ onto the **Tangent Space** $\mathcal{T}_X \mathcal{M}$ to obtain the **Riemannian Gradient**.
    3.  **Retraction**: Move along the tangent direction and map the point back onto the manifold to ensure the constraint is maintained.

---

## 24.3 Geodesics and Distance

!!! theorem "Theorem 24.1 (Geodesic Equation)"
    The shortest path between two points on a manifold is called a **Geodesic**. On the manifold of positive definite matrices, the geodesic has an explicit formula:
    $$\gamma(t) = A^{1/2} (A^{-1/2} B A^{-1/2})^t A^{1/2}$$
    At $t=0.5$, this is exactly the matrix geometric mean.

---

## Exercises

**1. [Basics] Describe the geometry of the Stiefel manifold $St(1, 3)$.**

??? success "Solution"
    **Analysis:**
    1. Here $n=3, k=1$.
    2. Matrix $X$ is a $3 \times 1$ vector $\mathbf{v}$.
    3. The constraint $X^T X = I_1$ implies $\mathbf{v}^T \mathbf{v} = 1$.
    **Conclusion**: This is the unit sphere $S^2$ in $\mathbb{R}^3$.

**2. [Dimension] Find the dimension of the Orthogonal Group $O(n)$ as a matrix manifold.**

??? success "Solution"
    **Calculation:**
    1. An $n \times n$ matrix has $n^2$ degrees of freedom.
    2. The constraint $Q^T Q = I$ is symmetric, providing $n(n+1)/2$ independent equations.
    3. Dimension $d = n^2 - n(n+1)/2 = n(n-1)/2$.
    **Note**: This matches the dimension of the space of skew-symmetric matrices (Lie algebra $\mathfrak{so}(n)$).

**3. [Tangent Space] Find the condition for a vector to be tangent to the unit circle $S^1$ at $(1, 0)$.**

??? success "Solution"
    **Derivation:**
    1. The constraint is $x^2 + y^2 = 1$.
    2. Differentiating with respect to time: $2x \dot{x} + 2y \dot{y} = 0$.
    3. At $(1, 0)$, we get $2(1)\dot{x} + 2(0)\dot{y} = 0 \implies \dot{x} = 0$.
    **Conclusion**: Tangent vectors must be of the form $(0, v_y)^T$, i.e., perpendicular to the radial direction.

**4. [Property] Prove: The tangent space of $O(n)$ at $Q$ consists of matrices $Q\Omega$ where $\Omega$ is skew-symmetric.**

??? success "Solution"
    **Proof:**
    1. Let $\gamma(t) = Q(t)$ be a curve on the manifold with $Q(0) = Q$.
    2. $Q(t)^T Q(t) = I \implies \dot{Q}^T Q + Q^T \dot{Q} = 0$.
    3. Let $\Omega = Q^T \dot{Q}$. Then $\Omega^T + \Omega = (Q^T \dot{Q})^T + Q^T \dot{Q} = \dot{Q}^T Q + Q^T \dot{Q} = 0$.
    **Conclusion**: $\Omega$ is skew-symmetric and $\dot{Q} = Q\Omega$.

**5. [Retraction] Why is a "Retraction" mapping needed in manifold optimization?**

??? success "Solution"
    **Reasoning:**
    The tangent space is a flat linear approximation. When we take a step $X + \eta$ in the tangent direction, the curvature of the manifold causes the new point to move "off" the surface. A retraction maps this point back onto the manifold (e.g., via the $Q$ factor of a QR decomposition), ensuring every iteration is a valid manifold element.

**6. [Grassmann] Briefly explain the relationship between Grassmann manifolds $Gr(k, n)$ and Stiefel manifolds.**

??? success "Solution"
    **Connection:**
    The Grassmann manifold is the **quotient space** of the Stiefel manifold: $Gr(k, n) = St(k, n) / O(k)$. In a Grassmannian, we do not care about the specific choice of basis, only the subspace they span. Two matrices in $St(k, n)$ whose columns span the same subspace are treated as the same point in $Gr(k, n)$.

**7. [Riemannian Gradient] Given an Euclidean gradient $G$, find the Riemannian gradient on the orthogonal group.**

??? success "Solution"
    **Formula:**
    $\operatorname{grad} f(X) = G - X G^T X$ (for numerator layout) or a similar symmetrized projection. The core idea is to strip away the components of the gradient that are perpendicular to the tangent plane.

**8. [Application] How are matrix manifolds used in Matrix Completion?**

??? success "Solution"
    The set of matrices with a fixed rank $r$ forms a manifold. By optimizing directly on the "rank-$r$ manifold," the rank constraint is strictly maintained throughout the process, which is often more efficient than nuclear-norm regularizations in the full space.

**9. [Distance] How does the Riemannian distance on $\mathcal{P}_n$ differ from Euclidean distance?**

??? success "Solution"
    **Comparison:**
    - **Euclidean**: $\|A - B\|_F$. It can lead to negative eigenvalues if the path crosses the boundary.
    - **Riemannian**: $\|\ln(A^{-1/2} B A^{-1/2})\|_F$. It accounts for the curvature, placing the boundary (singular matrices) at infinity. This distance is invariant under scaling transformations.

**10. [Application] Briefly state the role of manifold optimization in Camera Pose Estimation.**

??? success "Solution"
    The rotation of a camera belongs to the $SO(3)$ manifold. Using manifold algorithms allows us to optimize camera parameters directly on $SO(3)$, avoiding Gimbal Lock (from Euler angles) and the need for constant re-normalization (from quaternions), leading to the most efficient geometric registration.

## Chapter Summary

Matrix manifolds are the bridge from linear algebra to geometric optimization:

1.  **Internalizing Constraints**: They transform tedious non-linear equality constraints into intrinsic geometric features, establishing the "freedom within constraints" algorithmic logic.
2.  **Measuring Curvature**: Via Riemannian metrics and geodesics, matrix space evolves from a flat grid to a deep geometric entity, providing the only tool for defining natural distances between operators.
3.  **Algorithmic Elegance**: Manifold optimization proves that complex matrix constraint problems can be solved via projections and retractions, serving as the mathematical engine for modern large-scale industrial vision, robotics, and signal recovery.
