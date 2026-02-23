# Chapter 24: Matrix Manifolds

<div class="context-flow" markdown>

**Prerequisites**: Matrix Groups (Ch55) · Matrix Calculus (Ch47) · Unitary Matrices (Ch8) · Optimization (Ch25)

**Chapter Outline**: Differential Geometry Basics (Manifolds, Tangent Spaces) → Definition of Matrix Manifolds → Common Manifolds (Stiefel, Grassmann, Positive Definite Matrices) → Riemannian Metrics and Geodesics → Riemannian Gradients and Hessians → Retractions → Optimization Algorithms on Manifolds → Applications (Low-rank Matrix Completion, ICA, Computational Anatomy)

**Extension**: Matrix manifolds bring linear algebra into curved spaces; they provide a natural geometric framework for optimization problems with structural constraints (like orthogonality or fixed rank).

</div>

Matrix manifolds study sets of matrices with a geometric structure. Often, our search space is not the flat $\mathbb{R}^{n 	imes n}$, but a curved surface such as "all orthogonal matrices" or "all matrices of rank $k$." By optimizing directly on these manifolds, we can automatically satisfy complex matrix constraints.

---

## 24.1 Core Manifolds and Geometry

!!! definition "Definition 24.1 (Stiefel Manifold)"
    The Stiefel manifold $St(k, n)$ is the set of all $n 	imes k$ matrices satisfying $U^T U = I_k$. It is a generalization of the concept of orthogonal matrices.

!!! definition "Definition 24.3 (Grassmann Manifold)"
    The Grassmann manifold $Gr(k, n)$ is the set of all $k$-dimensional subspaces of $\mathbb{R}^n$. It differs from the Stiefel manifold because it treats bases that span the same space as the same point.

---

## Exercises

1. **[Tangent Space] Find the tangent space of the unit circle ($1 	imes 1$ orthogonal matrices) at point 1.**
   ??? success "Solution"
       The condition is $x^2 = 1$. Differentiating gives $2x \dot{x} = 0$. At $x=1$, $\dot{x}=0$. Thus the tangent space is the origin. Generalizing to $O(n)$, the tangent space at $I$ is the space of skew-symmetric matrices.

2. **[Dimension] Calculate the dimension of the $O(n)$ manifold of $n 	imes n$ orthogonal matrices.**
   ??? success "Solution"
       The dimension equals the dimension of its tangent space (skew-symmetric matrices): $n(n-1)/2$. This represents the degrees of freedom required to rotate an $n$-dimensional object.

3. **[Riemannian Gradient] Describe the relation between the Riemannian gradient $\operatorname{grad} f(X)$ and the Euclidean gradient $
abla f(X)$.**
   ??? success "Solution"
       The Riemannian gradient is the **orthogonal projection** of the Euclidean gradient onto the tangent space of the manifold. This ensures the search direction stays along the manifold surface and doesn't "fly off" the search space.

4. **[Geodesic] What is a geodesic? What role does it play in manifold optimization?**
   ??? success "Solution"
       A geodesic is the "shortest path" (generalized line) between two points on a manifold. In optimization, geodesics define the update direction for steps: $X_{k+1} = \operatorname{Exp}_{X_k}(\eta)$.

5. **[Retraction] Why are retractions often used instead of the exponential map in practical computations?**
   ??? success "Solution"
       The exponential map involves exact geodesic calculations, which are often computationally expensive (e.g., involving matrix exponentials). A retraction is a first-order approximation of the exponential map that maps tangent vectors back to the manifold at a lower cost while maintaining algorithm convergence.

6. **[PD Manifold] Describe the Riemannian metric on the space of positive definite matrices $S_n^{++}$.**
   ??? success "Solution"
       The common metric is $\langle \xi, \eta angle_X = \operatorname{tr}(X^{-1} \xi X^{-1} \eta)$. This metric makes the space of PD matrices a symmetric space with negative curvature, where the inversion operation becomes an isometry.

7. **[Stiefel vs Grassmann] If two matrices $U_1, U_2 \in St(k, n)$ satisfy $U_1 = U_2 Q$ ($Q$ is a rotation), do they represent the same point on the Grassmann manifold?**
   ??? success "Solution"
       Yes. The Grassmann manifold focuses on the subspace itself, not the choice of basis. Since the columns of $U_1$ and $U_2$ span the same space, they are equivalent in $Gr(k, n)$.

8. **[Application] How do matrix manifolds handle "low-rank" constraints?**
   ??? success "Solution"
       The set of fixed-rank matrices $\mathcal{M}_r = \{X \in \mathbb{R}^{m 	imes n} : \operatorname{rank}(X) = r\}$ is a smooth manifold. By defining gradient descent on this manifold, one can achieve exact low-rank matrix completion without explicitly using nuclear norm penalties.

9. **[Calculation] Write the equation satisfied by an element $\xi$ of the tangent space at $Q$ for the $O(n)$ manifold.**
   ??? success "Solution"
       From $Q^T Q = I$, differentiating yields $\dot{Q}^T Q + Q^T \dot{Q} = 0$. Letting $\xi = \dot{Q}$, then $\xi^T Q + Q^T \xi = 0$.

10. **[Riemannian Hessian] Why is the Riemannian Hessian more complex than the Euclidean Hessian?**
    ??? success "Solution"
        Because it must account for the curvature of the manifold. Its calculation involves the Levi-Civita connection (covariant derivative), which captures the rotation of the tangent plane as it moves.

## Chapter Summary

Matrix manifolds are the geometric destination for constrained optimization:

1. **Internalizing Structure**: Transforming matrix constraints from "external boundaries" into "internal properties" of the space itself.
2. **Curvature Wisdom**: Utilizing non-Euclidean metrics (like the affine-invariant metric for PD matrices) can significantly improve numerical stability.
3. **Algorithmic Unity**: Manifold optimization provides a unified calculus framework for solving problems with orthogonal, low-rank, or symplectic constraints.
