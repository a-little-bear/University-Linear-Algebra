# Chapter 24: Matrix Manifolds and Geometric Optimization

<div class="context-flow" markdown>

**Prerequisites**: Matrix Calculus (Ch47A) · Orthogonality (Ch7) · Matrix Groups (Ch55) · Optimization (Ch25)

**Chapter Outline**: Space of Matrices as a Manifold → Tangent and Normal Spaces → The Stiefel Manifold $St(k, n)$ → The Grassmann Manifold $Gr(k, n)$ → Orthogonal and Unitary Groups → Positive Definite Cone as a Manifold → Riemannian Metrics on Matrix Manifolds → Geodesics and Parallel Transport → Optimization on Manifolds (Retractions)

**Extension**: Matrix manifolds provide the geometric framework for constrained optimization, allowing for gradients and Newton steps to be taken while remaining on a curved surface like the set of orthogonal matrices.

</div>

Matrix manifolds treat subsets of matrices as smooth geometric surfaces. Instead of viewing an orthogonal matrix as a matrix with constraints, we view the **Orthogonal Group** $O(n)$ as a curved space (manifold) where each point is a matrix. This geometric perspective is essential for solving optimization problems where the solution must maintain a specific structure, such as being a projection or a rotation. This chapter formalizes the concepts of **Tangent Spaces** and **Retractions**, enabling the extension of standard optimization algorithms to curved matrix spaces.

---

## 24.1 Definitions and Core Manifolds

!!! definition "Definition 24.1 (Tangent Space)"
    The tangent space $T_A \mathcal{M}$ to a manifold $\mathcal{Cl}$ at point $A$ is the set of all possible velocities of curves passing through $A$. For $O(n)$, the tangent space at $I$ is the set of skew-symmetric matrices.

!!! theorem "Theorem 24.1 (The Stiefel Manifold)"
    The Stiefel manifold $St(k, n)$ is the set of all $n \times k$ matrices with orthonormal columns. Its dimension is $nk - \frac{1}{2}k(k+1)$.

---

## Exercises

1. **[Fundamentals] Find the tangent space to the Orthogonal Group $O(n)$ at $Q$.**
   ??? success "Solution"
       Let $Q(t)$ be a curve with $Q(0)=Q$ and $Q(t)^T Q(t) = I$. Differentiating at $t=0$ gives $\dot{Q}^T Q + Q^T \dot{Q} = 0$. Thus $\dot{Q} = Q \Omega$ where $\Omega$ is skew-symmetric. $T_Q O(n) = \{ Q\Omega : \Omega^T = -\Omega \}$.

2. **[Grassmannian] Define the Grassmann manifold $Gr(k, n)$ and its difference from the Stiefel manifold.**
   ??? success "Solution"
       The Stiefel manifold cares about the specific orthonormal basis (the columns). The Grassmannian only cares about the **subspace** spanned by those columns. $Gr(k, n)$ is the set of all $k$-dimensional subspaces of $\mathbb{R}^n$.

3. **[PSD Cone] Describe the geometry of the cone of positive definite matrices $S_n^{++}$.**
   ??? success "Solution"
       It is an open, convex cone. In Riemannian geometry, it is a symmetric space with a metric given by $ds^2 = \operatorname{tr}(A^{-1} dA A^{-1} dA)$, which ensures the boundary (singular matrices) is at infinite distance.

4. **[Retraction] What is a "Retraction" in manifold optimization?**
   ??? success "Solution"
       A mapping from the tangent space back to the manifold. For example, if we take a step in the tangent space of $O(n)$, the result is not orthogonal. The retraction (like the polar decomposition or QR) "pulls" the point back onto the manifold.

5. **[Projection] Relate the Grassmannian to projection matrices.**
   ??? success "Solution"
       $Gr(k, n)$ is isomorphic to the set of rank-$k$ symmetric projection matrices $P^2=P$. This provides a concrete matrix representation for abstract subspaces.

6. **[Geodesics] Define a geodesic on a matrix manifold.**
   ??? success "Solution"
       The "straightest possible" path between two matrices on the manifold. On $O(n)$, geodesics are of the form $Q(t) = Q e^{t\Omega}$ for skew-symmetric $\Omega$.

7. **[Dimension] Calculate the dimension of the Grassmannian $Gr(k, n)$.**
   ??? success "Solution"
       $\dim Gr(k, n) = k(n-k)$. This is the dimension of the Stiefel manifold minus the dimension of $O(k)$ (the group of basis changes).

8. **[Normal Space] What is the normal space to $O(n)$ at $I$?**
   ??? success "Solution"
       The set of all symmetric matrices. In the trace inner product, symmetric matrices are orthogonal to skew-symmetric matrices.

9. **[Eigenvalues] How do the eigenvalues of $A$ and $B$ relate to the distance between them on the PSD manifold?**
   ??? success "Solution"
       The Riemannian distance is $d(A, B) = \|\log(A^{-1/2} B A^{-1/2})\|_F$, which is the RMS of the logs of the generalized eigenvalues.

10. **[Applications] Where is optimization on the Stiefel manifold used?**
    ??? success "Solution"
        In electronic structure calculations (finding orthonormal orbitals that minimize energy) and in independent component analysis (ICA).

## Chapter Summary

This chapter establishes the geometric theory of constrained operators:

1. **Topological Mapping**: Defined the Stiefel and Grassmann manifolds as the geometric homes for orthonormal bases and subspaces.
2. **Tangent Calculus**: Developed the linearized local structure of matrix groups, enabling the use of derivatives on curved surfaces.
3. **Riemannian Metrics**: Formulated the concepts of distance and geodesics for comparing matrices based on their intrinsic geometry.
4. **Projective Optimization**: Integrated retractions and projections to enable robust optimization algorithms on structured matrix sets.
