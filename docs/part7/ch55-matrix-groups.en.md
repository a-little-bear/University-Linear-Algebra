# Chapter 55A: Matrix Groups and Classical Lie Groups

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Orthogonality (Ch07) · Basics of Group Theory

**Chapter Outline**: From Sets to Group Structures → The General Linear Group $GL(n, F)$ → The Special Linear Group $SL(n, F)$ → The Orthogonal Group $O(n)$ and Special Orthogonal Group $SO(n)$ → The Unitary Group $U(n)$ and Special Unitary Group $SU(n)$ → The Symplectic Group $Sp(2n)$ → Geometric Meaning of Matrix Groups as Lie Groups → Center, Connectivity, and Compactness → Applications: Symmetry in Physics (Rotations, Gauge Fields), and Motion Models in Computer Vision

**Extension**: Matrix groups are sets of transformations that preserve a specific structure; they combine algebraic operations with continuous manifolds, proving that symmetry can be quantified into specific subgroups. They form the universal algebraic foundation for theoretical physics, differential geometry, and modern image processing.

</div>

In linear algebra, we focus not only on individual matrices but on sets of matrices with common properties. If these sets form a group under multiplication, they are called **Matrix Groups**. **Classical Lie Groups** are the most important among these, describing rotations, reflections, and scalings in space. This chapter introduces the definitions, structures, and central role of these groups in describing the symmetries of nature.

---

## 55A.1 General and Special Linear Groups

!!! definition "Definition 55A.1 (General Linear Group $GL(n)$)"
    The group of all invertible $n \times n$ matrices over a field $F$:
    $$GL(n, F) = \{ A \in M_n(F) : \det(A) \neq 0 \}$$

!!! definition "Definition 55A.2 (Special Linear Group $SL(n)$)"
    The subgroup of matrices with determinant 1:
    $$SL(n, F) = \{ A \in M_n(F) : \det(A) = 1 \}$$
    **Geometric Insight**: $SL(n)$ contains all linear transformations that preserve volume and orientation.

---

## 55A.2 Orthogonal and Unitary Groups

!!! definition "Definition 55A.3 (Orthogonal Groups)"
    1.  **Orthogonal Group $O(n)$**: Matrices preserving the Euclidean inner product, $A^T A = I$.
    2.  **Special Orthogonal Group $SO(n)$**: Matrices in $O(n)$ with $\det(A) = 1$.
    **Physical Meaning**: $SO(n)$ describes **pure rotations** in $n$-dimensional space.

!!! definition "Definition 55A.4 (Unitary Group $U(n)$)"
    Matrices in complex space preserving the Hermitian inner product, $A^* A = I$.

---

## 55A.3 Geometric Features of Lie Groups

!!! technique "Intuition: Lie Groups"
    Matrix groups are not just algebraic objects but **continuous manifolds** in real or complex space.
    - **Compactness**: $O(n), SO(n), U(n), SU(n)$ are compact (bounded and closed), which ensures desirable representation properties.
    - **Connectivity**: $SO(n)$ is connected, while $O(n)$ consists of two disconnected components (rotations and reflections).

---

## Exercises

**1. [Basics] Determine if the set of all $n \times n$ singular matrices forms a group.**

??? success "Solution"
    **Determination:**
    1. A group must contain the identity matrix $I$.
    2. Singular matrices satisfy $\det(A) = 0$.
    3. Since $\det(I) = 1 \neq 0$, the identity is not in the set.
    **Conclusion**: No. Furthermore, the product of singular matrices is not always singular (though the product of non-singular ones is always non-singular).

**2. [Property] Prove: If $A, B \in SL(n)$, then $AB \in SL(n)$.**

??? success "Solution"
    **Proof:**
    1. By definition, $\det(A) = 1$ and $\det(B) = 1$.
    2. Using the multiplicative property: $\det(AB) = \det(A)\det(B)$.
    3. $\det(AB) = 1 \cdot 1 = 1$.
    **Conclusion**: Closure is satisfied. Since the identity and inverses also have determinant 1, $SL(n)$ is a group.

**3. [Calculation] Does the matrix $\begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$ in $O(2)$ belong to $SO(2)$?**

??? success "Solution"
    **Calculation:**
    $\det = \cos^2\theta - (-\sin^2\theta) = \cos^2\theta + \sin^2\theta = 1$.
    **Conclusion**: Yes. These matrices represent planar rotations and are the elements of $SO(2)$.

**4. [Counter-example] Give an example of a matrix in $O(2)$ that is not in $SO(2)$ and describe its meaning.**

??? success "Solution"
    **Example:** $A = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$.
    1. Orthogonality: $A^T A = I$. (Checks out).
    2. Determinant: $\det = -1$. (Not 1).
    **Geometric Insight**: This is a **reflection** across the $x$-axis. Reflection changes the orientation of the space.

**5. [Dimension] Find the dimension of $SO(3)$ as a manifold.**

??? success "Solution"
    **Conclusion: 3.**
    **Analysis**:
    1. A $3 \times 3$ matrix has 9 parameters.
    2. The constraint $A^T A = I$ provides 6 independent equations (3 for the diagonal, 3 for the upper triangle).
    3. Degrees of freedom $= 9 - 6 = 3$. These correspond to the three Euler angles.

**6. [Unitary] Prove that the eigenvalues of a unitary matrix $U \in U(n)$ have modulus 1.**

??? success "Solution"
    **Proof:**
    1. Let $U\mathbf{v} = \lambda \mathbf{v}$.
    2. Norm: $\|U\mathbf{v}\|^2 = (U\mathbf{v})^*(U\mathbf{v}) = \mathbf{v}^* U^* U \mathbf{v}$.
    3. Since $U \in U(n)$, $U^* U = I$, so $\|U\mathbf{v}\|^2 = \|\mathbf{v}\|^2$.
    4. Also, $\|U\mathbf{v}\|^2 = |\lambda|^2 \|\mathbf{v}\|^2$.
    5. Thus $|\lambda|^2 = 1 \implies |\lambda| = 1$.

**7. [Center] Find the center of $SL(n, \mathbb{C})$ (matrices that commute with all elements).**

??? success "Solution"
    **Conclusion:**
    The set of scalar matrices $\{ \omega I : \omega^n = 1 \}$, where $\omega$ is an $n$-th root of unity. A matrix must be scalar to commute with all others, and its determinant must be 1.

**8. [Application] Briefly state the relationship between $SU(2)$ and $SO(3)$.**

??? success "Solution"
    **Double Cover:**
    There is a 2:1 surjective homomorphism from $SU(2)$ to $SO(3)$. In physics, this explains why spinors (spin-1/2 particles) require a 720-degree rotation to return to their initial state—the mathematical core of rotation theory in quantum mechanics.

**9. [Comparison] To which group is $Sp(2, \mathbb{R})$ isomorphic?**

??? success "Solution"
    **Conclusion: $SL(2, \mathbb{R})$.**
    In 2D, a linear transformation preserving signed area (det=1) is identical to one preserving the symplectic form.

**10. [Physics] Why are physical laws often invariant under $SO(3)$?**

??? success "Solution"
    **Reasoning:**
    This represents the **Isotropy** of space. If physical laws are invariant under rotation, Noether's Theorem implies that the system satisfies **Conservation of Angular Momentum**. Matrix groups provide the precise micro-algebraic language for these macro-conservation laws.

## Chapter Summary

Matrix groups and Lie groups represent the high-level organization of linear algebra:

1.  **Quantification of Symmetry**: By defining groups like $SO(n)$ and $SU(n)$, we transform intuitive physical symmetries into rigorous algebraic constraint equations, establishing the standard paradigm for invariants.
2.  **Fusion of Continuous and Discrete**: Matrix groups possess both algebraic group properties and geometric manifold properties, serving as the junction between transformation theory and topology.
3.  **Skeleton of the Universe**: From atomic energy levels to galactic evolution, classical Lie groups provide the background for describing all conservation laws and transformation rules—one of the most universal structures in modern science.
