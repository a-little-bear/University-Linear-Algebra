# Chapter 55: Matrix Groups

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Determinants (Ch03) · Orthogonality and Unitarity (Ch07-08) · Symplectic Matrices (Ch53)

**Chapter Outline**: Group Theory Basics and the Definition of Matrix Groups → General Linear Group $GL(n)$ and Special Linear Group $SL(n)$ → Orthogonal Group $O(n)$ and Rotation Group $SO(n)$ → Unitary Group $U(n)$ and Special Unitary Group $SU(n)$ → Symplectic Group $Sp(2n)$ → Topological Properties (Connectedness, Compactness) → The Exponential Mapping and Introduction to Lie Algebras → Applications: Symmetries in Physics, Robot Attitude, and Homography in Computer Vision

**Extension**: Matrix groups are the intersection of linear algebra, topology, and group theory; they are not just mathematical carriers of symmetry but the backbone of modern physics, from crystal structures to gauge theory.

</div>

In previous chapters, we focused on the properties of individual matrices. Now, we shift our perspective to sets of matrices that share common properties—**Matrix Groups**. By studying the algebraic structure and geometric form of these groups, we can more deeply understand the symmetries of space. Each matrix group corresponds to the preservation of a specific geometric property (such as length, volume, or area).

---

## 55.1 Basic Group Classes

!!! definition "Definition 55.1 (General and Special Linear Groups)"
    1.  **General Linear Group $GL(n, F)$**: The group of all $n \times n$ non-singular matrices over field $F$.
    2.  **Special Linear Group $SL(n, F)$**: The subgroup of $GL(n, F)$ consisting of matrices with determinant 1. These preserve the **oriented volume** of space.

---

## 55.2 Classical Groups and Isometries

!!! definition "Definition 55.2 (Orthogonal and Rotation Groups)"
    1.  **Orthogonal Group $O(n)$**: The group of matrices satisfying $Q^T Q = I$. These preserve Euclidean distances and angles.
    2.  **Special Orthogonal Group $SO(n)$**: The subgroup of $O(n)$ with determinant 1, representing pure rotations.

!!! definition "Definition 55.3 (Unitary Groups $U(n)$ and $SU(n)$)"
    The **Unitary Group** consists of complex matrices satisfying $U^* U = I$. The **Special Unitary Group** $SU(n)$ further requires $\det(U) = 1$. These are the foundations of evolution in quantum mechanics.

---

## 55.3 Topological Properties and Lie Algebras

!!! technique "Compactness and Connectedness"
    - **Compactness**: $O(n), SO(n), U(n),$ and $SU(n)$ are compact (their entries are bounded and the sets are closed). $GL(n)$ and $SL(n)$ are non-compact.
    - **Connectedness**: $SO(n)$ is connected, while $O(n)$ consists of two disconnected components (rotations and reflections).

!!! theorem "Theorem 55.1 (The Exponential Map)"
    A matrix group is linked to its corresponding **Lie Algebra** (tangent space at the identity) via the exponential map $e^A$. For instance, the Lie algebra of $SO(n)$ consists of all skew-symmetric matrices.

---

## Exercises

1.  **[Basics] Prove: If $A, B \in SL(n)$, then $AB \in SL(n)$.**
    ??? success "Solution"
        $\det(AB) = \det(A)\det(B) = 1 \cdot 1 = 1$. Thus it satisfies the definition.

2.  **[Determinant] Prove that the determinant of any element in $O(n)$ is $\pm 1$.**
    ??? success "Solution"
        $Q^T Q = I \implies \det(Q)^2 = 1 \implies \det(Q) = \pm 1$.

3.  **[Inverse] Prove that if $A$ is an orthogonal matrix, its inverse $A^{-1}$ is also orthogonal.**
    ??? success "Solution"
        $A^{-1} = A^T$. Since $(A^T)^T A^T = A A^T = I$, $A^T$ is orthogonal.

4.  **[Unitary] How many real parameters are needed to describe $SU(2)$?**
    ??? success "Solution"
        3 parameters. It is isomorphic to the unit 3-sphere $S^3$, often represented using 3 Euler angles or Pauli matrices.

5.  **[Symplectic] What is the relationship between $Sp(2, \mathbb{R})$ and $SL(2, \mathbb{R})$?**
    ??? success "Solution"
        In 2 dimensions, the symplectic condition is equivalent to having a determinant of 1, so the groups are identical.

6.  **[Physics] Why is $U(n)$ critical in quantum mechanics?**
    ??? success "Solution"
        Because it preserves the norm (probability) of complex state vectors.

7.  **[Robotics] Which two parts comprise the Special Euclidean group $SE(3)$ for rigid body transforms?**
    ??? success "Solution"
        A rotation matrix from $SO(3)$ and a translation vector in $\mathbb{R}^3$.

8.  **[Topology] Is $GL(n, \mathbb{R})$ connected?**
    ??? success "Solution"
        No. it is partitioned into two components by the sign of the determinant ($\det A > 0$ and $\det A < 0$).

9.  **[Lie Algebra] Prove: If $X$ is a skew-symmetric matrix, then $e^X$ is an orthogonal matrix.**
    ??? success "Solution"
        $(e^X)^T = e^{X^T} = e^{-X} = (e^X)^{-1}$. Thus $(e^X)^T e^X = I$.

10. **[Application] What is a Haar measure?**

   ??? success "Solution"
        A shift-invariant integral measure defined on a compact matrix group, allowing for averaging over the group (e.g., averaging over all possible rotation orientations).

## Chapter Summary

Matrix groups establish the global geometric structure of linear transformations:

1.  **Systematization of Symmetry**: Matrix groups integrate isolated properties (length/volume preservation) into rigorous structures, providing standard algebraic language for physical laws.
2.  **Intertwining Algebra and Topology**: Studying compactness and connectedness reveals the large-scale geometry of operator space, essential for understanding gimbal lock, singularities, and stability.
3.  **Bridging Local and Global**: The exponential map demonstrates how infinitesimal transformations (Lie algebras) generate finite transformations (matrix groups), establishing the computational paradigm for continuous evolution.
