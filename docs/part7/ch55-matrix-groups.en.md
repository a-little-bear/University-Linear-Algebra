# Chapter 55: Matrix Groups and Lie Algebras

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Groups and Symmetry · Matrix Exponential (Ch13) · Differential Calculus

**Chapter Outline**: Matrix Lie Groups ($GL, SL, O, SO, U, SU, Sp$) → The Exponential Map $\exp: \mathfrak{g} \to G$ → Lie Algebras as Tangent Spaces → Commutators and the Lie Bracket → Adjoint Representation → Baker-Campbell-Hausdorff Formula → Relation between Group Topologies and Algebraic Structures

**Extension**: Matrix groups are the fundamental language of symmetry in physics (Standard Model) and geometry (Riemannian manifolds).

</div>

Matrix Lie groups are groups of matrices that also possess the structure of a smooth manifold. Every such group has an associated **Lie algebra**, which captures the infinitesimal structure of the group near the identity. The link between the group and its algebra is the **matrix exponential**, which maps "velocities" in the algebra to "rotations" or "displacements" in the group.

---

## 55.1 Core Groups and Algebras

!!! definition "Definition 55.1 (Matrix Lie Group)"
    A subgroup $G \subseteq GL(n, \mathbb{C})$ is a matrix Lie group if it is a closed subset of $GL(n, \mathbb{C})$.

!!! theorem "Theorem 55.1 (The Lie Algebra)"
    The Lie algebra $\mathfrak{g}$ of a group $G$ consists of all matrices $X$ such that $e^{tX} \in G$ for all $t \in \mathbb{R}$. The algebra is equipped with the **Lie bracket** $[X, Y] = XY - YX$.

---

## Exercises

1. **[Special Linear] Determine the Lie algebra $\mathfrak{sl}(n, \mathbb{C})$ of the group $SL(n, \mathbb{C})$.**
   ??? success "Solution"
       $A \in SL(n, \mathbb{C}) \implies \det A = 1$. Using $\det(e^{tX}) = e^{t \operatorname{tr}(X)}$, the condition $\det(e^{tX}) = 1$ for all $t$ implies $\operatorname{tr}(X) = 0$. Thus $\mathfrak{sl}(n, \mathbb{C})$ consists of all trace-zero matrices.

2. **[Orthogonal] Show that the Lie algebra $\mathfrak{so}(n, \mathbb{R})$ consists of skew-symmetric matrices.**
   ??? success "Solution"
       $R \in SO(n) \implies R^T R = I$. Let $R(t) = e^{tX}$. Then $(e^{tX})^T e^{tX} = e^{tX^T} e^{tX} = I$. Differentiating at $t=0$ gives $X^T + X = 0$, so $X^T = -X$.

3. **[Unitary] Characterize the Lie algebra $\mathfrak{su}(n)$ of the special unitary group $SU(n)$.**
   ??? success "Solution"
       It consists of traceless anti-Hermitian matrices ($X^* = -X$ and $\operatorname{tr}(X) = 0$). This algebra has dimension $n^2 - 1$.

4. **[BCH Formula] State the first few terms of the Baker-Campbell-Hausdorff formula for $e^X e^Y = e^Z$.**
   ??? success "Solution"
       $Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}([X, [X, Y]] + [Y, [Y, X]]) + \dots$. This shows that the group multiplication is locally determined by the Lie bracket.

5. **[Adjoint] Define the adjoint representation $\operatorname{Ad}_g(X)$.**
   ??? success "Solution"
       $\operatorname{Ad}_g(X) = g X g^{-1}$. This operator maps an element of the Lie algebra to another element of the same algebra via the group action, preserving the Lie bracket.

6. **[Commutator] Prove the Jacobi Identity for the Lie bracket: $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$.**
   ??? success "Solution"
       Expand using $[X, Y] = XY - YX$ and verify that all 12 terms cancel out in pairs. This identity is the fundamental axiom of Lie algebras.

7. **[Center] What is the center of the Lie algebra $\mathfrak{gl}(n, \mathbb{C})$?**
   ??? success "Solution"
       The center consists of matrices that commute with all other matrices in the algebra. For $\mathfrak{gl}(n)$, these are exactly the scalar matrices $cI$.

8. **[Spin Groups] Relate $SU(2)$ to $SO(3)$ via a homomorphism.**
   ??? success "Solution"
       There is a $2$-to-$1$ surjective homomorphism from $SU(2)$ to $SO(3)$. This implies that $SU(2)$ is the universal cover of $SO(3)$, a fact crucial for describing electron spin in quantum mechanics.

9. **[Exponential Map] Is the exponential map always surjective for $GL(n, \mathbb{C})$?**
   ??? success "Solution"
       Yes, for $GL(n, \mathbb{C})$, every invertible matrix has a matrix logarithm. However, for $SL(2, \mathbb{R})$, it is not surjective (e.g., matrices with negative eigenvalues and non-trivial Jordan blocks lack real logarithms).

10. **[Physics] Why are Lie groups central to the Standard Model of particle physics?**
    ??? success "Solution"
        The symmetries of the fundamental interactions are described by the gauge groups $U(1)$ (electromagnetism), $SU(2)$ (weak force), and $SU(3)$ (strong force). Particles correspond to irreducible representations of these matrix groups.

## Chapter Summary

This chapter explores the synthesis of algebra, topology, and geometry in matrix groups:

1. **Infinitesimal Analysis**: Used the matrix exponential to linearize groups into Lie algebras.
2. **Standard Symmetries**: Characterized the algebras of orthogonal, unitary, and special linear groups through trace and symmetry constraints.
3. **Algebraic Structure**: Defined the Lie bracket as the fundamental operation capturing non-commutativity.
4. **Global-Local Duality**: Demonstrated through the BCH formula and Lie's theorems how the local algebra structure dictates the global group behavior.
