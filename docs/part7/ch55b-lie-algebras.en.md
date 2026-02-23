# Chapter 55B: Lie Algebras and Infinitesimal Generators

<div class="context-flow" markdown>

**Prerequisites**: Matrix Groups (Ch55) · Matrix Exponential (Ch13) · Vector Spaces (Ch4) · Differential Geometry

**Chapter Outline**: Abstract Lie Algebra Definition → Commutator as Bracket → Structure Constants → Representations of Lie Algebras → Adjoint Representation $\operatorname{ad}$ → Solvable and Nilpotent Algebras → Engel's and Lie's Theorems → Killing Form → Cartan's Criterion → Semisimple Lie Algebras and Root Systems

**Extension**: Lie algebras are the "linearized" versions of Lie groups; they allow the study of continuous symmetries using purely linear algebraic tools (eigenvalues, roots, weights).

</div>

While Lie groups are manifolds, **Lie algebras** are vector spaces. By linearizing the group action at the identity, we transform non-linear geometric problems into linear algebraic ones. The bracket operation $[X, Y]$ replaces group multiplication, capturing the first-order non-commutativity of the symmetry group.

---

## 55B.1 Axioms and Basic Structures

!!! definition "Definition 55B.1 (Lie Algebra)"
    A Lie algebra $\mathfrak{g}$ over a field $K$ is a vector space equipped with a bilinear map $[\cdot, \cdot]: \mathfrak{g} 	imes \mathfrak{g} 	o \mathfrak{g}$ (the Lie bracket) satisfying:
    1. **Antisymmetry**: $[X, Y] = -[Y, X]$
    2. **Jacobi Identity**: $[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$

!!! theorem "Theorem 55B.3 (Engel's Theorem)"
    A finite-dimensional Lie algebra consists of nilpotent operators (in its adjoint representation) if and only if the algebra is **nilpotent**.

---

## Exercises

1. **[Fundamentals] Verify that $M_n(K)$ with the commutator bracket $[A, B] = AB - BA$ satisfies the Jacobi identity.**
   ??? success "Solution"
       Expand the bracket: $[A, [B, C]] = A(BC-CB) - (BC-CB)A = ABC - ACB - BCA + CBA$. Summing all cyclic permutations of $(A, B, C)$ shows that all 12 terms cancel in pairs (e.g., $ABC$ cancels with $-ABC$ from $[B, [C, A]]$).

2. **[Structure Constants] If $[e_i, e_j] = \sum c_{ij}^k e_k$, how does antisymmetry constrain the structure constants $c_{ij}^k$?**
   ??? success "Solution"
       $c_{ij}^k = -c_{ji}^k$. This reduces the number of independent parameters needed to define the Lie algebra structure.

3. **[Nilpotency] Show that the algebra of strictly upper triangular matrices is nilpotent.**
   ??? success "Solution"
       The bracket of two strictly upper triangular matrices increases the number of zero super-diagonals. After $n$ steps of nested brackets, the resulting matrix is zero. This is the prototypical example of a nilpotent Lie algebra.

4. **[Adjoint] Define the adjoint operator $\operatorname{ad}_X: \mathfrak{g} 	o \mathfrak{g}$.**
   ??? success "Solution"
       $\operatorname{ad}_X(Y) = [X, Y]$. The map $X \mapsto \operatorname{ad}_X$ is a representation of the Lie algebra into its own endomorphism space $\mathfrak{gl}(\mathfrak{g})$, called the **adjoint representation**.

5. **[Killing Form] What is the Killing form $B(X, Y)$?**
   ??? success "Solution"
       $B(X, Y) = \operatorname{tr}(\operatorname{ad}_X \circ \operatorname{ad}_Y)$. This symmetric bilinear form provides a geometric metric on the Lie algebra. Cartan proved that an algebra is semisimple iff its Killing form is non-degenerate.

6. **[Solvability] State Lie's Theorem for solvable Lie algebras over $\mathbb{C}$.**
   ??? success "Solution"
       Every finite-dimensional representation of a solvable Lie algebra has a common eigenvector. This implies that all matrices in the algebra can be simultaneously upper-triangularized in some basis.

7. **[Semisimplicity] Define a semisimple Lie algebra.**
   ??? success "Solution"
       An algebra is semisimple if it has no non-zero solvable ideals. Such algebras (like $\mathfrak{sl}_n$) are the "atoms" of Lie theory and can be completely classified using root systems and Dynkin diagrams.

8. **[Calculation] Compute $[H, E]$ in $\mathfrak{sl}_2(\mathbb{C})$ where $H = \begin{pmatrix} 1 & 0 \ 0 & -1 \end{pmatrix}, E = \begin{pmatrix} 0 & 1 \ 0 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $[H, E] = HE - EH = \begin{pmatrix} 0 & 1 \ 0 & 0 \end{pmatrix} - \begin{pmatrix} 0 & -1 \ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \ 0 & 0 \end{pmatrix} = 2E$. This shows $E$ is a root vector for $H$ with root (eigenvalue) 2.

9. **[Cartan Subalgebra] What is a Cartan subalgebra $\mathfrak{h} \subseteq \mathfrak{g}$?**
   ??? success "Solution"
       It is a maximal commutative subalgebra consisting of diagonalizable elements. It acts as the "diagonal" part of the algebra, allowing for the decomposition of $\mathfrak{g}$ into root spaces.

10. **[Representations] Relate Lie algebra representations to Lie group representations.**
    ??? success "Solution"
        If $\pi: G 	o GL(V)$ is a group representation, then its differential $d\pi: \mathfrak{g} 	o \mathfrak{gl}(V)$ is a Lie algebra representation. For simply connected groups, this correspondence is a bijection, allowing us to classify group symmetries using linear algebra.

## Chapter Summary

This chapter explores the linear framework of infinitesimal symmetry:

1. **Algebraic Foundation**: Defined Lie algebras via the bracket operation and Jacobi identity, linearizing the concept of non-commutativity.
2. **Structural Hierarchies**: Classified algebras into nilpotent, solvable, and semisimple types based on their derived series and ideals.
3. **Trace Metrics**: Introduced the Killing form as the definitive tool for assessing semisimplicity and global algebraic structure.
4. **Spectral Analysis**: Linked the study of representations to the geometry of root systems and common eigenvectors (Lie's Theorem).
