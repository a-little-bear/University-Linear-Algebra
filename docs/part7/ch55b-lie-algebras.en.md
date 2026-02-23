# Chapter 55B: Lie Algebras

<div class="context-flow" markdown>

**Prerequisites**: Matrix Groups (Ch55) · Matrix Analysis (Ch14) · Commutators (Ch02)

**Chapter Outline**: Abstract Definition of a Lie Algebra → The Lie Bracket: Bilinearity, Antisymmetry, and the Jacobi Identity → Matrix Lie Algebras as Tangent Spaces of Matrix Groups → Lie Algebras of Classical Groups ($\mathfrak{gl}, \mathfrak{sl}, \mathfrak{so}, \mathfrak{su}, \mathfrak{sp}$) → The Exponential and Logarithmic Mappings → Adjoint Representation $\operatorname{ad}_X(Y) = [X, Y]$ → Structure Constants → Physics Applications: Commutation Relations in Quantum Mechanics, Angular Momentum → Control Theory: Reachability in Non-holonomic Systems

**Extension**: A Lie algebra locally linearizes a "group" (a complex geometric object) into a "vector space"; it is the ultimate algebraic tool for studying continuous symmetry, supporting the mathematical formulation of everything from General Relativity to Gauge Theory.

</div>

If matrix groups describe "global transformations," then **Lie Algebras** describe "infinitesimal transformations." By studying the tangent space of a matrix group at the identity, we obtain a vector space equipped with a special multiplication operation (the Lie bracket). The superiority of Lie algebras lies in their ability to transform complex group multiplication into linear bracket operations, greatly simplifying the analysis of symmetry.

---

## 55B.1 Definition and Axioms

!!! definition "Definition 55B.1 (Lie Algebra)"
    A vector space $\mathfrak{g}$ over a field $F$ is a **Lie Algebra** if it is equipped with a binary operation $[\cdot, \cdot]: \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}$ (the **Lie Bracket**) satisfying:
    1.  **Bilinearity**: $[ax+by, z] = a[x, z] + b[y, z]$.
    2.  **Antisymmetry**: $[x, y] = -[y, x]$ (implies $[x, x] = 0$).
    3.  **Jacobi Identity**: $[x, [y, z]] + [y, [z, x]] + [z, [x, y]] = 0$.

!!! note "Matrix Lie Bracket"
    For square matrices, the **commutator** $[A, B] = AB - BA$ is the standard implementation satisfying these axioms.

---

## 55B.2 Lie Algebras of Classical Groups

!!! theorem "Theorem 55B.1 (Classification of Classical Lie Algebras)"
    1.  **$\mathfrak{gl}(n)$**: The Lie algebra of $GL(n)$, consisting of all $n \times n$ matrices.
    2.  **$\mathfrak{sl}(n)$**: The Lie algebra of $SL(n)$, consisting of **trace-zero** matrices ($\operatorname{tr}(X) = 0$).
    3.  **$\mathfrak{so}(n)$**: The Lie algebra of $O(n)$, consisting of **skew-symmetric** matrices ($X^T = -X$).
    4.  **$\mathfrak{su}(n)$**: The Lie algebra of $SU(n)$, consisting of **trace-zero skew-Hermitian** matrices ($X^* = -X, \operatorname{tr}(X) = 0$).

---

## 55B.3 Adjoint Representation and Structure Constants

!!! definition "Definition 55B.2 (Adjoint Representation $\operatorname{ad}$)"
    For $X \in \mathfrak{g}$, the operator $\operatorname{ad}_X: \mathfrak{g} \to \mathfrak{g}$ is defined as $\operatorname{ad}_X(Y) = [X, Y]$.
    This is a linear representation of the Lie algebra on itself. The trace of its matrix product is related to the **Killing Form**, used to identify semisimple algebras.

!!! technique "Structure Constants"
    Given a basis $\{E_i\}$, the bracket is $[E_i, E_j] = \sum C_{ij}^k E_k$. The constants $C_{ij}^k$ encode all local algebraic information of the Lie algebra.

---

## Exercises

1.  **[Calculation] Calculate the Lie bracket of $\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$ and $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**
    ??? success "Solution"
        $AB - BA = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} - \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ -2 & 0 \end{pmatrix}$.

2.  **[Jacobi] Verify that $[x, [y, y]] = 0$ is consistent with the Jacobi Identity.**
    ??? success "Solution"
        Since $[y, y] = 0$, the left side is 0. The Jacobi Identity trivially holds whenever one of the nested brackets is zero.

3.  **[Trace] Prove: if $A, B \in \mathfrak{sl}(n)$, then $[A, B] \in \mathfrak{sl}(n)$.**
    ??? success "Solution"
        $\operatorname{tr}(AB - BA) = \operatorname{tr}(AB) - \operatorname{tr}(BA) = 0$. Thus, the commutator is trace-zero.

4.  **[Dimension] Prove that the dimension of $\mathfrak{so}(3)$ is 3.**
    ??? success "Solution"
        A $3 \times 3$ skew-symmetric matrix has the form $\begin{pmatrix} 0 & a & b \\ -a & 0 & c \\ -b & -c & 0 \end{pmatrix}$, which is determined by the 3 independent parameters $a, b, c$.

5.  **[Quantum] Which Lie algebra corresponds to the commutation relations of the Pauli matrices?**
    ??? success "Solution"
        They correspond to $\mathfrak{su}(2)$ (which is isomorphic to $\mathfrak{so}(3)$). They represent the infinitesimal generators of angular momentum.

6.  **[Exponential] Describe the geometric meaning of the map $\exp: \mathfrak{so}(n) \to SO(n)$.**
    ??? success "Solution"
        It maps a rotation velocity (a skew-symmetric matrix) to a finite rotation matrix.

7.  **[Center] What is the center of a Lie algebra?**
    ??? success "Solution"
        The set of elements that commute with everything: $Z(\mathfrak{g}) = \{x \in \mathfrak{g} : [x, y]=0, \forall y \in \mathfrak{g}\}$.

8.  **[Adjoint] Prove $\operatorname{ad}_{[X, Y]} = [\operatorname{ad}_X, \operatorname{ad}_Y]$.**
    ??? success "Solution"
        This is equivalent to the Jacobi Identity: $[[X, Y], Z] = [X, [Y, Z]] - [Y, [X, Z]]$.

9.  **[Abelian] What type of Lie algebra has all structure constants equal to 0?**
    ??? success "Solution"
        An **Abelian Lie algebra**. In this case, the Lie bracket is always zero, corresponding to a commutative group.

****

??? success "Solution"
    

## Chapter Summary

Lie algebras are the highest form of linear algebra in the study of symmetry:

1.  **Analytic Power of the Infinitesimal**: Through tangent space techniques, they reduce complex non-linear group manifolds to linear vector spaces, making the classification of symmetries computable.
2.  **Unity of Operation**: The Lie bracket abstracts the essence of the matrix commutator, establishing algebraic rules for the interference between physical operators (such as momentum and spin).
3.  **Structural Mapping**: Through the exponential map and adjoint representation, Lie algebras establish a "dynamic-replaces-static" analysis mode, proving that local infinitesimal changes are sufficient to determine global transformation laws.
