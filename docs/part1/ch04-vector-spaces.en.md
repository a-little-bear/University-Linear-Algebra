# Chapter 04: Vector Spaces

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Matrix Algebra (Ch02)

**Chapter Outline**: Definition of Abstract Vector Spaces → The 8 Axioms → Examples (Coordinate Spaces, Matrix Spaces, Polynomials, Function Spaces) → Subspace Criteria → Linear Combinations and Span → Linear Independence/Dependence → Basis and Dimension → Coordinate Vectors → Four Fundamental Subspaces of a Matrix → Rank-Nullity Theorem → Change of Basis and Transition Matrices

**Extension**: Vector spaces abstract the concept of "arrows," allowing us to treat physical quantities, signals, and polynomials under a unified framework; it is the stage upon which Linear Transformations (Ch05) act.

</div>

If matrices are the skeleton of linear algebra, then vector spaces are its soul. By abstracting specific calculation rules into axioms, we gain a universal method for handling any objects that satisfy the principle of linear superposition. This chapter begins with rigorous definitions and systematically constructs the "scaffold" of a space: basis and dimension.

---

## 04.1 Vector Spaces and Axioms

!!! definition "Definition 04.1 (Vector Space)"
    A set $V$ is called a **vector space** over a field $F$ if it is equipped with addition and scalar multiplication satisfying the 8 axioms:
    1.  **Commutativity and Associativity of Addition**.
    2.  **Existence of Zero Vector $\mathbf{0}$**.
    3.  **Existence of Additive Inverses**.
    4.  **Associativity of Scalar Multiplication**.
    5.  **Unit Scalar Identity**: $1\mathbf{v} = \mathbf{v}$.
    6.  **Distributivity over Scalars**: $c(\mathbf{u}+\mathbf{v}) = c\mathbf{u} + c\mathbf{v}$.
    7.  **Distributivity over Vectors**: $(a+b)\mathbf{v} = a\mathbf{v} + b\mathbf{v}$.

---

## 04.2 Subspaces and Span

!!! definition "Definition 04.2 (Subspace)"
    A non-empty subset $W$ of $V$ is a **subspace** if $W$ is **closed** under addition and scalar multiplication:
    - If $\mathbf{u}, \mathbf{v} \in W$, then $\mathbf{u} + \mathbf{v} \in W$.
    - If $\mathbf{u} \in W$ and $c \in F$, then $c\mathbf{u} \in W$.

!!! definition "Definition 04.3 (Span)"
    The set of all linear combinations of a set of vectors $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ is called their **span**:
    $$\operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_k\} = \{c_1\mathbf{v}_1 + \cdots + c_k\mathbf{v}_k : c_i \in F\}$$

---

## 04.3 Basis and Dimension

!!! definition "Definition 04.4 (Linear Independence)"
    A set of vectors is **linearly independent** if $c_1\mathbf{v}_1 + \cdots + c_k\mathbf{v}_k = \mathbf{0}$ implies $c_i = 0$ for all $i$.

!!! definition "Definition 04.5 (Basis and Dimension)"
    If a set $B = \{\mathbf{b}_1, \ldots, \mathbf{b}_n\}$ is linearly independent and spans $V$, then $B$ is a **basis** for $V$. The number of vectors in a basis is the **dimension** $\dim(V)$.

---

## 04.4 Four Fundamental Subspaces

!!! theorem "Theorem 04.1 (Four Fundamental Subspaces)"
    For an $m \times n$ matrix $A$:
    1.  **Column Space $C(A)$**: Spanned by the columns of $A$, in $\mathbb{R}^m$.
    2.  **Nullspace $N(A)$**: All $x$ such that $Ax = \mathbf{0}$, in $\mathbb{R}^n$.
    3.  **Row Space $C(A^T)$**: Spanned by the rows of $A$, in $\mathbb{R}^n$.
    4.  **Left Nullspace $N(A^T)$**: All $y$ such that $A^T y = \mathbf{0}$, in $\mathbb{R}^m$.

!!! theorem "Theorem 04.2 (Rank-Nullity Theorem)"
    $\dim C(A) + \dim N(A) = n$ (number of columns). That is: **Rank + Nullity = Total Columns**.

---

## Exercises

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

Vector spaces establish the geometric and algebraic foundation of linear algebra:

1.  **Power of Abstraction**: Axiomatization allows us to treat disparate mathematical objects (polynomials, signals) as the same type of entity (vectors).
2.  **Structural Analysis**: The concepts of basis and dimension quantify the "capacity" and "orientation" of a space, establishing the minimal coordinate system for description.
3.  **Matrix Depth**: The classification of the four fundamental subspaces reveals the deepest internal construction of a linear operator $A$ and the necessity of dimension conservation between input and output.
