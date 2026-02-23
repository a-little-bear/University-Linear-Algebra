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

1. **[Axioms] Verify if the set of all $2 \times 2$ matrices $M_{2,2}$ is a vector space.**
   ??? success "Solution"
       Yes. Matrix addition is commutative and associative, there is a zero matrix, every matrix has an inverse, and scalar multiplication distributes.

2. **[Subspace] Determine if $W = \{(x, y) : x \ge 0\}$ is a subspace of $\mathbb{R}^2$.**
   ??? success "Solution"
       No. Subspaces must be closed under scalar multiplication. If $(1, 0) \in W$ and $c = -1$, then $c(1, 0) = (-1, 0) \notin W$.

3. **[Independence] Are $(1, 0), (0, 1), (1, 1)$ linearly independent?**
   ??? success "Solution"
       No. Because $(1, 1) = (1, 0) + (0, 1)$, they are linearly dependent.

4. **[Basis] Find the dimension of the subspace spanned by $(1, 0, 0)$ and $(0, 1, 0)$ in $\mathbb{R}^3$.**
   ??? success "Solution"
       The dimension is 2. These vectors are independent and form a basis for that plane.

5. **[Nullspace] Find the dimension of the nullspace $N(A)$ for $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $\operatorname{rank}(A) = 1$. By the rank-nullity theorem: $\dim N(A) = 2 - 1 = 1$.

6. **[Row Space] Does the column space of $A$ have the same dimension as its row space?**
   ??? success "Solution"
       Yes. This is the rank theorem: $\operatorname{rank}(A) = \operatorname{rank}(A^T)$.

7. **[Coordinates] Find the coordinates of $(2, 0)$ relative to the basis $B = \{(1, 1), (1, -1)\}$.**
   ??? success "Solution"
       Set $(2, 0) = c_1(1, 1) + c_2(1, -1)$. Solving gives $c_1=1, c_2=1$. Coordinates are $(1, 1)_B$.

8. **[Polynomials] What is the standard basis for $P_2$ (polynomials of degree < 2)?**
   ??? success "Solution"
       The standard basis is $\{1, x\}$.

9. **[Rank] If $A$ is $3 \times 5$ and $\operatorname{rank}(A)=3$, what is the nullity?**
   ??? success "Solution"
       $\dim N(A) = 5 - 3 = 2$.

10. **[Intersection] Prove the intersection of two subspaces is a subspace.**
    ??? success "Solution"
        If $\mathbf{u}, \mathbf{v} \in W_1 \cap W_2$, they are in both $W_1$ and $W_2$. Since both are subspaces, their sum and scalar multiples are in both, hence in the intersection.

## Chapter Summary

Vector spaces establish the geometric and algebraic foundation of linear algebra:

1.  **Power of Abstraction**: Axiomatization allows us to treat disparate mathematical objects (polynomials, signals) as the same type of entity (vectors).
2.  **Structural Analysis**: The concepts of basis and dimension quantify the "capacity" and "orientation" of a space, establishing the minimal coordinate system for description.
3.  **Matrix Depth**: The classification of the four fundamental subspaces reveals the deepest internal construction of a linear operator $A$ and the necessity of dimension conservation between input and output.
