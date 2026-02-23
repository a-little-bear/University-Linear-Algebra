# Chapter 04: Vector Spaces

<div class="context-flow" markdown>

**Prerequisites**: Systems of Linear Equations (Ch1)

**Chapter Outline**: Definition of Abstract Vector Spaces → Subspaces → Linear Combinations and Span → Linear Independence → Basis and Dimension → Change of Basis → Four Fundamental Subspaces (Rank-Nullity Theorem)

**Extension**: Vector spaces abstract the concept of "arrows," allowing us to treat functions, signals, and matrices just like geometric vectors.

</div>

Vector spaces are the stage for linear algebra. Starting from this chapter, we no longer focus solely on arrays of numbers, but move toward studying sets with specific algebraic structures. The concepts of linear independence and basis establish the "skeleton" of a space.

---

## 04.1 Definitions and Core Concepts

!!! definition "Definition 04.1 (Vector Space)"
    A set $V$ equipped with addition and scalar multiplication is a vector space over a field $F$ if it satisfies 8 axioms (associativity, commutativity, distributivity, etc.).

!!! theorem "Theorem 04.1 (Rank-Nullity Theorem)"
    For an $m \times n$ matrix $A$:
    $$\operatorname{rank}(A) + \operatorname{nullity}(A) = n$$
    The dimension of the column space plus the dimension of the nullspace equals the number of columns.

---

## Exercises

1. **[Axioms] Verify if the set of all $2 \times 2$ matrices $M_{2,2}$ is a vector space.**
   ??? success "Solution"
       Yes. Matrix addition is commutative and associative, there is a zero matrix, every matrix has an additive inverse, and scalar multiplication distributes. All vector space axioms hold.

2. **[Subspace Criteria] Determine if $W = \{(x, y) : x \ge 0\}$ is a subspace of $\mathbb{R}^2$.**
   ??? success "Solution"
       No. A subspace must be closed under scalar multiplication. If we take $(1, 0) \in W$ and scalar $c = -1$, then $c(1, 0) = (-1, 0) \notin W$. Thus $W$ is not a subspace.

3. **[Linear Independence] Determine if $v_1 = (1, 0), v_2 = (0, 1), v_3 = (1, 1)$ are linearly independent.**
   ??? success "Solution"
       Linearly dependent. Notice $v_3 = v_1 + v_2$. There exists a non-trivial linear combination that equals the zero vector.

4. **[Basis and Dimension] Find the dimension of the subspace of $\mathbb{R}^3$ spanned by $(1, 0, 0)$ and $(0, 1, 0)$.**
   ??? success "Solution"
       The dimension is 2. These two vectors are linearly independent and form a basis for that plane.

5. **[Nullspace] Find the nullspace $N(A)$ for $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Solve $Ax=0$, which is $x_1 + x_2 = 0 \implies x_1 = -x_2$.
       The nullspace is $\{ c(1, -1)^T : c \in \mathbb{R} \}$. Its dimension (nullity) is 1.

6. **[Column Space] How is the column space $C(A)$ related to the consistency of $Ax=b$?**
   ??? success "Solution"
       The system $Ax=b$ has a solution if and only if the vector $b$ belongs to the column space $C(A)$.

7. **[Coordinates] What are the coordinates of $v = (2, 0)$ in the basis $B = \{(1, 1), (1, -1)\}$?**
   ??? success "Solution"
       Solve $c_1(1, 1) + c_2(1, -1) = (2, 0)$.
       Adding: $2c_1=2 \implies c_1=1$. Subtracting: $2c_2=2 \implies c_2=1$.
       The coordinates are $[1, 1]_B^T$.

8. **[Rank] Prove: If $A$ is a $3 \times 5$ matrix, what is its maximum rank? What is its minimum nullity?**
   ??? success "Solution"
       Rank is limited by rows and columns, so $\operatorname{rank}(A) \le \min(3, 5) = 3$.
       By the rank-nullity theorem, $\operatorname{nullity}(A) = 5 - \operatorname{rank}(A) \ge 5 - 3 = 2$.

9. **[Polynomial Spaces] Prove that the dimension of $P_n$ (polynomials of degree less than $n$) is $n$.**
   ??? success "Solution"
       The standard basis is $\{1, x, x^2, \dots, x^{n-1}\}$. Since there are $n$ elements in the basis, the dimension is $n$.

10. **[Sum and Union] If $U, W$ are subspaces, is $U \cup W$ always a subspace?**
    ??? success "Solution"
        No. For example, the $x$-axis and $y$-axis in $\mathbb{R}^2$. Their union does not contain $(1, 1)$ (since $(1, 0) + (0, 1) = (1, 1)$ is not on the axes). A union is a subspace only if one contains the other.

## Chapter Summary

Vector spaces upgrade geometric intuition into algebraic rigor:

1. **Structural Abstraction**: Any objects satisfying linear operation rules (functions, matrices) are vectors.
2. **Basis and Dimension**: Dimension is the "degrees of freedom," and a basis provides the "navigation coordinates."
3. **Subspace Analysis**: The four fundamental subspaces reveal the underlying logic of linear systems' inputs and outputs.
