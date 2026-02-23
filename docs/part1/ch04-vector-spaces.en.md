# Chapter 04: Vector Spaces

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Fields (Ch00) · Set Theory

**Chapter Outline**: Definition of Vector Spaces and Subspaces → Axioms → Linear Combinations and Span → Linear Independence → Bases and Dimension → Change of Basis Matrix → The Four Fundamental Subspaces → Direct Sums and Complements → Coordinates and Isomorphisms

**Extension**: Vector spaces are the stage upon which all of linear algebra is performed; the concept of a basis allows us to assign coordinates to abstract objects like functions or polynomials.

</div>

Linear algebra is the study of **vector spaces** and the maps between them. A vector space is more than just a collection of arrows; it is any set of objects that can be added and scaled according to specific axioms. This abstraction allows us to treat functions, matrices, and signals as "vectors." The most important concept in this chapter is the **basis**, a minimal set of building blocks that uniquely describes every element in the space.

---

## 04.1 Definitions and Fundamental Concepts

!!! definition "Definition 04.1 (Vector Space)"
    A vector space $V$ over a field $F$ is a set equipped with addition ($+$) and scalar multiplication ($\cdot$) such that $(V, +)$ is an Abelian group and scalar multiplication distributes over addition and scaling.

!!! theorem "Theorem 04.1 (Basis Uniqueness)"
    Every basis of a finite-dimensional vector space has the same number of vectors. This number is the **dimension** ($\dim V$) of the space.

---

## Exercises

1. **[Subspaces] Show that the set of all symmetric $n \times n$ matrices is a subspace of $M_n(\mathbb{R})$.**
   ??? success "Solution"
       If $A, B$ are symmetric, $(A+B)^T = A^T + B^T = A+B$, and $(kA)^T = k A^T = kA$. The zero matrix is symmetric. Since it is closed under addition and scaling, it is a subspace.

2. **[Linear Independence] Are the vectors $\{(1, 0), (0, 1), (1, 1)\}$ linearly independent in $\mathbb{R}^2$?**
   ??? success "Solution"
       No. The third vector is the sum of the first two. Since $v_3 = 1 v_1 + 1 v_2$, the set is linearly dependent.

3. **[Basis] Find a basis for the space of polynomials of degree at most 2.**
   ??? success "Solution"
       The standard basis is $\{1, x, x^2\}$. Any polynomial $ax^2 + bx + c$ is a unique linear combination of these elements.

4. **[Dimension] What is the dimension of the subspace of $\mathbb{R}^3$ defined by $x+y+z=0$?**
   ??? success "Solution"
       The equation represents a plane through the origin. Since there is one constraint on 3 variables, the dimension is $3 - 1 = 2$.

5. **[Span] Describe the span of $\{(1, 0, 0), (0, 1, 0)\}$.**
   ??? success "Solution"
       The span is the $xy$-plane in $\mathbb{R}^3$, i.e., the set of all vectors $(x, y, 0)$.

6. **[Four Subspaces] Define the null space $\ker(A)$ and the column space $\operatorname{Im}(A)$.**
   ??? success "Solution"
       $\ker(A) = \{x : Ax = 0\}$. $\operatorname{Im}(A) = \{Ax : x \in \mathbb{R}^n\}$. These characterize the solution set and the range of the operator $A$.

7. **[Change of Basis] If $[v]_B = (1, 2)^T$ and the basis change matrix from $B$ to $C$ is $\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$, find $[v]_C$.**
   ??? success "Solution"
       $[v]_C = P_{C \leftarrow B} [v]_B = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 2 \end{pmatrix} = \begin{pmatrix} 3 \\ 2 \end{pmatrix}$.

8. **[Direct Sum] When is a vector space the direct sum $V = U \oplus W$?**
   ??? success "Solution"
       When every $v \in V$ can be uniquely written as $u + w$ ($u \in U, w \in W$). This requires $U + W = V$ and $U \cap W = \{0\}$.

9. **[Coordinate Mapping] Explain why every $n$-dimensional vector space over $F$ is isomorphic to $F^n$.**
   ??? success "Solution"
       Fixing a basis $\{b_1, \dots, b_n\}$ allows us to map each vector $v = \sum x_i b_i$ to its coordinate vector $(x_1, \dots, x_n) \in F^n$. This mapping is a linear bijection.

10. **[Rank-Nullity] If $A$ is a $4 \times 7$ matrix with $\operatorname{rank}(A) = 3$, what is $\dim \ker(A)$?**
    ??? success "Solution"
        By the Rank-Nullity Theorem, $\operatorname{rank}(A) + \dim \ker(A) = n$. Thus $3 + \dim \ker(A) = 7$, so $\dim \ker(A) = 4$.

## Chapter Summary

This chapter provides the axiomatic framework for linear structures:

1. **Axiomatic Abstraction**: Defined vector spaces through operational laws, extending linear analysis to matrices and functions.
2. **Structural Composition**: Developed the concepts of span and linear independence to describe the reach and redundancy of vector sets.
3. **Dimensional Scaling**: Established the basis as the unique coordinate system for finite-dimensional spaces.
4. **Subspace Analysis**: Linked matrix rank to the dimensions of the four fundamental subspaces, providing a geometric view of solvability.
