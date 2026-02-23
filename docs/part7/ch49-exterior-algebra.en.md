# Chapter 49: Exterior Algebra and Grassmannians

<div class="context-flow" markdown>

**Prerequisites**: Multilinear Algebra (Ch21) · Determinants (Ch03) · Subspaces and Bases (Ch04)

**Chapter Outline**: From Multilinearity to Antisymmetry → Definition and Axioms of the Exterior Product ($\wedge$) → Exterior Algebra (Grassmann Algebra) $\Lambda(V)$ → Basis Construction and Dimension $\binom{n}{k}$ → The Exterior Origin of Determinants → Decomposable Tensors and their Correspondence to Subspaces → Algebraic Description of the Grassmannian Manifold $Gr(k, V)$ → Plücker Coordinates and Plücker Relations → Applications: Differential Forms, Electromagnetic Field Tensors, Subspace Clustering, and Projective Geometry

**Extension**: Exterior algebra is the natural language for describing "oriented volumes"; by introducing antisymmetry, it treats subspaces themselves as algebraic objects. It is the core of modern differential geometry, general relativity, and network flow analysis in theoretical computer science.

</div>

In classical algebra, we study combinations of vectors. In **Exterior Algebra**, we elevate our perspective to "parallelepipeds spanned by vectors." By introducing the **Exterior Product** (wedge product) $\wedge$, we can handle areas, volumes, and their generalizations in a purely algebraic manner. This theory ultimately leads to the **Grassmannian**—the manifold of all $k$-dimensional subspaces—providing the ultimate framework for "subspace operations" in operator spaces.

---

## 49.1 The Exterior Product and Algebra

!!! definition "Definition 49.1 (Exterior Product)"
    The **exterior product** $v \wedge w$ of vectors $v, w \in V$ satisfies:
    1.  **Bilinearity**: $(av_1+bv_2) \wedge w = a(v_1 \wedge w) + b(v_2 \wedge w)$.
    2.  **Antisymmetry**: $v \wedge w = -w \wedge v$.
    3.  **Nilpotency**: $v \wedge v = 0$.

!!! note "Exterior Algebra $\Lambda(V)$"
    The direct sum of all $k$-th exterior powers $\Lambda^k(V)$ forms the exterior algebra. Its total dimension is $2^n$ for an $n$-dimensional space.

---

## 49.2 Geometric Meaning and Subspaces

!!! theorem "Theorem 49.1 (Linear Independence Criterion)"
    The set $\{v_1, \ldots, v_k\}$ is linearly independent iff $v_1 \wedge v_2 \wedge \cdots \wedge v_k \neq 0$.
    **Geometric Intuition**: A non-zero $k$-th order wedge product represents a "subspace fragment" with a specific orientation and $k$-dimensional volume.

---

## 49.3 Grassmannians and Plücker Coordinates

!!! definition "Definition 49.2 (Grassmannian)"
    $Gr(k, V)$ is the set of all $k$-dimensional subspaces of $V$.
    Each $k$-dimensional subspace can be uniquely represented (up to a scalar) by the exterior product of its basis vectors $v_1 \wedge \cdots \wedge v_k$ (a decomposable tensor).

!!! technique "Plücker Coordinates"
    Expanding $v_1 \wedge \cdots \wedge v_k$ in the basis of $\Lambda^k(V)$ yields coefficients called **Plücker Coordinates**. They satisfy a set of quadratic equations known as **Plücker Relations**.

---

## Exercises

**1. [Basics] Calculate $(e_1 + e_2) \wedge (e_1 - e_2)$.**

??? success "Solution"
    **Steps:**
    1. Apply distributivity: $e_1 \wedge e_1 - e_1 \wedge e_2 + e_2 \wedge e_1 - e_2 \wedge e_2$.
    2. Apply nilpotency: $e_1 \wedge e_1 = 0$ and $e_2 \wedge e_2 = 0$.
    3. Apply antisymmetry: $e_2 \wedge e_1 = -e_1 \wedge e_2$.
    4. Substitute: $0 - e_1 \wedge e_2 - e_1 \wedge e_2 - 0 = -2 e_1 \wedge e_2$.
    **Conclusion**: The result is $-2 e_1 \wedge e_2$.

**2. [Dimension] If $\dim V = 4$, find the dimension of $\Lambda^2(V)$ and list the basis.**

??? success "Solution"
    **Calculation:**
    1. Dimension is $\binom{4}{2} = \frac{4 \times 3}{2 \times 1} = 6$.
    2. Let the basis of $V$ be $\{e_1, e_2, e_3, e_4\}$.
    **Basis**: $\{e_1 \wedge e_2, e_1 \wedge e_3, e_1 \wedge e_4, e_2 \wedge e_3, e_2 \wedge e_4, e_3 \wedge e_4\}$.

**3. [Determinant] Prove that the determinant of an $n \times n$ matrix is the wedge product of its columns.**

??? success "Solution"
    **Algebraic Mapping:**
    1. Consider $v_1 \wedge v_2 \wedge \cdots \wedge v_n$.
    2. Substitute $v_j = \sum a_{ij} e_i$ and expand using antisymmetry.
    3. Terms are non-zero only when the indices are a permutation of $\{1, \ldots, n\}$.
    4. The coefficient of each term is the sign of the permutation $\operatorname{sgn}(\sigma)$.
    **Conclusion**: $v_1 \wedge \cdots \wedge v_n = \det(A) (e_1 \wedge \cdots \wedge e_n)$. This reveals the most fundamental definition of the determinant: it is the scaling factor of the top-level exterior power.

**4. [Decomposability] Determine if $v = e_1 \wedge e_2 + e_3 \wedge e_4$ is decomposable (i.e., can be written as $u \wedge w$).**

??? success "Solution"
    **Plücker Criterion:**
    1. in 4D, a 2nd-order exterior element is decomposable iff $v \wedge v = 0$.
    2. Compute $v \wedge v = (e_1 \wedge e_2 + e_3 \wedge e_4) \wedge (e_1 \wedge e_2 + e_3 \wedge e_4)$.
    3. Expand: $e_1 \wedge e_2 \wedge e_1 \wedge e_2 + e_1 \wedge e_2 \wedge e_3 \wedge e_4 + e_3 \wedge e_4 \wedge e_1 \wedge e_2 + e_3 \wedge e_4 \wedge e_3 \wedge e_4$.
    4. Repeated terms are zero, leaving $2 e_1 \wedge e_2 \wedge e_3 \wedge e_4$.
    **Conclusion**: Since $v \wedge v \neq 0$, the element is not decomposable. It does not represent a single 2D subspace but a superposition of them.

**5. [Property] Prove: If $v_1, \ldots, v_k$ are linearly dependent, then $v_1 \wedge \cdots \wedge v_k = 0$.**

??? success "Solution"
    **Proof:**
    1. If dependent, at least one vector is a combination of others, say $v_1 = \sum_{i=2}^k c_i v_i$.
    2. Substitute into the wedge product: $(\sum c_i v_i) \wedge v_2 \wedge \cdots \wedge v_k$.
    3. Every term in the sum contains a repeated vector (e.g., $c_2 v_2 \wedge v_2 \wedge \cdots$).
    4. By nilpotency, every term is zero.

**6. [Grassmannian] To what geometric object does $Gr(1, V)$ correspond?**

??? success "Solution"
    **Conclusion: The Projective Space $P(V)$.**
    $Gr(1, V)$ represents all lines through the origin. Geometrically, this is the definition of a projective space.

**7. [Plücker Relation] Write the single Plücker relation for $Gr(2, 4)$.**

??? success "Solution"
    **Formula:**
    Let the coordinates be $p_{ij}$.
    $p_{12}p_{34} - p_{13}p_{24} + p_{14}p_{23} = 0$.
    This is the algebraic equation identifying which 6D vectors represent 2D subspaces.

**8. [Duality] What is the relationship between $\Lambda^k(V)$ and $\Lambda^{n-k}(V)$?**

??? success "Solution"
    **Conclusion: They are isomorphic.**
    **Reasoning**: Their dimensions are equal $\binom{n}{k} = \binom{n}{n-k}$. In a space with an inner product, this correspondence is given by the **Hodge Star operator** ($\star$).

**9. [Basics] In $\mathbb{R}^3$, how does $v \wedge w$ relate to the cross product $v \times w$?**

??? success "Solution"
    **Connection:**
    $v \wedge w$ is a 2nd-order tensor (representing an oriented area element) in $\Lambda^2(\mathbb{R}^3)$. $v \times w$ is a vector in $\mathbb{R}^3$. Through the Hodge dual in 3D, area elements map uniquely to their normal vectors. Thus, the cross product is the dual manifestation of the wedge product in 3D.

**10. [Application] Briefly state the application of exterior algebra in electromagnetism.**

??? success "Solution"
    In relativistic electrodynamics, the electric and magnetic fields are unified into a single 2nd-order antisymmetric tensor (the Faraday tensor $F$), which is an element of the exterior algebra of 4D spacetime. Maxwell's equations can be written as $dF = 0$ and $d{\star F} = J$, beautifully demonstrating the natural advantage of exterior algebra in describing flux and circulation.

## Chapter Summary

Exterior algebra is the algebraic culmination of geometric intuition:

1.  **Algebraization of Subspaces**: Through the wedge product, we transform abstract "subspaces" into concrete "algebraic elements," realizing the jump from studying points to studying fragments of space.
2.  **Operatorization of Volume**: The origin and properties of determinants are most thoroughly explained within the exterior algebra framework, proving that antisymmetry is the algebraic essence of multi-dimensional measure theory.
3.  **Foundations of Manifolds**: The Grassmannian manifold and its Plücker coordinates provide fine-grained local characterization tools for modern geometry and topology, serving as a vital link between pure algebra and modern physics (e.g., string theory).
