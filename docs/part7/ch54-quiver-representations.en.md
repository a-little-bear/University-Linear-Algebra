# Chapter 54: Quiver Representations

<div class="context-flow" markdown>

**Prerequisites**: Graph Theory (Ch27) · Linear Transformations (Ch05) · Direct Sums of Vector Spaces (Ch04)

**Chapter Outline**: Definition of a Quiver (Directed Graph) → Definition of a Quiver Representation (Spaces on Vertices, Maps on Edges) → Morphisms and Isomorphisms → Subrepresentations and Direct Sums → Indecomposable Representations → Gabriel's Theorem: Finite Type Quivers and ADE Dynkin Diagrams → Dimension Vectors → Path Algebras → Applications: Representation Theory and Persistent Homology in Topological Data Analysis (TDA)

**Extension**: Quiver representation theory is a masterpiece of modern algebra; it translates abstract representation problems into intuitive graph-theoretic and linear-algebraic problems, revealing deep combinatorial classification laws behind symmetry.

</div>

What happens when we generalize a linear transformation (a mapping from a space to itself) to a complex network structure? **Quiver Representations** provide the framework for this study. By placing vector spaces at the vertices of a directed graph and linear mappings on the edges, we obtain an algebraic object that is both powerful and intuitive. This chapter will prove that these complex network structures can ultimately be reduced to several perfect symmetric forms.

---

## 54.1 Definition of Quivers and Representations

!!! definition "Definition 54.1 (Quiver)"
    A **Quiver** $Q = (V, E)$ is simply a directed graph, where $V$ is the set of vertices and $E$ is the set of directed edges (arrows). Self-loops and multiple edges are allowed.

!!! definition "Definition 54.2 (Representation of a Quiver)"
    A **representation** $V$ of a quiver $Q$ consists of:
    1.  For each vertex $i \in V$, a vector space $V_i$.
    2.  For each arrow $a: i \to j$, a linear transformation $V_a: V_i \to V_j$.

---

## 54.2 Indecomposability and Direct Sums

!!! definition "Definition 54.3 (Direct Sum and Indecomposability)"
    - The **direct sum** of two representations is obtained by taking the direct sum of the corresponding spaces and mappings.
    - A representation is **indecomposable** if it cannot be written as the direct sum of two non-zero subrepresentations.
    **Goal**: Classify all possible indecomposable representations.

---

## 54.3 Gabriel's Theorem

!!! theorem "Theorem 54.1 (Gabriel's Theorem)"
    A connected quiver $Q$ has finitely many indecomposable representations (is of **finite type**) if and only if its underlying graph (ignoring direction) is a **Dynkin diagram of type A, D, or E**:
    - **$A_n$**: A simple chain.
    - **$D_n$**: A chain with one branch.
    - **$E_6, E_7, E_8$**: Specific star-shaped diagrams.
    In this case, the dimension vectors of indecomposable representations correspond bijectively to the positive roots of the associated root system.

---

## 54.4 Path Algebras

!!! technique "Algebraic Perspective"
    The **Path Algebra** $kQ$ of a quiver $Q$ is the algebra formed by all possible paths in the graph. The representation theory of $Q$ is equivalent to the module theory of the path algebra $kQ$. This provides a bridge between geometric structures and pure algebra.

---

## Exercises

1.  **[Basics] Provide a non-zero representation for the quiver $A_2$ ($1 \to 2$).**
    ??? success "Solution"
        Assign $V_1 = \mathbb{C}$, $V_2 = \mathbb{C}$, and the edge map $V_a = I$ (the identity map).

2.  **[Dimension] What is the dimension vector of a representation?**
    ??? success "Solution"
        It is the vector $\mathbf{d} = (\dim V_1, \dim V_2, \ldots, \dim V_n)$ consisting of the dimensions of the vector spaces at each vertex.

3.  **[Indecomposable] Prove that the representation $V: \mathbb{C} \xrightarrow{1} \mathbb{C}$ of $A_2$ is indecomposable.**
    ??? success "Solution"
        This is a 2-dimensional representation. If it were decomposable, it would be the sum of two 1-dimensional representations. However, 1D representations are either $\mathbb{C} \to 0$ or $0 \to \mathbb{C}$. Their sum must have a zero map, which contradicts the original representation.

4.  **[Classification] How many finite-type quivers correspond to the Dynkin diagram $A_3$ ($1-2-3$)?**
    ??? success "Solution"
        There are $2^2 = 4$ choices of directions: $1 \to 2 \to 3$, $1 \leftarrow 2 \leftarrow 3$, $1 \leftarrow 2 \to 3$, and $1 \to 2 \leftarrow 3$.

5.  **[Gabriel] Why is the representation theory of $A_n$ type quivers relatively simple?**
    ??? success "Solution"
        Because it corresponds to the decomposition of a sequence of linear transformations, essentially a generalization of matrix chains or the SVD problem.

6.  **[Morphism] What commutative diagram must a morphism $\phi: V \to W$ satisfy?**
    ??? success "Solution"
        For every edge $a: i \to j$, the diagram must commute: $W_a \phi_i = \phi_j V_a$.

7.  **[Zero] What is the zero representation?**
    ??? success "Solution"
        A representation where all vertex spaces are the zero space $\{0\}$.

8.  **[Application] Briefly describe the role of quiver representations in TDA.**
    ??? success "Solution"
        Persistent homology analyzes a sequence of nested complexes (an $A_n$ type quiver representation) and uses the "barcodes" of indecomposable representations to identify topological features of data.

9.  **[Loops] Prove that a quiver with a loop (self-loop or cycle) cannot be of finite type.**
    ??? success "Solution"
        By placing maps with different eigenvalues $\lambda$ on the loop, one can construct infinitely many non-isomorphic indecomposable representations.

10. **[Dynkin] Which Dynkin diagram corresponds to a set of mutually orthogonal vectors?**

   ??? success "Solution"
        Type $A_1$ (isolated vertices), because there is no mapping coupling between the vertices.

## Chapter Summary

Quiver representation theory is a perfect symphony of combinatorics and linear algebra:

1.  **Networked Mappings**: It extends isolated linear transformations into complex topological networks, establishing a universal algebraic framework for describing multi-object interactions.
2.  **Taxonomy of Symmetry**: Gabriel's theorem reveals the geometric "harmony" behind algebraic structures, reducing complex classifications to classic Dynkin diagrams—a hallmark of unified thought in mathematics.
3.  **Topologizing Computation**: Through path algebras and dimension vectors, quiver theory provides solid theoretical support for modern topological data analysis, proving that "algebraic stability" is "geometric reality."
