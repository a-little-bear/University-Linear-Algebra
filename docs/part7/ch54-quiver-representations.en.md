# Chapter 54: Quiver Representations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Linear Transformations (Ch05) · Graph Theory Basics (Ch27)

**Chapter Outline**: From Graphs to Algebra → Definition of a Quiver (Directed Graph) → Representation of a Quiver: Spaces at Vertices and Maps on Edges → Homomorphisms and Isomorphisms of Representations → Path Algebra $kQ$ → Direct Sums and Irreducible Representations → Key Theorem: Gabriel’s Theorem (Relationship between Finite Representation Type and Dynkin Diagrams) → Projective and Injective Representations → Applications: Root Systems, Gauge Theory in Physics (String Theory), and the Category Theoretic Perspective

**Extension**: Quiver representation theory is a cutting-edge branch of modern algebra; it chains isolated linear transformations into complex network structures, proving that the topological shape of a graph directly dictates the "classification difficulty" of its linear operators. It is a perfect bridge between discrete graph theory and continuous transformation theory.

</div>

In traditional linear algebra, we typically study a single transformation $T: V \to V$ on a single space (this corresponds to a graph with one vertex and one loop). **Quiver Representations** study more complex systems: a network of multiple vector spaces and a set of linear maps between them. This theory not only unifies classic problems like matrix similarity and congruence but also reveals astonishing links between graph topology (such as Dynkin diagrams) and the classification of linear operators.

---

## 54.1 Definition of Quivers and Representations

!!! definition "Definition 54.1 (Quiver)"
    A **Quiver** $Q$ is a directed graph consisting of a set of vertices $Q_0$ and a set of edges $Q_1$.

!!! definition "Definition 54.2 (Representation of a Quiver)"
    A representation $V$ of $Q$ over a field $k$ consists of:
    1.  A vector space $V_i$ for each vertex $i \in Q_0$.
    2.  A linear map $V_a: V_i \to V_j$ for each arrow $a: i \to j$.

---

## 54.2 Core Theorem: Gabriel's Theorem

!!! theorem "Theorem 54.1 (Gabriel’s Theorem)"
    A quiver has finitely many indecomposable representations (is of finite type) iff its underlying undirected graph is a **Dynkin Diagram** ($A_n, D_n, E_6, E_7, E_8$).
    **Significance**: This profound result links the discrete geometry of graphs to the classification complexity of linear algebra.

---

## 54.3 Path Algebras $kQ$

!!! technique "Algebraic Transformation"
    The category of representations of a quiver $Q$ is equivalent to the category of left modules over its **Path Algebra** $kQ$. This allows the use of ring and module theory to study complex mapping networks.

---

## Exercises

**1. [Basics] Construct a non-trivial representation of the quiver $1 \to 2$.**

??? success "Solution"
    **Construction:**
    1. Assign space $V_1 = k$ to vertex 1.
    2. Assign space $V_2 = k$ to vertex 2.
    3. Assign the identity map $V_a = [1]$ to arrow $a: 1 \to 2$.
    **Conclusion**: This is the simplest indecomposable representation, denoted $P_1$ (a projective representation).

**2. [Calculation] Find all indecomposable representations of the quiver $1 \to 2$.**

??? success "Solution"
    **Analysis:**
    1. The underlying graph is $A_2$ (a Dynkin diagram).
    2. By Gabriel's theorem, it has finite type.
    3. The dimension vectors $(d_1, d_2)$ must be positive roots of $A_2$: $(1, 0), (0, 1), (1, 1)$.
    **Conclusion**: There are 3 indecomposables: $(k \to 0), (0 \to k)$, and $(k \xrightarrow{1} k)$.

**3. [Isomorphism] Determine if the representations $(k \xrightarrow{1} k)$ and $(k \xrightarrow{2} k)$ are isomorphic.**

??? success "Solution"
    **Determination:**
    1. Isomorphism requires non-zero scalars $c_1, c_2$ such that $c_2 \cdot 1 = 2 \cdot c_1$.
    2. Setting $c_1 = 1$ and $c_2 = 2$ satisfies this.
    **Conclusion**: Yes, they are isomorphic. Intuition: In quiver theory, pure basis scaling doesn't change the underlying structure.

**4. [Path Algebra] List the basis elements of the path algebra for $1 \xrightarrow{a} 2 \xrightarrow{b} 3$.**

??? success "Solution"
    **Enumeration:**
    1. Paths of length 0 (vertices): $e_1, e_2, e_3$.
    2. Paths of length 1 (arrows): $a, b$.
    3. Paths of length 2: $ba$.
    **Conclusion**: The basis for $kQ$ is $\{e_1, e_2, e_3, a, b, ba\}$.

**5. [Dynkin] Why is the quiver $1 \to 2 \gets 3$ of finite type?**

??? success "Solution"
    **Reasoning:**
    Its underlying undirected graph is $A_3$. According to Gabriel's theorem, all quivers of type $A_n$ are finite type, regardless of the orientation of the arrows.

**6. [Calculation] Find the sub-representations of $k \xrightarrow{0} k$.**

??? success "Solution"
    **Analysis:**
    A sub-representation $U \subset V$ requires $V_a(U_i) \subseteq U_j$. Since $V_a = 0$, any pair of subspaces $(U_1, U_2)$ works.
    **Conclusion**: Sub-representations are $(0 \to 0), (k \to 0), (0 \to k)$, and $(k \to k)$.

**7. [Classification] Briefly describe "Representation Types": Finite, Tame, and Wild.**

??? success "Solution"
    **Comparison:**
    - **Finite**: Finitely many indecomposables (Dynkin graphs).
    - **Tame**: Indecomposables come in continuous families (Euclidean graphs $\tilde{A}, \tilde{D}, \tilde{E}$).
    - **Wild**: Classification contains the problem of classifying all pairs of matrices; impossible to fully resolve.

**8. [Application] Which quiver corresponds to the matrix similarity problem?**

??? success "Solution"
    **Conclusion: The Jordan Quiver (one vertex, one loop).**
    A representation consists of a space $V$ and an operator $T: V \to V$. Since it is type $\tilde{A}_0$ (tame), it can be completely classified via the Jordan Canonical Form.

**9. [Dimension Vector] Define the dimension vector of a representation.**

??? success "Solution"
    **Definition:**
    An integer vector $\mathbf{d} = (\dim V_1, \dim V_2, \ldots, \dim V_n)$ that records the local dimensions at each node of the network.

**10. [Physics] Why do quivers appear in Superstring Theory?**

??? success "Solution"
    When studying D-branes at singularities, physicists found that the gauge symmetries and interaction terms between branes can be encoded as a quiver. Vertices represent brane types and arrows represent fermion fields. This mapping (Quiver Gauge Theory) transforms deep geometric data into intuitive linear algebraic networks.

## Chapter Summary

Quiver representation theory is the algebraic symphony of graph theory and linear operators:

1.  **Topologizing Maps**: It proves that the complexity of linear systems is defined by the "connectivity pattern" of mapping chains, binding operator properties to discrete graph theory.
2.  **Boundaries of Classification**: Gabriel's theorem establishes the ultimate boundary of "fully resolvable" linear systems (Dynkin graphs), revealing deep combinatorial symmetries in algebraic structures.
3.  **Modularity of Structure**: Through the theory of indecomposable representations, the quiver framework disassembles complex mapping networks into fundamental, irreducible components, providing powerful tools for modern algebraic geometry and physical system modeling.
