# Chapter 54: Quiver Representations

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch4) · Linear Transformations (Ch5) · Simultaneous Triangularization (Ch63A) · Graph Theory (Ch27)

**Chapter Outline**: Definition of Quivers (Directed Graphs) → Path Algebra → Quiver Representations as Vector Spaces and Maps → Morphisms between Representations → Indecomposable Representations → Gabriel's Theorem → Relation to Dynkin Diagrams → Root Systems → Application in Representation Theory

**Extension**: Quiver representation theory provides a geometric and diagrammatic way to study systems of linear maps, unifying diverse problems in linear algebra under a single visual framework.

</div>

A **quiver** is simply a directed graph. A **representation** of a quiver assigns a vector space to each vertex and a linear map to each arrow. This setup generalizes the study of a single linear operator (a quiver with one vertex and one loop) to entire configurations of interdependent spaces and maps. Gabriel's Theorem spectacularly links the complexity of these representations to the classical Dynkin diagrams of Lie theory.

---

## 54.1 Quivers and Their Representations

!!! definition "Definition 54.1 (Quiver Representation)"
    A representation $V$ of a quiver $Q = (V_0, E)$ over a field $K$ is a collection:
    - $\{V_i\}_{i \in V_0}$ of $K$-vector spaces;
    - $\{\phi_\alpha: V_i 	o V_j\}_{\alpha: i 	o j \in E}$ of linear maps.

!!! theorem "Theorem 54.1 (Gabriel's Theorem, 1972)"
    A quiver $Q$ has only finitely many indecomposable representations (up to isomorphism) if and only if its underlying undirected graph is a **Dynkin diagram** of type $A_n, D_n, E_6, E_7, E_8$.

---

## Exercises

1. **[Single Loop] Describe the representations of the quiver with one vertex and one loop.**
   ??? success "Solution"
       This corresponds to a vector space $V$ and an endomorphism $T: V 	o V$. Classifying these representations is equivalent to the Jordan Canonical Form problem. The indecomposable representations are exactly the Jordan blocks.

2. **[Path Algebra] What is the path algebra $KQ$ of a quiver?**
   ??? success "Solution"
       $KQ$ is the $K$-vector space with a basis consisting of all paths in $Q$ (including paths of length 0). Multiplication is defined by concatenation of paths. A quiver representation is equivalent to a module over its path algebra.

3. **[Dimension Vector] Define the dimension vector of a representation.**
   ??? success "Solution"
       The dimension vector $\mathbf{d} = (\dim V_1, \dots, \dim V_n)$ tracks the size of the vector space at each vertex. For finite type quivers, the dimension vectors of indecomposable representations are precisely the positive roots of the associated root system.

4. **[Morphism] When are two representations $V$ and $W$ isomorphic?**
   ??? success "Solution"
       $V \cong W$ if there exists a set of invertible linear maps $f_i: V_i 	o W_i$ such that for every arrow $\alpha: i 	o j$, the diagram commutes: $f_j \circ \phi_\alpha = \psi_\alpha \circ f_i$.

5. **[A2 Quiver] Find all indecomposable representations of the quiver $1 	o 2$ (Type $A_2$).**
   ??? success "Solution"
       There are three: $(K 	o 0)$, $(0 	o K)$, and $(K \xrightarrow{	ext{id}} K)$. These correspond to the roots $(1,0), (0,1), (1,1)$.

6. **[Dynkin Diagrams] Draw the Dynkin diagram $A_3$ and write its corresponding quiver.**
   ??? success "Solution"
       $A_3$ is $1-2-3$. A corresponding quiver could be $1 	o 2 	o 3$ or $1 \leftarrow 2 	o 3$. Gabriel's theorem implies the number of indecomposable representations is independent of the arrow orientation.

7. **[Kronecker Quiver] Analyze the Kronecker quiver (two vertices, two parallel arrows $1 ightrightarrows 2$).**
   ??? success "Solution"
       This is not a Dynkin diagram (it is type $	ilde{A}_1$ or Euclidean). It has infinitely many indecomposable representations, corresponding to the theory of matrix pencils $(A, B)$ and the Kronecker canonical form.

8. **[Direct Sum] Define the direct sum of two quiver representations.**
   ??? success "Solution"
       $(V \oplus W)_i = V_i \oplus W_i$ at each vertex, and for each arrow, the map is the block diagonal matrix $\begin{pmatrix} \phi_\alpha & 0 \ 0 & \psi_\alpha \end{pmatrix}$. This allows decomposing complex systems into simpler indecomposable units.

9. **[Subrepresentations] What is a subrepresentation of $V$?**
   ??? success "Solution"
       A subrepresentation $U \subseteq V$ is a collection of subspaces $U_i \subseteq V_i$ such that for every arrow $\alpha: i 	o j$, the map $\phi_\alpha$ maps $U_i$ into $U_j$.

10. **[Reflection Functors] Briefly describe the role of Bernstein-Gelfand-Ponomarev (BGP) reflection functors.**
    ??? success "Solution"
        Reflection functors provide a way to map representations of a quiver to representations of a quiver with reversed arrows at a vertex. This allows for an inductive proof of Gabriel's theorem by traversing the root system.

## Chapter Summary

This chapter unifies linear algebra through graph-based representations:

1. **Diagrammatic Synthesis**: Formulated configurations of spaces and maps as quiver representations, extending operator theory to network structures.
2. **Gabriel's Symmetry**: Discovered the profound link between finite representation types and Lie-theoretic Dynkin diagrams.
3. **Path Calculus**: Linked representations to path algebras, providing an associative algebra framework for graph-based maps.
4. **Root Correspondence**: Established the bijection between indecomposable representations and the geometry of root systems.
