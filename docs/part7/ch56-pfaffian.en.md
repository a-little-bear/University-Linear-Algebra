# Chapter 56: The Pfaffian and Antisymmetric Matrices

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Determinants (Ch3) · Exterior Algebra (Ch49) · Graph Theory (Ch27)

**Chapter Outline**: Definition of the Pfaffian → Relation to the Determinant $\det A = (\operatorname{pf} A)^2$ → Properties of Antisymmetric Matrices → Perfect Matchings in Graphs → FKT Algorithm for Planar Graphs → Pfaffian Orientation → Application in Statistical Mechanics (Ising Model)

**Extension**: The Pfaffian is the central tool for counting perfect matchings in planar graphs and describes the states of certain topological phases of matter.

</div>

For a skew-symmetric matrix $A$ (where $A^T = -A$), the determinant is always the square of a polynomial in its entries. This polynomial is called the **Pfaffian**, denoted $\operatorname{pf}(A)$. While the determinant involves all permutations, the Pfaffian only sums over partitions of the set into disjoint pairs, creating a deep link between linear algebra and the combinatorial theory of matchings.

---

## 56.1 Definition and Fundamental Relation

!!! definition "Definition 56.1 (The Pfaffian)"
    For a $2n 	imes 2n$ skew-symmetric matrix $A$, the Pfaffian is defined as:
    $$\operatorname{pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    For odd dimensions, the Pfaffian is defined to be zero.

!!! theorem "Theorem 56.1 (Pfaffian-Determinant Identity)"
    For any skew-symmetric matrix $A$:
    $$\det(A) = [\operatorname{pf}(A)]^2$$

---

## Exercises

1. **[Calculation] Compute the Pfaffian of $A = \begin{pmatrix} 0 & a \ -a & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $\operatorname{pf}(A) = a$. Note that $\det(A) = a^2$, consistent with the theorem.

2. **[Expansion] Compute the Pfaffian of a $4 	imes 4$ skew-symmetric matrix.**
   ??? success "Solution"
       $\operatorname{pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$. This corresponds to the three possible perfect matchings of 4 vertices.

3. **[Graph Theory] Relate the Pfaffian of the adjacency matrix of a graph to its matchings.**
   ??? success "Solution"
       If $A$ is the adjacency matrix of a graph with a Pfaffian orientation, then $\operatorname{pf}(A)$ counts the number of perfect matchings in the graph. This allows matching problems to be solved via linear algebra.

4. **[Properties] How does $\operatorname{pf}(P A P^T)$ relate to $\operatorname{pf}(A)$?**
   ??? success "Solution"
       $\operatorname{pf}(P A P^T) = \det(P) \operatorname{pf}(A)$. If $P$ is a permutation matrix, this shows the Pfaffian is invariant under relabeling of vertices up to the sign of the permutation.

5. **[Exterior Algebra] Express the Pfaffian using the wedge product.**
   ??? success "Solution"
       Let $\omega = \sum_{i<j} a_{ij} e_i \wedge e_j$. Then the $n$-th exterior power satisfies $\frac{1}{n!} \omega^n = \operatorname{pf}(A) e_1 \wedge \dots \wedge e_{2n}$.

6. **[Odd Dimension] Prove that the determinant of an $n 	imes n$ skew-symmetric matrix is zero if $n$ is odd.**
   ??? success "Solution"
       $\det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)$. For odd $n$, $\det(A) = -\det(A)$, so $\det(A) = 0$.

7. **[Eigenvalues] Describe the eigenvalues of a real skew-symmetric matrix.**
   ??? success "Solution"
       The eigenvalues are purely imaginary and come in conjugate pairs $\pm i	heta$. This is consistent with the fact that $\det A = \prod (i	heta_j)(-i	heta_j) = \prod 	heta_j^2 \ge 0$.

8. **[FKT Algorithm] What is the significance of the FKT algorithm?**
   ??? success "Solution"
       The Fisher-Kasteleyn-Temperley (FKT) algorithm computes the number of perfect matchings in a planar graph in polynomial time by finding a Pfaffian orientation and computing the Pfaffian of the resulting matrix.

9. **[Block Form] What is the canonical form of a skew-symmetric matrix under congruence?**
   ??? success "Solution"
       There exists a non-singular $P$ such that $P A P^T = \operatorname{diag}(J, J, \dots, J, 0, \dots, 0)$, where $J = \begin{pmatrix} 0 & 1 \ -1 & 0 \end{pmatrix}$.

10. **[Physics] In what context does the Pfaffian appear in the study of Fermions?**
    ??? success "Solution"
        The partition function of a free Fermion system can be expressed as a Pfaffian of the interaction matrix. It also appears in the description of the Pfaffian state in the fractional quantum Hall effect (a non-Abelian topological phase).

## Chapter Summary

This chapter explores the unique polynomial structure of antisymmetric matrices:

1. **Polynomial Square Root**: Defined the Pfaffian as the canonical "square root" of the determinant for skew-symmetric matrices.
2. **Combinatorial Matching**: Established the bijection between Pfaffian terms and perfect matchings in graphs.
3. **Topological Algorithms**: Detailed the FKT algorithm for counting matchings in planar graphs using linear algebra.
4. **Invariant Theory**: Formulated the Pfaffian in terms of exterior powers and congruence transformations.
