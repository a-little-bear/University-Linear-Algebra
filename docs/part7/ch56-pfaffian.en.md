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
    For a $2n \times 2n$ skew-symmetric matrix $A$, the Pfaffian is defined as:
    $$\operatorname{pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    For odd dimensions, the Pfaffian is defined to be zero.

!!! theorem "Theorem 56.1 (Pfaffian-Determinant Identity)"
    For any skew-symmetric matrix $A$:
    $$\det(A) = [\operatorname{pf}(A)]^2$$

---

## Exercises

1. **[Calculation] Compute the Pfaffian of $A = \begin{pmatrix} 0 & a \\ -a & 0 \end{pmatrix}$.**
   ??? success "Solution"
       $\operatorname{pf}(A) = a$. Note that $\det(A) = a^2$, which is consistent with the square root identity.

2. **[Expansion] Write the explicit expansion of the Pfaffian for a $4 \times 4$ skew-symmetric matrix.**
   ??? success "Solution"
       $\operatorname{pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$. This corresponds to the three possible ways to pair 4 vertices into perfect matchings.

3. **[Graph Theory] Explain how the Pfaffian relates to perfect matchings in a graph.**
   ??? success "Solution"
       If $A$ is the adjacency matrix of a graph with a Pfaffian orientation, then $\operatorname{pf}(A)$ counts the number of perfect matchings. This provides a bridge between linear algebra and combinatorial counting.

4. **[Invariance] How does $\operatorname{pf}(P A P^T)$ behave for a permutation matrix $P$?**
   ??? success "Solution"
       $\operatorname{pf}(P A P^T) = \det(P) \operatorname{pf}(A)$. This means the Pfaffian is invariant under vertex relabeling up to the sign of the permutation.

5. **[Exterior Algebra] Express the Pfaffian using the wedge product of a bivector.**
   ??? success "Solution"
       Let $\omega = \sum_{i<j} a_{ij} e_i \wedge e_j$. The $n$-th exterior power satisfies $\frac{1}{n!} \omega^n = \operatorname{pf}(A) e_1 \wedge \dots \wedge e_{2n}$. This gives a geometric foundation to the Pfaffian.

6. **[Odd Dimension] Prove that the determinant of an $n \times n$ skew-symmetric matrix is zero if $n$ is odd.**
   ??? success "Solution"
       $\det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)$. For odd $n$, this implies $\det(A) = -\det(A)$, hence $\det(A) = 0$.

7. **[Eigenvalues] Describe the spectrum of a real skew-symmetric matrix.**
   ??? success "Solution"
       The eigenvalues are purely imaginary and appear in conjugate pairs $\pm i\theta_j$. Since $\det A = \prod (i\theta_j)(-i\theta_j) = \prod \theta_j^2$, the determinant is always non-negative.

8. **[FKT Algorithm] What is the core idea of the FKT algorithm?**
   ??? success "Solution"
       The Fisher-Kasteleyn-Temperley algorithm counts perfect matchings in planar graphs by constructing a specific orientation (Pfaffian orientation) such that the Pfaffian of the oriented adjacency matrix correctly sums the matchings without cancellations.

9. **[Canonical Form] State the canonical form of a skew-symmetric matrix under congruence transformations.**
   ??? success "Solution"
       For every skew-symmetric $A$, there exists a non-singular $P$ such that $P A P^T$ is block diagonal with $2 \times 2$ blocks $\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ and a zero block.

10. **[Physics] Where does the Pfaffian appear in statistical mechanics?**
    ??? success "Solution"
        In the Ising model, the partition function on a planar lattice can be expressed as a Pfaffian. It also characterizes the ground state of topological superconductors (e.g., Majorana fermions).

## Chapter Summary

This chapter explores the unique polynomial structure of antisymmetric matrices:

1. **Polynomial Square Root**: Defined the Pfaffian as the canonical "square root" of the determinant for skew-symmetric matrices.
2. **Combinatorial Matching**: Established the bijection between Pfaffian terms and perfect matchings in graphs.
3. **Topological Algorithms**: Detailed the FKT algorithm for counting matchings in planar graphs using linear algebra.
4. **Invariant Theory**: Formulated the Pfaffian in terms of exterior powers and congruence transformations.
