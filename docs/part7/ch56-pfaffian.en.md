# Chapter 56: The Pfaffian and Skew-symmetric Matrices

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch03) · Matrix Groups (Ch55) · Symplectic Matrices (Ch53) · Basic Combinatorics

**Chapter Outline**: Special Properties of Skew-symmetric Matrices → Definition of the Pfaffian → Relationship with the Determinant $\operatorname{Pf}(A)^2 = \det(A)$ → Algebraic Properties (Invariance under Congruence) → Explicit Formula via Perfect Matchings → Pfaffian version of Binet-Cauchy → Applications: Counting Perfect Matchings in Planar Graphs (FKT Algorithm), the Ising Model in Statistical Physics, and Majorana Fermions

**Extension**: The Pfaffian is the "square root" of the determinant for skew-symmetric matrices; it provides a more refined description than the determinant for physical systems with parity symmetry and combinatorial counting problems.

</div>

For symmetric matrices, we have theories of positive definiteness and eigenvalues. For **Skew-symmetric Matrices**, however, the most profound scalar function is not the determinant, but the **Pfaffian**. Because the determinant of a skew-symmetric matrix is always the square of a polynomial, the Pfaffian is exactly that polynomial. This chapter reveals how this "algebraic square root of the determinant" connects linear algebra, graph theory, and quantum physics.

---

## 56.1 Definition of the Pfaffian

!!! definition "Definition 56.1 (The Pfaffian)"
    Let $A$ be a $2n \times 2n$ skew-symmetric matrix. Its **Pfaffian** $\operatorname{Pf}(A)$ is defined as:
    $$\operatorname{Pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    For a skew-symmetric matrix of odd order, the Pfaffian is defined to be 0.

!!! example "Example 56.1"
    For a $2 \times 2$ matrix $A = \begin{pmatrix} 0 & a \\ -a & 0 \end{pmatrix}$, the Pfaffian is $\operatorname{Pf}(A) = a$.
    Note that $\det(A) = a^2 = \operatorname{Pf}(A)^2$.

---

## 56.2 Core Properties

!!! theorem "Theorem 56.1 (Relation to Determinant)"
    For any $2n \times 2n$ skew-symmetric matrix $A$:
    $$\det(A) = [\operatorname{Pf}(A)]^2$$

!!! theorem "Theorem 56.2 (Congruence Invariance)"
    For any $2n \times 2n$ matrix $M$:
    $$\operatorname{Pf}(M A M^T) = \det(M) \operatorname{Pf}(A)$$
    This implies that the Pfaffian is invariant under rotations (when $\det(M)=1$).

---

## 56.3 Combinatorial Significance: FKT Algorithm

!!! technique "Counting Perfect Matchings"
    Let $G$ be a planar graph and $A$ its associated Pfaffian orientation matrix. The number of perfect matchings in $G$ is exactly $|\operatorname{Pf}(A)|$. This is the **FKT Algorithm**, which allows for the counting of perfect matchings in planar graphs in polynomial time, whereas the general problem is #P-complete.

---

## Exercises


****
??? success "Solution"
     $\operatorname{Pf} = 3$.


****
??? success "Solution"
     $\pm 4$. The sign depends on the specific arrangement of the matrix entries.


****
??? success "Solution"
     The determinant of an odd-order skew-symmetric matrix is always 0. Thus, its Pfaffian (as the square root) must also be 0.


****
??? success "Solution"
     For a block-diagonal skew-symmetric matrix, $\operatorname{Pf}(A \oplus B) = \operatorname{Pf}(A)\operatorname{Pf}(B)$. Thus $\operatorname{Pf}(J) = 1 \cdot 1 = 1$.


****
??? success "Solution"
     $\operatorname{Pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$. This is obtained by enumerating all perfect matchings of 4 elements.


****
??? success "Solution"
     It corresponds to Majorana bound states with non-Abelian statistics, which are candidates for topological quantum computation.


****
??? success "Solution"
     They appear in conjugate pairs of purely imaginary numbers: $\pm i\lambda_1, \pm i\lambda_2, \ldots$. The Pfaffian is the product of these imaginary parts: $\prod \lambda_j$.


****
??? success "Solution"
     From the determinant property $\det(kA) = k^{2n} \det(A)$, taking the square root yields the result.


****
??? success "Solution"
     Because general graphs do not possess a uniform Pfaffian orientation (it is impossible to assign signs to edges such that all cycles contribute consistently), reflecting the algebraic uniqueness of planarity.

****
??? success "Solution"
    ## Chapter Summary

The Pfaffian is an exquisite structure within skew-symmetric algebra:


****: It provides a perfect algebraic explanation for why the determinant of a skew-symmetric matrix must be a non-negative square, establishing the fundamental scalar under this specific symmetry.

****: Via the FKT algorithm, the Pfaffian pulls counting problems that are inherently exponential back into the safety of polynomial complexity, demonstrating algebra's immense power for dimensionality reduction in combinatorics.

****: From classical thermodynamics to cutting-edge topological superconductivity, the Pfaffian serves as the natural language for describing pairing phenomena, proving the depth of matrix theory in revealing the microscopic order of matter.
