# Chapter 40B: Immanants and Generalized Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Determinants (Ch3) · Permanents (Ch40A) · Group Representation Theory

**Chapter Outline**: $S_n$ Representation Theory (Young Diagrams, Irreducible Characters) → Immanant Definition → Merris Generalized Matrix Functions → Schur Inequalities → Lieb's Theorem & Permanent Dominance Conjecture → Stembridge Conjecture → Relation to Schur Functions → Latin Squares → Kasteleyn's Theorem & FKT Algorithm → VP vs. VNP

**Extension**: Immanants unify determinants and permanents through the representation theory of the symmetric group; they are central to algebraic combinatorics.

</div>

Determinants and permanents correspond to the two "extreme" characters of the symmetric group $S_n$: the sign character $\operatorname{sgn}$ and the trivial character $\mathbf{1}$. **Immanants** provide a unified framework by allowing the use of *any* irreducible character $\chi^\lambda$ of $S_n$. This chapter formalizes these functions and explores their profound connections to symmetric functions and computational complexity.

---

## 40B.1 Representation Theory of $S_n$

!!! definition "Definition 40B.1 (Partitions and Young Diagrams)"
    A partition $\lambda \vdash n$ is a non-increasing sequence $(\lambda_1, \dots, \lambda_k)$ of positive integers summing to $n$. Each partition corresponds to a unique irreducible character $\chi^\lambda$ of $S_n$.

!!! theorem "Theorem 40B.1 (Hook Length Formula)"
    The dimension of the irreducible representation $S^\lambda$ (denoted $f^\lambda$) is given by $f^\lambda = n! / \prod h(i,j)$, where $h(i,j)$ are the hook lengths of the Young diagram.

---

## 40B.2 Immanants

!!! definition "Definition 40B.6 (Immanant)"
    For an $n 	imes n$ matrix $A$ and a character $\chi$ of $S_n$, the $\chi$-immanant is:
    $$d_\chi(A) = \sum_{\sigma \in S_n} \chi(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}$$
    - If $\chi = \operatorname{sgn}$, $d_\chi(A) = \det(A)$.
    - If $\chi = \mathbf{1}$, $d_\chi(A) = \operatorname{perm}(A)$.

---

## Exercises

1. **[Fundamentals] Compute the $d_{(2,1)}$ immanant for a $3 	imes 3$ matrix $A$ using its character table.**

   ??? success "Solution"
       For $S_3$, $\chi^{(2,1)}$ values are: $\chi(	ext{id})=2, \chi(	ext{transposition})=0, \chi(	ext{3-cycle})=-1$.
       $d_{(2,1)}(A) = 2 a_{11}a_{22}a_{33} - (a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32})$.

2. **[Diagonal] Find $d_\lambda(D)$ for $D = \operatorname{diag}(x_1, \dots, x_n)$.**

   ??? success "Solution"
       $d_\lambda(D) = \chi^\lambda(	ext{id}) \prod x_i = f^\lambda \prod x_i$. Only the identity permutation contributes to the sum.

3. **[Schur Inequality] Prove that for $A \succeq 0$, $d_\lambda(A) \ge f^\lambda \det(A)$.**

   ??? success "Solution"
       This is Schur's theorem (1918). It follows from the fact that the immanant corresponds to the trace of a positive operator on the symmetry-adapted subspace $S^\lambda$.

4. **[Permanent Dominance] State Lieb's Theorem regarding immanants.**

   ??? success "Solution"
       For $A \succeq 0$, $|d_\chi(A)| \le \chi(	ext{id}) \operatorname{perm}(A)$. This proves the permanent is the "largest" among all immanants for positive semi-definite matrices.

5. **[Stembridge] What is the Stembridge Conjecture?**

   ??? success "Solution"
       It states that the immanants of certain Jacobi-Trudi matrices, when expressed as symmetric functions, are Schur-positive (non-negative combinations of Schur functions).

6. **[Latin Squares] How do permanents count the number of Latin squares?**

   ??? success "Solution"
       The number of ways to add a row to a $r 	imes n$ Latin rectangle is given by the permanent of a specific $(0,1)$-matrix representing allowed entries.

7. **[Kasteleyn] Why can the number of perfect matchings in a planar graph be computed in $O(n^3)$?**

   ??? success "Solution"
       Kasteleyn's theorem shows that for planar graphs, we can orient edges such that the number of matchings equals the Pfaffian of the adjacency matrix, which is the square root of the determinant.

8. **[Complexity] Define the classes VP and VNP.**

   ??? success "Solution"
       VP is the class of polynomials computable by polynomial-sized circuits (like the determinant). VNP is the class of polynomials where each coefficient is "easy" to compute (like the permanent).

9. **[Valiant] Is computing the permanent VNP-complete?**

   ??? success "Solution"
       Yes. Valiant (1979) proved that the permanent is the hardest polynomial in VNP. Proving it is not in VP is equivalent to the algebraic version of $P 
eq NP$.

10. **[Symmetry] Prove $d_\chi(A^T) = d_\chi(A)$ for real characters.**

   ??? success "Solution"
        $d_\chi(A^T) = \sum \chi(\sigma) \prod a_{\sigma(i),i} = \sum \chi(\sigma) \prod a_{i,\sigma^{-1}(i)}$. Since $\chi(\sigma) = \chi(\sigma^{-1})$ for $S_n$, the result follows.

## Chapter Summary

This chapter explores matrix invariants through the lens of group theory:

1. **Representation Framework**: Unified determinants and permanents as specific instances of immanants.
2. **Spectral Bounds**: Established the hierarchy of matrix functions for positive definite operators (Schur and Lieb theorems).
3. **Combinatorial Counting**: Linked immanants to Latin squares and planar matching theory.
4. **Complexity Barrier**: Introduced the VP vs VNP problem, highlighting the computational gap between different immanants.
