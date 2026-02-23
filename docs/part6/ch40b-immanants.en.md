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

****
??? success "Solution"
         $d_{(2,1)}(A) = 2 a_{11}a_{22}a_{33} - (a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32})$.

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    
eq NP$.

****
??? success "Solution"
    ## Chapter Summary

This chapter explores matrix invariants through the lens of group theory:


****: Unified determinants and permanents as specific instances of immanants.

****: Established the hierarchy of matrix functions for positive definite operators (Schur and Lieb theorems).

****: Linked immanants to Latin squares and planar matching theory.

****: Introduced the VP vs VNP problem, highlighting the computational gap between different immanants.
