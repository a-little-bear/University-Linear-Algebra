# Chapter 45A: Completely Positive Matrices

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Non-negative Matrices (Ch17) · Convex Cones (Ch25)

**Chapter Outline**: Definition of Completely Positive (CP) Matrices → cp-rank and Upper Bounds → Extreme Rays of the CP Cone → Doubly Non-negative (DNN) Matrices → CP=DNN for Small Dimensions → Graph-Theoretic Characterization → CP Semidefinite Matrices (Quantum Generalization)

**Extension**: The relation between CP matrices and non-negative rank is a core problem in communication complexity; CP semidefinite matrices are basic tools in the study of Bell inequalities.

</div>

A symmetric matrix $A$ is **completely positive** if it can be decomposed as $A = BB^T$, where $B$ is a non-negative matrix. This condition combines positive semi-definiteness with element-wise non-negativity, linking algebraic properties to combinatorial structures. The study of the CP cone $\mathcal{CP}_n$ is central to convex optimization and graph theory.

---

## 45A.1 Core Concepts

!!! definition "Definition 45A.1 (Completely Positive Matrix)"
    A symmetric matrix $A \in \mathbb{R}^{n 	imes n}$ is **completely positive** (CP) if there exists a matrix $B \in \mathbb{R}^{n 	imes k}_{\ge 0}$ such that $A = BB^T$.

!!! definition "Definition 45A.2 (cp-rank)"
    The **cp-rank** of $A \in \mathcal{CP}_n$ is the smallest number of columns $k$ required for a non-negative decomposition $A = BB^T$.

---

## Exercises

1. **[Fundamentals] Prove that every CP matrix is doubly non-negative (DNN).**
   ??? success "Solution"
       If $A = BB^T$ with $B \ge 0$, then:
       1. For any $x$, $x^T A x = \|B^T x\|^2 \ge 0$, so $A \succeq 0$.
       2. $a_{ij} = \sum_l b_{il} b_{jl} \ge 0$ since $b_{il}, b_{jl} \ge 0$.
       Thus, $A$ is both PSD and non-negative (DNN).

2. **[cp-rank] Find the cp-rank of the identity matrix $I_n$.**
   ??? success "Solution"
       $I_n = \sum_{i=1}^n e_i e_i^T$. Since the standard basis vectors $e_i$ are non-negative and the decomposition uses $n$ terms, and $\operatorname{rank}(I_n)=n$, the cp-rank is $n$.

3. **[Dimensions] For which dimensions $n$ does $\mathcal{CP}_n = \mathcal{DNN}_n$ hold?**
   ??? success "Solution"
       The equality holds if and only if $n \le 4$ (Maxfield-Minc Theorem). For $n \ge 5$, there exist DNN matrices that are not CP.

4. **[Caratheodory] State the Caratheodory upper bound for the cp-rank.**
   ??? success "Solution"
       $\operatorname{cp-rank}(A) \le \frac{n(n+1)}{2}$. This is the dimension of the space of symmetric matrices.

5. **[DJL Conjecture] What is the DJL conjecture regarding the cp-rank?**
   ??? success "Solution"
       It conjectures that for $n \ge 4$, $\operatorname{cp-rank}(A) \le \lfloor n^2/4 floor$.

6. **[Extreme Rays] What are the extreme rays of the CP cone?**
   ??? success "Solution"
       The extreme rays are the rank-1 matrices $bb^T$ where $b \ge 0$.

7. **[Graph Theory] Relate the comparison graph $G(A)$ to the CP property.**
   ??? success "Solution"
       According to the Berman-Xu theorem, a DNN matrix is CP if its comparison graph contains no odd holes (odd cycles of length $\ge 5$).

8. **[Horn Matrix] Why is the Horn matrix significant in CP theory?**
   ??? success "Solution"
       The Horn matrix is a $5 	imes 5$ matrix that is copositive but not in the dual of the CP cone, effectively proving that $\mathcal{CP}_5 
eq \mathcal{DNN}_5$.

9. **[Hadamard Product] Is the Hadamard product of two CP matrices always CP?**
   ??? success "Solution"
       Yes. If $A = \sum a_i a_i^T$ and $B = \sum b_j b_j^T$, then $A \circ B = \sum_{i,j} (a_i \circ b_j)(a_i \circ b_j)^T$, which is a sum of non-negative rank-1 matrices.

10. **[Submatrices] Prove that every principal submatrix of a CP matrix is also CP.**
    ??? success "Solution"
        If $A = BB^T$, a principal submatrix $A_S$ is formed by $B_S B_S^T$ where $B_S$ contains the rows of $B$ indexed by $S$. Since $B \ge 0$ implies $B_S \ge 0$, $A_S$ is CP.

## Chapter Summary

This chapter explores the boundary between non-negativity and positive semi-definiteness:

1. **Structured Positivity**: Defined CP matrices as those admitting a non-negative factorization.
2. **Rank Complexity**: Distinguished the cp-rank from the standard rank, highlighting its role as a combinatorial invariant.
3. **Geometric Limits**: Identified the separation between CP and DNN cones as a phenomenon arising from odd cycles in the underlying graph.
4. **Quantum Extension**: Introduced CPSD matrices, linking classical matrix positivity to quantum correlation theory.
