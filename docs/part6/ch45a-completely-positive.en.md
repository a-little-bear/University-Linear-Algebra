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

****

??? success "Solution"
    
       1. For any $x$, $x^T A x = \|B^T x\|^2 \ge 0$, so $A \succeq 0$.
       2. $a_{ij} = \sum_l b_{il} b_{jl} \ge 0$ since $b_{il}, b_{jl} \ge 0$.
       Thus, $A$ is both PSD and non-negative (DNN).

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
    
eq \mathcal{DNN}_5$.

****

??? success "Solution"
    

****

??? success "Solution"
    

## Chapter Summary

This chapter explores the boundary between non-negativity and positive semi-definiteness:

1. **Structured Positivity**: Defined CP matrices as those admitting a non-negative factorization.
2. **Rank Complexity**: Distinguished the cp-rank from the standard rank, highlighting its role as a combinatorial invariant.
3. **Geometric Limits**: Identified the separation between CP and DNN cones as a phenomenon arising from odd cycles in the underlying graph.
4. **Quantum Extension**: Introduced CPSD matrices, linking classical matrix positivity to quantum correlation theory.
