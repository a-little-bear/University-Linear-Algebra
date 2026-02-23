# Chapter 38B: P-Matrices, H-Matrices and Related Classes

<div class="context-flow" markdown>

**Prerequisites**: Non-negative Matrices (Ch17) · M-Matrices (Ch38A) · Eigenvalues (Ch6) · LCP Basics

**Chapter Outline**: P-Matrices & Principal Minors → LCP Unique Solvability → P0-Matrices → Comparison Matrices & H-Matrices → Diagonal Dominance Hierarchy → N-Matrices → Semipositive Matrices → Ostrowski-Reich Theorem → Iterative Stability

**Extension**: P-matrices are central to mathematical programming (LCP theory); H-matrices generalize M-matrix theory to complex matrices and are vital for the convergence analysis of iterative methods.

</div>

This chapter explores broader matrix classes than Z-matrices. **P-matrices** require all principal minors to be positive, with no restriction on sign patterns. **H-matrices** utilize the comparison matrix to extend M-matrix properties to general complex matrices. These classes provide the analytical foundation for Linear Complementarity Problems (LCP) and numerical stability.

---

## 38B.1 P-Matrices

!!! definition "Definition 38B.1 (P-Matrix)"
    A matrix $A \in M_n(\mathbb{R})$ is a **P-matrix** if all its principal minors are positive:
    $$\det(A[\alpha, \alpha]) > 0, \quad \forall\, \alpha \subseteq \{1, \dots, n\}, \alpha 
e \emptyset.$$

!!! theorem "Theorem 38B.3 (LCP Solvability)"
    $A$ is a P-matrix if and only if the Linear Complementarity Problem $\operatorname{LCP}(q, A)$ has a unique solution for every $q \in \mathbb{R}^n$.

---

## Exercises

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
         $\mathcal{M}(A) = \begin{pmatrix} 4 & -\sqrt{2} \ -\sqrt{5} & 5 \end{pmatrix}$.

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    
eq 0$) such that $Ax > 0$. Every non-singular M-matrix is semipositive.

****
??? success "Solution"
    ## Chapter Summary

This chapter classifies matrices based on the positivity of their sub-structures:


****: Defined P-matrices through principal minors, establishing their role in LCP solvability.

****: Developed H-matrices as the complex generalization of M-matrices via magnitude-based dominance.

****: Linked these classes to the numerical stability of iterative solvers and error bounds in interval arithmetic.
