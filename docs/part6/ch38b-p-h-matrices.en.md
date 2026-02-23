# Chapter 38B: P-matrices and H-matrices

<div class="context-flow" markdown>

**Prerequisites**: M-matrices (Ch38A) · Positive Definite Matrices (Ch16) · Determinants (Ch03)

**Chapter Outline**: Generalization of M-matrices → Definition of P-matrices (All Principal Minors Positive) → Solving Linear Complementarity Problems (LCP) → H-matrices (Comparison Matrix is an M-matrix) → Ostrowski Criterion → Hierarchical Relationships: Positive Definite $\subset$ M-matrices $\subset$ H-matrices → Applications in Stability and Numerical Computation: Circuit Analysis and Linearization of Non-linear Systems

**Extension**: P-matrices and H-matrices are the advanced algebraic language for "non-symmetric stability"; they preserve many key properties of positive definite matrices without requiring symmetry, serving as the core for solving complementarity problems in operations research.

</div>

If we relax the "non-positive off-diagonals" restriction of M-matrices, we obtain the more general **H-matrices**; if we only retain the "all principal minors positive" characteristic, we get **P-matrices**. These two classes of matrices are critical in non-linear programming and circuit simulation, as they guarantee the existence and uniqueness of equilibrium points in complex systems. This chapter reveals the hierarchical connections between these structures.

---

## 38B.1 P-matrices and the LCP

!!! definition "Definition 38B.1 (P-matrix)"
    A square matrix $A$ is a **P-matrix** if every principal minor (not just leading ones) is strictly positive.
    **Properties**: Symmetric positive definite matrices and non-singular M-matrices are always P-matrices.

!!! technique "Application: Linear Complementarity Problem (LCP)"
    P-matrices provide the necessary and sufficient condition for the LCP ($w - Az = q, w \ge 0, z \ge 0, w^T z = 0$) to have a unique solution for any $q$. This is invaluable in game theory and contact mechanics.

---

## 38B.2 H-matrices and Comparison Matrices

!!! definition "Definition 38B.2 (H-matrix)"
    A square matrix $A$ is an **H-matrix** if its **Comparison Matrix** $\mathcal{M}(A)$ is a non-singular M-matrix.
    The comparison matrix is defined as:
    - Diagonal: $(\mathcal{M}(A))_{ii} = |a_{ii}|$
    - Off-diagonal: $(\mathcal{M}(A))_{ij} = -|a_{ij}|$ ($i \neq j$)

---

## 38B.3 Criteria and Hierarchies

!!! theorem "Theorem 38B.1 (Ostrowski Criterion)"
    If $A$ is generalized strictly diagonally dominant (i.e., there exists a positive vector $d$ such that $|a_{ii}| d_i > \sum_{j \neq i} |a_{ij}| d_j$), then $A$ is a non-singular H-matrix.

---

## Exercises

**1. [Basics] Determine if $A = \begin{pmatrix} 1 & -2 \\ 0 & 1 \end{pmatrix}$ is a P-matrix.**

??? success "Solution"
    **Calculate Principal Minors:**
    1. $1 \times 1$ minors: $|1|=1, |1|=1$ (both positive).
    2. $2 \times 2$ minor: $\det(A) = 1 - 0 = 1 > 0$.
    **Conclusion**: Since all principal minors are positive, it is a P-matrix. Note: It is neither symmetric nor positive definite.

**2. [Comparison] Find the comparison matrix for $A = \begin{pmatrix} 2 & i \\ -1 & 3 \end{pmatrix}$.**

??? success "Solution"
    **Construction:**
    1. Absolute values of diagonals: $2, 3$.
    2. Negative absolute values of off-diagonals: $-|i| = -1, -|-1| = -1$.
    **Result**: $\mathcal{M}(A) = \begin{pmatrix} 2 & -1 \\ -1 & 3 \end{pmatrix}$.

**3. [H-matrix Check] Is the matrix $A$ from the previous problem an H-matrix?**

??? success "Solution"
    **Determination:**
    1. Check if $\mathcal{M}(A)$ is an M-matrix.
    2. It is a Z-matrix, and its principal minors are $2 > 0$ and $6-1=5 > 0$.
    **Conclusion**: Since the comparison matrix is a non-singular M-matrix, $A$ is a non-singular H-matrix.

**4. [Properties] Can the eigenvalues of a P-matrix have negative real parts?**

??? success "Solution"
    **Conclusion: Yes.**
    **Reasoning**: P-matrices only guarantee positive minors, not necessarily positive real parts for eigenvalues (which is a stronger property held by M-matrices or PD matrices).
    Example: $A = \begin{pmatrix} 1 & 2 \\ -2 & 1 \end{pmatrix}$ is a P-matrix (minors 1, 1, 5), but its eigenvalues are $1 \pm 2i$. While these have positive real parts, in higher dimensions, more complex distributions are possible.

**5. [Hierarchy] Briefly state the containment relationship between PD, M, and P matrices.**

??? success "Solution"
    **Inclusion Chain:**
    (Symmetric PD) $\cup$ (Non-singular M-matrices) $\subset$ (P-matrices).
    P-matrices represent the broadest generalization of "determinantal positivity."

**6. [Inversion] Is the inverse of an H-matrix always non-negative?**

??? success "Solution"
    **Conclusion: Not necessarily.**
    **Analysis**: Only M-matrices guarantee $A^{-1} \ge 0$. For a general H-matrix (e.g., with complex entries), the inverse is usually complex. However, H-matrices satisfy $\|A^{-1}\| \le \|\mathcal{M}(A)^{-1}\|$, which is useful for error bounds.

**7. [Calculation] Determine if $\begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ is a P-matrix.**

??? success "Solution"
    **Calculation:**
    Principal minors: 1, 1. But $\det = 1 - 4 = -3 < 0$.
    **Conclusion**: It is not a P-matrix.

**8. [Application] In circuit analysis, what type of matrix is the conductance matrix formed by resistors?**

??? success "Solution"
    Typically an **M-matrix** (due to diagonal dominance from Kirchhoff's Current Law and negative conductances on off-diagonals). If dependent sources are present, it may generalize to an **H-matrix**.

**9. [Stability] Why are H-matrices important in neural network stability?**

??? success "Solution"
    Because non-linear activation functions are often bounded in slope. Using the Ostrowski criterion for H-matrices allows one to determine if the weight matrix guarantees a unique global equilibrium.

**10. [Limit] As the degree of diagonal dominance increases, what does an H-matrix approach?**

??? success "Solution"
    It approaches a **diagonal matrix**. Diagonal dominance is the core numerical trait of H-matrices; the stronger the dominance, the more the inverse resembles that of a diagonal matrix, increasing system stability.

## Chapter Summary

P-matrices and H-matrices define the boundaries of generalized stability:

1.  **Dominance of Minors**: P-matrices prove that even without symmetry, the overall positivity of principal minors ensures unique solutions to linear complementarity problems—an algebraic cornerstone of operations research.
2.  **Abstraction of Magnitudes**: Through the comparison matrix technique, H-matrices simplify complex numerical (even complex-valued) operations into magnitude analysis using M-matrices, providing unified bounds for error propagation.
3.  **Evolution of Structure**: From Positive Definite to M to H, this hierarchy demonstrates how linear algebra generalizes stability logic to non-linear and non-symmetric domains by step-wise relaxation of constraints.
