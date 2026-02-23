# Chapter 60: Linear Preserver Problems (LPP)

<div class="context-flow" markdown>

**Prerequisites**: Linear Transformations (Ch05) · Determinants (Ch03) · Matrix Groups (Ch55) · Basics of Operator Theory

**Chapter Outline**: Motivation for Linear Preserver Problems → Preserving the Determinant (Frobenius's Theorem) → Preserving Rank (Dieudonné's Theorem) → Preserving the Spectrum and Eigenvalues → Preserving Invertibility and Singularity → Preserving Norms (Isometries) → Preserving Positive Definiteness and Commutativity → Applications: Operator Evolution in Quantum Information and Symmetry Analysis in Matrix Theory

**Extension**: Linear Preserver Problems investigate "which transformations leave the soul (essential properties) of a matrix intact"; they are the ultimate tool for understanding symmetries in matrix spaces, revealing the deepest rigid structures in operator algebras.

</div>

In studying linear operators $L: M_n \to M_n$ on matrix spaces, a natural question arises: if $L$ preserves a specific matrix property (such as the determinant, rank, or eigenvalues), what form must $L$ take? This is known as the **Linear Preserver Problem** (LPP). Such studies not only reveal the algebraic rigidity of matrix spaces but also provide the theoretical foundation for symmetry transformations in quantum mechanics and functional analysis.

---

## 60.1 Preserving the Determinant: Frobenius's Theorem

!!! theorem "Theorem 60.1 (Frobenius's Theorem)"
    Let $L: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ be a linear transformation. If $\det(L(A)) = \det(A)$ for all matrices $A$, then there exist non-singular matrices $M, N$ with $\det(MN) = 1$ such that $L$ has one of the following two forms:
    1.  **Equivalent Form**: $L(A) = MAN$
    2.  **Transpose Equivalent Form**: $L(A) = MA^T N$
    **Significance**: This shows that the determinant is not just a value; it strictly constrains the possible forms of transformations on the matrix space.

---

## 60.2 Preserving Rank and Singularity

!!! theorem "Theorem 60.2 (Dieudonné's Theorem)"
    A linear transformation $L$ preserves the set of all rank-1 matrices (or preserves the set of singular matrices) if and only if it has the form $L(A) = MAN$ or $L(A) = MA^T N$ for non-singular $M, N$.
    **Insight**: Rank is the most fundamental topological feature of a matrix; transformations that preserve it are almost inevitably a reorganization of the coordinate system.

---

## 60.3 Preserving the Spectrum and Eigenvalues

!!! theorem "Theorem 60.3 (Preserving the Spectrum)"
    If $L$ preserves the spectrum of a matrix (all eigenvalues and their multiplicities), then $L$ must be a **similarity transformation** or its transpose:
    $$L(A) = P^{-1} A P \quad \text{or} \quad L(A) = P^{-1} A^T P$$
    This proves that similarity transformations are the only "legitimate" symmetric operations in eigenvalue theory.

---

## 60.4 Preserving Positive Definiteness and Norms

!!! note "Other Preserved Properties"
    - **Isometries**: Transformations preserving the Frobenius norm must be of the form $L(A) = UAV$ where $U, V$ are unitary.
    - **Positive Definiteness**: Transformations preserving the set of symmetric positive definite matrices must have the form $L(A) = P A P^T$ for some non-singular $P$.

---

## Exercises


****
??? success "Solution"
     $\det(P^{-1} A P) = \det(P^{-1})\det(A)\det(P) = \det(A)$.


****
??? success "Solution"
     $L(A) = MAN$ where $M, N$ are not inverses of each other. For example, $M=2I, N=0.5I$. While $\det(MN)=1$ ensures the determinant is preserved, the individual eigenvalues will change because the transformation is not a similarity transform.


****
??? success "Solution"
     Not necessarily. It might scale the determinant by a constant factor.


****
??? success "Solution"
     The form is very broad and not limited to $MAN$. Any linear transformation satisfying $\operatorname{tr}(L(E_{ij})) = \delta_{ij}$ works.


****
??? success "Solution"
     If a linear transformation maps the set of non-singular matrices to itself, by dimensionality and topological continuity, it must also map the boundary (singular matrices) to itself.


****
??? success "Solution"
     Yes. Since $\det(A-\lambda I) = \det((A-\lambda I)^T) = \det(A^T - \lambda I)$, the eigenvalues are identical.


****
??? success "Solution"
     $\|UAV\|_F^2 = \operatorname{tr}((UAV)^*(UAV)) = \operatorname{tr}(V^* A^* U^* U A V) = \operatorname{tr}(V^* A^* A V) = \operatorname{tr}(A^* A) = \|A\|_F^2$.


****
??? success "Solution"
     Because a physical quantum channel must map a valid density matrix (positive definite with trace 1) to another valid density matrix.


****
??? success "Solution"
     Usually a combination of scalar multiplication and a similarity transformation.

****
??? success "Solution"
    ## Chapter Summary

Linear preserver problems establish the algebraic boundaries of matrix properties:


****: Theorems by Frobenius and Dieudonné prove that the core attributes of a matrix (determinant, rank) are extremely "solid," and any linear attempt to preserve them eventually returns to fundamental coordinate rotations and transpositions.

****: From preserving the spectrum (strong constraint) to preserving the trace (weak constraint), LPP demonstrates how different mathematical features exert varying degrees of compression on the transformation space.

****: By investigating preservation properties, we can redefine what constitutes "physically consistent" evolution from a purely algebraic perspective, providing rigorous criteria for modern quantum operator theory and matrix manifold analysis.
