# Chapter 46B: Matrix Means and Geometric Means

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Operator Monotone Functions (Ch46A) · Löwner Partial Order

**Chapter Outline**: Kubo-Ando Axioms → Kubo-Ando Theorem → Operator Connections → Classical Means (Arithmetic, Harmonic, Geometric) → Equivalent Definitions of Geometric Mean (Riccati Equation, Variational Characterization) → AM-GM-HM Inequality → Multivariate Means (Ando-Li-Mathias, Karcher) → Riemannian Geometry of the Positive Definite Manifold

**Extension**: Matrix geometric means are applied in medical imaging (diffusion tensor MRI), signal processing, and quantum fidelity.

</div>

What is the "middle point" between two positive definite matrices $A$ and $B$? Unlike scalars, the simple square root $\sqrt{AB}$ is not generally Hermitian. The correct definition, $A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$, emerges from the **Kubo-Ando theory** (1980), which identifies matrix means with normalized operator monotone functions. This chapter develops the axiomatic theory of means and explores the rich Riemannian geometry of the positive definite manifold.

---

## 46B.1 Core Concepts

!!! definition "Definition 46B.1 (Kubo-Ando Axioms)"
    A matrix mean $\sigma$ is a binary operation on $A, B \succ 0$ satisfying:
    1. **Monotonicity**: $A \preceq A', B \preceq B' \Rightarrow A \sigma B \preceq A' \sigma B'$.
    2. **Transformer Inequality**: $C^*(A \sigma B)C \preceq (C^*AC) \sigma (C^*BC)$.
    3. **Normalization**: $I \sigma I = I$.

!!! theorem "Theorem 46B.1 (Kubo-Ando Theorem)"
    There is a 1-to-1 correspondence between matrix means and normalized operator monotone functions $f$:
    $$A \sigma B = A^{1/2} f(A^{-1/2} B A^{-1/2}) A^{1/2}.$$

---

## Exercises

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
         $2(A^{-1} + B^{-1})^{-1} \preceq A \# B \preceq \frac{A+B}{2}$.

****
??? success "Solution"
    ****
??? success "Solution"
         Multiply by $C^*$ on the left and $C$ on the right: $C^* X A^{-1} X C = C^* B C$.
    Inserting $C^{-1}(C^*)^{-1}$: $(C^* X C) (C^* A C)^{-1} (C^* X C) = C^* B C$.
    Thus $C^* X C$ is the geometric mean of $C^* A C$ and $C^* B C$.

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
    ## Chapter Summary

This chapter explores the synthesis of algebraic means and differential geometry:


****: Identified the Kubo-Ando framework as the rigorous path to define matrix means.

****: Linked matrix geometric means to Riccati equations and the theory of operator monotone functions.

****: Developed the Riemannian geometry of the PSD manifold, identifying the geometric mean as the geodesic midpoint.

****: Extended the two-variable mean to $k$ variables via Karcher's variational approach.
