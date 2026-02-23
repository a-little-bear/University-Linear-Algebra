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

1. **[Fundamentals] Compute the geometric mean $A \# B$ for $A = \operatorname{diag}(4, 1)$ and $B = \operatorname{diag}(1, 9)$.**
   ??? success "Solution"
       Since the matrices are diagonal, they commute. $A \# B = (AB)^{1/2} = \operatorname{diag}(\sqrt{4 \cdot 1}, \sqrt{1 \cdot 9}) = \operatorname{diag}(2, 3)$.

2. **[Riccati] Verify that $X = A \# B$ satisfies the Riccati equation $X A^{-1} X = B$.**
   ??? success "Solution"
       $A^{1/2}(A^{-1/2} B A^{-1/2})^{1/2} A^{1/2} \cdot A^{-1} \cdot A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2} = A^{1/2} (A^{-1/2} B A^{-1/2}) A^{1/2} = B$.

3. **[AM-GM-HM] Arrange the three classical means in increasing order under the Löwner partial order.**
   ??? success "Solution"
       Harmonic Mean $\preceq$ Geometric Mean $\preceq$ Arithmetic Mean.
       $2(A^{-1} + B^{-1})^{-1} \preceq A \# B \preceq \frac{A+B}{2}$.

4. **[Determinant] Prove that $\det(A \# B) = \sqrt{\det A \cdot \det B}$.**
   ??? success "Solution"
       $\det(A \# B) = \det(A^{1/2}) \det((A^{-1/2} B A^{-1/2})^{1/2}) \det(A^{1/2}) = \det(A) \sqrt{\det(A^{-1}B)} = \det(A) \sqrt{\det(A)^{-1} \det B} = \sqrt{\det A \det B}$.

5. **[Transformer] Show that if $C$ is invertible, $C^*(A \# B)C = (C^*AC) \# (C^*BC)$.**
   ??? success "Solution"
       This is the transformer equality axiom. Let $X = A \# B$. Then $X A^{-1} X = B$.
       Multiply by $C^*$ on the left and $C$ on the right: $C^* X A^{-1} X C = C^* B C$.
       Inserting $C^{-1}(C^*)^{-1}$: $(C^* X C) (C^* A C)^{-1} (C^* X C) = C^* B C$.
       Thus $C^* X C$ is the geometric mean of $C^* A C$ and $C^* B C$.

6. **[Multivariate] Define the Karcher mean for $k$ matrices $A_1, \dots, A_k$.**
   ??? success "Solution"
       The Karcher mean is the unique positive definite minimizer of the sum of squared Riemannian distances: $G = \arg\min_X \sum w_i \|\log(X^{-1/2} A_i X^{-1/2})\|_F^2$.

7. **[Riemannian Metric] What is the Riemannian metric on the manifold of positive definite matrices?**
   ??? success "Solution"
       At a point $X$, the inner product of two tangent vectors (Hermitian matrices) $H, K$ is $\langle H, K angle_X = \operatorname{tr}(X^{-1} H X^{-1} K)$.

8. **[Geodesic] What is the geodesic path connecting $A$ and $B$?**
   ??? success "Solution"
       The geodesic is the weighted geometric mean $\gamma(t) = A \#_t B = A^{1/2}(A^{-1/2} B A^{-1/2})^t A^{1/2}$.

9. **[Commutativity] If $A$ and $B$ commute, how does $A \# B$ simplify?**
   ??? success "Solution"
       $A \# B = A^{1/2} B^{1/2} = (AB)^{1/2}$.

10. **[Singular Case] Calculate $A \# B$ for $A = \operatorname{diag}(1, 0)$ and $B = \operatorname{diag}(0, 1)$.**
    ??? success "Solution"
        Using the limit $\epsilon 	o 0$: $(A+\epsilon I) \# (B+\epsilon I) = \operatorname{diag}(\sqrt{\epsilon(1+\epsilon)}, \sqrt{\epsilon(1+\epsilon)})$. As $\epsilon 	o 0$, the mean becomes the zero matrix.

## Chapter Summary

This chapter explores the synthesis of algebraic means and differential geometry:

1. **Axiomatic Foundation**: Identified the Kubo-Ando framework as the rigorous path to define matrix means.
2. **Spectral Interplay**: Linked matrix geometric means to Riccati equations and the theory of operator monotone functions.
3. **Geometric Synthesis**: Developed the Riemannian geometry of the PSD manifold, identifying the geometric mean as the geodesic midpoint.
4. **Multivariate Generalization**: Extended the two-variable mean to $k$ variables via Karcher's variational approach.
