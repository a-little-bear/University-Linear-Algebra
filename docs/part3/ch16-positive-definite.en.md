# Chapter 16: Positive Definite Matrices

<div class="context-flow" markdown>

**Prerequisites**: Quadratic Forms (Ch09) · Matrix Decompositions (Ch10) · Eigenvalues (Ch06)

**Chapter Outline**: Definition of Positive Definite (PD) and Positive Semi-definite (PSD) Matrices → The Five Equivalent Criteria (Eigenvalues, Minors, Quadratic Forms, Cholesky, Gram Matrices) → Algebraic Properties → Schur Complement and Definiteness → Löwner Partial Order ($\succeq$) → Matrix Square Roots → Covariance Matrices in Statistics → Convexity in Optimization

**Extension**: Positive definite matrices are the geometric core of convex optimization; they are not just the extension of "positive numbers" to the matrix dimension, but the bedrock of modern Control Theory, Finance, and Loss Functions in Machine Learning.

</div>

Positive definite matrices are among the most celebrated classes of matrices in linear algebra. They are structurally symmetric, spectrally pure (all eigenvalues are positive), and geometrically stable. In physics, they represent the energy ground state; in statistics, they characterize covariance; in mathematics, they are the natural vehicle for defining measures and distances.

---

## 16.1 Definition and Equivalent Criteria

!!! definition "Definition 16.1 (PD and PSD)"
    For a real symmetric matrix $A$:
    1.  **Positive Definite (PD)**: $\mathbf{x}^T A \mathbf{x} > 0$ for all non-zero vectors $\mathbf{x}$. Denoted $A \succ 0$.
    2.  **Positive Semi-definite (PSD)**: $\mathbf{x}^T A \mathbf{x} \ge 0$ for all vectors $\mathbf{x}$. Denoted $A \succeq 0$.

!!! theorem "Theorem 16.1 (The Five Criteria)"
    For a real symmetric matrix $A$, the following are equivalent to $A \succ 0$:
    1.  **Eigenvalues**: All eigenvalues $\lambda_i > 0$.
    2.  **Principal Minors (Sylvester)**: All leading principal minors are positive.
    3.  **Cholesky**: There exists a unique lower triangular matrix $L$ with positive diagonal such that $A = LL^T$.
    4.  **Gram Matrix**: There exists a matrix $B$ with full column rank such that $A = B^T B$.
    5.  **Energy**: The graph of the quadratic form $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$ is an upward-opening hyper-paraboloid.

---

## 16.2 Schur Complement and Definiteness

!!! theorem "Theorem 16.2 (Block Matrix Definiteness)"
    Let $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$. Then:
    $M \succ 0 \iff A \succ 0$ and the **Schur complement** $S = C - B^T A^{-1} B \succ 0$.
    **Application**: This is key for constrained optimization and analyzing large power systems.

---

## 16.3 Löwner Partial Order

!!! definition "Definition 16.2 (Löwner Order)"
    For symmetric matrices $A, B$, we define $A \succeq B \iff A - B \succeq 0$.
    **Properties**:
    1.  If $A \succeq B$ and $C \succeq 0$, then $A + C \succeq B$.
    2.  If $A \succeq B \succ 0$, then $B^{-1} \succeq A^{-1} \succ 0$ (operator inversion reverses order).

---

## Exercises

1. **[Criteria] Determine if $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ is positive definite.**
   ??? success "Solution"
       Leading principal minors: $D_1 = 2 > 0, D_2 = 4 - 1 = 3 > 0$. Thus, it is PD.

2. **[Determinant] If $A \succ 0$, must its determinant be positive?**
   ??? success "Solution"
       Yes. $\det(A) = \prod \lambda_i$, and since all $\lambda_i > 0$, the product is positive.

3. **[Property] Prove: If $A \succ 0$ and $B \succ 0$, then $A + B \succ 0$.**
   ??? success "Solution"
       $\mathbf{x}^T (A+B) \mathbf{x} = \mathbf{x}^T A \mathbf{x} + \mathbf{x}^T B \mathbf{x} > 0 + 0 = 0$ for $\mathbf{x} \neq \mathbf{0}$.

4. **[Inverse] If $A \succ 0$, prove $A^{-1} \succ 0$.**
   ??? success "Solution"
       If eigenvalues of $A$ are $\lambda_i > 0$, then eigenvalues of $A^{-1}$ are $1/\lambda_i > 0$.

5. **[Gram] If $A = B^T B$ but $B$ is not full column rank, what kind of matrix is $A$?**
   ??? success "Solution"
       Positive semi-definite (PSD). There exists non-zero $x$ such that $Bx = 0 \implies x^T A x = \|Bx\|^2 = 0$.

6. **[Diagonal] Prove: Diagonal entries of a PD matrix must be positive.**
   ??? success "Solution"
       Let $\mathbf{x} = \mathbf{e}_i$, then $\mathbf{e}_i^T A \mathbf{e}_i = a_{ii} > 0$.

7. **[Schur] Determine if $\begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$ is PD using Schur complement.**
   ??? success "Solution"
       $A = (1) \succ 0$. Schur complement $S = 5 - 2(1)^{-1}2 = 1 > 0$. Thus, the matrix is PD.

8. **[Root] If $A \succeq 0$, prove there exists a unique $B \succeq 0$ such that $B^2 = A$.**
   ??? success "Solution"
       Using diagonalization $A = Q \Lambda Q^T$, take $B = Q \Lambda^{1/2} Q^T$.

9. **[Order] If $A \succeq B \succ 0$, prove $\operatorname{tr}(A) \ge \operatorname{tr}(B)$.**
   ??? success "Solution"
       $A - B \succeq 0 \implies \operatorname{tr}(A - B) = \sum \lambda_i(A-B) \ge 0 \implies \operatorname{tr}(A) \ge \operatorname{tr}(B)$.

10. **[Statistics] Why is a covariance matrix always semi-definite?**
    ??? success "Solution"
        It is $E[(\mathbf{x}-\mu)(\mathbf{x}-\mu)^T]$. For any $v$, $v^T \Sigma v = E[(v^T(\mathbf{x}-\mu))^2] \ge 0$, representing non-negative variance.

## Chapter Summary

Positive definite matrices construct the convex geometry of high-dimensional space:

1.  **Multi-dimensional Positivity**: PD-ness is the perfect extension of "greater than zero" to operators, establishing the legitimacy of energy, probability, and distance.
2.  **Diverse Characterization**: From micro-level entry patterns (minors) to macro-level energy behavior (quadratic forms), and internal structural decomposition (Cholesky), these equivalent perspectives provide flexible tools for different fields.
3.  **Computational Supremacy**: PD matrices possess natural "monotonic stability" in numerical computing. The Löwner order allows us to treat matrix functions like scalar inequalities, forming the core of operator analysis.
