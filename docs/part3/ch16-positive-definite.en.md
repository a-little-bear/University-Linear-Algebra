# Chapter 16: Positive Definite Matrices

<div class="context-flow" markdown>

**Prerequisites**: Quadratic Forms (Ch09) · Matrix Decompositions (Ch10) · Eigenvalues (Ch06)

**Chapter Outline**: Definition of Positive Definite (PD) and Positive Semi-definite (PSD) Matrices → The Five Equivalent Criteria (Eigenvalues, Principal Minors, Quadratic Forms, Cholesky, Gram Matrices) → Algebraic Properties of PD Matrices → Matrix Square Roots → Schur Complement and PD Criteria → Inequalities for Matrix Variables: The Löwner Partial Order ($\succeq$) → Applications: Covariance Matrices in Statistics, Convexity in Optimization, and Stiffness Matrices in Engineering Mechanics

**Extension**: Positive definite matrices are the geometric core of convex optimization; they are not just the extension of "positive numbers" to matrix dimensions but the bedrock of modern Control Theory, Finance, and Loss Functions in Machine Learning.

</div>

Positive definite matrices are among the most celebrated classes of matrices in linear algebra. They are structurally symmetric, spectrally pure (all eigenvalues are positive), and geometrically stable. In physics, they represent the energy ground state; in statistics, they characterize covariance. This chapter establishes multi-dimensional criteria for positive definiteness and introduces partial orders to describe matrix "magnitude."

---

## 16.1 Definition and the Five Criteria

!!! definition "Definition 16.1 (PD and PSD)"
    For a real symmetric matrix $A \in S_n$:
    1.  **Positive Definite (PD)**: $\mathbf{x}^T A \mathbf{x} > 0$ for all $\mathbf{x} \neq \mathbf{0}$. Denoted $A \succ 0$.
    2.  **Positive Semi-definite (PSD)**: $\mathbf{x}^T A \mathbf{x} \ge 0$ for all $\mathbf{x}$. Denoted $A \succeq 0$.

!!! theorem "Theorem 16.1 (Equivalence Criteria)"
    The following are equivalent for a real symmetric matrix $A$:
    1.  **Eigenvalues**: All eigenvalues $\lambda_i > 0$.
    2.  **Leading Principal Minors**: All $k$-th order leading principal minors are positive.
    3.  **Cholesky Factorization**: There exists a unique lower triangular matrix $L$ (with positive diagonals) such that $A = LL^T$.
    4.  **Gram Matrix**: There exists a matrix $B$ with full column rank such that $A = B^T B$.
    5.  **Energy**: The surface defined by the quadratic form is an upward-opening hyper-paraboloid.

---

## 16.2 Matrix Square Roots and Schur Complements

!!! theorem "Theorem 16.2 (Matrix Square Root)"
    If $A \succeq 0$, there exists a unique positive semi-definite matrix $B$ such that $B^2 = A$. We denote this as $B = A^{1/2}$.

!!! technique "Technique: Schur Complement Criterion"
    A partitioned matrix $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$ is positive definite iff $A \succ 0$ and the **Schur Complement** $S = C - B^T A^{-1} B \succ 0$.

---

## 16.3 The Löwner Partial Order

!!! definition "Definition 16.2 (Löwner Order)"
    For symmetric matrices $A, B$, we define $A \succeq B \iff A - B \succeq 0$.
    **Properties**:
    1.  If $A \succeq B$ and $C \succeq 0$, then $A + C \succeq B$.
    2.  If $A \succeq B \succ 0$, then $B^{-1} \succeq A^{-1} \succ 0$ (operator inversion reverses the order).

---

## Exercises

**1. [Criteria] Determine if $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ is positive definite.**

??? success "Solution"
    **Calculate Principal Minors:**
    1. 1st order: $D_1 = 2 > 0$.
    2. 2nd order: $D_2 = 2 \cdot 2 - (-1) \cdot (-1) = 3 > 0$.
    **Conclusion**: Since all leading principal minors are positive, the matrix is **positive definite**.

**2. [Determinant] If $A$ is positive definite, must its determinant $\det(A)$ be positive?**

??? success "Solution"
    **Analysis:**
    1. A positive definite matrix has all eigenvalues $\lambda_i > 0$.
    2. The determinant $\det(A) = \prod_{i=1}^n \lambda_i$.
    3. The product of multiple positive numbers is always positive.
    **Conclusion**: Yes, $\det(A) > 0$. Note: A positive determinant is a necessary but not sufficient condition (e.g., consider $\operatorname{diag}(-1, -1)$).

**3. [Properties] Prove: If $A \succ 0$ and $B \succ 0$, then $A + B \succ 0$.**

??? success "Solution"
    **Proof:**
    1. For any non-zero vector $\mathbf{x}$, consider the quadratic form $\mathbf{x}^T (A+B) \mathbf{x}$.
    2. By linearity: $= \mathbf{x}^T A \mathbf{x} + \mathbf{x}^T B \mathbf{x}$.
    3. Since $A, B \succ 0$, we have $\mathbf{x}^T A \mathbf{x} > 0$ and $\mathbf{x}^T B \mathbf{x} > 0$.
    4. The sum of two positive numbers is positive.
    **Conclusion**: The set of positive definite matrices is closed under addition.

**4. [Inverse] If $A \succ 0$, prove its inverse $A^{-1}$ exists and $A^{-1} \succ 0$.**

??? success "Solution"
    **Proof:**
    1. Since $\lambda_i > 0$, the determinant is non-zero, so the inverse exists.
    2. The eigenvalues of $A^{-1}$ are $1/\lambda_i$.
    3. Since $\lambda_i > 0$, their reciprocals $1/\lambda_i$ are also positive.
    4. By the eigenvalue criterion, $A^{-1}$ is positive definite.

**5. [Gram Matrix] Given $A = B^T B$ where $B$ has rank $r$. Under what condition on $r$ is $A$ positive definite?**

??? success "Solution"
    **Determination:**
    1. $A$ is PD iff its kernel is empty, i.e., $Ax=0 \implies x=0$.
    2. $\mathbf{x}^T A \mathbf{x} = \|B\mathbf{x}\|^2$. If $B\mathbf{x}=0$, then $\mathbf{x}^T A \mathbf{x} = 0$.
    3. To ensure $\mathbf{x}^T A \mathbf{x} > 0$ for all non-zero $\mathbf{x}$, we must have $B\mathbf{x} \neq 0$.
    4. This requires $B$ to have **full column rank**.
    **Conclusion**: If $B$ is $m \times n$, we require $r = n$.

**6. [Diagonal Entries] Prove: The diagonal entries of a PD matrix must be positive.**

??? success "Solution"
    **Proof:**
    1. Take the standard basis vector $\mathbf{e}_i$.
    2. By the definition of positive definiteness, $\mathbf{e}_i^T A \mathbf{e}_i > 0$.
    3. Expanding gives $\mathbf{e}_i^T A \mathbf{e}_i = a_{ii}$.
    **Conclusion**: $a_{ii} > 0$ for all $i$.

**7. [Schur Complement] Determine if $\begin{pmatrix} 1 & 2 \\ 2 & 5 \end{pmatrix}$ is PD using the Schur complement.**

??? success "Solution"
    **Steps:**
    1. Top-left $A = (1)$. Since $1 > 0$, the initial condition is met.
    2. Calculate Schur complement $S = C - B^T A^{-1} B = 5 - 2(1)^{-1}2 = 5 - 4 = 1$.
    3. Since $S = 1 > 0$, the Schur complement criterion is satisfied.
    **Conclusion**: The matrix is positive definite.

**8. [Square Root] Find the square root of $A = \begin{pmatrix} 4 & 0 \\ 0 & 9 \end{pmatrix}$.**

??? success "Solution"
    **Conclusion:**
    $A^{1/2} = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$.
    **Analysis**: For diagonal matrices, the square root is simply the square root of each diagonal entry.

**9. [Order] If $A \succeq B \succ 0$, prove $\operatorname{tr}(A) \ge \operatorname{tr}(B)$.**

??? success "Solution"
    **Proof:**
    1. $A - B \succeq 0$ implies all eigenvalues $\mu_i$ of $A-B$ are non-negative.
    2. Trace is the sum of eigenvalues: $\operatorname{tr}(A - B) = \sum \mu_i \ge 0$.
    3. By linearity: $\operatorname{tr}(A) - \operatorname{tr}(B) \ge 0 \implies \operatorname{tr}(A) \ge \operatorname{tr}(B)$.

**10. [Application] Why is a covariance matrix always semi-definite?**

??? success "Solution"
    **Statistical Reasoning:**
    1. The quadratic form $v^T \Sigma v$ of a covariance matrix $\Sigma$ represents the **variance** of the linear combination $v^T \mathbf{X}$.
    2. Variance is non-negative by definition: $\operatorname{Var}(Y) = E[(Y-E[Y])^2] \ge 0$.
    **Algebraic Conclusion**: Since the variance of any linear combination cannot be negative, the covariance matrix must be positive semi-definite.

## Chapter Summary

Positive definite matrices construct the convex geometry of high-dimensional space:

1.  **Multi-dimensional Positivity**: PD-ness is the perfect extension of "greater than zero" to operators, establishing the legitimacy of energy, probability, and distance.
2.  **Diverse Characterization**: From micro-level entry patterns (minors) to macro-level energy behavior (quadratic forms), and internal structural decomposition (Cholesky), these equivalent perspectives provide tools for different fields.
3.  **Computational Supremacy**: PD matrices possess natural "stability" in numerical computing. The Löwner order allows us to treat matrix functions like scalar inequalities, forming the core of operator analysis.
