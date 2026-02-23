# Chapter 20: Matrix Equations

<div class="context-flow" markdown>

**Prerequisites**: Kronecker Product and Vec Operator (Ch19) · Matrix Analysis (Ch14) · Eigenvalues (Ch06)

**Chapter Outline**: From Scalar to Matrix Variable Equations → Basic Linear Equations $AX=B$ and $AXB=C$ → Existence and Uniqueness of the Sylvester Equation ($AX+XB=C$) via Spectral Separation → Lyapunov Equation ($AX+XA^T=Q$) and Stability Criteria → Discrete vs. Continuous Cases → Non-linear Matrix Equations: The Algebraic Riccati Equation (ARE) → Computational Techniques: Bartels-Stewart Algorithm → Applications: Optimal Control (LQR), Kalman Filtering, and Model Reduction

**Extension**: Matrix equations are the link between modern control theory and numerical linear algebra; they solve the problem of linear interference between operators and provide the analytic path for optimal feedback gains through Riccati equations.

</div>

When the unknown itself is a matrix, we call the expression a **Matrix Equation**. Such equations are not just natural extensions of linear operator theory but direct products of dynamical system equilibrium analysis, control gain calculation, and complex numerical simulations. Starting from simple linear coupled equations and utilizing the vectorization tools from Ch19, this chapter progresses into the non-linear Riccati equations central to control theory.

---

## 20.1 Linear Matrix Equations and Vectorization

!!! definition "Definition 20.1 (Basic Linear Matrix Equations)"
    1.  **Left-multiplication $AX = B$**: Solvable iff every column of $B$ lies in the column space of $A$.
    2.  **Two-sided $AXB = C$**: Can be transformed into $(B^T \otimes A) \operatorname{vec}(X) = \operatorname{vec}(C)$. Its unique solvability depends on whether $B^T \otimes A$ is non-singular.

---

## 20.2 Sylvester and Lyapunov Equations

!!! definition "Definition 20.2 (Sylvester Equation)"
    The equation takes the form: $AX + XB = C$, where $A \in M_n, B \in M_m, C \in M_{n \times m}$.

!!! theorem "Theorem 20.1 (Uniqueness Criterion)"
    The Sylvester equation $AX + XB = C$ has a unique solution for every $C$ iff $A$ and $-B$ have **no common eigenvalues**:
    $$\sigma(A) \cap \sigma(-B) = \emptyset$$
    **Significance**: This spectral separation condition ensures that the Kronecker sum operator $A \oplus B^T$ is invertible.

---

## 20.3 Algebraic Riccati Equation (ARE)

!!! definition "Definition 20.3 (Continuous-time ARE)"
    $$A^T X + XA - X B R^{-1} B^T X + Q = 0$$
    This is a quadratic non-linear matrix equation.
    **Application**: Used to solve for the optimal state feedback gain in Linear Quadratic Regulator (LQR) problems.

---

## Exercises

**1. [Sylvester] Determine if $\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} X + X \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = C$ has a unique solution.**

??? success "Solution"
    **Determination Steps:**
    1. Identify matrices: $A = I, B = I$.
    2. Calculate spectra: $\sigma(A) = \{1\}, \sigma(-B) = \{-1\}$.
    3. Check intersection: $\{1\} \cap \{-1\} = \emptyset$.
    **Conclusion**: The spectral separation condition is satisfied, so there is a unique solution for any $C$.

**2. [Uniqueness] If $A$ and $B$ are both diagonal and share an eigenvalue, is the Sylvester equation necessarily unsolvable?**

??? success "Solution"
    **Analysis:**
    1. Not necessarily.
    2. If the spectra overlap, the vectorized matrix $M = I \otimes A + B^T \otimes I$ is singular (it has a zero eigenvalue).
    3. This means the equation either has infinitely many solutions (if $\operatorname{vec}(C)$ is in the column space of $M$) or no solution.
    **Conclusion**: We can only say the solution is not unique; we cannot immediately conclude there is no solution.

**3. [Lyapunov] Why is the Lyapunov equation $AX + XA^T = Q$ commonly used for stability analysis?**

??? success "Solution"
    **Physical Logic:**
    1. The solution $X$ forms a generalized energy function $V(x) = x^T X x$.
    2. Differentiating this yields $\dot{V} = x^T(AX + XA^T)x = x^T Q x$.
    3. If we can find $X \succ 0$ such that $\dot{V} < 0$ (meaning $Q \prec 0$), then the system is stable.
    **Algebraic Mapping**: The Lyapunov equation transforms "stability of differential equations" into a "positivity problem of linear matrix equations."

**4. [Vectorization] Transform $AXB + CXD = F$ into vectorized form.**

??? success "Solution"
    **Steps:**
    1. Apply the vectorization formula $\operatorname{vec}(MXN) = (N^T \otimes M) \operatorname{vec}(X)$ to each term.
    2. First term: $(B^T \otimes A) \operatorname{vec}(X)$.
    3. Second term: $(D^T \otimes C) \operatorname{vec}(X)$.
    **Conclusion**: $(B^T \otimes A + D^T \otimes C) \operatorname{vec}(X) = \operatorname{vec}(F)$.

**5. [Riccati] Prove: If $X$ is a solution to the Riccati equation, then its transpose $X^T$ is also a solution (assuming $Q, R$ are symmetric).**

??? success "Solution"
    **Proof:**
    1. Transpose both sides of $A^T X + XA - XBR^{-1}B^T X + Q = 0$.
    2. $(A^T X)^T + (XA)^T - (XBR^{-1}B^T X)^T + Q^T = 0$.
    3. $X^T A + A^T X^T - X^T B R^{-1} B^T X^T + Q = 0$.
    4. Observe that the resulting equation for $X^T$ is identical in form to the original equation for $X$.
    **Conclusion**: $X^T$ satisfies the same equation.

**6. [Calculation] Solve $(1)X + X(2) = (6)$ (scalar case).**

??? success "Solution"
    **Calculation:**
    $1 \cdot X + X \cdot 2 = 6 \implies 3X = 6 \implies X = 2$.
    This shows how matrix equations converge to standard linear equations in one dimension.

**7. [Property] If $A$ and $B$ are upper triangular, is the solution $X$ to the Sylvester equation necessarily upper triangular?**

??? success "Solution"
    **Conclusion**: Generally no.
    **Reasoning**: The structure of $X$ depends not only on $A$ and $B$ but heavily on the structure of the constant term $C$. If $C$ is a dense matrix, $X$ will typically be dense even if $A$ and $B$ are diagonal.

**8. [Discrete Case] Write the form of the Discrete Lyapunov Equation.**

??? success "Solution"
    **Formula:**
    $$A X A^T - X = Q$$
    or $A X A^T - X + Q = 0$. It describes the stability of the difference equation $x_{k+1} = Ax_k$.

**9. [Numerical Algorithm] Briefly explain the core idea of the Bartels-Stewart algorithm for solving $AX+XB=C$.**

??? success "Solution"
    **Idea:**
    1. Use Schur decomposition to transform $A$ into upper triangular and $B$ into lower triangular form.
    2. In the transformed basis, the equation becomes entry-wise coupled.
    3. We can solve for the elements of $X$ one by one, starting from one corner of the matrix, using a method similar to back-substitution.
    **Benefit**: Complexity is $O(n^3)$, far lower than the $O(n^6)$ required for direct vectorization.

**10. [Application] In Kalman filtering, is the prediction step covariance update a matrix equation?**

??? success "Solution"
    **Yes.**
    The covariance update $P_{k|k-1} = A P_{k-1|k-1} A^T + Q$ is essentially a linear matrix mapping. Finding the steady-state Kalman gain requires solving an algebraic Riccati equation.

## Chapter Summary

Matrix equations are the algebraic language of advanced linear systems:

1.  **Solution Coupling**: The Sylvester equation characterizes the linear interference between two different operators; its spectral separation condition reveals the essence of system resonance.
2.  **Energy and Stability**: The Lyapunov equation builds a bridge between algebra and analytic stability through quadratic matrices, forming the basis for feedback design in control engineering.
3.  **Optimal Search**: The Riccati equation demonstrates how non-linear structures naturally emerge in variational and optimal control problems, marking the transition from linear system analysis to linear system synthesis.
