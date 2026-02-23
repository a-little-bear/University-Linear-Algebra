# Chapter 20: Matrix Equations

<div class="context-flow" markdown>

**Prerequisites**: Kronecker Product (Ch19) · Matrix Analysis (Ch14) · Eigenvalues (Ch06)

**Chapter Outline**: Overview of Matrix Equations → Basic Linear Equations $AX=B$ and $AXB=C$ → Sylvester Equation ($AX+XB=C$) → Lyapunov Equation ($AX+XA^T=Q$) → Criteria for Existence and Uniqueness (Spectral Separation) → Algebraic Riccati Equation (ARE) → Continuous vs. Discrete Cases → Numerical Algorithms (Bartels-Stewart Algorithm) → Applications in Control Theory (LQR, Stability)

**Extension**: Matrix equations are the link between modern control theory and numerical linear algebra; they are the core mathematical tools for analyzing system observability, controllability, and solving optimal control strategies.

</div>

When the unknown itself is a matrix, we call it a matrix equation. Matrix equations are not just an extension of linear operator theory but also the direct product of dynamical system equilibrium analysis, control gain calculation, and numerical simulation. This chapter begins with simple linear coupled equations and progresses into complex non-linear Riccati equations.

---

## 20.1 Linear Matrix Equations

!!! definition "Definition 20.1 (Basic Equations)"
    1.  **Left-multiplication $AX = B$**: Solvable if and only if the columns of $B$ belong to the column space of $A$.
    2.  **Two-sided $AXB = C$**: Can be transformed via the Kronecker product into $(B^T \otimes A) \operatorname{vec}(X) = \operatorname{vec}(C)$.

---

## 20.2 Sylvester and Lyapunov Equations

!!! definition "Definition 20.2 (Sylvester Equation)"
    $$AX + XB = C$$
    where $A, B, C$ are given square matrices.

!!! theorem "Theorem 20.1 (Uniqueness Criterion)"
    The Sylvester equation $AX + XB = C$ has a unique solution for every $C$ iff $\sigma(A) \cap \sigma(-B) = \emptyset$ (i.e., $A$ and $-B$ have no common eigenvalues).

!!! definition "Definition 20.3 (Lyapunov Equation)"
    $$AX + XA^T = Q$$
    This is a special symmetric form of the Sylvester equation used to determine the stability of dynamical systems. If $A$ is stable (spectrum in the left half-plane), then for any $Q \prec 0$, there is a unique positive definite solution $X \succ 0$.

---

## 20.3 Algebraic Riccati Equation (ARE)

!!! definition "Definition 20.4 (Continuous-time ARE)"
    $$A^T X + XA - X B R^{-1} B^T X + Q = 0$$
    This is a quadratic non-linear matrix equation.
    **Application**: Used to solve for the optimal feedback gain $K = R^{-1} B^T X$ in Linear-Quadratic Regulator (LQR) problems.

---

## 20.4 Numerical Algorithms

!!! algorithm "Algorithm 20.1 (Bartels-Stewart Algorithm)"
    Used for efficient solution of linear matrix equations:
    1.  Perform Schur decomposition (triangularization) on $A$ and $B$.
    2.  Solve the transformed equation via forward/backward substitution.
    3.  Back-transform to obtain the solution to the original equation.
    **Complexity**: $O(n^3)$, much faster than direct vectorization ($O(n^6)$).

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

Matrix equations are the algebraic language of advanced linear systems:


****: The Sylvester equation characterizes the linear interference between two different operators; its spectral separation condition reveals the essence of system resonance.

****: The Lyapunov equation builds a bridge between algebra and analytical stability through quadratic matrices, serving as the operator expression of Lyapunov's second method in control engineering.

****: The Riccati equation demonstrates how non-linear structures naturally arise in variational and optimal control problems, marking the transition from linear system analysis to linear system synthesis.
