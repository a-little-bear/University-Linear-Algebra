# Chapter 20: Matrix Equations

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Kronecker Product (Ch19) · Stability (Ch36) · Eigenvalues (Ch6)

**Chapter Outline**: Sylvester Equation ($AX - XB = C$) → Lyapunov Equation ($AX + XA^T = Q$) → Existence and Uniqueness of Solutions → Solvability Conditions (Spectral Disjointness) → Solving via Kronecker Products → Solving via Schur Decomposition (Bartels-Stewart Algorithm) → Algebraic Riccati Equation (ARE) Overview → Applications in Control Theory and Stability

**Extension**: Matrix equations are the "algebraic laws" of systems theory; solving them is the primary way to construct Lyapunov functions and design optimal controllers.

</div>

Matrix equations are algebraic relationships where the unknown is a matrix $X$. Unlike the standard $Ax=b$, these equations involve transformations acting on $X$ from both sides. The **Sylvester equation** $AX - XB = C$ is the most fundamental, governing the existence of similarity transformations and the decoupling of blocks. Its special case, the **Lyapunov equation**, is the primary tool for verifying the stability of linear systems. This chapter details the solvability conditions for these equations and introduces the **Bartels-Stewart algorithm**, the standard numerical method for their solution.

---

## 20.1 Sylvester and Lyapunov Equations

!!! definition "Definition 20.1 (Sylvester Equation)"
    The Sylvester equation is $AX - XB = C$, where $A \in M_m, B \in M_n, C \in M_{m \times n}$ are given and $X$ is unknown.

!!! theorem "Theorem 20.1 (Existence and Uniqueness)"
    The Sylvester equation has a unique solution $X$ if and only if $A$ and $B$ have no common eigenvalues ($\sigma(A) \cap \sigma(B) = \emptyset$).

---

## Exercises

1. **[Fundamentals] Rewrite $AX - XB = C$ as a standard linear system $Mx = c$.**
   ??? success "Solution"
       Using the vec-Kronecker identity (Ch19): $(I \otimes A - B^T \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$. The matrix $M = I \otimes A - B^T \otimes I$ has eigenvalues $\{\lambda_i(A) - \mu_j(B)\}$.

2. **[Solvability] Why is $\sigma(A) \cap \sigma(B) = \emptyset$ required for uniqueness?**
   ??? success "Solution"
       The eigenvalues of the Kronecker sum $M$ are $\lambda_i - \mu_j$. If $A$ and $B$ share an eigenvalue, then $\lambda_i - \mu_j = 0$ for some pair, making $M$ singular. A singular $M$ means the solution is either non-unique or does not exist.

3. **[Lyapunov] State the Continuous Lyapunov Equation and its role in stability.**
   ??? success "Solution"
       $A^T P + PA = -Q$. For a Hurwitz stable matrix $A$, there exists a unique $P \succ 0$ for any $Q \succ 0$. $P$ defines the Lyapunov function $V(x) = x^T P x$.

4. **[Calculation] Solve $AX - XA = 0$.**
   ??? success "Solution"
       The solutions are exactly the matrices that commute with $A$. This set forms the **centralizer** of $A$ and includes all polynomials in $A$.

5. **[Schur Method] Describe the Bartels-Stewart Algorithm.**
   ??? success "Solution"
       1. Reduce $A$ and $B$ to Schur form ($T_A, T_B$) using unitary transformations. 2. Transform $C$ accordingly. 3. Solve the resulting triangular system $T_A \tilde{X} - \tilde{X} T_B = \tilde{C}$ using back-substitution.

6. **[Riccati] How does the Algebraic Riccati Equation (ARE) differ from the Sylvester equation?**
   ??? success "Solution"
       The ARE, $A^T P + PA - PBR^{-1}B^T P + Q = 0$, is **quadratic** in the unknown $P$. It arises in optimal control (LQR) and requires more sophisticated solvers (like the Schur method on Hamiltonian matrices).

7. **[Block Decoupling] Show how the Sylvester equation is used to block-diagonalize $\begin{pmatrix} A & C \\ 0 & B \end{pmatrix}$.**
   ??? success "Solution"
       A similarity transform with $\begin{pmatrix} I & X \\ 0 & I \end{pmatrix}$ results in $\begin{pmatrix} A & AX-XB+C \\ 0 & B \end{pmatrix}$. If $X$ solves $AX-XB = -C$, the off-diagonal block becomes zero.

8. **[Inertia] State the relationship between the Lyapunov equation and the inertia of $A$.**
   ??? success "Solution"
       This is the Lyapunov Inertia Theorem: if $A^T P + PA = Q \succ 0$, then $A$ and $P$ have the same number of eigenvalues in the right and left half-planes.

9. **[Iterative] Mention an iterative method for solving matrix equations.**
   ??? success "Solution"
       The Smith iteration: $X_{k+1} = A^T X_k A + Q$ for the discrete Lyapunov equation. It converges if $\rho(A) < 1$.

10. **[Trace] If $A^T P + PA = -I$, express $\operatorname{tr}(P)$ in terms of $A$.**
    ??? success "Solution"
        $\operatorname{tr}(P) = \int_0^\infty \operatorname{tr}(e^{A^T t} e^{At}) dt = \int_0^\infty \|e^{At}\|_F^2 dt$. The trace of $P$ measures the total energy dissipation of the system.

## Chapter Summary

This chapter establishes the algebraic theory of matrix-valued unknowns:

1. **Spectral Disjointness**: Defined the universal condition for the existence and uniqueness of solutions to Sylvester-type equations.
2. **Stability Engine**: Formulated the Lyapunov equation as the definitive algebraic test for dynamical equilibrium.
3. **Triangular Reduction**: Introduced the Bartels-Stewart algorithm as the numerically stable path for solving high-dimensional matrix systems.
4. **Structural Decoupling**: Demonstrated the role of matrix equations in the canonical reduction and block-diagonalization of operators.
