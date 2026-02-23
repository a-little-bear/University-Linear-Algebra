# Chapter 20: Matrix Equations

<div class="context-flow" markdown>

**Prerequisites**: Kronecker Product (Ch19) · Eigenvalues (Ch6) · Matrix Analysis (Ch14) · Matrix Stability (Ch36)

**Chapter Outline**: Sylvester Equation $AX - XB = C$ → Existence and Uniqueness Conditions → Lyapunov Equation $AX + XA^* = Q$ → Stability Determination and Inertia Theorem → Algebraic Riccati Equation (ARE) → Iterative Solution Methods → Applications (Controllability and Observability in Modern Control Theory)

**Extension**: Matrix equations are the language for transforming static algebra into dynamic characteristics; the solution to the Lyapunov equation directly maps the "energy" distribution of a dynamical system.

</div>

Matrix equations study algebraic equations where the unknowns are matrices. They typically appear in control system design, equilibrium analysis, and numerical analysis. The most famous, the Sylvester and Lyapunov equations, provide the mathematical framework for understanding the interactions between operators.

---

## 20.1 Sylvester and Lyapunov Equations

!!! definition "Definition 20.1 (Sylvester Equation)"
    The matrix equation $AX - XB = C$ is called the Sylvester equation.

!!! theorem "Theorem 20.3 (Existence and Uniqueness)"
    The Sylvester equation $AX - XB = C$ has a unique solution if and only if $A$ and $B$ have no common eigenvalues: $\sigma(A) \cap \sigma(B) = \emptyset$.

---

## Exercises

1. **[Existence] Does the equation $AX - XA = 0$ always have non-zero solutions?**
   ??? success "Solution"
       Yes. Any matrix that commutes with $A$ is a solution. In particular, $X=I$ and $X=A$ are always solutions. Since $\sigma(A) \cap \sigma(A) \neq \emptyset$ (unless the set is empty), the equation does not have a unique solution by Theorem 20.3.

2. **[Lyapunov] Assume $A$ is Hurwitz stable. Prove that for any $Q$, the solution to $AX + XA^* = Q$ can be expressed in integral form.**
   ??? success "Solution"
       The solution is $X = -\int_0^\infty e^{At} Q e^{A^*t} dt$. Since $A$ is stable, the exponential terms decay over time, guaranteeing the convergence of the integral.

3. **[Calculation] Solve $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} X - X \begin{pmatrix} 0 & 0 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Let $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$.
       Substitute: $\begin{pmatrix} 1x_{11} & 1x_{12} \\ 2x_{21} & 2x_{22} \end{pmatrix} - \begin{pmatrix} 0 & 3x_{12} \\ 0 & 3x_{22} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.
       Solve the system: $x_{11}=1, -2x_{12}=1, 2x_{21}=1, -x_{22}=1$.
       Thus $X = \begin{pmatrix} 1 & -0.5 \\ 0.5 & -1 \end{pmatrix}$.

4. **[Stability] Prove: If there exists $P \succ 0$ such that $A^T P + PA \prec 0$, then $A$ is stable.**
   ??? success "Solution"
       This is the Lyapunov stability theorem. Consider the energy function $V(x) = x^T P x$.
       Its derivative is $\dot{V} = x^T (A^T P + PA) x$. Since $A^T P + PA$ is negative definite, $\dot{V} < 0$, which means energy decreases over time and the system converges.

5. **[Trace Application] In the Lyapunov equation $AX + XA^T + BB^T = 0$, what does the trace of $X$ represent?**
   ??? success "Solution"
       If the equation describes the controllability Gramian of a control system, the trace of $X$ (sum of eigenvalues) represents the "average gain" or overall degree of controllability for control energy input from all directions.

6. **[Riccati Intro] What is the Algebraic Riccati Equation (ARE)? How does it differ from the Sylvester equation?**
   ??? success "Solution"
       ARE has the form $A^T P + PA - PBR^{-1}B^T P + Q = 0$. Unlike the linear Sylvester equation, ARE is a **quadratic** matrix equation. It plays a central role in optimal control (LQR).

7. **[Controllability] If the system controllability matrix is $W_c$, what matrix equation does it satisfy?**
   ??? success "Solution"
       It satisfies the Lyapunov equation $A W_c + W_c A^T + BB^T = 0$. The positive definiteness of the solution $W_c$ directly corresponds to the system's controllability.

8. **[Symmetry] If $A$ is stable and $Q$ is symmetric, is the solution $X$ to $AX + XA^T = Q$ always symmetric?**
   ??? success "Solution"
       Yes. Taking the transpose of the original equation gives $XA^T + AX^T = Q^T = Q$. Since the solution is unique and both $X$ and $X^T$ satisfy the same equation, $X = X^T$.

9. **[Elementary Operator] The Sylvester equation can be viewed as an inverse problem for which linear operator?**
   ??? success "Solution"
       It can be viewed as inverting the elementary operator $\mathcal{L}(X) = AX - XB$ acting on the space of matrices. Its eigenvalues are $\lambda_i(A) - \mu_j(B)$.

10. **[Application] How do matrix equations model degradation in image restoration?**
    ??? success "Solution"
        Blurring is often modeled as $Y = AXB + N$, where $A$ and $B$ represent horizontal and vertical blurring operators. Restoring the image involves solving this matrix equation (usually with a regularization term).

## Chapter Summary

Matrix equations are the dynamic logic of linear algebra:

1. **Spectral Interaction**: The existence of a solution depends on the separation of the spectra of the two operators.
2. **Energy Mapping**: Lyapunov equations build a direct bridge between algebra and dynamical stability.
3. **Optimization Core**: From linear to quadratic equations (Riccati), matrix equations form the computational base of modern control theory.
