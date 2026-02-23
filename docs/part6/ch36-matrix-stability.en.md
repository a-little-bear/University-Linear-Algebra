# Chapter 36: Matrix Stability

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Matrix Exponential (Ch13) · Differential Equations (Ch26) · Lyapunov Equations (Ch20)

**Chapter Outline**: Stability of Continuous-time Systems (Hurwitz) → Stability of Discrete-time Systems (Schur) → Lyapunov Stability Theorem → Lyapunov Equation $A^T P + PA = -Q$ → Routh-Hurwitz Criterion → Jury Stability Criterion → Stability Radius → Robust Stability

**Extension**: Matrix stability is the mathematical prerequisite for control theory (Ch66) and the analysis of equilibrium points in economic and biological systems.

</div>

Matrix stability concerns the behavior of trajectories in linear dynamical systems. A matrix is **Hurwitz stable** if all its eigenvalues lie in the open left-half complex plane, ensuring that solutions to $\dot{x} = Ax$ decay to zero. It is **Schur stable** if all eigenvalues lie within the unit circle, ensuring that $x_{k+1} = Ax_k$ is stable. The **Lyapunov equation** provides a way to verify stability using positive definite matrices without explicitly calculating eigenvalues.

---

## 36.1 Hurwitz and Schur Stability

!!! definition "Definition 36.1 (Hurwitz Matrix)"
    A matrix $A \in M_n(\mathbb{C})$ is **Hurwitz** if $\operatorname{Re}(\lambda_i) < 0$ for all eigenvalues $\lambda_i \in \sigma(A)$.

!!! definition "Definition 36.2 (Schur Matrix)"
    A matrix $A \in M_n(\mathbb{C})$ is **Schur** if $|\lambda_i| < 1$ for all eigenvalues $\lambda_i \in \sigma(A)$.

!!! theorem "Theorem 36.1 (Lyapunov Stability Theorem)"
    $A$ is Hurwitz stable if and only if for any $Q \succ 0$, there exists a unique $P \succ 0$ such that:
    $$A^T P + PA = -Q$$

---

## Exercises

1. **[Fundamentals] Is $A = \begin{pmatrix} -1 & 10 \ 0 & -2 \end{pmatrix}$ Hurwitz stable?**
   ??? success "Solution"
       Yes. The eigenvalues of an upper triangular matrix are its diagonal entries: $\{-1, -2\}$. Since both have negative real parts, $A$ is Hurwitz.

2. **[Schur Stability] Check if $A = \begin{pmatrix} 0.5 & 0.5 \ 0.5 & 0.5 \end{pmatrix}$ is Schur stable.**
   ??? success "Solution"
       The eigenvalues are $\{1, 0\}$. Since one eigenvalue lies on the unit circle ($|1|=1$), $A$ is not Schur stable (it is marginally stable).

3. **[Lyapunov Equation] Why must $P$ be positive definite in the Lyapunov equation?**
   ??? success "Solution"
       $V(x) = x^T P x$ serves as a energy-like function (Lyapunov function). $P \succ 0$ ensures $V(x) \ge 0$ and $V(x) 	o \infty$ as $\|x\| 	o \infty$. The equation $\dot{V} = -x^T Q x < 0$ ensures that energy always decreases until the origin is reached.

4. **[Trace and Det] What do $\operatorname{tr}(A)$ and $\det(A)$ tell you about Hurwitz stability?**
   ??? success "Solution"
       For a Hurwitz matrix, $\operatorname{tr}(A) = \sum \operatorname{Re}(\lambda_i) < 0$. Also, $(-1)^n \det A > 0$. These are necessary but not sufficient conditions.

5. **[Routh-Hurwitz] For a $2 	imes 2$ real matrix, state the necessary and sufficient conditions for Hurwitz stability in terms of its coefficients.**
   ??? success "Solution"
       For $p(\lambda) = \lambda^2 + a_1 \lambda + a_0$, stability requires $a_1 > 0$ and $a_0 > 0$. In matrix terms: $\operatorname{tr}(A) < 0$ and $\det(A) > 0$.

6. **[Bilinear Transform] How can you convert a Hurwitz stability problem into a Schur stability problem?**
   ??? success "Solution"
       Using the Cayley transform (bilinear transform) $z = \frac{1+s}{1-s}$. The left-half plane in $s$ is mapped to the unit disk in $z$.

7. **[Positive Matrices] Prove that a Metzler matrix (off-diagonals $\ge 0$) is Hurwitz iff there exists $d > 0$ such that $Ad < 0$.**
   ??? success "Solution"
       This is a property of M-matrices. Stability for positive systems can be verified using a linear (rather than quadratic) Lyapunov function.

8. **[Stability Radius] Define the complex stability radius $r_{\mathbb{C}}(A)$.**
   ??? success "Solution"
       $r_{\mathbb{C}}(A) = \min \{ \|\Delta\| : A+\Delta 	ext{ is unstable} \}$. It is equal to $1 / \sup_{\omega} \|(i\omega I - A)^{-1}\|$.

9. **[Commuting Families] If $A$ and $B$ are Hurwitz and $AB=BA$, is $A+B$ Hurwitz?**
   ??? success "Solution"
       Yes. If they commute, they can be simultaneously triangularized. The eigenvalues of $A+B$ are $\lambda_i(A) + \lambda_i(B)$, and the sum of two numbers with negative real parts has a negative real part.

10. **[Discrete Lyapunov] State the Lyapunov equation for Schur stability.**
    ??? success "Solution"
        $A^T P A - P = -Q$. Here $Q \succ 0$ implies $P \succ 0$ iff $A$ is Schur stable. This corresponds to $V(x_{k+1}) - V(x_k) = -x_k^T Q x_k$.

## Chapter Summary

This chapter establishes the criteria for the asymptotic convergence of linear systems:

1. **Spectral Domains**: Defined Hurwitz and Schur stability based on the location of eigenvalues relative to the imaginary axis and the unit circle.
2. **Lyapunov's Method**: Formulated stability as the existence of a positive definite solution to a linear matrix equation, avoiding eigenvalue computation.
3. **Polynomial Criteria**: Explored the Routh-Hurwitz and Jury tests for assessing stability directly from the characteristic equation.
4. **Robustness Metrics**: Introduced the stability radius to quantify the distance to the boundary of instability under model perturbations.
