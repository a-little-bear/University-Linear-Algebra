# Chapter 36: Matrix Stability and Inertia

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Matrix Equations (Ch20) · Positive Definite Matrices (Ch16)

**Chapter Outline**: Physical Motivation for Stability → Hurwitz Stability (Continuous Systems) → Schur Stability (Discrete Systems) → Lyapunov Stability Theorem (Positive Definite Criteria) → Definition of Matrix Inertia → Sylvester's Inertia Theorem Extensions → Routh-Hurwitz Criterion → D-Stability and P-Stability → Applications: Control Systems and Population Dynamics

**Extension**: Matrix stability is the algebraic criterion for determining whether the equilibrium points of a dynamical system (from weather models to economic cycles) possess "self-recovery" capabilities; it is the soul of Control Theory (Ch66).

</div>

If a physical system returns to its equilibrium state after a small perturbation, we call the system stable. In linear models, this physical property is completely encoded into the spectral distribution of the coefficient matrix. This chapter establishes the algebraic criteria for determining matrix stability and introduces "Inertia," a profound topological invariant that describes the distribution of the spectrum in the complex plane.

---

## 36.1 Hurwitz and Schur Stability

!!! definition "Definition 36.1 (Hurwitz Stability)"
    A square matrix $A \in M_n(\mathbb{C})$ is **Hurwitz Stable** (or asymptotically stable) if the real parts of all its eigenvalues are negative:
    $$\operatorname{Re}(\lambda_i) < 0, \quad \forall i$$
    **Physical Meaning**: The solution to the continuous system $\dot{\mathbf{x}} = A\mathbf{x}$ tends to zero as $t \to \infty$.

!!! definition "Definition 36.2 (Schur Stability)"
    A square matrix $A$ is **Schur Stable** if the moduli of all its eigenvalues are less than 1:
    $$|\lambda_i| < 1, \quad \forall i$$
    **Physical Meaning**: The solution to the discrete system $\mathbf{x}_{k+1} = A\mathbf{x}_k$ tends to zero as $k \to \infty$.

---

## 36.2 Lyapunov Stability Theorem

!!! theorem "Theorem 36.1 (Lyapunov Criterion)"
    A square matrix $A$ is Hurwitz stable if and only if for every positive definite matrix $Q \succ 0$, the following **Lyapunov Equation** has a unique positive definite solution $P \succ 0$:
    $$A^T P + PA = -Q$$
    **Significance**: This result transforms the global information of "spectral distribution" into the algebraic computation of a "matrix equation," avoiding the need to explicitly solve for eigenvalues.

---

## 36.3 Matrix Inertia

!!! definition "Definition 36.3 (Matrix Inertia)"
    The **Inertia** of a matrix $A$ is a triple $\operatorname{In}(A) = (\pi, \nu, \delta)$:
    - $\pi$: The number of eigenvalues with a positive real part.
    - $\nu$: The number of eigenvalues with a negative real part.
    - $\delta$: The number of eigenvalues with a zero real part.
    **Property**: $A$ is Hurwitz stable $\iff \operatorname{In}(A) = (0, n, 0)$.

---

## 36.4 Routh-Hurwitz Criterion

!!! technique "Technique: Routh Table"
    For a given characteristic polynomial $p(\lambda) = \sum a_k \lambda^k$, one can determine the number of eigenvalues with positive real parts without solving the equation by constructing a Routh table and checking the number of sign changes in its first column.

---

## Exercises

1.  **[Hurwitz] Determine if $A = \begin{pmatrix} -1 & 10 \\ 0 & -2 \end{pmatrix}$ is Hurwitz stable.**
    ??? success "Solution"
        Yes. The eigenvalues are -1 and -2, both of which have negative real parts.

2.  **[Schur] Determine if $\begin{pmatrix} 0.5 & 0.5 \\ 0 & 0.5 \end{pmatrix}$ is Schur stable.**
    ??? success "Solution"
        Yes. The eigenvalues are 0.5, both of which have moduli less than 1.

3.  **[Lyapunov] If the solution to $A^T P + PA = -I$ is $P = \operatorname{diag}(1, 2)$, is $A$ stable?**
    ??? success "Solution"
        Yes. Since $Q=I \succ 0$ and $P \succ 0$, by the Lyapunov theorem, $A$ is Hurwitz stable.

4.  **[Inertia] What is the inertia of the identity matrix $I_n$?**
    ??? success "Solution"
        $\operatorname{In}(I_n) = (n, 0, 0)$. All eigenvalues are 1.

5.  **[Trace] Prove: If $A$ is Hurwitz stable, then $\operatorname{tr}(A) < 0$.**
    ??? success "Solution"
        $\operatorname{tr}(A) = \sum \lambda_i$. Since each $\operatorname{Re}(\lambda_i) < 0$, the real part of their sum must also be negative. For a real matrix, the trace is real, so it must be less than 0.

6.  **[Skew-symmetric] Prove that a purely skew-symmetric matrix ($A^T = -A$) cannot be Hurwitz stable.**
    ??? success "Solution"
        The eigenvalues of a skew-symmetric matrix are purely imaginary (real part is 0), which does not satisfy the condition that the real part must be strictly less than 0. Its inertia is $(0, 0, n)$.

7.  **[Determinant] Prove: If an $n \times n$ real matrix $A$ is Hurwitz stable, then $(-1)^n \det(A) > 0$.**
    ??? success "Solution"
        $\det(A) = \prod \lambda_i$. Each real eigenvalue is negative, and complex eigenvalues appear in pairs with a positive product. Thus the sign is determined by the number of real negative roots, which is $(-1)^n$.

8.  **[D-Stability] What is D-stability?**
    ??? success "Solution"
        A matrix $A$ is D-stable if $DA$ is Hurwitz stable for every positive diagonal matrix $D$. This is crucial in the stability analysis of ecology and neural networks.

9.  **[Jury] What type of stability is the Jury criterion used to determine?**
    ??? success "Solution"
        It is used to determine the Schur stability of discrete-time systems (eigenvalues within the unit circle).

10. **[Application] In financial modeling, why does an eigenvalue near the imaginary axis signify risk?**

   ??? success "Solution"
        An eigenvalue with a real part near 0 means the system lacks damping power; perturbations will persist for a long time or even trigger resonance (Hopf bifurcation), causing the system to lose control.

## Chapter Summary

Matrix stability is the core criterion for analyzing dynamical systems:

1.  **Topological Classification of Spectra**: Hurwitz and Schur stability define the algebraic boundaries where "order" triumphs over "chaos" under continuous and discrete evolution, respectively.
2.  **Monotonicity of Energy**: Lyapunov's theorem proves that stability is essentially the monotonic decrease of a generalized energy (quadratic form) over time, transforming problems of analysis into problems of positive definite matrix algebra.
3.  **Structural Robustness**: Through tools like inertia and Routh tables, we can predict a system's survival space under parameter fluctuations without solving for exact eigenvalues, establishing a benchmark for robustness in control design.
