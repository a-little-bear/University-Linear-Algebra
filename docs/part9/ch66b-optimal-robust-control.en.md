# Chapter 66B: Optimal and Robust Control

<div class="context-flow" markdown>

**Prerequisites**: State-Space (Ch66A) · Positive Definite Matrices (Ch16) · Matrix Equations (Ch20) · Singular Values (Ch7)

**Chapter Outline**: Linear Quadratic Regulator (LQR) → Algebraic Riccati Equation (ARE) → State Estimation and Kalman Filtering → LQG and Separation Principle → $H_2$ and $H_\infty$ Norms → Robust Stability (Small Gain Theorem) → Linear Matrix Inequalities (LMIs) in Control

**Extension**: Optimal control theory is the core of aerospace, autonomous driving, and precision manufacturing; $H_\infty$ theory provides rigorous mathematical bounds for handling dynamical uncertainties.

</div>

Optimal and robust control study optimization problems of operators under constraints. LQR links dynamical stability with the solution of Algebraic Riccati Equations by minimizing variational functionals; robust control seeks hard guarantees for systems against perturbations within the sense of the Operator Norm.

---

## 66B.1 Linear Quadratic Regulator (LQR)

!!! definition "Definition 66B.1 (LQR Variational Problem)"
    Given a linear system $\dot{x} = Ax + Bu$, solve for the optimal control law $u(t)$ to minimize a quadratic performance index:
    $$J = \int_0^\infty (x^T Q x + u^T Ru) \, dt, \quad Q \succeq 0, R \succ 0$$
    The solution to this problem, derived via Hamiltonian system analysis, reduces to solving a nonlinear algebraic equation.

!!! theorem "Theorem 66B.1 (Algebraic Riccati Equation)"
    The optimal control gain for the LQR problem is $K = R^{-1} B^T P$, where $P$ is the unique symmetric positive definite solution to the following **Algebraic Riccati Equation (ARE)**:
    $$A^T P + PA - PBR^{-1}B^T P + Q = 0$$

---

## 66B.2 Robustness and the Small Gain Theorem

!!! theorem "Theorem 66B.3 (Small Gain Theorem)"
    Let a system $G$ and feedback perturbation $\Delta$ form a loop. The closed-loop system remains robustly stable if and only if the product of the $H_\infty$ norms of the loop components is less than 1:
    $$\|G(s)\|_\infty \cdot \|\Delta(s)\|_\infty < 1$$
    This reflects the stability constraints on the operator spectrum's position in the complex plane.

---

## Exercises

1. **[ARE Properties] Analyze the term $-PBR^{-1}B^TP$ in the Algebraic Riccati Equation and explain its role in maintaining closed-loop stability.**
   ??? success "Solution"
       This term represents energy dissipation introduced by the control input. In Lyapunov stability analysis, it provides a negative definite term that counteracts the energy growth brought by an unstable matrix $A$, ensuring closed-loop eigenvalues lie in the open left-half complex plane.

2. **[Hamiltonian Operator] Prove the relationship between the solution to the ARE and the spectral structure of the Hamiltonian matrix $H = \begin{pmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{pmatrix}$.**
   ??? success "Solution"
       The eigenvalues of a Hamiltonian matrix are symmetric with respect to the imaginary axis. The positive definite solution $P$ to the ARE is determined by the stable subspace of $H$ corresponding to the left-half complex plane. $P = X_2 X_1^{-1}$, where $[X_1; X_2]$ is the matrix formed by the stable eigenvectors of $H$.

3. **[Separation Principle] Prove that the transfer function matrix and pole distribution of an LQG control system satisfy the separation principle.**
   ??? success "Solution"
       The closed-loop eigenvalues are the union of $A-BK$ (controller poles) and $A-LC$ (observer poles). This is because the dynamics of the state estimation error $\tilde{x} = x - \hat{x}$ and the state evolution exhibit a block-triangular structure in the augmented space.

4. **[Calculation] Given a scalar system $\dot{x} = x + u$ with $Q=3, R=1$. Solve the ARE and determine the optimal feedback gain $K$.**
   ??? success "Solution"
       ARE: $1 \cdot P + P \cdot 1 - P \cdot 1 \cdot 1^{-1} \cdot 1 \cdot P + 3 = 0 \implies P^2 - 2P - 3 = 0$. The positive root is $P = 3$. The optimal gain is $K = 1^{-1} \cdot 1 \cdot 3 = 3$. Closed-loop dynamics: $\dot{x} = (1-3)x = -2x$.

5. **[H-infinity] Define the $H_\infty$ norm of a system's transfer matrix and explain its relationship to the maximum singular value curve $\bar{\sigma}(G(j\omega))$.**
   ??? success "Solution"
       $\|G(s)\|_\infty = \sup_{\omega \in \mathbb{R}} \bar{\sigma}(G(j\omega))$. It represents the maximum energy gain (peak gain) across all input frequencies and is the core metric for robust performance analysis.

6. **[Kalman Filtering] Describe how the optimal Kalman filter gain $L$ is determined algebraically and its consistency with the LQR problem in a dual sense.**
   ??? success "Solution"
       Kalman filtering is the dynamical version of least-squares estimation. Its optimal gain solution similarly reduces to an ARE, where the process noise covariance $W$ corresponds to $Q$ and the measurement noise covariance $V$ corresponds to $R$.

7. **[Robust Stability] Prove: If $\|G(s)\|_\infty < \gamma$, then for all stable perturbations satisfying $\|\Delta(s)\|_\infty \le 1/\gamma$, the closed-loop system remains stable.**
   ??? success "Solution"
       This is a direct application of the Small Gain Theorem. By the Nyquist stability criterion, the contraction of the loop gain within the unit circle ensures the characteristic locus does not encircle the critical point $-1$.

8. **[LMI] Explain the computational advantages of Linear Matrix Inequalities (LMIs) in multi-objective optimization for control systems.**
   ??? success "Solution"
       LMIs transform non-convex control constraints (e.g., simultaneous pole placement and $H_\infty$ metrics) into convex optimization problems over the PSD cone. This allows for finding global or near-optimal solutions efficiently using interior-point methods.

9. **[Definiteness] Prove: If $Q \succ 0$ and $(A, B)$ is stabilizable, then the symmetric positive definite solution $P$ to the ARE exists and is unique.**
   ??? success "Solution"
       Given stabilizability and observability (guaranteed by $Q$), the Hamiltonian matrix has no eigenvalues on the imaginary axis and possesses a unique stable invariant subspace, which induces a unique positive definite symmetric matrix $P$.

10. **[Control Cost] Analyze the convergence conditions for the LQR cost function $J$ over an infinite time horizon.**
    ??? success "Solution"
        The convergence of $J$ requires asymptotic stability of the closed-loop system. According to Lyapunov theory, if there exists $P \succ 0$ satisfying the ARE, then $x(t)^T P x(t)$ is a descending function, thereby guaranteeing the boundedness of the performance index.

## Chapter Summary

This chapter discusses the optimization and robustness criteria in linear control theory:

1. **Variational Optimization**: Established the algebraic balance between performance and control effort through quadratic cost functions and Algebraic Riccati Equations.
2. **State Reconstruction**: Formulated the algebraic framework for Kalman filtering, achieving optimal state estimation under stochastic noise.
3. **Frequency-Domain Robustness**: Utilized the $H_\infty$ norm and the Small Gain Theorem to quantify hard boundaries for systems against model perturbations.
4. **Unified Computational Platform**: Demonstrated Linear Matrix Inequalities (LMIs) as a general algebraic tool for solving complex constrained control problems.
