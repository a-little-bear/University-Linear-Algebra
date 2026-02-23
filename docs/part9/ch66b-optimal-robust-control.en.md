# Chapter 66B: Optimal and Robust Control

<div class="context-flow" markdown>

**Prerequisites**: State-Space and System Realization (Ch66A) · Matrix Equations (Ch20) · Matrix Inequalities (Ch18) · Convex Optimization (Ch25)

**Chapter Outline**: From Feedback to Performance Optimality → Mathematical Model of the Linear Quadratic Regulator (LQR) → The Core Operator: Algebraic Riccati Equation (ARE) → Optimal Feedback Gain $K = R^{-1} B^T P$ → System Stability: Matrix Form of Lyapunov Stability Theorem → Motivation for Robustness → $H_\infty$ Control and the Small Gain Theorem → Modern Methods: Dominance of Linear Matrix Inequalities (LMI) in Control Synthesis → Applications: UAV Balance Control, Robust Disturbance Rejection in Chemical Processes, and Satellite Orbit Precision Optimization

**Extension**: Optimal control is the "pinnacle of performance" in linear algebra; it transforms controller design into finding roots of specific matrix equations. It proves that system stability and disturbance rejection can be guaranteed entirely by the spectral distribution of operators and the feasible regions of matrix inequalities—the mathematical core of modern precision industry.

</div>

After solving the problem of whether a system is controllable, engineers focus on how to control it "best." **Optimal Control** seeks a control law that minimizes a cost function (such as energy consumption or error accumulation). **Robust Control** further considers model uncertainties. This chapter demonstrates how matrix equations and inequalities in linear algebra serve as the ultimate blueprint for designing stable, efficient, and robust control systems.

---

## 66B.1 Linear Quadratic Regulator (LQR)

!!! definition "Definition 66B.1 (LQR Problem)"
    Given a system $\dot{x} = Ax + Bu$, minimize the quadratic cost function:
    $$J = \int_0^\infty (x^T Q x + u^T R u) dt$$
    where $Q \succeq 0$ penalizes state deviations and $R \succ 0$ penalizes control energy.

!!! theorem "Theorem 66B.1 (LQR Optimal Solution)"
    The optimal control law is the linear feedback $\mathbf{u}(t) = -K \mathbf{x}(t)$, where the gain matrix is:
    $$K = R^{-1} B^T P$$
    and $P$ is the unique positive definite solution to the **Algebraic Riccati Equation (ARE)**:
    $$A^T P + PA - P B R^{-1} B^T P + Q = 0$$

---

## 66B.2 $H_\infty$ Control and Robustness

!!! definition "Definition 66B.2 ($H_\infty$ Norm)"
    The $H_\infty$ norm of a system represents the maximum amplification the system can exert on an external disturbance.
    **Small Gain Theorem**: If the $H_\infty$ norm of a closed-loop system is less than 1, then the system is robustly stable against any non-linear perturbation with a norm less than 1.

---

## 66B.3 LMI and Control Synthesis

!!! technique "Technique: Linear Matrix Inequalities (LMI)"
    Modern control design often solves a set of LMIs rather than explicit equations. For example, determining stability is equivalent to finding $P \succ 0$ such that:
    $$A^T P + PA \prec 0$$
    Using convex optimization tools, we can satisfy multiple constraints simultaneously, such as stability, decay rates, and actuator saturation.

---

## Exercises

**1. [Basics] In an LQR problem, what happens to the controller's behavior if the $R$ matrix becomes very large?**

??? success "Solution"
    **Algebraic Analysis:**
    1. $R$ represents the penalty on control effort $u$.
    2. As $R \to \infty$, the term $u^T Ru$ in the cost function explodes unless $u \to 0$.
    **Conclusion**: The controller becomes very "conservative" (lazy), leading to slower response times as it prioritizes saving control energy over rapid state regulation.

**2. [Riccati] For a scalar system $A=1, B=1, Q=3, R=1$, solve the Riccati equation for $P$.**

??? success "Solution"
    **Steps:**
    1. Write the ARE: $1P + P(1) - P(1)(1^{-1})(1)P + 3 = 0$.
    2. Simplify: $2P - P^2 + 3 = 0$.
    3. Rearrange: $P^2 - 2P - 3 = 0$.
    4. Factor: $(P-3)(P+1) = 0$.
    **Conclusion**: The positive definite solution is $P = 3$.

**3. [Calculation] Using the result above, find the optimal gain $K$ and the closed-loop system matrix.**

??? success "Solution"
    **Calculation:**
    1. $K = R^{-1} B^T P = (1)^{-1}(1)(3) = 3$.
    2. Closed-loop: $\dot{x} = (A - BK)x = (1 - 1 \cdot 3)x = -2x$.
    **Conclusion**: The closed-loop system is stable since the eigenvalue is $-2 < 0$.

**4. [Stability] Prove: If there exists $P \succ 0$ such that $A^T P + PA \prec 0$, then $A$ is Hurwitz stable.**

??? success "Solution"
    **Using Lyapunov Functions:**
    1. Let $V(x) = x^T P x$. Since $P \succ 0$, $V(x)$ is a valid energy function.
    2. Time derivative: $\dot{V} = \dot{x}^T P x + x^T P \dot{x} = (Ax)^T P x + x^T P (Ax) = x^T(A^T P + PA)x$.
    3. Since $A^T P + PA \prec 0$, then $\dot{V} < 0$ for all non-zero $x$.
    **Conclusion**: Energy decreases monotonically and is bounded below, so the system must converge to the origin.

**5. [LMI] How is $A^T P + PA + Q \prec 0$ transformed into standard LMI form?**

??? success "Solution"
    Since the expression is linear in the unknown matrix $P$ (involving only addition and constant multiplication) and satisfies a semi-definite cone constraint, it is already a standard Linear Matrix Inequality. In solvers like YALMIP or CVX, it is entered directly as a constraint.

**6. [Application] How are weights in the $Q$ matrix typically set in autonomous driving?**

??? success "Solution"
    **Explanation:**
    Usually, $Q$ is diagonal. High weights are assigned to components like "deviation from lane center" to enforce strict tracking. Lower weights are assigned to components like "steering wheel jitter" to allow some flexibility. This algebraically prioritizes safety over smoothness.

**7. [Robustness] What is a system's "Uncertainty Set"?**

??? success "Solution"
    It is the range of possible values for the actual system matrix $A$, often represented as $A_{real} = A_{nom} + \Delta$ where $\|\Delta\| \le \epsilon$. Robust control aims to find a single $K$ that stabilizes the system for **all** $A$ in the set.

**8. [Stability] If a system's eigenvalues are $\{-1, -2\}$ and a perturbation $E = \begin{pmatrix} 0 & 0.1 \\ 0.1 & 0 \end{pmatrix}$ is applied, is the system necessarily still stable?**

??? success "Solution"
    **Determination:**
    Using eigenvalue perturbation theory (Ch15/42), since the eigenvalues are far from the imaginary axis and the perturbation is small, the shift in eigenvalues is not enough to cross the axis.
    **Conclusion**: Due to the "stability margin," the system remains stable under small perturbations.

**9. [Duality] How does LQG control combine LQR and Kalman filtering?**

??? success "Solution"
    **Separation Principle:**
    1. Use a Kalman filter to estimate the optimal state $\hat{x}$ from noisy data.
    2. Apply the LQR gain $K$ to the estimated state: $u = -K\hat{x}$.
    Linear algebra proves these two processes can be designed independently without affecting overall optimality.

**10. [Application] Briefly state the role of $H_\infty$ control in suppressing wind disturbance for UAVs.**

??? success "Solution"
    Wind is modeled as an external energy source with bounded $L_2$ norm. The $H_\infty$ controller designs a gain matrix that minimizes the maximum singular value (gain) of the transfer operator from "wind" to "pose error." This guarantees that even in strong winds, pose deviations are strictly limited within algebraic bounds.

## Chapter Summary

Optimal and Robust control represent the "supreme command" of linear algebra in complex systems:

1.  **Algebraic Performance**: Through Riccati equations, it transforms the philosophical question of "what is best" into the arithmetic task of finding specific positive definite matrices, enabling hardcore quantification of performance.
2.  **Geometric Bounds of Robustness**: Via norms and LMIs, it establishes a "resilience radius" for systems, proving that robustness is essentially the tolerance of spectral structures to perturbation spaces.
3.  **Synthesis of Control**: This theory unifies observation (Kalman), decision-making (LQR), and defense (H-infinity) within matrix analysis, supporting the modern industrial backbone from deep-space exploration to precision manufacturing.
