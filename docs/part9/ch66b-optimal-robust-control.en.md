# Chapter 66B: Optimal and Robust Control

<div class="context-flow" markdown>

**Prerequisites**: State-Space Control (Ch66A) · Positive Definite Matrices (Ch16) · Matrix Equations (Ch20) · Matrix Inequalities (Ch18)

**Chapter Outline**: From Stabilization to Optimization → The Linear Quadratic Regulator (LQR) Model → Solving the Algebraic Riccati Equation (ARE) → Linear Quadratic Gaussian (LQG) Control and the Separation Principle → Motivation for Robustness (System Uncertainty) → Definition of $H_2$ and $H_\infty$ Norms → The Small Gain Theorem → Robust Control Design via LMIs → Applications: High-Performance Aircraft, Active Suspension, and Precision Manufacturing

**Extension**: Optimal control seeks the "least-cost" path of evolution, while robust control seeks the "most-resilient" defense strategy against disturbances; they elevate controller design from simple pole placement to multi-objective, matrix-norm optimization.

</div>

After achieving basic stabilization, the next engineering goal is often "optimality." We want the system to reach its target with minimum energy consumption or maximum speed. **Optimal Control** utilizes quadratic cost functions to provide standard algebraic solutions. Meanwhile, facing real-world fluctuations in model parameters, **Robust Control** ensures system safety by minimizing operator norms. This chapter explores these advanced algebraic topics in modern control theory.

---

## 66B.1 Linear Quadratic Regulator (LQR)

!!! definition "Definition 66B.1 (The LQR Problem)"
    Given a system $\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$, find the control sequence $\mathbf{u}(t)$ that minimizes the cost function $J$:
    $$J = \int_0^\infty (\mathbf{x}^T Q \mathbf{x} + \mathbf{u}^T R \mathbf{u}) dt$$
    where $Q \succeq 0$ penalizes state deviation and $R \succ 0$ penalizes energy consumption.

!!! theorem "Theorem 66B.1 (Optimal LQR Solution)"
    The optimal control law is a state feedback $\mathbf{u} = -K\mathbf{x}$, where:
    $$K = R^{-1} B^T P$$
    and $P$ is the unique positive definite solution to the **Algebraic Riccati Equation (ARE)**:
    $$A^T P + PA - P B R^{-1} B^T P + Q = 0$$

---

## 66B.2 Linear Quadratic Gaussian (LQG) and Separation

!!! technique "LQG Control"
    In the presence of measurement and process noise, LQG combines optimal estimation (the Kalman filter) with optimal control (LQR).
    **Separation Principle**: Proves that the gains for the optimal estimator and the optimal controller can be designed independently. This greatly simplifies the synthesis of complex systems.

---

## 66B.3 Robust Control and the $H_\infty$ Norm

!!! definition "Definition 66B.2 (The $H_\infty$ Norm)"
    For a system transfer function $G(s)$, the $H_\infty$ norm is the maximum singular value of its frequency response:
    $$\|G\|_\infty = \sup_{\omega} \sigma_{\max}(G(i\omega))$$
    It represents the maximum amplification factor of external disturbances.

!!! theorem "Theorem 66B.2 (Small Gain Theorem)"
    A closed-loop system with uncertainty feedback $\Delta$ is stable if the product of the loop gains satisfies:
    $$\|G\|_\infty \|\Delta\|_\infty < 1$$
    This provides a rigorous mathematical criterion for evaluating model simplification errors.

---

## 66B.4 Controller Design via LMIs

!!! technique "LMI-based Design"
    Modern control design often avoids explicit equation solving, instead framing stability, decay rate, and $H_\infty$ performance as a set of **Linear Matrix Inequalities (LMI)**.
    $$A^T P + PA + P B R^{-1} B^T P + Q \prec 0 \quad (\text{solved for } P \text{ via interior point methods})$$

---

## Exercises

1.  **[Basics] In the LQR cost function, does increasing the weight of matrix $R$ increase or decrease control energy?**
    ??? success "Solution"
        It decreases control energy. A larger $R$ penalizes the control input $u$ more heavily, causing the system to favor gentler actions.

2.  **[Riccati] Write the Riccati equation for the scalar system $\dot{x} = x + u$ with $Q=3, R=1$.**
    ??? success "Solution"
        $1P + P1 - P(1)(1)^{-1}(1)P + 3 = 0 \implies 2P - P^2 + 3 = 0$.
        Solving gives $P = 3$ (taking the positive root).

3.  **[Gain] Using the result from the previous problem, find the optimal feedback gain $K$.**
    ??? success "Solution"
        $K = R^{-1} B^T P = 1^{-1} \cdot 1 \cdot 3 = 3$. The optimal closed-loop system is $\dot{x} = (1-3)x = -2x$.

4.  **[LQG] In an LQG system, what determines the observer gain $L$?**
    ??? success "Solution"
        It is determined by the process noise covariance $W$ and the measurement noise covariance $V$ (by solving a dual Riccati equation, yielding the Kalman gain).

5.  **[Norm] Calculate the $H_\infty$ norm of the scalar system $G(s) = \frac{1}{s+1}$.**
    ??? success "Solution"
        $|G(i\omega)| = 1/\sqrt{\omega^2+1}$. The maximum occurs at $\omega=0$, so $\|G\|_\infty = 1$.

6.  **[Robustness] If a system has 20% parameter uncertainty, the Small Gain Theorem requires the $H_\infty$ norm of the nominal system to be less than what value?**
    ??? success "Solution"
        Less than $1/0.2 = 5$.

7.  **[Stability] Prove: if $Q$ is positive definite, the LQR closed-loop system is always asymptotically stable.**
    ??? success "Solution"
        The boundedness of $J$ and $Q \succ 0$ ensures that energy dissipates over time. Using the Lyapunov function $V(x) = x^T P x$, one can show $\dot{V} = -(x^T Q x + u^T R u) < 0$.

8.  **[Separation] In what scenarios does the separation principle fail?**
    ??? success "Solution"
        It often fails for non-linear systems or robust control systems with specific types of structured uncertainty where estimation and control become coupled.

9.  **[Implementation] Why prefer LMIs over Riccati equations in numerical implementation?**
    ??? success "Solution"
        LMIs can handle more flexible constraints (e.g., bounds on control gains, pole region restrictions) and benefit from the universality of convex optimization solvers.

10. **[Application] Briefly describe the role of $H_\infty$ control in wind-resistant drones.**

   ??? success "Solution"
        By minimizing the $H_\infty$ norm from wind disturbance to attitude deviation, we ensure that attitude fluctuations remain within safe limits even under strong gusts.

## Chapter Summary

Optimal and robust control theories define the performance limits of engineering control:

1.  **Algebraicized Trade-offs**: LQR transforms subjective human preferences for "performance" and "cost" into rigorous quadratic optimization problems through the configuration of $Q$ and $R$ matrices.
2.  **Rationality Under Noise**: LQG and the separation principle prove that optimal information extraction and optimal execution strategies are perfectly compatible under statistical uncertainty, establishing the logical pillars of feedback control.
3.  **Metrics of Uncertainty**: The $H_\infty$ norm and the small gain theorem provide quantitative geometric measures for system "safety," proving that the operator norm theory of linear algebra is the mathematical floodgate preventing real-world catastrophes.
