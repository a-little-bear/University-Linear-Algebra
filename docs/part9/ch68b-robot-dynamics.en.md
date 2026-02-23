# Chapter 68B: Robot Dynamics

<div class="context-flow" markdown>

**Prerequisites**: Kinematics (Ch68A) · Matrix Differentiation (Ch47) · Eigenvalues (Ch6)

**Chapter Outline**: Lagrange Equations → Inertia Matrix ($M$) → Coriolis and Centrifugal Matrix ($C$) → Gravity Term ($G$) → Linear Parameterization → Least Squares Identification → Impedance Control → Principle of Virtual Work

**Extension**: The positive definiteness of the inertia matrix is a mathematical prerequisite for control law robustness; the linear parameterization property enables robot self-calibration through data.

</div>

Robot dynamics studies the causal relationship between joint torques and motion trajectories. Through Lagrangian mechanics or Newton-Euler methods, the laws of motion for multi-rigid-body systems are expressed as systems of non-linear differential equations parameterized by joint positions and velocities.

---

## 68B.1 Dynamic Equations and Algebraic Properties

!!! definition "Definition 68B.1 (Standard Second-Order Model)"
    The equations of motion for an $n$-degree-of-freedom robot follow this matrix form:
    $$M(q) \ddot{q} + C(q, \dot{q}) \dot{q} + G(q) = \tau$$
    where:
    - $M(q) \in \mathbb{R}^{n \times n}$ is the **Inertia Matrix**;
    - $C(q, \dot{q}) \dot{q}$ represents the Coriolis and centrifugal terms;
    - $G(q)$ is the gravity vector;
    - $\tau$ is the vector of actuator torques.

!!! theorem "Theorem 68B.1 (Positive Definiteness of the Inertia Matrix)"
    For any configuration $q$, the inertia matrix $M(q)$ is always **symmetric and positive definite** ($M^T = M, M \succ 0$). This stems from the physical positive definiteness of the system's total kinetic energy $T = \frac{1}{2} \dot{q}^T M(q) \dot{q}$.

---

## Exercises

1. **[Kinetic Energy] Prove: If the total kinetic energy of the system is positive for any non-zero velocity, then the inertia matrix $M(q)$ must have only positive eigenvalues.**
   ??? success "Solution"
       Kinetic energy $T = \frac{1}{2}\dot{q}^T M(q) \dot{q}$ is a quadratic form. Since $T > 0$ for all $\dot{q} \neq 0$, the quadratic form is positive definite. According to the properties of real symmetric matrices, this implies all eigenvalues of $M(q)$ are strictly positive.

2. **[Energy Conservation] Prove: The matrix $\dot{M}(q) - 2C(q, \dot{q})$ is skew-symmetric. What is the implication for control design?**
   ??? success "Solution"
       From energy conservation $\dot{T} = \dot{q}^T \tau$. Differentiating $T$ and substituting the dynamic equation leads to $\dot{q}^T (\dot{M}-2C) \dot{q} = 0$. This skew-symmetry allows for the cancellation of $M$ derivative terms in Lyapunov-based adaptive control laws.

3. **[Linear Parameterization] Show that the dynamic equation is linear with respect to the physical parameter vector $\Phi$ (e.g., link masses, inertia tensors).**
   ??? success "Solution"
       Dynamic terms like $m \ddot{x}$ or $I \dot{\omega}$ are linear in the mass or inertia coefficients. By rearranging the equation into $\tau = Y(q, \dot{q}, \ddot{q}) \Phi$, where $Y$ is the regression matrix, we transform a non-linear state problem into a linear parameter estimation problem.

4. **[System Identification] Describe the algebraic process of identifying robotic physical parameters via least squares.**
   ??? success "Solution"
       Given experimental data pairs $(\tau_k, q_k, \dot{q}_k, \ddot{q}_k)$, we stack them into an overdetermined system $\mathcal{T} = \mathcal{Y} \Phi$. The optimal estimate is found using the pseudoinverse: $\hat{\Phi} = (\mathcal{Y}^T \mathcal{Y})^{-1} \mathcal{Y}^T \mathcal{T}$.

5. **[Impedance Control] Analyze the requirement for choosing target matrices $M_d, B_d, K_d$ as positive definite in impedance control.**
   ??? success "Solution"
       The target closed-loop system behaves like a virtual mass-spring-damper system. Positive definiteness ensures that the resulting linear system is asymptotically stable and dissipative, preventing energy gain during human-robot interaction.

6. **[Condition Number] What are the effects of a high condition number $\kappa(M)$ on robot servo performance?**
   ??? success "Solution"
       A large condition number indicates that the inertia is much larger in some directions than others. This amplifies numerical noise in the control loop and can lead to oscillations in the directions corresponding to small eigenvalues.

7. **[Modal Analysis] How do the natural frequencies of a robot relate to the eigenvalues of $M^{-1}K$ (where $K$ is the joint stiffness matrix)?**
   ??? success "Solution"
       The linearized homogeneous equation is $M \ddot{x} + K x = 0$. Assuming $x(t) = v e^{j\omega t}$ leads to the generalized eigenvalue problem $(K - \omega^2 M)v = 0$. The natural frequencies $\omega_i$ are the square roots of the eigenvalues of $M^{-1}K$.

8. **[Projected Dynamics] Explain the use of the projection matrix $P = I - J^\dagger J$ in constrained robotic motion.**
   ??? success "Solution"
       The matrix $P$ projects the dynamics onto the nullspace of the Jacobian. This allows for controlling internal forces (which do not cause motion) separately from external task-space motion.

9. **[Virtual Work] Use the Jacobian to derive the mapping from task-space forces $f$ to joint torques $\tau$.**
   ??? success "Solution"
       By the principle of virtual work, $f^T \delta x = \tau^T \delta \theta$. Substituting $\delta x = J \delta \theta$ gives $f^T J \delta \theta = \tau^T \delta \theta$. Since $\delta \theta$ is arbitrary, we have $\tau = J^T f$.

10. **[Inverse Dynamics] Why is the Newton-Euler recursive algorithm computationally more efficient than the Lagrangian expansion for large $n$?**
    ??? success "Solution"
        Lagrangian methods often result in $O(n^4)$ complexity if expanded symbolically due to redundant cross-terms. The recursive Newton-Euler algorithm exploits the local connectivity of the chain to achieve $O(n)$ complexity by propagating velocities and forces sequentially.

## Chapter Summary

This chapter discusses the dynamic framework describing the energy evolution and force characteristics of multi-degree-of-freedom mechanical systems:

1. **Matrix Equation Systems**: Established standard matrix expressions composed of inertia, Coriolis, and gravity terms.
2. **Definiteness Criteria**: Revealed physical constraints on the inertia matrix through quadratic forms, providing the basis for stability analysis.
3. **Identification Logic**: Demonstrated the linear parameterization of dynamics, establishing a linear algebraic path for parameter estimation.
4. **Interaction Control**: Utilized virtual work and projection operators to formulate impedance control and constrained motion laws.
