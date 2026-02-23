# Chapter 68B: Robot Dynamics

<div class="context-flow" markdown>

**Prerequisites**: Robot Kinematics (Ch68A) · Positive Definite Matrices (Ch16) · Matrix Calculus (Ch47A) · Matrix Equations (Ch20)

**Chapter Outline**: Matrix Form of the Robot Dynamic Equations → Properties of the Mass Matrix $M(q)$ (Symmetry & Positive Definiteness) → The Coriolis and Centrifugal Matrix $C(q, \dot{q})$ → The Gravity Vector $G(q)$ → The Fundamental Identity: Skew-symmetry of $\dot{M}-2C$ and Energy Conservation → Operational Space Dynamics → Recursive Newton-Euler Algorithm (RNEA) → Applications: Model Predictive Control (MPC), Force/Position Control, and Parameter Identification

**Extension**: Robot dynamics is "linear algebra with mass"; it extends static geometric transforms into second-order matrix differential equations involving inertia, damping, and potential energy, forming the theoretical basis for high-performance motion control.

</div>

Kinematics tells us *where* a robot is, but **Robot Dynamics** tells us *how much force* is needed to get there. The dynamic equations link joint accelerations to the torques produced by the motors. in high-dimensional state space, these equations manifest as a highly non-linear matrix system. This chapter exploits the structural properties of linear algebra (e.g., positive definiteness, skew-symmetry) to simplify this complex physical description.

---

## 68B.1 Matrix Form of the Dynamic Equations

!!! definition "Definition 68B.1 (Standard Dynamic Model)"
    The dynamic equation for a robot with $n$ joints is written as:
    $$M(q)\ddot{q} + C(q, \dot{q})\dot{q} + G(q) = \tau$$
    - $q, \dot{q}, \ddot{q}$: Joint positions, velocities, and accelerations.
    - $M(q)$: $n \times n$ **Inertia Matrix** (or Mass Matrix).
    - $C(q, \dot{q})$: $n \times n$ **Coriolis and Centrifugal Matrix**.
    - $G(q)$: $n \times 1$ **Gravity Vector**.
    - $\tau$: Vector of joint torques.

---

## 68B.2 Algebraic Properties of the Inertia Matrix

!!! theorem "Theorem 68B.1 (Positive Definiteness of the Mass Matrix)"
    The inertia matrix $M(q)$ is always **symmetric and positive definite** ($M \succ 0$).
    **Physical Meaning**: The kinetic energy $K = \frac{1}{2} \dot{q}^T M(q) \dot{q}$ is always positive for any non-zero velocity, and the matrix is invertible, ensuring a unique acceleration for any given torque: $\ddot{q} = M^{-1}(\tau - C\dot{q} - G)$.

---

## 68B.3 Energy Conservation and Skew-Symmetry

!!! theorem "Theorem 68B.2 (Skew-Symmetry Identity)"
    The matrix $N = \dot{M}(q) - 2C(q, \dot{q})$ is a **skew-symmetric matrix**, satisfying $x^T N x = 0$ for any vector $x$.
    **Significance**: This property reflects the cancellation of internal work within the system and is a core technique for proving the stability of controllers, especially in Lyapunov-based adaptive control.

---

## 68B.4 Linear Structure of Parameter Identification

!!! technique "Linear Parameterization"
    The dynamic equations are **linear** with respect to the physical parameters (mass, center of mass, inertia tensor):
    $$\tau = Y(q, \dot{q}, \ddot{q}) \boldsymbol{\phi}$$
    where $Y$ is the **Regressor Matrix** and $\boldsymbol{\phi}$ is the vector of parameters to be identified. This allows for the use of **Least Squares** (Ch07) to accurately reconstruct the robot's dynamic model.

---

## Exercises


****
??? success "Solution"
     Kinetic energy $K = \frac{1}{2} m (l\dot{q})^2 = \frac{1}{2} (ml^2) \dot{q}^2$. Thus $M(q) = ml^2$ (a scalar inertia matrix).


****
??? success "Solution"
     Because kinetic energy is a positive definite quadratic form. If $M$ had non-positive eigenvalues, it would mean that moving in certain directions would result in zero or negative kinetic energy, which contradicts physical laws of mass distribution.


****
??? success "Solution"
     In a Lyapunov stability proof, the derivative of the kinetic energy term $\frac{1}{2} \frac{d}{dt}(\dot{q}^T M \dot{q})$ generates a $\frac{1}{2} \dot{q}^T \dot{M} \dot{q}$ term. By combining this with the $C$ matrix term, parts of the expression cancel out, simplifying the stability analysis.


****
??? success "Solution"
     $\dot{M} = \begin{pmatrix} -b\sin(q_2)\dot{q}_2 & 0 \\ 0 & 0 \end{pmatrix}$. Note the use of the chain rule.


****
??? success "Solution"
     $\Lambda(x) = (J M^{-1} J^T)^{-1}$. It describes the effective mass perceived at the end-effector in Cartesian space.


****
??? success "Solution"
     It utilizes the inverse of the dynamic equations to algebraically cancel the non-linear terms $C\dot{q}$ and $G$, forcing the system to behave as a linear double integrator: $\ddot{q} = u$.


****
??? success "Solution"
     This follows directly from the definition of the conservative force term in the Euler-Lagrange equations.


****
??? success "Solution"
     No, because of the non-linear velocity-squared terms $\dot{q}$ (Coriolis forces) and the constant gravity term $G$. It is only approximately linear in quasi-static conditions or when velocity is negligible.


****
??? success "Solution"
     $n^3$ symbols. Each $C_{ij} = \sum c_{ijk} \dot{q}_k$.

****
??? success "Solution"
    ## Chapter Summary

Robot dynamics demonstrates the depth of linear algebra in handling non-linear physical laws:


****: Through the three operators $M, C$, and $G$, dynamics condenses complex Newtonian derivations into standard matrix equations, establishing a universal format for multi-body simulation.

****: The positive definiteness of the mass matrix and the skew-symmetry of $\dot{M}-2C$ are not just mathematical tricks but projections of the Law of Conservation of Energy, providing a theoretical moat for robust control design.

****: The parameter-linearization property proves that even extremely complex motions are governed by underlying physical constants that can be isolated via linear regression, showcasing linear algebra as a powerful tool for scientific analysis.
