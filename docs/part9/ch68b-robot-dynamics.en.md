# Chapter 68B: Robot Dynamics

<div class="context-flow" markdown>

**Prerequisites**: Kinematics (Ch68A) · Matrix Differentiation (Ch47) · Eigenvalues (Ch6)

**Chapter Outline**: Lagrange Equations → Inertia Matrix → Coriolis and Centrifugal Matrix → Gravity Term → Linear Parameterization of Dynamic Equations → Least Squares Identification → Impedance Control → Principle of Virtual Work

**Extension**: The positive definiteness of the inertia matrix is a mathematical prerequisite for control law robustness; the linear parameterization property enables robot self-calibration through motion data.

</div>

Robot dynamics studies the causal relationship between joint torques and motion trajectories. Through Lagrangian mechanics or Newton-Euler methods, the laws of motion for multi-rigid-body systems are expressed as systems of non-linear differential equations parameterized by joint positions and velocities.

---

## 68B.1 Dynamic Equations and Algebraic Properties

!!! definition "Definition 68B.1 (Standard Second-Order Dynamic Model)"
    The equations of motion for an $n$-degree-of-freedom robot follow this matrix form:
    $$M(q) \ddot{q} + C(q, \dot{q}) \dot{q} + G(q) = 	au$$
    where:
    - $M(q) \in \mathbb{R}^{n 	imes n}$ is the **Inertia Matrix**;
    - $C(q, \dot{q}) \dot{q}$ represents the Coriolis and centrifugal terms;
    - $G(q)$ is the gravity vector;
    - $	au$ is the vector of actuator torques.

!!! theorem "Theorem 68B.1 (Positive Definiteness of the Inertia Matrix)"
    For any configuration $q$, the inertia matrix $M(q)$ is always **symmetric and positive definite** ($M^T = M, M \succ 0$). This stems from the physical positive definiteness of the system's total kinetic energy $T = \frac{1}{2} \dot{q}^T M(q) \dot{q}$.

---

## Exercises

1. **[Kinetic Energy Quadratic Form] Prove: If the total kinetic energy of the system is positive, then the inertia matrix $M(q)$ must have no zero or negative eigenvalues.**
   ??? success "Solution"
       Kinetic energy $T = \frac{1}{2}\dot{q}^T M \dot{q}$ is a quadratic form in $\dot{q}$. Since the kinetic energy of a physical system is strictly greater than 0 for any non-zero velocity, the quadratic form is positive definite. According to the matrix criterion for positive definite quadratic forms, $M$ must be a symmetric positive definite matrix.

2. **[Energy Conservation Property] Prove: The matrix $\dot{M}(q) - 2C(q, \dot{q})$ is skew-symmetric. What is the use of this property in control law design?**
   ??? success "Solution"
       Derived from energy conservation $\dot{T} = \dot{q}^T 	au$. Expanding the derivative terms and substituting the equations of motion proves $\dot{q}^T (\dot{M}-2C) \dot{q} = 0$. This skew-symmetry allows for canceling the derivative term of $M$ when constructing Lyapunov function terms like $\frac{1}{2}\dot{q}^T M \dot{q}$ in adaptive control.

3. **[Linear Parameterization] Prove that the dynamic equations are linear with respect to the set of physical parameters $\Phi$ (e.g., masses and moments of inertia of each link), i.e., $	au = Y(q, \dot{q}, \ddot{q}) \Phi$.**
   ??? success "Solution"
       Dynamic terms such as $m \ddot{x}$ or $I \dot{\omega}$ are linear multiples of the mass or inertia of each link. By extracting these parameters, one can construct a Regressor Matrix $Y$, thereby transforming the non-linear equation into a linear regression problem for parameter identification.

4. **[System Identification] Detail the algebraic process of identifying robotic physical parameters using the least squares method.**
   ??? success "Solution"
       Collect a series of data pairs $(	au_k, q_k, \dot{q}_k, \ddot{q}_k)$ and construct an overdetermined equation $\mathbf{T} = \mathbf{Y} \Phi$. Solve for $\Phi$ using the Moore-Penrose pseudoinverse: $\Phi = (\mathbf{Y}^T \mathbf{Y})^{-1} \mathbf{Y}^T \mathbf{T}$.

5. **[Impedance Control] Analyze the necessity of choosing the target inertia matrix $M_d$, damping matrix $B_d$, and stiffness matrix $K_d$ as positive definite matrices in impedance control.**
   ??? success "Solution"
       This ensures that the target closed-loop system $M_d \ddot{e} + B_d \dot{e} + K_d e = 0$ is asymptotically stable. Any negative or zero eigenvalue would lead to divergence or static bias during force interactions.

6. **[Calculation] Given a two-axis robotic arm, if its inertia matrix $M$ has an extremely large condition number $\kappa(M)$, what adverse effects does this have on trajectory tracking control?**
   ??? success "Solution"
       A large condition number means that in certain singular directions, the system exhibits extremely high sensitivity or massive inertial resistance. Numerical truncation errors in control commands will be amplified, leading to torque oscillations or a sharp drop in servo precision.

7. **[Modal Analysis] Analyze the mathematical connection between the natural frequencies of a linearized dynamic system and the eigenvalues of the matrix $M^{-1}K$.**
   ??? success "Solution"
       Near a linearized equilibrium point, $M \ddot{x} + K x = 0$. Assuming $x = v e^{j\omega t}$ yields $(\omega^2 M - K)v = 0$, or $(M^{-1}K - \omega^2 I)v = 0$. The squares of the natural frequencies correspond to the eigenvalues of $M^{-1}K$.

8. **[Projected Dynamics] Explain the projected representation of dynamic equations on the constraint subspace $N(J)$ under constrained motion.**
   ??? success "Solution"
       Utilize the projection matrix $P = I - J^\dagger J$ derived from the Jacobian matrix $J$. The dynamics of the constrained part can be decoupled from the full-degree-of-freedom equations using the projection matrix to analyze the separation of internal forces and motion.

9. **[Virtual Work Balance] Prove: The equivalent joint torques produced by an end-effector external force $f$ satisfy $	au_{ext} = J^T f$.**
   ??? success "Solution"
       From the principle of virtual work $\delta W = f^T \delta x = 	au^T \delta 	heta$. Using the mapping $\delta x = J \delta 	heta$ and substituting gives $f^T J \delta 	heta = 	au^T \delta 	heta$. Since $\delta 	heta$ is arbitrary, $	au = J^T f$.

10. **[Inverse Dynamics] Explain the algebraic reason why the Newton-Euler recursive algorithm is computationally superior to directly expanding the Lagrangian equations.**
    ??? success "Solution"
        The Newton-Euler algorithm exploits the local coupling of the chained structure, with a computational complexity of $O(n)$. Directly expanding the Lagrangian equations involves numerous partial derivative operations, generating redundant terms whose symbolic computation cost grows exponentially or as a high-order polynomial with $n$.

## Chapter Summary

This chapter discusses the dynamic framework describing the energy evolution and force characteristics of multi-degree-of-freedom mechanical systems:

1. **Matrix Equation Description**: Established standard matrix expressions composed of inertia, Coriolis, and gravity terms, unifying the evolution laws of rigid-body systems.
2. **Positive Definiteness Criteria**: Revealed physical constraints on the inertia matrix through kinetic energy quadratic forms, providing algebraic prerequisites for control stability analysis.
3. **Identification and Estimation**: Demonstrated the linear parameterization property of dynamic equations, establishing a linear algebraic path from motion data to physical models.
4. **Force Interaction**: Utilized the principle of virtual work to establish mapping rules for the Jacobian matrix in force space, forming the theoretical core of impedance control and force feedback.
