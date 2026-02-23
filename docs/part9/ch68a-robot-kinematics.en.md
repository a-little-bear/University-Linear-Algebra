# Chapter 68A: Robot Kinematics

<div class="context-flow" markdown>

**Prerequisites**: Homogeneous Coordinates (Ch67) · Rotation Representations (Ch67) · Matrix Differentiation (Ch47)

**Chapter Outline**: Joint Space and Task Space → Cascaded Transformation Matrices → Denavit-Hartenberg (D-H) Parameterization → Forward Kinematics → Inverse Kinematics → Jacobian Matrix → Velocity Mapping → Singularity Analysis → Manipulability Ellipsoid

**Extension**: The Jacobian matrix is the computational cornerstone for robot force control and visual servoing; manipulability metrics provide algebraic criteria for robot mechanism design optimization.

</div>

Robot kinematics studies the mapping between the configuration space of multi-rigid-body systems and the Cartesian space. This field utilizes the $SE(3)$ Lie group and its algebra to describe complex chained geometric constraints, converting mechanical motion laws into systems of non-linear algebraic equations of transformation matrices.

---

## 68A.1 Kinematic Modeling Theory

!!! definition "Definition 68A.1 (Forward Kinematics Cascade)"
    The end-effector pose $T_e$ of a robotic arm is formed by the sequential composition of joint transformation matrices $T_i(\theta_i)$:
    $$T_e(\theta) = T_1(\theta_1) T_2(\theta_2) \dots T_n(\theta_n) \in SE(3)$$
    where $\theta$ is a vector in the joint variable space (configuration space).

!!! theorem "Theorem 68A.3 (Geometric Jacobian Operator)"
    The Jacobian matrix $J(\theta)$ establishes a linear mapping from the joint velocity vector $\dot{\theta}$ to the end-effector spatial velocity (linear and angular) $v = [\omega; v]^T$:
    $$v = J(\theta) \dot{\theta}$$
    The nullspace and image of this matrix characterize the local motion constraints of the mechanism.

---

## Exercises

1. **[Transformation Cascade] Explain why forward kinematics modeling for serial robots can be entirely reduced to a matrix chain multiplication problem.**
   ??? success "Solution"
       Each link rotates or translates relative to the previous link through a joint. According to the properties of homogeneous transformation matrices, the composite transformation of the end-effector relative to the base equals the product of all intermediate relative transformation matrices.

2. **[D-H Parameters] Analyze the algebraic structure of each component matrix in the Denavit-Hartenberg parameterization.**
   ??? success "Solution"
       A D-H transformation is composed of a rotation/translation about the $z$-axis and a rotation/translation about the $x$-axis. It requires the common normal between adjacent $z$-axes to serve as the $x$-axis, reducing a general 6-parameter transformation to 4 parameters and standardizing kinematic modeling.

3. **[Inverse Kinematics] Explain the non-uniqueness of inverse kinematics (IK) solutions and its algebraic correspondence to multi-solution phenomena.**
   ??? success "Solution"
       Solving IK is equivalent to solving the matrix equation $T(\theta) = T_{target}$. Due to the periodic nature of $\sin$ and $\cos$, this non-linear system typically yields multiple isolated solutions (e.g., "elbow-up" vs. "elbow-down" configurations) or no solution if the target is outside the workspace.

4. **[Jacobian Derivation] Prove: For an $n$-axis manipulator, the $i$-th column of its Jacobian matrix corresponds to the velocity component at the end-effector produced by the $i$-th joint.**
   ??? success "Solution"
       From the total differential relation $v = \sum \frac{\partial T_e}{\partial \theta_i} \dot{\theta_i}$. For a revolute joint, the $i$-th column is $[z_{i-1}; z_{i-1} \times (p_e - p_{i-1})]$, reflecting the geometric contribution of that specific joint's motion to the overall end-effector velocity.

5. **[Singularity Analysis] Analyze the algebraic significance of $\det(J J^T) = 0$ and its impact on control systems.**
   ??? success "Solution"
       A singularity point means $J$ is rank-deficient, and the system loses one or more degrees of freedom in the task space. At this point, the inverse kinematics velocity solution $\dot{\theta} = J^{-1} v$ can produce infinite joint velocities, leading to physical actuator failure.

6. **[Calculation] Given 2-DOF planar linkage (lengths $l_1, l_2$), derive its Jacobian matrix and identify all singular configurations.**
   ??? success "Solution"
       $J = \begin{pmatrix} -l_1 s_1 - l_2 s_{12} & -l_2 s_{12} \\ l_1 c_1 + l_2 c_{12} & l_2 c_{12} \end{pmatrix}$.
       The determinant is $\det(J) = l_1 l_2 s_2$. Singularities occur when $\theta_2 = 0$ or $\pi$ (the linkage is fully extended or folded).

7. **[Manipulability] Define the Yoshikawa manipulability index $w = \sqrt{\det(J J^T)}$ and explain its role.**
   ??? success "Solution"
       This index measures the "volume" of motion capability of the end-effector in all directions. In trajectory planning, one should avoid regions where $w \to 0$ to ensure the robot has sufficient motion margin to compensate for disturbances.

8. **[Principle of Virtual Work] Use linear algebra to prove that joint torques $\tau$ and end-effector forces $f$ satisfy the equilibrium equation $\tau = J^T f$.**
   ??? success "Solution"
       From power conservation: $f^T v = \tau^T \dot{\theta}$. Substituting $v = J \dot{\theta}$ yields $f^T J \dot{\theta} = \tau^T \dot{\theta}$. Since $\dot{\theta}$ is arbitrary, $\tau^T = f^T J \implies \tau = J^T f$. This defines the dual mapping between force and motion spaces.

9. **[Redundancy Resolution] For redundant robots ($n > 6$), explain how to utilize the generalized inverse of $J$ and nullspace projection.**
   ??? success "Solution"
       The general solution is $\dot{\theta} = J^\dagger v + (I - J^\dagger J) z$. Here $z$ is an arbitrary vector, and the term $(I - J^\dagger J) z$ belongs to the nullspace of $J$. This "internal motion" does not change the end-effector pose and can be used for obstacle avoidance.

10. **[Screw Theory] Briefly describe the advantages of the Product of Exponentials (PoE) formula over traditional D-H methods.**
    ??? success "Solution"
        The PoE formula is based on the exponential mapping $e^{\hat{\xi}\theta}$ of Lie algebras. It treats the entire kinematic chain as a series of screw motions, avoiding ambiguities in local frame definitions and providing more global consistency in differential analysis.

## Chapter Summary

This chapter discusses the algebraic systems describing the geometric evolution of multi-rigid-body robots:

1. **Configuration Mapping**: Established an analytical mapping model from joint space to task space through transformation matrix cascades.
2. **Differential Kinematics**: Introduced the Jacobian matrix as the core operator, revealing the velocity and torque conversion relations between spaces.
3. **Performance Evaluation**: Formulated algebraic metrics such as singularity and manipulability, providing quantitative criteria for robot configuration optimization.
4. **Redundancy Handling**: Demonstrated the critical role of linear subspace theory in task planning for systems with redundant degrees of freedom.
