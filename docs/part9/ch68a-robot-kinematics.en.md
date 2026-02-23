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
    The end-effector pose $T_e$ of a robotic arm is formed by the sequential composition of joint transformation matrices $T_i(	heta_i)$:
    $$T_e(	heta) = T_1(	heta_1) T_2(	heta_2) \dots T_n(	heta_n) \in SE(3)$$
    where $	heta$ is a vector in the joint variable space (configuration space).

!!! theorem "Theorem 68A.3 (Geometric Jacobian Operator)"
    The Jacobian matrix $J(	heta)$ establishes a linear mapping from the joint velocity vector $\dot{	heta}$ to the end-effector spatial velocity (linear and angular) $v = [\omega; v]^T$:
    $$v = J(	heta) \dot{	heta}$$
    The nullspace and image of this matrix characterize the local motion constraints of the mechanism.

---

## Exercises

1. **[Transformation Cascade] Explain why forward kinematics modeling for serial robots can be entirely reduced to a matrix chain multiplication problem.**
   ??? success "Solution"
       Each link rotates or translates relative to the previous link through a joint. According to the properties of homogeneous transformation matrices, the composite transformation of the end-effector relative to the base equals the product of all intermediate relative transformation matrices.

2. **[D-H Parameters] Analyze the algebraic structure of each component matrix in the Denavit-Hartenberg parameterization and its constraints on coordinate frame selection.**
   ??? success "Solution"
       A D-H transformation is composed of a rotation/translation about the $z$-axis and a rotation/translation about the $x$-axis. It requires the common normal between adjacent $z$-axes to serve as the $x$-axis, thereby reducing a general 6-parameter transformation to 4 parameters and standardizing kinematic modeling.

3. **[Inverse Kinematics] Explain the non-uniqueness of inverse kinematics (IK) solutions and its algebraic correspondence to multi-solution phenomena in transcendental equations.**
   ??? success "Solution"
       Solving IK is equivalent to solving $T(	heta) = T_{target}$. Due to the presence of numerous $\sin$ and $\cos$ terms, this non-linear system typically yields multiple isolated solutions (e.g., "elbow-up" vs. "elbow-down" configurations) or no solution due to workspace boundaries.

4. **[Jacobian Derivation] Prove: For an $n$-axis manipulator, the $i$-th column of its Jacobian matrix corresponds to the velocity component at the end-effector produced by the $i$-th joint.**
   ??? success "Solution"
       From the total differential relation $v = \sum \frac{\partial T_e}{\partial 	heta_i} \dot{	heta_i}$. For a revolute joint, the $i$-th column is $[z_{i-1}; z_{i-1} 	imes (p_e - p_{i-1})]$, reflecting the geometric contribution of the instantaneous center of rotation.

5. **[Singularity Analysis] Analyze the algebraic significance of $\det(J J^T) = 0$ and its impact on control systems.**
   ??? success "Solution"
       A singularity point means $J$ is rank-deficient, and the system loses one or more degrees of freedom in the task space. At this point, the inverse kinematics velocity solution $\dot{	heta} = J^{-1} v$ can produce infinite joint velocities, leading to controller failure.

6. **[Calculation] Given a 2-DOF planar linkage (lengths $l_1, l_2$), derive its Jacobian matrix and identify all singular configurations.**
   ??? success "Solution"
       $J = \begin{pmatrix} -l_1 s_1 - l_2 s_{12} & -l_2 s_{12} \ l_1 c_1 + l_2 c_{12} & l_2 c_{12} \end{pmatrix}$.
       $\det(J) = l_1 l_2 s_2$. Singularities occur at $	heta_2 = 0$ or $\pi$ (the linkage is fully extended or folded).

7. **[Manipulability] Define the Yoshikawa manipulability index $w = \sqrt{\det(J J^T)}$ and explain its role in trajectory planning.**
   ??? success "Solution"
       This index measures the "volume" of motion capability of the end-effector in all directions. In trajectory planning, one should avoid regions where $w 	o 0$ to ensure sufficient motion margin against external disturbances.

8. **[Principle of Virtual Work] Use linear algebra to prove that joint torques $	au$ and end-effector forces $f$ satisfy the equilibrium equation $	au = J^T f$.**
   ??? success "Solution"
       From power conservation (virtual work equality): $f^T v = 	au^T \dot{	heta}$. Substituting $v = J \dot{	heta}$ yields $f^T J \dot{	heta} = 	au^T \dot{	heta}$. Since $\dot{	heta}$ is arbitrary, $	au^T = f^T J \implies 	au = J^T f$.

9. **[Redundancy Resolution] For redundant robots ($n > 6$), explain how to utilize the generalized inverse of $J$ and nullspace projection for task priority control.**
   ??? success "Solution"
       The general solution is $\dot{	heta} = J^\dagger v + (I - J^\dagger J) z$. Here $z$ is an arbitrary vector, and the projection term $(I - J^\dagger J) z$ belongs to the nullspace of $J$, not affecting the end-effector motion. It can be used to satisfy secondary goals like obstacle avoidance.

10. **[Screw Theory] Briefly describe the advantages of the Product of Exponentials (PoE) formula over traditional D-H methods in handling differential motion.**
    ??? success "Solution"
        The PoE formula is based on the exponential mapping $e^{\hat{\xi}	heta}$ of Lie algebras. It treats the entire kinematic chain as a series of screw motions, avoiding ambiguities in local coordinate frame definitions and providing more global consistency in Jacobian calculations.

## Chapter Summary

This chapter discusses the algebraic systems describing the geometric evolution of multi-rigid-body robots:

1. **Configuration Mapping**: Established an analytical mapping model from joint space to task space through transformation matrix cascades.
2. **Differential Kinematics**: Introduced the Jacobian matrix as the core operator, revealing the velocity and torque conversion relations between configuration and task spaces.
3. **Performance Evaluation**: Formulated algebraic metrics such as singularity and manipulability, providing quantitative criteria for robot configuration optimization and control stability.
4. **Redundancy Handling**: Demonstrated the critical role of linear subspace theory in task planning for systems with redundant degrees of freedom.
