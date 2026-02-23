# Chapter 68A: Robot Kinematics

<div class="context-flow" markdown>

**Prerequisites**: Homogeneous Coordinates (Ch67) · Rotation Matrices (Ch05) · Geometric Algebra (Ch50)

**Chapter Outline**: From Degrees of Freedom to Pose Description → Representations of Rotation: Axis-angle, Euler Angles, and Quaternions → Core Description: Denavit-Hartenberg (DH) Parameterization → Forward Kinematics (FK): Concatenated Multiplication of Transformation Matrices → Algebraic Challenges of Inverse Kinematics (IK) → The Robot Jacobian Matrix ($J$) and Singular Poses → Differential Motion: From Velocity Space to Force Space → Applications: Industrial Manipulator Path Planning, Surgical Robots, and Humanoid Bipedal Locomotion

**Extension**: Robot kinematics is the "dynamic topology" of linear algebra; it transforms the physical constraints of mechanical structures into continuous motion on matrix manifolds. It proves that precise control of an end-effector is essentially the solving of a series of highly coupled non-linear matrix equations—the mathematical foundation for embodied intelligence.

</div>

In robotics, every joint movement corresponds to a coordinate transformation. The core task of **Robot Kinematics** is to establish a precise mathematical mapping between joint angles (internal variables) and the end-effector pose (external variables). By cascading **Homogeneous Transformation Matrices**, we can describe extremely complex mechanical chains like building blocks. This chapter introduces the motion geometry that powers industrial automation and precision manufacturing.

---

## 68A.1 Pose Description and DH Parameters

!!! definition "Definition 68A.1 (Pose Matrix)"
    The pose of a robot end-effector relative to the base frame is described by a $4 \times 4$ matrix $T$:
    $$T = \begin{pmatrix} R & \mathbf{p} \\ 0 & 1 \end{pmatrix}$$
    where $R$ represents orientation (rotation) and $\mathbf{p}$ represents spatial position.

!!! technique "Technique: DH Parameterization"
    Using four parameters (link length $a$, link twist $\alpha$, link offset $d$, and joint angle $\theta$), the transformation between adjacent joints is uniquely represented by a matrix $A_i$.

---

## 68A.2 Forward and Inverse Kinematics

!!! theorem "Theorem 68A.1 (Fundamental Kinematic Equations)"
    - **Forward Kinematics (FK)**: Uniquely determines the end-effector pose by computing the product $T = A_1 A_2 \cdots A_n$.
    - **Inverse Kinematics (IK)**: Given a target $T$, solve for the joint vector $\mathbf{q} = (\theta_1, \ldots, \theta_n)$.
    - **Property**: FK is a direct mapping, while IK is an inverse mapping that often involves multiple solutions or none at all.

---

## 68A.3 The Robot Jacobian and Velocity

!!! definition "Definition 68A.2 (Robot Jacobian $J$)"
    Establishes the linear mapping between joint velocities $\dot{\mathbf{q}}$ and end-effector velocities $\mathbf{v}$:
    $$\mathbf{v} = J(\mathbf{q}) \dot{\mathbf{q}}$$
    **Singular Pose**: A configuration where $\det(J) = 0$, causing the robot to lose mobility in certain directions.

---

## Exercises

**1. [Basics] Write the DH transformation matrix for a rotation $\alpha$ around the $x$-axis followed by a translation $a$ along the $x$-axis.**

??? success "Solution"
    **Construction:**
    According to DH rules, this is a composite transform: $A = \operatorname{Rot}(x, \alpha) \operatorname{Trans}(x, a)$.
    **Result**: $\begin{pmatrix} 1 & 0 & 0 & a \\ 0 & \cos\alpha & -\sin\alpha & 0 \\ 0 & \sin\alpha & \cos\alpha & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}$.

**2. [Calculation] For a 2-link planar arm with link lengths $L_1=L_2=1$ and angles $\theta_1=0, \theta_2=90^\circ$, find the end-effector coordinates.**

??? success "Solution"
    **Geometric Derivation:**
    1. First segment extends 1 unit along $x$.
    2. Second segment rotates $90^\circ$ relative to the first, extending 1 unit along $y$.
    3. Final coordinates: $x = 1+0=1, y=0+1=1$.
    **Matrix Verification**: $T = \operatorname{Rot}_z(0)\operatorname{Trans}_x(1) \cdot \operatorname{Rot}_z(90^\circ)\operatorname{Trans}_x(1)$ yields the same $(1, 1)$ result.

**3. [Inverse Kinematics] Explain why a 6-DOF manipulator typically has 8 inverse solutions.**

??? success "Solution"
    **Algebraic Context:**
    1. IK involves solving a system of coupled trigonometric equations.
    2. Each major configuration ("Shoulder," "Elbow," and "Wrist") has two symmetric possibilities (e.g., Left/Right, Elbow up/down).
    3. $2 \times 2 \times 2 = 8$.
    **Conclusion**: This reflects the global multi-valuedness of non-linear mappings.

**4. [Jacobian] What do the columns of the Jacobian matrix represent geometrically?**

??? success "Solution"
    **Conclusion:**
    The $i$-th column represents the velocity vector at the end-effector generated if only the $i$-th joint moves. For revolute joints, it is the cross product of the joint axis and the vector from the joint to the tip.

**5. [Singularity] Determine: If a manipulator is fully extended (straight line), is it in a singular pose?**

??? success "Solution"
    **Yes.**
    **Reasoning**: In a straight configuration, the arm cannot move any further outward radially (the velocity component in that direction is always 0). Algebraically, the Jacobian loses rank, and a zero appears in its singular value spectrum.

**6. [Application] Briefly state the role of "Pseudoinverse Control" in redundant robots.**

??? success "Solution"
    If a robot has >6 DOF (redundant), the equation $\mathbf{v} = J\dot{\mathbf{q}}$ has infinitely many solutions. Using the **Moore-Penrose Pseudoinverse** $J^+$, we get $\dot{\mathbf{q}} = J^+ \mathbf{v}$. This guarantees the end-effector task is met while minimizing joint energy consumption (the norm of $\dot{\mathbf{q}}$).

**7. [Quaternions] Why is Spherical Linear Interpolation (SLERP) preferred over linear interpolation of rotation matrices?**

??? success "Solution"
    **Reasoning:**
    Linear combinations of rotation matrices are generally not orthogonal (causing distortion). SLERP moves along a geodesic on the unit 4-sphere, ensuring constant angular velocity and maintaining valid rotation properties throughout the transition (see Ch51).

**8. [Calculation] A point $P^B$ is observed in frame $B$. If frame $B$ relative to $A$ is $T_B^A$, find $P^A$.**

??? success "Solution"
    **Formula:**
    $P^A = T_B^A P^B$ (using homogeneous coordinates).
    This demonstrates how linear algebra unifies observations across different reference frames.

**9. [Duality] What role does the transpose $J^T$ play in robotics?**

??? success "Solution"
    **Conclusion: Statics Mapping.**
    The formula $\boldsymbol{\mu} = J^T \mathbf{f}$ maps the external force $\mathbf{f}$ at the end-effector to the required joint torques $\boldsymbol{\mu}$. This embodies the duality between velocity space and force space (Principle of Virtual Work).

**10. [Application] Briefly state the meaning of the Image Jacobian in Visual Servoing.**

??? success "Solution"
    The Image Jacobian relates camera motion to the movement of feature points in the image plane. By inverting this matrix, a robot can automatically calculate motor corrections based on image errors to achieve closed-loop target tracking.

## Chapter Summary

Robot kinematics is the "embodiment" of linear algebra in physical space:

1.  **Chained Logic**: Through matrix cascading, complex mechanical topologies are simplified into continuous operator compositions, establishing a universal algebraic paradigm for multi-link systems.
2.  **Differential View of Velocity**: The Jacobian matrix proves that local linearization is the key to understanding global motion, providing a unified framework for velocity, torque, and singularity analysis.
3.  **Inverse Challenges**: The multi-valuedness and non-linearity of inverse mappings reveal the inherent difficulty of moving from "Task Space" back to "Joint Space," driving continuous breakthroughs in geometric algebra and numerical optimization.
