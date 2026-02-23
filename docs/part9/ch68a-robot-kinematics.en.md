# Chapter 68A: Robot Kinematics

<div class="context-flow" markdown>

**Prerequisites**: Matrix Transformations (Ch67) · Matrix Groups (Ch55) · Generalized Inverses (Ch33) · Matrix Calculus (Ch47A)

**Chapter Outline**: Algebraic Abstraction of Linkages → Forward Kinematics (FK) & the D-H Parameter Method → Algebraic Challenges of Inverse Kinematics (IK) → The Robot Jacobian Matrix ($J$) → Velocity Transformations and Torque Propagation → Singularity Analysis ($\det(J)=0$) → Product of Exponentials (PoE) Formula → Control of Redundant Robots via Pseudoinverses → Applications: Industrial Manipulators, Humanoid Motion Planning, and Surgical Robotics

**Extension**: Robot kinematics is the chained composition of linear transformations along a physical skeleton; it maps the rotation/extension of joint spaces to the pose of end-effectors in Cartesian space, serving as the geometric heart of modern automation.

</div>

How can we make a robotic arm precisely grasp a cup? This is essentially a complex coordinate transformation problem. **Robot Kinematics** utilizes matrix theory to describe the geometric relationships between links and joints. By representing the transformation of each link as a 4x4 matrix, we can calculate the end-effector's pose through simple matrix concatenation. This chapter explores this algebraic mapping from joint angles to spatial coordinates.

---

## 68A.1 Forward Kinematics and D-H Parameters

!!! definition "Definition 68A.1 (Forward Kinematics - FK)"
    FK computes the mapping from joint variable vector $\theta$ to the end-effector pose $T$: $T = f(\theta)$.

!!! technique "The Denavit-Hartenberg (D-H) Method"
    A standard way to describe transformations between adjacent joints using four parameters: link length $a$, link twist $\alpha$, link offset $d$, and joint angle $\theta$.
    $$A_i = \operatorname{Rot}_z(\theta_i) \operatorname{Trans}_z(d_i) \operatorname{Trans}_x(a_i) \operatorname{Rot}_x(\alpha_i)$$
    The total transformation is $T_n^0 = A_1 A_2 \cdots A_n$.

---

## 68A.2 The Robot Jacobian Matrix

!!! definition "Definition 68A.2 (Jacobian Matrix)"
    The Jacobian matrix $J(\theta)$ describes the linear relationship between joint velocities $\dot{\theta}$ and the linear/angular velocity $v$ of the end-effector:
    $$\mathbf{v} = J(\theta) \dot{\boldsymbol{\theta}}$$
    **Significance**: It is the differential form of kinematics, revealing the system's local flexibility in its current pose.

---

## 68A.3 Singularity Analysis

!!! theorem "Theorem 68A.1 (Kinematic Singularity)"
    When $\det(J(\theta)) = 0$, the robot is in a **singular configuration**.
    - In this pose, the robot loses the ability to move in certain directions (rank deficiency).
    - Attempting to move in these directions would result in joint velocities tending to infinity, potentially causing mechanical failure.

---

## 68A.4 Redundancy and Pseudoinverse Control

!!! technique "Algebraic Solution to IK"
    For a target velocity $v$, the minimum-norm joint velocity is given by the **Moore-Penrose Pseudoinverse**:
    $$\dot{\theta} = J^+ v + (I - J^+ J)z$$
    where $(I - J^+ J)z$ represents "null-space motion," which adjusts the joint configuration without affecting the end-effector's pose.

---

## Exercises

1.  **[Basics] For a 2D arm with two revolute joints, find the relation between end-effector coordinates $(x, y)$ and joint angles $(\theta_1, \theta_2)$.**
    ??? success "Solution"
        $x = l_1 \cos\theta_1 + l_2 \cos(\theta_1+\theta_2)$
        $y = l_1 \sin\theta_1 + l_2 \sin(\theta_1+\theta_2)$

2.  **[Jacobian] Calculate the Jacobian matrix for the 2D arm above.**
    ??? success "Solution"
        $J = \begin{pmatrix} \frac{\partial x}{\partial \theta_1} & \frac{\partial x}{\partial \theta_2} \\ \frac{\partial y}{\partial \theta_1} & \frac{\partial y}{\partial \theta_2} \end{pmatrix} = \begin{pmatrix} -y & -l_2 \sin(\theta_1+\theta_2) \\ x & l_2 \cos(\theta_1+\theta_2) \end{pmatrix}$.

3.  **[Singularity] For the 2D arm, what happens when $l_2 \sin\theta_2 = 0$?**
    ??? success "Solution"
        Then $\det(J) = 0$. The robot is singular (fully extended or folded), and it cannot produce any radial velocity component along the arm.

4.  **[D-H Parameters] Why does the D-H method require only 4 parameters per link instead of 6?**
    ??? success "Solution"
        Because D-H frame placement follows specific constraints (the $x$-axis intersects and is perpendicular to the previous $z$-axis), eliminating two degrees of freedom and simplifying the model.

5.  **[IK Challenge] Why is inverse kinematics typically harder than forward kinematics?**
    ??? success "Solution"
        FK is a unique matrix product, whereas IK involves solving non-linear equations, which may have multiple solutions (configurations), one solution, or no solution (if the target is outside the workspace).

6.  **[Workspace] What is the "reachable workspace" of a robot?**
    ??? success "Solution"
        The set of all points the end-effector can reach. In linear algebra terms, it is the image of the joint space manifold under the operator $f$.

7.  **[Statics] Relate joint torques $\tau$ to end-effector forces $F$.**
    ??? success "Solution"
        $\tau = J^T F$. This demonstrates the importance of the Jacobian transpose via the Principle of Virtual Work.

8.  **[Redundancy] If a 7-DOF arm operates in 3D space, what are the dimensions of its Jacobian?**
    ??? success "Solution"
        $6 \times 7$ (6 velocity components, 7 joint variables). Since there are more columns than rows, the system is redundant.

9.  **[PoE] What is the Product of Exponentials (PoE) formula?**
    ??? success "Solution"
        An FK representation based on Screw Theory: $T = e^{S_1 \theta_1} \cdots e^{S_n \theta_n} M$. It is more geometrically intuitive and easier to differentiate than D-H.

10. **[Application] Why are high-precision IK algorithms vital in surgical robots?**

   ??? success "Solution"
        Because surgery requires precise control in narrow spaces; small algebraic errors could lead to accidents. Redundant DOFs are used to avoid colliding with delicate tissues.

## Chapter Summary

Robot kinematics is the chained integration of linear algebra and rigid body geometry:

1.  **Accumulation of Transforms**: Proved that complex mechanical motion can be deconstructed into a cascade of local 4x4 matrices, establishing standard algebraic protocols for motion modeling.
2.  **Tangent Space of Velocity**: The Jacobian matrix maps non-linear spatial configurations to linear velocity vectors, providing precise algebraic criteria for analyzing system flexibility and singularity.
3.  **Freedom under Constraint**: Through null-space and pseudoinverse techniques, kinematics demonstrates how to optimize secondary goals (obstacle avoidance, energy saving) while satisfying the primary task (pose), showcasing linear algebra's superiority in handling redundant constraints.
