# Chapter 68B: Robot Dynamics and Perception

<div class="context-flow" markdown>

**Prerequisites**: Robot Kinematics (Ch68A) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16) · Probability and Statistics

**Chapter Outline**: From Geometry to Mechanics → Matrix Formulation of Lagrangian Dynamics → The Core Operator: The **Mass Matrix** $M(q)$ → Linear Representation of Centrifugal, Coriolis, and Gravity Terms → Forward Dynamics (Finding Acceleration) vs. Inverse Dynamics (Finding Torques) → Robot Perception: Sensor Calibration and Coordinate Transforms → Data Fusion: Introduction to Linear Kalman Filtering → Covariance Evolution in State Estimation → Applications: Dynamic Compensation Control, Collision Detection, and Mobile Robot Localization (SLAM)

**Extension**: Robot dynamics is the "energy expression" of linear algebra; it proves that complex mechanical behavior can be encapsulated as a positive-definite matrix operator evolving with configuration. Perception is the process of seeking truth through the propagation of matrix uncertainty—the physical core of building closed-loop intelligent systems.

</div>

If kinematics studies "where" a robot is, **Dynamics** studies "how" it moves. **Robot Dynamics** describes the causal relationship between joint torques and the resulting accelerations. By defining the **Mass Matrix** $M(q)$, we transform the distribution of mechanical energy into algebraic properties of a matrix. Simultaneously, **Robot Perception** uses matrix transformations and statistical filtering (such as the Kalman Filter) to turn noisy sensor data into precise estimates of the true state. This chapter introduces the algebraic system serving as the link between a robot's "brain" and its "muscles."

---

## 68B.1 Matrix Form of the Equations of Motion

!!! definition "Definition 68B.1 (Standard Dynamics Equation)"
    $$M(q) \ddot{q} + C(q, \dot{q})\dot{q} + g(q) = \tau$$
    - $M(q)$: **Inertia (Mass) Matrix** (Symmetric and Positive Definite).
    - $C(q, \dot{q})$: **Coriolis and Centrifugal Matrix**.
    - $g(q)$: **Gravity Vector**.
    - $\tau$: **Joint Torque Inputs**.

---

## 68B.2 Perception and the Kalman Filter

!!! technique "Technique: Linear State Estimation"
    In perception, system state updates follow linear equations: $x_k = Ax_{k-1} + Bu_k + w$. The **Kalman Filter** uses the evolution of the covariance matrix $P$ to weight the reliability of "prediction" versus "observation." The update step involves matrix inversion, which is essentially solving a weighted least-squares problem.

---

## 68B.3 Coordinate Transforms and Calibration

!!! definition "Definition 68B.2 (Sensor Extrinsics)"
    The transformation matrix $T_{sensor}^{body}$ from a sensor (like LiDAR) to the robot body is called the extrinsic calibration. It can be determined by solving matrix equations of the form $AX=XB$ (Hand-Eye Calibration).

---

## Exercises

**1. [Basics] Prove that the robot's inertia matrix $M(q)$ is always positive definite.**

??? success "Solution"
    **Physical Proof:**
    1. The kinetic energy of the system is $K = \frac{1}{2} \dot{q}^T M(q) \dot{q}$.
    2. According to physics, as long as there is motion ($\dot{q} \neq 0$), kinetic energy must be strictly positive.
    3. The quadratic form $v^T M v > 0$ is exactly the definition of a **Positive Definite Matrix** (Ch16).
    **Significance**: This ensures the numerical stability of the dynamics equations during integration.

**2. [Calculation] Let $M = \operatorname{diag}(2, 1)$. To generate acceleration $\ddot{q} = (1, 5)^T$, neglecting other terms, what torques are required?**

??? success "Solution"
    **Steps:**
    $\tau = M \ddot{q} = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 1 \\ 5 \end{pmatrix} = \begin{pmatrix} 2 \\ 5 \end{pmatrix}$.
    **Conclusion**: Joint 1 requires 2 units of torque, and Joint 2 requires 5 units.

**3. [Property] What is the "Skew-Symmetry Property" in robot dynamics?**

??? success "Solution"
    **Conclusion: $\dot{M} - 2C$ is a skew-symmetric matrix.**
    **Significance**: This property reflects **Conservation of Energy** in the physical system. It plays a core role in designing robust non-linear controllers (like adaptive control), allowing algebraic cancellation of terms to prove closed-loop stability.

**4. [Perception] Sensor data has variance $\sigma^2$. Where does this reside in a matrix?**

??? success "Solution"
    **Conclusion: The diagonal.**
    The diagonal entries $R_{ii}$ of the measurement noise covariance matrix $R$ represent the variance of the $i$-th sensor. Off-diagonal entries represent noise correlation between different sensors (usually assumed to be zero).

**5. [Kalman Filter] In the Kalman Gain $K = PH^T(HPH^T + R)^{-1}$, what is the physical meaning of the inversion?**

??? success "Solution"
    **Algebraic Interpretation:**
    The inversion corresponds to the **inverse of information** (precision). When measurement noise $R$ is massive, the inverse matrix becomes small, leading to a smaller gain $K$. This means the algorithm trusts its prediction more and ignores the unreliable observation data.

**6. [Calculation] Two sensors measure the same quantity with variances $P_1=1$ and $P_2=4$. What is the fused variance?**

??? success "Solution"
    **Formula:**
    $P_{combined} = (P_1^{-1} + P_2^{-1})^{-1} = (1 + 0.25)^{-1} = 1/1.25 = 0.8$.
    This corresponds to the **operator harmonic mean** (Ch46B). By fusing more information, the total uncertainty is reduced.

**7. [Calibration] Explain the meaning of the $AX=XB$ equation in robot calibration.**

??? success "Solution"
    This is the standard form for **Hand-Eye Calibration**.
    - $A$: Relative transform between two poses of the end-effector (known).
    - $B$: Relative transform between two observations of a target by the camera (known).
    - $X$: The unknown mounting pose of the camera relative to the end-effector.
    Solving this equation is essentially finding the unique similarity transformation that closes the two motion chains.

**8. [Dynamics] What is "Computed Torque" control?**

??? success "Solution"
    **Explanation:**
    High-speed computers calculate the required torque at each instant: $\tau_{calc} = M\ddot{q}_d + C\dot{q} + g$, and send this directly as a motor command. Linear algebra proves that through this algebraic compensation, a complex non-linear robot can be "linearized" into a set of simple second-order decoupled systems.

**9. [Mobile Robots] In SLAM, why are map point coordinates stored in the state vector?**

??? success "Solution"
    Because map point positions are uncertain. By merging landmark coordinates and robot pose into one massive state vector $\mathbf{X}$, the Kalman filter uses the **off-diagonal entries of the covariance matrix** to record the correlation between the "robot pose" and "landmark positions." When the robot sees a known landmark, this correlation helps correct accumulated errors in the entire historical trajectory.

**10. [Application] Briefly state the role of the "Mass Matrix" in collision detection.**

??? success "Solution"
    By comparing the actual torque $\tau_{act}$ inferred from current measurements with the theoretical $\tau_{calc}$, we compute the difference $\Delta \tau$. If the energy term $\Delta \tau^T M^{-1} \Delta \tau$ exceeds a threshold, an unexpected collision is detected. Here $M^{-1}$ acts as a converter from torque space to acceleration (kinetic energy) space.

## Chapter Summary

Robot dynamics and perception are the "closed-loop logic" of linear algebra in the physical world:

1.  **Algebraization of Mechanics**: The mass matrix $M(q)$ transforms complex mass distributions into positive-definite operators, establishing a rigorous causal chain between control torques and motion response.
2.  **Probabilistic Measure of Perception**: The Kalman filter, through the linear evolution of covariance matrices, proves that "truth" can be found by optimal weighted projection of residuals, providing mathematical criteria for handling real-world uncertainty.
3.  **Synthesis of Systems**: From geometric calibration via $AX=XB$ to dynamic compensation control, linear algebra achieves complete mathematical encapsulation from mechanical structure to perceptual information—the underlying operating system for modern industrial robots and autonomous vehicles.
