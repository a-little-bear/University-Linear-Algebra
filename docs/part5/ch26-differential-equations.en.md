# Chapter 26: Linear Algebra in Differential Equations

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Exponential (Ch13) · Linear Equations (Ch01)

**Chapter Outline**: From Single Equations to Systems → Matrix Representation of Linear Ordinary Differential Equations (ODEs) → Homogeneous Systems and Fundamental Solutions → The Central Role of the Matrix Exponential $e^{At}$ → Non-homogeneous Systems and the Variation of Parameters → Phase Portraits and Stability Analysis → Classification of Equilibrium Points: Sinks, Sources, Saddles, and Centers → Reduction of Higher-order Equations → Applications: Coupled Oscillators, Circuit Analysis, and Dynamical Modeling

**Extension**: Linear algebra provides the "geometric life" to differential equations; it decomposes continuous time-evolution into scaling and rotation along characteristic directions. It proves that the essence of complex dynamics is pre-determined by the spectral structure of the operator—a magnificent bridge between algebra and analysis.

</div>

In calculus, we learn to solve single equations like $\dot{x} = ax$. In reality, however, variables are often coupled (like two interacting objects). **Systems of Linear ODEs** deconstruct these couplings via matrix forms. By utilizing **Eigen-decomposition** and the **Matrix Exponential**, we transform complex dynamic evolutions into simple geometric projections. This chapter shows how linear algebra serves as the ultimate crystal ball for predicting the long-term behavior of systems.

---

## 26.1 Matrix Formulation of ODE Systems

!!! definition "Definition 26.1 (Linear ODE System)"
    A first-order linear system can be written as:
    $$\dot{\mathbf{x}}(t) = A\mathbf{x}(t)$$
    where $\mathbf{x}(t)$ is the state vector and $A$ is the coefficient matrix.

!!! theorem "Theorem 26.1 (Analytic Solution)"
    Given the initial condition $\mathbf{x}(0) = \mathbf{x}_0$, the unique solution is:
    $$\mathbf{x}(t) = e^{At} \mathbf{x}_0$$
    If $A$ is diagonalizable as $PDP^{-1}$, the solution simplifies to $\mathbf{x}(t) = P e^{Dt} P^{-1} \mathbf{x}_0$.

---

## 26.2 Stability and Phase Portraits

!!! technique "Geometry: Classification of Equilibria"
    The properties of the origin (equilibrium point) are determined by the signs of the eigenvalues $\lambda_i$:
    1.  **Sink**: Real parts are all negative. All trajectories converge to the origin.
    2.  **Source**: Real parts are all positive. Trajectories move away from the origin.
    3.  **Saddle**: Eigenvalues have mixed signs.
    4.  **Center**: Purely imaginary eigenvalues. Trajectories form closed loops.

---

## 26.3 Non-homogeneous Systems

!!! algorithm "Algorithm 26.1 (Variation of Parameters)"
    For $\dot{\mathbf{x}} = A\mathbf{x} + \mathbf{f}(t)$, the general solution is:
    $$\mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)} \mathbf{f}(s) ds$$
    This reflects the convolution of the external force with the system's impulse response.

---

## Exercises

**1. [Basics] Transform the second-order equation $\ddot{y} + 3\dot{y} + 2y = 0$ into a first-order matrix system.**

??? success "Solution"
    **Steps:**
    1. Define state variables: $x_1 = y, x_2 = \dot{y}$.
    2. Find the derivatives:
       - $\dot{x}_1 = x_2$.
       - $\dot{x}_2 = \ddot{y} = -2y - 3\dot{y} = -2x_1 - 3x_2$.
    **Conclusion**: The matrix form is $\begin{pmatrix} \dot{x}_1 \\ \dot{x}_2 \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix}$.

**2. [Calculation] Find the eigenvalues of the system above and determine its stability.**

??? success "Solution"
    **Steps:**
    1. Characteristic equation: $\lambda(\lambda+3) + 2 = \lambda^2 + 3\lambda + 2 = 0$.
    2. Factorize: $(\lambda+1)(\lambda+2) = 0$.
    3. Eigenvalues: $\lambda_1 = -1, \lambda_2 = -2$.
    **Conclusion**: Both eigenvalues are negative real numbers. The system is **asymptotically stable** (the origin is a sink).

**3. [Analytic Solution] Find $\mathbf{x}(t)$ if $A = \operatorname{diag}(1, -1)$.**

??? success "Solution"
    **Calculation:**
    1. Matrix exponential $e^{At} = \operatorname{diag}(e^t, e^{-t})$.
    2. Multiply by initial state: $\mathbf{x}(t) = \begin{pmatrix} e^t & 0 \\ 0 & e^{-t} \end{pmatrix} \begin{pmatrix} x_1(0) \\ x_2(0) \end{pmatrix}$.
    **Conclusion**: $x_1(t) = x_1(0)e^t$ and $x_2(t) = x_2(0)e^{-t}$.

**4. [Phase Portrait] If the eigenvalues are $\pm i\omega$, what is the shape of the trajectories in the phase plane?**

??? success "Solution"
    **Conclusion: Ellipses or Circles.**
    **Reasoning**: Purely imaginary eigenvalues imply solutions involving $\sin(\omega t)$ and $\cos(\omega t)$. Energy is conserved; the system neither dissipates nor explodes, but oscillates in a closed loop around the equilibrium.

**5. [Jordan Terms] If $A = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$, will terms like $te^{\lambda t}$ appear in the solution?**

??? success "Solution"
    **Yes.**
    **Reasoning**: $e^{At} = e^{\lambda t} \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$ (see Ch13). The $t$ term in the upper right causes a linear time growth factor. Physically, this corresponds to "critical damping" or resonance phenomena.

**6. [Non-homogeneous] What is the "Impulse Response" of a linear system?**

??? success "Solution"
    **Conclusion: It is the matrix $e^{At}$.**
    In the convolution integral, $e^{A(t-s)}$ describes how a unit impulse input at time $s$ affects the state at time $t$. It is the fundamental transfer function linking external forces to internal dynamics.

**7. [Stability Criterion] Prove: If $\operatorname{tr}(A) > 0$, the system is unstable.**

??? success "Solution"
    **Proof:**
    1. The trace equals the sum of eigenvalues: $\operatorname{tr}(A) = \sum \lambda_i$.
    2. If $\operatorname{tr}(A) > 0$, at least one eigenvalue must have a positive real part.
    3. The component in that direction will grow exponentially as $e^{\lambda t}$.
    **Conclusion**: The system diverges at the origin.

**8. [Calculation] Find the general solution to $\dot{x} = y, \dot{y} = x$.**

??? success "Solution"
    **Steps:**
    1. Matrix $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.
    2. Eigenvalues $\pm 1$ with vectors $(1, 1)^T$ and $(1, -1)^T$.
    3. General solution: $\mathbf{x}(t) = c_1 e^t \begin{pmatrix} 1 \\ 1 \end{pmatrix} + c_2 e^{-t} \begin{pmatrix} 1 \\ -1 \end{pmatrix}$. This is a classic **Saddle Point**.

**9. [Decoupling] Briefly state the significance of "Modal Coordinate Transformation."**

??? success "Solution"
    Using $\mathbf{x} = P\mathbf{z}$ transforms a coupled system $\dot{\mathbf{x}} = A\mathbf{x}$ into a diagonalized system $\dot{\mathbf{z}} = D\mathbf{z}$. This allows us to solve for each mode (eigen-component) independently, decomposing complex system evolution into separate univariate exponential growths or decays.

**10. [Application] How is linear algebra used in the linearization of epidemic models (like SIR)?**

??? success "Solution"
    **Explanation:**
    Non-linear equations can be approximated by a first-order Taylor expansion near an equilibrium point. The resulting **Jacobian matrix** determines the initial trend of the epidemic. If its dominant eigenvalue has a positive real part (corresponding to $R_0 > 1$), the outbreak explodes; otherwise, it naturally fizzles out.

## Chapter Summary

Linear algebra is the "internal decoder" for differential dynamics:

1.  **Algebraization of Dynamics**: The matrix exponential $e^{At}$ compresses continuous time-evolution into the action of a static operator, unifying the logic of discrete mappings and continuous flows.
2.  **Determinism of Morphology**: The algebraic properties of eigenvalues (real/imaginary parts) translate directly into the geometric destiny of the system (stable, divergent, oscillatory), establishing the standard framework for qualitative analysis.
3.  **Spectral Perspective**: Through eigenvector decomposition, complex coupled interactions are reduced to a superposition of independent modes, revealing the orderly algebraic hierarchy behind seemingly chaotic dynamical systems.
