# Chapter 26: Linear Algebra in Differential Equations

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues and Diagonalization (Ch06) · Jordan Canonical Form (Ch12) · Matrix Functions (Ch13)

**Chapter Outline**: First-order Linear ODE Systems → The Matrix Exponential $e^{At}$ as a Fundamental Solution → Stability Analysis (Hurwitz Criteria) → Phase Plane Geometry (Classification of Singularities) → Converting High-order ODEs to Systems (Companion Matrices) → Non-homogeneous Systems and Duhamel's Principle → Spectral Methods for PDEs (Sturm-Liouville Theory) → Periodic Systems & Floquet Theory

**Extension**: Linear differential equations are the classical application of linear algebra; the long-term behavior of a dynamic system (stability, oscillation, decay) is entirely dictated by the spectrum of its coefficient matrix.

</div>

The study of change is the domain of differential equations, and the study of structured change is the domain of linear algebra. Linear ordinary differential equations (ODEs) can be completely solved using the matrix exponential $e^{At}$, and their stability—whether they converge to zero or explode—is determined by the real parts of their eigenvalues. This chapter builds the bridge between discrete operator theory and continuous dynamical evolution.

---

## 26.1 Linear Systems of ODEs

<div class="context-flow" markdown>

**The Core Structure**: The solution space of $\mathbf{x}' = A\mathbf{x}$ is an $n$-dimensional vector space. The fundamental matrix $\Phi(t)$ provides a basis for this space.

</div>

!!! definition "Definition 26.1 (Linear Constant-Coefficient ODE System)"
    A system of $n$ first-order linear ODEs is written as:
    $$\mathbf{x}'(t) = A\mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0$$
    where $A$ is an $n \times n$ constant matrix and $\mathbf{x}(t) \in \mathbb{R}^n$.

!!! theorem "Theorem 26.1 (The Fundamental Solution)"
    The unique solution to the system $\mathbf{x}' = A\mathbf{x}$ with initial condition $\mathbf{x}(0) = \mathbf{x}_0$ is:
    $$\mathbf{x}(t) = e^{At} \mathbf{x}_0$$
    where $e^{At}$ is the **matrix exponential** (see Ch13).

---

## 26.2 Stability Analysis

<div class="context-flow" markdown>

**Crucial Insight**: A system is stable if its state does not grow without bound over time. This depends solely on the eigenvalues of $A$.

</div>

!!! definition "Definition 26.2 (Stability Classifications)"
    1.  **Asymptotically Stable**: $\mathbf{x}(t) \to \mathbf{0}$ as $t \to \infty$ for all $\mathbf{x}_0$.
    2.  **Stable (Lyapunov Stable)**: $\|\mathbf{x}(t)\|$ remains bounded for all $t \ge 0$.
    3.  **Unstable**: $\|\mathbf{x}(t)\| \to \infty$ for some $\mathbf{x}_0$.

!!! theorem "Theorem 26.2 (Eigenvalue Stability Criterion)"
    Let $\sigma(A) = \{\lambda_1, \ldots, \lambda_n\}$ be the eigenvalues of $A$.
    - The system is **asymptotically stable** iff $\operatorname{Re}(\lambda_i) < 0$ for all $i$ ($A$ is a **Hurwitz matrix**).
    - The system is **stable** if $\operatorname{Re}(\lambda_i) \le 0$ for all $i$, and for any $\lambda_i$ with $\operatorname{Re}(\lambda_i) = 0$, its algebraic multiplicity equals its geometric multiplicity.
    - The system is **unstable** if any $\operatorname{Re}(\lambda_i) > 0$.

---

## 26.3 Phase Plane Analysis (2D Systems)

!!! technique "Classification of 2D Equilibrium Points"
    For a $2 \times 2$ system, the origin is classified based on $\lambda_1, \lambda_2$:
    - **Node**: Real eigenvalues of the same sign.
    - **Saddle**: Real eigenvalues of opposite signs.
    - **Spiral (Focus)**: Complex conjugate eigenvalues with $\operatorname{Re}(\lambda) \neq 0$.
    - **Center**: Purely imaginary eigenvalues ($\operatorname{Re}(\lambda) = 0$).

---

## 26.4 High-order ODEs and Companion Matrices

!!! technique "Reduction of Order"
    An $n$-th order scalar ODE $y^{(n)} + a_{n-1}y^{(n-1)} + \cdots + a_0 y = 0$ can be converted into a first-order system $\mathbf{x}' = C\mathbf{x}$ where $C$ is the **companion matrix**:
    $$C = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1} \end{pmatrix}$$
    The eigenvalues of $C$ are exactly the roots of the original characteristic equation.

---

## 26.5 Non-homogeneous Systems: Duhamel's Principle

!!! theorem "Theorem 26.3 (Variation of Parameters)"
    The solution to $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$ is:
    $$\mathbf{x}(t) = e^{At} \mathbf{x}_0 + \int_0^t e^{A(t-s)} \mathbf{f}(s) \, ds$$
    The integral term represents the "forced response" of the system.

---

## 26.6 Spectral Methods for PDEs (Sturm-Liouville)

!!! theorem "Theorem 26.4 (Sturm-Liouville Eigenfunctions)"
    Many linear PDEs (like the heat or wave equation) can be solved by separating variables, which leads to an eigenvalue problem for a differential operator. Sturm-Liouville theory guarantees that such operators have a complete set of orthogonal eigenfunctions, acting as a "Fourier basis" for the solution space.

---

## Exercises


****
??? success "Solution"
     Eigenvalues: $\lambda^2 + 3\lambda + 2 = 0 \implies \lambda_1 = -1, \lambda_2 = -2$.
     Eigenvectors: $\mathbf{v}_1 = (1, -1)^T, \mathbf{v}_2 = (1, -2)^T$.
     $\mathbf{x}(t) = c_1 e^{-t} \begin{pmatrix} 1 \\ -1 \end{pmatrix} + c_2 e^{-2t} \begin{pmatrix} 1 \\ -2 \end{pmatrix}$.
     $\mathbf{x}(0) = \begin{pmatrix} 1 \\ 0 \end{pmatrix} \implies c_1 = 2, c_2 = -1$.
     $\mathbf{x}(t) = \begin{pmatrix} 2e^{-t} - e^{-2t} \\ -2e^{-t} + 2e^{-2t} \end{pmatrix}$.


****
??? success "Solution"
     Yes, asymptotically stable. The eigenvalues are -1 and -2 (diagonal entries of a triangular matrix), both of which have negative real parts.


****
??? success "Solution"
     Eigenvalues: $\lambda^2 + 1 = 0 \implies \lambda = \pm i$. The origin is a **Center**, representing periodic oscillatory motion (closed orbits).


****
??? success "Solution"
     $\mathbf{x}(t) = e^{At} \mathbf{x}_0 + (e^{At} - I)A^{-1} \mathbf{b}$ (if $A$ is invertible).


****
??? success "Solution"
     $\begin{pmatrix} 0 & 1 \\ -6 & -5 \end{pmatrix}$.


****
??? success "Solution"
     $A^2 = O$, so $e^{At} = I + At = \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$.


****
??? success "Solution"
     If there exists a positive definite $P$ for a given $Q \succ 0$, then $A$ is Hurwitz (asymptotically stable). $V(x) = x^T P x$ acts as an energy function that strictly decreases over time.


****
??? success "Solution"
     The solution will contain terms like $t e^{\lambda t}$ in addition to $e^{\lambda t}$, representing slower-than-exponential decay or growth (or linear growth at resonance).


****
??? success "Solution"
     Because the Sturm-Liouville operator is **self-adjoint** under a specific inner product, and the spectral theorem for self-adjoint operators guarantees orthogonal eigenspaces.

****
??? success "Solution"
    ## Chapter Summary

Linear algebra provides the universal solution template for continuous change:


****: Identified $e^{At}$ as the fundamental operator that maps initial states to future states, unifying all linear ODE systems.

****: Linked the geometric position of eigenvalues in the complex plane to the qualitative physical behavior of systems (convergence vs. explosion).

****: Showed how high-order physical laws (like Newton's second law) can be linearized into first-order matrix systems.

****: Extended the theory to PDEs via Sturm-Liouville eigenfunctions, demonstrating that complex waves are just linear combinations of independent modal vibrations.
