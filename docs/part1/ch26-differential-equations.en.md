# Chapter 26: Applications of Linear Algebra in Differential Equations

<div class="context-flow" markdown>

**Prerequisites**: Matrix Exponential (Ch13) · Eigenvalues (Ch6) · Stability (Ch36)

**Chapter Outline**: First-order Linear Systems $\dot{x} = Ax$ → Fundamental Solution Matrix → State Transition Matrix → Homogeneous and Non-homogeneous Solutions → Phase Portrait Analysis (Nodes, Saddles, Foci) → Reduction of Order for High-order Equations → Operator Methods → Introduction to PDE Discretization

**Extension**: Differential equations describe "change" in the world, while linear algebra provides the algebraic framework to capture the "steady states" and "modes" of those changes.

</div>

Linear systems of differential equations are at the heart of dynamical systems research. Through the matrix exponential, we transform complex calculus evolution into pure algebraic multiplication. The signs and magnitudes of eigenvalues directly determine whether a system settles into tranquility or drifts into chaos.

---

## 26.1 Core Solution Structures

!!! definition "Definition 26.1 (Fundamental Solution)"
    The general solution to the first-order homogeneous linear system $\dot{x}(t) = Ax(t)$ is $x(t) = e^{At}x(0)$.

!!! theorem "Theorem 26.3 (Stability Criterion)"
    The system $\dot{x} = Ax$ is asymptotically stable if and only if all eigenvalues of $A$ have strictly negative real parts.

---

## Exercises

1. **[Basic Calculation] Solve $\dot{x} = 2x, \dot{y} = 3y$ with initial values $(1, 1)$.**
   ??? success "Solution"
       The system matrix is $A = \begin{pmatrix} 2 & 0 \ 0 & 3 \end{pmatrix}$.
       The solution is $x(t) = e^{2t}, y(t) = e^{3t}$.

2. **[Phase Portrait] If $A = \begin{pmatrix} 1 & 0 \ 0 & -1 \end{pmatrix}$, describe the type of the equilibrium point $(0,0)$.**
   ??? success "Solution"
       The eigenvalues are $1$ and $-1$. One direction expands, while the other contracts. The equilibrium is a **Saddle Point**.

3. **[Reduction of Order] How do you convert the second-order equation $\ddot{y} + 3\dot{y} + 2y = 0$ into a first-order matrix equation?**
   ??? success "Solution"
       Let $x_1 = y$ and $x_2 = \dot{y}$.
       Then $\dot{x}_1 = x_2$ and $\dot{x}_2 = -2x_1 - 3x_2$.
       In matrix form: $\dot{\mathbf{x}} = \begin{pmatrix} 0 & 1 \ -2 & -3 \end{pmatrix} \mathbf{x}$.

4. **[Non-homogeneous] Write the Variation of Parameters formula for $\dot{x} = Ax + f(t)$.**
   ??? success "Solution"
       $x(t) = e^{At}x(0) + \int_0^t e^{A(t-	au)}f(	au)d	au$.

5. **[Stability Test] Is the system corresponding to matrix $\begin{pmatrix} -1 & 100 \ 0 & -2 \end{pmatrix}$ stable?**
   ??? success "Solution"
       Yes. The eigenvalues are $-1, -2$. Although there is a large off-diagonal term (which might cause short-term transient growth), the real parts are negative, so the system eventually converges to zero.

6. **[Solution Structure] If $A$ has a pair of conjugate complex eigenvalues $\pm i \beta$, what shape do the trajectories take?**
   ??? success "Solution"
       They form a **Center**, where trajectories are closed ellipses, representing undamped periodic oscillation.

7. **[Matrix Exponential] Calculate the solution matrix for $\dot{x} = \begin{pmatrix} 0 & 1 \ -1 & 0 \end{pmatrix} x$.**
   ??? success "Solution"
       The exponential matrix is the rotation matrix $e^{At} = \begin{pmatrix} \cos t & \sin t \ - \sin t & \cos t \end{pmatrix}$.

8. **[Initial Value Problem] If $x(0)$ is an eigenvector of $A$ for $\lambda$, prove $x(t) = e^{\lambda t} x(0)$.**
   ??? success "Solution"
       $x(t) = e^{At} x(0) = (I + At + \frac{1}{2}A^2t^2 + \dots) x(0)$.
       Since $A^k x(0) = \lambda^k x(0)$, substituting gives $x(t) = (1 + \lambda t + \frac{1}{2}\lambda^2 t^2 + \dots) x(0) = e^{\lambda t} x(0)$. This shows that eigenvectors define single-mode evolution.

9. **[Discretization] If using Euler discretization with step size $\Delta t$, what is the corresponding difference equation matrix?**
   ??? success "Solution"
       $x_{k+1} = x_k + \Delta t A x_k = (I + \Delta t A) x_k$.

10. **[Application] Why is "Modal Superposition" needed in studying large-scale vibrations?**
    ??? success "Solution"
        Modal superposition is essentially a coordinate transformation using the eigenvector matrix $P$, which decouples the coupled dynamical equations into a set of independent scalar second-order equations. This simplifies computation and allows us to focus only on dominant low-frequency modes.

## Chapter Summary

Differential equations are the dynamic laboratory of linear algebra:

1. **Modal Perspective**: Eigenvalues and eigenvectors define the fundamental frequencies and shapes of system evolution.
2. **Operator Bridge**: Matrix exponential transforms static structural parameters into dynamic time flows.
3. **Stability Bedrock**: Spectral theory provides the rigorous algebraic criteria for determining if a system will collapse or reach equilibrium.
