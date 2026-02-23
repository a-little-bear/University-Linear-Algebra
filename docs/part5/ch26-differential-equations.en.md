# Chapter 26: Linear Systems of Differential Equations

<div class="context-flow" markdown>

**Prerequisites**: Matrix Exponential (Ch13) · Jordan Form (Ch12) · Eigenvalues (Ch6) · Diagonalization (Ch6)

**Chapter Outline**: First-order Linear Systems $\dot{x} = Ax$ → Homogeneous and Non-homogeneous Systems → Solution via Matrix Exponential → Decoupling via Diagonalization → Stability Analysis of Equilibrium Points → Phase Portraits → Variation of Parameters → Higher-order Linear ODEs as First-order Systems → Fundamental Matrix Solution

**Extension**: Systems of differential equations describe the time-evolution of almost all physical phenomena; linear algebra provides the exact solution for the "small oscillation" limit of these systems.

</div>

Systems of linear differential equations are the dynamic manifestation of linear algebra. The fundamental equation $\dot{x} = Ax$ describes a velocity field where the rate of change of the state $x$ is a linear transformation of the state itself. The solution is the **matrix exponential** $e^{At}$, which acts as a "propagator" mapping initial conditions to future states. This chapter shows how diagonalization and the Jordan form allow us to decouple complex interactions into independent modes, and how the eigenvalues of $A$ determine whether the system is stable, oscillating, or exploding.

---

## 26.1 The Matrix Propagator

!!! definition "Definition 26.1 (Homogeneous Linear System)"
    A system of the form $\dot{x}(t) = Ax(t)$, where $x(0) = x_0$. The solution is $x(t) = e^{At}x_0$.

!!! theorem "Theorem 26.1 (Stability of the Origin)"
    The equilibrium point $x=0$ is:
    1. **Asymptotically Stable** if all $\operatorname{Re}(\lambda_i) < 0$.
    2. **Unstable** if any $\operatorname{Re}(\lambda_i) > 0$.
    3. **Marginally Stable** if all $\operatorname{Re}(\lambda_i) \le 0$ and those with $\operatorname{Re}(\lambda_i) = 0$ are non-defective.

---

## Exercises

1. **[Fundamentals] Solve $\dot{x}=3x, \dot{y}=2y$ with $x(0)=1, y(0)=4$.**
   ??? success "Solution"
       The system is decoupled. $x(t) = e^{3t}$ and $y(t) = 4e^{2t}$.

2. **[Diagonalization] Solve $\dot{x} = Ax$ where $A = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$.**
   ??? success "Solution"
       Eigenvalues are 1 and 2. $A = PDP^{-1}$ with $D = \operatorname{diag}(1, 2)$ and $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$. Solution is $x(t) = P e^{Dt} P^{-1} x_0 = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} e^t & 0 \\ 0 & e^{2t} \end{pmatrix} \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix} x_0$.

3. **[Stability] Is the system with $A = \begin{pmatrix} -1 & 100 \\ 0 & -2 \end{pmatrix}$ stable?**
   ??? success "Solution"
       Yes, asymptotically stable. The eigenvalues are $\{-1, -2\}$, both of which have negative real parts. Note that the large off-diagonal element may cause transient growth (Ch43).

4. **[Higher Order] Convert $\ddot{y} + 3\dot{y} + 2y = 0$ into a first-order system.**
   ??? success "Solution"
       Let $x_1 = y, x_2 = \dot{y}$. Then $\dot{x}_1 = x_2$ and $\dot{x}_2 = -2x_1 - 3x_2$. Matrix $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$.

5. **[Non-homogeneous] State the Variation of Parameters formula for $\dot{x} = Ax + f(t)$.**
   ??? success "Solution"
       $x(t) = e^{At}x_0 + \int_0^t e^{A(t-\tau)} f(\tau) d\tau$. This is the convolution of the impulse response with the input.

6. **[Jordan Form] Solve $\dot{x} = Jx$ for $J = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$.**
   ??? success "Solution"
       $e^{Jt} = e^{\lambda t} \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$. Thus $x_1(t) = e^{\lambda t}(x_1(0) + t x_2(0))$ and $x_2(t) = e^{\lambda t} x_2(0)$.

7. **[Phase Portrait] Describe the trajectories of a system with eigenvalues $\pm i\omega$.**
   ??? success "Solution"
       They are concentric ellipses centered at the origin. The system is a harmonic oscillator (center).

8. **[Invariance] If $v$ is an eigenvector of $A$, what happens to the trajectory starting at $x_0 = v$?**
   ??? success "Solution"
       $x(t) = e^{At}v = e^{\lambda t}v$. The state stays on the line spanned by $v$, merely scaling by $e^{\lambda t}$.

9. **[Fundamental Matrix] What is the Fundamental Matrix $\Phi(t)$?**
   ??? success "Solution"
       A matrix whose columns are linearly independent solutions. For LTI systems, $\Phi(t) = e^{At}$. It satisfy $\dot{\Phi} = A\Phi$.

10. **[Symmetry] If $A$ is skew-symmetric, what property does $e^{At}$ have?**
    ??? success "Solution"
        $e^{At}$ is an orthogonal matrix. Trajectories preserve the Euclidean norm $\|x(t)\|_2 = \|x_0\|_2$, corresponding to pure rotation in phase space.

## Chapter Summary

This chapter establishes the analytic solution to linear evolution:

1. **Propagator Calculus**: Defined the matrix exponential as the universal solution to linear time-invariant systems.
2. **Modal Decoupling**: Utilized spectral decomposition to transform coupled equations into independent scalar growth laws.
3. **Equilibrium Topology**: Formulated the stability criteria based on the real parts of eigenvalues, categorizing phase space behaviors.
4. **Convolution Logic**: Extended the theory to non-homogeneous systems, providing the mathematical basis for control and signal response.
