# Chapter 66A: State-Space Control Theory

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Matrix Equations (Ch20)

**Chapter Outline**: Dynamical Background of Control Theory → State-Space Representation $(A, B, C, D)$ → Controllability Definition and the Controllability Matrix → Observability Definition and the Observability Matrix → Pole Placement & State Feedback → Observer Design (Luenberger Observers) → Transfer Functions from State-Space → Realization Theory → Introduction to Sampled-Data Systems

**Extension**: The state-space method is the cornerstone of modern control theory; it transforms external input-output descriptions into internal state evolutions, enabling the use of matrix theory (e.g., rank, subspace decomposition) to solve complex problems in industrial automation and aerospace guidance.

</div>

In classical control theory, the Laplace transform is used to handle Single-Input Single-Output (SISO) systems. **State-Space Methods**, however, utilize linear algebra to handle Multi-Input Multi-Output (MIMO) systems. It compresses all historical information of a system into a "state vector" and describes the system's evolution as a first-order matrix differential equation in state space. This chapter demonstrates how the profound algebraic properties of controllability and observability dictate the control limits of a system.

---

## 66A.1 State-Space Representation

!!! definition "Definition 66A.1 (Linear Time-Invariant - LTI System)"
    The state-space model of an LTI system is defined as:
    $$\begin{cases} \dot{\mathbf{x}}(t) = A\mathbf{x}(t) + B\mathbf{u}(t) \\ \mathbf{y}(t) = C\mathbf{x}(t) + D\mathbf{u}(t) \end{cases}$$
    - $\mathbf{x} \in \mathbb{R}^n$: State vector.
    - $\mathbf{u} \in \mathbb{R}^m$: Input vector.
    - $\mathbf{y} \in \mathbb{R}^p$: Output vector.

---

## 66A.2 Controllability and Observability

!!! definition "Definition 66A.2 (Controllability)"
    A system $(A, B)$ is **controllable** if it is possible to steer the system from any initial state to any target state in finite time using an appropriate input.
    **Criterion**: The **Controllability Matrix** $\mathcal{C} = [B \ AB \ A^2B \ \cdots \ A^{n-1}B]$ must have **full rank** ($n$).

!!! definition "Definition 66A.3 (Observability)"
    A system $(A, C)$ is **observable** if the initial state can be uniquely determined from the observation of input and output over a finite time interval.
    **Criterion**: The **Observability Matrix** $\mathcal{O} = [C^T \ (CA)^T \ \cdots \ (CA^{n-1})^T]^T$ must have **full rank** ($n$).

---

## 66A.3 Pole Placement and Feedback

!!! technique "State Feedback"
    By choosing a control law $u = -Kx$, the closed-loop system dynamics become $\dot{x} = (A-BK)x$.
    If the system is controllable, the closed-loop poles (eigenvalues of $A-BK$) can be placed at any desired locations in the complex plane by selecting an appropriate gain matrix $K$. This provides an algebraic guarantee for stabilizing systems.

---

## Exercises

1.  **[Basics] Write the state-space representation for $y'' + 3y' + 2y = u$.**
    ??? success "Solution"
        Let $x_1 = y, x_2 = y'$. Then $\dot{x}_1 = x_2$ and $\dot{x}_2 = -2x_1 - 3x_2 + u$.
        $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}, B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}, C = \begin{pmatrix} 1 & 0 \end{pmatrix}$.

2.  **[Controllability] Determine if $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}, B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ is controllable.**
    ??? success "Solution"
        $\mathcal{C} = [B \ AB] = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$. The rank is 1 < 2, so the system is uncontrollable.

3.  **[Observability] If $C = \begin{pmatrix} 1 & 0 \end{pmatrix}$ and $A$ is upper triangular, is the system always observable?**
    ??? success "Solution"
        Not necessarily. One must check the rank of $\mathcal{O}$. If $A = \operatorname{diag}(1, 1)$, the second state never appears in the output, making it unobservable.

4.  **[Transfer Function] How is the transfer function $G(s)$ of an LTI system represented in matrix form?**
    ??? success "Solution"
        $G(s) = C(sI - A)^{-1}B + D$.

5.  **[Stability] What is the condition for the stability of the closed-loop system $(A-BK)$?**
    ??? success "Solution"
        $A-BK$ must be Hurwitz stable (all eigenvalues have negative real parts).

6.  **[Observer] What is the error equation for a Luenberger observer?**
    ??? success "Solution"
        $\dot{e} = (A-LC)e$. We design $L$ such that $A-LC$ is stable, ensuring the estimate converges to the true state.

7.  **[Cayley-Hamilton] Why does the controllability matrix only go up to $A^{n-1}B$?**
    ??? success "Solution"
        By the Cayley-Hamilton theorem, $A^n$ and higher powers are linear combinations of lower powers and will not increase the rank of the span.

8.  **[Calculation] Find the poles of $\dot{x} = x+u$ under the feedback $u=-2x$.**
    ??? success "Solution"
        $\dot{x} = (1-2)x = -x$. The pole is at -1 (the system is stabilized).

9.  **[Realization] What is a minimal realization?**
    ??? success "Solution"
        A state-space realization of minimum dimension. A realization is minimal iff it is both controllable and observable.

10. **[Application] Why are state-space methods mandatory in aerospace?**

   ??? success "Solution"
        Because aircraft have multiple highly coupled control surfaces (MIMO) and involve numerous internal sensors (states). State-space methods provide a unified framework for optimal control and filtering (e.g., Kalman filtering).

## Chapter Summary

State-space control theory is the dynamical pinnacle of linear algebra:

1.  **Discovery of Internal States**: It makes the hidden logic of "black box" systems explicit, proving that the entire future evolution of a system is stored in the interaction between the current vector and the operator.
2.  **Physical Meaning of Rank**: Controllability and observability transform dry matrix rank definitions into physical measures of "intervention capability" and "depth of insight," defining the algebraic boundaries of human control over the physical world.
3.  **Pole Design**: Through eigenvalue placement, control theory demonstrates how to use algebraic feedback to reshape an operator's spectrum, forcing chaotic or unstable systems to exhibit desired smoothness and response speeds.
