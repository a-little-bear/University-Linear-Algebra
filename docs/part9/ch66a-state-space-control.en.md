# Chapter 66A: State-Space and System Realization

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch01) · Eigenvalues (Ch06) · Matrix Functions (Ch13) · Matrix Pencils (Ch41A)

**Chapter Outline**: From External Descriptions to Internal Representations → Definition of State-Space Equations → System, Input, Output, and Feedthrough Matrices → Structure of Solutions: Natural and Forced Responses → Core Criteria: Controllability and Observability → Kalman Decomposition Theorem → Transfer Function Matrices and Algebraic Equivalence → System Realization: Mapping from Transfer Functions to State-Space → Minimal Realization → Applications: Multivariable Feedback Control, Aerospace Guidance, and State Modeling of Circuit Systems

**Extension**: State-space methods are the soul of modern control theory; they encapsulate the physical state of a system as a vector and map causality to matrix operators. They prove that the evolution of complex systems depends not just on inputs but on the geometric evolution of internal energy states—the main artery connecting linear algebra to engineering automation.

</div>

In classical control, we use transfer functions to describe the relationship between inputs and outputs. However, this "black box" view fails to reveal the internal dynamics of the system. **State-Space Representation** introduces a "state vector" $\mathbf{x}(t)$ to transform complex dynamics into a set of linear matrix equations. It allows us to precisely analyze which internal states can be manipulated by external forces (**Controllability**) and which can be detected via sensors (**Observability**). This chapter introduces the algebraic framework at the heart of automatic control and systems science.

---

## 66A.1 State-Space Equations

!!! definition "Definition 66A.1 (LTI State-Space)"
    $$\begin{cases} \dot{\mathbf{x}}(t) = A\mathbf{x}(t) + B\mathbf{u}(t) \\ \mathbf{y}(t) = C\mathbf{x}(t) + D\mathbf{u}(t) \end{cases}$$
    - $A$: **System Matrix** (intrinsic dynamics).
    - $B$: **Input Matrix** (control action).
    - $C$: **Output Matrix** (observation method).
    - $D$: **Feedthrough Matrix**.

---

## 66A.2 Controllability and Observability

!!! theorem "Theorem 66A.1 (Kalman Criteria)"
    1.  **Controllability**: A system is completely controllable iff the **Controllability Matrix** $\mathcal{C} = [B, AB, A^2B, \ldots, A^{n-1}B]$ has full row rank.
    2.  **Observability**: A system is completely observable iff the **Observability Matrix** $\mathcal{O} = [C^T, (CA)^T, \ldots, (CA^{n-1})^T]^T$ has full column rank.

---

## 66A.3 System Realization and Transfer Functions

!!! technique "Technique: Transfer Function Matrix"
    The input-output behavior under Laplace transformation is given by:
    $$G(s) = C(sI - A)^{-1}B + D$$
    **Minimal Realization**: A realization $(A, B, C, D)$ is minimal (lowest possible dimension) iff it is both completely controllable and completely observable.

---

## Exercises

**1. [Basics] Given $A = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix}$ and $B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$, calculate the controllability matrix $\mathcal{C}$.**

??? success "Solution"
    **Steps:**
    1. Compute $AB = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix} \begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \\ -3 \end{pmatrix}$.
    2. Construct $\mathcal{C} = [B, AB] = \begin{pmatrix} 0 & 1 \\ 1 & -3 \end{pmatrix}$.
    **Conclusion**: $\det(\mathcal{C}) = -1 \neq 0$, so the matrix has full rank. The system is **completely controllable**.

**2. [Observability] If $C = \begin{pmatrix} 1 & 0 \end{pmatrix}$, determine if the system above is observable.**

??? success "Solution"
    **Steps:**
    1. Compute $CA = \begin{pmatrix} 1 & 0 \end{pmatrix} \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix} = \begin{pmatrix} 0 & 1 \end{pmatrix}$.
    2. Construct $\mathcal{O} = \begin{pmatrix} C \\ CA \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.
    **Conclusion**: $\mathcal{O} = I$, which has full rank. The system is **completely observable**.

**3. [Response] What is the solution to the system under zero input $\mathbf{u}(t)=0$?**

??? success "Solution"
    **Conclusion: $\mathbf{x}(t) = e^{At}\mathbf{x}_0$.**
    This is known as the **zero-input response** or free response. It is entirely determined by the eigenvalues of the system matrix $A$.

**4. [Transfer Function] Prove: If $D=0, B=b, C=c^T$ (SISO), the transfer function is a scalar fraction.**

??? success "Solution"
    **Proof:**
    $G(s) = c^T(sI - A)^{-1}b$. Since $(sI-A)^{-1}$ is an $n \times n$ matrix, left-multiplication by row vector $c^T$ and right-multiplication by column vector $b$ yields a $1 \times 1$ scalar. Using the adjugate formula, $G(s) = \frac{c^T \operatorname{adj}(sI-A) b}{\det(sI-A)}$, where the denominator is the characteristic polynomial.

**5. [PBH Test] What is the PBH test for controllability?**

??? success "Solution"
    **Description:**
    The pair $(A, B)$ is controllable iff the matrix pencil $[ \lambda I - A \ | \ B ]$ has full row rank for all eigenvalues $\lambda$ of $A$.
    **Significance**: This implies that the control action $B$ must be able to affect every eigenmode of the system.

**6. [Calculation] For $A = \operatorname{diag}(1, -1)$ and $B = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, is the system controllable?**

??? success "Solution"
    **Determination:**
    1. $\mathcal{C} = [B, AB] = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$.
    2. Rank is 1, which is less than the dimension 2.
    **Conclusion**: Not controllable.
    **Intuition**: Because the weight of the input on the second eigenvalue (-1) is 0 and $A$ is diagonal (no coupling), the second state remains untouched.

**7. [Realization] Does a single transfer function $G(s)$ correspond to a unique state-space realization?**

??? success "Solution"
    **Conclusion: No.**
    **Reasoning**: If we perform a coordinate transformation $\mathbf{z} = P\mathbf{x}$, the new system $(P A P^{-1}, PB, CP^{-1}, D)$ yields the exact same transfer function. Physical implementations have the freedom of basis change.

**8. [Minimal] If a realization is minimal, what can be said about its dimension $n$?**

??? success "Solution"
    **Conclusion: $n$ is the degree of the denominator of the transfer function.**
    Minimal realization implies there are no redundant, uncontrollable, or unobservable states.

**9. [Duality] Briefly state the duality between controllability and observability.**

??? success "Solution"
    **Theorem:**
    The pair $(A, B)$ is controllable iff the dual pair $(A^T, C^T)$ is observable. This symmetry allows us to transform controller design problems into observer (state estimation) design problems.

**10. [Application] Why is state-space preferred over transfer functions in aerospace guidance?**

??? success "Solution"
    **Reasoning:**
    1. **Multivariable Coupling**: Aircraft have 6 degrees of freedom; inputs and outputs are highly coupled, making transfer function matrices extremely messy.
    2. **Time-variance**: Mass (fuel consumption) and aerodynamic parameters change during flight; state-space handles $A(t)$ naturally.
    3. **Modern Optimization**: State-space allows the direct application of optimal control algorithms like LQR (based on eigenvalues) for precision guidance.

## Chapter Summary

State-space and system realization are the "brute force output" of linear algebra into modern engineering:

1.  **Transparency of the Internal**: It deconstructs complex dynamics into trajectories of state vectors in multi-dimensional space, proving that system solvability is essentially the integration of linear operators over time.
2.  **Boundaries of Ability**: Controllability and observability criteria establish the algebraic limits of human intervention in natural systems, defining the mathematical borders of the "knowable" and "controllable."
3.  **Flexibility of Realization**: By revealing the many-to-one mapping from state-space to transfer functions, the theory provides engineering designs with immense coordinate freedom, supporting the transition from classical feedback to modern intelligent control.
