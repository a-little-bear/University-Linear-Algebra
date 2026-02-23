# Chapter 25: Linear Algebra in Optimization

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Matrix Decompositions (Ch10) · Positive Definite Matrices (Ch16) · SVD (Ch11) · Vector Calculus

**Chapter Outline**: The Geometry of Optimization → Linear Programming (LP) and the Simplex Method → Quadratic Programming (QP) & Hessian Matrices → Semidefinite Programming (SDP) → Matrix Completion & Nuclear Norm Minimization → Compressed Sensing & Sparsity → First-order Methods (Gradient Descent, Acceleration) → Second-order Methods (Newton, Quasi-Newton, BFGS) → Applications: Finance, Machine Learning, and Robotics

**Extension**: Linear Algebra is the "engine" of optimization; every step of an optimization algorithm involves solving a linear system or a spectral problem. Modern AI is essentially large-scale matrix optimization.

</div>

Optimization is the science of finding the "best" solution under constraints. Whether it's minimizing cost in a supply chain, maximizing return in a portfolio, or training a billion-parameter neural network, the underlying problem always reduces to linear algebra. This chapter explores the algebraic structures that enable efficient optimization, from the polyhedral geometry of Linear Programming to the cone-based theory of Semidefinite Programming.

---

## 25.1 Linear and Quadratic Programming

<div class="context-flow" markdown>

**Linear Programming (LP)**: $\min c^T x$ subject to $Ax = b, x \ge 0$.
**Quadratic Programming (QP)**: $\min \frac{1}{2} x^T Q x + c^T x$ subject to linear constraints.

</div>

!!! theorem "Theorem 25.1 (Optimality Conditions for QP)"
    For a QP problem with a symmetric positive definite matrix $Q$, the global minimum is achieved when the gradient vanishes:
    $$\nabla f(x) = Qx + c = 0 \implies x^* = -Q^{-1}c$$
    If $Q$ is not positive definite, the problem may be unbounded or have no local minimum.

!!! technique "Interior Point Methods"
    Unlike the Simplex method which walks along the edges of a polytope, Interior Point methods travel through the center of the feasible region. This requires solving a sequence of linear systems involving the matrix $A D A^T$, where $D$ is a diagonal scaling matrix.

---

## 25.2 Semidefinite Programming (SDP)

<div class="context-flow" markdown>

**The Frontier of Optimization**: SDP is the most powerful generalization of LP, where variables are matrices and constraints involve the **Positive Definite Cone** ($\mathcal{S}_n^+$).

</div>

!!! definition "Definition 25.1 (Semidefinite Program)"
    An SDP problem takes the form:
    $$\min_{X \in \mathcal{S}_n} \langle C, X \rangle \quad \text{subject to } \langle A_i, X \rangle = b_i, \quad X \succeq 0$$
    where $\langle A, B \rangle = \operatorname{tr}(A^T B)$ is the Frobenius inner product.

!!! theorem "Theorem 25.2 (SDP Duality)"
    Every SDP has a dual problem:
    $$\max_{y \in \mathbb{R}^m} b^T y \quad \text{subject to } \sum y_i A_i \preceq C$$
    Under strong duality (Slater's condition), the optimal values of the primal and dual are equal.

---

## 25.3 Matrix Completion and Nuclear Norm

!!! technique "The Nuclear Norm Trick"
    In problems like the "Netflix Challenge" (recommender systems), we want to fill in a large matrix $M$ using only a few entries. This is a rank-minimization problem, which is NP-hard. We relax it to a convex problem using the **Nuclear Norm** (sum of singular values):
    $$\min \|X\|_* \quad \text{subject to } X_{ij} = M_{ij} \text{ for observed } (i,j)$$
    This is equivalent to an SDP and can be solved efficiently.

---

## 25.4 Compressed Sensing and Sparsity

!!! theorem "Theorem 25.3 (L1-Relaxation)"
    Finding the sparsest solution to $Ax=b$ (minimizing $\|x\|_0$) is NP-hard. However, if $A$ satisfies the **Restricted Isometry Property** (RIP), the sparsest solution is exactly the same as the one that minimizes the $L_1$ norm ($\sum |x_i|$). $L_1$ minimization is a Linear Program.

---

## 25.5 Second-order Methods: Newton and BFGS

<div class="context-flow" markdown>

**Core Idea**: Second-order methods use the curvature of the function (the Hessian $H$) to take faster, smarter steps.

</div>

!!! algorithm "Algorithm 25.1 (Newton's Method)"
    1.  Approximate $f(x)$ by a quadratic Taylor expansion at $x_k$.
    2.  The update is $x_{k+1} = x_k - H(x_k)^{-1} \nabla f(x_k)$.
    3.  **Linear Algebra cost**: Solving the linear system $H \Delta x = -\nabla f$ in each step.

!!! technique "Quasi-Newton (BFGS)"
    In many large-scale problems, computing the full Hessian $H$ is too expensive ($O(n^2)$ storage, $O(n^3)$ inversion). BFGS approximates $H^{-1}$ using only first-order gradients through a sequence of rank-2 updates, maintaining positive definiteness and requiring only $O(n^2)$ work.

---

## Exercises

1.  **[Quadratic] Find the minimum of $f(x, y) = x^2 + 2y^2 - 4x + 4y + 7$.**
    ??? success "Solution"
        The Hessian is $H = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}$, which is PD.
        The gradient is $\nabla f = (2x-4, 4y+4)^T$.
        Setting $\nabla f = 0$ gives $x=2, y=-1$.
        The minimum value is $f(2, -1) = 4 + 2 - 8 - 4 + 7 = 1$.

2.  **[KKT] State the Karush-Kuhn-Tucker (KKT) conditions for equality-constrained optimization.**
    ??? success "Solution"
        1. Primal feasibility: $h(x) = 0$.
        2. Stationarity: $\nabla f(x) + \sum \lambda_i \nabla h_i(x) = 0$.

3.  **[SDP] Express the constraint $\lambda_{\max}(X) \le 1$ as an LMI (Linear Matrix Inequality).**
    ??? success "Solution"
        $X \preceq I$, or $I - X \succeq 0$.

4.  **[SVD in Opt] Why is SVD used in low-rank matrix recovery?**
    ??? success "Solution"
        The proximal operator of the nuclear norm involves a "singular value thresholding" step: computing the SVD and shrinking all singular values towards zero.

5.  **[Nuclear Norm] Prove $\|X\|_* = \operatorname{tr}(\sqrt{X^T X})$.**
    ??? success "Solution"
        $\sqrt{X^T X}$ has eigenvalues equal to the singular values $\sigma_i$ of $X$. Since $\|X\|_* = \sum \sigma_i$ and the trace is the sum of eigenvalues, the equality holds.

6.  **[Convergence] Why does Newton's method converge faster than Gradient Descent?**
    ??? success "Solution"
        Gradient Descent only uses the first derivative (linear approximation), while Newton's uses the second derivative (quadratic approximation). Newton's method has quadratic convergence near the optimum, whereas Gradient Descent is linear.

7.  **[Conditioning] How does the condition number $\kappa(H)$ affect Gradient Descent?**
    ??? success "Solution"
        High condition numbers result in "elongated" contours (ellipses). Gradient Descent "zig-zags" in these narrow valleys, leading to very slow convergence.

8.  **[Dual] Write the dual of the LP $\min c^T x$ s.t. $Ax \ge b, x \ge 0$.**
    ??? success "Solution"
        $\max b^T y$ s.t. $A^T y \le c, y \ge 0$.

9.  **[SVM] Why is the Support Vector Machine (SVM) a Quadratic Program?**
    ??? success "Solution"
        Training an SVM involves maximizing the margin $2/\|w\|$, which is equivalent to minimizing $\frac{1}{2} \|w\|^2$ (a quadratic objective) subject to linear constraints $y_i(w^T x_i + b) \ge 1$.

10. **[StQP] What is a Standard Quadratic Program?**
    ??? success "Solution"
        Minimizing a quadratic form $x^T Q x$ over the standard simplex $\{x \ge 0, \sum x_i = 1\}$. This is NP-hard if $Q$ is indefinite but can be solved via copositive programming.

## Chapter Summary

Optimization is the ultimate application of Matrix Analysis:

1.  **Algebraic Search**: Established that solving optimization problems is synonymous with solving linear systems or finding spectral bounds.
2.  **Cone Geometry**: Generalized linear constraints to the positive definite cone, creating the powerful framework of SDP.
3.  **Curvature Utilization**: Developed second-order methods that exploit the Hessian matrix to achieve superior convergence rates.
4.  **Information Exploitation**: Demonstrated how linear algebra tricks like the nuclear norm and L1-minimization can extract low-dimensional signals from high-dimensional noise.
