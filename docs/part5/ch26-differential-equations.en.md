# Chapter 26  Linear Algebra in Differential Equations

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues and diagonalization (Ch6) · Jordan canonical form (Ch12) · Matrix functions (Ch13) · **Chapter arc**: $\mathbf{x}' = A\mathbf{x}$ (homogeneous systems) → $e^{At}$ (matrix exponential) → real parts of eigenvalues (stability) → phase plane (geometric classification) → companion matrix (higher-order to systems) → variation of parameters (nonhomogeneous) → Sturm-Liouville (PDEs) → Floquet (periodic systems)
**Essence**: Both qualitative and quantitative theory of differential equations rely on linear algebra at their core — eigenvalues determine stability, matrix exponentials give solutions, Jordan form handles repeated roots

</div>

Linear differential equations represent one of the most classical and important application domains of linear algebra. The solutions of linear constant-coefficient ODE systems can be completely described by the matrix exponential $e^{At}$, while the long-term behavior of the system (stability, oscillation, decay) is entirely determined by the eigenvalues of the coefficient matrix $A$. This chapter begins with first-order constant-coefficient linear systems and progresses through the matrix exponential, stability analysis, phase plane classification, higher-order equations, nonhomogeneous equations, eigenvalue problems in PDEs, and Floquet theory for periodic-coefficient systems.

---

## 26.1 Linear ODE Systems

<div class="context-flow" markdown>

**Core structure**: The solution space of $\mathbf{x}' = A\mathbf{x}$ is an $n$-dimensional linear space → each column of the fundamental matrix $\Phi(t)$ is a solution → $\Phi(t) = e^{At}$ when $\Phi(0) = I$
**Link**: The dimension theorem from Ch4 (linear spaces) is directly manifested here

</div>

Linear constant-coefficient ODE systems are fundamental models in dynamical systems, control theory, and physics.

!!! definition "Definition 26.1 (Linear constant-coefficient ODE system)"
    Let $A \in \mathbb{R}^{n \times n}$ be a constant matrix. A **linear constant-coefficient ordinary differential equation system** is

    $$
    \mathbf{x}'(t) = A\mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$

    where $\mathbf{x}(t) \in \mathbb{R}^n$ is the state vector, $t \ge 0$.

!!! definition "Definition 26.2 (Fundamental matrix)"
    An $n \times n$ matrix-valued function $\Phi(t)$ is called a **fundamental matrix** of the equation $\mathbf{x}' = A\mathbf{x}$ if

    $$
    \Phi'(t) = A\Phi(t), \quad \det \Phi(t) \neq 0, \quad \forall\, t.
    $$

    If furthermore $\Phi(0) = I$, then $\Phi(t)$ is called the **principal fundamental matrix** or **state transition matrix**.

!!! theorem "Theorem 26.1 (Solution space structure)"
    The solution set of the homogeneous equation $\mathbf{x}' = A\mathbf{x}$ forms an $n$-dimensional linear space over $\mathbb{R}^n$. If $\Phi(t)$ is any fundamental matrix, the general solution is

    $$
    \mathbf{x}(t) = \Phi(t)\mathbf{c}, \quad \mathbf{c} \in \mathbb{R}^n.
    $$

    The unique solution satisfying the initial condition $\mathbf{x}(0) = \mathbf{x}_0$ is $\mathbf{x}(t) = \Phi(t)\Phi(0)^{-1}\mathbf{x}_0$.

??? proof "Proof"
    **Linearity**: If $\mathbf{x}_1(t)$ and $\mathbf{x}_2(t)$ are solutions, then $\alpha \mathbf{x}_1(t) + \beta \mathbf{x}_2(t)$ is also a solution, since

    $$
    (\alpha \mathbf{x}_1 + \beta \mathbf{x}_2)' = \alpha A\mathbf{x}_1 + \beta A\mathbf{x}_2 = A(\alpha \mathbf{x}_1 + \beta \mathbf{x}_2).
    $$

    **Dimension is $n$**: The ODE existence and uniqueness theorem guarantees a unique solution for each initial value $\mathbf{x}_0 \in \mathbb{R}^n$. The map from initial values to solutions $\mathbf{x}_0 \mapsto \mathbf{x}(t)$ is a linear bijection, so the solution space is isomorphic to $\mathbb{R}^n$ and has dimension $n$.

    The $n$ columns of the fundamental matrix $\Phi(t)$ are linearly independent (since $\det \Phi(t) \neq 0$), forming a basis for the solution space. $\blacksquare$

!!! example "Example 26.1"
    **Solution of a two-dimensional linear system.** Consider

    $$
    \mathbf{x}' = \begin{pmatrix} 0 & 1 \\ -2 & -3 \end{pmatrix} \mathbf{x}.
    $$

    The characteristic equation $\lambda^2 + 3\lambda + 2 = 0$ gives $\lambda_1 = -1$, $\lambda_2 = -2$. The corresponding eigenvectors are $\mathbf{v}_1 = (1, -1)^T$, $\mathbf{v}_2 = (1, -2)^T$. The general solution is

    $$
    \mathbf{x}(t) = c_1 e^{-t} \begin{pmatrix} 1 \\ -1 \end{pmatrix} + c_2 e^{-2t} \begin{pmatrix} 1 \\ -2 \end{pmatrix}.
    $$

    The fundamental matrix is $\Phi(t) = \begin{pmatrix} e^{-t} & e^{-2t} \\ -e^{-t} & -2e^{-2t} \end{pmatrix}$, satisfying $\Phi'(t) = A\Phi(t)$.

!!! theorem "Theorem 26.2 (Liouville's formula)"
    Let $\Phi(t)$ be a fundamental matrix of $\mathbf{x}' = A(t)\mathbf{x}$ (where $A(t)$ may be time-varying). Then

    $$
    \det \Phi(t) = \det \Phi(t_0) \cdot \exp\!\left(\int_{t_0}^{t} \operatorname{tr} A(s)\, ds\right).
    $$

    In particular, for constant-coefficient systems $A(t) \equiv A$, we have $\det \Phi(t) = \det \Phi(0) \cdot e^{t\, \operatorname{tr}(A)}$.

??? proof "Proof"
    Let $W(t) = \det \Phi(t)$ (the Wronskian). Using the formula for differentiating a determinant with respect to a matrix,

    $$
    W'(t) = \sum_{i=1}^{n} \det \Phi_i(t),
    $$

    where $\Phi_i(t)$ is obtained by replacing the $i$-th row of $\Phi(t)$ with its derivative. Since $\Phi' = A\Phi$, the derivative of the $i$-th row equals the product of the $i$-th row of $A$ with $\Phi$. After expansion, only the diagonal entries $a_{ii}$ contribute to the determinant, yielding $W'(t) = \operatorname{tr}(A(t)) W(t)$. Solving this scalar ODE gives the result. $\blacksquare$

---

## 26.2 Matrix Exponential and Solution Representation

<div class="context-flow" markdown>

**Core formula**: $e^{At} = \sum_{k=0}^{\infty} \frac{(At)^k}{k!}$ → when diagonalizable, $e^{At} = P\, \text{diag}(e^{\lambda_i t})\, P^{-1}$ → for Jordan blocks, terms involve $t^k e^{\lambda t}$
**Link**: The most important special case of Ch13 matrix function theory

</div>

The matrix exponential is the central tool in linear ODE theory, perfectly unifying matrix algebra and differential equations.

!!! definition "Definition 26.3 (Matrix exponential)"
    For $A \in \mathbb{R}^{n \times n}$, the **matrix exponential** is defined as

    $$
    e^{A} = \sum_{k=0}^{\infty} \frac{A^k}{k!} = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots
    $$

    This series converges absolutely for all $A$.

!!! theorem "Theorem 26.3 (Matrix exponential and ODEs)"
    The unique solution of $\mathbf{x}' = A\mathbf{x}$, $\mathbf{x}(0) = \mathbf{x}_0$ is

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0.
    $$

    $e^{At}$ is the principal fundamental matrix: $\frac{d}{dt}e^{At} = Ae^{At}$, $e^{A \cdot 0} = I$.

??? proof "Proof"
    Differentiating term by term:

    $$
    \frac{d}{dt}e^{At} = \frac{d}{dt}\sum_{k=0}^{\infty} \frac{(At)^k}{k!} = \sum_{k=1}^{\infty} \frac{A^k t^{k-1}}{(k-1)!} = A\sum_{k=1}^{\infty} \frac{(At)^{k-1}}{(k-1)!} = Ae^{At}.
    $$

    Therefore $\mathbf{x}(t) = e^{At}\mathbf{x}_0$ satisfies $\mathbf{x}' = Ae^{At}\mathbf{x}_0 = A\mathbf{x}(t)$ and $\mathbf{x}(0) = I\mathbf{x}_0 = \mathbf{x}_0$. By the uniqueness theorem, this is the unique solution. $\blacksquare$

!!! theorem "Theorem 26.4 (Computing the matrix exponential)"
    Let $A = PJP^{-1}$ be the Jordan decomposition, $J = \operatorname{diag}(J_1, \ldots, J_s)$, where each $J_i$ is a Jordan block. Then

    $$
    e^{At} = P\, \operatorname{diag}(e^{J_1 t}, \ldots, e^{J_s t})\, P^{-1}.
    $$

    For a $k \times k$ Jordan block $J_i = \lambda_i I + N$ ($N$ nilpotent),

    $$
    e^{J_i t} = e^{\lambda_i t} \begin{pmatrix} 1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\ 0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\ \vdots & & \ddots & & \vdots \\ 0 & 0 & \cdots & 1 & t \\ 0 & 0 & \cdots & 0 & 1 \end{pmatrix}.
    $$

??? proof "Proof"
    From $A = PJP^{-1}$ we get $A^k = PJ^k P^{-1}$, hence

    $$
    e^{At} = \sum_{k=0}^{\infty} \frac{(PJP^{-1})^k t^k}{k!} = P\left(\sum_{k=0}^{\infty} \frac{J^k t^k}{k!}\right)P^{-1} = Pe^{Jt}P^{-1}.
    $$

    Since $J$ is block diagonal, $e^{Jt} = \operatorname{diag}(e^{J_1 t}, \ldots, e^{J_s t})$. For a Jordan block $J_i = \lambda_i I + N$, since $\lambda_i I$ and $N$ commute, $e^{J_i t} = e^{\lambda_i t} e^{Nt}$, where $N^k = 0$ (for $k \ge$ block size), so $e^{Nt}$ is a finite sum, yielding the upper triangular matrix form. $\blacksquare$

!!! example "Example 26.2"
    **Matrix exponential of a system with repeated eigenvalues.** Let

    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}.
    $$

    $A$ has a single eigenvalue $\lambda = 2$ with Jordan block $J = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ ($P = I$).

    $$
    e^{At} = e^{2t}\begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}.
    $$

    For initial value $\mathbf{x}_0 = (1, 3)^T$, the solution is

    $$
    \mathbf{x}(t) = e^{2t}\begin{pmatrix} 1 + 3t \\ 3 \end{pmatrix}.
    $$

    The presence of the $te^{2t}$ term is the hallmark behavior of Jordan blocks (non-diagonalizable matrices).

!!! definition "Definition 26.4 (Fundamental properties of the matrix exponential)"
    The matrix exponential satisfies the following properties:

    1. $e^{O} = I$ ($O$ is the zero matrix).
    2. $(e^{A})^{-1} = e^{-A}$.
    3. $e^{(s+t)A} = e^{sA}e^{tA}$ (semigroup property).
    4. If $AB = BA$, then $e^{A+B} = e^{A}e^{B}$.
    5. $\det(e^{A}) = e^{\operatorname{tr}(A)}$.
    6. In general, $e^{A+B} \neq e^{A}e^{B}$.

!!! example "Example 26.3"
    **Non-commutativity of the matrix exponential.** Let $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$, $B = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$.

    $e^A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$, $e^B = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$, $e^A e^B = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$.

    But $A + B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ has eigenvalues $\pm 1$,

    $$
    e^{A+B} = \begin{pmatrix} \cosh 1 & \sinh 1 \\ \sinh 1 & \cosh 1 \end{pmatrix} \neq e^A e^B.
    $$

    Since $AB \neq BA$, the Baker-Campbell-Hausdorff formula introduces nontrivial correction terms.

---

## 26.3 Stability Analysis

<div class="context-flow" markdown>

**Criterion chain**: All $\operatorname{Re}(\lambda_i) < 0$ $\iff$ asymptotically stable → Lyapunov equation $A^T P + PA = -Q$ has positive definite solution $\iff$ Hurwitz matrix
**Applications**: Convergence of control systems, stability of filters in signal processing

</div>

The stability of linear systems is entirely determined by the eigenvalues of the coefficient matrix — one of the most profound applications of linear algebra in dynamical systems.

!!! definition "Definition 26.5 (Stability classification)"
    For the system $\mathbf{x}' = A\mathbf{x}$, the equilibrium $\mathbf{x} = \mathbf{0}$ is called:

    - **Asymptotically stable**: if $\lim_{t \to \infty} \mathbf{x}(t) = \mathbf{0}$ for all initial values $\mathbf{x}_0$;
    - **Stable** (Lyapunov stable): if for every $\varepsilon > 0$, there exists $\delta > 0$ such that $\|\mathbf{x}_0\| < \delta$ implies $\|\mathbf{x}(t)\| < \varepsilon$ for all $t \ge 0$;
    - **Unstable**: if it is not stable.

!!! theorem "Theorem 26.5 (Eigenvalue stability criterion)"
    For the constant-coefficient system $\mathbf{x}' = A\mathbf{x}$:

    1. **Asymptotically stable** $\iff$ all eigenvalues of $A$ satisfy $\operatorname{Re}(\lambda_i) < 0$ ($A$ is a **Hurwitz matrix**).
    2. **Stable** $\iff$ all eigenvalues satisfy $\operatorname{Re}(\lambda_i) \le 0$, and eigenvalues with $\operatorname{Re}(\lambda_i) = 0$ have Jordan blocks of size $1$.
    3. **Unstable** $\iff$ there exists $\operatorname{Re}(\lambda_i) > 0$, or there exists $\operatorname{Re}(\lambda_i) = 0$ with Jordan block size $> 1$.

??? proof "Proof"
    From $\mathbf{x}(t) = e^{At}\mathbf{x}_0$ and the Jordan decomposition $e^{At} = Pe^{Jt}P^{-1}$, each Jordan block contributes terms of the form $t^j e^{\lambda_i t}$ ($0 \le j <$ block size).

    - If $\operatorname{Re}(\lambda_i) < 0$, then $|t^j e^{\lambda_i t}| = t^j e^{\operatorname{Re}(\lambda_i) t} \to 0$ ($t \to \infty$), regardless of $j$.
    - If $\operatorname{Re}(\lambda_i) = 0$ and $j = 0$, then $|e^{\lambda_i t}| = 1$, bounded but not tending to zero.
    - If $\operatorname{Re}(\lambda_i) = 0$ and $j \ge 1$, then $|t^j e^{\lambda_i t}| = t^j \to \infty$.
    - If $\operatorname{Re}(\lambda_i) > 0$, then $|e^{\lambda_i t}| \to \infty$. $\blacksquare$

!!! theorem "Theorem 26.6 (Lyapunov stability theorem)"
    $A$ is a Hurwitz matrix if and only if for any positive definite matrix $Q \succ 0$, the **Lyapunov equation**

    $$
    A^T P + PA = -Q
    $$

    has a unique positive definite solution $P \succ 0$. In this case, $V(\mathbf{x}) = \mathbf{x}^T P \mathbf{x}$ is a Lyapunov function for the system.

??? proof "Proof"
    **Sufficiency**: If $P \succ 0$ satisfies $A^T P + PA = -Q \prec 0$, define $V(\mathbf{x}) = \mathbf{x}^T P \mathbf{x} > 0$. Differentiating along trajectories:

    $$
    \dot{V} = \mathbf{x}'^T P \mathbf{x} + \mathbf{x}^T P \mathbf{x}' = \mathbf{x}^T(A^T P + PA)\mathbf{x} = -\mathbf{x}^T Q \mathbf{x} < 0.
    $$

    **Necessity**: If $A$ is Hurwitz, define $P = \int_0^{\infty} e^{A^T t} Q e^{At}\, dt$ (convergent since $\operatorname{Re}(\lambda_i) < 0$). Direct verification shows $A^T P + PA = -Q$ and $P \succ 0$. Uniqueness follows from the linearity of the Lyapunov equation and the Hurwitz condition. $\blacksquare$

!!! example "Example 26.4"
    **Stability analysis of a circuit system.** An RLC circuit has the state equation

    $$
    A = \begin{pmatrix} 0 & 1 \\ -\frac{1}{LC} & -\frac{R}{L} \end{pmatrix}.
    $$

    The characteristic equation is $\lambda^2 + \frac{R}{L}\lambda + \frac{1}{LC} = 0$, with eigenvalues

    $$
    \lambda_{1,2} = -\frac{R}{2L} \pm \sqrt{\frac{R^2}{4L^2} - \frac{1}{LC}}.
    $$

    Since $R, L, C > 0$, we have $\operatorname{Re}(\lambda_{1,2}) < 0$, so the system is asymptotically stable. Physical meaning: the resistor dissipates energy, and current and voltage eventually decay to zero. When $R^2 > 4L/C$, the system is overdamped (two real negative eigenvalues); when $R^2 < 4L/C$, it is underdamped (complex conjugate eigenvalues, decaying oscillation); when $R^2 = 4L/C$, it is critically damped.

---

## 26.4 Phase Plane Analysis

<div class="context-flow" markdown>

**Classification**: Complete classification of $2 \times 2$ systems depends on $\operatorname{tr}(A)$ and $\det(A)$ → the sign and real/imaginary nature of eigenvalues determine nodes/saddles/spirals/centers
**Geometry**: The qualitative shape of trajectories is completely determined by eigenvalues and eigenvectors

</div>

Phase plane analysis of two-dimensional linear systems $\mathbf{x}' = A\mathbf{x}$ serves as the starting point for understanding higher-dimensional systems.

!!! definition "Definition 26.6 (Classification of critical points)"
    For the two-dimensional system $\mathbf{x}' = A\mathbf{x}$ ($A \in \mathbb{R}^{2 \times 2}$, $\det A \neq 0$), let the eigenvalues be $\lambda_1, \lambda_2$. The type of the critical point at the origin is:

    | Eigenvalue condition | Type |
    |:---:|:---:|
    | $\lambda_1, \lambda_2 \in \mathbb{R}$, same sign | **Node** |
    | $\lambda_1, \lambda_2 \in \mathbb{R}$, opposite signs | **Saddle point** |
    | $\lambda_{1,2} = \alpha \pm \beta i$, $\alpha \neq 0$ | **Spiral/focus** |
    | $\lambda_{1,2} = \pm \beta i$, $\alpha = 0$ | **Center** |
    | $\lambda_1 = \lambda_2 \in \mathbb{R}$, non-diagonalizable | **Degenerate/improper node** |

!!! theorem "Theorem 26.7 (Trace-determinant classification)"
    Let $\tau = \operatorname{tr}(A)$, $\delta = \det(A)$, $\Delta = \tau^2 - 4\delta$. Then:

    - $\delta < 0$: **Saddle point** (unstable).
    - $\delta > 0$, $\Delta > 0$, $\tau < 0$: **Stable node**.
    - $\delta > 0$, $\Delta > 0$, $\tau > 0$: **Unstable node**.
    - $\delta > 0$, $\Delta < 0$, $\tau < 0$: **Stable spiral**.
    - $\delta > 0$, $\Delta < 0$, $\tau > 0$: **Unstable spiral**.
    - $\delta > 0$, $\tau = 0$: **Center**.

??? proof "Proof"
    The eigenvalues are $\lambda_{1,2} = \frac{\tau \pm \sqrt{\Delta}}{2}$, with $\lambda_1 \lambda_2 = \delta$ and $\lambda_1 + \lambda_2 = \tau$.

    - $\delta < 0$ means $\lambda_1, \lambda_2$ have opposite signs.
    - $\delta > 0$, $\Delta > 0$: two real roots of the same sign, with the sign of $\tau$ determining stability.
    - $\delta > 0$, $\Delta < 0$: complex conjugate roots $\alpha \pm \beta i$ ($\alpha = \tau/2$), with the sign of $\alpha$ determining stability.
    - $\tau = 0$ implies $\alpha = 0$, giving purely imaginary eigenvalues. $\blacksquare$

!!! example "Example 26.5"
    **Phase plane of a spring-damper system.** The spring-mass system $mx'' + cx' + kx = 0$ becomes

    $$
    \begin{pmatrix} x' \\ v' \end{pmatrix} = \begin{pmatrix} 0 & 1 \\ -k/m & -c/m \end{pmatrix} \begin{pmatrix} x \\ v \end{pmatrix}.
    $$

    Here $\tau = -c/m < 0$, $\delta = k/m > 0$, $\Delta = c^2/m^2 - 4k/m$.

    - **Underdamped** ($c^2 < 4mk$): $\Delta < 0$, stable spiral. Trajectories spiral toward the origin.
    - **Overdamped** ($c^2 > 4mk$): $\Delta > 0$, stable node. Trajectories approach the origin monotonically along eigenvector directions.
    - **Undamped** ($c = 0$): $\tau = 0$, center. Trajectories are ellipses centered at the origin, corresponding to energy conservation.

---

## 26.5 Higher-Order Linear ODEs and the Companion Matrix

<div class="context-flow" markdown>

**Reduction**: $n$-th order linear ODE → first-order $n \times n$ system → the companion matrix $C$'s characteristic polynomial = the original equation's characteristic polynomial
**Link**: The characteristic polynomial theory of Ch6 is directly applied here

</div>

Higher-order linear ODEs can be converted to first-order systems by introducing state variables, with the corresponding coefficient matrix having a special companion matrix structure.

!!! definition "Definition 26.7 (Higher-order linear ODE)"
    An $n$-th order linear constant-coefficient ODE is

    $$
    y^{(n)} + a_{n-1}y^{(n-1)} + \cdots + a_1 y' + a_0 y = 0.
    $$

!!! definition "Definition 26.8 (Companion matrix)"
    The **companion matrix** corresponding to the above equation is

    $$
    C = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1} \end{pmatrix} \in \mathbb{R}^{n \times n}.
    $$

    The characteristic polynomial of $C$ is exactly $\lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0$.

!!! theorem "Theorem 26.8 (Equivalence of higher-order ODEs and first-order systems)"
    Let $x_1 = y$, $x_2 = y'$, $\ldots$, $x_n = y^{(n-1)}$. Then the higher-order ODE is equivalent to the first-order system

    $$
    \mathbf{x}' = C\mathbf{x},
    $$

    where $\mathbf{x} = (x_1, \ldots, x_n)^T$ and $C$ is the companion matrix. The eigenvalues of $C$ are exactly the characteristic roots of the original equation, and the Jordan structure of $C$ determines the form of the solution.

??? proof "Proof"
    By definition, $x_i' = x_{i+1}$ ($1 \le i \le n-1$), and

    $$
    x_n' = y^{(n)} = -a_{n-1}y^{(n-1)} - \cdots - a_0 y = -a_{n-1}x_n - \cdots - a_0 x_1.
    $$

    This is precisely the form $\mathbf{x}' = C\mathbf{x}$. For $C$'s characteristic polynomial: expanding the determinant $\det(\lambda I - C)$ along the last row and using the structure of $1$'s above the diagonal, one obtains $\lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_0$ through successive simplification. $\blacksquare$

!!! example "Example 26.6"
    **A fourth-order equation as a system.** The equation $y^{(4)} + 2y''' + 3y'' + 2y' + y = 0$ corresponds to the companion matrix

    $$
    C = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & -2 & -3 & -2 \end{pmatrix}.
    $$

    The characteristic polynomial is $\lambda^4 + 2\lambda^3 + 3\lambda^2 + 2\lambda + 1 = (\lambda^2 + \lambda + 1)^2$. The eigenvalues are $\lambda = \frac{-1 \pm i\sqrt{3}}{2}$ (each with multiplicity 2). The solutions include $e^{-t/2}\cos(\frac{\sqrt{3}}{2}t)$, $e^{-t/2}\sin(\frac{\sqrt{3}}{2}t)$, and their products with $t$.

---

## 26.6 Nonhomogeneous Equations and Variation of Parameters

<div class="context-flow" markdown>

**General solution of nonhomogeneous** = general homogeneous solution + particular solution → particular solution given by convolution integral with $e^{At}$ (variation of parameters) → essentially a Green's function / impulse response
**Link**: A direct manifestation of Ch4 quotient spaces and affine subspaces

</div>

The solution of nonhomogeneous linear systems $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$ relies on the matrix exponential and variation of parameters.

!!! definition "Definition 26.9 (Nonhomogeneous linear ODE system)"
    A **nonhomogeneous linear ODE system** is

    $$
    \mathbf{x}'(t) = A\mathbf{x}(t) + \mathbf{f}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$

    where $\mathbf{f}(t) \in \mathbb{R}^n$ is a given continuous forcing term.

!!! theorem "Theorem 26.9 (Variation of parameters formula)"
    The solution of the nonhomogeneous system $\mathbf{x}' = A\mathbf{x} + \mathbf{f}(t)$, $\mathbf{x}(0) = \mathbf{x}_0$ is

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\, ds.
    $$

    The first term is the homogeneous solution (free response), and the second is a particular solution (forced response) — the convolution of the matrix exponential with the forcing term.

??? proof "Proof"
    Set $\mathbf{x}(t) = e^{At}\mathbf{c}(t)$ (the core Ansatz of variation of parameters). Substituting into the equation:

    $$
    Ae^{At}\mathbf{c}(t) + e^{At}\mathbf{c}'(t) = Ae^{At}\mathbf{c}(t) + \mathbf{f}(t).
    $$

    Canceling $Ae^{At}\mathbf{c}(t)$ yields $e^{At}\mathbf{c}'(t) = \mathbf{f}(t)$, i.e., $\mathbf{c}'(t) = e^{-At}\mathbf{f}(t)$. Integrating:

    $$
    \mathbf{c}(t) = \mathbf{x}_0 + \int_0^t e^{-As}\mathbf{f}(s)\, ds.
    $$

    Therefore $\mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\, ds$. $\blacksquare$

!!! example "Example 26.7"
    **A forced oscillation system.** Consider the system

    $$
    \mathbf{x}' = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}\mathbf{x} + \begin{pmatrix} 0 \\ \cos t \end{pmatrix}, \quad \mathbf{x}(0) = \mathbf{0}.
    $$

    The homogeneous part has eigenvalues $\pm i$ (center), with $e^{At} = \begin{pmatrix} \cos t & \sin t \\ -\sin t & \cos t \end{pmatrix}$.

    The forced response is

    $$
    \int_0^t \begin{pmatrix} \cos(t-s) & \sin(t-s) \\ -\sin(t-s) & \cos(t-s) \end{pmatrix} \begin{pmatrix} 0 \\ \cos s \end{pmatrix} ds.
    $$

    Computing the first component: $\int_0^t \sin(t-s)\cos s\, ds = \frac{t}{2}\sin t$ (resonance). Since the driving frequency matches the natural frequency, the amplitude grows linearly — this is the linear algebra explanation of the resonance phenomenon.

---

## 26.7 Eigenvalue Problems in Linear PDEs

<div class="context-flow" markdown>

**Separation of variables** → spatial part becomes eigenvalue problem $\mathcal{L}u = \lambda u$ → Sturm-Liouville theory guarantees real eigenvalues, orthogonal eigenfunctions → spectral theorem on infinite-dimensional inner product spaces
**Link**: An infinite-dimensional generalization of the Ch8 spectral theorem for inner product spaces

</div>

The method of separation of variables in partial differential equations converts PDEs into ODE eigenvalue problems. Sturm-Liouville theory is the infinite-dimensional generalization of the finite-dimensional spectral theorem.

!!! definition "Definition 26.10 (Sturm-Liouville problem)"
    The **Sturm-Liouville problem** asks for $\lambda$ and nonzero functions $y(x)$ ($x \in [a, b]$) such that

    $$
    -(p(x)y')' + q(x)y = \lambda w(x) y,
    $$

    subject to boundary conditions (e.g., $y(a) = y(b) = 0$), where $p(x) > 0$ and $w(x) > 0$ is the weight function.

!!! theorem "Theorem 26.10 (Sturm-Liouville spectral theorem)"
    Under appropriate regularity conditions, the Sturm-Liouville problem has the following properties:

    1. There exist countably infinitely many real eigenvalues $\lambda_1 < \lambda_2 < \cdots$, with $\lambda_n \to \infty$.
    2. The corresponding eigenfunctions $\{y_n(x)\}$ are orthogonal with respect to the weight function $w(x)$:

        $$
        \int_a^b y_m(x) y_n(x) w(x)\, dx = 0, \quad m \neq n.
        $$

    3. $\{y_n\}$ forms a complete orthogonal basis for $L^2_w([a, b])$ (analogous to the orthogonal eigenbasis of a finite-dimensional symmetric matrix).

??? proof "Proof"
    **Self-adjointness**: Define the operator $\mathcal{L}y = \frac{1}{w}[-(py')' + qy]$. Through integration by parts and the boundary conditions, one can verify that $\langle \mathcal{L}y, z \rangle_w = \langle y, \mathcal{L}z \rangle_w$, i.e., $\mathcal{L}$ is self-adjoint with respect to the $L^2_w$ inner product.

    **Orthogonality**: If $\mathcal{L}y_m = \lambda_m y_m$ and $\mathcal{L}y_n = \lambda_n y_n$ ($\lambda_m \neq \lambda_n$), then

    $$
    \lambda_m \langle y_m, y_n \rangle_w = \langle \mathcal{L}y_m, y_n \rangle_w = \langle y_m, \mathcal{L}y_n \rangle_w = \lambda_n \langle y_m, y_n \rangle_w.
    $$

    Therefore $(\lambda_m - \lambda_n)\langle y_m, y_n \rangle_w = 0$, i.e., $\langle y_m, y_n \rangle_w = 0$. This proof is exactly parallel to showing that eigenvectors of a symmetric matrix corresponding to distinct eigenvalues are orthogonal. $\blacksquare$

!!! example "Example 26.8"
    **Separation of variables for the heat equation.** Consider the one-dimensional heat equation

    $$
    u_t = u_{xx}, \quad u(0, t) = u(\pi, t) = 0, \quad u(x, 0) = f(x).
    $$

    Setting $u(x, t) = X(x)T(t)$, separation of variables gives $T'/T = X''/X = -\lambda$.

    **Spatial problem**: $X'' = -\lambda X$, $X(0) = X(\pi) = 0$. This is a Sturm-Liouville problem with $p = w = 1$, $q = 0$. Eigenvalues $\lambda_n = n^2$, eigenfunctions $X_n = \sin(nx)$ ($n = 1, 2, 3, \ldots$).

    **Temporal problem**: $T_n' = -n^2 T_n$, solution $T_n = e^{-n^2 t}$.

    **Complete solution**:

    $$
    u(x, t) = \sum_{n=1}^{\infty} b_n e^{-n^2 t} \sin(nx), \quad b_n = \frac{2}{\pi}\int_0^{\pi} f(x)\sin(nx)\, dx.
    $$

    Each mode $n$ decays at rate $n^2$ — higher-frequency components decay faster, explaining the smoothing effect of heat conduction. The eigenvalues $n^2$ completely determine the time scales of diffusion.

---

## 26.8 Floquet Theory

<div class="context-flow" markdown>

**Problem**: $\mathbf{x}' = A(t)\mathbf{x}$, $A(t+T) = A(t)$ → solutions are no longer $e^{At}$ → Floquet's theorem: $\Phi(t) = P(t)e^{Bt}$ (periodic part $\times$ exponential part) → stability determined by eigenvalues of $B$
**Applications**: Parametric resonance (Mathieu equation), Bloch's theorem in crystal lattices

</div>

When the coefficient matrix has periodicity $A(t+T) = A(t)$, Floquet theory reduces the periodic system to a constant-coefficient problem.

!!! definition "Definition 26.11 (Monodromy matrix)"
    For the periodic system $\mathbf{x}' = A(t)\mathbf{x}$ ($A(t+T) = A(t)$), let $\Phi(t)$ be the fundamental matrix satisfying $\Phi(0) = I$. The **monodromy matrix** is defined as

    $$
    M = \Phi(T).
    $$

    The eigenvalues of $M$ are called **Floquet multipliers**.

!!! definition "Definition 26.12 (Floquet exponents)"
    If $\rho$ is a Floquet multiplier, then $\mu = \frac{1}{T}\ln \rho$ (with appropriate branch) is called a **Floquet exponent**. These are the "effective eigenvalues" of the equivalent constant-coefficient system.

!!! theorem "Theorem 26.11 (Floquet's theorem)"
    Let $A(t)$ be a continuous $T$-periodic matrix. The fundamental matrix $\Phi(t)$ ($\Phi(0) = I$) can be written as

    $$
    \Phi(t) = P(t) e^{Bt},
    $$

    where $P(t+T) = P(t)$ is a $T$-periodic nonsingular matrix and $B$ is a constant matrix satisfying $e^{BT} = M = \Phi(T)$.

??? proof "Proof"
    **Key observation**: $\Psi(t) = \Phi(t+T)$ is also a fundamental matrix, since

    $$
    \Psi'(t) = \Phi'(t+T) = A(t+T)\Phi(t+T) = A(t)\Psi(t).
    $$

    By uniqueness of fundamental matrices (up to right multiplication by a constant matrix), $\Phi(t+T) = \Phi(t)M$.

    Choose $B$ such that $e^{BT} = M$ (this always exists — take the matrix logarithm). Define $P(t) = \Phi(t)e^{-Bt}$. Then

    $$
    P(t+T) = \Phi(t+T)e^{-B(t+T)} = \Phi(t)Me^{-BT}e^{-Bt} = \Phi(t)e^{-Bt} = P(t).
    $$

    Therefore $P(t)$ is $T$-periodic. $\blacksquare$

!!! theorem "Theorem 26.12 (Floquet stability criterion)"
    The periodic system $\mathbf{x}' = A(t)\mathbf{x}$ is asymptotically stable at the equilibrium if and only if all Floquet multipliers satisfy $|\rho_i| < 1$ (equivalently, all Floquet exponents satisfy $\operatorname{Re}(\mu_i) < 0$).

??? proof "Proof"
    From $\Phi(nT) = M^n$ and $\Phi(t) = P(t)e^{Bt}$, we have $\|\Phi(nT)\| = \|M^n\|$. $M^n \to 0$ if and only if all eigenvalues (Floquet multipliers) of $M$ have modulus less than $1$. Between periods, $P(t)$ is bounded, so $\|\Phi(t)\| \to 0$ if and only if $|\rho_i| < 1$. $\blacksquare$

!!! example "Example 26.9"
    **The Mathieu equation and parametric resonance.** The Mathieu equation

    $$
    y'' + (\delta + \varepsilon \cos 2t) y = 0
    $$

    becomes the system $\mathbf{x}' = A(t)\mathbf{x}$ with $A(t) = \begin{pmatrix} 0 & 1 \\ -\delta - \varepsilon\cos 2t & 0 \end{pmatrix}$, period $T = \pi$.

    In the parameter plane $(\delta, \varepsilon)$, there exist **stability regions** ($|\rho_i| \le 1$) and **instability regions/resonance tongues** ($|\rho_i| > 1$). When $\varepsilon = 0$, the instability regions degenerate to points $\delta = n^2$ ($n = 1, 2, \ldots$). When $\varepsilon > 0$, the instability regions expand into "tongue"-shaped domains.

    Physical application: An inverted pendulum can be stabilized under appropriate frequency and amplitude of vertical oscillation (the Kapitza pendulum) — Floquet analysis gives the precise stability conditions.

!!! example "Example 26.10"
    **Floquet theory and Bloch's theorem.** In solid-state physics, the Schrodinger equation for an electron in a lattice

    $$
    -\psi''(x) + V(x)\psi(x) = E\psi(x), \quad V(x+a) = V(x)
    $$

    is a periodic-coefficient ODE. Floquet's theorem guarantees the existence of solutions in the Bloch wave form

    $$
    \psi(x) = e^{ikx} u(x), \quad u(x+a) = u(x),
    $$

    where $k$ is the wave number (imaginary part of the Floquet exponent). The energies $E(k)$ for different values of $k$ form the **band structure**, with band gaps corresponding to boundaries where the Floquet multipliers have unit modulus.
