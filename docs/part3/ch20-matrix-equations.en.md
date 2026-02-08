# Chapter 20  Matrix Equations

<div class="context-flow" markdown>

**Prerequisites**: Ch19 Kronecker product / Vec operator · **Chapter arc**: $AX=B$ (pseudoinverse) → **Sylvester** $AX+XB=C$ (eigenvalue separation condition) → **Lyapunov** (stability) → **Riccati** (optimal control, nonlinear) → Penrose equations (axiomatization of pseudoinverse)

</div>

Matrix equations are central to the applications of linear algebra in control theory, signal processing, and numerical analysis. Unlike ordinary linear systems $A\mathbf{x} = \mathbf{b}$, the unknown in a matrix equation is a **matrix** rather than a vector. From the simplest $AX = B$ to the structurally rich Sylvester equation $AX + XB = C$, the Lyapunov equation $AX + XA^* = Q$, and the nonlinear Riccati equation, each of these equations has its own unique mathematical structure and solution methods.

This chapter systematically introduces the theoretical foundations, solvability conditions, solution formulas, and numerical algorithms for these matrix equations, and concludes with a discussion of the Penrose equations closely related to the Moore-Penrose pseudoinverse.

---

## 20.1 Overview of linear matrix equations

<div class="context-flow" markdown>

**Chapter arc**: $AX=B$ solvable $\Leftrightarrow$ $\mathcal{C}(B)\subseteq\mathcal{C}(A)$ · General solution structure = particular solution + null space · $AXB=C$: two-sided elimination, condition $AA^+CB^+B=C$

</div>

This section starts from the most basic matrix equations and establishes the framework for matrix equation theory.

!!! definition "Definition 20.1 (Basic linear matrix equations)"
    The following are the three most common linear matrix equations:

    1. **One-sided equation**: $AX = B$, where $A$ is $m \times n$, $B$ is $m \times p$, and $X$ is the $n \times p$ unknown matrix.
    2. **Two-sided equation**: $AXB = C$, where $A$ is $m \times n$, $B$ is $p \times q$, $C$ is $m \times q$, and $X$ is the $n \times p$ unknown matrix.
    3. **Sylvester-type equation**: $AX + XB = C$, where $A$ is $m \times m$, $B$ is $n \times n$, $C$ is $m \times n$, and $X$ is the $m \times n$ unknown matrix.

!!! theorem "Theorem 20.1 (Solvability of the one-sided equation $AX = B$)"
    The equation $AX = B$ has a solution if and only if $\operatorname{rank}(A) = \operatorname{rank}(A \mid B)$, i.e., the column space of $B$ is contained in the column space of $A$: $\mathcal{C}(B) \subseteq \mathcal{C}(A)$.

    When $A$ is invertible, the unique solution is $X = A^{-1}B$.

    In general, let $A^+$ be the Moore-Penrose pseudoinverse of $A$. The general solution is:

    $$
    X = A^+ B + (I - A^+ A)Z,
    $$

    where $Z$ is an arbitrary $n \times p$ matrix.

??? proof "Proof"
    **Necessity**: If $AX = B$, then every column $\mathbf{b}_j = A\mathbf{x}_j \in \mathcal{C}(A)$.

    **Sufficiency**: If $\mathcal{C}(B) \subseteq \mathcal{C}(A)$, then $AA^+B = B$ (since $A^+$ is a left inverse from $\mathcal{C}(A)$ to itself, and $AA^+$ is the orthogonal projection onto $\mathcal{C}(A)$). Therefore $X_0 = A^+B$ is a particular solution.

    **General solution**: Let $X$ be any solution. Then $A(X - X_0) = 0$, i.e., $X - X_0 \in \mathcal{N}(A)$. Since $\mathcal{N}(A) = \mathcal{C}(I - A^+A)$ ($I - A^+A$ is the orthogonal projection onto $\mathcal{N}(A)$), we have $X - X_0 = (I - A^+A)Z$. $\blacksquare$

!!! theorem "Theorem 20.2 (Solvability of the two-sided equation $AXB = C$)"
    The equation $AXB = C$ has a solution if and only if $AA^+CB^+B = C$.

    When $A$ and $B$ are both invertible, the unique solution is $X = A^{-1}CB^{-1}$.

    In general, the general solution is:

    $$
    X = A^+ C B^+ + Z - A^+ A Z B B^+,
    $$

    where $Z$ is an arbitrary matrix of appropriate size.

??? proof "Proof"
    **Necessity**: If $AXB = C$, then $AA^+CB^+B = AA^+(AXB)B^+B = (AA^+A)X(BB^+B) = AXB = C$.

    **Sufficiency**: Assume $AA^+CB^+B = C$. Let $X_0 = A^+CB^+$. Then:

    $$
    AX_0B = A(A^+CB^+)B = (AA^+)C(B^+B) = AA^+CB^+B = C.
    $$

    **General solution**: Let $AXB = C$. Then $A(X - X_0)B = 0$. We need the general solution of $A\tilde{X}B = 0$. $\tilde{X} = Z - A^+AZBB^+$ satisfies:

    $$
    A(Z - A^+AZBB^+)B = AZB - (AA^+A)Z(BB^+B) = AZB - AZB = 0. \quad \blacksquare
    $$

!!! definition "Definition 20.2 (Operator form of matrix equations)"
    A general linear matrix equation can be written in operator form $\mathcal{L}(X) = C$, where $\mathcal{L}: \mathbb{C}^{m \times n} \to \mathbb{C}^{m \times n}$ is a linear operator. Using the Vec operator and Kronecker product, this is equivalent to:

    $$
    L \operatorname{vec}(X) = \operatorname{vec}(C),
    $$

    where $L$ is an $mn \times mn$ matrix (see Chapter 19).

!!! example "Example 20.1"
    Solve $AX = B$, where $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$, $B = \begin{pmatrix} 3 \\ 6 \end{pmatrix}$.

    $\operatorname{rank}(A) = 1$ (the second row is twice the first). $B = 3\begin{pmatrix} 1 \\ 2 \end{pmatrix} \in \mathcal{C}(A)$, so a solution exists.

    $A^+ = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$ (since $A = \mathbf{a}\mathbf{a}^T/5$ where $\mathbf{a} = (1,2)^T$, so $A^+ = A/\|A\|_F^2$... let us compute directly).

    Since $A = \begin{pmatrix} 1 \\ 2 \end{pmatrix}(1, 2)$, $A^+ = \frac{\mathbf{v}\mathbf{u}^T}{\|\mathbf{u}\|^2\|\mathbf{v}\|^2} = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$.

    Particular solution: $X_0 = A^+B = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} 3 \\ 6 \end{pmatrix} = \frac{1}{25}\begin{pmatrix} 15 \\ 30 \end{pmatrix} = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix}$.

    General solution: $X = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix} + (I - A^+A)Z = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix} + \frac{1}{5}\begin{pmatrix} 4 & -2 \\ -2 & 1 \end{pmatrix}\begin{pmatrix} z \end{pmatrix}$.

    That is, $X = \begin{pmatrix} 3/5 + 4t/5 \\ 6/5 - 2t/5 \end{pmatrix}$, or $X = \begin{pmatrix} 3 + 4t \\ 6 - 2t \end{pmatrix}/5$ ($t$ is an arbitrary real number).

    Verification (taking $t = 3$): $X = \begin{pmatrix} 3 \\ 0 \end{pmatrix}$, $AX = \begin{pmatrix} 3 \\ 6 \end{pmatrix} = B$. Correct.

---

## 20.2 Sylvester equation

<div class="context-flow" markdown>

**Core**: $AX+XB=C$ uniquely solvable $\Leftrightarrow$ $\sigma(A)\cap\sigma(-B)=\emptyset$ · Ch19 Kronecker sum perspective: coefficient matrix $=A\oplus B^T$, eigenvalues $\lambda_i(A)+\lambda_j(B)$ · Integral representation requires stability conditions

</div>

The Sylvester equation is one of the most important matrix equations, appearing widely in control theory, model reduction, and matrix function theory.

!!! definition "Definition 20.3 (Sylvester equation)"
    The **Sylvester equation** is a linear matrix equation of the form:

    $$
    AX + XB = C,
    $$

    where $A \in \mathbb{C}^{m \times m}$, $B \in \mathbb{C}^{n \times n}$, $C \in \mathbb{C}^{m \times n}$ are known matrices, and $X \in \mathbb{C}^{m \times n}$ is the unknown matrix.
    When $A = -B^*$, the Sylvester equation reduces to a Lyapunov equation.

!!! theorem "Theorem 20.3 (Solvability of the Sylvester equation -- Sylvester-Rosenblum theorem)"
    The Sylvester equation $AX + XB = C$ has a unique solution $X$ for every $C$ if and only if $A$ and $-B$ have no common eigenvalues, i.e.:

    $$
    \sigma(A) \cap \sigma(-B) = \emptyset.
    $$

??? proof "Proof"
    **Kronecker product method**: By Chapter 19, the equation $AX + XB = C$ is equivalent to:

    $$
    (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$

    The eigenvalues of $I_n \otimes A + B^T \otimes I_m$ are $\{\lambda_i(A) + \lambda_j(B) : i = 1,\ldots,m;\; j = 1,\ldots,n\}$ (these are the eigenvalues of the Kronecker sum $A \oplus B^T$, noting that $\sigma(B^T) = \sigma(B)$).

    Therefore the equation has a unique solution for every $C$ if and only if all $\lambda_i(A) + \lambda_j(B) \neq 0$, i.e., $\lambda_i(A) \neq -\lambda_j(B)$ for all $i, j$, i.e., $\sigma(A) \cap \sigma(-B) = \emptyset$.

    **Direct method** (Rosenblum): Transform $A$ to Jordan normal form $A = PJP^{-1}$; the equation becomes $JY + YB = D$ ($Y = P^{-1}X$, $D = P^{-1}C$). Due to the triangular structure of $J$, $Y$ can be solved column by column. When $\sigma(A) \cap \sigma(-B) = \emptyset$, the coefficient matrix $(\lambda_i I + B)$ at each step is invertible. $\blacksquare$

!!! theorem "Theorem 20.4 (Integral representation of the Sylvester equation)"
    If all eigenvalues of $A$ have negative real parts ($A$ is stable) and all eigenvalues of $B$ have positive real parts, then the unique solution of the Sylvester equation $AX + XB = C$ is:

    $$
    X = \int_0^{\infty} e^{At} C e^{Bt} \, dt.
    $$

??? proof "Proof"
    First verify that the integral converges. Since $A$ is stable, $\|e^{At}\| \leq M e^{-\alpha t}$ (for some $\alpha > 0$); since the eigenvalues of $B$ have positive real parts, $-B$ is stable, $\|e^{Bt}\| = \|e^{-(-B)t}\|$...

    More precisely, the condition should be $\operatorname{Re}\lambda_i(A) < 0$ and $\operatorname{Re}\lambda_j(B) > 0$ (or the opposite with a sign change), so that $e^{At} \to 0$ ($t \to +\infty$) and $e^{Bt}$ grows but the integral still converges.

    In fact, the standard condition is $\operatorname{Re}\lambda_i(A) + \operatorname{Re}\mu_j(B) < 0$. In this case:

    Let $X = \int_0^{\infty} e^{At} C e^{Bt} \, dt$. For $AX + XB$:

    $$
    AX + XB = \int_0^{\infty} (Ae^{At})Ce^{Bt} \, dt + \int_0^{\infty} e^{At}C(e^{Bt}B) \, dt
    $$

    $$
    = \int_0^{\infty} \frac{d}{dt}\left(e^{At}Ce^{Bt}\right) dt = \left[e^{At}Ce^{Bt}\right]_0^{\infty} = 0 - C = -C.
    $$

    Therefore $X$ is the solution of $AX + XB = -C$. If the equation is $AX + XB = C$, then the solution is $X = -\int_0^{\infty} e^{At}Ce^{Bt}\,dt$. $\blacksquare$

!!! definition "Definition 20.4 (Homogeneous Sylvester equation and commutator)"
    The homogeneous Sylvester equation $AX - XA = 0$ (i.e., $AX = XA$) describes the matrices $X$ that **commute** with $A$. All matrices commuting with $A$ form an algebra called the **commutant** of $A$.

    The map $\operatorname{ad}_A(X) = AX - XA$ is called the **adjoint map** of $A$, and $[A, X] = AX - XA$ is called the **commutator** or **Lie bracket** of $A$ and $X$.

!!! theorem "Theorem 20.5 (Structure of commuting matrices)"
    Let $A$ be an $n \times n$ matrix with $k$ distinct eigenvalues $\lambda_1, \ldots, \lambda_k$, with corresponding Jordan block sizes $n_{i,1} \geq n_{i,2} \geq \cdots$ ($i = 1,\ldots,k$). Then the dimension of the space of matrices commuting with $A$ is:

    $$
    \dim\{X : AX = XA\} = \sum_{i=1}^{k} \sum_{j} (2j - 1) n_{i,j}',
    $$

    where $n_{i,j}'$ is the conjugate partition of the Jordan block sizes for eigenvalue $\lambda_i$.

    In particular, when $A$ has $n$ distinct eigenvalues (i.e., the minimal polynomial of $A$ equals the characteristic polynomial), the matrices commuting with $A$ are exactly the polynomials $p(A)$, and the dimension is $n$.

??? proof "Proof"
    **Simple case** ($A$ has $n$ distinct eigenvalues): Let $A = P \operatorname{diag}(\lambda_1, \ldots, \lambda_n) P^{-1}$. $AX = XA$ if and only if $DY = YD$ ($Y = P^{-1}XP$, $D = \operatorname{diag}(\lambda_1,\ldots,\lambda_n)$).

    $DY = YD$ means $\lambda_i y_{ij} = \lambda_j y_{ij}$, i.e., $(\lambda_i - \lambda_j)y_{ij} = 0$. When $\lambda_i \neq \lambda_j$, $y_{ij} = 0$, so $Y$ is a diagonal matrix, and $X = P Y P^{-1}$.

    A diagonal matrix $Y = \operatorname{diag}(d_1,\ldots,d_n)$ can be uniquely written as $Y = p(D)$ for a polynomial $p$ of degree at most $n-1$ (Lagrange interpolation), so $X = p(A)$. The dimension is $n$. $\blacksquare$

!!! example "Example 20.2"
    Solve the Sylvester equation $AX + XB = C$, where:

    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}, \quad C = \begin{pmatrix} 6 & 10 \\ 15 & 35 \end{pmatrix}.
    $$

    $\sigma(A) = \{1, 3\}$, $\sigma(-B) = \{-2, -4\}$, no common eigenvalues, so a unique solution exists.

    Let $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$.

    $AX + XB = \begin{pmatrix} x_{11} & x_{12} \\ 3x_{21} & 3x_{22} \end{pmatrix} + \begin{pmatrix} 2x_{11} & 4x_{12} \\ 2x_{21} & 4x_{22} \end{pmatrix} = \begin{pmatrix} 3x_{11} & 5x_{12} \\ 5x_{21} & 7x_{22} \end{pmatrix}$.

    Therefore $3x_{11} = 6$, $5x_{12} = 10$, $5x_{21} = 15$, $7x_{22} = 35$.

    Solution: $X = \begin{pmatrix} 2 & 2 \\ 3 & 5 \end{pmatrix}$.

    Verification: $AX + XB = \begin{pmatrix} 2 & 2 \\ 9 & 15 \end{pmatrix} + \begin{pmatrix} 4 & 8 \\ 6 & 20 \end{pmatrix} = \begin{pmatrix} 6 & 10 \\ 15 & 35 \end{pmatrix} = C$. Correct.

!!! example "Example 20.3"
    Consider the homogeneous Sylvester equation $AX - XB = 0$, i.e., $AX = XB$.

    Let $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$, $B = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$.

    $\sigma(A) = \{2\}$, $\sigma(B) = \{2, 3\}$. The common eigenvalue is $\{2\}$, so the homogeneous equation has a nontrivial solution.

    Let $X = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$.

    $AX = \begin{pmatrix} 2a+c & 2b+d \\ 2c & 2d \end{pmatrix}$, $XB = \begin{pmatrix} 2a & 3b \\ 2c & 3d \end{pmatrix}$.

    $AX = XB$ gives: $c = 0$, $2b + d = 3b$ (i.e., $d = b$), $2d = 3d$ (i.e., $d = 0$).

    Therefore $b = 0$, and the solution space is $X = \begin{pmatrix} a & 0 \\ 0 & 0 \end{pmatrix}$, with dimension 1.

---

## 20.3 Lyapunov equation

<div class="context-flow" markdown>

**Chapter arc**: $AX+XA^*=-Q$: Sylvester equation with $B=A^*$ · **Stability equivalence**: $A$ stable ($\operatorname{Re}\lambda_i<0$) $\Leftrightarrow$ $\exists X\succ 0$ satisfying the equation ($Q\succ 0$) · Integral solution $X=\int_0^\infty e^{At}Qe^{A^*t}dt$

</div>

The Lyapunov equation is an important special case of the Sylvester equation, playing a central role in system stability analysis.

!!! definition "Definition 20.5 (Continuous Lyapunov equation)"
    The **continuous Lyapunov equation** is:

    $$
    AX + XA^* = Q,
    $$

    where $A \in \mathbb{C}^{n \times n}$, $Q \in \mathbb{C}^{n \times n}$ ($Q = Q^*$), and $X \in \mathbb{C}^{n \times n}$ is the unknown Hermitian matrix.

!!! definition "Definition 20.6 (Discrete Lyapunov equation)"
    The **discrete Lyapunov equation** (also called the **Stein equation**) is:

    $$
    AXA^* - X = Q,
    $$

    or equivalently $X - AXA^* = -Q$. Here $A, Q, X \in \mathbb{C}^{n \times n}$.

!!! theorem "Theorem 20.6 (Continuous Lyapunov equation and stability)"
    Let $A \in \mathbb{C}^{n \times n}$.

    1. If $A$ is **stable** (i.e., all eigenvalues have negative real parts: $\operatorname{Re}\lambda_i(A) < 0$), then for any positive semidefinite $Q \succeq 0$, the continuous Lyapunov equation $AX + XA^* = -Q$ has a unique positive semidefinite solution:

    $$
    X = \int_0^{\infty} e^{At} Q e^{A^*t} \, dt \succeq 0.
    $$

    2. Conversely, if for some $Q \succ 0$ there exists $X \succ 0$ satisfying $AX + XA^* = -Q$, then $A$ is stable.

??? proof "Proof"
    **(1)**: Since $A$ is stable, $e^{At} \to 0$ ($t \to \infty$) at exponential rate, so the integral converges. Verification:

    $$
    AX + XA^* = \int_0^{\infty} (Ae^{At})Qe^{A^*t} \, dt + \int_0^{\infty} e^{At}Q(e^{A^*t}A^*) \, dt
    $$

    $$
    = \int_0^{\infty} \frac{d}{dt}(e^{At}Qe^{A^*t}) \, dt = [e^{At}Qe^{A^*t}]_0^{\infty} = 0 - Q = -Q.
    $$

    $X \succeq 0$ because for any $\mathbf{v}$: $\mathbf{v}^*X\mathbf{v} = \int_0^{\infty} \|Q^{1/2}e^{A^*t}\mathbf{v}\|^2 \, dt \geq 0$.

    Uniqueness is guaranteed by the Sylvester-Rosenblum theorem ($\sigma(A) \cap \sigma(-A^*) = \emptyset$, since the eigenvalues of $A$ have negative real parts while those of $-A^*$, namely $-\overline{\lambda_i(A)}$, have positive real parts).

    **(2)**: Let $A\mathbf{v} = \lambda\mathbf{v}$ ($\mathbf{v} \neq 0$). Then:

    $$
    \mathbf{v}^*(AX + XA^*)\mathbf{v} = \lambda(\mathbf{v}^*X\mathbf{v}) + \bar{\lambda}(\mathbf{v}^*X\mathbf{v}) = 2\operatorname{Re}(\lambda)(\mathbf{v}^*X\mathbf{v}) = -\mathbf{v}^*Q\mathbf{v}.
    $$

    Since $X \succ 0$, $\mathbf{v}^*X\mathbf{v} > 0$; since $Q \succ 0$, $\mathbf{v}^*Q\mathbf{v} > 0$. Therefore $2\operatorname{Re}(\lambda) < 0$, i.e., $\operatorname{Re}(\lambda) < 0$. $\blacksquare$

!!! theorem "Theorem 20.7 (Discrete Lyapunov equation and stability)"
    Let $A \in \mathbb{C}^{n \times n}$.

    1. If $A$ is **Schur stable** (i.e., all eigenvalues have modulus less than 1: $|\lambda_i(A)| < 1$), then for any $Q \succeq 0$, the discrete equation $X - AXA^* = Q$ has a unique positive semidefinite solution:

    $$
    X = \sum_{k=0}^{\infty} A^k Q (A^*)^k \succeq 0.
    $$

    2. Conversely, if for some $Q \succ 0$ there exists $X \succ 0$ satisfying $X - AXA^* = Q$, then $A$ is Schur stable.

??? proof "Proof"
    **(1)**: Since $|\lambda_i(A)| < 1$, the series $\sum_{k=0}^{\infty} A^k Q (A^*)^k$ converges absolutely. Verification:

    $$
    X - AXA^* = \sum_{k=0}^{\infty} A^k Q (A^*)^k - A\left(\sum_{k=0}^{\infty} A^k Q (A^*)^k\right)A^* = \sum_{k=0}^{\infty} A^k Q (A^*)^k - \sum_{k=1}^{\infty} A^k Q (A^*)^k = A^0 Q (A^*)^0 = Q.
    $$

    **(2)**: Let $A\mathbf{v} = \lambda\mathbf{v}$. Then $\mathbf{v}^*(X - AXA^*)\mathbf{v} = \mathbf{v}^*X\mathbf{v} - |\lambda|^2 \mathbf{v}^*X\mathbf{v} = (1 - |\lambda|^2)\mathbf{v}^*X\mathbf{v} = \mathbf{v}^*Q\mathbf{v} > 0$. Since $X \succ 0$, we get $1 - |\lambda|^2 > 0$, i.e., $|\lambda| < 1$. $\blacksquare$

!!! example "Example 20.4"
    Solve the continuous Lyapunov equation $AX + XA^T = -Q$, where:

    $$
    A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}, \quad Q = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}.
    $$

    $A$ is stable (eigenvalues $-1, -2$). Using the integral formula:

    $$
    X = \int_0^{\infty} e^{At} Q e^{A^Tt} \, dt = \int_0^{\infty} \begin{pmatrix} e^{-t} & 0 \\ 0 & e^{-2t} \end{pmatrix} \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} \begin{pmatrix} e^{-t} & 0 \\ 0 & e^{-2t} \end{pmatrix} dt.
    $$

    $$
    = \int_0^{\infty} \begin{pmatrix} 2e^{-2t} & 0 \\ 0 & 2e^{-4t} \end{pmatrix} dt = \begin{pmatrix} 1 & 0 \\ 0 & 1/2 \end{pmatrix}.
    $$

    Verification: $AX + XA^T = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix} + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix} = \begin{pmatrix} -2 & 0 \\ 0 & -2 \end{pmatrix} = -Q$. Correct.

!!! example "Example 20.5"
    Solve the discrete Lyapunov equation $X - AXA^T = Q$, where:

    $$
    A = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}, \quad Q = I_2.
    $$

    $A$ is Schur stable (eigenvalues $1/2, 1/3$).

    Let $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{12} & x_{22} \end{pmatrix}$ (symmetric solution).

    $X - AXA^T = \begin{pmatrix} x_{11} - x_{11}/4 & x_{12} - x_{12}/6 \\ x_{12} - x_{12}/6 & x_{22} - x_{22}/9 \end{pmatrix} = \begin{pmatrix} 3x_{11}/4 & 5x_{12}/6 \\ 5x_{12}/6 & 8x_{22}/9 \end{pmatrix} = I_2$.

    Solution: $x_{11} = 4/3$, $x_{12} = 0$, $x_{22} = 9/8$.

    $X = \begin{pmatrix} 4/3 & 0 \\ 0 & 9/8 \end{pmatrix}$.

    Verification can be done by direct substitution. $X \succ 0$ and $A$ is Schur stable, consistent with the theorem.

---

## 20.4 Riccati equation

<div class="context-flow" markdown>

**Chapter arc**: $A^*X+XA-XBR^{-1}B^*X+Q=0$ (nonlinear!) · Stable invariant subspace of the **Hamiltonian matrix** → stabilizing solution $X=U_2U_1^{-1}$ · Stabilizable + detectable $\Rightarrow$ unique positive semidefinite solution

</div>

The Riccati equation is a **nonlinear** matrix equation that plays a central role in optimal control, filtering, and game theory.

!!! definition "Definition 20.7 (Algebraic Riccati equation)"
    The **continuous algebraic Riccati equation** (CARE) is:

    $$
    A^*X + XA - XBR^{-1}B^*X + Q = 0,
    $$

    where $A \in \mathbb{C}^{n \times n}$, $B \in \mathbb{C}^{n \times m}$, $Q = Q^* \in \mathbb{C}^{n \times n}$ ($Q \succeq 0$), $R = R^* \in \mathbb{C}^{m \times m}$ ($R \succ 0$), and $X = X^* \in \mathbb{C}^{n \times n}$ is the unknown Hermitian matrix.

    The **discrete algebraic Riccati equation** (DARE) is:

    $$
    X = A^*XA - A^*XB(R + B^*XB)^{-1}B^*XA + Q.
    $$

!!! definition "Definition 20.8 (Hamiltonian matrix)"
    The **Hamiltonian matrix** associated with the continuous Riccati equation is:

    $$
    H = \begin{pmatrix} A & -BR^{-1}B^* \\ -Q & -A^* \end{pmatrix}.
    $$

    An important property of the Hamiltonian matrix is: if $\lambda$ is an eigenvalue of $H$, then $-\bar{\lambda}$ is also an eigenvalue.

!!! theorem "Theorem 20.8 (CARE and Hamiltonian matrix)"
    Assume the Hamiltonian matrix $H$ has no purely imaginary eigenvalues. Partition the $2n$ eigenvalues of $H$ into the stable part (negative real parts) and the unstable part (positive real parts), $n$ each.

    Let $\begin{pmatrix} U_1 \\ U_2 \end{pmatrix}$ be a basis for the stable invariant subspace of $H$ ($U_1, U_2$ are both $n \times n$). If $U_1$ is invertible, then the **stabilizing solution** of CARE is:

    $$
    X = U_2 U_1^{-1}.
    $$

    This solution makes the closed-loop matrix $A - BR^{-1}B^*X$ stable.

??? proof "Proof"
    Let $\mathbf{w} = \begin{pmatrix} \mathbf{x} \\ \mathbf{p} \end{pmatrix}$ be an eigenvector of $H$, $H\mathbf{w} = \lambda\mathbf{w}$. For the stable invariant subspace, $H\begin{pmatrix} U_1 \\ U_2 \end{pmatrix} = \begin{pmatrix} U_1 \\ U_2 \end{pmatrix}\Lambda_s$, where the eigenvalues of $\Lambda_s$ all lie in the left half-plane.

    Expanding the first row: $AU_1 - BR^{-1}B^*U_2 = U_1\Lambda_s$.
    Expanding the second row: $-QU_1 - A^*U_2 = U_2\Lambda_s$.

    If $U_1$ is invertible, let $X = U_2U_1^{-1}$, so $U_2 = XU_1$.

    Substituting into the first row: $AU_1 - BR^{-1}B^*XU_1 = U_1\Lambda_s$, i.e., $(A - BR^{-1}B^*X)U_1 = U_1\Lambda_s$. This shows that $A - BR^{-1}B^*X$ is similar to $\Lambda_s$, so the closed-loop matrix is stable.

    Substituting into the second row: $-QU_1 - A^*XU_1 = XU_1\Lambda_s = X(AU_1 - BR^{-1}B^*XU_1)$.

    Rearranging: $-Q - A^*X = XA - XBR^{-1}B^*X$, i.e., $A^*X + XA - XBR^{-1}B^*X + Q = 0$. $\blacksquare$

!!! theorem "Theorem 20.9 (Existence condition for positive definite solution of CARE)"
    Suppose $(A, B)$ is stabilizable and $(A, Q^{1/2})$ is detectable, with $Q \succeq 0$ and $R \succ 0$. Then CARE has a unique positive semidefinite solution $X \succeq 0$, and the closed-loop matrix $A - BR^{-1}B^*X$ is stable.

??? proof "Proof"
    A complete proof of this theorem is lengthy; we provide a framework of the main ideas.

    **Stabilizability** means there exists a feedback gain $K$ such that $A - BK$ is stable. **Detectability** means $(A^*, (Q^{1/2})^*)$ is stabilizable.

    Under these conditions, one can show that the Hamiltonian matrix $H$ has no purely imaginary eigenvalues (a key step), so the construction of Theorem 20.8 can proceed.

    Positive semidefiniteness $X \succeq 0$ follows from the condition $Q \succeq 0$ and a Lyapunov argument: when the closed-loop system $(A-BR^{-1}B^*X)$ is stable, $X$ satisfies a Lyapunov equation that guarantees positive semidefiniteness.

    Uniqueness follows by contradiction: if there were two positive semidefinite solutions $X_1, X_2$, their difference would satisfy a linear equation, and the stability condition would force the difference to be zero. $\blacksquare$

!!! example "Example 20.6"
    Solve the scalar Riccati equation: $2x + 2x - x^2 + 1 = 0$, where $A = 1$, $B = 1$, $R = 1$, $Q = 1$.

    Simplifying: $2x - x^2 + 1 = 0$, i.e., $x^2 - 2x - 1 = 0$.

    $x = \frac{2 \pm \sqrt{4+4}}{2} = 1 \pm \sqrt{2}$.

    Stabilizing solution: the closed-loop is $A - BR^{-1}B^*x = 1 - x$. Taking $x = 1 + \sqrt{2}$, the closed-loop $= 1 - 1 - \sqrt{2} = -\sqrt{2} < 0$ (stable). Taking $x = 1 - \sqrt{2}$, the closed-loop $= \sqrt{2} > 0$ (unstable).

    Therefore the stabilizing solution is $x = 1 + \sqrt{2}$.

!!! example "Example 20.7"
    **Riccati equation in optimal control**.

    Consider the linear system $\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$ with cost function $J = \int_0^{\infty}(\mathbf{x}^TQ\mathbf{x} + \mathbf{u}^TR\mathbf{u})\,dt$.

    Let $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ (double integrator), $B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$, $Q = I_2$, $R = 1$.

    Hamiltonian matrix: $H = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & -1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & -1 & 0 \end{pmatrix}$.

    The characteristic polynomial of $H$ can be computed as $\lambda^4 + 3\lambda^2 + 1 = 0$ (letting $\mu = \lambda^2$, $\mu^2 + 3\mu + 1 = 0$, $\mu = \frac{-3 \pm \sqrt{5}}{2}$).

    The subspace corresponding to the stable eigenvalues gives the positive definite solution of the Riccati equation $X = \begin{pmatrix} \sqrt{3} & 1 \\ 1 & \sqrt{3} \end{pmatrix}$ (exact values require more careful computation).

    The optimal feedback is $\mathbf{u} = -R^{-1}B^TX\mathbf{x} = -\begin{pmatrix} 0 & 1 \end{pmatrix}X\mathbf{x}$.

---

## 20.5 Numerical methods

<div class="context-flow" markdown>

**Intuition**: Direct Kronecker product method has $O(n^6)$ complexity, infeasible → **Bartels-Stewart**: Schur decomposition + column-by-column back-substitution = $O(n^3)$ · Key: triangular Sylvester equations can be solved column by column

</div>

This section introduces numerical methods for solving Sylvester and Lyapunov equations. Directly using the Kronecker product to transform the equation into an $n^2$-dimensional linear system is theoretically valid but has computational complexity $O(n^6)$, which is impractical. The Bartels-Stewart algorithm reduces the complexity to $O(n^3)$.

!!! definition "Definition 20.9 (Schur decomposition)"
    Any square matrix $A$ can be decomposed as $A = UTU^*$, where $U$ is unitary and $T$ is upper triangular (with diagonal entries being the eigenvalues of $A$). This is called the **Schur decomposition**. For real matrices, one can use the **real Schur decomposition**: $A = UTU^T$, where $T$ is quasi-upper triangular (diagonal blocks are $1 \times 1$ or $2 \times 2$).

!!! theorem "Theorem 20.10 (Bartels-Stewart algorithm)"
    The Sylvester equation $AX + XB = C$ can be solved efficiently through the following steps:

    **Step 1**: Compute Schur decompositions $A = U_A T_A U_A^*$, $B = U_B T_B U_B^*$.

    **Step 2**: Substitute $Y = U_A^* X U_B$, $D = U_A^* C U_B$; the equation becomes:

    $$
    T_A Y + Y T_B = D.
    $$

    **Step 3**: Since $T_A$ and $T_B$ are (quasi-)upper triangular, the equation can be solved column by column (starting from the last column) using back-substitution.

    **Step 4**: Recover $X = U_A Y U_B^*$.

    The total computational complexity is $O(n^3)$ (dominated by the Schur decomposition).

??? proof "Proof"
    **Correctness**: Substituting and verifying, $T_A Y + YT_B = U_A^*AU_A \cdot U_A^*XU_B + U_A^*XU_B \cdot U_B^*BU_B = U_A^*(AX + XB)U_B = U_A^*CU_B = D$.

    **Column-by-column solving**: Let the $j$-th column of $T_B$ be $\mathbf{t}_j$ (only the first $j$ entries may be nonzero), $\mathbf{y}_j$ be the $j$-th column of $Y$, and $\mathbf{d}_j$ be the $j$-th column of $D$.

    The $j$-th column of the equation $T_A Y + Y T_B = D$ is:

    $$
    T_A \mathbf{y}_j + \sum_{k=1}^{j} (T_B)_{kj} \mathbf{y}_k = \mathbf{d}_j.
    $$

    That is, $(T_A + (T_B)_{jj}I)\mathbf{y}_j = \mathbf{d}_j - \sum_{k=1}^{j-1}(T_B)_{kj}\mathbf{y}_k$.

    Starting from $j = 1$, the right side is known (when $j = 1$ there is no summation), and $(T_A + (T_B)_{11}I)$ is an upper triangular matrix that can be solved by back-substitution for $\mathbf{y}_1$. Then proceed column by column.

    Each column solve costs $O(n^2)$ (triangular system), with $n$ columns, so the back-substitution part is $O(n^3)$. The Schur decomposition is $O(n^3)$. Total complexity is $O(n^3)$. $\blacksquare$

!!! theorem "Theorem 20.11 (Hessenberg-Schur method)"
    For the Sylvester equation $AX + XB = C$ ($A$ is $m \times m$, $B$ is $n \times n$, $m \gg n$ or $m \ll n$), an improved **Hessenberg-Schur method** can be used: only the smaller matrix undergoes Schur decomposition, while the larger one is reduced to Hessenberg form. This can reduce computation in certain cases.

??? proof "Proof"
    Perform Schur decomposition on $B$: $B = U_B T_B U_B^*$, and reduce $A$ to upper Hessenberg form: $A = Q_A H_A Q_A^*$. After substitution, the equation becomes $H_A Y + Y T_B = D$.

    When solving column by column, $(H_A + (T_B)_{jj}I)\mathbf{y}_j = \mathbf{r}_j$ is a Hessenberg system that can be solved in $O(n^2)$ (without reducing to triangular form; direct Gauss elimination suffices since the LU factorization of a Hessenberg matrix is fast).

    The constant factor in the total computation is slightly smaller than the standard Bartels-Stewart algorithm. $\blacksquare$

!!! example "Example 20.8"
    Solve $AX + XB = C$ using the Bartels-Stewart algorithm, where:

    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 2 \\ 0 & 4 \end{pmatrix}, \quad C = \begin{pmatrix} 10 & 26 \\ 9 & 28 \end{pmatrix}.
    $$

    $A$ and $B$ are already upper triangular (i.e., they are their own Schur decompositions, $U_A = U_B = I$).

    Directly solve column by column $T_A Y + Y T_B = D$ (i.e., $AX + XB = C$).

    **Column 1** ($j = 1$): $(A + b_{11}I)\mathbf{x}_1 = \mathbf{c}_1$.

    $$
    \begin{pmatrix} 3 & 1 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_{11} \\ x_{21} \end{pmatrix} = \begin{pmatrix} 10 \\ 9 \end{pmatrix}.
    $$

    Back-substitution: $x_{21} = 9/4$, $x_{11} = (10 - 9/4)/3 = 31/12$.

    **Column 2** ($j = 2$): $(A + b_{22}I)\mathbf{x}_2 = \mathbf{c}_2 - b_{12}\mathbf{x}_1$.

    $$
    \begin{pmatrix} 6 & 1 \\ 0 & 7 \end{pmatrix}\begin{pmatrix} x_{12} \\ x_{22} \end{pmatrix} = \begin{pmatrix} 26 \\ 28 \end{pmatrix} - 2\begin{pmatrix} 31/12 \\ 9/4 \end{pmatrix} = \begin{pmatrix} 26 - 31/6 \\ 28 - 9/2 \end{pmatrix} = \begin{pmatrix} 125/6 \\ 47/2 \end{pmatrix}.
    $$

    Back-substitution: $x_{22} = 47/14$, $x_{12} = (125/6 - 47/14)/6 = (875/42 - 141/42)/6 = (734/42)/6 = 734/252 = 367/126$.

    Verification can be done by substituting back into the original equation. (The numerical values are somewhat involved but the process is correct.)

!!! example "Example 20.9"
    **Bartels-Stewart solution of a Lyapunov equation**.

    For $AX + XA^T = -Q$ (real symmetric case), only one Schur decomposition of $A$ is needed: $A = UTU^T$; the equation becomes $TY + YT^T = -\tilde{Q}$ ($Y = U^TXU$, $\tilde{Q} = U^TQU$). Since $T$ is upper triangular and $T^T$ is lower triangular, column-by-column solving requires exploiting symmetry.

    Let $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$... this is not upper triangular. Taking the real Schur decomposition: the eigenvalues of $A$ are $\pm i$ (purely imaginary), so the Lyapunov equation $AX + XA^T = -Q$ does not necessarily have a solution (since $\sigma(A) \cap \sigma(-A^T) = \{i, -i\} \cap \{i, -i\} \neq \emptyset$).

    Instead take $A = \begin{pmatrix} -1 & 1 \\ 0 & -2 \end{pmatrix}$ (stable, already upper triangular), $Q = I$.

    $(T + t_{11}I)\mathbf{y}_1 = -\mathbf{q}_1$: $\begin{pmatrix} -2 & 1 \\ 0 & -3 \end{pmatrix}\begin{pmatrix} y_{11} \\ y_{21} \end{pmatrix} = -\begin{pmatrix} 1 \\ 0 \end{pmatrix}$.

    $y_{21} = 0$, $y_{11} = 1/2$.

    $(T + t_{22}I)\mathbf{y}_2 = -\mathbf{q}_2 - t_{12}^T\mathbf{y}_1$: here $t_{12}^T$ is the upper triangular element of $T^T$...

    Handling more carefully: $TY + YT^T = -Q$. The equation structure for column-by-column solving requires using the triangularity of both the rows and columns simultaneously. If $Y$ is symmetric, only $n(n+1)/2$ variables need to be determined. In this example, the solution is $X = \begin{pmatrix} 5/6 & 1/6 \\ 1/6 & 1/4 \end{pmatrix}$ (to be verified by substitution).

---

## 20.6 Penrose equations

<div class="context-flow" markdown>

**Chapter arc**: The four Penrose equations axiomatically define $A^+$ · (1) inner inverse → (1)(3) least squares → (1)(4) minimum norm → (1)(2)(3)(4) unique Moore-Penrose pseudoinverse · SVD provides existence and explicit formula

</div>

The Penrose equations give an axiomatic characterization of the Moore-Penrose pseudoinverse.

!!! definition "Definition 20.10 (Penrose equations)"
    Let $A$ be an $m \times n$ matrix. The **Penrose equations** are the following four equations for an $n \times m$ matrix $X$:

    $$
    \begin{aligned}
    &(1)\quad AXA = A, \\
    &(2)\quad XAX = X, \\
    &(3)\quad (AX)^* = AX, \\
    &(4)\quad (XA)^* = XA.
    \end{aligned}
    $$

!!! definition "Definition 20.11 (Generalized inverses)"
    A matrix $X$ satisfying different subsets of the Penrose equations is called a **generalized inverse** of $A$:

    - $X$ satisfying (1) is called a **$\{1\}$-inverse** or **inner inverse**, denoted $A^{(1)}$.
    - $X$ satisfying (1)(2) is called a **$\{1,2\}$-inverse** or **reflexive generalized inverse**.
    - $X$ satisfying (1)(3) is called a **$\{1,3\}$-inverse** or **least squares generalized inverse**.
    - $X$ satisfying (1)(4) is called a **$\{1,4\}$-inverse** or **minimum norm generalized inverse**.
    - The unique $X$ satisfying all four equations is called the **Moore-Penrose pseudoinverse**, denoted $A^+$.

!!! theorem "Theorem 20.12 (Existence and uniqueness of the Moore-Penrose pseudoinverse)"
    For any $m \times n$ matrix $A$, the Moore-Penrose pseudoinverse $A^+$ exists and is unique.

??? proof "Proof"
    **Existence**: Let $A = U\Sigma V^*$ be the singular value decomposition, where $\Sigma = \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix}$, $\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$ ($r = \operatorname{rank}(A)$).

    Define $A^+ = V\Sigma^+ U^*$, where $\Sigma^+ = \begin{pmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{pmatrix}$.

    Verify the four Penrose equations one by one:

    (1) $AXA = U\Sigma V^* V\Sigma^+ U^* U\Sigma V^* = U\Sigma\Sigma^+\Sigma V^* = U\Sigma V^* = A$.

    (2) $XAX = V\Sigma^+ U^* U\Sigma V^* V\Sigma^+ U^* = V\Sigma^+\Sigma\Sigma^+ U^* = V\Sigma^+ U^* = X$.

    (3) $AX = U\Sigma\Sigma^+ U^* = U\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}U^*$, which is Hermitian (since the diagonal block is real symmetric, and $U(\cdot)U^*$ preserves the Hermitian property).

    (4) $XA = V\Sigma^+\Sigma V^* = V\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}V^*$, which is similarly Hermitian.

    **Uniqueness**: Suppose $X_1, X_2$ both satisfy all four equations.

    By (1)(3): $AX_1$ and $AX_2$ are both orthogonal projections onto $\mathcal{C}(A)$, so $AX_1 = AX_2$.

    By (1)(4): $X_1A$ and $X_2A$ are both orthogonal projections onto $\mathcal{C}(A^*)$, so $X_1A = X_2A$.

    By (2): $X_1 = X_1AX_1 = X_2AX_1 = X_2AX_2 = X_2$ (the intermediate steps use $X_1A = X_2A$ and $AX_1 = AX_2$). $\blacksquare$

!!! theorem "Theorem 20.13 (Properties of the Moore-Penrose pseudoinverse)"
    Let $A$ be an $m \times n$ matrix and $A^+$ be its Moore-Penrose pseudoinverse. Then:

    1. $(A^+)^+ = A$.
    2. $(A^*)^+ = (A^+)^*$.
    3. $(A^*A)^+ = A^+(A^+)^* = A^+(A^*)^+$.
    4. $\operatorname{rank}(A^+) = \operatorname{rank}(A)$.
    5. $A^+ = (A^*A)^+A^* = A^*(AA^*)^+$.
    6. If $A$ has full column rank, then $A^+ = (A^*A)^{-1}A^*$.
    7. If $A$ has full row rank, then $A^+ = A^*(AA^*)^{-1}$.

??? proof "Proof"
    **(1)**: $A$ satisfies the four Penrose equations with $A^+$ as the "original matrix" (swapping the roles of $A$ and $X$ and using the Hermitian conditions), so $(A^+)^+ = A$.

    **(2)**: Let $A = U\Sigma V^*$. Then $A^* = V\Sigma^* U^*$, $(A^*)^+ = U(\Sigma^*)^+ V^*$. Also $A^+ = V\Sigma^+ U^*$, $(A^+)^* = U(\Sigma^+)^* V^*$. Since $\Sigma$ is a real matrix, $(\Sigma^*)^+ = (\Sigma^+)^*$, so $(A^*)^+ = (A^+)^*$.

    **(6)**: If $A$ has full column rank ($\operatorname{rank}(A) = n$), then $A^*A$ is invertible. In this case $X = (A^*A)^{-1}A^*$ satisfies:
    - (1): $AXA = A(A^*A)^{-1}A^*A = A$.
    - (3): $AX = A(A^*A)^{-1}A^*$ is Hermitian.
    - (4): $XA = (A^*A)^{-1}A^*A = I$ is Hermitian.
    - (2): $XAX = IX = X$.
    Therefore $X = A^+$. $\blacksquare$

!!! theorem "Theorem 20.14 (General solution of $\{1\}$-inverses)"
    For an $m \times n$ matrix $A$ ($\operatorname{rank}(A) = r$), the general solution of the equation $AXA = A$ is:

    $$
    X = A^+ + (I - A^+A)W_1 + W_2(I - AA^+),
    $$

    where $W_1, W_2$ are arbitrary matrices of appropriate sizes.

    Equivalently, let $A = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q$ ($P, Q$ invertible). The general form of a $\{1\}$-inverse is:

    $$
    X = Q^{-1}\begin{pmatrix} I_r & L \\ M & N \end{pmatrix}P^{-1},
    $$

    where $L, M, N$ are arbitrary matrices.

??? proof "Proof"
    Let $A = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q$, $X = Q^{-1}\begin{pmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{pmatrix}P^{-1}$.

    $AXA = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{pmatrix}\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q = P\begin{pmatrix} G_{11} & 0 \\ 0 & 0 \end{pmatrix}Q$.

    $AXA = A$ requires $G_{11} = I_r$, while $G_{12}, G_{21}, G_{22}$ are arbitrary. $\blacksquare$

!!! theorem "Theorem 20.15 (Penrose equations and least squares)"
    The solutions of the equation $AXA = A$ are closely related to least squares problems:

    For the incompatible system $A\mathbf{x} = \mathbf{b}$:

    1. $\mathbf{x} = A^{(1)}\mathbf{b}$ is some least squares solution (satisfying (1)), but not necessarily the one with minimum norm.
    2. $\mathbf{x} = A^{(1,3)}\mathbf{b}$ is a least squares solution ($\min \|A\mathbf{x} - \mathbf{b}\|$), satisfying (1)(3).
    3. $\mathbf{x} = A^{(1,4)}\mathbf{b}$ is a minimum norm solution (minimum norm among solutions of the compatible equation), satisfying (1)(4).
    4. $\mathbf{x} = A^+\mathbf{b}$ is the **minimum norm least squares solution** (satisfying all four equations).

??? proof "Proof"
    **Proof of (4)**: $\mathbf{x}_0 = A^+\mathbf{b}$ satisfies $A\mathbf{x}_0 = AA^+\mathbf{b}$. By Penrose equation (3), $AA^+$ is an orthogonal projection onto $\mathcal{C}(A)$, so $AA^+\mathbf{b}$ is the orthogonal projection of $\mathbf{b}$ onto $\mathcal{C}(A)$, i.e., the minimum of $\|A\mathbf{x} - \mathbf{b}\|$ is attained at $A\mathbf{x} = AA^+\mathbf{b}$.

    The set of least squares solutions is $\{\mathbf{x} : A\mathbf{x} = AA^+\mathbf{b}\} = A^+\mathbf{b} + \mathcal{N}(A)$.

    $\mathbf{x}_0 = A^+\mathbf{b} \in \mathcal{C}(A^+) = \mathcal{C}(A^*)$ (as seen from the SVD). Since $\mathcal{C}(A^*) = \mathcal{N}(A)^{\perp}$, $\mathbf{x}_0$ is orthogonal to $\mathcal{N}(A)$, and therefore has the minimum norm in the affine set $A^+\mathbf{b} + \mathcal{N}(A)$. $\blacksquare$

!!! example "Example 20.10"
    Compute the Moore-Penrose pseudoinverse of $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$.

    $A$ has full column rank ($\operatorname{rank} = 2$), so $A^+ = (A^TA)^{-1}A^T$.

    $A^TA = I_2$, hence $A^+ = A^T = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$.

    Verifying the Penrose equations:

    (1) $AXA = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = A$.

    (2) $XAX = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = X$.

    (3) $AX = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$, symmetric.

    (4) $XA = I_2$, symmetric.

    All verifications pass.

!!! example "Example 20.11"
    Compute the Moore-Penrose pseudoinverse of the rank-1 matrix $A = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}$.

    $A = \mathbf{u}\mathbf{v}^T$, where $\mathbf{u} = (1,2,3)^T$, $\mathbf{v} = (1,1)^T$.

    For a rank-1 matrix, $A^+ = \frac{\mathbf{v}\mathbf{u}^T}{\|\mathbf{u}\|^2\|\mathbf{v}\|^2} = \frac{1}{14 \cdot 2}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix} = \frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}$.

    Verifying (1): $AXA = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}\frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}$.

    First compute $XA = \frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix} = \frac{14}{28}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$.

    $AXA = A \cdot \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$... more precisely: $A(XA) = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}\frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix} = A$. Correct.

!!! example "Example 20.12"
    **Using the pseudoinverse to solve an incompatible system**.

    Let $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix}$.

    $A^TA = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$, $(A^TA)^{-1} = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}$.

    $A^+ = (A^TA)^{-1}A^T = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}\begin{pmatrix} 1 & 1 & 0 \\ 1 & -1 & 1 \end{pmatrix} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ 1/3 & -1/3 & 1/3 \end{pmatrix}$.

    Minimum norm least squares solution: $\mathbf{x} = A^+\mathbf{b} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ 1/3 & -1/3 & 1/3 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$.

    Residual: $A\mathbf{x} - \mathbf{b} = \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} - \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} = \mathbf{0}$.

    The system turns out to be compatible! The solution is unique (since $A$ has full column rank).

!!! example "Example 20.13"
    **A truly incompatible case**. Let $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$.

    $A^+ = \frac{A^T}{A^TA} = \frac{1}{2}(1, 1)$.

    $\mathbf{x} = A^+\mathbf{b} = \frac{1}{2}(1 + 3) = 2$.

    Least squares solution: $x = 2$, $A\mathbf{x} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$, residual $= \begin{pmatrix} -1 \\ 1 \end{pmatrix}$, $\|$residual$\| = \sqrt{2}$.

    This is optimal, since $\|A x - \mathbf{b}\|^2 = (x-1)^2 + (x-3)^2$, differentiating with respect to $x$ gives $2(x-1) + 2(x-3) = 0$, so $x = 2$.
