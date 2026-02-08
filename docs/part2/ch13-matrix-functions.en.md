# Chapter 13  Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Ch12 Jordan normal form · Ch10 Spectral decomposition · **Chapter arc**: $p(A)$ (polynomials) → Cayley-Hamilton → Power series/convergence → $e^A$ (matrix exponential) → $\log A$, $A^{1/2}$ → General $f(A)$ (Jordan/Cauchy integral)
Essence: $f(\lambda) \to f(A)$ — the Jordan form lets scalar functions **act on matrices block by block**, with derivative information filling the superdiagonal

</div>

In the preceding chapters we have become familiar with basic matrix operations such as addition, multiplication, and inversion. This chapter extends the concept of functions from scalars to matrices: given a function $f$ (such as the exponential, logarithm, or square root), we wish to define what $f(A)$ means. Matrix functions have important applications in differential equations, control theory, quantum mechanics, and other fields. Starting from matrix polynomials and proceeding through power series and the Jordan normal form, this chapter systematically builds the complete theory of matrix functions.

---

## 13.1 Matrix Polynomials

<div class="context-flow" markdown>

$p(A) = a_k A^k + \cdots + a_0 I$: the most basic matrix function → Preserves similarity $p(PBP^{-1}) = Pp(B)P^{-1}$ → **Cayley-Hamilton**: $p_A(A) = 0$

</div>

Matrix polynomials are the most basic way to define matrix functions and the starting point for understanding more general matrix functions.

!!! definition "Definition 13.1 (Matrix Polynomial)"
    Let $p(\lambda) = a_k \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_1\lambda + a_0$ be a scalar polynomial and $A$ an $n \times n$ matrix. The **matrix polynomial** is defined as
    $$
    p(A) = a_k A^k + a_{k-1}A^{k-1} + \cdots + a_1 A + a_0 I,
    $$
    where $A^0 = I$ (the identity matrix).

!!! theorem "Theorem 13.1 (Basic properties of matrix polynomials)"
    Let $p, q$ be scalar polynomials and $A$ an $n \times n$ matrix. Then:

    1. $(p + q)(A) = p(A) + q(A)$;
    2. $(pq)(A) = p(A)q(A)$;
    3. $p(A)$ and $q(A)$ commute: $p(A)q(A) = q(A)p(A)$;
    4. If $A = PBP^{-1}$, then $p(A) = Pp(B)P^{-1}$;
    5. The eigenvalues of $p(A)$ are $p(\lambda_1), \ldots, p(\lambda_n)$, where $\lambda_i$ are the eigenvalues of $A$.

??? proof "Proof"
    **(1)** and **(2)** follow directly from the distributive and associative laws of matrix multiplication.

    **(3)** By (2), $p(A)q(A) = (pq)(A) = (qp)(A) = q(A)p(A)$ (using commutativity of scalar polynomial multiplication).

    **(4)** $A^k = PB^kP^{-1}$ (easily proved by induction), so
    $$
    p(A) = \sum a_i A^i = \sum a_i PB^iP^{-1} = P\left(\sum a_i B^i\right)P^{-1} = Pp(B)P^{-1}.
    $$

    **(5)** If $A\mathbf{v} = \lambda\mathbf{v}$, then $A^k\mathbf{v} = \lambda^k\mathbf{v}$, so $p(A)\mathbf{v} = p(\lambda)\mathbf{v}$. That is, $p(\lambda)$ is an eigenvalue of $p(A)$. $\blacksquare$

<div class="context-flow" markdown>

**Insight**: Cayley-Hamilton means $A^n$ can be expressed as a linear combination of $I, A, \ldots, A^{n-1}$ — the matrix algebra $\mathbb{F}[A]$ has dimension $\le n$

</div>

!!! theorem "Theorem 13.2 (Cayley-Hamilton Theorem)"
    Let $A$ be an $n \times n$ matrix with characteristic polynomial $p_A(\lambda) = \det(\lambda I - A)$. Then
    $$
    p_A(A) = 0.
    $$
    That is, every matrix satisfies its own characteristic equation.

??? proof "Proof"
    Let $p_A(\lambda) = \lambda^n + c_{n-1}\lambda^{n-1} + \cdots + c_0$.

    **Method 1 (Adjugate matrix method):** Let $B(\lambda) = \operatorname{adj}(\lambda I - A)$ be the adjugate of $\lambda I - A$. Then
    $$
    (\lambda I - A)B(\lambda) = \det(\lambda I - A) \cdot I = p_A(\lambda) I.
    $$
    Each entry of $B(\lambda)$ is a polynomial in $\lambda$ of degree at most $(n-1)$, so we can write
    $$
    B(\lambda) = B_{n-1}\lambda^{n-1} + B_{n-2}\lambda^{n-2} + \cdots + B_0,
    $$
    where $B_i$ are constant matrices. Expanding $(\lambda I - A)B(\lambda) = p_A(\lambda) I$ and comparing coefficients of each power of $\lambda$ yields a system of matrix equations. Left-multiplying the $k$-th equation by $A^k$ and summing all equations, the cancellations give $p_A(A) = 0$. $\blacksquare$

!!! example "Example 13.1"
    Verify the Cayley-Hamilton theorem for $A = \begin{pmatrix}1&2\\3&4\end{pmatrix}$.

    **Solution:** Characteristic polynomial:
    $$
    p_A(\lambda) = \lambda^2 - 5\lambda - 2.
    $$

    Compute $p_A(A) = A^2 - 5A - 2I$:
    $$
    A^2 = \begin{pmatrix}7&10\\15&22\end{pmatrix}, \quad 5A = \begin{pmatrix}5&10\\15&20\end{pmatrix}, \quad 2I = \begin{pmatrix}2&0\\0&2\end{pmatrix}.
    $$
    $$
    p_A(A) = \begin{pmatrix}7&10\\15&22\end{pmatrix} - \begin{pmatrix}5&10\\15&20\end{pmatrix} - \begin{pmatrix}2&0\\0&2\end{pmatrix} = \begin{pmatrix}0&0\\0&0\end{pmatrix}. \quad \checkmark
    $$

---

## 13.2 Matrix Power Series

<div class="context-flow" markdown>

$f(A) = \sum c_k A^k$ converges ↔ **spectral radius** $\rho(A) < R$ (radius of convergence) → Neumann series $(I-A)^{-1} = \sum A^k$ ($\rho(A)<1$)

</div>

Extending polynomials to power series requires introducing the concept of convergence for matrix series.

!!! definition "Definition 13.2 (Convergence of Matrix Series)"
    Let $\{A_k\}$ be a sequence of $m \times n$ matrices. The series $\sum_{k=0}^{\infty} A_k$ is said to **converge** if the scalar series at each entry position converges, i.e., for all $1 \le i \le m$, $1 \le j \le n$, $\sum_{k=0}^{\infty} [A_k]_{ij}$ converges. In this case, we define
    $$
    \sum_{k=0}^{\infty} A_k = \left[\sum_{k=0}^{\infty} [A_k]_{ij}\right]_{m \times n}.
    $$

!!! definition "Definition 13.3 (Spectral Radius)"
    The **spectral radius** of a matrix $A$ is defined as
    $$
    \rho(A) = \max\{|\lambda| : \lambda \text{ is an eigenvalue of } A\}.
    $$

!!! theorem "Theorem 13.3 (Convergence criterion for matrix power series)"
    Let $f(z) = \sum_{k=0}^{\infty} c_k z^k$ be a power series with radius of convergence $R$, and let $A$ be an $n \times n$ matrix. A sufficient condition for the matrix power series
    $$
    f(A) = \sum_{k=0}^{\infty} c_k A^k
    $$
    to converge is $\rho(A) < R$.

??? proof "Proof"
    Let $A = PJP^{-1}$, where $J$ is the Jordan normal form. Then $A^k = PJ^kP^{-1}$, so
    $$
    f(A) = P\left(\sum_{k=0}^{\infty} c_k J^k\right)P^{-1} = Pf(J)P^{-1}.
    $$
    Since $J$ is block diagonal, $f(J) = \operatorname{diag}(f(J_{n_1}(\lambda_1)), \ldots)$.

    For a Jordan block $J_m(\lambda)$, the $(p,q)$ entry ($q \ge p$) of $f(J_m(\lambda))$ is
    $$
    \frac{f^{(q-p)}(\lambda)}{(q-p)!} = \sum_{k=q-p}^{\infty} c_k \binom{k}{q-p} \lambda^{k-q+p}.
    $$
    When $|\lambda| < R$, all derivatives of $f$ at $\lambda$ converge, so the expression above converges. $\rho(A) < R$ ensures all eigenvalues $\lambda$ satisfy $|\lambda| < R$. $\blacksquare$

!!! theorem "Theorem 13.4 (Neumann Series)"
    Let $A$ be an $n \times n$ matrix. If $\rho(A) < 1$, then $I - A$ is invertible and
    $$
    (I - A)^{-1} = \sum_{k=0}^{\infty} A^k = I + A + A^2 + \cdots.
    $$

??? proof "Proof"
    This is the matrix version of $f(z) = \frac{1}{1-z} = \sum_{k=0}^{\infty} z^k$ (radius of convergence $R = 1$).

    By Theorem 13.3, when $\rho(A) < 1$, $\sum A^k$ converges. Let $S_N = \sum_{k=0}^N A^k$. Then
    $$
    (I - A)S_N = I - A^{N+1}.
    $$
    Since $\rho(A) < 1$, one can show $A^{N+1} \to 0$ (entrywise), so $(I - A) \lim S_N = I$, i.e., $(I-A)^{-1} = \sum_{k=0}^{\infty} A^k$. $\blacksquare$

!!! example "Example 13.2"
    Let $A = \begin{pmatrix}0.5&0.1\\0&0.3\end{pmatrix}$. Compute $(I-A)^{-1}$.

    **Solution:** The eigenvalues of $A$ are $0.5$ and $0.3$, so $\rho(A) = 0.5 < 1$ and the Neumann series converges.

    Direct computation:
    $$
    I - A = \begin{pmatrix}0.5&-0.1\\0&0.7\end{pmatrix}, \quad (I-A)^{-1} = \begin{pmatrix}2&\frac{2}{7}\\0&\frac{10}{7}\end{pmatrix}.
    $$

    Verification (approximation by the first few terms of the Neumann series):
    $$
    I + A + A^2 + A^3 + \cdots \approx \begin{pmatrix}2&0.2857\\0&1.4286\end{pmatrix} \approx \begin{pmatrix}2&\frac{2}{7}\\0&\frac{10}{7}\end{pmatrix}. \quad \checkmark
    $$

---

## 13.3 The Matrix Exponential

<div class="context-flow" markdown>

$e^A = \sum \frac{A^k}{k!}$ converges for **any** $A$ → $\det(e^A) = e^{\operatorname{tr}(A)}$ → But $e^{A+B} = e^Ae^B$ only when $AB = BA$

</div>

The matrix exponential is the most important matrix function, playing a central role in the theory of linear differential equations.

!!! definition "Definition 13.4 (Matrix Exponential)"
    For an $n \times n$ matrix $A$, the **matrix exponential** is defined as
    $$
    e^A = \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!} = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots.
    $$
    This series converges absolutely for any matrix $A$ (since the radius of convergence of $e^z$ is $\infty$).

!!! theorem "Theorem 13.5 (Basic properties of the matrix exponential)"
    Let $A, B$ be $n \times n$ matrices. Then:

    1. $e^{0} = I$;
    2. $(e^A)^{-1} = e^{-A}$, i.e., $e^A$ is always invertible;
    3. $e^{(s+t)A} = e^{sA} e^{tA}$, for any scalars $s, t$;
    4. If $AB = BA$, then $e^{A+B} = e^A e^B$;
    5. $\det(e^A) = e^{\operatorname{tr}(A)}$;
    6. $e^{PAP^{-1}} = P e^A P^{-1}$, for any invertible matrix $P$.

??? proof "Proof"
    **(1)** $e^0 = I + 0 + 0 + \cdots = I$.

    **(2)** By (4) (taking $B = -A$; clearly $A$ and $-A$ commute), $e^A e^{-A} = e^{A+(-A)} = e^0 = I$.

    **(3)** $sA$ and $tA$ commute, so by (4), $e^{sA}e^{tA} = e^{sA+tA} = e^{(s+t)A}$.

    **(4)** When $AB = BA$, the binomial theorem applies:
    $$
    (A+B)^k = \sum_{j=0}^k \binom{k}{j}A^j B^{k-j}.
    $$
    Therefore
    $$
    e^{A+B} = \sum_{k=0}^{\infty}\frac{(A+B)^k}{k!} = \sum_{k=0}^{\infty}\sum_{j=0}^k \frac{A^j B^{k-j}}{j!(k-j)!} = \left(\sum_{j=0}^{\infty}\frac{A^j}{j!}\right)\left(\sum_{l=0}^{\infty}\frac{B^l}{l!}\right) = e^A e^B.
    $$
    (Cauchy product; absolute convergence justifies rearrangement.)

    **(5)** Let the Jordan normal form of $A$ be $J$, with $A = PJP^{-1}$. Then $e^A = Pe^JP^{-1}$ and $\det(e^A) = \det(e^J)$. $e^J = \operatorname{diag}(e^{J_{k_i}(\lambda_i)})$, and $\det(e^{J_k(\lambda)}) = (e^\lambda)^k = e^{k\lambda}$. Therefore $\det(e^A) = e^{\sum k_i\lambda_i} = e^{\operatorname{tr}(A)}$.

    **(6)** $(PAP^{-1})^k = PA^kP^{-1}$; substituting into the series gives the result. $\blacksquare$

!!! note "Note"
    **Caution:** When $AB \neq BA$, in general $e^{A+B} \neq e^A e^B$. This is an important difference between the matrix exponential and the scalar exponential. For example, taking $A = \begin{pmatrix}0&1\\0&0\end{pmatrix}$ and $B = \begin{pmatrix}0&0\\1&0\end{pmatrix}$, one can verify that $e^{A+B} \neq e^A e^B$.

!!! example "Example 13.3"
    Compute $e^A$, where $A = \begin{pmatrix}0&-\theta\\\theta&0\end{pmatrix}$.

    **Solution:** Note that
    $$
    A^2 = \begin{pmatrix}-\theta^2&0\\0&-\theta^2\end{pmatrix} = -\theta^2 I, \quad A^3 = -\theta^2 A, \quad A^4 = \theta^4 I, \ldots
    $$
    In general, $A^{2k} = (-1)^k \theta^{2k} I$, $A^{2k+1} = (-1)^k \theta^{2k} A$.

    $$
    e^A = \sum_{k=0}^{\infty} \frac{A^k}{k!} = \left(\sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k}}{(2k)!}\right) I + \left(\sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k}}{(2k+1)!}\right) \frac{A}{\theta}
    $$
    $$
    = \cos\theta \cdot I + \frac{\sin\theta}{\theta} \cdot A = \begin{pmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}.
    $$

    This is precisely a rotation matrix! $A$ is skew-symmetric, and $e^A$ is orthogonal.

---

## 13.4 Computing the Matrix Exponential

<div class="context-flow" markdown>

Three approaches: **Diagonalization** $e^A = Pe^{\Lambda}P^{-1}$ · **Jordan form** $e^{J_k(\lambda)t} = e^{\lambda t}\sum \frac{t^j}{j!}N^j$ · **Cayley-Hamilton** method using eigenvalue conditions to determine coefficients

</div>

Computing the matrix exponential is a central problem in applications. This section introduces several main methods.

### 13.4.1 Diagonal Matrices

!!! theorem "Theorem 13.6 (Exponential of a diagonal matrix)"
    If $A = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$, then
    $$
    e^A = \operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n}).
    $$
    More generally, if $A = P\operatorname{diag}(\lambda_1, \ldots, \lambda_n)P^{-1}$ ($A$ is diagonalizable), then
    $$
    e^A = P\operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n})P^{-1}.
    $$

??? proof "Proof"
    $A^k = \operatorname{diag}(\lambda_1^k, \ldots, \lambda_n^k)$, so
    $$
    e^A = \sum_{k=0}^{\infty}\frac{A^k}{k!} = \operatorname{diag}\left(\sum_{k=0}^{\infty}\frac{\lambda_1^k}{k!}, \ldots, \sum_{k=0}^{\infty}\frac{\lambda_n^k}{k!}\right) = \operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n}). \qquad \blacksquare
    $$

!!! example "Example 13.4"
    Compute $e^{At}$, where $A = \begin{pmatrix}1&0\\0&-2\end{pmatrix}$.

    **Solution:** $A$ is a diagonal matrix, so
    $$
    e^{At} = \begin{pmatrix}e^t&0\\0&e^{-2t}\end{pmatrix}.
    $$

### 13.4.2 Jordan Blocks

!!! theorem "Theorem 13.7 (Exponential of a Jordan block)"
    For a Jordan block $J_k(\lambda)$,
    $$
    e^{J_k(\lambda)t} = e^{\lambda t}\begin{pmatrix}
    1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\
    0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & 1 & t \\
    0 & \cdots & 0 & 0 & 1
    \end{pmatrix}.
    $$

??? proof "Proof"
    $J_k(\lambda)t = \lambda t I + N_k t$, where $\lambda t I$ and $N_k t$ commute. Therefore
    $$
    e^{J_k(\lambda)t} = e^{\lambda t I}e^{N_k t} = e^{\lambda t}\sum_{j=0}^{k-1}\frac{(N_k t)^j}{j!} = e^{\lambda t}\sum_{j=0}^{k-1}\frac{t^j}{j!}N_k^j.
    $$
    Since $N_k^j$ has ones on the $j$-th superdiagonal and zeros elsewhere, substituting gives the result. $\blacksquare$

!!! example "Example 13.5"
    Compute $e^{At}$, where $A = \begin{pmatrix}3&1&0\\0&3&1\\0&0&3\end{pmatrix} = J_3(3)$.

    **Solution:** By Theorem 13.7,
    $$
    e^{At} = e^{3t}\begin{pmatrix}1&t&\frac{t^2}{2}\\0&1&t\\0&0&1\end{pmatrix}.
    $$

### 13.4.3 The Cayley-Hamilton Method

!!! definition "Definition 13.5 (Cayley-Hamilton method for computing matrix functions)"
    By the Cayley-Hamilton theorem, $A^n$ can be expressed as a linear combination of $I, A, \ldots, A^{n-1}$. Therefore any matrix function $f(A)$ can be written as
    $$
    f(A) = \alpha_0 I + \alpha_1 A + \cdots + \alpha_{n-1} A^{n-1},
    $$
    where the coefficients $\alpha_0, \ldots, \alpha_{n-1}$ are determined by the following conditions: for each eigenvalue $\lambda_i$ of $A$ (with algebraic multiplicity $m_i$),
    $$
    f^{(j)}(\lambda_i) = \alpha_0^{(j)} + \alpha_1 \cdot j! + \cdots \quad (j = 0, 1, \ldots, m_i - 1),
    $$
    i.e., the function values and derivatives of the $\alpha$-polynomial at $\lambda_i$ match those of $f$.

!!! example "Example 13.6"
    Compute $e^{At}$ using the Cayley-Hamilton method, where $A = \begin{pmatrix}2&1\\0&2\end{pmatrix}$.

    **Solution:** Eigenvalue $\lambda = 2$ (algebraic multiplicity 2). $n = 2$, so let
    $$
    e^{At} = \alpha_0(t) I + \alpha_1(t) A.
    $$

    From the condition $f(\lambda) = e^{\lambda t}$ matching at $\lambda = 2$ up to first derivative:

    - $f(2) = e^{2t}$: $\alpha_0 + 2\alpha_1 = e^{2t}$;
    - $f'(2) = te^{2t}$: $\alpha_1 = te^{2t}$.

    Solving: $\alpha_1 = te^{2t}$, $\alpha_0 = e^{2t} - 2te^{2t}$.

    $$
    e^{At} = (e^{2t} - 2te^{2t})I + te^{2t}A = e^{2t}\begin{pmatrix}1-2t&0\\0&1-2t\end{pmatrix} + te^{2t}\begin{pmatrix}2&1\\0&2\end{pmatrix}
    $$
    $$
    = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}.
    $$

---

## 13.5 The Matrix Exponential and Differential Equations

<div class="context-flow" markdown>

$\mathbf{x}' = A\mathbf{x}$ → $\mathbf{x}(t) = e^{At}\mathbf{x}_0$ (unique solution) → Nonhomogeneous case via **variation of parameters** → Jordan form determines asymptotic behavior ($e^{\lambda t}$ times polynomial)

</div>

The most important application of the matrix exponential is solving linear constant-coefficient systems of differential equations.

<div class="context-flow" markdown>

**Insight**: The cleverness of the uniqueness proof — let $\mathbf{z}(t) = e^{-At}\mathbf{y}(t)$, use $\mathbf{z}'=0$ to conclude $\mathbf{y} = e^{At}\mathbf{x}_0$

</div>

!!! theorem "Theorem 13.8 (Solution of homogeneous linear ODE systems)"
    The ODE system
    $$
    \mathbf{x}'(t) = A\mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$
    has the unique solution
    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0.
    $$

??? proof "Proof"
    **Existence:** Let $\mathbf{x}(t) = e^{At}\mathbf{x}_0$. Then
    $$
    \mathbf{x}'(t) = \frac{d}{dt}e^{At}\mathbf{x}_0 = Ae^{At}\mathbf{x}_0 = A\mathbf{x}(t),
    $$
    where $\frac{d}{dt}e^{At} = \sum_{k=1}^{\infty}\frac{kA^k t^{k-1}}{k!} = A\sum_{k=1}^{\infty}\frac{A^{k-1}t^{k-1}}{(k-1)!} = Ae^{At}$.

    Also $\mathbf{x}(0) = e^{0}\mathbf{x}_0 = \mathbf{x}_0$.

    **Uniqueness:** Let $\mathbf{y}(t)$ be another solution. Define $\mathbf{z}(t) = e^{-At}\mathbf{y}(t)$. Then
    $$
    \mathbf{z}'(t) = -Ae^{-At}\mathbf{y}(t) + e^{-At}\mathbf{y}'(t) = -Ae^{-At}\mathbf{y}(t) + e^{-At}A\mathbf{y}(t) = 0.
    $$
    Therefore $\mathbf{z}(t) = \mathbf{z}(0) = \mathbf{x}_0$, i.e., $\mathbf{y}(t) = e^{At}\mathbf{x}_0$. $\blacksquare$

!!! theorem "Theorem 13.9 (Nonhomogeneous linear ODE systems)"
    The ODE system
    $$
    \mathbf{x}'(t) = A\mathbf{x}(t) + \mathbf{f}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$
    has the solution (variation of parameters)
    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds.
    $$

??? proof "Proof"
    Let $\mathbf{x}(t) = e^{At}\mathbf{c}(t)$ (variation of parameters). Substituting into the equation gives
    $$
    Ae^{At}\mathbf{c}(t) + e^{At}\mathbf{c}'(t) = Ae^{At}\mathbf{c}(t) + \mathbf{f}(t).
    $$
    Simplifying: $e^{At}\mathbf{c}'(t) = \mathbf{f}(t)$, i.e., $\mathbf{c}'(t) = e^{-At}\mathbf{f}(t)$. Integrating:
    $$
    \mathbf{c}(t) = \mathbf{x}_0 + \int_0^t e^{-As}\mathbf{f}(s)\,ds.
    $$
    Therefore $\mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds$. $\blacksquare$

!!! example "Example 13.7"
    Solve the ODE system
    $$
    \begin{cases} x_1' = 3x_1 + x_2, \\ x_2' = -x_1 + x_2, \end{cases} \quad \mathbf{x}(0) = \begin{pmatrix}1\\0\end{pmatrix}.
    $$

    **Solution:** $A = \begin{pmatrix}3&1\\-1&1\end{pmatrix}$. Eigenvalues: $\lambda^2 - 4\lambda + 4 = (\lambda-2)^2 = 0$, so $\lambda = 2$ (repeated root).

    $A - 2I = \begin{pmatrix}1&1\\-1&-1\end{pmatrix}$, $\operatorname{rank} = 1$, geometric multiplicity = 1.

    Jordan form $J = J_2(2)$. Eigenvector $\mathbf{v}_1 = \begin{pmatrix}1\\-1\end{pmatrix}$, generalized eigenvector $(A-2I)\mathbf{v}_2 = \mathbf{v}_1$:
    $$
    \begin{pmatrix}1&1\\-1&-1\end{pmatrix}\mathbf{v}_2 = \begin{pmatrix}1\\-1\end{pmatrix}, \quad \Rightarrow \quad \mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix}.
    $$

    $P = \begin{pmatrix}1&1\\-1&0\end{pmatrix}$, $P^{-1} = \begin{pmatrix}0&-1\\1&1\end{pmatrix}$.

    $$
    e^{At} = Pe^{Jt}P^{-1} = \begin{pmatrix}1&1\\-1&0\end{pmatrix}e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}\begin{pmatrix}0&-1\\1&1\end{pmatrix}
    $$
    $$
    = e^{2t}\begin{pmatrix}1&1\\-1&0\end{pmatrix}\begin{pmatrix}t&-1+t\\1&1\end{pmatrix} = e^{2t}\begin{pmatrix}1+t&t\\-t&1-t\end{pmatrix}.
    $$

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}(0) = e^{2t}\begin{pmatrix}1+t\\-t\end{pmatrix}.
    $$

---

## 13.6 Matrix Logarithm

<div class="context-flow" markdown>

Inverse problem of $e^X = A$ → Every invertible matrix has a logarithm → For Jordan blocks $\log J_k(\lambda) = (\log\lambda)I + \sum \frac{(-1)^{j+1}}{j}(\lambda^{-1}N)^j$ (finite sum)

</div>

The matrix logarithm is the inverse operation of the matrix exponential.

!!! definition "Definition 13.6 (Matrix Logarithm)"
    Let $A$ be an $n \times n$ invertible matrix. If there exists a matrix $X$ such that $e^X = A$, then $X$ is called a **matrix logarithm** of $A$, written $X = \log A$ or $X = \ln A$.

!!! theorem "Theorem 13.10 (Existence of the matrix logarithm)"
    Let $A$ be an $n \times n$ invertible complex matrix. Then a matrix logarithm of $A$ exists, i.e., there exists a matrix $X$ such that $e^X = A$.

    More precisely, if $A$ has no negative real eigenvalues, then there exists a unique matrix logarithm $X$ such that all eigenvalues of $X$ have imaginary parts in $(-\pi, \pi)$. This $X$ is called the **principal logarithm** of $A$.

??? proof "Proof"
    **Existence (constructive proof):**

    Let $A = PJP^{-1}$, $J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$. It suffices to define the logarithm for each Jordan block.

    For $J_k(\lambda)$ ($\lambda \neq 0$), write
    $$
    J_k(\lambda) = \lambda(I + \lambda^{-1}N_k) = \lambda(I + M),
    $$
    where $M = \lambda^{-1}N_k$ is nilpotent. Take
    $$
    \log J_k(\lambda) = (\log\lambda) I + \log(I + M) = (\log\lambda)I + \sum_{j=1}^{k-1}\frac{(-1)^{j+1}}{j}M^j.
    $$
    The series is finite (since $M$ is nilpotent), and $e^{\log J_k(\lambda)} = J_k(\lambda)$.

    Let $\log A = P \operatorname{diag}(\log J_{k_1}(\lambda_1), \ldots) P^{-1}$. $\blacksquare$

!!! example "Example 13.8"
    Find $\log A$, where $A = \begin{pmatrix}1&1\\0&1\end{pmatrix} = J_2(1)$.

    **Solution:** $A = I + N$, where $N = \begin{pmatrix}0&1\\0&0\end{pmatrix}$.

    $$
    \log A = \log(I + N) = N - \frac{N^2}{2} + \frac{N^3}{3} - \cdots = N = \begin{pmatrix}0&1\\0&0\end{pmatrix}.
    $$
    (Since $N^2 = 0$, the series has only one term.)

    Verification: $e^N = I + N + \frac{N^2}{2!} + \cdots = I + N = A$. $\checkmark$

---

## 13.7 Matrix Square Root

<div class="context-flow" markdown>

$X^2 = A$ → Positive definite matrices have a unique **positive definite square root** $A^{1/2} = Q\Lambda^{1/2}Q^T$ → Can also be defined via $e^{\frac{1}{2}\log A}$ → Appears in Ch10 polar decomposition $P = (A^HA)^{1/2}$

</div>

!!! definition "Definition 13.7 (Matrix Square Root)"
    Let $A$ be an $n \times n$ matrix. If there exists a matrix $X$ such that $X^2 = A$, then $X$ is called a **matrix square root** of $A$, written $X = A^{1/2}$.

!!! theorem "Theorem 13.11 (Unique positive definite square root of a positive definite matrix)"
    Let $A$ be an $n \times n$ real symmetric positive definite matrix. Then there exists a unique real symmetric positive definite matrix $B$ such that $B^2 = A$. $B$ is called the **positive definite square root** of $A$.

??? proof "Proof"
    **Existence:** $A$ is real symmetric positive definite, so by the spectral theorem, $A = Q\Lambda Q^T$, where $Q$ is orthogonal, $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$, and $\lambda_i > 0$.

    Let $B = Q\Lambda^{1/2}Q^T$, where $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$. Then
    $$
    B^2 = Q\Lambda^{1/2}Q^T Q\Lambda^{1/2}Q^T = Q\Lambda Q^T = A.
    $$
    $B$ is symmetric: $B^T = (Q\Lambda^{1/2}Q^T)^T = Q\Lambda^{1/2}Q^T = B$. The eigenvalues of $B$ are $\sqrt{\lambda_i} > 0$, so $B$ is positive definite.

    **Uniqueness:** Suppose $C$ is also symmetric positive definite with $C^2 = A$. Since $C$ is symmetric positive definite, let $C = R\Gamma R^T$, $\Gamma = \operatorname{diag}(\gamma_1, \ldots, \gamma_n)$, $\gamma_i > 0$. Then $A = C^2 = R\Gamma^2 R^T$.

    By the uniqueness of the spectral decomposition of $A$ (once eigenvalues are determined, eigenspaces are determined), $\Gamma^2 = \Lambda$ (after reordering eigenvalues). Since $\gamma_i > 0$, $\gamma_i = \sqrt{\lambda_i}$, so $C = B$. $\blacksquare$

!!! definition "Definition 13.8 (Square root of a positive semidefinite matrix)"
    For a real symmetric positive semidefinite matrix $A$ (eigenvalues $\lambda_i \ge 0$), the positive definite square root construction generalizes to: $A^{1/2} = Q\operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})Q^T$. In this case $A^{1/2}$ is positive semidefinite (the unique positive semidefinite square root).

!!! theorem "Theorem 13.12 (Existence of square root for invertible matrices)"
    Let $A$ be an $n \times n$ invertible complex matrix with no negative real eigenvalues. Then $A$ has a unique square root $A^{1/2}$ such that all eigenvalues of $A^{1/2}$ have positive real part.

??? proof "Proof"
    Using the principal logarithm: let $A^{1/2} = e^{\frac{1}{2}\log A}$, where $\log A$ is the principal logarithm. Then
    $$
    (A^{1/2})^2 = e^{\frac{1}{2}\log A} e^{\frac{1}{2}\log A} = e^{\log A} = A.
    $$
    The eigenvalues of $A^{1/2}$ are $e^{\frac{1}{2}\log\lambda_i}$, where $\log\lambda_i$ has imaginary part in $(-\pi, \pi)$, so $\frac{1}{2}\log\lambda_i$ has imaginary part in $(-\frac{\pi}{2}, \frac{\pi}{2})$, and $e^{\frac{1}{2}\log\lambda_i}$ has positive real part. Uniqueness follows from this condition. $\blacksquare$

!!! example "Example 13.9"
    Find the positive definite square root of $A = \begin{pmatrix}2&1\\1&2\end{pmatrix}$.

    **Solution:** Eigenvalues $\lambda_1 = 3$, $\lambda_2 = 1$. Orthogonal eigenvectors: $\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$, $\mathbf{q}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\-1\end{pmatrix}$.

    $$
    A^{1/2} = Q\begin{pmatrix}\sqrt{3}&0\\0&1\end{pmatrix}Q^T = \frac{1}{2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}\begin{pmatrix}\sqrt{3}&0\\0&1\end{pmatrix}\begin{pmatrix}1&1\\1&-1\end{pmatrix}
    $$
    $$
    = \frac{1}{2}\begin{pmatrix}\sqrt{3}+1&\sqrt{3}-1\\\sqrt{3}-1&\sqrt{3}+1\end{pmatrix}.
    $$

    Verification: $(A^{1/2})^2 = \frac{1}{4}\begin{pmatrix}\sqrt{3}+1&\sqrt{3}-1\\\sqrt{3}-1&\sqrt{3}+1\end{pmatrix}^2$.

    Diagonal entry: $\frac{1}{4}[(\sqrt{3}+1)^2 + (\sqrt{3}-1)^2] = \frac{1}{4}(4+2\sqrt{3}+4-2\sqrt{3}) = \frac{8}{4} = 2$.

    Off-diagonal entry: $\frac{1}{4}[(\sqrt{3}+1)(\sqrt{3}-1) + (\sqrt{3}-1)(\sqrt{3}+1)] = \frac{1}{4}(2 \times 2) = 1$.

    Therefore $(A^{1/2})^2 = \begin{pmatrix}2&1\\1&2\end{pmatrix} = A$. $\checkmark$

---

## 13.8 General Matrix Functions

<div class="context-flow" markdown>

Unified framework: the $(p,q)$ entry of $f(J_k(\lambda))$ is $\frac{f^{(q-p)}(\lambda)}{(q-p)!}$ → Cauchy integral definition $f(A) = \frac{1}{2\pi i}\oint f(z)(zI-A)^{-1}dz$ is equivalent → Spectral mapping $\sigma(f(A)) = f(\sigma(A))$

</div>

The preceding sections discussed specific matrix functions (exponential, logarithm, square root). This section gives the general definition framework for matrix functions.

### 13.8.1 Definition via Jordan Normal Form

<div class="context-flow" markdown>

**Insight**: $f(J_k(\lambda))$ involves $f, f', f'', \ldots, f^{(k-1)}$ — the size of the Jordan block determines how many orders of **differentiability** $f$ requires

</div>

!!! definition "Definition 13.9 (General matrix function — Jordan form definition)"
    Let $f$ be a function defined on the spectrum of $A$, with $A = PJP^{-1}$ and $J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$. Define
    $$
    f(A) = P f(J) P^{-1} = P \operatorname{diag}(f(J_{k_1}(\lambda_1)), \ldots, f(J_{k_s}(\lambda_s))) P^{-1},
    $$
    where for each Jordan block
    $$
    f(J_k(\lambda)) = \begin{pmatrix}
    f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\
    0 & f(\lambda) & f'(\lambda) & \cdots & \frac{f^{(k-2)}(\lambda)}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & f(\lambda) & f'(\lambda) \\
    0 & \cdots & 0 & 0 & f(\lambda)
    \end{pmatrix}.
    $$
    This requires $f$ to be at least $k_i - 1$ times differentiable at each eigenvalue $\lambda_i$.

!!! theorem "Theorem 13.13 (Well-definedness of the Jordan form definition)"
    The $f(A)$ in Definition 13.9 does not depend on the choice of Jordan decomposition (i.e., does not depend on the choice of $P$), so $f(A)$ is well-defined.

??? proof "Proof"
    Suppose $A = P_1 J P_1^{-1} = P_2 J P_2^{-1}$ (same Jordan form $J$, but different transition matrices). Then $P_2^{-1}P_1$ commutes with $J$, i.e., $P_2^{-1}P_1 J = J P_2^{-1}P_1$.

    Since $f(J)$ is the limit of polynomials in $J$ (when $f$ is analytic), $P_2^{-1}P_1$ also commutes with $f(J)$. Therefore
    $$
    P_1 f(J) P_1^{-1} = P_2 (P_2^{-1}P_1) f(J) (P_2^{-1}P_1)^{-1} P_2^{-1} = P_2 f(J) P_2^{-1}. \qquad \blacksquare
    $$

!!! example "Example 13.10"
    Compute $\sin(A)$, where $A = \begin{pmatrix}0&\pi\\0&0\end{pmatrix}$.

    **Solution:** $A = J_2(0)$, $\lambda = 0$, $k = 2$.

    $f(\lambda) = \sin(\lambda)$. $f(0) = 0$, $f'(0) = \cos(0) = 1$.

    $$
    \sin(A) = f(J_2(0)) = \begin{pmatrix}f(0)&f'(0)\\0&f(0)\end{pmatrix} = \begin{pmatrix}0&1\\0&0\end{pmatrix}.
    $$

    Note that this is not the entrywise $\sin$! $\sin\begin{pmatrix}0&\pi\\0&0\end{pmatrix} \neq \begin{pmatrix}\sin 0&\sin\pi\\\sin 0&\sin 0\end{pmatrix}$.

### 13.8.2 Cauchy Integral Definition

!!! definition "Definition 13.10 (Matrix function — Cauchy integral definition)"
    Let $f$ be analytic on an open set $\Omega$ containing all eigenvalues of $A$, and let $\Gamma$ be a simple closed curve in $\Omega$ enclosing all eigenvalues. Define
    $$
    f(A) = \frac{1}{2\pi i} \oint_\Gamma f(z)(zI - A)^{-1}\,dz.
    $$

!!! theorem "Theorem 13.14 (Equivalence of the Cauchy integral and Jordan form definitions)"
    When $f$ is analytic on some open neighborhood of the spectrum of $A$, Definitions 13.9 and 13.10 give the same $f(A)$.

??? proof "Proof"
    Let $A = PJP^{-1}$. Then $(zI - A)^{-1} = P(zI - J)^{-1}P^{-1}$, so
    $$
    \frac{1}{2\pi i}\oint_\Gamma f(z)(zI-A)^{-1}dz = P\left(\frac{1}{2\pi i}\oint_\Gamma f(z)(zI-J)^{-1}dz\right)P^{-1}.
    $$

    For a Jordan block $J_k(\lambda)$, the $(p,q)$ entry ($q \ge p$) of $(zI - J_k(\lambda))^{-1}$ is $\frac{1}{(z-\lambda)^{q-p+1}}$. By the Cauchy integral formula:
    $$
    \frac{1}{2\pi i}\oint \frac{f(z)}{(z-\lambda)^{q-p+1}}dz = \frac{f^{(q-p)}(\lambda)}{(q-p)!}.
    $$
    This is precisely the $(p,q)$ entry of $f(J_k(\lambda))$ from Definition 13.9. $\blacksquare$

!!! theorem "Theorem 13.15 (Spectral Mapping Theorem for matrix functions)"
    Let $f$ be an analytic function on the spectrum of matrix $A$. Then
    $$
    \sigma(f(A)) = f(\sigma(A)) = \{f(\lambda) : \lambda \in \sigma(A)\},
    $$
    i.e., the eigenvalues of $f(A)$ are exactly the images of the eigenvalues of $A$ under $f$.

??? proof "Proof"
    Let $A = PJP^{-1}$, so $f(A) = Pf(J)P^{-1}$. The diagonal entries of $f(J)$ are $f(\lambda_i)$ (the diagonal entries of each Jordan block), and these are precisely the eigenvalues of $f(A)$ (similarity transformations preserve the spectrum). $\blacksquare$

!!! example "Example 13.11"
    Suppose the eigenvalues of $A$ are $1, 2, 3$. Find the eigenvalues of $e^A$ and of $\cos(A)$.

    **Solution:** By the spectral mapping theorem:

    - The eigenvalues of $e^A$ are $e^1, e^2, e^3$, i.e., $e, e^2, e^3$.
    - The eigenvalues of $\cos(A)$ are $\cos 1, \cos 2, \cos 3$.

!!! example "Example 13.12"
    Let $A = \begin{pmatrix}1&0&0\\0&2&1\\0&0&2\end{pmatrix}$. Compute $\sqrt{A}$ (the principal square root).

    **Solution:** $A = J_1(1) \oplus J_2(2)$. $f(\lambda) = \sqrt{\lambda}$, $f'(\lambda) = \frac{1}{2\sqrt{\lambda}}$.

    $$
    f(J_1(1)) = (1) = (\sqrt{1}) = (1).
    $$
    $$
    f(J_2(2)) = \begin{pmatrix} f(2) & f'(2) \\ 0 & f(2) \end{pmatrix} = \begin{pmatrix} \sqrt{2} & \frac{1}{2\sqrt{2}} \\ 0 & \sqrt{2} \end{pmatrix} = \begin{pmatrix} \sqrt{2} & \frac{\sqrt{2}}{4} \\ 0 & \sqrt{2} \end{pmatrix}.
    $$

    Therefore
    $$
    \sqrt{A} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \sqrt{2} & \frac{\sqrt{2}}{4} \\ 0 & 0 & \sqrt{2} \end{pmatrix}.
    $$

    Verification: the $(2,3)$ entry of $(\sqrt{A})^2$ is $\sqrt{2} \cdot \frac{\sqrt{2}}{4} + \frac{\sqrt{2}}{4} \cdot \sqrt{2} = \frac{2}{4} + \frac{2}{4} = 1 = A_{23}$. $\checkmark$

---

## Chapter Summary

This chapter systematically introduced the theory of matrix functions, including:

1. **Matrix polynomials** $p(A)$ are the most basic matrix functions, with the Cayley-Hamilton theorem as the core result;
2. **Matrix power series** convergence is controlled by the spectral radius $\rho(A)$;
3. The **matrix exponential** $e^A$ is defined for any matrix, satisfies $\det(e^A) = e^{\operatorname{tr}(A)}$, but note that $e^{A+B} = e^A e^B$ holds only when $AB = BA$;
4. **Matrix exponential and differential equations**: the solution of $\mathbf{x}'(t) = A\mathbf{x}(t)$ is $\mathbf{x}(t) = e^{At}\mathbf{x}_0$;
5. The **matrix logarithm** exists for invertible matrices; positive definite matrices have a unique positive semidefinite logarithm;
6. **Matrix square root**: positive definite matrices have a unique positive definite square root $A^{1/2} = Q\Lambda^{1/2}Q^T$;
7. **General matrix functions** can be defined via the Jordan normal form or the Cauchy integral; the spectral mapping theorem $\sigma(f(A)) = f(\sigma(A))$ is an important property.

Matrix function theory deeply fuses calculus with linear algebra and is an important toolbox of modern applied mathematics.
