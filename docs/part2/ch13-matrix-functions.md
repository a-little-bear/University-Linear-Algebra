# 第 13 章 矩阵函数

<div class="context-flow" markdown>

**前置**：Ch12 Jordan 标准形 · Ch10 谱分解 · **本章脉络**：$p(A)$（多项式） → Cayley-Hamilton → 幂级数/收敛 → $e^A$（矩阵指数） → $\log A$, $A^{1/2}$ → 一般 $f(A)$（Jordan/Cauchy 积分）
本质：$f(\lambda) \to f(A)$ —— Jordan 形让标量函数**逐块作用于矩阵**，导数信息填充超对角线

</div>

在前面的章节中，我们已经熟悉了矩阵的加法、乘法和求逆等基本运算。本章将把函数的概念从标量推广到矩阵：给定一个函数 $f$（如指数函数、对数函数、平方根等），我们希望定义 $f(A)$ 的含义。矩阵函数（matrix function）在微分方程、控制理论、量子力学等领域有着重要应用。本章从矩阵多项式出发，经由幂级数和 Jordan 标准形，逐步建立矩阵函数的完整理论。

---

## 13.1 矩阵多项式

<div class="context-flow" markdown>

$p(A) = a_k A^k + \cdots + a_0 I$：最基础的矩阵函数 → 保相似性 $p(PBP^{-1}) = Pp(B)P^{-1}$ → **Cayley-Hamilton**：$p_A(A) = 0$

</div>

矩阵多项式是定义矩阵函数的最基础方式，也是理解更一般矩阵函数的起点。

!!! definition "定义 13.1 (矩阵多项式 Matrix Polynomial)"
    设 $p(\lambda) = a_k \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_1\lambda + a_0$ 为一个标量多项式，$A$ 为 $n \times n$ 矩阵。定义**矩阵多项式**
    $$
    p(A) = a_k A^k + a_{k-1}A^{k-1} + \cdots + a_1 A + a_0 I,
    $$
    其中 $A^0 = I$（单位矩阵）。

!!! theorem "定理 13.1 (矩阵多项式的基本性质)"
    设 $p, q$ 为标量多项式，$A$ 为 $n \times n$ 矩阵，则：

    1. $(p + q)(A) = p(A) + q(A)$；
    2. $(pq)(A) = p(A)q(A)$；
    3. $p(A)$ 与 $q(A)$ 可交换：$p(A)q(A) = q(A)p(A)$；
    4. 若 $A = PBP^{-1}$，则 $p(A) = Pp(B)P^{-1}$；
    5. $p(A)$ 的特征值为 $p(\lambda_1), \ldots, p(\lambda_n)$，其中 $\lambda_i$ 为 $A$ 的特征值。

??? proof "证明"
    **(1)** 和 **(2)** 由矩阵乘法的分配律和结合律直接得出。

    **(3)** 由 (2)，$p(A)q(A) = (pq)(A) = (qp)(A) = q(A)p(A)$（利用标量多项式乘法的交换性）。

    **(4)** $A^k = PB^kP^{-1}$（数学归纳法易证），因此
    $$
    p(A) = \sum a_i A^i = \sum a_i PB^iP^{-1} = P\left(\sum a_i B^i\right)P^{-1} = Pp(B)P^{-1}.
    $$

    **(5)** 设 $A\mathbf{v} = \lambda\mathbf{v}$，则 $A^k\mathbf{v} = \lambda^k\mathbf{v}$，因此 $p(A)\mathbf{v} = p(\lambda)\mathbf{v}$。即 $p(\lambda)$ 为 $p(A)$ 的特征值。$\blacksquare$

<div class="context-flow" markdown>

**洞察**：Cayley-Hamilton 意味着 $A^n$ 可表示为 $I, A, \ldots, A^{n-1}$ 的线性组合——矩阵代数 $\mathbb{F}[A]$ 的维数 $\le n$

</div>

!!! theorem "定理 13.2 (Cayley-Hamilton 定理)"
    设 $A$ 为 $n \times n$ 矩阵，$p_A(\lambda) = \det(\lambda I - A)$ 为其特征多项式，则
    $$
    p_A(A) = 0.
    $$
    即每个矩阵都满足自己的特征方程。

??? proof "证明"
    设 $p_A(\lambda) = \lambda^n + c_{n-1}\lambda^{n-1} + \cdots + c_0$。

    **方法一（伴随矩阵法）：** 设 $B(\lambda) = \operatorname{adj}(\lambda I - A)$ 为 $\lambda I - A$ 的伴随矩阵（adjugate），则
    $$
    (\lambda I - A)B(\lambda) = \det(\lambda I - A) \cdot I = p_A(\lambda) I.
    $$
    $B(\lambda)$ 的每个元素是 $\lambda$ 的至多 $(n-1)$ 次多项式，因此可以写为
    $$
    B(\lambda) = B_{n-1}\lambda^{n-1} + B_{n-2}\lambda^{n-2} + \cdots + B_0,
    $$
    其中 $B_i$ 为常数矩阵。将 $(\lambda I - A)B(\lambda) = p_A(\lambda) I$ 展开并比较 $\lambda$ 的各次幂系数，可得到一组矩阵等式。将第 $k$ 个等式左乘 $A^k$ 后全部相加，利用抵消可得 $p_A(A) = 0$。$\blacksquare$

!!! example "例 13.1"
    验证 Cayley-Hamilton 定理对 $A = \begin{pmatrix}1&2\\3&4\end{pmatrix}$ 成立。

    **解：** 特征多项式：
    $$
    p_A(\lambda) = \lambda^2 - 5\lambda - 2.
    $$

    计算 $p_A(A) = A^2 - 5A - 2I$：
    $$
    A^2 = \begin{pmatrix}7&10\\15&22\end{pmatrix}, \quad 5A = \begin{pmatrix}5&10\\15&20\end{pmatrix}, \quad 2I = \begin{pmatrix}2&0\\0&2\end{pmatrix}.
    $$
    $$
    p_A(A) = \begin{pmatrix}7&10\\15&22\end{pmatrix} - \begin{pmatrix}5&10\\15&20\end{pmatrix} - \begin{pmatrix}2&0\\0&2\end{pmatrix} = \begin{pmatrix}0&0\\0&0\end{pmatrix}. \quad \checkmark
    $$

---

## 13.2 矩阵幂级数

<div class="context-flow" markdown>

$f(A) = \sum c_k A^k$ 收敛 ↔ **谱半径** $\rho(A) < R$（收敛半径） → Neumann 级数 $(I-A)^{-1} = \sum A^k$（$\rho(A)<1$）

</div>

将多项式推广到幂级数，需要引入矩阵级数的收敛性概念。

!!! definition "定义 13.2 (矩阵级数的收敛 Convergence of Matrix Series)"
    设 $\{A_k\}$ 为 $m \times n$ 矩阵序列。称 $\sum_{k=0}^{\infty} A_k$ **收敛**，若每个元素位置的标量级数都收敛，即对所有 $1 \le i \le m$、$1 \le j \le n$，$\sum_{k=0}^{\infty} [A_k]_{ij}$ 收敛。此时定义
    $$
    \sum_{k=0}^{\infty} A_k = \left[\sum_{k=0}^{\infty} [A_k]_{ij}\right]_{m \times n}.
    $$

!!! definition "定义 13.3 (谱半径 Spectral Radius)"
    矩阵 $A$ 的**谱半径**（spectral radius）定义为
    $$
    \rho(A) = \max\{|\lambda| : \lambda \text{ 为 } A \text{ 的特征值}\}.
    $$

!!! theorem "定理 13.3 (矩阵幂级数的收敛判定)"
    设 $f(z) = \sum_{k=0}^{\infty} c_k z^k$ 为收敛半径为 $R$ 的幂级数，$A$ 为 $n \times n$ 矩阵。则矩阵幂级数
    $$
    f(A) = \sum_{k=0}^{\infty} c_k A^k
    $$
    收敛的充分条件是 $\rho(A) < R$。

??? proof "证明"
    设 $A = PJP^{-1}$，其中 $J$ 为 Jordan 标准形。则 $A^k = PJ^kP^{-1}$，故
    $$
    f(A) = P\left(\sum_{k=0}^{\infty} c_k J^k\right)P^{-1} = Pf(J)P^{-1}.
    $$
    由于 $J$ 为分块对角矩阵，$f(J) = \operatorname{diag}(f(J_{n_1}(\lambda_1)), \ldots)$。

    对 Jordan 块 $J_m(\lambda)$，$f(J_m(\lambda))$ 的 $(p,q)$ 元素（$q \ge p$）为
    $$
    \frac{f^{(q-p)}(\lambda)}{(q-p)!} = \sum_{k=q-p}^{\infty} c_k \binom{k}{q-p} \lambda^{k-q+p}.
    $$
    当 $|\lambda| < R$ 时，$f$ 在 $\lambda$ 处的各阶导数都收敛，因此上式收敛。$\rho(A) < R$ 保证所有特征值 $\lambda$ 满足 $|\lambda| < R$。$\blacksquare$

!!! theorem "定理 13.4 (Neumann 级数)"
    设 $A$ 为 $n \times n$ 矩阵。若 $\rho(A) < 1$，则 $I - A$ 可逆且
    $$
    (I - A)^{-1} = \sum_{k=0}^{\infty} A^k = I + A + A^2 + \cdots.
    $$

??? proof "证明"
    这是 $f(z) = \frac{1}{1-z} = \sum_{k=0}^{\infty} z^k$（收敛半径 $R = 1$）的矩阵版本。

    由定理 13.3，当 $\rho(A) < 1$ 时，$\sum A^k$ 收敛。设 $S_N = \sum_{k=0}^N A^k$，则
    $$
    (I - A)S_N = I - A^{N+1}.
    $$
    由 $\rho(A) < 1$，可证 $A^{N+1} \to 0$（逐元素），因此 $(I - A) \lim S_N = I$，即 $(I-A)^{-1} = \sum_{k=0}^{\infty} A^k$。$\blacksquare$

!!! example "例 13.2"
    设 $A = \begin{pmatrix}0.5&0.1\\0&0.3\end{pmatrix}$，计算 $(I-A)^{-1}$。

    **解：** $A$ 的特征值为 $0.5$ 和 $0.3$，$\rho(A) = 0.5 < 1$，因此 Neumann 级数收敛。

    直接计算：
    $$
    I - A = \begin{pmatrix}0.5&-0.1\\0&0.7\end{pmatrix}, \quad (I-A)^{-1} = \begin{pmatrix}2&\frac{2}{7}\\0&\frac{10}{7}\end{pmatrix}.
    $$

    验证（用 Neumann 级数前几项近似）：
    $$
    I + A + A^2 + A^3 + \cdots \approx \begin{pmatrix}2&0.2857\\0&1.4286\end{pmatrix} \approx \begin{pmatrix}2&\frac{2}{7}\\0&\frac{10}{7}\end{pmatrix}. \quad \checkmark
    $$

---

## 13.3 矩阵指数

<div class="context-flow" markdown>

$e^A = \sum \frac{A^k}{k!}$ 对**任意** $A$ 收敛 → $\det(e^A) = e^{\operatorname{tr}(A)}$ → 但 $e^{A+B} = e^Ae^B$ 仅当 $AB = BA$

</div>

矩阵指数（matrix exponential）是最重要的矩阵函数，在线性微分方程理论中起核心作用。

!!! definition "定义 13.4 (矩阵指数 Matrix Exponential)"
    对 $n \times n$ 矩阵 $A$，**矩阵指数**定义为
    $$
    e^A = \exp(A) = \sum_{k=0}^{\infty} \frac{A^k}{k!} = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots.
    $$
    该级数对任意矩阵 $A$ 都绝对收敛（因为 $e^z$ 的收敛半径为 $\infty$）。

!!! theorem "定理 13.5 (矩阵指数的基本性质)"
    设 $A, B$ 为 $n \times n$ 矩阵，则：

    1. $e^{0} = I$；
    2. $(e^A)^{-1} = e^{-A}$，即 $e^A$ 总是可逆的；
    3. $e^{(s+t)A} = e^{sA} e^{tA}$，对任意标量 $s, t$；
    4. 若 $AB = BA$，则 $e^{A+B} = e^A e^B$；
    5. $\det(e^A) = e^{\operatorname{tr}(A)}$；
    6. $e^{PAP^{-1}} = P e^A P^{-1}$，对任意可逆矩阵 $P$。

??? proof "证明"
    **(1)** $e^0 = I + 0 + 0 + \cdots = I$。

    **(2)** 由 (4)（取 $B = -A$，显然 $A$ 与 $-A$ 可交换），$e^A e^{-A} = e^{A+(-A)} = e^0 = I$。

    **(3)** $sA$ 与 $tA$ 可交换，由 (4) 得 $e^{sA}e^{tA} = e^{sA+tA} = e^{(s+t)A}$。

    **(4)** 当 $AB = BA$ 时，可用二项式定理：
    $$
    (A+B)^k = \sum_{j=0}^k \binom{k}{j}A^j B^{k-j}.
    $$
    因此
    $$
    e^{A+B} = \sum_{k=0}^{\infty}\frac{(A+B)^k}{k!} = \sum_{k=0}^{\infty}\sum_{j=0}^k \frac{A^j B^{k-j}}{j!(k-j)!} = \left(\sum_{j=0}^{\infty}\frac{A^j}{j!}\right)\left(\sum_{l=0}^{\infty}\frac{B^l}{l!}\right) = e^A e^B.
    $$
    （Cauchy 乘积，绝对收敛保证重排合法。）

    **(5)** 设 $A$ 的 Jordan 标准形为 $J$，$A = PJP^{-1}$。则 $e^A = Pe^JP^{-1}$，$\det(e^A) = \det(e^J)$。$e^J = \operatorname{diag}(e^{J_{k_i}(\lambda_i)})$，而 $\det(e^{J_k(\lambda)}) = (e^\lambda)^k = e^{k\lambda}$。因此 $\det(e^A) = e^{\sum k_i\lambda_i} = e^{\operatorname{tr}(A)}$。

    **(6)** $(PAP^{-1})^k = PA^kP^{-1}$，代入级数即得。$\blacksquare$

!!! note "注"
    **注意：** 当 $AB \neq BA$ 时，一般 $e^{A+B} \neq e^A e^B$。这是矩阵指数与标量指数的一个重要区别。例如取 $A = \begin{pmatrix}0&1\\0&0\end{pmatrix}$，$B = \begin{pmatrix}0&0\\1&0\end{pmatrix}$，可以验证 $e^{A+B} \neq e^A e^B$。

!!! example "例 13.3"
    计算 $e^A$，其中 $A = \begin{pmatrix}0&-\theta\\\theta&0\end{pmatrix}$。

    **解：** 注意到
    $$
    A^2 = \begin{pmatrix}-\theta^2&0\\0&-\theta^2\end{pmatrix} = -\theta^2 I, \quad A^3 = -\theta^2 A, \quad A^4 = \theta^4 I, \ldots
    $$
    一般地，$A^{2k} = (-1)^k \theta^{2k} I$，$A^{2k+1} = (-1)^k \theta^{2k} A$。

    $$
    e^A = \sum_{k=0}^{\infty} \frac{A^k}{k!} = \left(\sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k}}{(2k)!}\right) I + \left(\sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k}}{(2k+1)!}\right) \frac{A}{\theta}
    $$
    $$
    = \cos\theta \cdot I + \frac{\sin\theta}{\theta} \cdot A = \begin{pmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix}.
    $$

    这正是旋转矩阵！$A$ 是反对称矩阵，$e^A$ 是正交矩阵。

---

## 13.4 矩阵指数的计算

<div class="context-flow" markdown>

三条路线：**对角化** $e^A = Pe^{\Lambda}P^{-1}$ · **Jordan 形** $e^{J_k(\lambda)t} = e^{\lambda t}\sum \frac{t^j}{j!}N^j$ · **Cayley-Hamilton** 法用特征值条件定系数

</div>

矩阵指数的计算是应用中的核心问题。本节介绍几种主要方法。

### 13.4.1 对角矩阵

!!! theorem "定理 13.6 (对角矩阵的指数)"
    若 $A = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，则
    $$
    e^A = \operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n}).
    $$
    更一般地，若 $A = P\operatorname{diag}(\lambda_1, \ldots, \lambda_n)P^{-1}$（$A$ 可对角化），则
    $$
    e^A = P\operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n})P^{-1}.
    $$

??? proof "证明"
    $A^k = \operatorname{diag}(\lambda_1^k, \ldots, \lambda_n^k)$，因此
    $$
    e^A = \sum_{k=0}^{\infty}\frac{A^k}{k!} = \operatorname{diag}\left(\sum_{k=0}^{\infty}\frac{\lambda_1^k}{k!}, \ldots, \sum_{k=0}^{\infty}\frac{\lambda_n^k}{k!}\right) = \operatorname{diag}(e^{\lambda_1}, \ldots, e^{\lambda_n}). \qquad \blacksquare
    $$

!!! example "例 13.4"
    计算 $e^{At}$，其中 $A = \begin{pmatrix}1&0\\0&-2\end{pmatrix}$。

    **解：** $A$ 是对角矩阵，因此
    $$
    e^{At} = \begin{pmatrix}e^t&0\\0&e^{-2t}\end{pmatrix}.
    $$

### 13.4.2 Jordan 块

!!! theorem "定理 13.7 (Jordan 块的指数)"
    对 Jordan 块 $J_k(\lambda)$，
    $$
    e^{J_k(\lambda)t} = e^{\lambda t}\begin{pmatrix}
    1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\
    0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & 1 & t \\
    0 & \cdots & 0 & 0 & 1
    \end{pmatrix}.
    $$

??? proof "证明"
    $J_k(\lambda)t = \lambda t I + N_k t$，其中 $\lambda t I$ 和 $N_k t$ 可交换。因此
    $$
    e^{J_k(\lambda)t} = e^{\lambda t I}e^{N_k t} = e^{\lambda t}\sum_{j=0}^{k-1}\frac{(N_k t)^j}{j!} = e^{\lambda t}\sum_{j=0}^{k-1}\frac{t^j}{j!}N_k^j.
    $$
    而 $N_k^j$ 是第 $j$ 条超对角线全为 1 的矩阵，代入即得结论。$\blacksquare$

!!! example "例 13.5"
    计算 $e^{At}$，其中 $A = \begin{pmatrix}3&1&0\\0&3&1\\0&0&3\end{pmatrix} = J_3(3)$。

    **解：** 由定理 13.7，
    $$
    e^{At} = e^{3t}\begin{pmatrix}1&t&\frac{t^2}{2}\\0&1&t\\0&0&1\end{pmatrix}.
    $$

### 13.4.3 Cayley-Hamilton 方法

!!! definition "定义 13.5 (Cayley-Hamilton 方法计算矩阵函数)"
    由 Cayley-Hamilton 定理，$A^n$ 可以表示为 $I, A, \ldots, A^{n-1}$ 的线性组合。因此任意矩阵函数 $f(A)$ 可以表示为
    $$
    f(A) = \alpha_0 I + \alpha_1 A + \cdots + \alpha_{n-1} A^{n-1},
    $$
    其中系数 $\alpha_0, \ldots, \alpha_{n-1}$ 由以下条件确定：对 $A$ 的每个特征值 $\lambda_i$（代数重数 $m_i$），
    $$
    f^{(j)}(\lambda_i) = \alpha_0^{(j)} + \alpha_1 \cdot j! + \cdots \quad (j = 0, 1, \ldots, m_i - 1),
    $$
    即 $\alpha$ 多项式在 $\lambda_i$ 处的函数值及导数值与 $f$ 的一致。

!!! example "例 13.6"
    用 Cayley-Hamilton 方法计算 $e^{At}$，其中 $A = \begin{pmatrix}2&1\\0&2\end{pmatrix}$。

    **解：** 特征值 $\lambda = 2$（代数重数 2）。$n = 2$，设
    $$
    e^{At} = \alpha_0(t) I + \alpha_1(t) A.
    $$

    由条件 $f(\lambda) = e^{\lambda t}$ 在 $\lambda = 2$ 处匹配到一阶导数：

    - $f(2) = e^{2t}$：$\alpha_0 + 2\alpha_1 = e^{2t}$；
    - $f'(2) = te^{2t}$：$\alpha_1 = te^{2t}$。

    解得 $\alpha_1 = te^{2t}$，$\alpha_0 = e^{2t} - 2te^{2t}$。

    $$
    e^{At} = (e^{2t} - 2te^{2t})I + te^{2t}A = e^{2t}\begin{pmatrix}1-2t&0\\0&1-2t\end{pmatrix} + te^{2t}\begin{pmatrix}2&1\\0&2\end{pmatrix}
    $$
    $$
    = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}.
    $$

---

## 13.5 矩阵指数与微分方程

<div class="context-flow" markdown>

$\mathbf{x}' = A\mathbf{x}$ → $\mathbf{x}(t) = e^{At}\mathbf{x}_0$（唯一解） → 非齐次用**常数变易法** → Jordan 形决定解的渐近行为（$e^{\lambda t}$ × 多项式）

</div>

矩阵指数的最重要应用是求解线性常系数微分方程组。

<div class="context-flow" markdown>

**洞察**：唯一性证明的巧妙之处——令 $\mathbf{z}(t) = e^{-At}\mathbf{y}(t)$，利用 $\mathbf{z}'=0$ 即得 $\mathbf{y} = e^{At}\mathbf{x}_0$

</div>

!!! theorem "定理 13.8 (齐次线性微分方程组的解)"
    微分方程组
    $$
    \mathbf{x}'(t) = A\mathbf{x}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$
    的唯一解为
    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0.
    $$

??? proof "证明"
    **存在性：** 令 $\mathbf{x}(t) = e^{At}\mathbf{x}_0$。则
    $$
    \mathbf{x}'(t) = \frac{d}{dt}e^{At}\mathbf{x}_0 = Ae^{At}\mathbf{x}_0 = A\mathbf{x}(t),
    $$
    其中 $\frac{d}{dt}e^{At} = \sum_{k=1}^{\infty}\frac{kA^k t^{k-1}}{k!} = A\sum_{k=1}^{\infty}\frac{A^{k-1}t^{k-1}}{(k-1)!} = Ae^{At}$。

    且 $\mathbf{x}(0) = e^{0}\mathbf{x}_0 = \mathbf{x}_0$。

    **唯一性：** 设 $\mathbf{y}(t)$ 也是解。令 $\mathbf{z}(t) = e^{-At}\mathbf{y}(t)$，则
    $$
    \mathbf{z}'(t) = -Ae^{-At}\mathbf{y}(t) + e^{-At}\mathbf{y}'(t) = -Ae^{-At}\mathbf{y}(t) + e^{-At}A\mathbf{y}(t) = 0.
    $$
    因此 $\mathbf{z}(t) = \mathbf{z}(0) = \mathbf{x}_0$，即 $\mathbf{y}(t) = e^{At}\mathbf{x}_0$。$\blacksquare$

!!! theorem "定理 13.9 (非齐次线性微分方程组)"
    微分方程组
    $$
    \mathbf{x}'(t) = A\mathbf{x}(t) + \mathbf{f}(t), \quad \mathbf{x}(0) = \mathbf{x}_0,
    $$
    的解为（常数变易法 / variation of parameters）
    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds.
    $$

??? proof "证明"
    令 $\mathbf{x}(t) = e^{At}\mathbf{c}(t)$（常数变易法），代入方程得
    $$
    Ae^{At}\mathbf{c}(t) + e^{At}\mathbf{c}'(t) = Ae^{At}\mathbf{c}(t) + \mathbf{f}(t).
    $$
    化简得 $e^{At}\mathbf{c}'(t) = \mathbf{f}(t)$，即 $\mathbf{c}'(t) = e^{-At}\mathbf{f}(t)$。积分得
    $$
    \mathbf{c}(t) = \mathbf{x}_0 + \int_0^t e^{-As}\mathbf{f}(s)\,ds.
    $$
    因此 $\mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}\mathbf{f}(s)\,ds$。$\blacksquare$

!!! example "例 13.7"
    求微分方程组
    $$
    \begin{cases} x_1' = 3x_1 + x_2, \\ x_2' = -x_1 + x_2, \end{cases} \quad \mathbf{x}(0) = \begin{pmatrix}1\\0\end{pmatrix}
    $$
    的解。

    **解：** $A = \begin{pmatrix}3&1\\-1&1\end{pmatrix}$。特征值：$\lambda^2 - 4\lambda + 4 = (\lambda-2)^2 = 0$，$\lambda = 2$（重根）。

    $A - 2I = \begin{pmatrix}1&1\\-1&-1\end{pmatrix}$，$\operatorname{rank} = 1$，几何重数 = 1。

    Jordan 形 $J = J_2(2)$。特征向量 $\mathbf{v}_1 = \begin{pmatrix}1\\-1\end{pmatrix}$，广义特征向量 $(A-2I)\mathbf{v}_2 = \mathbf{v}_1$：
    $$
    \begin{pmatrix}1&1\\-1&-1\end{pmatrix}\mathbf{v}_2 = \begin{pmatrix}1\\-1\end{pmatrix}, \quad \Rightarrow \quad \mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix}.
    $$

    $P = \begin{pmatrix}1&1\\-1&0\end{pmatrix}$，$P^{-1} = \begin{pmatrix}0&-1\\1&1\end{pmatrix}$。

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

## 13.6 矩阵对数

<div class="context-flow" markdown>

$e^X = A$ 的逆问题 → 可逆矩阵必有对数 → 对 Jordan 块 $\log J_k(\lambda) = (\log\lambda)I + \sum \frac{(-1)^{j+1}}{j}(\lambda^{-1}N)^j$（有限和）

</div>

矩阵对数（matrix logarithm）是矩阵指数的逆运算。

!!! definition "定义 13.6 (矩阵对数 Matrix Logarithm)"
    设 $A$ 为 $n \times n$ 可逆矩阵。若存在矩阵 $X$ 使得 $e^X = A$，则称 $X$ 为 $A$ 的**矩阵对数**，记为 $X = \log A$ 或 $X = \ln A$。

!!! theorem "定理 13.10 (矩阵对数的存在性)"
    设 $A$ 为 $n \times n$ 可逆复数矩阵，则 $A$ 的矩阵对数存在。即存在矩阵 $X$ 使得 $e^X = A$。

    更精确地，若 $A$ 没有负实数特征值，则存在唯一的矩阵对数 $X$，使得 $X$ 的所有特征值的虚部都在 $(-\pi, \pi)$ 内。此 $X$ 称为 $A$ 的**主对数**（principal logarithm）。

??? proof "证明"
    **存在性（构造性证明）：**

    设 $A = PJP^{-1}$，$J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$。只需对每个 Jordan 块定义对数。

    对 $J_k(\lambda)$（$\lambda \neq 0$），写
    $$
    J_k(\lambda) = \lambda(I + \lambda^{-1}N_k) = \lambda(I + M),
    $$
    其中 $M = \lambda^{-1}N_k$ 是幂零矩阵。取
    $$
    \log J_k(\lambda) = (\log\lambda) I + \log(I + M) = (\log\lambda)I + \sum_{j=1}^{k-1}\frac{(-1)^{j+1}}{j}M^j.
    $$
    级数是有限的（因为 $M$ 幂零），且 $e^{\log J_k(\lambda)} = J_k(\lambda)$。

    令 $\log A = P \operatorname{diag}(\log J_{k_1}(\lambda_1), \ldots) P^{-1}$。$\blacksquare$

!!! example "例 13.8"
    求 $\log A$，其中 $A = \begin{pmatrix}1&1\\0&1\end{pmatrix} = J_2(1)$。

    **解：** $A = I + N$，其中 $N = \begin{pmatrix}0&1\\0&0\end{pmatrix}$。

    $$
    \log A = \log(I + N) = N - \frac{N^2}{2} + \frac{N^3}{3} - \cdots = N = \begin{pmatrix}0&1\\0&0\end{pmatrix}.
    $$
    （因为 $N^2 = 0$，级数只有一项。）

    验证：$e^N = I + N + \frac{N^2}{2!} + \cdots = I + N = A$。 $\checkmark$

---

## 13.7 矩阵平方根

<div class="context-flow" markdown>

$X^2 = A$ → 正定矩阵有唯一**正定平方根** $A^{1/2} = Q\Lambda^{1/2}Q^T$ → 也可通过 $e^{\frac{1}{2}\log A}$ 定义 → 出现在 Ch10 极分解 $P = (A^HA)^{1/2}$

</div>

!!! definition "定义 13.7 (矩阵平方根 Matrix Square Root)"
    设 $A$ 为 $n \times n$ 矩阵。若存在矩阵 $X$ 使得 $X^2 = A$，则称 $X$ 为 $A$ 的**矩阵平方根**，记为 $X = A^{1/2}$。

!!! theorem "定理 13.11 (正定矩阵的唯一正定平方根)"
    设 $A$ 为 $n \times n$ 实对称正定矩阵。则存在唯一的实对称正定矩阵 $B$ 使得 $B^2 = A$。$B$ 称为 $A$ 的**正定平方根**。

??? proof "证明"
    **存在性：** $A$ 是实对称正定矩阵，由谱定理，$A = Q\Lambda Q^T$，其中 $Q$ 是正交矩阵，$\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$\lambda_i > 0$。

    令 $B = Q\Lambda^{1/2}Q^T$，其中 $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$。则
    $$
    B^2 = Q\Lambda^{1/2}Q^T Q\Lambda^{1/2}Q^T = Q\Lambda Q^T = A.
    $$
    $B$ 是对称的：$B^T = (Q\Lambda^{1/2}Q^T)^T = Q\Lambda^{1/2}Q^T = B$。$B$ 的特征值为 $\sqrt{\lambda_i} > 0$，因此 $B$ 正定。

    **唯一性：** 设 $C$ 也是对称正定矩阵且 $C^2 = A$。由 $C$ 对称正定，设 $C = R\Gamma R^T$，$\Gamma = \operatorname{diag}(\gamma_1, \ldots, \gamma_n)$，$\gamma_i > 0$。则 $A = C^2 = R\Gamma^2 R^T$。

    由 $A$ 的谱分解的唯一性（特征值确定后特征空间确定），$\Gamma^2 = \Lambda$（特征值重新排列后）。因为 $\gamma_i > 0$，$\gamma_i = \sqrt{\lambda_i}$，故 $C = B$。$\blacksquare$

!!! definition "定义 13.8 (半正定矩阵的平方根)"
    对实对称半正定矩阵 $A$（特征值 $\lambda_i \ge 0$），正定平方根的构造推广为：$A^{1/2} = Q\operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})Q^T$。此时 $A^{1/2}$ 是半正定的（唯一的半正定平方根）。

!!! theorem "定理 13.12 (可逆矩阵平方根的存在性)"
    设 $A$ 为 $n \times n$ 可逆复数矩阵，且 $A$ 没有负实数特征值。则 $A$ 存在唯一的平方根 $A^{1/2}$，使得 $A^{1/2}$ 的所有特征值具有正实部。

??? proof "证明"
    利用主对数：令 $A^{1/2} = e^{\frac{1}{2}\log A}$，其中 $\log A$ 为主对数。则
    $$
    (A^{1/2})^2 = e^{\frac{1}{2}\log A} e^{\frac{1}{2}\log A} = e^{\log A} = A.
    $$
    $A^{1/2}$ 的特征值为 $e^{\frac{1}{2}\log\lambda_i}$，其中 $\log\lambda_i$ 的虚部在 $(-\pi, \pi)$ 内，因此 $\frac{1}{2}\log\lambda_i$ 的虚部在 $(-\frac{\pi}{2}, \frac{\pi}{2})$ 内，$e^{\frac{1}{2}\log\lambda_i}$ 的实部为正。唯一性可由此条件推出。$\blacksquare$

!!! example "例 13.9"
    求 $A = \begin{pmatrix}2&1\\1&2\end{pmatrix}$ 的正定平方根。

    **解：** 特征值 $\lambda_1 = 3$，$\lambda_2 = 1$。正交特征向量：$\mathbf{q}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$，$\mathbf{q}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\-1\end{pmatrix}$。

    $$
    A^{1/2} = Q\begin{pmatrix}\sqrt{3}&0\\0&1\end{pmatrix}Q^T = \frac{1}{2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}\begin{pmatrix}\sqrt{3}&0\\0&1\end{pmatrix}\begin{pmatrix}1&1\\1&-1\end{pmatrix}
    $$
    $$
    = \frac{1}{2}\begin{pmatrix}\sqrt{3}+1&\sqrt{3}-1\\\sqrt{3}-1&\sqrt{3}+1\end{pmatrix}.
    $$

    验证：$(A^{1/2})^2 = \frac{1}{4}\begin{pmatrix}\sqrt{3}+1&\sqrt{3}-1\\\sqrt{3}-1&\sqrt{3}+1\end{pmatrix}^2$。

    计算对角元：$\frac{1}{4}[(\sqrt{3}+1)^2 + (\sqrt{3}-1)^2] = \frac{1}{4}(4+2\sqrt{3}+4-2\sqrt{3}) = \frac{8}{4} = 2$。

    计算非对角元：$\frac{1}{4}[(\sqrt{3}+1)(\sqrt{3}-1) + (\sqrt{3}-1)(\sqrt{3}+1)] = \frac{1}{4}(2 \times 2) = 1$。

    因此 $(A^{1/2})^2 = \begin{pmatrix}2&1\\1&2\end{pmatrix} = A$。$\checkmark$

---

## 13.8 一般矩阵函数

<div class="context-flow" markdown>

统一框架：$f(J_k(\lambda))$ 的 $(p,q)$ 元素 = $\frac{f^{(q-p)}(\lambda)}{(q-p)!}$ → Cauchy 积分定义 $f(A) = \frac{1}{2\pi i}\oint f(z)(zI-A)^{-1}dz$ 与之等价 → 谱映射 $\sigma(f(A)) = f(\sigma(A))$

</div>

前面几节讨论了特定的矩阵函数（指数、对数、平方根）。本节给出矩阵函数的一般定义框架。

### 13.8.1 通过 Jordan 标准形定义

<div class="context-flow" markdown>

**洞察**：$f(J_k(\lambda))$ 中出现 $f, f', f'', \ldots, f^{(k-1)}$ ——Jordan 块的大小决定了 $f$ 需要多少阶**可微性**

</div>

!!! definition "定义 13.9 (一般矩阵函数 — Jordan 形定义)"
    设 $f$ 为定义在 $A$ 的谱上的函数，$A = PJP^{-1}$，$J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$。定义
    $$
    f(A) = P f(J) P^{-1} = P \operatorname{diag}(f(J_{k_1}(\lambda_1)), \ldots, f(J_{k_s}(\lambda_s))) P^{-1},
    $$
    其中对每个 Jordan 块
    $$
    f(J_k(\lambda)) = \begin{pmatrix}
    f(\lambda) & f'(\lambda) & \frac{f''(\lambda)}{2!} & \cdots & \frac{f^{(k-1)}(\lambda)}{(k-1)!} \\
    0 & f(\lambda) & f'(\lambda) & \cdots & \frac{f^{(k-2)}(\lambda)}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & f(\lambda) & f'(\lambda) \\
    0 & \cdots & 0 & 0 & f(\lambda)
    \end{pmatrix}.
    $$
    要求 $f$ 在每个特征值 $\lambda_i$ 处至少 $k_i - 1$ 次可微。

!!! theorem "定理 13.13 (Jordan 形定义的合理性)"
    定义 13.9 中的 $f(A)$ 不依赖于 Jordan 分解的选取（即不依赖 $P$ 的选取），因此 $f(A)$ 是良定义的。

??? proof "证明"
    设 $A = P_1 J P_1^{-1} = P_2 J P_2^{-1}$（同一个 Jordan 形 $J$，但不同的过渡矩阵）。则 $P_2^{-1}P_1$ 与 $J$ 可交换，即 $P_2^{-1}P_1 J = J P_2^{-1}P_1$。

    由于 $f(J)$ 是 $J$ 的多项式的极限（当 $f$ 是解析函数时），$P_2^{-1}P_1$ 也与 $f(J)$ 可交换。因此
    $$
    P_1 f(J) P_1^{-1} = P_2 (P_2^{-1}P_1) f(J) (P_2^{-1}P_1)^{-1} P_2^{-1} = P_2 f(J) P_2^{-1}. \qquad \blacksquare
    $$

!!! example "例 13.10"
    计算 $\sin(A)$，其中 $A = \begin{pmatrix}0&\pi\\0&0\end{pmatrix}$。

    **解：** $A = J_2(0)$，$\lambda = 0$，$k = 2$。

    $f(\lambda) = \sin(\lambda)$。$f(0) = 0$，$f'(0) = \cos(0) = 1$。

    $$
    \sin(A) = f(J_2(0)) = \begin{pmatrix}f(0)&f'(0)\\0&f(0)\end{pmatrix} = \begin{pmatrix}0&1\\0&0\end{pmatrix}.
    $$

    注意这不是逐元素的 $\sin$！$\sin\begin{pmatrix}0&\pi\\0&0\end{pmatrix} \neq \begin{pmatrix}\sin 0&\sin\pi\\\sin 0&\sin 0\end{pmatrix}$。

### 13.8.2 Cauchy 积分定义

!!! definition "定义 13.10 (矩阵函数 — Cauchy 积分定义)"
    设 $f$ 在包含 $A$ 的所有特征值的开集 $\Omega$ 上解析，$\Gamma$ 为 $\Omega$ 中包围所有特征值的简单闭曲线。定义
    $$
    f(A) = \frac{1}{2\pi i} \oint_\Gamma f(z)(zI - A)^{-1}\,dz.
    $$

!!! theorem "定理 13.14 (Cauchy 积分定义与 Jordan 形定义的等价性)"
    当 $f$ 在 $A$ 的谱的某个开邻域上解析时，定义 13.9 和定义 13.10 给出的 $f(A)$ 相同。

??? proof "证明"
    设 $A = PJP^{-1}$，$(zI - A)^{-1} = P(zI - J)^{-1}P^{-1}$。因此
    $$
    \frac{1}{2\pi i}\oint_\Gamma f(z)(zI-A)^{-1}dz = P\left(\frac{1}{2\pi i}\oint_\Gamma f(z)(zI-J)^{-1}dz\right)P^{-1}.
    $$

    对 Jordan 块 $J_k(\lambda)$，$(zI - J_k(\lambda))^{-1}$ 的 $(p,q)$ 元素（$q \ge p$）为 $\frac{1}{(z-\lambda)^{q-p+1}}$。由 Cauchy 积分公式：
    $$
    \frac{1}{2\pi i}\oint \frac{f(z)}{(z-\lambda)^{q-p+1}}dz = \frac{f^{(q-p)}(\lambda)}{(q-p)!}.
    $$
    这正是定义 13.9 中 $f(J_k(\lambda))$ 的 $(p,q)$ 元素。$\blacksquare$

!!! theorem "定理 13.15 (矩阵函数的谱映射定理)"
    设 $f$ 为矩阵 $A$ 的谱上的解析函数，则
    $$
    \sigma(f(A)) = f(\sigma(A)) = \{f(\lambda) : \lambda \in \sigma(A)\},
    $$
    即 $f(A)$ 的特征值恰好是 $A$ 的特征值经 $f$ 映射后的值。

??? proof "证明"
    设 $A = PJP^{-1}$，$f(A) = Pf(J)P^{-1}$。$f(J)$ 的对角元素为 $f(\lambda_i)$（各 Jordan 块的对角元素），这些正是 $f(A)$ 的特征值（相似变换保谱）。$\blacksquare$

!!! example "例 13.11"
    设 $A$ 的特征值为 $1, 2, 3$。求 $e^A$ 的特征值和 $\cos(A)$ 的特征值。

    **解：** 由谱映射定理：

    - $e^A$ 的特征值为 $e^1, e^2, e^3$，即 $e, e^2, e^3$。
    - $\cos(A)$ 的特征值为 $\cos 1, \cos 2, \cos 3$。

!!! example "例 13.12"
    设 $A = \begin{pmatrix}1&0&0\\0&2&1\\0&0&2\end{pmatrix}$，计算 $\sqrt{A}$（取主平方根）。

    **解：** $A = J_1(1) \oplus J_2(2)$。$f(\lambda) = \sqrt{\lambda}$，$f'(\lambda) = \frac{1}{2\sqrt{\lambda}}$。

    $$
    f(J_1(1)) = (1) = (\sqrt{1}) = (1).
    $$
    $$
    f(J_2(2)) = \begin{pmatrix} f(2) & f'(2) \\ 0 & f(2) \end{pmatrix} = \begin{pmatrix} \sqrt{2} & \frac{1}{2\sqrt{2}} \\ 0 & \sqrt{2} \end{pmatrix} = \begin{pmatrix} \sqrt{2} & \frac{\sqrt{2}}{4} \\ 0 & \sqrt{2} \end{pmatrix}.
    $$

    因此
    $$
    \sqrt{A} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & \sqrt{2} & \frac{\sqrt{2}}{4} \\ 0 & 0 & \sqrt{2} \end{pmatrix}.
    $$

    验证：$(\sqrt{A})^2$ 的 $(2,3)$ 元素 = $\sqrt{2} \cdot \frac{\sqrt{2}}{4} + \frac{\sqrt{2}}{4} \cdot \sqrt{2} = \frac{2}{4} + \frac{2}{4} = 1 = A_{23}$。$\checkmark$

---

## 本章小结

本章系统介绍了矩阵函数理论，包括：

1. **矩阵多项式** $p(A)$ 是最基础的矩阵函数，Cayley-Hamilton 定理是核心结果；
2. **矩阵幂级数**的收敛由谱半径 $\rho(A)$ 控制；
3. **矩阵指数** $e^A$ 对任意矩阵都有定义，满足 $\det(e^A) = e^{\operatorname{tr}(A)}$，但注意 $e^{A+B} = e^A e^B$ 仅在 $AB = BA$ 时成立；
4. **矩阵指数与微分方程**：$\mathbf{x}'(t) = A\mathbf{x}(t)$ 的解为 $\mathbf{x}(t) = e^{At}\mathbf{x}_0$；
5. **矩阵对数**在可逆矩阵上存在，正定矩阵有唯一的半正定对数；
6. **矩阵平方根**：正定矩阵有唯一的正定平方根 $A^{1/2} = Q\Lambda^{1/2}Q^T$；
7. **一般矩阵函数**可通过 Jordan 标准形或 Cauchy 积分定义，谱映射定理 $\sigma(f(A)) = f(\sigma(A))$ 是重要性质。

矩阵函数理论将微积分与线性代数深度融合，是现代应用数学的重要工具箱。
