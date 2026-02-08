# 第 24 章 矩阵流形

<div class="context-flow" markdown>

**前置**：正交矩阵(Ch7) · SVD(Ch8) · 正定矩阵(Ch7)

**脉络**：矩阵约束集 → 光滑流形 → Riemannian 几何(测地线/梯度) → 流形优化

**核心对象**：$O(n)$ / $\text{St}(k,n)$ / $\text{Gr}(k,n)$ / $\mathcal{P}(n)$ ——将约束优化转化为流形上的无约束优化 → 链接 Ch25

**延伸**：矩阵流形在计算机视觉（本质矩阵估计、形状分析）、机器人学（姿态估计 $SO(3)$）、医学影像（扩散张量 MRI）、机器学习（低秩优化、子空间跟踪）中有广泛应用

</div>

矩阵流形（Matrix Manifold）是指由满足特定约束条件的矩阵所构成的光滑流形。这些流形在控制理论、信号处理、计算机视觉、机器学习等领域中无处不在。本章将系统介绍若干重要的矩阵流形，包括一般线性群、正交群、Stiefel 流形、Grassmann 流形和正定矩阵流形，以及它们上的微分几何结构与优化方法。

---

## 24.1 矩阵流形基本概念

<div class="context-flow" markdown>

**基础工具箱**：嵌入子流形($h^{-1}(0)$, 隐函数定理) → 切空间($\ker Dh$) → Riemannian 度量(切空间上的内积) → 测地线(局部最短路) → 指数映射

</div>

我们首先回顾流形论的基本语言，并说明矩阵集合如何自然地获得流形结构。

!!! definition "定义 24.1 (光滑流形)"
    一个 **光滑流形**（smooth manifold）$\mathcal{M}$ 是一个 Hausdorff 的第二可数拓扑空间，配备一个极大光滑图册 $\{(U_\alpha, \varphi_\alpha)\}$，其中每个 $U_\alpha$ 是 $\mathcal{M}$ 的开子集，$\varphi_\alpha: U_\alpha \to \mathbb{R}^d$ 为同胚，且转移映射 $\varphi_\beta \circ \varphi_\alpha^{-1}$ 为光滑映射。$d$ 称为 $\mathcal{M}$ 的维数。

!!! definition "定义 24.2 (嵌入子流形)"
    设 $\mathcal{M}$ 为 $\mathbb{R}^{n \times p}$ 的子集。若 $\mathcal{M}$ 可以表示为

    $$
    \mathcal{M} = \{ X \in \mathbb{R}^{n \times p} : h(X) = 0 \},
    $$

    其中 $h: \mathbb{R}^{n \times p} \to \mathbb{R}^m$ 为光滑映射，且对所有 $X \in \mathcal{M}$，$Dh(X)$ 为满秩映射（即 $\operatorname{rank}(Dh(X)) = m$），则 $\mathcal{M}$ 为 $\mathbb{R}^{n \times p}$ 的一个 **嵌入子流形**（embedded submanifold），维数为 $np - m$。

!!! definition "定义 24.3 (切空间)"
    设 $\mathcal{M}$ 为光滑流形，$X \in \mathcal{M}$。$\mathcal{M}$ 在 $X$ 处的 **切空间**（tangent space）$T_X \mathcal{M}$ 是所有经过 $X$ 的光滑曲线在 $X$ 处的速度向量的集合：

    $$
    T_X \mathcal{M} = \left\{ \gamma'(0) : \gamma: (-\epsilon, \epsilon) \to \mathcal{M} \text{ 光滑}, \, \gamma(0) = X \right\}.
    $$

    若 $\mathcal{M} = h^{-1}(0)$ 为嵌入子流形，则 $T_X \mathcal{M} = \ker(Dh(X))$。

!!! definition "定义 24.4 (Riemannian 度量)"
    流形 $\mathcal{M}$ 上的一个 **Riemannian 度量**（Riemannian metric）是一个光滑地依赖于底点的内积族 $\langle \cdot, \cdot \rangle_X: T_X\mathcal{M} \times T_X\mathcal{M} \to \mathbb{R}$，$X \in \mathcal{M}$。配备 Riemannian 度量的流形称为 Riemannian 流形。

!!! definition "定义 24.5 (测地线)"
    Riemannian 流形 $(\mathcal{M}, \langle \cdot, \cdot \rangle)$ 上的 **测地线**（geodesic）是局部最短曲线。形式上，$\gamma: [0, 1] \to \mathcal{M}$ 为测地线当且仅当 $\nabla_{\dot\gamma} \dot\gamma = 0$，其中 $\nabla$ 为 Levi-Civita 联络。**指数映射**（exponential map）$\operatorname{Exp}_X: T_X\mathcal{M} \to \mathcal{M}$ 将切向量 $\xi$ 映射到从 $X$ 出发、初始速度为 $\xi$ 的测地线在 $t = 1$ 时的端点。

!!! theorem "定理 24.1 (隐函数定理与子流形)"
    设 $h: \mathbb{R}^{n \times p} \to \mathbb{R}^m$ 为光滑映射，$\mathbf{0}$ 为 $h$ 的正则值（即对所有 $X \in h^{-1}(\mathbf{0})$，$Dh(X)$ 满秩）。则 $\mathcal{M} = h^{-1}(\mathbf{0})$ 为 $\mathbb{R}^{n \times p}$ 的光滑嵌入子流形，维数为 $np - m$。

??? proof "证明"
    这是隐函数定理的直接推论。对 $X_0 \in \mathcal{M}$，$Dh(X_0)$ 满秩意味着存在 $X_0$ 的邻域 $U$ 和光滑映射 $\psi$，使得 $\mathcal{M} \cap U$ 可以参数化为 $np - m$ 个自由变量的光滑函数的图像。这些局部参数化构成 $\mathcal{M}$ 的光滑图册。$\blacksquare$

!!! example "例 24.1"
    **单位球面作为子流形。**

    $S^{n-1} = \{\mathbf{x} \in \mathbb{R}^n : \|\mathbf{x}\|^2 = 1\} = h^{-1}(0)$，其中 $h(\mathbf{x}) = \mathbf{x}^T\mathbf{x} - 1$。$Dh(\mathbf{x}) = 2\mathbf{x}^T$，对 $\mathbf{x} \in S^{n-1}$ 显然满秩。切空间 $T_{\mathbf{x}} S^{n-1} = \{\mathbf{v} : \mathbf{x}^T \mathbf{v} = 0\}$，即与 $\mathbf{x}$ 正交的向量。

---

## 24.2 一般线性群 $GL(n)$

<div class="context-flow" markdown>

**母群**：$GL(n) = \{\det \ne 0\}$ 是 $\mathbb{R}^{n^2}$ 的开集 → Lie 代数 $\mathfrak{gl}(n) = \mathbb{R}^{n \times n}$ → 矩阵指数 $e^A$ 连接 Lie 代数与 Lie 群

**所有矩阵 Lie 群都是 $GL(n)$ 的闭子群**（Cartan 定理）

</div>

!!! definition "定义 24.6 (一般线性群)"
    **一般线性群**（General Linear Group）$GL(n, \mathbb{R})$（简记 $GL(n)$）是所有 $n \times n$ 可逆实矩阵的集合：

    $$
    GL(n) = \{ A \in \mathbb{R}^{n \times n} : \det(A) \ne 0 \}.
    $$

    $GL(n)$ 关于矩阵乘法构成群，且为 $\mathbb{R}^{n \times n}$ 的开子集，因此是 $n^2$ 维光滑流形。类似地，$GL(n, \mathbb{C})$ 为所有 $n \times n$ 可逆复矩阵之群。

!!! definition "定义 24.7 (Lie 代数 $\mathfrak{gl}(n)$)"
    $GL(n)$ 的 **Lie 代数**（Lie algebra）$\mathfrak{gl}(n)$ 是 $GL(n)$ 在单位元 $I$ 处的切空间，即全体 $n \times n$ 实矩阵：

    $$
    \mathfrak{gl}(n) = T_I GL(n) = \mathbb{R}^{n \times n}.
    $$

    $\mathfrak{gl}(n)$ 上的 Lie 括号定义为矩阵交换子 $[A, B] = AB - BA$。

!!! theorem "定理 24.2 (矩阵指数映射)"
    映射 $\exp: \mathfrak{gl}(n) \to GL(n)$，$A \mapsto e^A = \sum_{k=0}^{\infty} \frac{A^k}{k!}$ 是光滑映射，且在 $A = 0$ 处的微分为恒等映射。因此 $\exp$ 在 $0$ 的邻域内是微分同胚。

??? proof "证明"
    矩阵指数级数 $e^A = \sum_{k=0}^\infty \frac{A^k}{k!}$ 对所有 $A \in \mathbb{R}^{n \times n}$ 绝对收敛（因为 $\|A^k/k!\| \le \|A\|^k/k!$，级数的和不超过 $e^{\|A\|}$）。光滑性由幂级数的逐项微分得到。在 $A = 0$ 处，$D\exp(0)(H) = \frac{d}{dt}\Big|_{t=0} e^{tH} = H$，故微分为恒等映射。由反函数定理，$\exp$ 在 $0$ 的邻域内为微分同胚。$\blacksquare$

!!! theorem "定理 24.3 ($GL(n)$ 上的左不变向量场)"
    $GL(n)$ 上的每个左不变向量场由其在 $I$ 处的值唯一确定。具体地，对 $A \in \mathfrak{gl}(n)$，左不变向量场 $\tilde{A}$ 定义为 $\tilde{A}_g = gA$（$g \in GL(n)$），且 $[\tilde{A}, \tilde{B}] = \widetilde{[A, B]}$。

??? proof "证明"
    左平移 $L_g: GL(n) \to GL(n)$，$h \mapsto gh$ 的微分为 $(dL_g)_h(V) = gV$。向量场 $\tilde{A}$ 满足 $(dL_g)\tilde{A}_h = g(hA) = (gh)A = \tilde{A}_{gh}$，故左不变。Lie 括号 $[\tilde{A}, \tilde{B}]_I = \frac{d}{dt}\Big|_{t=0} \frac{d}{ds}\Big|_{s=0} (e^{tA} e^{sB} e^{-tA}) = AB - BA = [A, B]$。$\blacksquare$

!!! example "例 24.2"
    **$GL(2)$ 上的指数映射。**

    设 $A = \begin{pmatrix} 0 & -\theta \\ \theta & 0 \end{pmatrix}$，则

    $$
    e^A = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \in SO(2).
    $$

    这表明反对称矩阵通过指数映射生成旋转矩阵。

!!! example "例 24.3"
    **对角矩阵的指数。**

    设 $D = \operatorname{diag}(d_1, \ldots, d_n)$，则 $e^D = \operatorname{diag}(e^{d_1}, \ldots, e^{d_n})$。$e^D$ 的特征值都是正的，因此 $e^D$ 是正定对角矩阵。

---

## 24.3 正交群与酉群

<div class="context-flow" markdown>

**约束 → 流形**：$Q^TQ = I$ → $O(n)$，$\dim = \frac{n(n-1)}{2}$ · Lie 代数 = **反对称矩阵** $\text{Skew}(n)$ → 测地线 $\gamma(t) = Qe^{t\Omega}$ · Ch23 GOE 的不变性正是 $O(n)$-共轭不变

</div>

正交群和酉群是最重要的矩阵 Lie 群，在几何、物理和工程中有广泛应用。

!!! definition "定义 24.8 (正交群)"
    **正交群**（Orthogonal Group）$O(n)$ 是所有 $n \times n$ 正交矩阵的集合：

    $$
    O(n) = \{ Q \in \mathbb{R}^{n \times n} : Q^T Q = I_n \}.
    $$

    $O(n)$ 是紧 Lie 群，维数为 $\frac{n(n-1)}{2}$。**特殊正交群**（Special Orthogonal Group）$SO(n) = \{Q \in O(n) : \det(Q) = 1\}$ 是 $O(n)$ 的单位连通分量。

!!! definition "定义 24.9 (酉群)"
    **酉群**（Unitary Group）$U(n)$ 是所有 $n \times n$ 酉矩阵的集合：

    $$
    U(n) = \{ U \in \mathbb{C}^{n \times n} : U^* U = I_n \}.
    $$

    $U(n)$ 是紧 Lie 群，（实）维数为 $n^2$。**特殊酉群**（Special Unitary Group）$SU(n) = \{U \in U(n) : \det(U) = 1\}$，维数为 $n^2 - 1$。

!!! theorem "定理 24.4 ($O(n)$ 和 $SO(n)$ 的切空间和 Lie 代数)"
    $O(n)$ 在 $I$ 处的切空间（即 Lie 代数）为反对称矩阵空间：

    $$
    \mathfrak{o}(n) = T_I O(n) = \operatorname{Skew}(n) = \{ \Omega \in \mathbb{R}^{n \times n} : \Omega^T = -\Omega \}.
    $$

    更一般地，在 $Q \in O(n)$ 处，$T_Q O(n) = \{ Q\Omega : \Omega \in \operatorname{Skew}(n) \}$。$SO(n)$ 的 Lie 代数也是 $\operatorname{Skew}(n)$（因为 $\det$ 条件不影响切空间）。

??? proof "证明"
    $O(n)$ 由约束 $h(Q) = Q^T Q - I = 0$ 定义。$Dh(Q)(V) = V^T Q + Q^T V$。在 $Q = I$ 处，$Dh(I)(V) = V^T + V$。因此

    $$
    T_I O(n) = \ker(Dh(I)) = \{ V : V^T + V = 0 \} = \operatorname{Skew}(n).
    $$

    $\operatorname{Skew}(n)$ 的维数为 $\frac{n(n-1)}{2}$，与 $O(n)$ 的维数一致。

    在一般点 $Q$ 处，令 $\gamma(t) = Qe^{t\Omega}$（$\Omega \in \operatorname{Skew}(n)$），则 $\gamma(0) = Q$，$\gamma(t)^T\gamma(t) = e^{t\Omega^T}Q^TQe^{t\Omega} = e^{-t\Omega}e^{t\Omega} = I$，故 $\gamma(t) \in O(n)$，$\gamma'(0) = Q\Omega$。这说明 $T_Q O(n) \supseteq \{Q\Omega : \Omega \in \operatorname{Skew}(n)\}$。维数计数表明这是等式。$\blacksquare$

!!! theorem "定理 24.5 ($U(n)$ 的 Lie 代数)"
    $U(n)$ 的 Lie 代数为反 Hermite 矩阵空间：

    $$
    \mathfrak{u}(n) = T_I U(n) = \{ \Omega \in \mathbb{C}^{n \times n} : \Omega^* = -\Omega \}.
    $$

    $SU(n)$ 的 Lie 代数为 $\mathfrak{su}(n) = \{ \Omega \in \mathfrak{u}(n) : \operatorname{tr}(\Omega) = 0 \}$。

??? proof "证明"
    与 $O(n)$ 的证明完全类似。约束为 $U^* U = I$，微分得 $V^* U + U^* V = 0$，在 $U = I$ 时 $V^* + V = 0$，即 $V$ 为反 Hermite 矩阵。对于 $SU(n)$，额外条件 $\det(U) = 1$ 的微分为 $\operatorname{tr}(U^{-1}V) = \operatorname{tr}(V) = 0$（在 $U = I$ 处）。$\blacksquare$

!!! theorem "定理 24.6 ($O(n)$ 上的测地线)"
    赋予 $O(n)$ 双不变度量 $\langle \xi, \eta \rangle_Q = \operatorname{tr}(\xi^T \eta)$（$\xi, \eta \in T_Q O(n)$），则从 $Q$ 出发、初始速度 $Q\Omega$（$\Omega \in \operatorname{Skew}(n)$）的测地线为

    $$
    \gamma(t) = Q e^{t\Omega}.
    $$

    $O(n)$ 上两点 $Q_1, Q_2$ 之间的测地距离为

    $$
    d(Q_1, Q_2) = \|\log(Q_1^T Q_2)\|_F = \left(\sum_{k=1}^{\lfloor n/2 \rfloor} \theta_k^2\right)^{1/2},
    $$

    其中 $\theta_k$ 为 $Q_1^T Q_2$ 的旋转角。

??? proof "证明"
    在双不变度量下，$O(n)$ 上的 Levi-Civita 联络为 $\nabla_{\tilde{A}}\tilde{B} = \frac{1}{2}[\tilde{A}, \tilde{B}]$（其中 $\tilde{A}, \tilde{B}$ 为左不变向量场）。曲线 $\gamma(t) = Qe^{t\Omega}$ 的速度为 $\dot\gamma = Qe^{t\Omega}\Omega$，对应的左不变速度 $\gamma^{-1}\dot\gamma = \Omega$ 为常数，因此 $\nabla_{\dot\gamma}\dot\gamma = 0$，即 $\gamma$ 为测地线。距离公式由测地线长度 $L = \int_0^1 \|\dot\gamma\|dt = \|\Omega\|_F$ 得到，其中 $e^\Omega = Q_1^T Q_2$。$\blacksquare$

!!! example "例 24.4"
    **$SO(3)$ 中的旋转。**

    $SO(3)$ 的 Lie 代数 $\mathfrak{so}(3)$ 由 $3 \times 3$ 反对称矩阵组成，维数为 $3$。标准基为

    $$
    L_1 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & 1 & 0 \end{pmatrix}, \quad
    L_2 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ -1 & 0 & 0 \end{pmatrix}, \quad
    L_3 = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
    $$

    每个 $\Omega \in \mathfrak{so}(3)$ 可以写成 $\Omega = \theta_1 L_1 + \theta_2 L_2 + \theta_3 L_3$，对应绕轴 $(\theta_1, \theta_2, \theta_3)$ 旋转角度 $\|\boldsymbol{\theta}\|$。Rodrigues 公式给出

    $$
    e^\Omega = I + \frac{\sin\theta}{\theta}\Omega + \frac{1 - \cos\theta}{\theta^2}\Omega^2, \quad \theta = \|\boldsymbol{\theta}\|.
    $$

!!! example "例 24.5"
    **$O(n)$ 的维数计算。**

    $O(n)$ 由约束 $Q^T Q = I$ 定义，共 $n^2$ 个矩阵元素，$\frac{n(n+1)}{2}$ 个独立约束（因为 $Q^T Q$ 是对称的）。因此

    $$
    \dim O(n) = n^2 - \frac{n(n+1)}{2} = \frac{n(n-1)}{2}.
    $$

    例如 $\dim O(3) = 3$，$\dim O(4) = 6$。

---

## 24.4 Stiefel 流形

<div class="context-flow" markdown>

**推广链**：$S^{n-1} = \text{St}(1,n) \subset \text{St}(k,n) \subset O(n) = \text{St}(n,n)$ · 切空间：$X^TZ$ 反对称 → 投影 $\Pi_X(Z) = Z - X\text{sym}(X^TZ)$ → 链接 §24.8 Riemannian 梯度下降

</div>

Stiefel 流形推广了正交群，是约束优化中的重要对象。

!!! definition "定义 24.10 (Stiefel 流形)"
    **Stiefel 流形**（Stiefel manifold）$V_k(\mathbb{R}^n)$（$1 \le k \le n$）是 $\mathbb{R}^n$ 中所有 $k$-标准正交组的集合：

    $$
    V_k(\mathbb{R}^n) = \operatorname{St}(k, n) = \{ X \in \mathbb{R}^{n \times k} : X^T X = I_k \}.
    $$

    当 $k = n$ 时，$V_n(\mathbb{R}^n) = O(n)$；当 $k = 1$ 时，$V_1(\mathbb{R}^n) = S^{n-1}$。$\operatorname{St}(k, n)$ 是紧流形，维数为 $nk - \frac{k(k+1)}{2}$。

!!! theorem "定理 24.7 (Stiefel 流形的切空间)"
    $\operatorname{St}(k, n)$ 在 $X$ 处的切空间为

    $$
    T_X \operatorname{St}(k, n) = \{ Z \in \mathbb{R}^{n \times k} : X^T Z + Z^T X = 0 \}.
    $$

    等价地，$Z \in T_X \operatorname{St}(k, n)$ 当且仅当 $X^T Z$ 是反对称矩阵。

??? proof "证明"
    约束为 $h(X) = X^T X - I_k = 0$。微分：$Dh(X)(Z) = X^T Z + Z^T X$。因此

    $$
    T_X \operatorname{St}(k, n) = \ker(Dh(X)) = \{ Z : X^T Z + Z^T X = 0 \}.
    $$

    验证满秩条件：$Dh(X)$ 是从 $\mathbb{R}^{n \times k}$ 到 $\operatorname{Sym}(k)$（$k \times k$ 对称矩阵空间）的线性映射。给定任意对称矩阵 $S$，取 $Z = \frac{1}{2}XS$，则 $X^T Z + Z^T X = S$，故 $Dh(X)$ 满射。$\blacksquare$

!!! theorem "定理 24.8 (Stiefel 流形上的典范度量与测地线)"
    赋予 $\operatorname{St}(k, n)$ 从欧氏空间 $\mathbb{R}^{n \times k}$ 继承的度量（典范度量）$\langle Z_1, Z_2 \rangle = \operatorname{tr}(Z_1^T Z_2)$。正交投影到切空间 $T_X \operatorname{St}(k, n)$ 的公式为

    $$
    \Pi_X(Z) = Z - X \operatorname{sym}(X^T Z),
    $$

    其中 $\operatorname{sym}(A) = \frac{A + A^T}{2}$。

    在此度量下，从 $X$ 出发、初始速度 $Z \in T_X \operatorname{St}(k, n)$ 的测地线可以通过如下方法计算：设 $A = X^T Z$（反对称），$QR = (I - XX^T)Z$ 为 QR 分解（$Q \in \mathbb{R}^{n \times k}$，$R \in \mathbb{R}^{k \times k}$），则

    $$
    \gamma(t) = \begin{pmatrix} X & Q \end{pmatrix} \exp\!\left(t \begin{pmatrix} A & -R^T \\ R & 0 \end{pmatrix}\right) \begin{pmatrix} I_k \\ 0 \end{pmatrix}.
    $$

??? proof "证明"
    投影公式：对 $Z \in \mathbb{R}^{n \times k}$，分解 $Z = \Pi_X(Z) + X \cdot \operatorname{sym}(X^T Z)$。验证 $\Pi_X(Z) \in T_X\operatorname{St}(k,n)$：$X^T\Pi_X(Z) = X^TZ - \operatorname{sym}(X^TZ)$ 是反对称的。残差 $X \cdot \operatorname{sym}(X^T Z)$ 与切空间正交（在典范内积下）。

    测地线公式的推导需要解联络方程 $\nabla_{\dot\gamma}\dot\gamma = 0$。将速度分解为 $X$ 方向和 $X$ 的正交补方向，利用矩阵指数表示旋转，可以将测地线方程化为有限维矩阵 ODE，其解即为上述公式。详细推导见 Edelman-Arias-Smith (1998)。$\blacksquare$

!!! example "例 24.6"
    **$\operatorname{St}(2, 3)$ 的维数和切空间。**

    $\operatorname{St}(2, 3) = \{X \in \mathbb{R}^{3 \times 2} : X^TX = I_2\}$，维数 $= 3 \times 2 - \frac{2 \times 3}{2} = 3$。

    取 $X = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$，则切空间为满足 $X^TZ$ 反对称的 $Z \in \mathbb{R}^{3 \times 2}$：

    $$
    T_X \operatorname{St}(2, 3) = \left\{ \begin{pmatrix} 0 & -a \\ a & 0 \\ b & c \end{pmatrix} : a, b, c \in \mathbb{R} \right\},
    $$

    确实是 $3$ 维空间。

---

## 24.5 Grassmann 流形

<div class="context-flow" markdown>

**商空间**：$\text{Gr}(k,n) = \text{St}(k,n)/O(k)$——参数化 $k$ 维子空间而非基 · 距离 = **主角** $\theta_i = \arccos\sigma_i(X_1^TX_2)$（SVD, Ch8）

**应用**：PCA = $\text{Gr}(k,n)$ 上的优化（Ch25）· 子空间跟踪/计算机视觉

</div>

Grassmann 流形参数化固定维数的子空间，在 PCA、子空间跟踪等问题中自然出现。

!!! definition "定义 24.11 (Grassmann 流形)"
    **Grassmann 流形**（Grassmann manifold）$\operatorname{Gr}(k, n)$ 是 $\mathbb{R}^n$ 中所有 $k$ 维子空间的集合：

    $$
    \operatorname{Gr}(k, n) = \{ \mathcal{V} \subseteq \mathbb{R}^n : \dim \mathcal{V} = k \}.
    $$

    $\operatorname{Gr}(k, n)$ 可以视为 Stiefel 流形的商空间：

    $$
    \operatorname{Gr}(k, n) \cong \operatorname{St}(k, n) / O(k),
    $$

    其中等价关系为 $X \sim XQ$（$Q \in O(k)$），因为 $X$ 和 $XQ$ 张成相同的子空间。$\operatorname{Gr}(k, n)$ 的维数为 $k(n - k)$。

!!! definition "定义 24.12 (Grassmann 流形的投影表示)"
    每个 $k$ 维子空间 $\mathcal{V} \in \operatorname{Gr}(k, n)$ 可以用其正交投影矩阵 $P = XX^T$（$X \in \operatorname{St}(k, n)$ 为 $\mathcal{V}$ 的标准正交基）唯一表示。因此

    $$
    \operatorname{Gr}(k, n) \cong \{ P \in \mathbb{R}^{n \times n} : P^2 = P, \, P^T = P, \, \operatorname{tr}(P) = k \}.
    $$

!!! theorem "定理 24.9 (Grassmann 流形的切空间)"
    在投影表示 $P = XX^T$ 下，$\operatorname{Gr}(k, n)$ 在 $P$ 处的切空间为

    $$
    T_P \operatorname{Gr}(k, n) = \{ \Delta \in \mathbb{R}^{n \times n} : \Delta = \Delta^T, \, P\Delta + \Delta P = \Delta \}.
    $$

    等价地，在标准正交基表示下，$T_{[X]} \operatorname{Gr}(k, n)$ 中的水平切向量为

    $$
    \mathcal{H}_X = \{ Z \in \mathbb{R}^{n \times k} : X^T Z = 0 \},
    $$

    即列空间与 $X$ 的列空间正交的矩阵。

??? proof "证明"
    在商空间的框架下，Stiefel 流形上的切向量 $Z \in T_X \operatorname{St}(k, n)$ 分解为垂直分量（沿 $O(k)$ 轨道的方向，即 $XA$，$A$ 反对称）和水平分量（与轨道正交的方向）。在典范度量下，水平空间为

    $$
    \mathcal{H}_X = \{ Z \in T_X \operatorname{St}(k, n) : X^T Z = 0 \} = \{ Z \in \mathbb{R}^{n \times k} : X^T Z = 0 \}.
    $$

    维数为 $k(n - k) = \dim \operatorname{Gr}(k, n)$。$\blacksquare$

!!! theorem "定理 24.10 (Grassmann 流形上的测地线与距离)"
    赋予 $\operatorname{Gr}(k, n)$ 从 Stiefel 流形商结构继承的 Riemannian 度量。两个子空间 $\mathcal{V}_1, \mathcal{V}_2 \in \operatorname{Gr}(k, n)$ 之间的测地距离为

    $$
    d(\mathcal{V}_1, \mathcal{V}_2) = \|\boldsymbol{\theta}\|_2 = \left(\sum_{i=1}^{k} \theta_i^2\right)^{1/2},
    $$

    其中 $\theta_1, \ldots, \theta_k \in [0, \pi/2]$ 为 $\mathcal{V}_1$ 和 $\mathcal{V}_2$ 之间的 **主角**（principal angles），定义为 $\cos\theta_i = \sigma_i(X_1^T X_2)$，$\sigma_i$ 为奇异值。

    从 $[X_1]$ 沿水平方向 $Z$ 的测地线为：设 $Z = U\Sigma V^T$ 为紧 SVD，则

    $$
    \gamma(t) = [X_1 V \cos(\Sigma t) + U \sin(\Sigma t)].
    $$

??? proof "证明"
    设 $X_1, X_2 \in \operatorname{St}(k, n)$ 分别为 $\mathcal{V}_1, \mathcal{V}_2$ 的标准正交基。$X_1^T X_2$ 的 SVD 为 $X_1^T X_2 = P \operatorname{diag}(\cos\theta_1, \ldots, \cos\theta_k) Q^T$。取 $\tilde{X}_1 = X_1 P$，$\tilde{X}_2 = X_2 Q$，则 $\tilde{X}_1^T \tilde{X}_2 = \operatorname{diag}(\cos\theta_i)$。

    构造测地线 $\gamma(t) = \tilde{X}_1 \operatorname{diag}(\cos(t\theta_i)) + U_\perp \operatorname{diag}(\sin(t\theta_i))$，其中 $U_\perp = (\tilde{X}_2 - \tilde{X}_1 \operatorname{diag}(\cos\theta_i)) \operatorname{diag}(\sin\theta_i)^{-1}$。验证 $\gamma(0) = [\tilde{X}_1] = [X_1]$，$\gamma(1) = [X_2]$，且 $\gamma$ 满足测地线方程。长度为 $\int_0^1 \|\dot\gamma\| dt = \|\boldsymbol{\theta}\|_2$。$\blacksquare$

!!! example "例 24.7"
    **$\operatorname{Gr}(1, n)$——实射影空间。**

    $\operatorname{Gr}(1, n)$ 是 $\mathbb{R}^n$ 中所有一维子空间（过原点的直线）的集合，即实射影空间 $\mathbb{RP}^{n-1}$。维数为 $1 \times (n-1) = n - 1$。两条直线之间的测地距离就是它们的夹角 $\theta \in [0, \pi/2]$。

!!! example "例 24.8"
    **计算两个子空间的主角。**

    设 $X_1 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$，$X_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}$（已正交化）。计算

    $$
    X_1^T X_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix},
    $$

    奇异值为 $\sigma_1 = \sigma_2 = 1/\sqrt{2}$，主角 $\theta_1 = \theta_2 = \pi/4$。测地距离 $d = \sqrt{(\pi/4)^2 + (\pi/4)^2} = \frac{\pi}{2\sqrt{2}}$。

---

## 24.6 正定矩阵流形

<div class="context-flow" markdown>

**特殊几何**：$\mathcal{P}(n)$ 上仿射不变度量 $\langle \xi,\eta\rangle_P = \text{tr}(P^{-1}\xi P^{-1}\eta)$ → 测地线 $P^{1/2}e^{tP^{-1/2}\xi P^{-1/2}}P^{1/2}$ → **Hadamard 流形**（非正曲率，唯一 Frechet 均值）

**应用**：协方差估计、扩散张量成像中的几何均值 · $\mathcal{P}(1) = (0,\infty)$ 上距离 $= |\ln(p/q)|$

</div>

对称正定矩阵在统计学、扩散张量成像、协方差估计等领域中扮演重要角色。

!!! definition "定义 24.13 (正定矩阵流形)"
    **对称正定矩阵流形**（manifold of symmetric positive definite matrices）定义为

    $$
    \mathcal{P}(n) = \operatorname{Sym}^+(n) = \{ P \in \mathbb{R}^{n \times n} : P = P^T, \, P \succ 0 \}.
    $$

    $\mathcal{P}(n)$ 是 $\operatorname{Sym}(n)$ 的一个开子集，因此是 $\frac{n(n+1)}{2}$ 维光滑流形，切空间为 $T_P \mathcal{P}(n) = \operatorname{Sym}(n)$。

!!! definition "定义 24.14 (仿射不变 Riemannian 度量)"
    $\mathcal{P}(n)$ 上的 **仿射不变 Riemannian 度量**（affine-invariant Riemannian metric）定义为

    $$
    \langle \xi, \eta \rangle_P = \operatorname{tr}(P^{-1} \xi P^{-1} \eta), \quad \xi, \eta \in T_P \mathcal{P}(n) = \operatorname{Sym}(n).
    $$

    该度量在群作用 $P \mapsto APA^T$（$A \in GL(n)$）下不变。

<div class="context-flow" markdown>

**洞察**：仿射不变度量下 $d(P,Q) = \|\log(P^{-1/2}QP^{-1/2})\|_F$ ——距离由特征值的**对数比**决定，几何均值取代算术均值成为自然的"中心"概念

</div>

!!! theorem "定理 24.11 ($\mathcal{P}(n)$ 上的测地线)"
    在仿射不变度量下，$\mathcal{P}(n)$ 上从 $P$ 出发、初始速度 $\xi \in \operatorname{Sym}(n)$ 的测地线为

    $$
    \gamma(t) = P^{1/2} \exp(t P^{-1/2} \xi P^{-1/2}) P^{1/2}.
    $$

    两点 $P, Q \in \mathcal{P}(n)$ 之间的测地距离为

    $$
    d(P, Q) = \left\| \log(P^{-1/2} Q P^{-1/2}) \right\|_F = \left( \sum_{i=1}^{n} \ln^2 \lambda_i \right)^{1/2},
    $$

    其中 $\lambda_1, \ldots, \lambda_n$ 为 $P^{-1}Q$（或等价地 $P^{-1/2}QP^{-1/2}$）的特征值。

??? proof "证明"
    在 $P = I$ 处，度量简化为 $\langle \xi, \eta \rangle_I = \operatorname{tr}(\xi\eta)$，测地线为 $\gamma(t) = e^{t\xi}$。验证：$\gamma(0) = I$，$\dot\gamma(0) = \xi$，且 $\gamma(t)$ 满足 $\mathcal{P}(n)$ 上的测地线方程 $\ddot\gamma - \dot\gamma \gamma^{-1} \dot\gamma = 0$（这里的运动方程由 Levi-Civita 联络 $\nabla_\xi \eta = D_\xi \eta - \frac{1}{2}(\xi P^{-1}\eta + \eta P^{-1}\xi)$ 推导出）。

    对一般 $P$，利用等距变换 $\Phi_A(P) = APA^T$ 将 $P$ 映射到 $I$（取 $A = P^{-1/2}$），得到一般测地线公式。距离为

    $$
    d(P, Q) = d(I, P^{-1/2}QP^{-1/2}) = \|\log(P^{-1/2}QP^{-1/2})\|_F.
    $$

    $\blacksquare$

!!! theorem "定理 24.12 ($\mathcal{P}(n)$ 上的 Frechet 均值)"
    设 $P_1, \ldots, P_m \in \mathcal{P}(n)$。其关于仿射不变度量的 **Frechet 均值**（Frechet mean）为

    $$
    \bar{P} = \arg\min_{P \in \mathcal{P}(n)} \sum_{i=1}^{m} d(P, P_i)^2.
    $$

    Frechet 均值存在且唯一（因为 $\mathcal{P}(n)$ 在仿射不变度量下是非正曲率完备 Riemannian 流形，即 Hadamard 流形）。

??? proof "证明"
    $\mathcal{P}(n)$ 在仿射不变度量下是完备的（测地线存在于任意时间），且截面曲率 $K \le 0$（这可以通过计算曲率张量验证）。Hadamard 流形上的 Frechet 均值存在且唯一，这是 Cartan-Hadamard 定理和 Karcher 定理的推论。具体地，严格负曲率保证目标函数 $f(P) = \sum_i d(P, P_i)^2$ 是严格凸的（在测地线意义下），因此有唯一极小值点。$\blacksquare$

!!! example "例 24.9"
    **$\mathcal{P}(1)$ 即正实数。**

    $\mathcal{P}(1) = (0, \infty)$，度量为 $ds^2 = dp^2/p^2$（双曲度量）。测地线为 $\gamma(t) = p_0 e^{vt}$，距离为 $d(p, q) = |\ln(p/q)|$。两点 $p, q > 0$ 的 Frechet 均值为几何均值 $\bar{p} = (p \cdot q)^{1/2}$（对两点情形）。

!!! example "例 24.10"
    **计算 $2 \times 2$ 正定矩阵间的距离。**

    设 $P = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$，$Q = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$。则

    $$
    P^{-1}Q = \begin{pmatrix} 1/2 & 0 \\ 0 & 2 \end{pmatrix}, \quad \log(P^{-1}Q) = \begin{pmatrix} -\ln 2 & 0 \\ 0 & \ln 2 \end{pmatrix}.
    $$

    距离 $d(P, Q) = \sqrt{(\ln 2)^2 + (\ln 2)^2} = \ln 2 \cdot \sqrt{2} \approx 0.980$。

---

## 24.7 矩阵 Lie 群与 Lie 代数

<div class="context-flow" markdown>

**统一理论**：BCH 公式 $e^Xe^Y = e^{X+Y+\frac{1}{2}[X,Y]+\cdots}$ 将群乘法编码为 Lie 括号 → Lie 群同态 ↔ Lie 代数同态（单连通时完全对应）

**实用**：$SO(3)$ 中 Rodrigues 公式 · 机器人学/计算机视觉中的小旋转合成

</div>

本节系统讨论矩阵 Lie 群的一般理论。

!!! definition "定义 24.15 (矩阵 Lie 群)"
    **矩阵 Lie 群**（matrix Lie group）是 $GL(n, \mathbb{C})$ 的闭子群 $G$。由 Cartan 闭子群定理，$G$ 自动是光滑流形，因此是 Lie 群。常见的矩阵 Lie 群包括 $GL(n), SL(n), O(n), SO(n), U(n), SU(n), Sp(2n)$ 等。

!!! theorem "定理 24.13 (Baker-Campbell-Hausdorff 公式)"
    设 $X, Y \in \mathfrak{g}$（某 Lie 代数），且 $\|X\|, \|Y\|$ 足够小。则存在 $Z \in \mathfrak{g}$ 使得 $e^X e^Y = e^Z$，且

    $$
    Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}\big([X, [X, Y]] - [Y, [X, Y]]\big) + \cdots
    $$

    该级数称为 **Baker-Campbell-Hausdorff (BCH) 公式**，其中每一项都是 $X, Y$ 的嵌套 Lie 括号。特别地，当 $[X, Y] = 0$ 时，$Z = X + Y$。

??? proof "证明"
    **证明思路。** 在 $e^X e^Y = e^Z$ 两边取对数（在形式幂级数意义下），利用矩阵指数和对数的级数展开

    $$
    e^X = \sum_k \frac{X^k}{k!}, \quad \log(I + W) = \sum_k \frac{(-1)^{k+1}}{k} W^k,
    $$

    将 $Z = \log(e^X e^Y)$ 展开为 $X, Y$ 的幂级数。通过 Dynkin 的显式公式，每一项可以表示为嵌套交换子。

    前几项的具体验证：

    $$
    e^X e^Y = (I + X + \tfrac{X^2}{2} + \cdots)(I + Y + \tfrac{Y^2}{2} + \cdots) = I + (X+Y) + (XY + \tfrac{X^2 + Y^2}{2}) + \cdots
    $$

    $$
    \log(e^X e^Y) = (X+Y) + \frac{XY - YX}{2} + \cdots = X + Y + \frac{1}{2}[X, Y] + \cdots
    $$

    收敛性在 $\|X\| + \|Y\|$ 足够小时由 Lie 代数的完备性保证。$\blacksquare$

!!! theorem "定理 24.14 (Lie 群同态与 Lie 代数同态的对应)"
    设 $G, H$ 为矩阵 Lie 群，$\mathfrak{g}, \mathfrak{h}$ 为其 Lie 代数。若 $\Phi: G \to H$ 为 Lie 群同态（光滑群同态），则 $d\Phi_I: \mathfrak{g} \to \mathfrak{h}$ 为 Lie 代数同态，即保持 Lie 括号：

    $$
    d\Phi_I([X, Y]) = [d\Phi_I(X), d\Phi_I(Y)], \quad \forall X, Y \in \mathfrak{g}.
    $$

    反之，若 $G$ 单连通，则每个 Lie 代数同态 $\phi: \mathfrak{g} \to \mathfrak{h}$ 都可提升为唯一的 Lie 群同态 $\Phi: G \to H$ 使得 $d\Phi_I = \phi$。

??? proof "证明"
    设 $\Phi: G \to H$ 为 Lie 群同态。对 $X, Y \in \mathfrak{g}$，

    $$
    d\Phi_I([X, Y]) = d\Phi_I\!\left(\frac{d}{dt}\Big|_{t=0} e^{tX} Y e^{-tX}\right) = \frac{d}{dt}\Big|_{t=0} \Phi(e^{tX}) d\Phi_I(Y) \Phi(e^{-tX}).
    $$

    由 $\Phi(e^{tX}) = e^{t \, d\Phi_I(X)}$，上式等于 $[d\Phi_I(X), d\Phi_I(Y)]$。

    反方向需要用到 $G$ 的单连通性和指数映射的局部同胚性来逐步扩展 $\phi$ 到整个群。$\blacksquare$

!!! example "例 24.11"
    **$\det: GL(n) \to \mathbb{R}^*$ 的 Lie 代数映射。**

    行列式 $\det: GL(n) \to GL(1) = \mathbb{R}^*$ 是 Lie 群同态。其在 $I$ 处的微分为

    $$
    d(\det)_I(X) = \operatorname{tr}(X).
    $$

    这是因为 $\det(I + tX) = 1 + t\operatorname{tr}(X) + O(t^2)$。对应的 Lie 代数同态为 $\operatorname{tr}: \mathfrak{gl}(n) \to \mathfrak{gl}(1) = \mathbb{R}$，保持 Lie 括号：$\operatorname{tr}([A, B]) = 0 = [\operatorname{tr}(A), \operatorname{tr}(B)]$。

!!! example "例 24.12"
    **BCH 公式的应用：近似乘积。**

    对于小角度旋转 $R_1 = e^{\Omega_1}, R_2 = e^{\Omega_2} \in SO(3)$，BCH 公式给出

    $$
    R_1 R_2 \approx \exp\!\left(\Omega_1 + \Omega_2 + \frac{1}{2}[\Omega_1, \Omega_2]\right).
    $$

    这在机器人学和计算机视觉中用于小旋转的合成。当 $\Omega_1, \Omega_2$ 足够小时，可以进一步近似为 $R_1 R_2 \approx e^{\Omega_1 + \Omega_2}$。

---

## 24.8 矩阵流形上的优化

<div class="context-flow" markdown>

**核心转化**：约束优化($X^TX=I$ 等) → 流形上无约束优化

**Riemannian 梯度** = 欧氏梯度投影到切空间 → **收回映射**(retraction)替代昂贵的测地线步

**链接**：特征值问题 = $\text{St}(k,n)$ 上 $\max \text{tr}(X^TAX)$（Ch25）· PCA = $\text{Gr}(k,n)$ 上优化

</div>

许多实际问题可以建模为矩阵流形上的优化问题。流形优化将约束优化转化为无约束的 Riemannian 优化。

!!! definition "定义 24.16 (Riemannian 梯度)"
    设 $f: \mathcal{M} \to \mathbb{R}$ 为 Riemannian 流形 $(\mathcal{M}, \langle \cdot, \cdot \rangle)$ 上的光滑函数。$f$ 在 $X \in \mathcal{M}$ 处的 **Riemannian 梯度**（Riemannian gradient）$\operatorname{grad} f(X) \in T_X \mathcal{M}$ 定义为满足

    $$
    \langle \operatorname{grad} f(X), \xi \rangle_X = Df(X)[\xi], \quad \forall \xi \in T_X \mathcal{M}
    $$

    的唯一切向量。

!!! definition "定义 24.17 (收回映射)"
    **收回映射**（retraction）$R_X: T_X \mathcal{M} \to \mathcal{M}$ 是指数映射的一阶近似，满足：

    1. $R_X(0) = X$；
    2. $\frac{d}{dt}\Big|_{t=0} R_X(t\xi) = \xi$，$\forall \xi \in T_X \mathcal{M}$。

    收回映射比指数映射计算代价低，在优化算法中常用于替代测地线步。

!!! theorem "定理 24.15 (Riemannian 梯度下降收敛性)"
    设 $f: \mathcal{M} \to \mathbb{R}$ 为紧 Riemannian 流形上的光滑函数，$R$ 为收回映射。Riemannian 梯度下降迭代

    $$
    X_{k+1} = R_{X_k}(-\alpha_k \operatorname{grad} f(X_k))
    $$

    满足：在适当的步长选择（如 Armijo 线搜索）下，$\|\operatorname{grad} f(X_k)\| \to 0$，即迭代点趋向临界点。

??? proof "证明"
    **证明思路。** 利用收回映射的性质和 Taylor 展开：

    $$
    f(R_X(-\alpha \operatorname{grad} f)) = f(X) - \alpha \|\operatorname{grad} f(X)\|^2 + O(\alpha^2).
    $$

    在 Armijo 条件 $f(X_{k+1}) \le f(X_k) - c \alpha_k \|\operatorname{grad} f(X_k)\|^2$（$0 < c < 1$）下，

    $$
    \sum_{k=0}^{\infty} \alpha_k \|\operatorname{grad} f(X_k)\|^2 \le f(X_0) - \inf f < \infty.
    $$

    由步长的下界（回溯线搜索保证 $\alpha_k \ge \alpha_{\min} > 0$），可得 $\|\operatorname{grad} f(X_k)\| \to 0$。紧性保证了下界和 Lipschitz 常数的存在。$\blacksquare$

!!! theorem "定理 24.16 (Stiefel 流形上的 Riemannian 梯度)"
    设 $f: \mathbb{R}^{n \times k} \to \mathbb{R}$ 为光滑函数，限制到 $\operatorname{St}(k, n)$ 上。在典范度量下，Riemannian 梯度为

    $$
    \operatorname{grad} f(X) = \nabla f(X) - X \operatorname{sym}(X^T \nabla f(X)),
    $$

    其中 $\nabla f(X)$ 为 $f$ 在欧氏空间中的梯度，$\operatorname{sym}(A) = (A + A^T)/2$。

??? proof "证明"
    Riemannian 梯度是欧氏梯度到切空间的正交投影：$\operatorname{grad} f(X) = \Pi_X(\nabla f(X))$。由定理 24.8 的投影公式，$\Pi_X(Z) = Z - X\operatorname{sym}(X^TZ)$。代入 $Z = \nabla f(X)$ 即得。$\blacksquare$

!!! example "例 24.13"
    **Stiefel 流形上的特征值问题。**

    求 $A \in \operatorname{Sym}(n)$ 的前 $k$ 个最大特征值对应的特征向量，等价于

    $$
    \max_{X \in \operatorname{St}(k, n)} \operatorname{tr}(X^T A X).
    $$

    欧氏梯度 $\nabla f(X) = 2AX$，Riemannian 梯度为

    $$
    \operatorname{grad} f(X) = 2AX - X \operatorname{sym}(X^T \cdot 2AX) = 2AX - X(X^TAX + (X^TAX)^T)/1 = 2(I - XX^T)AX,
    $$

    其中利用了 $X^TAX$ 的对称性。Riemannian 梯度上升迭代配合收回映射（如 QR 分解收回 $R_X(Z) = \operatorname{qf}(X + Z)$，其中 $\operatorname{qf}$ 表示 QR 分解的 Q 因子）可以求解此问题。

!!! example "例 24.14"
    **Grassmann 流形上的子空间拟合。**

    给定数据矩阵 $Y \in \mathbb{R}^{n \times m}$，寻找最佳 $k$ 维子空间以最大化投影方差：

    $$
    \max_{[X] \in \operatorname{Gr}(k, n)} \|X^T Y\|_F^2 = \max_{[X] \in \operatorname{Gr}(k, n)} \operatorname{tr}(X^T Y Y^T X).
    $$

    这是 PCA 问题的流形优化表述。Riemannian 梯度为 $\operatorname{grad} f([X]) = 2(I - XX^T)YY^TX$，收回映射可以取极分解收回 $R_X(Z) = (X + Z)(I + Z^TZ)^{-1/2}$。

!!! example "例 24.15"
    **正定矩阵流形上的几何均值。**

    给定 $P_1, \ldots, P_m \in \mathcal{P}(n)$，几何均值（Frechet 均值）可以通过 Riemannian 梯度下降求解：

    $$
    \operatorname{grad} f(P) = -\sum_{i=1}^{m} \log(P^{-1} P_i),
    $$

    其中 $f(P) = \frac{1}{2}\sum_i d(P, P_i)^2$，$\log$ 为矩阵对数。迭代

    $$
    P_{k+1} = P_k^{1/2} \exp\!\left(\frac{\alpha}{m} P_k^{-1/2} \sum_{i=1}^{m} \log(P_k^{-1/2} P_i P_k^{-1/2}) P_k^{-1/2}\right) P_k^{1/2}
    $$

    在适当步长 $\alpha$ 下收敛到几何均值。

---

## 本章小结

本章系统介绍了矩阵流形的理论框架和计算方法：

1. **基本概念**：嵌入子流形、切空间、Riemannian 度量为矩阵约束问题提供了几何语言。
2. **一般线性群** $GL(n)$ 是所有矩阵 Lie 群的母群，其 Lie 代数 $\mathfrak{gl}(n)$ 即全矩阵空间。
3. **正交群与酉群** $O(n), U(n)$ 及其特殊子群是最基本的紧 Lie 群，Lie 代数分别为反对称和反 Hermite 矩阵空间。
4. **Stiefel 流形** $\operatorname{St}(k, n)$ 参数化标准正交 $k$-组，是正交约束优化的自然空间。
5. **Grassmann 流形** $\operatorname{Gr}(k, n)$ 参数化 $k$ 维子空间，主角提供了子空间间的自然距离。
6. **正定矩阵流形** $\mathcal{P}(n)$ 在仿射不变度量下是非正曲率空间，具有唯一的 Frechet 均值。
7. **矩阵 Lie 群与 Lie 代数**通过指数映射和 BCH 公式联系。
8. **流形优化**将约束优化转化为 Riemannian 无约束优化，Riemannian 梯度和收回映射是核心工具。
