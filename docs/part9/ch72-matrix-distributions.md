# 第 72 章 矩阵值分布

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 行列式(Ch3) · 矩阵函数(Ch13) · 随机矩阵(Ch23)

**本章脉络**：矩阵变量的概率 $\to$ 矩阵变换的 Jacobi 行列式 $\to$ 矩阵正态分布 $\to$ Wishart 分布 $\to$ 逆 Wishart 分布 $\to$ 矩阵 t 分布 $\to$ 矩阵 Beta 分布 $\to$ 多元分析中的检验统计量

**延伸**：Wishart 分布是多元统计分析的基石；逆 Wishart 分布是贝叶斯统计中协方差矩阵的共轭先验；矩阵值分布在金融风险管理（协方差矩阵建模）和脑成像（扩散张量 MRI）中有直接应用

</div>

经典概率论主要研究标量或向量值的随机变量。然而，在多元统计分析、信号处理和机器学习等领域，我们经常遇到以矩阵为值的随机变量。例如，从多元正态总体中抽取的样本协方差矩阵服从 Wishart 分布；在贝叶斯统计中，协方差矩阵的先验分布通常取逆 Wishart 分布。

矩阵值分布的理论核心是**矩阵变换的 Jacobi 行列式**——它在矩阵空间上的角色类似于一元微积分中的换元法。掌握了 Jacobi 行列式计算，就能推导出所有矩阵值分布的密度函数。

本章将从线性代数的角度系统地发展矩阵值分布理论，强调行列式、正定矩阵和 Kronecker 积在其中的核心作用。

---

## 72.1 矩阵变量的概率基础

<div class="context-flow" markdown>

**核心问题**：如何为矩阵值随机变量定义概率密度？矩阵随机变量的期望和协方差结构如何刻画？

</div>

!!! definition "定义 72.1 (矩阵值随机变量)"
    一个 $p \times n$ **矩阵值随机变量** $X = (X_{ij})$ 是定义在概率空间上的、取值于 $\mathbb{R}^{p \times n}$ 的可测映射。$X$ 的概率密度函数（若存在）是关于 $\mathbb{R}^{p \times n} \cong \mathbb{R}^{pn}$ 上 Lebesgue 测度的 Radon-Nikodym 导数。

    Lebesgue 测度的体积元为：

    $$(dX) = \prod_{i=1}^{p}\prod_{j=1}^{n} dX_{ij}$$

!!! definition "定义 72.2 (矩阵期望)"
    矩阵值随机变量 $X$ 的**期望**定义为逐元素取期望：

    $$E[X] = (E[X_{ij}])_{p \times n}$$

!!! definition "定义 72.3 (矩阵随机变量的协方差结构)"
    $p \times n$ 矩阵值随机变量 $X$ 的协方差结构可以通过 $\text{vec}$ 算子来描述。令 $\mathbf{x} = \text{vec}(X) \in \mathbb{R}^{pn}$（将 $X$ 的各列依次堆叠），则 $\mathbf{x}$ 的协方差矩阵为：

    $$\text{Cov}(\text{vec}(X)) = E[(\text{vec}(X) - \text{vec}(E[X]))(\text{vec}(X) - \text{vec}(E[X]))^T] \in \mathbb{R}^{pn \times pn}$$

    当这个 $pn \times pn$ 协方差矩阵具有 Kronecker 积结构 $\Omega \otimes \Sigma$（$\Omega \in \mathbb{R}^{n \times n}$，$\Sigma \in \mathbb{R}^{p \times p}$）时，我们说 $X$ 的行之间的协方差由 $\Sigma$ 控制，列之间的协方差由 $\Omega$ 控制。

!!! theorem "定理 72.1 (vec 算子与 Kronecker 积)"
    对矩阵乘积 $Y = AXB$（$A \in \mathbb{R}^{m \times p}$，$B \in \mathbb{R}^{n \times q}$）：

    $$\text{vec}(AXB) = (B^T \otimes A)\text{vec}(X)$$

    因此若 $\text{Cov}(\text{vec}(X)) = \Omega \otimes \Sigma$，则：

    $$\text{Cov}(\text{vec}(AXB)) = (B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T) = (B^T\Omega B) \otimes (A\Sigma A^T)$$

??? proof "证明"
    设 $X$ 的第 $j$ 列为 $\mathbf{x}_j$，则 $Y = AXB$ 的第 $k$ 列为 $\mathbf{y}_k = A \sum_j B_{jk}\mathbf{x}_j = \sum_j B_{jk} A\mathbf{x}_j$。

    因此 $\text{vec}(Y) = \text{vec}(AXB)$，其第 $k$ 块（对应 $Y$ 的第 $k$ 列）为 $\sum_j B_{jk} A\mathbf{x}_j$。

    写成矩阵形式：

    $$\text{vec}(Y) = \begin{pmatrix} B_{1,1}A & B_{2,1}A & \cdots \\ B_{1,2}A & B_{2,2}A & \cdots \\ \vdots & & \ddots \end{pmatrix} \text{vec}(X) = (B^T \otimes A)\text{vec}(X)$$

    协方差的变换直接由线性变换的协方差公式得出。

    $\blacksquare$

!!! example "例 72.1"
    设 $X \in \mathbb{R}^{2 \times 3}$，$E[X] = 0$，$\text{Cov}(\text{vec}(X)) = \Omega \otimes \Sigma$，其中：

    $$\Sigma = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}, \quad \Omega = \begin{pmatrix} 1 & 0.5 & 0 \\ 0.5 & 1 & 0.5 \\ 0 & 0.5 & 1 \end{pmatrix}$$

    则 $\text{Cov}(\text{vec}(X))$ 是 $6 \times 6$ 矩阵：

    $$\Omega \otimes \Sigma = \begin{pmatrix} 1\cdot\Sigma & 0.5\cdot\Sigma & 0\cdot\Sigma \\ 0.5\cdot\Sigma & 1\cdot\Sigma & 0.5\cdot\Sigma \\ 0\cdot\Sigma & 0.5\cdot\Sigma & 1\cdot\Sigma \end{pmatrix}$$

    $\Sigma$ 描述同一列中两行的相关性（行协方差），$\Omega$ 描述同一行中不同列的相关性（列协方差）。

---

## 72.2 矩阵变换的 Jacobi 行列式

<div class="context-flow" markdown>

**核心问题**：当对矩阵随机变量进行变换（如求逆、Cholesky 分解）时，密度函数如何变化？

</div>

矩阵变换的 Jacobi 行列式是推导矩阵值分布的核心工具。

!!! theorem "定理 72.2 (线性变换 $X \mapsto AXB$ 的 Jacobi 行列式)"
    设 $A \in \mathbb{R}^{p \times p}$ 和 $B \in \mathbb{R}^{n \times n}$ 为非奇异矩阵，$X \in \mathbb{R}^{p \times n}$，$Y = AXB$。则：

    $$(dY) = |\det A|^n \cdot |\det B|^p \cdot (dX)$$

??? proof "证明"
    $\text{vec}(Y) = (B^T \otimes A)\text{vec}(X)$。线性变换的 Jacobi 行列式等于变换矩阵的行列式的绝对值：

    $$\left|\det(B^T \otimes A)\right| = |\det B^T|^p \cdot |\det A|^n = |\det A|^n \cdot |\det B|^p$$

    其中用到了 Kronecker 积的行列式公式 $\det(B^T \otimes A) = (\det B^T)^p (\det A)^n$（对 $p \times p$ 矩阵 $A$ 和 $n \times n$ 矩阵 $B^T$）。

    $\blacksquare$

!!! theorem "定理 72.3 (矩阵求逆 $X \mapsto X^{-1}$ 的 Jacobi 行列式)"
    设 $X \in \mathbb{R}^{p \times p}$ 为非奇异矩阵，$Y = X^{-1}$。则：

    $$(dY) = |\det X|^{-(p+1)} \cdot (dX)$$

    若限制在对称正定矩阵上（$X = X^T > 0$），则对上三角部分的 $p(p+1)/2$ 个独立变量：

    $$(dY) = |\det X|^{-(p+1)} \cdot (dX)$$

??? proof "证明"
    对 $Y = X^{-1}$，微分得 $dY = -X^{-1}(dX)X^{-1}$。

    将 $dX$ 视为 $X$ 的一个微小扰动矩阵，$dY$ 是对应的 $Y$ 的扰动。映射 $dX \mapsto dY = -X^{-1}(dX)X^{-1}$ 是线性映射。

    由定理 72.2，线性映射 $Z \mapsto -X^{-1}ZX^{-1}$ 的 Jacobi 行列式为：

    $$|\det(-X^{-1})|^p \cdot |\det(X^{-1})|^p = |\det X|^{-p} \cdot |\det X|^{-p} = |\det X|^{-2p}$$

    但这是对 $p^2$ 个独立变量的计算。对一般（非对称）矩阵，实际的 Jacobi 行列式需更仔细的计算。

    通过直接计算（利用 $\partial Y_{kl} / \partial X_{ij} = -(Y)_{ki}(Y)_{jl}$），完整的 Jacobi 矩阵为 $-(Y \otimes Y^T)$，其行列式为 $(-1)^{p^2} (\det Y)^{2p} = (\det X)^{-2p}$。

    对对称正定矩阵，独立变量只有 $p(p+1)/2$ 个（上三角部分），通过类似但更精细的计算，Jacobi 行列式为 $|\det X|^{-(p+1)}$。

    $\blacksquare$

!!! theorem "定理 72.4 (Cholesky 分解 $X = T^TT$ 的 Jacobi 行列式)"
    设 $X$ 为 $p \times p$ 对称正定矩阵，$T$ 为上三角矩阵且对角线元素为正，$X = T^T T$（Cholesky 分解）。则：

    $$(dX) = 2^p \prod_{i=1}^{p} t_{ii}^{p-i+1} \cdot (dT)$$

    其中 $(dX) = \prod_{i \leq j} dX_{ij}$（对称矩阵的上三角部分），$(dT) = \prod_{i \leq j} dT_{ij}$（上三角矩阵）。

??? proof "证明"
    $X = T^T T$，故 $X_{ij} = \sum_{k=1}^{\min(i,j)} T_{ki} T_{kj}$。微分：

    $$dX = (dT)^T T + T^T (dT)$$

    需要计算从 $(dT)$ 到 $(dX)$ 的 Jacobi 行列式。

    对角元素 $X_{ii} = \sum_{k=1}^{i} T_{ki}^2$，故 $dX_{ii}$ 包含 $2T_{ii} dT_{ii}$ 项。

    非对角元素 $X_{ij}$（$i < j$）包含 $T_{ii} dT_{ij}$ 项（以及涉及其他 $dT$ 元素的项）。

    通过归纳法或仔细的行列式计算，Jacobi 行列式为：

    $$2^p \prod_{i=1}^p t_{ii}^{p-i+1}$$

    因子 $2^p$ 来自对角元素的偏导数中的系数 2（$\partial X_{ii}/\partial T_{ii} = 2T_{ii}$），指数 $p - i + 1$ 来自 $T_{ii}$ 在第 $i$ 行及以下的参与度。

    $\blacksquare$

!!! example "例 72.2"
    $p = 2$ 的情形。$X = \begin{pmatrix} X_{11} & X_{12} \\ X_{12} & X_{22} \end{pmatrix} = T^T T$，$T = \begin{pmatrix} t_{11} & t_{12} \\ 0 & t_{22} \end{pmatrix}$。

    $$X_{11} = t_{11}^2, \quad X_{12} = t_{11} t_{12}, \quad X_{22} = t_{12}^2 + t_{22}^2$$

    Jacobi 矩阵：

    $$J = \frac{\partial(X_{11}, X_{12}, X_{22})}{\partial(t_{11}, t_{12}, t_{22})} = \begin{pmatrix} 2t_{11} & 0 & 0 \\ t_{12} & t_{11} & 0 \\ 0 & 2t_{12} & 2t_{22} \end{pmatrix}$$

    $$|\det J| = 2t_{11} \cdot t_{11} \cdot 2t_{22} = 4t_{11}^2 t_{22} = 2^2 \cdot t_{11}^{2} \cdot t_{22}^{1}$$

    与公式 $2^p \prod_{i=1}^p t_{ii}^{p-i+1} = 2^2 t_{11}^2 t_{22}^1$ 一致。

---

## 72.3 矩阵正态分布

<div class="context-flow" markdown>

**核心问题**：如何将多元正态分布推广到矩阵值随机变量？矩阵正态分布的参数有什么结构？

</div>

!!! definition "定义 72.4 (矩阵正态分布)"
    $p \times n$ 矩阵值随机变量 $X$ 服从**矩阵正态分布** $\text{MN}_{p,n}(M, \Sigma, \Omega)$，若 $\text{vec}(X) \sim N_{pn}(\text{vec}(M), \Omega \otimes \Sigma)$，即：

    $$f(X) = \frac{1}{(2\pi)^{pn/2} |\Sigma|^{n/2} |\Omega|^{p/2}} \exp\left(-\frac{1}{2}\text{tr}\left[\Sigma^{-1}(X - M)\Omega^{-1}(X - M)^T\right]\right)$$

    其中：

    - $M \in \mathbb{R}^{p \times n}$：均值矩阵
    - $\Sigma \in \mathbb{R}^{p \times p}$：行协方差矩阵（$\Sigma > 0$）
    - $\Omega \in \mathbb{R}^{n \times n}$：列协方差矩阵（$\Omega > 0$）

!!! theorem "定理 72.5 (矩阵正态分布的性质)"
    设 $X \sim \text{MN}_{p,n}(M, \Sigma, \Omega)$，$A \in \mathbb{R}^{q \times p}$，$B \in \mathbb{R}^{n \times m}$，$C \in \mathbb{R}^{q \times m}$。则：

    1. **仿射变换**：$AXB + C \sim \text{MN}_{q,m}(AMB + C, A\Sigma A^T, B^T\Omega B)$
    2. **行边缘分布**：$X$ 的第 $i$ 行 $\sim N_n(M_{i\cdot}, \Sigma_{ii}\Omega)$
    3. **列边缘分布**：$X$ 的第 $j$ 列 $\sim N_p(M_{\cdot j}, \Omega_{jj}\Sigma)$
    4. **条件分布**：对 $X$ 进行行分块 $X = \begin{pmatrix} X_1 \\ X_2 \end{pmatrix}$，$X_1 \mid X_2$ 也是矩阵正态

??? proof "证明"
    **(1)** $\text{vec}(AXB+C) = (B^T \otimes A)\text{vec}(X) + \text{vec}(C)$。由多元正态的仿射变换性质，

    $$\text{vec}(AXB+C) \sim N\left(\text{vec}(AMB+C), (B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T)\right)$$

    由 Kronecker 积的混合乘积性质：

    $$(B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T) = (B^T\Omega B) \otimes (A\Sigma A^T)$$

    故 $AXB + C \sim \text{MN}_{q,m}(AMB+C, A\Sigma A^T, B^T\Omega B)$。

    **(2)** 取 $A = \mathbf{e}_i^T$（第 $i$ 个标准基向量的转置），$B = I_n$。

    **(3)** 取 $A = I_p$，$B = \mathbf{e}_j$。

    **(4)** 这是多元正态条件分布的矩阵推广，利用 Kronecker 积结构可以直接验证。

    $\blacksquare$

!!! example "例 72.3"
    设 $X \sim \text{MN}_{2,3}(0, I_2, I_3)$，即 $X$ 的 6 个元素是独立标准正态。

    对 $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，$AX \sim \text{MN}_{2,3}(0, AA^T, I_3) = \text{MN}_{2,3}\left(0, \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}, I_3\right)$。

    变换后的行协方差变为 $AA^T = 2I_2$，但列之间仍然独立。

---

## 72.4 Wishart 分布

<div class="context-flow" markdown>

**核心问题**：从多元正态总体中抽取的样本协方差矩阵服从什么分布？

</div>

Wishart 分布是 $\chi^2$ 分布向矩阵的推广，是多元统计分析中最重要的分布。

!!! definition "定义 72.5 (Wishart 分布)"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_n \overset{\text{iid}}{\sim} N_p(\mathbf{0}, \Sigma)$，$X = (\mathbf{x}_1, \ldots, \mathbf{x}_n)^T \in \mathbb{R}^{n \times p}$。则

    $$W = X^T X = \sum_{i=1}^{n} \mathbf{x}_i \mathbf{x}_i^T$$

    服从 **Wishart 分布** $W_p(n, \Sigma)$，参数为自由度 $n \geq p$ 和尺度矩阵 $\Sigma$。

!!! theorem "定理 72.6 (Wishart 分布的密度函数)"
    当 $n \geq p$ 时，$W \sim W_p(n, \Sigma)$ 的密度函数（关于对称正定矩阵空间上的 Lebesgue 测度 $(dW) = \prod_{i \leq j} dW_{ij}$）为：

    $$f(W) = \frac{1}{2^{np/2} \Gamma_p(n/2) |\Sigma|^{n/2}} |W|^{(n-p-1)/2} \exp\left(-\frac{1}{2}\text{tr}(\Sigma^{-1}W)\right)$$

    其中 $W > 0$（正定），$\Gamma_p$ 为**多元 Gamma 函数**：

    $$\Gamma_p(a) = \pi^{p(p-1)/4} \prod_{i=1}^{p} \Gamma\left(a - \frac{i-1}{2}\right)$$

??? proof "证明"
    **步骤 1**：先证 $\Sigma = I_p$ 的情形。$X \in \mathbb{R}^{n \times p}$ 的元素 $X_{ij} \overset{\text{iid}}{\sim} N(0,1)$，联合密度为：

    $$f_X(X) = (2\pi)^{-np/2} \exp\left(-\frac{1}{2}\text{tr}(X^T X)\right)$$

    作变换 $W = X^T X$。利用 $X$ 的 QR 分解 $X = QT$（$Q \in \mathbb{R}^{n \times p}$ 半正交，$T$ 上三角正对角），$W = T^T T$。

    Jacobi 行列式分两步：(a) $X \mapsto (Q, T)$ 的 Jacobi 行列式涉及 Stiefel 流形的体积元；(b) $T \mapsto W = T^T T$ 的 Jacobi 行列式由定理 72.4 给出。

    对 $Q$ 积分（在 Stiefel 流形 $V_{p,n}$ 上积分，得 Stiefel 流形的体积 $\frac{2^p \pi^{np/2}}{\Gamma_p(n/2)}$），得到 $W$ 的密度。

    **步骤 2**：对一般 $\Sigma$，令 $\Sigma = LL^T$（Cholesky 分解），$\mathbf{z}_i = L^{-1}\mathbf{x}_i \sim N(0, I_p)$。$Z = XL^{-T}$，$W = X^TX = L(Z^TZ)L^T$。利用变换公式得一般情形的密度。

    $\blacksquare$

!!! theorem "定理 72.7 (Wishart 分布的性质)"
    设 $W \sim W_p(n, \Sigma)$。

    1. **期望**：$E[W] = n\Sigma$
    2. **可加性**：若 $W_1 \sim W_p(n_1, \Sigma)$，$W_2 \sim W_p(n_2, \Sigma)$ 独立，则 $W_1 + W_2 \sim W_p(n_1 + n_2, \Sigma)$
    3. **线性变换**：$AWA^T \sim W_q(n, A\Sigma A^T)$（$A \in \mathbb{R}^{q \times p}$，$q \leq p$，$A$ 满秩）
    4. **$p = 1$ 的特例**：$W_1(n, \sigma^2) = \sigma^2 \chi^2_n$

??? proof "证明"
    **(1)** $E[W] = E[\sum_i \mathbf{x}_i \mathbf{x}_i^T] = \sum_i E[\mathbf{x}_i \mathbf{x}_i^T] = n\Sigma$。

    **(2)** $W_1 = \sum_{i=1}^{n_1} \mathbf{x}_i \mathbf{x}_i^T$，$W_2 = \sum_{j=1}^{n_2} \mathbf{y}_j \mathbf{y}_j^T$，独立性保证 $W_1 + W_2 = \sum_{k=1}^{n_1+n_2} \mathbf{z}_k \mathbf{z}_k^T$，其中 $\mathbf{z}_k \overset{\text{iid}}{\sim} N(0, \Sigma)$。

    **(3)** $A\mathbf{x}_i \sim N_q(0, A\Sigma A^T)$，故 $AWA^T = \sum_i (A\mathbf{x}_i)(A\mathbf{x}_i)^T \sim W_q(n, A\Sigma A^T)$。

    **(4)** $p = 1$：$W = \sum_i x_i^2$，$x_i \sim N(0, \sigma^2)$，故 $W/\sigma^2 \sim \chi^2_n$。

    $\blacksquare$

!!! definition "定义 72.6 (Bartlett 分解)"
    若 $W \sim W_p(n, I_p)$，则 $W$ 的 Cholesky 因子 $T$（$W = T^TT$）的元素分布为：

    - $t_{ii}^2 \sim \chi^2_{n-i+1}$（$i = 1, \ldots, p$），独立
    - $t_{ij} \sim N(0, 1)$（$i < j$），独立
    - 所有元素相互独立

!!! example "例 72.4"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_{20} \overset{\text{iid}}{\sim} N_3(\boldsymbol{\mu}, \Sigma)$，样本协方差矩阵为：

    $$S = \frac{1}{n-1}\sum_{i=1}^{n}(\mathbf{x}_i - \bar{\mathbf{x}})(\mathbf{x}_i - \bar{\mathbf{x}})^T$$

    则 $(n-1)S \sim W_3(n-1, \Sigma) = W_3(19, \Sigma)$。$E[(n-1)S] = 19\Sigma$，故 $E[S] = \Sigma$（无偏估计）。

---

## 72.5 逆 Wishart 分布

<div class="context-flow" markdown>

**核心问题**：Wishart 矩阵的逆服从什么分布？这个分布为什么在贝叶斯统计中如此重要？

</div>

!!! definition "定义 72.7 (逆 Wishart 分布)"
    若 $W \sim W_p(n, \Sigma)$（$n > p + 1$），则 $W^{-1}$ 服从**逆 Wishart 分布** $\text{IW}_p(n, \Psi)$，其中 $\Psi = \Sigma^{-1}$。密度函数为：

    $$f(X) = \frac{|\Psi|^{n/2}}{2^{np/2}\Gamma_p(n/2)} |X|^{-(n+p+1)/2} \exp\left(-\frac{1}{2}\text{tr}(\Psi X^{-1})\right)$$

    其中 $X > 0$。

??? proof "证明"
    由 $W \sim W_p(n, \Sigma)$ 的密度和变换 $X = W^{-1}$，Jacobi 行列式由定理 72.3 给出：

    $$(dW) = |X|^{-(p+1)} (dX)$$

    代入 $W$ 的密度并令 $W = X^{-1}$：

    $$f(X) = \frac{|X^{-1}|^{(n-p-1)/2}}{2^{np/2}\Gamma_p(n/2)|\Sigma|^{n/2}} \exp\left(-\frac{1}{2}\text{tr}(\Sigma^{-1}X^{-1})\right) \cdot |X|^{-(p+1)}$$

    $$= \frac{|X|^{-(n+p+1)/2}}{2^{np/2}\Gamma_p(n/2)|\Sigma|^{n/2}} \exp\left(-\frac{1}{2}\text{tr}(\Psi X^{-1})\right)$$

    其中用了 $|X^{-1}|^{(n-p-1)/2} \cdot |X|^{-(p+1)} = |X|^{-(n-p-1)/2-(p+1)} = |X|^{-(n+p+1)/2+p/2}$——更仔细地计算：

    $|X^{-1}|^{(n-p-1)/2} = |X|^{-(n-p-1)/2}$，乘以 $|X|^{-(p+1)}$ 得 $|X|^{-(n-p-1)/2-(p+1)} = |X|^{-(n+p+1)/2}$。

    $\blacksquare$

!!! theorem "定理 72.8 (逆 Wishart 的矩)"
    设 $X \sim \text{IW}_p(n, \Psi)$（$n > p + 1$）。

    1. **期望**：$E[X] = \frac{\Psi}{n - p - 1}$
    2. **众数**：$\text{Mode}(X) = \frac{\Psi}{n + p + 1}$

!!! theorem "定理 72.9 (逆 Wishart 作为共轭先验)"
    在贝叶斯框架中，设数据 $\mathbf{x}_1, \ldots, \mathbf{x}_m \overset{\text{iid}}{\sim} N_p(\boldsymbol{\mu}, \Sigma)$（$\boldsymbol{\mu}$ 已知），先验 $\Sigma \sim \text{IW}_p(\nu_0, \Psi_0)$。则后验分布：

    $$\Sigma \mid \mathbf{x}_1, \ldots, \mathbf{x}_m \sim \text{IW}_p\left(\nu_0 + m, \Psi_0 + \sum_{i=1}^{m}(\mathbf{x}_i - \boldsymbol{\mu})(\mathbf{x}_i - \boldsymbol{\mu})^T\right)$$

??? proof "证明"
    似然函数：

    $$L(\Sigma) \propto |\Sigma|^{-m/2} \exp\left(-\frac{1}{2}\text{tr}\left(\Sigma^{-1}\sum_i(\mathbf{x}_i-\boldsymbol{\mu})(\mathbf{x}_i-\boldsymbol{\mu})^T\right)\right)$$

    先验：

    $$\pi(\Sigma) \propto |\Sigma|^{-(\nu_0+p+1)/2} \exp\left(-\frac{1}{2}\text{tr}(\Psi_0\Sigma^{-1})\right)$$

    后验 $\propto$ 似然 $\times$ 先验：

    $$p(\Sigma | \text{data}) \propto |\Sigma|^{-(\nu_0+m+p+1)/2} \exp\left(-\frac{1}{2}\text{tr}\left((\Psi_0 + S_0)\Sigma^{-1}\right)\right)$$

    其中 $S_0 = \sum_i (\mathbf{x}_i - \boldsymbol{\mu})(\mathbf{x}_i - \boldsymbol{\mu})^T$。这是参数为 $(\nu_0 + m, \Psi_0 + S_0)$ 的逆 Wishart 密度。

    $\blacksquare$

!!! example "例 72.5"
    **贝叶斯协方差估计**。设 $p = 2$，先验 $\Sigma \sim \text{IW}_2(5, \Psi_0)$，其中 $\Psi_0 = 3I_2$。观测到 $m = 10$ 个样本，$S_0 = \sum(\mathbf{x}_i - \boldsymbol{\mu})(\mathbf{x}_i - \boldsymbol{\mu})^T = \begin{pmatrix} 8 & 3 \\ 3 & 12 \end{pmatrix}$。

    后验：$\Sigma \mid \text{data} \sim \text{IW}_2\left(15, \begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix}\right)$

    后验期望：$E[\Sigma \mid \text{data}] = \frac{1}{15 - 2 - 1}\begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix} = \frac{1}{12}\begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix} \approx \begin{pmatrix} 0.917 & 0.250 \\ 0.250 & 1.250 \end{pmatrix}$

    与最大似然估计 $S_0/m = \begin{pmatrix} 0.8 & 0.3 \\ 0.3 & 1.2 \end{pmatrix}$ 相比，贝叶斯估计向先验方向收缩。

---

## 72.6 矩阵 Beta 和 F 分布

<div class="context-flow" markdown>

**核心问题**：如何将标量的 Beta 分布和 F 分布推广到矩阵值？这些分布在多元假设检验中起什么作用？

</div>

!!! definition "定义 72.8 (矩阵 Beta 分布)"
    设 $W_1 \sim W_p(n_1, I_p)$ 和 $W_2 \sim W_p(n_2, I_p)$ 独立（$n_1, n_2 \geq p$）。定义 **I 型矩阵 Beta 变量**：

    $$B = (W_1 + W_2)^{-1/2} W_1 (W_1 + W_2)^{-1/2}$$

    或等价地（在适当意义下）：

    $$B = W_1(W_1 + W_2)^{-1}$$

    $B$ 服从参数为 $(n_1/2, n_2/2)$ 的**矩阵 Beta 分布** $\text{Beta}_p(n_1/2, n_2/2)$。$B$ 的特征值 $\beta_1, \ldots, \beta_p \in (0, 1)$。

!!! definition "定义 72.9 (矩阵 F 分布)"
    **矩阵 F 变量**定义为：

    $$F = W_2^{-1/2} W_1 W_2^{-1/2}$$

    或 $F = W_1 W_2^{-1}$（非对称版本）。$F$ 的特征值为 $f_i = \beta_i / (1 - \beta_i)$。

!!! definition "定义 72.10 (Wilks' Lambda)"
    **Wilks' Lambda** 统计量定义为：

    $$\Lambda = \frac{|W_E|}{|W_E + W_H|} = \prod_{i=1}^{p}(1 - \beta_i) = \det(I - B)$$

    其中 $W_E$ 为组内（误差）平方和矩阵，$W_H$ 为组间（假设）平方和矩阵。

    $\Lambda \in (0, 1)$：$\Lambda$ 接近 0 表示组间差异大（拒绝原假设），$\Lambda$ 接近 1 表示组间差异小。

!!! theorem "定理 72.10 (Wilks' Lambda 的分布)"
    在原假设下（各组均值相同），$\Lambda$ 的分布由 $W_E \sim W_p(n_E, \Sigma)$ 和 $W_H \sim W_p(n_H, \Sigma)$ 的独立性决定。$\Lambda$ 的精确分布与矩阵 Beta 分布的行列式有关。

    特殊情形下 $\Lambda$ 有精确的 F 分布变换：

    - $p = 1$：$\frac{1-\Lambda}{\Lambda} \cdot \frac{n_E}{n_H} \sim F(n_H, n_E)$
    - $p = 2$：$\frac{1-\sqrt{\Lambda}}{\sqrt{\Lambda}} \cdot \frac{n_E - 1}{n_H} \sim F(2n_H, 2(n_E-1))$

!!! example "例 72.6"
    **单因素 MANOVA**。三组样本，每组 20 个观测，$p = 2$ 维响应变量。

    $W_E = \begin{pmatrix} 100 & 30 \\ 30 & 80 \end{pmatrix}$，$W_H = \begin{pmatrix} 25 & 10 \\ 10 & 15 \end{pmatrix}$

    $\Lambda = \frac{|W_E|}{|W_E + W_H|} = \frac{100 \times 80 - 30^2}{125 \times 95 - 40^2} = \frac{7100}{10275} \approx 0.691$

    自由度 $n_H = 2$（组数 - 1），$n_E = 57$（总样本数 - 组数）。利用 $p = 2$ 的精确变换检验。

---

## 72.7 多元分析中的检验统计量

<div class="context-flow" markdown>

**核心问题**：除了 Wilks' Lambda，还有哪些基于矩阵特征值的检验统计量？它们之间有什么关系？

</div>

所有 MANOVA 检验统计量都可以表示为 $W_H W_E^{-1}$ 的特征值 $\theta_1 \geq \theta_2 \geq \cdots \geq \theta_s$ 的函数（$s = \min(p, n_H)$）。

!!! definition "定义 72.11 (四大 MANOVA 检验统计量)"
    设 $\theta_1, \ldots, \theta_s$ 为 $W_H W_E^{-1}$ 的非零特征值。

    1. **Wilks' Lambda**：

    $$\Lambda = \prod_{i=1}^{s} \frac{1}{1 + \theta_i} = \frac{|W_E|}{|W_E + W_H|}$$

    2. **Pillai 迹**：

    $$V = \sum_{i=1}^{s} \frac{\theta_i}{1 + \theta_i} = \text{tr}(W_H(W_E + W_H)^{-1})$$

    3. **Hotelling-Lawley 迹**：

    $$U = \sum_{i=1}^{s} \theta_i = \text{tr}(W_H W_E^{-1})$$

    4. **Roy 最大根**：

    $$\Theta = \frac{\theta_1}{1 + \theta_1} = \max_{\mathbf{a}} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T(W_E + W_H)\mathbf{a}}$$

!!! theorem "定理 72.11 (检验统计量的关系)"
    这四个统计量在以下意义下等价：当且仅当 $s = 1$（$p = 1$ 或 $n_H = 1$）时，它们给出相同的检验。当 $s > 1$ 时，不同统计量可能给出不同的检验结论。

    在原假设 $H_0: \boldsymbol{\mu}_1 = \cdots = \boldsymbol{\mu}_g$ 下：

    - **Wilks' Lambda** 对所有备择假设方向有平衡的检验功效
    - **Pillai 迹** 对违背正态性最为稳健
    - **Hotelling-Lawley 迹** 在备择假设集中在一个方向时功效最高
    - **Roy 最大根** 在单方向备择假设下最优，但对多方向备择功效低

!!! theorem "定理 72.12 (Roy 最大根的变分刻画)"
    $$\theta_1 = \max_{\mathbf{a} \neq 0} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T W_E \mathbf{a}}$$

    这是广义特征值问题 $W_H \mathbf{a} = \theta W_E \mathbf{a}$ 的最大特征值。

??? proof "证明"
    由 Rayleigh 商理论（Ch6），对称矩阵 $W_E^{-1/2}W_H W_E^{-1/2}$ 的最大特征值为：

    $$\lambda_{\max}(W_E^{-1/2}W_H W_E^{-1/2}) = \max_{\mathbf{b} \neq 0} \frac{\mathbf{b}^T W_E^{-1/2}W_H W_E^{-1/2}\mathbf{b}}{\mathbf{b}^T\mathbf{b}}$$

    令 $\mathbf{a} = W_E^{-1/2}\mathbf{b}$：

    $$= \max_{\mathbf{a} \neq 0} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T W_E \mathbf{a}}$$

    而 $W_E^{-1/2}W_H W_E^{-1/2}$ 与 $W_H W_E^{-1}$ 有相同的特征值，故 $\theta_1 = \lambda_{\max}(W_H W_E^{-1})$。

    $\blacksquare$

!!! example "例 72.7"
    续例 72.6。

    $$W_H W_E^{-1} = \begin{pmatrix} 25 & 10 \\ 10 & 15 \end{pmatrix}\begin{pmatrix} 100 & 30 \\ 30 & 80 \end{pmatrix}^{-1}$$

    先求 $W_E^{-1}$：$|W_E| = 7100$，$W_E^{-1} = \frac{1}{7100}\begin{pmatrix} 80 & -30 \\ -30 & 100 \end{pmatrix}$

    $$W_H W_E^{-1} = \frac{1}{7100}\begin{pmatrix} 25 & 10 \\ 10 & 15 \end{pmatrix}\begin{pmatrix} 80 & -30 \\ -30 & 100 \end{pmatrix} = \frac{1}{7100}\begin{pmatrix} 1700 & 250 \\ 350 & 1200 \end{pmatrix}$$

    特征值：$\text{tr} = 2900/7100 \approx 0.4085$，$\det = (1700 \times 1200 - 250 \times 350)/7100^2 = (2040000 - 87500)/50410000 \approx 0.03873$。

    $\theta_1 + \theta_2 \approx 0.4085$，$\theta_1 \theta_2 \approx 0.03873$。

    四个统计量：

    - Wilks' $\Lambda \approx 0.691$
    - Pillai $V = \sum \frac{\theta_i}{1+\theta_i} \approx 0.297$
    - Hotelling-Lawley $U = \sum \theta_i \approx 0.409$
    - Roy $\Theta = \theta_1/(1+\theta_1)$

---

## 72.8 应用

<div class="context-flow" markdown>

**核心问题**：矩阵值分布在实际问题中如何应用？

</div>

!!! example "例 72.8 (贝叶斯协方差估计)"
    在投资组合优化中，资产收益率的协方差矩阵 $\Sigma$ 的估计至关重要。经典的样本协方差矩阵在 $p$ 接近 $n$ 时表现不佳（特征值过度分散）。

    贝叶斯方法使用逆 Wishart 先验 $\Sigma \sim \text{IW}_p(\nu_0, \Psi_0)$，其中先验可以编码对协方差结构的先验信念（如 $\Psi_0 = \nu_0 I_p$ 表示先验认为资产间不相关）。

    后验估计 $\hat{\Sigma}_{\text{Bayes}} = E[\Sigma \mid \text{data}]$ 自动实现了向先验的收缩（shrinkage），在高维情形下比样本协方差更稳定。

!!! example "例 72.9 (扩散张量 MRI)"
    在扩散张量 MRI（DTI）中，每个体素（三维像素）对应一个 $3 \times 3$ 对称正定矩阵 $D$——**扩散张量**，描述水分子在该位置的各向异性扩散。

    扩散张量的统计建模自然涉及矩阵值分布。Wishart 分布用于估计扩散张量的不确定性，其特征值（本征值）给出三个主扩散方向的扩散系数：

    - $\lambda_1 \gg \lambda_2 \approx \lambda_3$：各向异性扩散，表示白质纤维束（轴索平行排列）
    - $\lambda_1 \approx \lambda_2 \approx \lambda_3$：各向同性扩散，表示灰质或脑脊液

    **分数各向异性**（Fractional Anisotropy，FA）定义为：

    $$\text{FA} = \sqrt{\frac{3}{2}} \cdot \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$

    其中 $\bar{\lambda} = (\lambda_1 + \lambda_2 + \lambda_3)/3$。FA $\in [0, 1]$，FA 越大表示扩散越各向异性。

!!! example "例 72.10 (Wishart 过程在金融中的应用)"
    协方差矩阵的时变性在金融中非常重要（波动率聚集现象）。**Wishart 过程** $\{W_t\}$ 是 Wishart 分布到连续时间的推广，满足矩阵值随机微分方程：

    $$dW_t = (\nu Q^T Q + M W_t + W_t M^T) dt + \sqrt{W_t} dB_t Q + Q^T dB_t^T \sqrt{W_t}$$

    其中 $B_t$ 是矩阵值 Brown 运动，$M$ 为均值回复矩阵，$Q$ 为波动率矩阵，$\nu$ 为自由度参数。

    Wishart 过程的关键性质：

    - $W_t > 0$（正定性自动保持，只要 $\nu$ 足够大）
    - 边际分布为（非中心）Wishart
    - 可以模拟协方差矩阵的动态变化

!!! theorem "定理 72.13 (矩阵值分布的信息几何)"
    Wishart 分布族 $\{W_p(n, \Sigma) : \Sigma > 0\}$ 构成一个统计流形，其 Fisher 信息矩阵为：

    $$g_{ij,kl} = \frac{n}{2}(\Sigma^{-1})_{ik}(\Sigma^{-1})_{jl} + \frac{n}{2}(\Sigma^{-1})_{il}(\Sigma^{-1})_{jk}$$

    对应的测地距离（Rao 距离）为：

    $$d(\Sigma_1, \Sigma_2) = \sqrt{\frac{n}{2}} \left\|\log(\Sigma_1^{-1/2}\Sigma_2\Sigma_1^{-1/2})\right\|_F$$

    即与正定矩阵流形上的仿射不变 Riemannian 距离成比例。

??? proof "证明"
    Wishart 分布的对数密度为：

    $$\ell(\Sigma) = -\frac{n}{2}\log|\Sigma| - \frac{1}{2}\text{tr}(\Sigma^{-1}W) + \text{const}$$

    Fisher 信息：

    $$I(\Sigma)_{ij,kl} = -E\left[\frac{\partial^2 \ell}{\partial \Sigma_{ij}\partial \Sigma_{kl}}\right]$$

    利用 $\frac{\partial}{\partial \Sigma_{ij}}\log|\Sigma| = (\Sigma^{-1})_{ji}$ 和 $\frac{\partial}{\partial \Sigma_{ij}}\text{tr}(\Sigma^{-1}W) = -(\Sigma^{-1}W\Sigma^{-1})_{ji}$，经过计算得到上述结果。

    测地距离的推导涉及正定矩阵流形的微分几何，与 $\text{GL}(p)/O(p)$ 对称空间的结构有关。

    $\blacksquare$

---

## 本章小结

本章发展了矩阵值分布的系统理论，核心内容包括：

1. **Jacobi 行列式**：矩阵变换（线性变换、求逆、Cholesky 分解）的 Jacobi 行列式是推导所有矩阵值分布密度的关键工具。线性变换 $X \mapsto AXB$ 的 Jacobi 行列式为 $|\det A|^n |\det B|^p$，求逆变换的 Jacobi 行列式为 $|\det X|^{-(p+1)}$。

2. **矩阵正态分布**：$\text{MN}_{p,n}(M, \Sigma, \Omega)$ 通过 Kronecker 积 $\Omega \otimes \Sigma$ 编码行、列之间的协方差结构，是多元正态向矩阵的自然推广。

3. **Wishart 分布**：$W_p(n, \Sigma)$ 是样本协方差矩阵的分布，是 $\chi^2$ 分布的矩阵推广。其密度涉及行列式 $|W|^{(n-p-1)/2}$ 和矩阵迹 $\text{tr}(\Sigma^{-1}W)$。

4. **逆 Wishart 分布**：$\text{IW}_p(n, \Psi)$ 是 Wishart 逆的分布，是贝叶斯统计中协方差矩阵的共轭先验。

5. **MANOVA 检验**：Wilks' Lambda、Pillai 迹、Hotelling-Lawley 迹和 Roy 最大根都可以表示为 $W_H W_E^{-1}$ 的特征值的函数，通过广义特征值问题的 Rayleigh 商刻画。

6. **应用**：从贝叶斯协方差估计到扩散张量 MRI，从金融风险管理到信息几何，矩阵值分布理论将正定矩阵的代数结构与统计推断紧密联系在一起。

贯穿本章的核心线性代数工具是行列式、正定矩阵、Kronecker 积和特征值分析。这些工具在矩阵值概率论中的角色，正如它们在确定性矩阵理论中一样基础且不可或缺。
