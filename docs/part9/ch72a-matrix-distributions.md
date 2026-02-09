# 第 72A 章 矩阵值分布

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 行列式(Ch3) · 矩阵函数(Ch13) · Kronecker 积(Ch10) · 随机矩阵(Ch23) · 多元正态分布(Ch27)

**本章脉络**：矩阵变量的概率基础与 vec-Kronecker 框架 $\to$ 矩阵变换的 Jacobi 行列式(线性变换、求逆、Cholesky、特征值分解) $\to$ 矩阵正态分布 $\to$ Wishart 分布(密度、Bartlett 分解、性质) $\to$ 奇异 Wishart 分布 $\to$ 逆 Wishart 分布 $\to$ 矩阵 $t$ 分布 $\to$ LKJ 分布

**延伸**：Wishart 分布是多元统计分析的基石；逆 Wishart 分布是贝叶斯统计中协方差矩阵的经典共轭先验；矩阵 $t$ 分布联系了贝叶斯推断与稳健统计；LKJ 分布是现代贝叶斯建模中相关矩阵先验的标准选择；矩阵值分布在金融风险管理、脑成像(DTI)和空间统计中有直接应用

</div>

经典概率论主要研究标量或向量值的随机变量。然而，在多元统计分析、信号处理和机器学习等领域，我们经常遇到以矩阵为值的随机变量。例如，从多元正态总体中抽取的样本协方差矩阵服从 Wishart 分布；在贝叶斯统计中，协方差矩阵的先验分布通常取逆 Wishart 分布。

矩阵值分布的理论核心是**矩阵变换的 Jacobi 行列式**——它在矩阵空间上的角色类似于一元微积分中的换元法。掌握了 Jacobi 行列式计算，就能推导出所有矩阵值分布的密度函数。

本章将从线性代数的角度系统地发展矩阵值分布理论，强调行列式、正定矩阵和 Kronecker 积在其中的核心作用。

---

## 72A.1 矩阵变量的概率基础

<div class="context-flow" markdown>

**核心问题**：如何为矩阵值随机变量定义概率密度？矩阵随机变量的期望和协方差结构如何刻画？

</div>

!!! definition "定义 72A.1 (矩阵值随机变量)"
    一个 $p \times n$ **矩阵值随机变量** $X = (X_{ij})$ 是定义在概率空间上的、取值于 $\mathbb{R}^{p \times n}$ 的可测映射。$X$ 的概率密度函数（若存在）是关于 $\mathbb{R}^{p \times n} \cong \mathbb{R}^{pn}$ 上 Lebesgue 测度的 Radon-Nikodym 导数。

    Lebesgue 测度的体积元为：

    $$(dX) = \prod_{i=1}^{p}\prod_{j=1}^{n} dX_{ij}$$

!!! definition "定义 72A.2 (矩阵期望)"
    矩阵值随机变量 $X$ 的**期望**定义为逐元素取期望：

    $$E[X] = (E[X_{ij}])_{p \times n}$$

!!! definition "定义 72A.3 (矩阵随机变量的协方差结构)"
    $p \times n$ 矩阵值随机变量 $X$ 的协方差结构可以通过 $\text{vec}$ 算子来描述。令 $\mathbf{x} = \text{vec}(X) \in \mathbb{R}^{pn}$（将 $X$ 的各列依次堆叠），则 $\mathbf{x}$ 的协方差矩阵为：

    $$\text{Cov}(\text{vec}(X)) = E[(\text{vec}(X) - \text{vec}(E[X]))(\text{vec}(X) - \text{vec}(E[X]))^T] \in \mathbb{R}^{pn \times pn}$$

    当这个 $pn \times pn$ 协方差矩阵具有 Kronecker 积结构 $\Omega \otimes \Sigma$（$\Omega \in \mathbb{R}^{n \times n}$，$\Sigma \in \mathbb{R}^{p \times p}$）时，我们说 $X$ 的行之间的协方差由 $\Sigma$ 控制，列之间的协方差由 $\Omega$ 控制。这种**可分离**（separable）协方差结构大大降低了参数量（从 $O(p^2 n^2)$ 到 $O(p^2 + n^2)$）。

!!! theorem "定理 72A.1 (vec 算子与 Kronecker 积)"
    对矩阵乘积 $Y = AXB$（$A \in \mathbb{R}^{m \times p}$，$B \in \mathbb{R}^{n \times q}$）：

    $$\text{vec}(AXB) = (B^T \otimes A)\text{vec}(X)$$

    因此若 $\text{Cov}(\text{vec}(X)) = \Omega \otimes \Sigma$，则：

    $$\text{Cov}(\text{vec}(AXB)) = (B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T) = (B^T\Omega B) \otimes (A\Sigma A^T)$$

??? proof "证明"
    设 $X$ 的第 $j$ 列为 $\mathbf{x}_j$，则 $Y = AXB$ 的第 $k$ 列为：

    $$\mathbf{y}_k = A \sum_{j=1}^{n} B_{jk}\mathbf{x}_j = \sum_{j=1}^{n} B_{jk} A\mathbf{x}_j$$

    因此 $\text{vec}(Y)$ 的第 $k$ 块（对应 $Y$ 的第 $k$ 列）为 $\sum_j B_{jk} A\mathbf{x}_j$。

    写成矩阵形式：

    $$\text{vec}(Y) = \begin{pmatrix} B_{11}A & B_{21}A & \cdots & B_{n1}A \\ B_{12}A & B_{22}A & \cdots & B_{n2}A \\ \vdots & & & \vdots \\ B_{1q}A & B_{2q}A & \cdots & B_{nq}A \end{pmatrix} \text{vec}(X) = (B^T \otimes A)\text{vec}(X)$$

    协方差变换直接由线性变换的协方差公式得出：

    $$\text{Cov}(\text{vec}(Y)) = (B^T \otimes A) \cdot (\Omega \otimes \Sigma) \cdot (B^T \otimes A)^T$$

    利用 Kronecker 积的混合乘积性质 $(P \otimes Q)(R \otimes S) = (PR) \otimes (QS)$ 和 $(P \otimes Q)^T = P^T \otimes Q^T$：

    $$(B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T) = (B^T \Omega B) \otimes (A \Sigma A^T)$$

    $\blacksquare$

!!! example "例 72A.1"
    设 $X \in \mathbb{R}^{2 \times 3}$，$E[X] = 0$，$\text{Cov}(\text{vec}(X)) = \Omega \otimes \Sigma$，其中：

    $$\Sigma = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}, \quad \Omega = \begin{pmatrix} 1 & 0.5 & 0 \\ 0.5 & 1 & 0.5 \\ 0 & 0.5 & 1 \end{pmatrix}$$

    则 $\text{Cov}(\text{vec}(X))$ 是 $6 \times 6$ 矩阵：

    $$\Omega \otimes \Sigma = \begin{pmatrix} 1\cdot\Sigma & 0.5\cdot\Sigma & 0\cdot\Sigma \\ 0.5\cdot\Sigma & 1\cdot\Sigma & 0.5\cdot\Sigma \\ 0\cdot\Sigma & 0.5\cdot\Sigma & 1\cdot\Sigma \end{pmatrix}$$

    $\Sigma$ 描述同一列中两行的相关性（行协方差），$\Omega$ 描述同一行中不同列的相关性（列协方差）。

---

## 72A.2 矩阵变换的 Jacobi 行列式

<div class="context-flow" markdown>

**核心问题**：当对矩阵随机变量进行变换（如求逆、Cholesky 分解、特征值分解）时，密度函数如何变化？

</div>

矩阵变换的 Jacobi 行列式是推导矩阵值分布的核心工具。

!!! theorem "定理 72A.2 (线性变换 $X \mapsto AXB$ 的 Jacobi 行列式)"
    设 $A \in \mathbb{R}^{p \times p}$ 和 $B \in \mathbb{R}^{n \times n}$ 为非奇异矩阵，$X \in \mathbb{R}^{p \times n}$，$Y = AXB$。则：

    $$(dY) = |\det A|^n \cdot |\det B|^p \cdot (dX)$$

??? proof "证明"
    $\text{vec}(Y) = (B^T \otimes A)\text{vec}(X)$。线性变换 $\mathbf{y} = C\mathbf{x}$ 的 Jacobi 行列式等于变换矩阵的行列式的绝对值 $|\det C|$。

    因此 Jacobi 行列式为：

    $$\left|\det(B^T \otimes A)\right|$$

    由 Kronecker 积的行列式公式：对 $p \times p$ 矩阵 $A$ 和 $n \times n$ 矩阵 $B^T$，

    $$\det(B^T \otimes A) = (\det B^T)^p \cdot (\det A)^n$$

    因此：

    $$|\det(B^T \otimes A)| = |\det A|^n \cdot |\det B|^p$$

    $\blacksquare$

!!! theorem "定理 72A.3 (矩阵求逆 $X \mapsto X^{-1}$ 的 Jacobi 行列式)"
    设 $X \in \mathbb{R}^{p \times p}$ 为非奇异矩阵，$Y = X^{-1}$。

    **一般（非对称）矩阵**：$(dY) = |\det X|^{-2p} \cdot (dX)$

    **对称正定矩阵**（限于上三角 $p(p+1)/2$ 个独立变量）：$(dY) = |\det X|^{-(p+1)} \cdot (dX)$

??? proof "证明"
    对 $Y = X^{-1}$，微分得 $dY = -X^{-1}(dX)X^{-1}$。

    **一般矩阵情形**：映射 $dX \mapsto dY = -X^{-1}(dX)X^{-1}$ 是 $\mathbb{R}^{p \times p}$ 到自身的线性映射。利用 vec 化：

    $$\text{vec}(dY) = -(X^{-T} \otimes X^{-1})\text{vec}(dX)$$

    Jacobi 行列式为：

    $$|\det(-(X^{-T} \otimes X^{-1}))| = |\det X^{-T}|^p \cdot |\det X^{-1}|^p = |\det X|^{-2p}$$

    **对称正定矩阵情形**：$X = X^T > 0$，独立变量为 $X_{ij}$（$i \leq j$）。

    设 $Y = X^{-1}$，对独立变量的偏导数为 $\partial Y_{kl}/\partial X_{ij}$。通过对 $XY = I$ 微分并利用对称性约束，可以证明 Jacobi 矩阵（将 $(X_{ij})_{i \leq j}$ 映射到 $(Y_{kl})_{k \leq l}$）的行列式绝对值为 $|\det X|^{-(p+1)}$。

    直观理解：一般矩阵有 $p^2$ 个独立变量，给出 $|\det X|^{-2p}$；对称矩阵有 $p(p+1)/2$ 个独立变量，"指数减半再加一"给出 $|\det X|^{-(p+1)}$。严格证明需要仔细处理对称约束下的微分形式。

    $\blacksquare$

!!! theorem "定理 72A.4 (Cholesky 分解 $X = T^TT$ 的 Jacobi 行列式)"
    设 $X$ 为 $p \times p$ 对称正定矩阵，$T$ 为上三角矩阵且对角线元素为正，$X = T^T T$（Cholesky 分解）。则：

    $$(dX) = 2^p \prod_{i=1}^{p} t_{ii}^{p-i+1} \cdot (dT)$$

    其中 $(dX) = \prod_{i \leq j} dX_{ij}$，$(dT) = \prod_{i \leq j} dT_{ij}$。

??? proof "证明"
    $X = T^T T$ 逐元素展开：$X_{ij} = \sum_{k=1}^{\min(i,j)} T_{ki} T_{kj}$。微分：

    $$dX = (dT)^T T + T^T (dT)$$

    需要计算从 $p(p+1)/2$ 个变量 $(T_{ij})_{i \leq j}$ 到 $p(p+1)/2$ 个变量 $(X_{ij})_{i \leq j}$ 的 Jacobi 行列式。

    **关键观察**：按照适当的变量排序（先对角元素，再按行排列非对角元素），Jacobi 矩阵是三角形的。

    对角元素：$X_{ii} = \sum_{k=1}^{i} T_{ki}^2$，故 $\partial X_{ii}/\partial T_{ii} = 2T_{ii}$。

    非对角元素（$i < j$）：$X_{ij} = \sum_{k=1}^{i} T_{ki}T_{kj}$，其中涉及 $T_{ij}$ 的项为 $T_{ii}T_{ij}$，故 $\partial X_{ij}/\partial T_{ij} = T_{ii}$。

    采用归纳法：设结论对 $p-1$ 阶成立。将 $X$ 和 $T$ 分块：

    $$X = \begin{pmatrix} X_{11} & \mathbf{x}_{12}^T \\ \mathbf{x}_{12} & X_{22} \end{pmatrix}, \quad T = \begin{pmatrix} t_{11} & \mathbf{t}_{12}^T \\ \mathbf{0} & T_{22} \end{pmatrix}$$

    则 $X_{11} = t_{11}^2$，$\mathbf{x}_{12} = t_{11}\mathbf{t}_{12}$（注意这里用的是 $T^TT$ 约定），$X_{22} = \mathbf{t}_{12}\mathbf{t}_{12}^T + T_{22}^T T_{22}$。

    变换 $t_{11} \mapsto X_{11}$ 的 Jacobi 为 $2t_{11}$；$\mathbf{t}_{12} \mapsto \mathbf{x}_{12}$ 的 Jacobi 为 $t_{11}^{p-1}$（因为 $\mathbf{x}_{12} = t_{11}\mathbf{t}_{12}$，是 $(p-1)$ 维标量乘法）；剩余变量 $T_{22} \mapsto X_{22} - \mathbf{t}_{12}\mathbf{t}_{12}^T$ 由归纳假设给出 $2^{p-1}\prod_{i=2}^{p} t_{ii}^{p-i+1}$。

    合计：$2t_{11} \cdot t_{11}^{p-1} \cdot 2^{p-1}\prod_{i=2}^p t_{ii}^{p-i+1} = 2^p \prod_{i=1}^p t_{ii}^{p-i+1}$。

    $\blacksquare$

!!! theorem "定理 72A.5 (特征值分解的 Jacobi 行列式)"
    设 $X$ 为 $p \times p$ 实对称矩阵，特征值分解 $X = H \Lambda H^T$，其中 $\Lambda = \text{diag}(\lambda_1, \ldots, \lambda_p)$（$\lambda_1 > \lambda_2 > \cdots > \lambda_p$），$H \in O(p)$ 为正交矩阵。则：

    $$(dX) = \prod_{i < j} |\lambda_i - \lambda_j| \cdot (d\Lambda)(dH)$$

    其中 $(dX) = \prod_{i \leq j} dX_{ij}$，$(d\Lambda) = \prod_i d\lambda_i$，$(dH)$ 为 $O(p)/(\mathbb{Z}_2)^p$ 上的 Haar 测度元（因为 $H$ 的列的符号不唯一）。

    因子 $\prod_{i < j}|\lambda_i - \lambda_j|$ 称为 **Vandermonde 行列式**的绝对值，它反映了特征值之间的"排斥"效应。

??? proof "证明"
    $X = H\Lambda H^T$，微分：$dX = (dH)\Lambda H^T + H(d\Lambda)H^T + H\Lambda(dH)^T$。

    由于 $H^TH = I$，$(dH)^TH + H^T(dH) = 0$，令 $\Omega = H^T(dH)$，则 $\Omega$ 是反对称矩阵（$\Omega = -\Omega^T$）。

    左乘 $H^T$、右乘 $H$：$H^T(dX)H = \Omega\Lambda - \Lambda\Omega + d\Lambda$。

    令 $Z = H^T(dX)H$，则 $Z$ 是对称矩阵（因为 $dX$ 对称）。

    对角元素：$Z_{ii} = d\lambda_i$（$\Omega$ 的对角为 0）。

    非对角元素（$i < j$）：$Z_{ij} = \Omega_{ij}(\lambda_j - \lambda_i)$，故 $\Omega_{ij} = Z_{ij}/(\lambda_j - \lambda_i)$。

    变换 $(dX)_{i \leq j} \mapsto (Z_{ij})_{i \leq j}$ 是正交变换（$Z = H^T(dX)H$），Jacobi 行列式为 1。

    变换 $(Z_{ij})_{i \leq j} \mapsto (d\lambda_i, \Omega_{ij}|_{i < j})$ 中，对角元素直接对应 $d\lambda_i$，非对角元素 $Z_{ij} = (\lambda_j - \lambda_i)\Omega_{ij}$。

    因此总 Jacobi 行列式为 $\prod_{i < j}|\lambda_j - \lambda_i|$。

    $\blacksquare$

!!! example "例 72A.2"
    $p = 2$ 的 Cholesky 分解情形。$X = \begin{pmatrix} X_{11} & X_{12} \\ X_{12} & X_{22} \end{pmatrix} = T^T T$，$T = \begin{pmatrix} t_{11} & t_{12} \\ 0 & t_{22} \end{pmatrix}$。

    $$X_{11} = t_{11}^2, \quad X_{12} = t_{11} t_{12}, \quad X_{22} = t_{12}^2 + t_{22}^2$$

    Jacobi 矩阵：

    $$J = \frac{\partial(X_{11}, X_{12}, X_{22})}{\partial(t_{11}, t_{12}, t_{22})} = \begin{pmatrix} 2t_{11} & 0 & 0 \\ t_{12} & t_{11} & 0 \\ 0 & 2t_{12} & 2t_{22} \end{pmatrix}$$

    $$|\det J| = 2t_{11} \cdot t_{11} \cdot 2t_{22} = 4t_{11}^2 t_{22} = 2^2 \cdot t_{11}^{2} \cdot t_{22}^{1}$$

    与公式 $2^p \prod_{i=1}^p t_{ii}^{p-i+1} = 2^2 t_{11}^2 t_{22}^1$ 一致。

!!! example "例 72A.3"
    $p = 2$ 的特征值分解。对 $2 \times 2$ 对称矩阵 $X = \begin{pmatrix} a & b \\ b & c \end{pmatrix}$，特征值 $\lambda_1, \lambda_2$ 和旋转角 $\theta$ 满足 $X = H\Lambda H^T$，其中 $H = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix}$。

    Jacobi 行列式为 $|\lambda_1 - \lambda_2|$。这意味着当两个特征值接近时，$(a, b, c)$ 空间中的"体积"被压缩到 $(\lambda_1, \lambda_2, \theta)$ 空间中——特征值"排斥"效应的几何解释。

---

## 72A.3 矩阵正态分布

<div class="context-flow" markdown>

**核心问题**：如何将多元正态分布推广到矩阵值随机变量？矩阵正态分布的参数有什么结构？

</div>

!!! definition "定义 72A.4 (矩阵正态分布)"
    $p \times n$ 矩阵值随机变量 $X$ 服从**矩阵正态分布** $\text{MN}_{p,n}(M, \Sigma, \Omega)$，若 $\text{vec}(X) \sim N_{pn}(\text{vec}(M), \Omega \otimes \Sigma)$，即：

    $$f(X) = \frac{1}{(2\pi)^{pn/2} |\Sigma|^{n/2} |\Omega|^{p/2}} \exp\left(-\frac{1}{2}\text{tr}\left[\Sigma^{-1}(X - M)\Omega^{-1}(X - M)^T\right]\right)$$

    其中：

    - $M \in \mathbb{R}^{p \times n}$：均值矩阵
    - $\Sigma \in \mathbb{R}^{p \times p}$：行协方差矩阵（$\Sigma > 0$）
    - $\Omega \in \mathbb{R}^{n \times n}$：列协方差矩阵（$\Omega > 0$）

!!! theorem "定理 72A.6 (矩阵正态分布的性质)"
    设 $X \sim \text{MN}_{p,n}(M, \Sigma, \Omega)$，$A \in \mathbb{R}^{q \times p}$，$B \in \mathbb{R}^{n \times m}$，$C \in \mathbb{R}^{q \times m}$。则：

    1. **仿射变换**：$AXB + C \sim \text{MN}_{q,m}(AMB + C, A\Sigma A^T, B^T\Omega B)$
    2. **行边缘分布**：$X$ 的第 $i$ 行 $\sim N_n(M_{i\cdot}, \Sigma_{ii}\Omega)$
    3. **列边缘分布**：$X$ 的第 $j$ 列 $\sim N_p(M_{\cdot j}, \Omega_{jj}\Sigma)$
    4. **条件分布**：对 $X$ 进行行分块 $X = \begin{pmatrix} X_1 \\ X_2 \end{pmatrix}$，$X_1 \mid X_2$ 也是矩阵正态
    5. **最大似然估计**：给定 $X_1, \ldots, X_N \overset{\text{iid}}{\sim} \text{MN}_{p,n}(M, \Sigma, \Omega)$，MLE 满足翻转（flip-flop）方程，需要迭代求解

??? proof "证明"
    **(1)** $\text{vec}(AXB+C) = (B^T \otimes A)\text{vec}(X) + \text{vec}(C)$。由多元正态的仿射变换性质，

    $$\text{vec}(AXB+C) \sim N\left(\text{vec}(AMB+C), (B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T)\right)$$

    由 Kronecker 积的混合乘积性质：

    $$(B^T \otimes A)(\Omega \otimes \Sigma)(B \otimes A^T) = (B^T\Omega B) \otimes (A\Sigma A^T)$$

    故 $AXB + C \sim \text{MN}_{q,m}(AMB+C, A\Sigma A^T, B^T\Omega B)$。

    **(2)** 取 $A = \mathbf{e}_i^T$（第 $i$ 个标准基向量的转置），$B = I_n$。则 $A\Sigma A^T = \Sigma_{ii}$（标量），$B^T\Omega B = \Omega$。

    **(3)** 取 $A = I_p$，$B = \mathbf{e}_j$。则 $A\Sigma A^T = \Sigma$，$B^T\Omega B = \Omega_{jj}$。

    **(4)** 这是多元正态条件分布的矩阵推广。利用 Kronecker 积结构和 Schur 补公式可以直接验证。

    **(5)** MLE 的对数似然为 $\ell(M, \Sigma, \Omega) = -\frac{Nn}{2}\log|\Sigma| - \frac{Np}{2}\log|\Omega| - \frac{1}{2}\sum_{k=1}^N \text{tr}[\Sigma^{-1}(X_k - M)\Omega^{-1}(X_k - M)^T]$。关于 $M$ 最大化得 $\hat{M} = \bar{X}$。关于 $\Sigma$ 和 $\Omega$ 的最大化是耦合的（每个的 MLE 依赖于另一个），需要交替优化（flip-flop 算法）。

    $\blacksquare$

!!! example "例 72A.4"
    设 $X \sim \text{MN}_{2,3}(0, I_2, I_3)$，即 $X$ 的 6 个元素是独立标准正态。

    对 $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，$AX \sim \text{MN}_{2,3}(0, AA^T, I_3) = \text{MN}_{2,3}\left(0, \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}, I_3\right)$。

    变换后的行协方差变为 $AA^T = 2I_2$，但列之间仍然独立。

---

## 72A.4 Wishart 分布

<div class="context-flow" markdown>

**核心问题**：从多元正态总体中抽取的样本协方差矩阵服从什么分布？

</div>

Wishart 分布是 $\chi^2$ 分布向矩阵的推广，是多元统计分析中最重要的分布。

!!! definition "定义 72A.5 (Wishart 分布)"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_n \overset{\text{iid}}{\sim} N_p(\mathbf{0}, \Sigma)$，$X = (\mathbf{x}_1, \ldots, \mathbf{x}_n)^T \in \mathbb{R}^{n \times p}$。则

    $$W = X^T X = \sum_{i=1}^{n} \mathbf{x}_i \mathbf{x}_i^T$$

    服从 **Wishart 分布** $W_p(n, \Sigma)$，参数为自由度 $n \geq p$ 和尺度矩阵 $\Sigma$。

!!! theorem "定理 72A.7 (Wishart 分布的密度函数)"
    当 $n \geq p$ 时，$W \sim W_p(n, \Sigma)$ 的密度函数为：

    $$f(W) = \frac{1}{2^{np/2} \Gamma_p(n/2) |\Sigma|^{n/2}} |W|^{(n-p-1)/2} \exp\left(-\frac{1}{2}\text{tr}(\Sigma^{-1}W)\right)$$

    其中 $W > 0$（正定），$\Gamma_p$ 为**多元 Gamma 函数**：

    $$\Gamma_p(a) = \pi^{p(p-1)/4} \prod_{i=1}^{p} \Gamma\left(a - \frac{i-1}{2}\right)$$

??? proof "证明"
    **步骤 1：$\Sigma = I_p$ 的情形。**

    $X \in \mathbb{R}^{n \times p}$ 的元素 $X_{ij} \overset{\text{iid}}{\sim} N(0,1)$，联合密度为：

    $$f_X(X) = (2\pi)^{-np/2} \exp\left(-\frac{1}{2}\text{tr}(X^T X)\right)$$

    作变换 $W = X^T X$。利用 $X$ 的 QR 分解 $X = QT$（$Q \in V_{p,n}$ 为 $n \times p$ 半正交矩阵——Stiefel 流形的元素，$T$ 为 $p \times p$ 上三角正对角矩阵），$W = T^T T$。

    **步骤 1a**：$X \mapsto (Q, T)$ 的 Jacobi 行列式。$X$ 的密度关于 Stiefel 流形 $V_{p,n}$ 上的 $Q$ 和上三角矩阵 $T$ 分解为：

    $$(dX) = \prod_{i=1}^{p} t_{ii}^{n-i} \cdot (dQ)(dT)$$

    **步骤 1b**：$T \mapsto W = T^T T$ 的 Jacobi 行列式由定理 72A.4 给出：$(dW) = 2^p \prod_{i=1}^p t_{ii}^{p-i+1} (dT)$，即 $(dT) = 2^{-p}\prod_i t_{ii}^{-(p-i+1)} (dW)$。

    代入联合密度，并利用 $\text{tr}(X^TX) = \text{tr}(W)$，$|W| = \prod_i t_{ii}^2$：

    $$f_X dX = (2\pi)^{-np/2} e^{-\text{tr}(W)/2} \prod_i t_{ii}^{n-i} \cdot 2^{-p}\prod_i t_{ii}^{-(p-i+1)} (dQ)(dW)$$

    $$= (2\pi)^{-np/2} 2^{-p} e^{-\text{tr}(W)/2} \prod_i t_{ii}^{n-p-1} (dQ)(dW)$$

    由于 $\prod_i t_{ii}^{n-p-1} = |W|^{(n-p-1)/2}$，对 $Q$ 在 Stiefel 流形上积分得体积 $\text{Vol}(V_{p,n}) = \frac{2^p \pi^{np/2}}{\Gamma_p(n/2)}$。

    最终得 $W \sim W_p(n, I_p)$ 的密度。

    **步骤 2：一般 $\Sigma$。** 令 $\Sigma = LL^T$（Cholesky），$\mathbf{z}_i = L^{-1}\mathbf{x}_i \sim N(0, I_p)$。$W_0 = \sum \mathbf{z}_i\mathbf{z}_i^T \sim W_p(n, I_p)$，$W = LW_0L^T$。由线性变换的 Jacobi 行列式和密度变换公式，得一般情形的密度。

    $\blacksquare$

!!! definition "定义 72A.6 (Bartlett 分解)"
    若 $W \sim W_p(n, I_p)$，则 $W$ 的 Cholesky 因子 $T$（$W = T^TT$）的元素分布为：

    - $t_{ii}^2 \sim \chi^2_{n-i+1}$（$i = 1, \ldots, p$），独立
    - $t_{ij} \sim N(0, 1)$（$i < j$），独立
    - 所有元素相互独立

!!! theorem "定理 72A.8 (Bartlett 分解的证明)"
    Bartlett 分解成立，即 Wishart 矩阵 $W \sim W_p(n, I_p)$ 的 Cholesky 因子 $T$ 的元素独立且分布如上所述。

??? proof "证明"
    由步骤 1 的推导，$W = X^TX$ 且 $X = QT$（QR 分解）。$X$ 的元素独立标准正态。

    $T$ 的元素分布可以从 QR 分解的构造中推导。Gram-Schmidt 正交化过程给出：

    - $T$ 的第一列：$t_{11} = \|\mathbf{x}_{\cdot 1}\|$，其中 $\mathbf{x}_{\cdot 1}$ 是 $X$ 的第一列（$n$ 个独立标准正态），故 $t_{11}^2 = \|\mathbf{x}_{\cdot 1}\|^2 \sim \chi^2_n$。
    - $t_{12} = \mathbf{q}_1^T \mathbf{x}_{\cdot 2}$，其中 $\mathbf{q}_1 = \mathbf{x}_{\cdot 1}/\|\mathbf{x}_{\cdot 1}\|$ 与 $\mathbf{x}_{\cdot 2}$ 独立（给定 $\mathbf{q}_1$ 后，$\mathbf{q}_1^T\mathbf{x}_{\cdot 2} \sim N(0,1)$）。
    - $t_{22}^2 = \|\mathbf{x}_{\cdot 2} - t_{12}\mathbf{q}_1\|^2$，这是 $\mathbf{x}_{\cdot 2}$ 在 $\mathbf{q}_1$ 的正交补上的投影的模方，$\sim \chi^2_{n-1}$。

    一般地，$t_{ii}^2$ 是 $\mathbf{x}_{\cdot i}$ 在前 $i-1$ 个正交基向量的正交补（$n-(i-1)$ 维子空间）上投影的模方，$\sim \chi^2_{n-i+1}$。$t_{ij}$（$i < j$）是 $\mathbf{x}_{\cdot j}$ 在第 $i$ 个正交基向量上的投影系数，$\sim N(0,1)$。

    Gram-Schmidt 过程的逐步构造保证了所有这些随机变量的独立性。

    $\blacksquare$

!!! theorem "定理 72A.9 (Wishart 分布的性质)"
    设 $W \sim W_p(n, \Sigma)$。

    1. **期望**：$E[W] = n\Sigma$
    2. **可加性**：若 $W_1 \sim W_p(n_1, \Sigma)$，$W_2 \sim W_p(n_2, \Sigma)$ 独立，则 $W_1 + W_2 \sim W_p(n_1 + n_2, \Sigma)$
    3. **线性变换**：$AWA^T \sim W_q(n, A\Sigma A^T)$（$A \in \mathbb{R}^{q \times p}$，$q \leq p$，$A$ 满秩）
    4. **$p = 1$ 的特例**：$W_1(n, \sigma^2) = \sigma^2 \chi^2_n$

??? proof "证明"
    **(1)** $E[W] = E[\sum_i \mathbf{x}_i \mathbf{x}_i^T] = \sum_i E[\mathbf{x}_i \mathbf{x}_i^T] = n\Sigma$。

    **(2)** $W_1 = \sum_{i=1}^{n_1} \mathbf{x}_i \mathbf{x}_i^T$，$W_2 = \sum_{j=1}^{n_2} \mathbf{y}_j \mathbf{y}_j^T$，独立性保证 $W_1 + W_2 = \sum_{k=1}^{n_1+n_2} \mathbf{z}_k \mathbf{z}_k^T$，其中 $\mathbf{z}_k \overset{\text{iid}}{\sim} N(0, \Sigma)$。

    **(3)** $A\mathbf{x}_i \sim N_q(0, A\Sigma A^T)$，故 $AWA^T = \sum_i (A\mathbf{x}_i)(A\mathbf{x}_i)^T \sim W_q(n, A\Sigma A^T)$。

    **(4)** $p = 1$：$W = \sum_i x_i^2$，$x_i \sim N(0, \sigma^2)$，故 $W/\sigma^2 = \sum_i (x_i/\sigma)^2 \sim \chi^2_n$。

    $\blacksquare$

!!! example "例 72A.5"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_{20} \overset{\text{iid}}{\sim} N_3(\boldsymbol{\mu}, \Sigma)$，样本协方差矩阵为：

    $$S = \frac{1}{n-1}\sum_{i=1}^{n}(\mathbf{x}_i - \bar{\mathbf{x}})(\mathbf{x}_i - \bar{\mathbf{x}})^T$$

    则 $(n-1)S \sim W_3(n-1, \Sigma) = W_3(19, \Sigma)$。$E[(n-1)S] = 19\Sigma$，故 $E[S] = \Sigma$（无偏估计）。

---

## 72A.5 奇异 Wishart 分布

<div class="context-flow" markdown>

**核心问题**：当样本量 $n$ 小于维度 $p$ 时，样本协方差矩阵是奇异的。这种情况下的分布理论是什么？

</div>

在高维统计中（$p > n$），样本协方差矩阵 $W = X^TX$ 的秩至多为 $n < p$，因此 $W$ 不是正定的，而是半正定的。

!!! definition "定义 72A.7 (奇异 Wishart 分布)"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_n \overset{\text{iid}}{\sim} N_p(\mathbf{0}, \Sigma)$，$n < p$。则 $W = \sum_{i=1}^n \mathbf{x}_i\mathbf{x}_i^T$ 服从**奇异 Wishart 分布**（singular Wishart distribution）$W_p(n, \Sigma)$（$n < p$）。

    $W$ 是 $p \times p$ 对称半正定矩阵，秩为 $\min(n, \text{rank}(\Sigma)) = n$（几乎必然，当 $\Sigma > 0$ 时）。$W$ 没有关于 $\mathbb{R}^{p(p+1)/2}$ 上 Lebesgue 测度的密度。

!!! theorem "定理 72A.10 (奇异 Wishart 的密度——在支撑流形上)"
    当 $n < p$ 时，$W \sim W_p(n, \Sigma)$ 的密度关于秩 $n$ 正半定矩阵流形 $\mathcal{S}_n^+ = \{W \in \mathbb{R}^{p \times p} : W \geq 0, \text{rank}(W) = n\}$ 上的诱导测度存在。

    通过参数化 $W = B^TB$（$B \in \mathbb{R}^{n \times p}$）或特征值分解 $W = H\Lambda H^T$（$\Lambda \in \mathbb{R}^{n \times n}$，$H \in V_{n,p}$），可以写出密度。

    特别地，$W$ 的非零特征值 $\lambda_1 > \cdots > \lambda_n > 0$ 的联合密度为：

    $$f(\lambda_1, \ldots, \lambda_n) \propto \prod_{i=1}^n \lambda_i^{(p-n-1)/2} \cdot \prod_{i<j}(\lambda_i - \lambda_j) \cdot \exp\left(-\frac{1}{2}\text{tr}(\Sigma^{-1}W)\right)$$

    （当 $\Sigma = I_p$ 时的简化形式。）

!!! theorem "定理 72A.11 (奇异 Wishart 的性质)"
    奇异 Wishart 分布 $W \sim W_p(n, \Sigma)$（$n < p$）仍满足：

    1. $E[W] = n\Sigma$
    2. 可加性：若 $W_1 \sim W_p(n_1, \Sigma)$，$W_2 \sim W_p(n_2, \Sigma)$ 独立，则 $W_1 + W_2 \sim W_p(n_1 + n_2, \Sigma)$（即使 $n_1, n_2 < p$，但 $n_1 + n_2$ 可以 $\geq p$）
    3. 线性变换：$AWA^T \sim W_q(n, A\Sigma A^T)$

??? proof "证明"
    **(1)** 和 **(3)** 的证明与非奇异情形完全相同（不依赖正定性）。

    **(2)** $W_1 + W_2 = \sum_{k=1}^{n_1+n_2} \mathbf{z}_k\mathbf{z}_k^T$，其中 $\mathbf{z}_k \overset{\text{iid}}{\sim} N_p(0, \Sigma)$。当 $n_1 + n_2 \geq p$ 时，$W_1 + W_2$ 几乎必然正定，服从（非奇异）$W_p(n_1+n_2, \Sigma)$。

    这个性质很有用：即使单个样本协方差矩阵是奇异的，合并多个独立数据集的协方差矩阵可以得到非奇异的 Wishart。

    $\blacksquare$

!!! example "例 72A.6"
    **高维基因表达数据**。$p = 1000$ 个基因，$n = 50$ 个样本。样本协方差矩阵 $W = X^TX$ 的秩为 50，有 950 个零特征值。

    $W$ 不可逆，无法直接估计精度矩阵 $\Sigma^{-1}$。这是高维统计中的核心挑战，解决方法包括：

    - 正则化（如 Graphical LASSO，第 70B 章）
    - 收缩估计（Ledoit-Wolf）
    - 贝叶斯方法（逆 Wishart 先验提供正则化）

---

## 72A.6 逆 Wishart 分布

<div class="context-flow" markdown>

**核心问题**：Wishart 矩阵的逆服从什么分布？这个分布为什么在贝叶斯统计中如此重要？

</div>

!!! definition "定义 72A.8 (逆 Wishart 分布)"
    若 $W \sim W_p(n, \Sigma)$（$n > p + 1$），则 $W^{-1}$ 服从**逆 Wishart 分布** $\text{IW}_p(n, \Psi)$，其中 $\Psi = \Sigma^{-1}$。密度函数为：

    $$f(X) = \frac{|\Psi|^{n/2}}{2^{np/2}\Gamma_p(n/2)} |X|^{-(n+p+1)/2} \exp\left(-\frac{1}{2}\text{tr}(\Psi X^{-1})\right)$$

    其中 $X > 0$。

??? proof "证明"
    由 $W \sim W_p(n, \Sigma)$ 的密度和变换 $X = W^{-1}$，Jacobi 行列式由定理 72A.3 给出：

    $$(dW) = |X|^{-(p+1)} (dX)$$

    代入 $W$ 的密度并令 $W = X^{-1}$：

    $$f_W(X^{-1}) \cdot |X|^{-(p+1)} = \frac{|X^{-1}|^{(n-p-1)/2}}{2^{np/2}\Gamma_p(n/2)|\Sigma|^{n/2}} \exp\left(-\frac{1}{2}\text{tr}(\Sigma^{-1}X^{-1})\right) \cdot |X|^{-(p+1)}$$

    由 $|X^{-1}|^{(n-p-1)/2} = |X|^{-(n-p-1)/2}$：

    $$f_X(X) = \frac{|X|^{-(n-p-1)/2-(p+1)}}{2^{np/2}\Gamma_p(n/2)|\Sigma|^{n/2}} \exp\left(-\frac{1}{2}\text{tr}(\Psi X^{-1})\right)$$

    指数为 $-(n-p-1)/2 - (p+1) = -(n+p+1)/2$，并将 $|\Sigma|^{n/2} = |\Psi|^{-n/2}$ 代入得结果。

    $\blacksquare$

!!! theorem "定理 72A.12 (逆 Wishart 的矩)"
    设 $X \sim \text{IW}_p(n, \Psi)$（$n > p + 1$）。

    1. **期望**：$E[X] = \frac{\Psi}{n - p - 1}$
    2. **众数**：$\text{Mode}(X) = \frac{\Psi}{n + p + 1}$

!!! theorem "定理 72A.13 (逆 Wishart 作为共轭先验)"
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

!!! example "例 72A.7"
    **贝叶斯协方差估计**。设 $p = 2$，先验 $\Sigma \sim \text{IW}_2(5, \Psi_0)$，其中 $\Psi_0 = 3I_2$。观测到 $m = 10$ 个样本，$S_0 = \sum(\mathbf{x}_i - \boldsymbol{\mu})(\mathbf{x}_i - \boldsymbol{\mu})^T = \begin{pmatrix} 8 & 3 \\ 3 & 12 \end{pmatrix}$。

    后验：$\Sigma \mid \text{data} \sim \text{IW}_2\left(15, \begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix}\right)$

    后验期望：$E[\Sigma \mid \text{data}] = \frac{1}{15 - 2 - 1}\begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix} = \frac{1}{12}\begin{pmatrix} 11 & 3 \\ 3 & 15 \end{pmatrix} \approx \begin{pmatrix} 0.917 & 0.250 \\ 0.250 & 1.250 \end{pmatrix}$

    与最大似然估计 $S_0/m = \begin{pmatrix} 0.8 & 0.3 \\ 0.3 & 1.2 \end{pmatrix}$ 相比，贝叶斯估计向先验方向收缩（对角线增大，非对角线缩小）。

---

## 72A.7 矩阵 $t$ 分布

<div class="context-flow" markdown>

**核心问题**：如何将标量 $t$ 分布推广到矩阵值？矩阵 $t$ 分布与贝叶斯推断有什么联系？

</div>

矩阵 $t$ 分布是矩阵正态分布的"重尾"推广，在贝叶斯统计和稳健估计中有重要应用。

!!! definition "定义 72A.9 (矩阵 $t$ 分布)"
    $p \times n$ 矩阵值随机变量 $X$ 服从**矩阵 $t$ 分布** $\text{MT}_{p,n}(\nu, M, \Sigma, \Omega)$，若其密度为：

    $$f(X) = \frac{\Gamma_p((\nu + n + p - 1)/2)}{(\pi)^{np/2}\Gamma_p((\nu + p - 1)/2)} \frac{|\Sigma|^{-n/2}|\Omega|^{-p/2}}{|I_n + \Omega^{-1}(X-M)^T\Sigma^{-1}(X-M)|^{(\nu+n+p-1)/2}}$$

    其中 $\nu > 0$ 为自由度参数，$M \in \mathbb{R}^{p \times n}$ 为位置参数，$\Sigma \in \mathbb{R}^{p \times p}$ 和 $\Omega \in \mathbb{R}^{n \times n}$ 为正定尺度矩阵。

!!! theorem "定理 72A.14 (矩阵 $t$ 分布的层级表示)"
    矩阵 $t$ 分布可以表示为矩阵正态分布关于逆 Wishart 分布的混合：

    若 $\Sigma_0 \sim \text{IW}_p(\nu, \Psi)$，$X \mid \Sigma_0 \sim \text{MN}_{p,n}(M, \Sigma_0, \Omega)$，

    则 $X$ 的边缘分布为矩阵 $t$ 分布 $\text{MT}_{p,n}(\nu, M, \Psi/(\nu - p - 1), \Omega)$（在适当的参数化下）。

??? proof "证明"
    边缘密度为：

    $$f(X) = \int_{\Sigma_0 > 0} f(X \mid \Sigma_0) \cdot \pi(\Sigma_0) \, d\Sigma_0$$

    $$= \int_{\Sigma_0 > 0} \frac{e^{-\frac{1}{2}\text{tr}[\Sigma_0^{-1}(X-M)\Omega^{-1}(X-M)^T]}}{(2\pi)^{np/2}|\Sigma_0|^{n/2}|\Omega|^{p/2}} \cdot \frac{|\Psi|^{\nu/2}|\Sigma_0|^{-(\nu+p+1)/2}e^{-\frac{1}{2}\text{tr}(\Psi\Sigma_0^{-1})}}{2^{\nu p/2}\Gamma_p(\nu/2)} \, d\Sigma_0$$

    合并 $\Sigma_0$ 的幂次：$|\Sigma_0|^{-n/2 - (\nu+p+1)/2} = |\Sigma_0|^{-(\nu+n+p+1)/2}$。

    合并迹中的项：$\text{tr}[\Sigma_0^{-1}(\Psi + (X-M)\Omega^{-1}(X-M)^T)]$。

    令 $\Psi' = \Psi + (X-M)\Omega^{-1}(X-M)^T$，被积函数正比于 $\text{IW}_p(\nu + n, \Psi')$ 的密度。对 $\Sigma_0$ 积分得归一化常数：

    $$\int \propto \frac{\Gamma_p((\nu+n)/2)}{|\Psi'|^{(\nu+n)/2}}$$

    代入整理后得到矩阵 $t$ 分布的密度。

    $\blacksquare$

!!! theorem "定理 72A.15 (矩阵 $t$ 分布的性质)"
    设 $X \sim \text{MT}_{p,n}(\nu, M, \Sigma, \Omega)$。

    1. **期望**：$E[X] = M$（当 $\nu > 2$ 时存在）
    2. **协方差**：$\text{Cov}(\text{vec}(X)) = \frac{\nu}{\nu - 2}(\Omega \otimes \Sigma)$（当 $\nu > 2$ 时存在）
    3. **极限行为**：当 $\nu \to \infty$ 时，$\text{MT}_{p,n}(\nu, M, \Sigma, \Omega) \to \text{MN}_{p,n}(M, \Sigma, \Omega)$
    4. **仿射变换**：$AXB + C \sim \text{MT}_{q,m}(\nu, AMB+C, A\Sigma A^T, B^T\Omega B)$
    5. **边缘分布**：$X$ 的任意子矩阵也服从矩阵 $t$ 分布

!!! example "例 72A.8"
    **贝叶斯多元回归的预测分布**。设 $Y = XB + E$，$E$ 的行独立同分布 $N_p(0, \Sigma)$。若使用共轭先验 $B \sim \text{MN}(B_0, \Sigma, \Lambda_0^{-1})$ 和 $\Sigma \sim \text{IW}_p(\nu_0, \Psi_0)$，则 $B$ 的后验边缘分布为矩阵 $t$ 分布。

    对新观测 $\mathbf{x}_{\text{new}}$ 的预测分布 $Y_{\text{new}} \mid \text{data}$ 也是（多元）$t$ 分布，比正态分布有更重的尾部，自然地反映了参数不确定性。

---

## 72A.8 LKJ 分布

<div class="context-flow" markdown>

**核心问题**：如何为**相关矩阵**（而非协方差矩阵）指定先验分布？

</div>

逆 Wishart 分布是**协方差矩阵**的先验，但它对相关矩阵和方差参数施加了耦合的先验结构，不够灵活。LKJ 分布直接定义在**相关矩阵**空间上，是现代贝叶斯建模的标准选择。

!!! definition "定义 72A.10 (LKJ 分布)"
    **LKJ 分布**（Lewandowski-Kurowicka-Joe, 2009）是定义在 $p \times p$ 相关矩阵空间上的分布。相关矩阵 $R$ 是对称正定矩阵且对角线全为 1（$R_{ii} = 1$），非对角元素 $R_{ij} \in (-1, 1)$。

    LKJ 分布的密度（关于相关矩阵空间上的自然测度）为：

    $$f(R \mid \eta) \propto |\det R|^{\eta - 1}$$

    其中 $\eta > 0$ 为形状参数。

    - $\eta = 1$：均匀分布——所有有效相关矩阵等概率
    - $\eta > 1$：偏好接近单位阵的相关矩阵（弱相关）
    - $0 < \eta < 1$：偏好强相关的矩阵

!!! theorem "定理 72A.16 (LKJ 分布的性质)"
    设 $R \sim \text{LKJ}(\eta)$。

    1. **边缘分布**：$R$ 的任意 $(i,j)$ 非对角元素的边缘分布以 0 为中心，方差随 $\eta$ 增大而减小。具体地，$R_{ij}$ 的边缘分布近似为 Beta 分布的变换。

    2. **与 Cholesky 因子的关系**：$R$ 的 Cholesky 因子 $L$（$R = LL^T$）的元素有更简单的独立先验。LKJ 分布可以通过对 $L$ 的非对角元素施加 Beta 先验来诱导。

    3. **归一化常数**：

    $$\int |\det R|^{\eta-1} (dR) = \prod_{k=2}^{p} B\left(\eta + \frac{p-k}{2}, \eta + \frac{p-k}{2}\right)^{-1} \cdot C_p$$

    其中 $B$ 为 Beta 函数，$C_p$ 为相关矩阵空间的体积。

??? proof "证明"
    **(1)** 对 $p = 2$ 的情形，$R = \begin{pmatrix} 1 & r \\ r & 1 \end{pmatrix}$，$\det R = 1 - r^2$。

    $f(r \mid \eta) \propto (1 - r^2)^{\eta - 1}$，$r \in (-1, 1)$。

    这正是参数为 $(\eta, \eta)$ 的 Beta 分布在 $(-1, 1)$ 上的变换：令 $u = (1+r)/2 \in (0,1)$，则 $u \sim \text{Beta}(\eta, \eta)$。

    当 $\eta = 1$ 时，$r$ 均匀分布在 $(-1, 1)$ 上。当 $\eta$ 增大时，$r$ 集中在 0 附近。

    **(2)** Lewandowski, Kurowicka 和 Joe 的原始论文通过 Cholesky 因子的参数化推导了 LKJ 分布。令 $R = LL^T$，$L$ 为下三角且对角线为正。定义 $L_{ij} = z_{ij} \prod_{k=1}^{j-1}\sqrt{1-z_{ik}^2}$（$i > j$），其中 $z_{ij} \in (-1,1)$。在适当的 Beta 先验下对 $z_{ij}$ 赋权，可以诱导 $|\det R|^{\eta-1}$ 的密度。

    $\blacksquare$

!!! example "例 72A.9"
    **贝叶斯层级模型中的协方差先验**。在实际建模中，常将协方差矩阵分解为：

    $$\Sigma = \text{diag}(\boldsymbol{\sigma}) \cdot R \cdot \text{diag}(\boldsymbol{\sigma})$$

    其中 $\boldsymbol{\sigma} = (\sigma_1, \ldots, \sigma_p)^T$ 为标准差向量，$R$ 为相关矩阵。

    先验设置：$\sigma_i \sim \text{Half-Cauchy}(\tau)$（独立），$R \sim \text{LKJ}(\eta)$。

    与逆 Wishart 先验相比，这种分解允许对方差和相关结构施加**独立的**先验信息，更加灵活。特别是在 Stan、PyMC 等概率编程语言中，LKJ 先验已成为默认选择。

---

## 习题

!!! exercise "习题 72A.1"
    证明 Kronecker 积的行列式公式：$\det(A \otimes B) = (\det A)^q (\det B)^p$，其中 $A \in \mathbb{R}^{p \times p}$，$B \in \mathbb{R}^{q \times q}$。

    提示：利用 $A \otimes B = (A \otimes I_q)(I_p \otimes B)$。

!!! exercise "习题 72A.2"
    对 $p = 3$ 的对称正定矩阵，Cholesky 分解 $X = T^TT$ 有 6 个独立变量。写出 Jacobi 矩阵（$6 \times 6$）并直接验证其行列式等于 $2^3 t_{11}^3 t_{22}^2 t_{33}^1 = 8 t_{11}^3 t_{22}^2 t_{33}$。

!!! exercise "习题 72A.3"
    证明矩阵正态分布 $X \sim \text{MN}_{p,n}(0, \Sigma, \Omega)$ 的矩母函数为：

    $$E[e^{\text{tr}(T^TX)}] = \exp\left(\frac{1}{2}\text{tr}(\Sigma T \Omega T^T)\right)$$

    其中 $T \in \mathbb{R}^{p \times n}$ 为参数矩阵。

!!! exercise "习题 72A.4"
    利用 Bartlett 分解证明：若 $W \sim W_p(n, I_p)$，则 $\det W = \prod_{i=1}^p t_{ii}^2$，其中 $t_{ii}^2 \overset{\text{indep}}{\sim} \chi^2_{n-i+1}$。由此推导 $\log\det W$ 的均值和方差。

!!! exercise "习题 72A.5"
    证明逆 Wishart 分布的期望公式 $E[X] = \Psi/(n-p-1)$。

    提示：利用 $E[W^{-1}]$ 其中 $W \sim W_p(n, \Psi^{-1})$，对 Wishart 密度作适当的变量替换。

!!! exercise "习题 72A.6"
    设 $W_1 \sim W_p(n_1, \Sigma)$ 和 $W_2 \sim W_p(n_2, \Sigma)$ 独立。证明 $W_1 + W_2 \sim W_p(n_1+n_2, \Sigma)$ 而不需要使用矩母函数，只需回到定义（即 $W_i$ 表示为正态向量外积之和）。

!!! exercise "习题 72A.7"
    对 $p = 2$ 的特征值分解 Jacobi 行列式进行直接验证。设 $X = \begin{pmatrix} a & b \\ b & c \end{pmatrix}$，$\lambda_1, \lambda_2$ 为特征值，$\theta$ 为旋转角。

    (a) 写出 $(a, b, c) \leftrightarrow (\lambda_1, \lambda_2, \theta)$ 的显式映射关系。

    (b) 计算 $3 \times 3$ Jacobi 矩阵并验证其行列式绝对值为 $|\lambda_1 - \lambda_2|$。

!!! exercise "习题 72A.8"
    **矩阵 $t$ 分布的退化**。证明当 $\nu \to \infty$ 时，$\text{MT}_{p,n}(\nu, M, \Sigma, \Omega)$ 的密度逐点收敛于 $\text{MN}_{p,n}(M, \Sigma, \Omega)$ 的密度。

    提示：利用 $(1 + x/\nu)^{-\nu/2} \to e^{-x/2}$（$\nu \to \infty$）的推广。

!!! exercise "习题 72A.9"
    对 LKJ($\eta$) 分布，当 $p = 3$ 时相关矩阵 $R = \begin{pmatrix} 1 & r_{12} & r_{13} \\ r_{12} & 1 & r_{23} \\ r_{13} & r_{23} & 1 \end{pmatrix}$：

    (a) 写出 $\det R$ 关于 $r_{12}, r_{13}, r_{23}$ 的表达式。

    (b) 画出 $\eta = 1, 2, 5$ 时 $r_{12}$ 的边缘密度的定性形状（假设其他相关系数被积分掉）。

    (c) 说明为什么 $\eta > 1$ 在实践中是合理的默认选择。

!!! exercise "习题 72A.10"
    **奇异 Wishart 的模拟**。

    (a) 写出 $p = 5, n = 3, \Sigma = I_5$ 时奇异 Wishart 矩阵 $W$ 的生成步骤。

    (b) 验证 $W$ 的秩为 3（几乎必然）。

    (c) 计算 $W$ 的非零特征值，并与 Marchenko-Pastur 律（第 72B 章）的预测进行比较。
