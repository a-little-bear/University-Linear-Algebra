# 第 23 章 随机矩阵初步

<div class="context-flow" markdown>

**前置**：特征值/SVD(Ch6-8) · 正定矩阵(Ch7)

**脉络**：$n \to \infty$ 时特征值的统计行为——半圆律(Wigner) → MP 律(Wishart/样本协方差) → Tracy-Widom(边缘涨落) → 普适性

**延伸**：随机矩阵在无线通信（MIMO 信道容量）、量子混沌、数论（Riemann zeta 函数零点的 Montgomery-Odlyzko 统计）、金融数学（资产相关性建模）中有深刻应用

</div>

随机矩阵理论（Random Matrix Theory, RMT）研究矩阵元素为随机变量时，矩阵特征值和特征向量的统计性质。该理论起源于 20 世纪 50 年代 Wigner 对原子核能级统计的研究，此后在数学物理、数论、无线通信、高维统计等众多领域产生了深远影响。本章将系统介绍随机矩阵的基本概念、核心极限定理以及若干前沿应用。

---

## 23.1 随机矩阵基本概念

<div class="context-flow" markdown>

**三大系综**：GOE($\beta=1$, 实对称) / GUE($\beta=2$, Hermite) / GSE($\beta=4$, 四元数) → 对称性决定 $\beta$ → 联合密度 $\propto \prod_{i<j}|\lambda_i - \lambda_j|^\beta \cdot e^{-\sum \lambda_i^2}$

**Wishart 矩阵** $W = \frac{1}{n}X^TX$：样本协方差矩阵的原型 → 链接 Ch25 PCA

</div>

随机矩阵是指其元素由随机变量构成的矩阵。我们关注的核心问题是：当矩阵维数 $n \to \infty$ 时，特征值的经验分布呈现何种确定性极限？

!!! definition "定义 23.1 (随机矩阵)"
    设 $(\Omega, \mathcal{F}, P)$ 为概率空间。一个 **随机矩阵**（random matrix）是一个可测映射 $M: \Omega \to \mathbb{K}^{n \times n}$，其中 $\mathbb{K} = \mathbb{R}$ 或 $\mathbb{C}$，即 $M$ 的每个元素 $M_{ij}(\omega)$ 都是定义在该概率空间上的随机变量。

!!! definition "定义 23.2 (高斯正交系综 GOE)"
    **高斯正交系综**（Gaussian Orthogonal Ensemble, GOE）是 $n \times n$ 实对称随机矩阵 $M$ 的概率分布，其密度函数为

    $$
    f(M) = C_n \exp\!\left( -\frac{n}{4} \operatorname{tr}(M^2) \right),
    $$

    其中 $C_n$ 为归一化常数。等价地，$M$ 的上三角元素独立，对角元素 $M_{ii} \sim N(0, 2/n)$，非对角元素 $M_{ij} \sim N(0, 1/n)$（$i < j$），且 $M_{ji} = M_{ij}$。GOE 的分布在正交共轭 $M \mapsto O^T M O$（$O \in O(n)$）下不变。

!!! definition "定义 23.3 (高斯酉系综 GUE)"
    **高斯酉系综**（Gaussian Unitary Ensemble, GUE）是 $n \times n$ Hermite 随机矩阵 $M$ 的概率分布，其密度函数为

    $$
    f(M) = \widetilde{C}_n \exp\!\left( -\frac{n}{2} \operatorname{tr}(M^2) \right).
    $$

    等价地，对角元素 $M_{ii} \sim N(0, 1/n)$ 为实随机变量；非对角元素（$i < j$）的实部和虚部独立且均服从 $N(0, 1/(2n))$，并令 $M_{ji} = \overline{M_{ij}}$。GUE 的分布在酉共轭 $M \mapsto U^* M U$（$U \in U(n)$）下不变。

!!! definition "定义 23.4 (高斯辛系综 GSE)"
    **高斯辛系综**（Gaussian Symplectic Ensemble, GSE）是 $2n \times 2n$ 自对偶四元数 Hermite 矩阵的概率分布。其分布在辛共轭 $M \mapsto S^* M S$（$S \in Sp(2n)$）下不变，密度函数为

    $$
    f(M) = \hat{C}_n \exp\!\left( -n \operatorname{tr}(M^2) \right).
    $$

    GOE、GUE、GSE 分别对应 Dyson 指标 $\beta = 1, 2, 4$。

!!! definition "定义 23.5 (Wishart 矩阵)"
    设 $X$ 为 $n \times p$ 矩阵，其行向量独立同分布于 $N(\mathbf{0}, \Sigma)$，则 **Wishart 矩阵**（Wishart matrix）定义为

    $$
    W = \frac{1}{n} X^T X.
    $$

    当 $\Sigma = I_p$ 时，$W$ 称为白 Wishart 矩阵。Wishart 分布记为 $W \sim \mathcal{W}_p(\Sigma, n)$。

!!! definition "定义 23.6 (经验谱分布)"
    设 $M$ 为 $n \times n$ Hermite 矩阵，特征值为 $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$。**经验谱分布**（Empirical Spectral Distribution, ESD）定义为

    $$
    F_n(x) = \frac{1}{n} \#\{ i : \lambda_i \le x \} = \frac{1}{n} \sum_{i=1}^{n} \mathbf{1}_{\{\lambda_i \le x\}}.
    $$

    对应的经验谱测度为 $\mu_n = \frac{1}{n} \sum_{i=1}^{n} \delta_{\lambda_i}$。

!!! note "注"
    三大高斯系综的统一框架可以通过 **$\beta$-系综**（$\beta$-ensemble）给出：对 $\beta > 0$，联合特征值密度为

    $$
    p(\lambda_1, \ldots, \lambda_n) = Z_{n,\beta}^{-1} \prod_{i < j} |\lambda_i - \lambda_j|^\beta \prod_{i=1}^{n} e^{-\frac{n\beta}{4} \lambda_i^2}.
    $$

!!! theorem "定理 23.1 (GOE/GUE 联合特征值密度)"
    设 $M$ 为 GOE（$\beta=1$）或 GUE（$\beta=2$）矩阵，则其特征值 $\lambda_1, \ldots, \lambda_n$ 的联合概率密度函数为

    $$
    p(\lambda_1, \ldots, \lambda_n) = Z_{n,\beta}^{-1} \prod_{1 \le i < j \le n} |\lambda_i - \lambda_j|^\beta \cdot \prod_{i=1}^{n} e^{-\frac{n\beta}{4} \lambda_i^2},
    $$

    其中 $Z_{n,\beta}$ 为归一化常数。

??? proof "证明"
    以 GUE（$\beta = 2$）为例。设 $M = U \Lambda U^*$，其中 $U \in U(n)$，$\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$。对密度 $f(M) \propto e^{-\frac{n}{2}\operatorname{tr}(M^2)}$ 做变量替换，由 Jacobi 公式，变换的 Jacobi 行列式为

    $$
    |J| = \prod_{i < j} (\lambda_i - \lambda_j)^2.
    $$

    注意 $\operatorname{tr}(M^2) = \operatorname{tr}(\Lambda^2) = \sum_i \lambda_i^2$，对 $U$ 上的 Haar 测度积分后得到

    $$
    p(\lambda_1, \ldots, \lambda_n) \propto \prod_{i < j} (\lambda_i - \lambda_j)^2 \cdot \prod_i e^{-\frac{n}{2}\lambda_i^2}.
    $$

    此即 $\beta = 2$ 的情形。GOE（$\beta=1$）的推导类似，区别在于 Jacobi 行列式中指数为 $1$。$\blacksquare$

!!! example "例 23.1"
    **验证 $2 \times 2$ GUE 特征值密度。**

    设 $M = \begin{pmatrix} a & z \\ \bar{z} & b \end{pmatrix}$，其中 $a, b \sim N(0, 1/n)$，$z = x + iy$，$x, y \sim N(0, 1/(2n))$。特征值为

    $$
    \lambda_\pm = \frac{a+b}{2} \pm \sqrt{\left(\frac{a-b}{2}\right)^2 + |z|^2}.
    $$

    做变量替换 $(a, b, x, y) \to (\lambda_+, \lambda_-, \theta)$，其中 $\theta$ 为特征向量参数，可以验证联合密度中出现因子 $|\lambda_+ - \lambda_-|^2$，与定理 23.1 一致。

---

## 23.2 Wigner 矩阵与半圆律

<div class="context-flow" markdown>

**核心定理**：Wigner 矩阵（对称+独立+零均值+方差1）的经验谱分布 → $\rho_{sc}(x) = \frac{1}{2\pi}\sqrt{4-x^2}$（半圆）

**证明路线**：矩量法——$\frac{1}{n}\mathbb{E}[\text{tr}(W^{2m})]$ → 闭合路径计数 → **Catalan 数** $C_m$ → 唯一确定半圆分布

</div>

Wigner 半圆律是随机矩阵理论中最基本的极限定理，它描述了 Wigner 矩阵经验谱分布的极限行为。

!!! definition "定义 23.7 (Wigner 矩阵)"
    **Wigner 矩阵**（Wigner matrix）是 $n \times n$ Hermite（或实对称）随机矩阵 $W_n = \frac{1}{\sqrt{n}}(X_{ij})_{1 \le i,j \le n}$，其中：

    1. $\{X_{ij} : i \le j\}$ 独立；
    2. 对角元素 $X_{ii}$ 独立同分布，$\mathbb{E}[X_{ii}] = 0$，$\mathbb{E}[X_{ii}^2] < \infty$；
    3. 非对角元素 $X_{ij}$（$i < j$）独立同分布，$\mathbb{E}[X_{ij}] = 0$，$\mathbb{E}[|X_{ij}|^2] = 1$；
    4. $X_{ji} = \overline{X_{ij}}$。

!!! theorem "定理 23.2 (Wigner 半圆律)"
    设 $W_n$ 为 Wigner 矩阵，令 $\mu_n = \frac{1}{n}\sum_{i=1}^{n} \delta_{\lambda_i}$ 为其经验谱测度。则当 $n \to \infty$ 时，$\mu_n$ 几乎必然弱收敛到 **半圆分布**（semicircle distribution）$\mu_{sc}$，其密度为

    $$
    \rho_{sc}(x) = \frac{1}{2\pi} \sqrt{4 - x^2} \cdot \mathbf{1}_{[-2, 2]}(x).
    $$

??? proof "证明"
    **矩量法**（Method of Moments）的核心思路如下。

    **第一步：计算矩量。** 需要证明对任意正整数 $k$，

    $$
    \frac{1}{n} \mathbb{E}\!\left[\operatorname{tr}(W_n^k)\right] \to m_k = \int x^k \, \rho_{sc}(x) \, dx.
    $$

    展开 $\operatorname{tr}(W_n^k) = \frac{1}{n^{k/2}} \sum_{i_1, \ldots, i_k} X_{i_1 i_2} X_{i_2 i_3} \cdots X_{i_k i_1}$，每一项对应一条闭合路径 $(i_1, i_2, \ldots, i_k, i_1)$。

    **第二步：图论组合分析。** 由于 $\mathbb{E}[X_{ij}] = 0$，只有每条边至少被经过两次的路径才有非零贡献。对 $k$ 为奇数，满足条件的路径不存在，故 $m_k = 0$。对 $k = 2m$，主要贡献来自每条边恰好被经过两次的路径，这类路径与 **Catalan 数**（Catalan number）$C_m$ 一一对应。

    **第三步：确认 Catalan 数。** 可以证明

    $$
    m_{2m} = C_m = \frac{1}{m+1}\binom{2m}{m},
    $$

    而半圆分布的 $2m$ 阶矩恰好等于 Catalan 数。

    **第四步：矩量唯一确定分布。** 由于半圆分布的支撑有界（$[-2, 2]$），其矩量序列唯一确定该分布。

    **第五步：几乎必然收敛。** 利用方差估计 $\operatorname{Var}\!\left(\frac{1}{n}\operatorname{tr}(W_n^k)\right) = O(n^{-2})$，由 Borel-Cantelli 引理可将期望收敛提升为几乎必然收敛。$\blacksquare$

!!! theorem "定理 23.3 (半圆律的 Stieltjes 变换刻画)"
    半圆分布 $\mu_{sc}$ 的 Stieltjes 变换为

    $$
    s(z) = \int \frac{\rho_{sc}(x)}{x - z} \, dx = \frac{-z + \sqrt{z^2 - 4}}{2}, \quad z \in \mathbb{C}^+,
    $$

    其中取使得 $\operatorname{Im}(s(z)) > 0$（当 $\operatorname{Im}(z) > 0$）的分支。等价地，$s(z)$ 满足方程

    $$
    s(z)^2 + z \, s(z) + 1 = 0.
    $$

??? proof "证明"
    直接计算。设 $z \in \mathbb{C}^+$，则

    $$
    s(z) = \frac{1}{2\pi} \int_{-2}^{2} \frac{\sqrt{4 - x^2}}{x - z} \, dx.
    $$

    做替换 $x = 2\cos\theta$，$dx = -2\sin\theta \, d\theta$，$\sqrt{4 - x^2} = 2\sin\theta$，则

    $$
    s(z) = \frac{1}{2\pi} \int_0^{\pi} \frac{4\sin^2\theta}{2\cos\theta - z} \, (-2) \, d\theta \cdot \frac{1}{-1} = \frac{2}{\pi} \int_0^{\pi} \frac{\sin^2\theta}{z - 2\cos\theta} \, d\theta.
    $$

    利用留数定理（令 $w = e^{i\theta}$，化为围道积分），可以计算得

    $$
    s(z) = \frac{-z + \sqrt{z^2 - 4}}{2}.
    $$

    验证：$s^2 + zs + 1 = \frac{z^2 - 2z\sqrt{z^2-4} + z^2 - 4}{4} + \frac{-z^2 + z\sqrt{z^2-4}}{2} + 1 = 0$。$\blacksquare$

!!! example "例 23.2"
    **数值验证半圆律。**

    取 $n = 1000$，生成 GOE 矩阵 $M = \frac{1}{\sqrt{n}} A$，其中 $A$ 为对称矩阵，上三角元素独立标准正态。计算特征值并绘制直方图，观察直方图与 $\rho_{sc}(x) = \frac{1}{2\pi}\sqrt{4 - x^2}$ 的拟合。实验表明即使在 $n = 1000$ 时，经验谱分布已非常接近半圆分布。

!!! example "例 23.3"
    **计算半圆分布的矩量。**

    半圆分布的奇数阶矩为零（对称性）。偶数阶矩：

    $$
    m_2 = \int_{-2}^{2} x^2 \cdot \frac{\sqrt{4 - x^2}}{2\pi} \, dx = 1, \quad m_4 = \int_{-2}^{2} x^4 \cdot \frac{\sqrt{4 - x^2}}{2\pi} \, dx = 2.
    $$

    一般地，$m_{2k} = C_k = \frac{1}{k+1}\binom{2k}{k}$，其中 $C_k$ 为第 $k$ 个 Catalan 数。前几项：$C_0=1, C_1=1, C_2=2, C_3=5, C_4=14$。

---

## 23.3 样本协方差矩阵与 Marchenko-Pastur 律

<div class="context-flow" markdown>

**从 Wigner 到 Wishart**：半圆律 = 对称矩阵的谱极限 → **MP 律** = $\frac{1}{n}X^TX$ 的谱极限，支撑 $[(1-\sqrt{y})^2, (1+\sqrt{y})^2]$，$y = p/n$

**高维统计的基石**：$p/n \to y > 0$ 时样本协方差矩阵与真实协方差偏差巨大 → 经典统计理论失效 → 链接 Ch25 PCA

</div>

当我们从高维总体中抽取样本时，样本协方差矩阵的谱行为由 Marchenko-Pastur 律刻画。

!!! definition "定义 23.8 (样本协方差矩阵)"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_n \in \mathbb{R}^p$ 为独立同分布随机向量，$\mathbb{E}[\mathbf{x}_i] = \mathbf{0}$，$\operatorname{Cov}(\mathbf{x}_i) = \Sigma$。**样本协方差矩阵**（sample covariance matrix）定义为

    $$
    S_n = \frac{1}{n} \sum_{i=1}^{n} \mathbf{x}_i \mathbf{x}_i^T = \frac{1}{n} X^T X,
    $$

    其中 $X = (\mathbf{x}_1, \ldots, \mathbf{x}_n)^T$ 为 $n \times p$ 数据矩阵。

!!! theorem "定理 23.4 (Marchenko-Pastur 律)"
    设 $X$ 为 $n \times p$ 矩阵，元素 $X_{ij}$ 独立同分布，$\mathbb{E}[X_{ij}] = 0$，$\mathbb{E}[X_{ij}^2] = 1$。令 $S_n = \frac{1}{n} X^T X$，$\gamma = p/n \to y \in (0, \infty)$。则 $S_n$ 的经验谱分布几乎必然弱收敛到 **Marchenko-Pastur 分布**（Marchenko-Pastur distribution）$\mu_{MP}$，其密度为

    $$
    \rho_{MP}(x) = \frac{1}{2\pi xy} \sqrt{(\lambda_+ - x)(x - \lambda_-)} \cdot \mathbf{1}_{[\lambda_-, \lambda_+]}(x),
    $$

    其中 $\lambda_\pm = (1 \pm \sqrt{y})^2$。当 $y > 1$ 时，$\mu_{MP}$ 在 $x = 0$ 处有质量 $1 - 1/y$ 的点质量。

??? proof "证明"
    **Stieltjes 变换方法。** 令 $s_n(z)$ 为 $S_n$ 经验谱分布的 Stieltjes 变换。利用恒等式

    $$
    s_n(z) = \frac{1}{p} \operatorname{tr}\!\left( S_n - zI \right)^{-1},
    $$

    以及矩阵恒等式

    $$
    \frac{1}{n} X^T X - zI_p = -z \left( I_p - \frac{1}{nz} X^T X \right),
    $$

    结合 Sherman-Morrison-Woodbury 公式逐列删除技巧，可以证明 $s_n(z)$ 收敛到满足如下方程的 $s(z)$：

    $$
    s(z) = \frac{1}{-z + y/(1 + s(z) \cdot z)} \cdot \frac{1}{1},
    $$

    化简得

    $$
    z^2 s(z) + z(1 - y - z \cdot 1) \cdot s(z) \text{ 的隐式方程为 } s(z) = \frac{1}{(1-y) - z - y z s(z)}.
    $$

    更精确地，Marchenko-Pastur 分布的 Stieltjes 变换满足

    $$
    y z s^2(z) + (z - 1 + y) s(z) + 1 = 0.
    $$

    解此方程并验证 Stieltjes 反转公式 $\rho(x) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \operatorname{Im} s(x + i\eta)$ 即可还原密度 $\rho_{MP}$。$\blacksquare$

!!! theorem "定理 23.5 (一般总体的 Marchenko-Pastur 律)"
    设 $X$ 为 $n \times p$ 矩阵，元素独立同分布，$\mathbb{E}[X_{ij}] = 0$，$\mathbb{E}[X_{ij}^2] = 1$。设总体协方差 $\Sigma$ 的经验谱分布收敛到 $H$。令 $S_n = \frac{1}{n} X \Sigma X^T$，$p/n \to y$。则 $S_n$ 的极限谱分布 $F$ 的 Stieltjes 变换 $s(z)$ 满足

    $$
    s(z) = \int \frac{1}{\tau(1 - y - y z s(z)) - z} \, dH(\tau).
    $$

??? proof "证明"
    证明思路与定理 23.4 类似，但在逐列删除时需要考虑 $\Sigma$ 的结构。利用 $S_n = \frac{1}{n} X \Sigma X^T$ 和预解矩阵的秩一扰动公式，经过集中不等式和截断论证，可以证明 $s_n(z)$ 满足的近似方程在 $n \to \infty$ 时收敛到上述确定性方程。详细证明见 Silverstein-Bai (1995)。$\blacksquare$

!!! example "例 23.4"
    **$y = 1$ 时的 Marchenko-Pastur 分布。**

    当 $p = n$（即 $y = 1$）时，$\lambda_- = 0$，$\lambda_+ = 4$，密度为

    $$
    \rho_{MP}(x) = \frac{1}{2\pi x} \sqrt{(4 - x) \cdot x} = \frac{1}{2\pi} \sqrt{\frac{4 - x}{x}}, \quad x \in (0, 4].
    $$

    注意此时 $\rho_{MP}(x) \to \infty$ 当 $x \to 0^+$，即在零点附近特征值密度趋于无穷。

!!! example "例 23.5"
    **比较不同 $y$ 值下的 MP 分布。**

    - $y = 0.2$：$\lambda_- = (1 - \sqrt{0.2})^2 \approx 0.106$，$\lambda_+ = (1 + \sqrt{0.2})^2 \approx 2.294$，分布集中在 $1$ 附近。
    - $y = 1$：$\lambda_- = 0$，$\lambda_+ = 4$，分布扩展到 $[0, 4]$。
    - $y = 5$：$\lambda_- = 0$，$\lambda_+ = (1 + \sqrt{5})^2 \approx 10.47$，在 $0$ 处有点质量 $1 - 1/5 = 0.8$。

    随着 $y$ 增大（样本量相对于维数减小），谱分布越来越展宽，反映了高维噪声的放大效应。

---

## 23.4 经验谱分布与 Stieltjes 变换方法

<div class="context-flow" markdown>

**方法论**：矩量法适合证明存在性 → **Stieltjes 变换** $s(z) = \int \frac{d\mu(x)}{x-z}$ 才是计算利器 → 反转公式 $\rho(x) = \frac{1}{\pi}\text{Im}\,s(x+i0^+)$ 从变换还原密度

</div>

Stieltjes 变换是研究随机矩阵极限谱分布的核心分析工具。

!!! definition "定义 23.9 (Stieltjes 变换)"
    设 $\mu$ 为 $\mathbb{R}$ 上的概率测度。$\mu$ 的 **Stieltjes 变换**（Stieltjes transform）定义为

    $$
    s_\mu(z) = \int_{\mathbb{R}} \frac{1}{x - z} \, d\mu(x), \quad z \in \mathbb{C}^+ = \{z : \operatorname{Im}(z) > 0\}.
    $$

    $s_\mu$ 是 $\mathbb{C}^+$ 上的全纯函数，且 $\operatorname{Im}(s_\mu(z)) > 0$。

!!! theorem "定理 23.6 (Stieltjes 反转公式)"
    设 $\mu$ 为概率测度，$s(z)$ 为其 Stieltjes 变换。若 $\mu$ 在 $(a, b)$ 上有连续密度 $\rho$，则

    $$
    \rho(x) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \operatorname{Im}\, s(x + i\eta), \quad x \in (a, b).
    $$

    更一般地，对 $\mu$ 的连续点 $a < b$，

    $$
    \mu\!\left((a, b)\right) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \int_a^b \operatorname{Im}\, s(x + i\eta) \, dx.
    $$

??? proof "证明"
    由 $s(x + i\eta) = \int \frac{1}{t - x - i\eta} \, d\mu(t)$，取虚部得

    $$
    \operatorname{Im}\, s(x + i\eta) = \int \frac{\eta}{(t - x)^2 + \eta^2} \, d\mu(t).
    $$

    注意 $\frac{\eta}{\pi((t-x)^2 + \eta^2)}$ 是以 $x$ 为中心、半宽 $\eta$ 的 Cauchy（Poisson）核，当 $\eta \to 0$ 时趋于 $\delta(t - x)$。因此

    $$
    \frac{1}{\pi} \operatorname{Im}\, s(x + i\eta) = \int \frac{1}{\pi} \frac{\eta}{(t - x)^2 + \eta^2} \, d\mu(t) \to \rho(x),
    $$

    其中收敛在 $\rho$ 的连续点成立。$\blacksquare$

!!! theorem "定理 23.7 (Stieltjes 变换的连续性定理)"
    设 $\{\mu_n\}$ 为概率测度序列，$s_n(z)$ 为对应的 Stieltjes 变换。若对所有 $z \in \mathbb{C}^+$，$s_n(z) \to s(z)$，且 $s(z)$ 是某个概率测度 $\mu$ 的 Stieltjes 变换，则 $\mu_n \xrightarrow{w} \mu$（弱收敛）。

??? proof "证明"
    Stieltjes 变换与分布函数之间存在一一对应关系（在适当条件下）。$s_n(z) \to s(z)$ 逐点收敛蕴含了矩量的收敛（通过 Laurent 展开），进而由矩量问题的唯一性得到弱收敛。严格证明需要用到紧性论证（Helly 选择定理）以及 Stieltjes 变换唯一确定测度的性质。$\blacksquare$

!!! example "例 23.6"
    **用 Stieltjes 变换验证半圆律。**

    对 Wigner 矩阵 $W_n$，其预解矩阵的归一化迹 $s_n(z) = \frac{1}{n}\operatorname{tr}(W_n - zI)^{-1}$ 满足近似方程

    $$
    s_n(z) \approx \frac{1}{-z - s_n(z)},
    $$

    令 $n \to \infty$，$s(z)$ 满足 $s = \frac{1}{-z - s}$，即 $s^2 + zs + 1 = 0$，解为 $s(z) = \frac{-z + \sqrt{z^2 - 4}}{2}$，恰为半圆分布的 Stieltjes 变换。

---

## 23.5 特征值间距与排斥现象

<div class="context-flow" markdown>

**微观行为**：半圆律/MP律 = 宏观（密度）→ 间距统计 = 微观

**排斥**：$p(s) \sim s^\beta$（$s \to 0$）vs Poisson $p(s) = e^{-s}$ → 特征值"互相推开"，这是随机矩阵与独立随机变量的本质区别

</div>

随机矩阵的一个显著特征是特征值之间的排斥效应：特征值倾向于互相远离，其间距统计与独立随机变量有本质区别。

!!! definition "定义 23.10 (标准化间距)"
    设 $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$ 为随机矩阵的有序特征值。在谱内部点 $E$ 处，局部特征值密度为 $\rho(E)$。**标准化间距**（normalized spacing）定义为

    $$
    s_i = n \rho(E) (\lambda_{i+1} - \lambda_i),
    $$

    使得标准化后间距的均值为 $1$。

!!! theorem "定理 23.8 (Wigner 特征值间距统计)"
    对于 GUE（$\beta = 2$）矩阵，在谱内部，标准化间距的分布趋近于 **Gaudin 分布**。对于小间距 $s \to 0$，间距概率密度

    $$
    p(s) \sim c_\beta \, s^\beta, \quad s \to 0,
    $$

    其中 $\beta$ 为 Dyson 指标（GOE: $\beta=1$，GUE: $\beta=2$，GSE: $\beta=4$）。这表明小间距出现的概率极小——特征值相互排斥。

??? proof "证明"
    以 GUE 为例。利用行列式点过程的结构：GUE 特征值构成以 sine 核

    $$
    K(x, y) = \frac{\sin(\pi(x - y))}{\pi(x - y)}
    $$

    为关联核的行列式点过程。间距分布函数为

    $$
    P(s > t) = \det(I - K_t),
    $$

    其中 $K_t$ 是 sine 核在区间 $[0, t]$ 上的限制算子。通过 Fredholm 行列式展开，当 $t \to 0$ 时，

    $$
    P(s \le t) = 1 - \det(I - K_t) \sim \frac{(\pi t)^2}{2} + O(t^4),
    $$

    故 $p(t) \sim \pi^2 t$，即 $p(s) \propto s^2$（$\beta = 2$）。$\blacksquare$

!!! theorem "定理 23.9 (Wigner-Dyson-Mehta 间距分布近似)"
    在实际应用中，常用以下近似公式（Wigner surmise）来近似间距分布：

    - GOE（$\beta = 1$）：$p(s) = \frac{\pi}{2} s \, e^{-\pi s^2 / 4}$；
    - GUE（$\beta = 2$）：$p(s) = \frac{32}{\pi^2} s^2 \, e^{-4s^2 / \pi}$；
    - GSE（$\beta = 4$）：$p(s) = \frac{2^{18}}{3^6 \pi^3} s^4 \, e^{-64 s^2 / (9\pi)}$。

    这些公式精确描述了 $2 \times 2$ 矩阵的间距分布，对大矩阵也是极好的近似。

??? proof "证明"
    以 GOE 的 $2 \times 2$ 情形为例。设 $M = \begin{pmatrix} a & b \\ b & c \end{pmatrix}$，$a, c \sim N(0,1)$，$b \sim N(0, 1/2)$。特征值间距 $s = \lambda_+ - \lambda_- = \sqrt{(a-c)^2 + 4b^2}$。令 $u = a - c$，$v = 2b$，则 $u \sim N(0, 2)$，$v \sim N(0, 2)$，$s = \sqrt{u^2 + v^2}$。转化为极坐标：$p(s) = \frac{s}{2} e^{-s^2/4}$。经适当归一化使 $\langle s \rangle = 1$，得到 Wigner surmise $p(s) = \frac{\pi}{2} s \, e^{-\pi s^2/4}$。$\blacksquare$

!!! example "例 23.7"
    **Poisson 与 GUE 间距统计的比较。**

    - 独立随机特征值（如对角随机矩阵）的间距服从指数分布 $p(s) = e^{-s}$（Poisson 统计），$p(0) = 1$。
    - GUE 的间距为 $p(s) \approx \frac{32}{\pi^2} s^2 e^{-4s^2/\pi}$，$p(0) = 0$。

    GUE 在 $s = 0$ 处密度为零，体现了特征值排斥；而 Poisson 统计在零间距处密度最大，说明独立特征值可以任意接近。

---

## 23.6 Tracy-Widom 分布

<div class="context-flow" markdown>

**精细尺度**：半圆律(宏观) → 间距(微观) → **最大特征值涨落**(边缘) · $\lambda_{\max} \approx 2 + n^{-2/3} \cdot F_\beta$ → Airy 核 + Painleve II 方程

**普适性**：$F_\beta$ 不依赖于矩阵元素的具体分布，只依赖对称性类 $\beta$ → 链接 §23.7 统计检验

</div>

半圆律描述了特征值的整体分布，而 Tracy-Widom 分布刻画了最大特征值的涨落。

!!! definition "定义 23.11 (Tracy-Widom 分布)"
    **Tracy-Widom 分布**（Tracy-Widom distribution）$F_\beta$（$\beta = 1, 2, 4$）描述了随机矩阵最大特征值经适当中心化和缩放后的极限分布。对于 $\beta = 2$（GUE），分布函数为

    $$
    F_2(s) = \exp\!\left( -\int_s^{\infty} (x - s) \, q(x)^2 \, dx \right),
    $$

    其中 $q(x)$ 为 Painleve II 方程 $q''(x) = xq(x) + 2q(x)^3$ 满足 Airy 衰减条件 $q(x) \sim \operatorname{Ai}(x)$（$x \to +\infty$）的唯一解。

!!! theorem "定理 23.10 (GUE 最大特征值的 Tracy-Widom 极限)"
    设 $M_n$ 为 $n \times n$ GUE 矩阵，$\lambda_{\max}$ 为其最大特征值。则

    $$
    \frac{\lambda_{\max} - 2}{n^{-2/3}} \xrightarrow{d} F_2,
    $$

    即 $P\!\left( n^{2/3}(\lambda_{\max} - 2) \le s \right) \to F_2(s)$，其中 $F_2$ 为 Tracy-Widom 分布。

??? proof "证明"
    **证明思路。** GUE 特征值构成行列式点过程，关联核为

    $$
    K_n(x, y) = \sum_{k=0}^{n-1} \varphi_k(x) \varphi_k(y),
    $$

    其中 $\varphi_k$ 为标准化 Hermite 函数。最大特征值的分布为

    $$
    P(\lambda_{\max} \le t) = \det(I - K_n)|_{L^2(t, \infty)}.
    $$

    在谱边缘 $t = 2 + s n^{-2/3}$ 处，利用 Plancherel-Rotach 渐近公式，$K_n$ 在适当缩放下收敛到 **Airy 核**

    $$
    K_{\text{Airy}}(x, y) = \frac{\operatorname{Ai}(x)\operatorname{Ai}'(y) - \operatorname{Ai}'(x)\operatorname{Ai}(y)}{x - y}.
    $$

    因此 $P(\lambda_{\max} \le 2 + sn^{-2/3}) \to \det(I - K_{\text{Airy}})|_{L^2(s, \infty)} = F_2(s)$。$\blacksquare$

!!! theorem "定理 23.11 (GOE 最大特征值的 Tracy-Widom 极限)"
    设 $M_n$ 为 $n \times n$ GOE 矩阵，则

    $$
    \frac{\lambda_{\max} - 2}{n^{-2/3}} \xrightarrow{d} F_1,
    $$

    其中 $F_1(s) = \exp\!\left( -\frac{1}{2}\int_s^\infty q(x) \, dx \right) \cdot F_2(s)^{1/2}$，$q$ 为上述 Painleve II 解。

??? proof "证明"
    GOE 特征值构成 Pfaffian 点过程（而非行列式点过程）。在谱边缘缩放后，关联核收敛到 Airy 核的对称化版本。最大特征值的分布可以写成 Fredholm Pfaffian，其极限给出 $F_1$。$\blacksquare$

!!! example "例 23.8"
    **Tracy-Widom 分布的数值特征。**

    $F_2$ 分布的数值特征：
    - 均值 $\approx -1.771$
    - 标准差 $\approx 0.813$
    - 偏度 $\approx 0.224$
    - 峰度 $\approx 0.093$

    该分布是左偏的，这意味着最大特征值倾向于略低于其均值 $2$。$F_1$ 分布比 $F_2$ 的涨落更大（标准差 $\approx 1.268$），因为实对称矩阵的自由度更少。

---

## 23.7 随机矩阵在统计中的应用

<div class="context-flow" markdown>

**落地**：**BBP 相变**——信号强度 $\theta > \sqrt{p/n}$ 时最大特征值从 MP 边缘"弹出" → 信号检测的理论阈值 · Tracy-Widom 分布替代经典 $F$/$\chi^2$ 分布用于高维假设检验

</div>

随机矩阵理论为高维统计提供了理论基础和实用工具。

!!! definition "定义 23.12 (高维渐近框架)"
    在 **高维渐近框架**（high-dimensional asymptotic framework）中，数据维数 $p$ 和样本量 $n$ 同时趋于无穷，且比值 $p/n \to y \in (0, \infty)$。这与经典统计中 $p$ 固定、$n \to \infty$ 的框架根本不同。

<div class="context-flow" markdown>

**洞察**：BBP 相变揭示了一个深刻的"信息论极限"——信号强度 $\theta$ 低于 $\sqrt{y}$ 时，**任何方法**都无法从样本协方差的特征值中检测到信号

</div>

!!! theorem "定理 23.12 (Baik-Ben Arous-Peche 相变)"
    **BBP 相变**（BBP phase transition）：设总体协方差矩阵 $\Sigma = I + \theta \mathbf{v}\mathbf{v}^T$（秩一扰动），$p/n \to y$。令 $\ell_1$ 为样本协方差矩阵 $S_n$ 的最大特征值。则：

    - 若 $\theta < \sqrt{y}$，则 $\ell_1 \to (1 + \sqrt{y})^2$（与 $\theta = 0$ 时相同）；
    - 若 $\theta > \sqrt{y}$，则 $\ell_1 \to (1 + \theta)(1 + y/\theta) > (1 + \sqrt{y})^2$。

    临界值 $\theta_c = \sqrt{y}$ 标志着信号是否能被检测到的相变。

??? proof "证明"
    **证明思路。** 利用 Stieltjes 变换方法。在 $\Sigma = I + \theta \mathbf{v}\mathbf{v}^T$ 下，$S_n$ 的极限谱分布仍由 Marchenko-Pastur 律描述（秩一扰动不影响极限谱分布），但最大特征值的行为取决于 $\theta$ 的大小。

    关键步骤是分析预解矩阵 $(S_n - zI)^{-1}$ 在谱边缘外的行为。当 $z$ 在 $\lambda_+ = (1+\sqrt{y})^2$ 外时，$\frac{1}{n}\operatorname{tr}(S_n - zI)^{-1} \to s(z)$。利用矩阵摄动公式，样本中最大特征值满足

    $$
    1 + \theta \cdot \frac{1}{p}\operatorname{tr}\!\left(\frac{S_n^{(0)}}{S_n^{(0)} - \ell_1 I}\right) \approx 0,
    $$

    其中 $S_n^{(0)}$ 为去掉信号后的矩阵。当 $\theta > \sqrt{y}$ 时此方程在 $(1+\sqrt{y})^2$ 外有解，即最大特征值从 Marchenko-Pastur 支撑中"弹出"。$\blacksquare$

!!! example "例 23.9"
    **信号检测问题。**

    在无线通信中，接收信号模型为 $\mathbf{x} = \sqrt{\theta} \, s \, \mathbf{a} + \mathbf{n}$，其中 $s$ 为信号，$\mathbf{a}$ 为方向向量，$\mathbf{n} \sim N(\mathbf{0}, I)$。$n$ 次观测后样本协方差矩阵为 $S_n = \frac{1}{n}XX^T$。由 BBP 相变，当信噪比 $\theta > \sqrt{p/n}$ 时，$S_n$ 的最大特征值将显著偏离 Marchenko-Pastur 分布的上边缘，从而可以检测到信号。

!!! example "例 23.10"
    **Roy 最大根检验的修正。**

    经典的 Roy 最大根检验统计量为 $\ell_1(S_1 S_2^{-1})$。在 $p, n \to \infty$、$p/n \to y$ 的高维框架下，该统计量在零假设下的极限分布不再是经典的 Roy 分布，而是 Tracy-Widom 分布 $F_1$。因此检验的拒绝域应基于 Tracy-Widom 分位数而非传统表格。

---

## 23.8 自由概率简介

<div class="context-flow" markdown>

**代数化**：随机矩阵 $A_n, B_n$ 独立 → 渐近**自由**（Voiculescu）→ $A+B$ 的极限谱由 **$R$-变换**可加性计算 → 自由 CLT：极限是半圆分布（对比经典 CLT → 正态分布）

**统一视角**：半圆律 = 自由概率的中心极限定理

</div>

自由概率论（Free Probability）是研究非交换随机变量的数学理论，它为随机矩阵的渐近谱行为提供了代数化框架。

!!! definition "定义 23.13 (非交换概率空间)"
    **非交换概率空间**（noncommutative probability space）是一对 $(\mathcal{A}, \varphi)$，其中 $\mathcal{A}$ 为一个含单位元的代数（不一定交换），$\varphi: \mathcal{A} \to \mathbb{C}$ 为一个线性泛函，满足 $\varphi(1) = 1$（称为迹态）。$\mathcal{A}$ 中的元素称为非交换随机变量。

!!! definition "定义 23.14 (自由独立性)"
    在非交换概率空间 $(\mathcal{A}, \varphi)$ 中，子代数 $\mathcal{A}_1, \ldots, \mathcal{A}_k$ 称为 **自由独立的**（freely independent），若对于任意 $a_j \in \mathcal{A}_{i_j}$（$j = 1, \ldots, m$），满足 $\varphi(a_j) = 0$ 且相邻元素来自不同子代数（$i_1 \ne i_2 \ne \cdots \ne i_m$）时，有

    $$
    \varphi(a_1 a_2 \cdots a_m) = 0.
    $$

    自由独立性是经典独立性在非交换设定下的类比，但两者有本质区别。

!!! theorem "定理 23.13 (Voiculescu 渐近自由性定理)"
    设 $A_n, B_n$ 为 $n \times n$ 独立随机矩阵，$A_n$ 为 Wigner 矩阵，$B_n = U_n D_n U_n^*$，其中 $D_n$ 为确定性对角矩阵（经验谱分布收敛到 $\nu$），$U_n$ 为 Haar 酉矩阵且与 $A_n$ 独立。则当 $n \to \infty$ 时，$A_n$ 和 $B_n$ **渐近自由**，即关于归一化迹 $\varphi(\cdot) = \frac{1}{n}\operatorname{tr}(\cdot)$，它们的混合矩量满足自由独立性的代数关系。

??? proof "证明"
    **证明思路。** 需要验证对于中心化后的交替乘积，归一化迹趋于零。即对 $p(A_n) = A_n^k - \varphi(A_n^k)I$ 和 $q(B_n) = B_n^l - \varphi(B_n^l)I$，需证

    $$
    \frac{1}{n}\operatorname{tr}\!\left(p_1(A_n) q_1(B_n) p_2(A_n) q_2(B_n) \cdots \right) \to 0.
    $$

    关键工具是 Weingarten 积分公式，它给出了 Haar 酉矩阵元素乘积对 $U_n$ 的积分。利用该公式，可以将上述迹展开为关于置换的求和，主要项之间的消去（由中心化保证）使得整体趋于零。$\blacksquare$

!!! definition "定义 23.15 (自由卷积)"
    设 $\mu, \nu$ 为 $\mathbb{R}$ 上的概率测度。若 $a, b$ 为非交换概率空间中的自伴元素，分布分别为 $\mu, \nu$，且 $a, b$ 自由独立，则 $a + b$ 的分布称为 $\mu$ 和 $\nu$ 的 **自由（加法）卷积**（free additive convolution），记为 $\mu \boxplus \nu$。

!!! theorem "定理 23.14 (自由卷积的 $R$-变换)"
    设 $\mu, \nu$ 为概率测度，$G_\mu(z) = \int \frac{d\mu(x)}{z - x}$ 为 Cauchy 变换（注意与 Stieltjes 变换差一个符号），$K_\mu$ 为其函数逆（$G_\mu(K_\mu(w)) = w$），$R_\mu(w) = K_\mu(w) - 1/w$。则

    $$
    R_{\mu \boxplus \nu}(w) = R_\mu(w) + R_\nu(w).
    $$

    即自由卷积下 $R$-变换是可加的，这类似于经典独立性下特征函数的对数可加性。

??? proof "证明"
    利用组合自由概率论（自由累积量理论）。定义自由累积量 $\kappa_n(\mu)$ 通过矩量-累积量公式

    $$
    m_n = \sum_{\pi \in NC(n)} \prod_{V \in \pi} \kappa_{|V|},
    $$

    其中求和遍历 $\{1, \ldots, n\}$ 的所有非交叉分割 $NC(n)$。可以证明自由独立性等价于混合自由累积量为零。因此 $\kappa_n(\mu \boxplus \nu) = \kappa_n(\mu) + \kappa_n(\nu)$。而 $R$-变换的 Laurent 展开系数恰好是自由累积量：$R_\mu(w) = \sum_{n=1}^\infty \kappa_n(\mu) w^{n-1}$。$\blacksquare$

!!! theorem "定理 23.15 (半圆律的自由中心极限定理)"
    设 $a_1, a_2, \ldots$ 为自由独立同分布的自伴非交换随机变量，$\varphi(a_i) = 0$，$\varphi(a_i^2) = 1$。则

    $$
    \frac{a_1 + a_2 + \cdots + a_n}{\sqrt{n}} \xrightarrow{d} \mu_{sc},
    $$

    即部分和的归一化极限分布为半圆分布。这是经典中心极限定理（极限为正态分布）的自由概率类比。

??? proof "证明"
    由自由累积量的可加性，$S_n = \frac{1}{\sqrt{n}}(a_1 + \cdots + a_n)$ 的自由累积量为 $\kappa_k(S_n) = n^{1 - k/2} \kappa_k(a_1)$。当 $n \to \infty$ 时，$\kappa_1(S_n) = 0$，$\kappa_2(S_n) = 1$，$\kappa_k(S_n) \to 0$（$k \ge 3$）。而半圆分布的自由累积量恰好是 $\kappa_2 = 1$、$\kappa_k = 0$（$k \ne 2$），这是因为半圆分布的 $R$-变换为 $R(w) = w$。$\blacksquare$

!!! example "例 23.11"
    **两个半圆分布的自由卷积。**

    设 $\mu = \mu_{sc}(0, 1)$（标准半圆分布，$R_\mu(w) = w$），$\nu = \mu_{sc}(0, 1)$。则

    $$
    R_{\mu \boxplus \nu}(w) = 2w,
    $$

    对应的分布为 $\mu_{sc}(0, \sqrt{2})$，即半径为 $2\sqrt{2}$、方差为 $2$ 的半圆分布。更一般地，方差为 $\sigma_1^2$ 和 $\sigma_2^2$ 的半圆分布的自由卷积仍为半圆分布，方差为 $\sigma_1^2 + \sigma_2^2$。

!!! example "例 23.12"
    **利用自由概率计算 $A + B$ 的极限谱。**

    设 $A_n$ 为 Wigner 矩阵（极限谱为半圆分布 $\mu_{sc}$），$B_n$ 为独立的确定性矩阵经 Haar 酉共轭后的矩阵，$B_n$ 的经验谱分布收敛到 Bernoulli 分布 $\nu = \frac{1}{2}\delta_{-1} + \frac{1}{2}\delta_1$。由渐近自由性，$A_n + B_n$ 的极限谱分布为 $\mu_{sc} \boxplus \nu$，可以通过 $R$-变换方法计算。$\nu$ 的 $R$-变换为 $R_\nu(w) = \frac{w}{1 - w^2}$（利用矩量-累积量关系），故 $R_{\mu_{sc} \boxplus \nu}(w) = w + \frac{w}{1 - w^2}$，再由反函数关系可以数值求解极限密度。

---

## 本章小结

本章介绍了随机矩阵理论的基本框架和核心结果：

1. **高斯系综**（GOE, GUE, GSE）作为随机矩阵的经典模型，其联合特征值密度具有优美的行列式/Pfaffian 结构。
2. **Wigner 半圆律**描述了 Wigner 矩阵经验谱分布的宏观极限。
3. **Marchenko-Pastur 律**刻画了高维样本协方差矩阵的谱行为，是高维统计的理论基石。
4. **Stieltjes 变换**是研究极限谱分布的核心分析工具。
5. 特征值**排斥现象**和 Wigner-Dyson 间距统计揭示了随机矩阵与独立随机变量的本质区别。
6. **Tracy-Widom 分布**描述了最大特征值的精细涨落。
7. **BBP 相变**为高维信号检测提供了理论阈值。
8. **自由概率论**为随机矩阵的渐近谱计算提供了代数化工具。
