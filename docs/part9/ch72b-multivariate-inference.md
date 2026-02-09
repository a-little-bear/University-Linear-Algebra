# 第 72B 章 多元统计推断

<div class="context-flow" markdown>

**前置**：Wishart 分布与逆 Wishart 分布(Ch72A) · 特征值与特征向量(Ch6) · 正定矩阵(Ch16) · 随机矩阵(Ch23) · 矩阵变换的 Jacobi 行列式(Ch72A)

**本章脉络**：矩阵 Beta/F 分布 $\to$ Hotelling $T^2$（多元 $t$ 检验）$\to$ Cochran 定理的矩阵版本 $\to$ Wilks' Lambda 分布 $\to$ MANOVA 四大检验统计量 $\to$ Roy 最大根的变分刻画 $\to$ Box M 检验 $\to$ 高维极限理论：Marchenko-Pastur 律 $\to$ Tracy-Widom 分布 $\to$ Spiked 协方差模型与 BBP 相变 $\to$ Fisher 信息度量 $\to$ 应用：金融、脑成像、Wishart 过程

**延伸**：Hotelling $T^2$ 是多元假设检验的基石；MANOVA 检验统计量将单变量 $F$ 检验推广到多响应变量情形；Marchenko-Pastur 律是高维统计和随机矩阵理论的里程碑；BBP 相变揭示了高维协方差估计中信号检测的基本极限；这些理论在金融风险管理、基因组学和无线通信中有直接应用

</div>

多元统计推断将经典的假设检验和参数估计从标量和向量推广到矩阵。其核心工具是矩阵值分布的特征值理论：从 Hotelling 的 $T^2$ 检验到 MANOVA 的四大检验统计量，都可以表示为特征值的函数。

当维度 $p$ 与样本量 $n$ 同时趋于无穷时，经典的有限维理论不再适用，取而代之的是随机矩阵理论的极限分布——Marchenko-Pastur 律描述了样本协方差矩阵特征值的整体分布，Tracy-Widom 分布刻画了最大特征值的波动，BBP 相变揭示了信号检测的基本极限。

本章将从经典多元统计推断出发，一直推进到现代高维随机矩阵理论。

---

## 72B.1 矩阵 Beta 和 F 分布

<div class="context-flow" markdown>

**核心问题**：如何将标量的 Beta 分布和 F 分布推广到矩阵值？这些分布在多元假设检验中起什么作用？

</div>

!!! definition "定义 72B.1 (矩阵 Beta 分布)"
    设 $W_1 \sim W_p(n_1, I_p)$ 和 $W_2 \sim W_p(n_2, I_p)$ 独立（$n_1, n_2 \geq p$）。定义 **I 型矩阵 Beta 变量**：

    $$B = (W_1 + W_2)^{-1/2} W_1 (W_1 + W_2)^{-1/2}$$

    或等价地（在适当意义下）：

    $$B = W_1(W_1 + W_2)^{-1}$$

    $B$ 服从参数为 $(n_1/2, n_2/2)$ 的**矩阵 Beta 分布** $\text{Beta}_p(n_1/2, n_2/2)$。$B$ 的特征值 $\beta_1, \ldots, \beta_p \in (0, 1)$。

!!! definition "定义 72B.2 (矩阵 F 分布)"
    **矩阵 F 变量**定义为：

    $$F = W_2^{-1/2} W_1 W_2^{-1/2}$$

    或 $F = W_1 W_2^{-1}$（非对称版本）。$F$ 的特征值为 $f_i = \beta_i / (1 - \beta_i)$。

    $p = 1$ 的特例：$F = (W_1/n_1)/(W_2/n_2) \cdot (n_2/n_1) \sim F(n_1, n_2)$。

!!! theorem "定理 72B.1 (矩阵 Beta 分布的密度)"
    $B \sim \text{Beta}_p(a, b)$（$a, b > (p-1)/2$）的密度（关于 $(0, I_p)$ 中正定矩阵的测度）为：

    $$f(B) = \frac{\Gamma_p(a+b)}{\Gamma_p(a)\Gamma_p(b)} |B|^{a-(p+1)/2} |I-B|^{b-(p+1)/2}$$

    其中 $0 < B < I_p$（即 $B$ 和 $I_p - B$ 均正定）。

??? proof "证明"
    设 $W = W_1 + W_2 \sim W_p(n_1+n_2, I_p)$（可加性），$B = W^{-1/2}W_1W^{-1/2}$。

    通过变量替换 $(W_1, W_2) \to (B, W)$——其中 $W_1 = W^{1/2}BW^{1/2}$，$W_2 = W^{1/2}(I-B)W^{1/2}$——可以证明 $B$ 和 $W$ 独立。

    Jacobi 行列式的计算涉及对称正定矩阵上的链式法则。最终结果是 $B$ 的边缘密度为矩阵 Beta 密度。

    与标量 Beta 密度 $f(x) \propto x^{a-1}(1-x)^{b-1}$ 对比，矩阵版本中的指数修正了 $(p+1)/2$ 以适应矩阵 Lebesgue 测度。

    $\blacksquare$

---

## 72B.2 Hotelling $T^2$ 分布

<div class="context-flow" markdown>

**核心问题**：如何检验多元正态总体的均值向量是否等于某个给定值？这是单变量 $t$ 检验的多元推广。

</div>

!!! definition "定义 72B.3 (Hotelling $T^2$ 统计量)"
    设 $\mathbf{x}_1, \ldots, \mathbf{x}_n \overset{\text{iid}}{\sim} N_p(\boldsymbol{\mu}, \Sigma)$，$\bar{\mathbf{x}} = \frac{1}{n}\sum_i \mathbf{x}_i$，$S = \frac{1}{n-1}\sum_i (\mathbf{x}_i - \bar{\mathbf{x}})(\mathbf{x}_i - \bar{\mathbf{x}})^T$。

    **Hotelling $T^2$ 统计量**检验 $H_0: \boldsymbol{\mu} = \boldsymbol{\mu}_0$：

    $$T^2 = n(\bar{\mathbf{x}} - \boldsymbol{\mu}_0)^T S^{-1} (\bar{\mathbf{x}} - \boldsymbol{\mu}_0)$$

!!! theorem "定理 72B.2 (Hotelling $T^2$ 的精确分布)"
    在原假设 $H_0: \boldsymbol{\mu} = \boldsymbol{\mu}_0$ 下：

    $$T^2 \sim \frac{(n-1)p}{n-p} F(p, n-p)$$

    即 $\frac{n-p}{(n-1)p} T^2 \sim F(p, n-p)$。

??? proof "证明"
    **第一步**：$\bar{\mathbf{x}} \sim N_p(\boldsymbol{\mu}_0, \Sigma/n)$（在 $H_0$ 下），且 $(n-1)S \sim W_p(n-1, \Sigma)$，并且 $\bar{\mathbf{x}}$ 和 $S$ 独立。

    独立性的证明：$\bar{\mathbf{x}} = n^{-1}\mathbf{1}^TX$，$(n-1)S = X^T(I - n^{-1}\mathbf{1}\mathbf{1}^T)X$。由于 $n^{-1}\mathbf{1}^T$ 和 $I - n^{-1}\mathbf{1}\mathbf{1}^T$ 正交（$n^{-1}\mathbf{1}^T(I - n^{-1}\mathbf{1}\mathbf{1}^T) = 0$），且 $X$ 的行独立正态，故二者独立。

    **第二步**：令 $\mathbf{z} = \sqrt{n}\Sigma^{-1/2}(\bar{\mathbf{x}} - \boldsymbol{\mu}_0) \sim N_p(0, I_p)$，$V = \Sigma^{-1/2}(n-1)S\Sigma^{-1/2} \sim W_p(n-1, I_p)$。

    $$T^2 = \mathbf{z}^T \left(\frac{V}{n-1}\right)^{-1} \mathbf{z} = (n-1)\mathbf{z}^T V^{-1}\mathbf{z}$$

    $\mathbf{z}^T\mathbf{z} \sim \chi^2_p$（因为 $\mathbf{z} \sim N_p(0, I_p)$），$V \sim W_p(n-1, I_p)$，且独立。

    **第三步**：利用 Wishart 分布的性质。$V = \sum_{j=1}^{n-1}\mathbf{u}_j\mathbf{u}_j^T$，$\mathbf{u}_j \overset{\text{iid}}{\sim} N_p(0, I_p)$。令 $A = (\mathbf{z}, \mathbf{u}_1, \ldots, \mathbf{u}_{n-1})^T \in \mathbb{R}^{n \times p}$。经过正交变换和条件分析：

    $$T^2 = (n-1) \cdot \frac{\mathbf{z}^TV^{-1}\mathbf{z}}{1}$$

    可以证明 $\frac{T^2/p}{(n-1-p+1)/(n-1)} = \frac{n-p}{(n-1)p}T^2 \sim F(p, n-p)$。

    这利用了以下事实：$\mathbf{z}^TV^{-1}\mathbf{z}$ 在适当正交变换下等价于 $\chi^2_p/(\chi^2_{n-p}/(...))$ 形式的比值。

    $\blacksquare$

!!! example "例 72B.1"
    **药物临床试验**。测量 $n = 30$ 名患者治疗前后 $p = 3$ 个生理指标的差值 $\mathbf{d}_i = \mathbf{x}_i^{\text{after}} - \mathbf{x}_i^{\text{before}}$。检验 $H_0: \boldsymbol{\mu}_d = 0$（药物无效）。

    $$T^2 = 30 \cdot \bar{\mathbf{d}}^T S_d^{-1} \bar{\mathbf{d}}$$

    在 $H_0$ 下，$\frac{n-p}{(n-1)p}T^2 = \frac{27}{87}T^2 \sim F(3, 27)$。

    若 $T^2 = 15.6$，则 $F = 27 \times 15.6/87 = 4.84$。$F(3, 27)$ 的 5% 临界值约为 2.96，拒绝 $H_0$。

    注意：$T^2$ 同时考虑了三个指标的变化及其相关性，避免了逐个检验带来的多重比较问题。

---

## 72B.3 Cochran 定理的矩阵版本

<div class="context-flow" markdown>

**核心问题**：如何将二次型的独立 $\chi^2$ 分解推广到矩阵情形？

</div>

!!! theorem "定理 72B.3 (矩阵 Cochran 定理)"
    设 $X \in \mathbb{R}^{n \times p}$ 的行独立同分布 $N_p(0, \Sigma)$。设 $A_1, A_2, \ldots, A_k$ 为 $n \times n$ 对称幂等矩阵（$A_i^2 = A_i = A_i^T$），秩分别为 $r_1, \ldots, r_k$。则以下条件等价：

    1. $A_i A_j = 0$ 对所有 $i \neq j$
    2. $\sum_{i=1}^k r_i = \text{rank}(\sum_{i=1}^k A_i)$
    3. $W_i = X^T A_i X \sim W_p(r_i, \Sigma)$ 相互独立

??? proof "证明"
    **(1) $\Rightarrow$ (3)**：设 $A_i A_j = 0$（$i \neq j$）。则 $A_i$ 的像空间相互正交。

    令 $P_i = A_i$（投影矩阵），$Y_i = P_i X$。由于 $P_i P_j = 0$，$Y_i$ 的行是 $X$ 的行在正交子空间上的投影，因此 $Y_i$ 和 $Y_j$ 独立。

    $W_i = X^T A_i X = X^T P_i X = Y_i^T Y_i$。由于 $P_i$ 秩为 $r_i$，$Y_i$ 可以视为 $r_i$ 个独立 $N_p(0, \Sigma)$ 向量的矩阵。故 $W_i \sim W_p(r_i, \Sigma)$。

    独立性：$W_i = Y_i^T Y_i$ 只依赖于 $Y_i$，而 $Y_i$ 相互独立，故 $W_i$ 相互独立。

    **(3) $\Rightarrow$ (1)**：若 $W_i$ 相互独立且各自为 Wishart，则利用 Wishart 分布的特征函数和独立性条件，可以推出投影矩阵必须正交。

    **(1) $\Leftrightarrow$ (2)**：$A_i A_j = 0$ 意味着 $\text{Im}(A_i) \perp \text{Im}(A_j)$，故 $\text{rank}(\sum A_i) = \sum \text{rank}(A_i) = \sum r_i$。

    $\blacksquare$

!!! example "例 72B.2"
    **方差分析的矩阵推广**。在单因素 MANOVA 中，总变异矩阵分解为：

    $$X^TX = X^T\left(\frac{1}{n}\mathbf{1}\mathbf{1}^T\right)X + X^T\left(P_G - \frac{1}{n}\mathbf{1}\mathbf{1}^T\right)X + X^T\left(I - P_G\right)X$$

    $$= W_{\text{mean}} + W_H + W_E$$

    其中 $P_G$ 是分组投影矩阵。矩阵 Cochran 定理保证 $W_H$ 和 $W_E$ 独立且各自为 Wishart。

---

## 72B.4 Wilks' Lambda 与 MANOVA

<div class="context-flow" markdown>

**核心问题**：如何检验多个多元正态总体的均值向量是否相同？

</div>

!!! definition "定义 72B.4 (Wilks' Lambda)"
    **Wilks' Lambda** 统计量定义为：

    $$\Lambda = \frac{|W_E|}{|W_E + W_H|} = \prod_{i=1}^{s}(1 - \beta_i) = \det(I - B)$$

    其中 $W_E$ 为组内（误差）平方和矩阵，$W_H$ 为组间（假设）平方和矩阵，$\beta_i$ 为矩阵 Beta 变量 $B = W_H(W_E + W_H)^{-1}$ 的特征值，$s = \min(p, n_H)$。

    $\Lambda \in (0, 1)$：$\Lambda$ 接近 0 表示组间差异大（拒绝 $H_0$），$\Lambda$ 接近 1 表示组间差异小（不拒绝）。

!!! theorem "定理 72B.4 (Wilks' Lambda 的精确分布)"
    在原假设下（各组均值相同），$\Lambda$ 的分布与矩阵 Beta 分布相关。在一些特殊情形下有精确的 $F$ 分布变换：

    - $p = 1$：$\frac{1-\Lambda}{\Lambda} \cdot \frac{n_E}{n_H} \sim F(n_H, n_E)$（即标准 ANOVA 的 $F$ 统计量）

    - $p = 2$：$\frac{1-\sqrt{\Lambda}}{\sqrt{\Lambda}} \cdot \frac{n_E - 1}{n_H} \sim F(2n_H, 2(n_E-1))$（Rao 的 $F$ 近似变为精确的）

    - $n_H = 1$（两组比较）：$\frac{1-\Lambda}{\Lambda} \cdot \frac{n_E - p + 1}{p} \sim F(p, n_E-p+1)$（即 Hotelling $T^2$ 检验）

    - $n_H = 2$：$\frac{1-\sqrt{\Lambda}}{\sqrt{\Lambda}} \cdot \frac{n_E - p + 1}{p} \sim F(2p, 2(n_E-p+1))$

    一般情形下，Rao 的 $F$ 近似为：

    $$F_{\text{approx}} = \frac{1 - \Lambda^{1/t}}{\Lambda^{1/t}} \cdot \frac{mt - 2u}{pn_H}$$

    其中 $t = \sqrt{\frac{p^2 n_H^2 - 4}{p^2 + n_H^2 - 5}}$，$m = n_E - (p - n_H + 1)/2$，$u = (pn_H - 2)/4$。

??? proof "证明"
    $p = 1$ 的情形直接验证：$\Lambda = W_E/(W_E + W_H)$，$W_E \sim \sigma^2\chi^2_{n_E}$，$W_H \sim \sigma^2\chi^2_{n_H}$（在 $H_0$ 下）且独立。

    $$\frac{1-\Lambda}{\Lambda} \cdot \frac{n_E}{n_H} = \frac{W_H}{W_E} \cdot \frac{n_E}{n_H} = \frac{W_H/n_H}{W_E/n_E} \sim F(n_H, n_E)$$

    $n_H = 1$ 的情形：$W_H = n\bar{\mathbf{d}}\bar{\mathbf{d}}^T$（秩 1），$\Lambda = |W_E|/|W_E + W_H|$。利用矩阵行列式引理 $|W_E + \mathbf{a}\mathbf{a}^T| = |W_E|(1 + \mathbf{a}^TW_E^{-1}\mathbf{a})$：

    $$\Lambda = \frac{1}{1 + n\bar{\mathbf{d}}^TW_E^{-1}\bar{\mathbf{d}}} = \frac{1}{1 + T^2/(n-1)}$$

    因此 $T^2 = (n-1)(1/\Lambda - 1)$，而 $T^2$ 的分布由 Hotelling 定理给出。

    $\blacksquare$

!!! definition "定义 72B.5 (四大 MANOVA 检验统计量)"
    设 $\theta_1 \geq \theta_2 \geq \cdots \geq \theta_s$ 为 $W_H W_E^{-1}$ 的非零特征值（$s = \min(p, n_H)$）。

    1. **Wilks' Lambda**：

    $$\Lambda = \prod_{i=1}^{s} \frac{1}{1 + \theta_i} = \frac{|W_E|}{|W_E + W_H|}$$

    2. **Pillai 迹**：

    $$V = \sum_{i=1}^{s} \frac{\theta_i}{1 + \theta_i} = \text{tr}(W_H(W_E + W_H)^{-1})$$

    3. **Hotelling-Lawley 迹**：

    $$U = \sum_{i=1}^{s} \theta_i = \text{tr}(W_H W_E^{-1})$$

    4. **Roy 最大根**：

    $$\Theta = \frac{\theta_1}{1 + \theta_1} = \max_{\mathbf{a}} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T(W_E + W_H)\mathbf{a}}$$

!!! theorem "定理 72B.5 (MANOVA 统计量的渐近分布)"
    在原假设 $H_0: \boldsymbol{\mu}_1 = \cdots = \boldsymbol{\mu}_g$ 下，当 $n_E \to \infty$：

    1. $-n_E \ln\Lambda \xrightarrow{d} \chi^2_{pn_H}$
    2. $n_E \cdot V \xrightarrow{d} \chi^2_{pn_H}$（当 $s = 1$）；一般情形需要更复杂的近似
    3. $n_E \cdot U / s \xrightarrow{d} \chi^2_{pn_H}/s$
    4. Roy 最大根没有简单的渐近 $\chi^2$ 形式，其临界值表需要单独查

    在实际应用中：

    - **Wilks' Lambda** 对所有备择假设方向有平衡的检验功效
    - **Pillai 迹** 对违背正态性和方差齐性最为稳健
    - **Hotelling-Lawley 迹** 在备择假设集中在一个方向时功效最高
    - **Roy 最大根** 在单方向备择假设下最优，但对多方向备择功效低

!!! theorem "定理 72B.6 (Roy 最大根的变分刻画)"
    $$\theta_1 = \max_{\mathbf{a} \neq 0} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T W_E \mathbf{a}}$$

    这是广义特征值问题 $W_H \mathbf{a} = \theta W_E \mathbf{a}$ 的最大特征值。等价地，$\theta_1 = \lambda_{\max}(W_E^{-1}W_H)$。

??? proof "证明"
    由 Rayleigh 商理论（Ch6），对称矩阵 $W_E^{-1/2}W_H W_E^{-1/2}$ 的最大特征值为：

    $$\lambda_{\max}(W_E^{-1/2}W_H W_E^{-1/2}) = \max_{\mathbf{b} \neq 0} \frac{\mathbf{b}^T W_E^{-1/2}W_H W_E^{-1/2}\mathbf{b}}{\mathbf{b}^T\mathbf{b}}$$

    令 $\mathbf{a} = W_E^{-1/2}\mathbf{b}$（可逆变换），$\mathbf{b} = W_E^{1/2}\mathbf{a}$：

    $$= \max_{\mathbf{a} \neq 0} \frac{\mathbf{a}^T W_E^{1/2} W_E^{-1/2}W_H W_E^{-1/2} W_E^{1/2}\mathbf{a}}{\mathbf{a}^T W_E \mathbf{a}} = \max_{\mathbf{a} \neq 0} \frac{\mathbf{a}^T W_H \mathbf{a}}{\mathbf{a}^T W_E \mathbf{a}}$$

    而 $W_E^{-1/2}W_H W_E^{-1/2}$ 与 $W_E^{-1}W_H = W_H W_E^{-1}$（在对称化意义下）有相同的特征值。

    最优解 $\mathbf{a}^*$ 是**判别方向**——在该方向上组间差异相对于组内变异最大。

    $\blacksquare$

!!! example "例 72B.3"
    **单因素 MANOVA**。三组样本，每组 20 个观测，$p = 2$ 维响应变量。

    $W_E = \begin{pmatrix} 100 & 30 \\ 30 & 80 \end{pmatrix}$，$W_H = \begin{pmatrix} 25 & 10 \\ 10 & 15 \end{pmatrix}$

    $\Lambda = \frac{|W_E|}{|W_E + W_H|} = \frac{100 \times 80 - 30^2}{125 \times 95 - 40^2} = \frac{7100}{10275} \approx 0.691$

    $$W_H W_E^{-1} = \frac{1}{7100}\begin{pmatrix} 25 & 10 \\ 10 & 15 \end{pmatrix}\begin{pmatrix} 80 & -30 \\ -30 & 100 \end{pmatrix} = \frac{1}{7100}\begin{pmatrix} 1700 & 250 \\ 350 & 1200 \end{pmatrix}$$

    特征值：$\text{tr} = 2900/7100 \approx 0.4085$，$\det = (1700 \times 1200 - 250 \times 350)/7100^2 \approx 0.03873$。

    四个统计量：Wilks' $\Lambda \approx 0.691$，Pillai $V \approx 0.297$，Hotelling-Lawley $U \approx 0.409$，Roy $\Theta = \theta_1/(1+\theta_1)$。

---

## 72B.5 Box M 检验

<div class="context-flow" markdown>

**核心问题**：MANOVA 假设各组的协方差矩阵相同（$\Sigma_1 = \Sigma_2 = \cdots = \Sigma_g$）。如何检验这个假设？

</div>

!!! definition "定义 72B.6 (Box M 统计量)"
    设 $g$ 组样本的组内协方差矩阵为 $S_1, \ldots, S_g$，各组样本量为 $n_1, \ldots, n_g$，合并协方差矩阵为 $S_{\text{pool}} = \sum_k (n_k-1)S_k / \sum_k(n_k-1)$。

    **Box M 统计量**为：

    $$M = \sum_{k=1}^{g}(n_k - 1)\ln|S_{\text{pool}}| - \sum_{k=1}^{g}(n_k-1)\ln|S_k|$$

    $$= \left(\sum_{k=1}^g (n_k-1)\right)\ln|S_{\text{pool}}| - \sum_{k=1}^g (n_k-1)\ln|S_k|$$

    在 $H_0: \Sigma_1 = \cdots = \Sigma_g$ 下，$M$ 渐近服从 $\chi^2$ 分布（带校正因子）。

!!! theorem "定理 72B.7 (Box M 的近似分布)"
    令 $N = \sum_{k=1}^g n_k$，$\nu = N - g$。Box 的校正因子为：

    $$c = \left(\sum_{k=1}^g \frac{1}{n_k-1} - \frac{1}{\nu}\right) \cdot \frac{2p^2 + 3p - 1}{6(p+1)(g-1)}$$

    在 $H_0$ 下：

    $$(1 - c) M \xrightarrow{d} \chi^2_f, \quad f = \frac{p(p+1)(g-1)}{2}$$

??? proof "证明"
    $M$ 是似然比统计量的 $-2\ln$ 倍的推广。对数似然比为：

    $$\Lambda_{LR} = \prod_{k=1}^g \left(\frac{|S_{\text{pool}}|}{|S_k|}\right)^{(n_k-1)/2}$$

    $$-2\ln\Lambda_{LR} = M$$

    由经典似然比理论，$-2\ln\Lambda_{LR}$ 在 $H_0$ 下渐近 $\chi^2_f$，自由度 $f$ 等于 $H_0$ 下的约束数：$g$ 个 $p \times p$ 对称矩阵的参数数（$gp(p+1)/2$）减去公共矩阵的参数数（$p(p+1)/2$），即 $f = (g-1)p(p+1)/2$。

    Box 的校正因子 $(1-c)$ 修正了有限样本偏差，使近似更准确。

    $\blacksquare$

!!! example "例 72B.4"
    三组（$g = 3$），$p = 2$，$n_1 = n_2 = n_3 = 20$。

    $$S_1 = \begin{pmatrix} 2.1 & 0.5 \\ 0.5 & 1.8 \end{pmatrix}, \quad S_2 = \begin{pmatrix} 1.9 & 0.3 \\ 0.3 & 2.2 \end{pmatrix}, \quad S_3 = \begin{pmatrix} 2.3 & 0.7 \\ 0.7 & 1.5 \end{pmatrix}$$

    $S_{\text{pool}} = (19S_1 + 19S_2 + 19S_3)/57$。

    $f = (3-1) \times 2 \times 3/2 = 6$。若 $(1-c)M < \chi^2_{6, 0.05} = 12.59$，则不拒绝方差齐性假设。

---

## 72B.6 Marchenko-Pastur 律

<div class="context-flow" markdown>

**核心问题**：当维度 $p$ 和样本量 $n$ 同时趋于无穷（$p/n \to \gamma > 0$）时，样本协方差矩阵的特征值分布趋向什么极限？

</div>

Marchenko-Pastur 律是高维统计学和随机矩阵理论中最基本的结果，揭示了样本协方差矩阵在高维情形下的根本不同行为。

!!! theorem "定理 72B.8 (Marchenko-Pastur 律)"
    设 $X \in \mathbb{R}^{n \times p}$ 的元素独立同分布，均值为 0，方差为 $\sigma^2/n$。令 $W = X^TX$（$p \times p$ 样本协方差矩阵）。当 $p, n \to \infty$ 且 $p/n \to \gamma \in (0, \infty)$ 时，$W$ 的经验谱分布（特征值的经验分布）几乎必然弱收敛到 **Marchenko-Pastur 分布** $F_\gamma$：

    $$dF_\gamma(\lambda) = \left(1 - \frac{1}{\gamma}\right)^+ \delta_0(\lambda) + \frac{\sqrt{(\lambda_+ - \lambda)(\lambda - \lambda_-)}}{2\pi\gamma\sigma^2\lambda} d\lambda$$

    其中 $\lambda_\pm = \sigma^2(1 \pm \sqrt{\gamma})^2$ 为谱的上下边界，$(x)^+ = \max(x, 0)$。

    关键特征：

    - **支撑**：$[\lambda_-, \lambda_+]$（当 $\gamma \leq 1$ 时）或 $\{0\} \cup [\lambda_-, \lambda_+]$（当 $\gamma > 1$ 时）
    - **相变**：$\gamma = 1$ 时，$\lambda_- = 0$（谱的下边界触及 0）
    - **当 $\gamma > 1$**：有 $p - n$ 个零特征值（$W$ 奇异），占比 $(1-1/\gamma)$

??? proof "证明"
    **证明思路**（Stieltjes 变换方法）：

    经验谱分布 $\hat{F}_p(\lambda) = \frac{1}{p}\sum_{i=1}^p \mathbf{1}(\lambda_i \leq \lambda)$ 的 Stieltjes 变换为：

    $$m_p(z) = \int \frac{1}{\lambda - z} d\hat{F}_p(\lambda) = \frac{1}{p}\text{tr}(W - zI)^{-1}, \quad z \in \mathbb{C}^+$$

    **第一步**：建立 $m(z)$（$m_p$ 的极限）的自洽方程。

    对 $W = X^TX = \frac{1}{n}\sum_{i=1}^n \mathbf{x}_i\mathbf{x}_i^T$（$\mathbf{x}_i$ 为 $X$ 的行），利用 Sherman-Morrison 公式逐一去除 $\mathbf{x}_i$ 的贡献，并在 $p, n \to \infty$ 的极限下取期望，得到：

    $$m(z) = \frac{1}{-z + \gamma\sigma^2/(1 + \sigma^2 m(z) \cdot (-z)^{-1} \cdot ...)}$$

    更精确地，极限 Stieltjes 变换满足**自洽方程**：

    $$m(z) = \frac{1}{\sigma^2(1 - \gamma) - z + \gamma\sigma^2 z m(z) \cdot ...}$$

    最终化简为关于 $m(z)$ 的二次方程：

    $$\gamma\sigma^2 z m(z)^2 - (\sigma^2(\gamma - 1) + z)m(z) - 1 = 0$$

    **第二步**：解此二次方程，选取满足 $\text{Im}(m(z)) > 0$（$z \in \mathbb{C}^+$）的根：

    $$m(z) = \frac{-(\sigma^2(\gamma-1)+z) + \sqrt{(\sigma^2(\gamma-1)+z)^2 + 4\gamma\sigma^2 z}}{2\gamma\sigma^2 z}$$

    **第三步**：由 Stieltjes 变换的反演公式：

    $$dF_\gamma(\lambda) = \lim_{\epsilon \to 0^+} \frac{1}{\pi} \text{Im}[m(\lambda + i\epsilon)] d\lambda$$

    代入 $m(z)$ 的表达式，当 $z = \lambda + i\epsilon$ 趋向实轴时，$\text{Im}(m)$ 非零当且仅当判别式的平方根为纯虚数，即当 $\lambda \in [\lambda_-, \lambda_+]$ 时。计算得到 Marchenko-Pastur 密度。

    **第四步**：当 $\gamma > 1$ 时，$\text{rank}(W) = \min(n, p) = n < p$，故有 $p - n$ 个零特征值，占比 $1 - n/p \to 1 - 1/\gamma$，对应 $\delta_0$ 的质量。

    $\blacksquare$

!!! example "例 72B.5"
    **特征值的散布**。设真实协方差为 $\Sigma = I_p$（所有变量独立同方差），$\sigma^2 = 1$。

    - **低维情形**（$p = 10, n = 1000, \gamma = 0.01$）：$\lambda_\pm = (1 \pm 0.1)^2$，特征值集中在 $[0.81, 1.21]$——接近真值 1。

    - **高维情形**（$p = 500, n = 1000, \gamma = 0.5$）：$\lambda_\pm = (1 \pm \sqrt{0.5})^2 \approx [0.086, 2.914]$——样本特征值散布极大，最大特征值约为真值的 3 倍！

    - **临界情形**（$p = n = 1000, \gamma = 1$）：$\lambda_\pm = (1 \pm 1)^2 = [0, 4]$——样本协方差矩阵奇异，特征值从 0 到 4。

    这解释了为什么高维协方差估计如此困难：即使真实协方差矩阵是单位阵，样本特征值也会严重偏离。

---

## 72B.7 Tracy-Widom 分布

<div class="context-flow" markdown>

**核心问题**：Marchenko-Pastur 律描述了特征值的整体分布。最大特征值的波动遵循什么分布？

</div>

!!! definition "定义 72B.7 (Tracy-Widom 分布)"
    **Tracy-Widom 分布** $\text{TW}_\beta$（$\beta = 1, 2, 4$ 分别对应实、复、四元数情形）描述了大随机矩阵最大特征值围绕其极限值的波动。

    对 $\beta = 1$（实情形），$\text{TW}_1$ 由 Painlev\'e II 方程定义：

    $$F_1(s) = \exp\left(-\frac{1}{2}\int_s^\infty \left[q(x) + (x-s)q^2(x)\right] dx\right)$$

    其中 $q(x)$ 是 Painlev\'e II 方程 $q''(x) = xq(x) + 2q^3(x)$ 满足 $q(x) \sim \text{Ai}(x)$（$x \to \infty$）的解。

    $\text{TW}_1$ 分布是**非对称**的，左偏（偏度约 $-0.29$），均值约 $-1.21$，方差约 $1.27$。

!!! theorem "定理 72B.9 (最大特征值的 Tracy-Widom 极限)"
    设 $W \sim W_p(n, I_p)$，$\lambda_1$ 为 $W$ 的最大特征值。当 $p, n \to \infty$，$p/n \to \gamma \in (0, 1]$ 时：

    $$\frac{\lambda_1 - \mu_{n,p}}{\sigma_{n,p}} \xrightarrow{d} \text{TW}_1$$

    其中中心化和标度参数为：

    $$\mu_{n,p} = (\sqrt{n} + \sqrt{p})^2, \quad \sigma_{n,p} = (\sqrt{n} + \sqrt{p})\left(\frac{1}{\sqrt{n}} + \frac{1}{\sqrt{p}}\right)^{1/3}$$

    注意 $\mu_{n,p}/n = (1 + \sqrt{\gamma})^2 = \lambda_+$（Marchenko-Pastur 上界），符合直觉。

!!! example "例 72B.6"
    **显著性检验**。$p = 100, n = 500$。零假设 $\Sigma = I_{100}$。

    $\mu = (\sqrt{500} + \sqrt{100})^2 = (22.36 + 10)^2 \approx 1047.5$

    $\sigma = (22.36 + 10)(1/22.36 + 1/10)^{1/3} \approx 32.36 \times 0.528 \approx 17.1$

    若观察到 $\lambda_1 = 1120$，标准化值 $(1120 - 1047.5)/17.1 \approx 4.24$。

    $\text{TW}_1$ 的 5% 右尾临界值约 $0.98$，$1\%$ 约 $2.02$。$4.24$ 远超临界值，拒绝 $H_0$——存在超出噪声水平的信号特征值。

---

## 72B.8 Spiked 协方差模型与 BBP 相变

<div class="context-flow" markdown>

**核心问题**：当真实协方差矩阵有少量大特征值（信号）叠加在噪声上时，样本特征值能否检测到信号？存在什么相变现象？

</div>

!!! definition "定义 72B.8 (Spiked 协方差模型)"
    **Spiked 协方差模型**（Johnstone, 2001）假设真实协方差矩阵为：

    $$\Sigma = I_p + \sum_{k=1}^{r} \ell_k \mathbf{u}_k \mathbf{u}_k^T$$

    即 $\Sigma$ 有 $r$ 个"spiked"特征值 $1 + \ell_1 \geq \cdots \geq 1 + \ell_r > 1$（信号），其余 $p - r$ 个特征值为 1（噪声）。$r$ 固定，$p, n \to \infty$。

!!! theorem "定理 72B.10 (BBP 相变——Baik-Ben Arous-P\'ech\'e)"
    在 spiked 协方差模型中，$p/n \to \gamma$，第 $k$ 个样本特征值 $\hat{\lambda}_k$ 的行为取决于信号强度 $\ell_k$：

    **相变阈值**：$\ell_{\text{crit}} = \sqrt{\gamma}$。

    1. **超临界情形**（$\ell_k > \sqrt{\gamma}$）：$\hat{\lambda}_k$ 几乎必然收敛到

    $$\hat{\lambda}_k \to (1 + \ell_k)\left(1 + \frac{\gamma}{\ell_k}\right) > \lambda_+ = (1 + \sqrt{\gamma})^2$$

    即样本特征值脱离 Marchenko-Pastur 谱的上边界，可以检测到信号。

    2. **亚临界情形**（$\ell_k < \sqrt{\gamma}$）：$\hat{\lambda}_k \to \lambda_+ = (1 + \sqrt{\gamma})^2$

    样本特征值"黏附"在 MP 上边界，**信号不可检测**。

    3. **临界情形**（$\ell_k = \sqrt{\gamma}$）：$\hat{\lambda}_k$ 围绕 $\lambda_+$ 以 Tracy-Widom 标度波动。

??? proof "证明"
    **证明思路**（基于 Stieltjes 变换的相变分析）：

    设 $\Sigma = I_p + \ell \mathbf{u}\mathbf{u}^T$（单个 spike），$\ell > 0$。样本协方差矩阵 $\hat{\Sigma} = \frac{1}{n}X^TX$，$X$ 的行独立 $\sim N_p(0, \Sigma)$。

    **关键思想**：$\hat{\Sigma}$ 可以分解为噪声部分和信号部分的扰动。利用 Sherman-Morrison 型论证：

    令 $\hat{\Sigma}_0 = \frac{1}{n}\tilde{X}^T\tilde{X}$（$\Sigma = I_p$ 时的样本协方差），则 $\hat{\Sigma}$ 与 $\hat{\Sigma}_0$ 的关系通过 $\Sigma^{1/2} = I + (\sqrt{1+\ell}-1)\mathbf{u}\mathbf{u}^T$ 建立。

    最大特征值 $\hat{\lambda}_1$ 是以下方程的解：

    $$\det(\hat{\Sigma} - \lambda I) = 0$$

    利用矩阵行列式引理和 Stieltjes 变换的渐近展开，可以证明 $\hat{\lambda}_1$ 满足：

    $$\frac{1}{\ell} = -m_{\gamma}(\hat{\lambda}_1)$$

    其中 $m_\gamma$ 是 Marchenko-Pastur 律的 Stieltjes 变换在实轴上 $\lambda > \lambda_+$ 的值。

    $m_\gamma$ 在 $(\lambda_+, \infty)$ 上从 $-\infty$ 单调递增到 $0$，且 $m_\gamma(\lambda_+) = -1/\sqrt{\gamma}$（右极限）。

    因此方程 $1/\ell = -m_\gamma(\lambda)$ 在 $(\lambda_+, \infty)$ 上有解，当且仅当 $1/\ell < 1/\sqrt{\gamma}$，即 $\ell > \sqrt{\gamma}$。

    当 $\ell \leq \sqrt{\gamma}$ 时无解——$\hat{\lambda}_1$ 无法脱离 Marchenko-Pastur 上边界。

    $\blacksquare$

!!! theorem "定理 72B.11 (信号特征值的偏差与特征向量的不一致性)"
    在超临界情形 $\ell > \sqrt{\gamma}$ 下：

    1. **特征值的渐近偏差**：$\hat{\lambda} / (1+\ell) \to 1 + \gamma/\ell > 1$，即样本特征值**高估**真实特征值。

    2. **特征向量的不一致性**：样本特征向量 $\hat{\mathbf{u}}$ 与真实特征向量 $\mathbf{u}$ 之间的内积平方满足：

    $$|\langle \hat{\mathbf{u}}, \mathbf{u} \rangle|^2 \to 1 - \frac{\gamma}{\ell^2} < 1$$

    即使信号可检测，样本特征向量也不一致地估计真实方向——**相关系数严格小于 1**。当 $\ell \to \sqrt{\gamma}^+$ 时，$|\langle \hat{\mathbf{u}}, \mathbf{u}\rangle|^2 \to 0$——特征向量完全失效。

!!! example "例 72B.7"
    **基因组学中的应用**。$p = 10000$ 个基因，$n = 500$ 个样本，$\gamma = 20$。

    BBP 阈值 $\ell_{\text{crit}} = \sqrt{20} \approx 4.47$。

    - 若真实信号强度 $\ell = 10 > 4.47$：信号可检测，但样本特征值 $(1+10)(1+20/10) = 33$ 远高于真值 $11$，且 $|\langle \hat{\mathbf{u}}, \mathbf{u}\rangle|^2 = 1 - 20/100 = 0.8$。

    - 若 $\ell = 3 < 4.47$：信号不可检测，样本最大特征值约为 $(1+\sqrt{20})^2 \approx 25.9$，完全由噪声主导。

    这对 PCA 在高维基因组数据中的应用有深刻影响：只有足够强的种群结构信号才能被 PCA 检测到。

---

## 72B.9 Fisher 信息度量与 Wishart 流形

<div class="context-flow" markdown>

**核心问题**：Wishart 分布族构成的统计流形具有什么几何结构？

</div>

!!! theorem "定理 72B.12 (Wishart 分布的 Fisher 信息)"
    Wishart 分布族 $\{W_p(n, \Sigma) : \Sigma > 0\}$ 构成一个统计流形，其 Fisher 信息矩阵为：

    $$g_{ij,kl} = \frac{n}{2}\left[(\Sigma^{-1})_{ik}(\Sigma^{-1})_{jl} + (\Sigma^{-1})_{il}(\Sigma^{-1})_{jk}\right]$$

    对应的测地距离（Rao 距离）为：

    $$d(\Sigma_1, \Sigma_2) = \sqrt{\frac{n}{2}} \left\|\log(\Sigma_1^{-1/2}\Sigma_2\Sigma_1^{-1/2})\right\|_F$$

    即与正定矩阵流形上的仿射不变 Riemannian 距离成比例。

??? proof "证明"
    Wishart 分布的对数似然（作为 $\Sigma$ 的函数）为：

    $$\ell(\Sigma) = -\frac{n}{2}\log|\Sigma| - \frac{1}{2}\text{tr}(\Sigma^{-1}W) + \text{const}$$

    Fisher 信息矩阵的 $(ij, kl)$ 元素为：

    $$I(\Sigma)_{ij,kl} = -E\left[\frac{\partial^2 \ell}{\partial \Sigma_{ij}\partial \Sigma_{kl}}\right]$$

    利用矩阵微分的标准结果：

    - $\frac{\partial}{\partial \Sigma_{ij}}\log|\Sigma| = (\Sigma^{-1})_{ji}$（对对称矩阵，$= (\Sigma^{-1})_{ij}$）
    - $\frac{\partial}{\partial \Sigma_{ij}}\text{tr}(\Sigma^{-1}W) = -(\Sigma^{-1}W\Sigma^{-1})_{ji}$
    - $\frac{\partial^2 \ell}{\partial \Sigma_{ij}\partial \Sigma_{kl}}$ 涉及 $\Sigma^{-1}$ 的二阶导数

    经过计算并取期望（利用 $E[W] = n\Sigma$），得到上述 Fisher 信息。

    测地距离的推导利用了正定矩阵流形 $\mathcal{P}_p = \text{GL}(p)/O(p)$ 是齐次空间的事实。Fisher 信息度量在此空间上诱导的距离与仿射不变 Riemannian 距离 $d(\Sigma_1, \Sigma_2) = \|\log(\Sigma_1^{-1/2}\Sigma_2\Sigma_1^{-1/2})\|_F$ 成比例。

    $\blacksquare$

---

## 72B.10 应用

<div class="context-flow" markdown>

**核心问题**：矩阵值分布和多元推断在实际问题中如何应用？

</div>

!!! example "例 72B.8 (金融：投资组合理论)"
    在 Markowitz 投资组合优化中，$p$ 种资产的收益率协方差矩阵 $\Sigma$ 的估计至关重要。经典的样本协方差矩阵在 $p$ 接近 $n$ 时表现不佳——Marchenko-Pastur 律告诉我们样本特征值会严重偏离真实值。

    **收缩估计**（Ledoit-Wolf）：$\hat{\Sigma}_{\text{shrink}} = \alpha S + (1-\alpha)\text{tr}(S)/p \cdot I_p$，其中 $\alpha$ 的最优选择由 Marchenko-Pastur 律的矩决定。

    **BBP 相变的应用**：在 $\gamma = p/n \approx 0.5$ 的金融数据中，只有前几个主成分（对应行业因子）信号强度超过 $\sqrt{\gamma}$ 阈值，其余主成分不可靠。

!!! example "例 72B.9 (脑成像：扩散张量 MRI)"
    在扩散张量 MRI（DTI）中，每个体素对应一个 $3 \times 3$ 对称正定矩阵 $D$——扩散张量。

    扩散张量的统计建模自然涉及矩阵值分布。Wishart 分布用于估计扩散张量的不确定性，其特征值（本征值）给出三个主扩散方向的扩散系数：

    - $\lambda_1 \gg \lambda_2 \approx \lambda_3$：各向异性扩散（白质纤维束）
    - $\lambda_1 \approx \lambda_2 \approx \lambda_3$：各向同性扩散（灰质或脑脊液）

    **分数各向异性**（FA）：

    $$\text{FA} = \sqrt{\frac{3}{2}} \cdot \frac{\sqrt{(\lambda_1 - \bar{\lambda})^2 + (\lambda_2 - \bar{\lambda})^2 + (\lambda_3 - \bar{\lambda})^2}}{\sqrt{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}$$

    Fisher 信息度量在 DTI 中用于定义张量之间的"距离"，比 Frobenius 范数更自然（具有仿射不变性）。

!!! example "例 72B.10 (Wishart 过程)"
    协方差矩阵的时变性在金融中非常重要（波动率聚集现象）。**Wishart 过程** $\{W_t\}$ 是 Wishart 分布到连续时间的推广，满足矩阵值 SDE：

    $$dW_t = (\nu Q^T Q + M W_t + W_t M^T) dt + \sqrt{W_t} dB_t Q + Q^T dB_t^T \sqrt{W_t}$$

    其中 $B_t$ 是矩阵值 Brown 运动，$M$ 为均值回复矩阵，$Q$ 为波动率矩阵，$\nu$ 为自由度。

    关键性质：$W_t > 0$（正定性自动保持，当 $\nu$ 足够大时），边际分布为（非中心）Wishart，可以模拟协方差矩阵的动态变化。

---

## 习题

!!! exercise "习题 72B.1"
    证明当 $p = 1$ 时，矩阵 Beta 分布退化为标量 Beta 分布：若 $W_1 \sim \chi^2_{n_1}$ 和 $W_2 \sim \chi^2_{n_2}$ 独立，则 $B = W_1/(W_1+W_2) \sim \text{Beta}(n_1/2, n_2/2)$。

!!! exercise "习题 72B.2"
    对 $n = 25$，$p = 4$ 的 Hotelling $T^2$ 检验：

    (a) $T^2$ 与 $F$ 统计量之间的精确关系是什么？

    (b) $H_0$ 下 $F$ 统计量的自由度是什么？

    (c) 若 $T^2 = 22.5$，是否在 5% 水平上拒绝 $H_0$？

!!! exercise "习题 72B.3"
    证明矩阵 Cochran 定理中条件 (1) $\Leftrightarrow$ (2)。即：对称幂等矩阵 $A_1, \ldots, A_k$ 满足 $A_iA_j = 0$（$i \neq j$）当且仅当 $\sum_i \text{rank}(A_i) = \text{rank}(\sum_i A_i)$。

    提示：利用幂等矩阵的像空间分解。

!!! exercise "习题 72B.4"
    在例 72B.3 中，完成以下计算：

    (a) 求 $W_HW_E^{-1}$ 的两个特征值 $\theta_1, \theta_2$。

    (b) 计算 Roy 最大根 $\Theta = \theta_1/(1+\theta_1)$。

    (c) 利用 $p = 2, n_H = 2$ 的精确 $F$ 变换，检验 $H_0$。

!!! exercise "习题 72B.5"
    **Marchenko-Pastur 律的矩**。证明 Marchenko-Pastur 分布的 $k$ 阶矩为 Catalan 数的推广：

    $$\int \lambda^k dF_\gamma(\lambda) = \sigma^{2k} \sum_{j=0}^{k-1} \frac{1}{j+1}\binom{k}{j}\binom{k-1}{j} \gamma^{j+1}$$

    特别地，一阶矩（均值）为 $\sigma^2$，二阶矩为 $\sigma^4(1 + \gamma)$。

!!! exercise "习题 72B.6"
    对 spiked 协方差模型 $\Sigma = I_p + \ell \mathbf{u}\mathbf{u}^T$：

    (a) 当 $\gamma = 1/4$（$n = 4p$）时，BBP 阈值 $\ell_{\text{crit}}$ 是多少？

    (b) 若 $\ell = 2$，计算渐近样本最大特征值 $\hat{\lambda}_1$ 和特征向量一致性 $|\langle\hat{\mathbf{u}}, \mathbf{u}\rangle|^2$。

    (c) 讨论：在什么意义下 PCA 在高维中"失效"？提出可能的修正方法。

!!! exercise "习题 72B.7"
    **Tracy-Widom 检验**。$p = 50, n = 200, \gamma = 0.25$。零假设 $\Sigma = I_{50}$。

    (a) 计算 $\mu_{n,p}$ 和 $\sigma_{n,p}$。

    (b) 若观察到最大特征值 $\lambda_1 = 2.8$，计算标准化值。

    (c) 若 $\text{TW}_1$ 的 95% 分位点约为 $0.98$，是否拒绝 $H_0$？

!!! exercise "习题 72B.8"
    证明 Box M 统计量可以写成 Kullback-Leibler 散度的加权和：

    $$M = 2\sum_{k=1}^g (n_k-1) D_{KL}(N(0, S_k) \| N(0, S_{\text{pool}}))$$

    其中 $D_{KL}(N(0, A) \| N(0, B)) = \frac{1}{2}[\text{tr}(B^{-1}A) - p - \ln|AB^{-1}|]$。

!!! exercise "习题 72B.9"
    **Marchenko-Pastur 律的模拟验证**。

    (a) 生成 $p = 200, n = 500$ 的随机矩阵 $X$（元素独立标准正态），计算 $W = X^TX/n$ 的特征值。

    (b) 画出特征值直方图并与 $\gamma = 0.4$ 的 Marchenko-Pastur 密度比较。

    (c) 验证最大特征值接近 $\lambda_+ = (1 + \sqrt{0.4})^2 \approx 2.97$。

!!! exercise "习题 72B.10"
    **Fisher 信息度量的不变性**。证明 Wishart 流形上的 Rao 距离 $d(\Sigma_1, \Sigma_2) = \sqrt{n/2}\|\log(\Sigma_1^{-1/2}\Sigma_2\Sigma_1^{-1/2})\|_F$ 满足仿射不变性：

    $$d(A\Sigma_1A^T, A\Sigma_2A^T) = d(\Sigma_1, \Sigma_2)$$

    对任意可逆矩阵 $A$。这为什么在 DTI 等应用中是重要的性质？
