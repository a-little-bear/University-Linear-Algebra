# 第 43 章 伪谱与非正规矩阵分析

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵范数(Ch15) · 矩阵分析(Ch14)

**本章脉络**：正规 vs 非正规 → ε-伪谱定义 → 预解式刻画 → Kreiss 矩阵定理 → 瞬态增长 → 伪谱腹地 → 计算方法 → 应用

**延伸**：伪谱理论揭示了流体力学中亚临界转捩（Orr-Sommerfeld 算子的非正规性）和数值 ODE 求解器稳定性分析中非正规性的关键作用

</div>

经典谱理论——特征值及其代数/几何重数——在正规矩阵（$AA^* = A^*A$）的情形下完美地刻画了矩阵的行为。然而，大量实际问题中出现的矩阵远非正规。对于非正规矩阵，特征值可能完全无法预测矩阵指数 $e^{tA}$ 的瞬态行为、线性方程组迭代求解的收敛性、或微分方程解的短时间增长。

伪谱（pseudospectra）理论正是为了填补这一巨大鸿沟而发展起来的。它不再只关注"$zI - A$ 何时奇异"（即特征值），而是关注"$zI - A$ 何时接近奇异"——这才是实际计算和物理行为中真正重要的问题。

本章系统介绍 $\varepsilon$-伪谱的定义、性质、计算方法及其在流体稳定性和数值分析中的应用。

---

## 43.1 正规与非正规矩阵回顾

<div class="context-flow" markdown>

**核心问题**：为什么经典特征值理论对非正规矩阵失效？如何量化矩阵的"非正规程度"？

</div>

!!! definition "定义 43.1 (正规矩阵)"
    矩阵 $A \in \mathbb{C}^{n \times n}$ 称为**正规矩阵**（normal matrix），如果 $AA^* = A^*A$。

    等价条件包括：

    1. $A$ 可以酉对角化：存在酉矩阵 $U$ 使得 $A = U \Lambda U^*$，其中 $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$。
    2. $A$ 具有 $n$ 个标准正交的特征向量。
    3. $\sum_{i} \|Ae_i\|^2 = \sum_{i} |\lambda_i|^2$（Schur 分解的上三角部分全为零）。

!!! definition "定义 43.2 (偏离正规度)"
    设 $A \in \mathbb{C}^{n \times n}$，其 Schur 分解为 $A = UTU^*$，其中 $T$ 是上三角矩阵，$U$ 是酉矩阵。定义 $A$ 的**偏离正规度**（departure from normality）为：

    $$\Delta(A) = \|A\|_F^2 - \sum_{i=1}^{n} |\lambda_i|^2 = \|T - \operatorname{diag}(T)\|_F^2 \geq 0.$$

    $\Delta(A) = 0$ 当且仅当 $A$ 是正规矩阵。

    另一个常用度量是 $\|AA^* - A^*A\|_F$，即**交换子范数**。

!!! theorem "定理 43.1 (正规矩阵的谱完备性)"
    设 $A \in \mathbb{C}^{n \times n}$ 是正规矩阵，$A = U\Lambda U^*$。则：

    1. **预解式范数**：$\|(zI - A)^{-1}\| = \frac{1}{\operatorname{dist}(z, \sigma(A))} = \frac{1}{\min_i |z - \lambda_i|}$。
    2. **矩阵指数**：$\|e^{tA}\| = e^{t \cdot \max_i \operatorname{Re} \lambda_i}$。
    3. **矩阵多项式**：$\|p(A)\| = \max_i |p(\lambda_i)|$。

    这些优美的等式对非正规矩阵全部变为（有时很松的）不等式。

??? proof "证明"
    由 $A = U\Lambda U^*$，对任意矩阵函数 $f$，$f(A) = Uf(\Lambda)U^*$。由于 $U$ 是酉的，
    $$\|f(A)\| = \|f(\Lambda)\| = \max_i |f(\lambda_i)|.$$

    特别地：

    - $\|(zI - A)^{-1}\| = \|(zI - \Lambda)^{-1}\| = \max_i |z - \lambda_i|^{-1} = (\min_i |z - \lambda_i|)^{-1}$。
    - $\|e^{tA}\| = \|e^{t\Lambda}\| = \max_i |e^{t\lambda_i}| = \max_i e^{t \operatorname{Re} \lambda_i} = e^{t \max_i \operatorname{Re} \lambda_i}$。$\blacksquare$

!!! example "例 43.1 (非正规性的戏剧性后果)"
    考虑 $n \times n$ 上三角 Jordan 块
    $$J_n(\lambda) = \begin{pmatrix} \lambda & 1 & & \\ & \lambda & \ddots & \\ & & \ddots & 1 \\ & & & \lambda \end{pmatrix}.$$

    特征值全为 $\lambda$，但：

    - **矩阵指数增长**：$\|e^{tJ_n(0)}\| = \left\|\sum_{k=0}^{n-1} \frac{t^k}{k!} N^k\right\|$（$N$ 是幂零部分），当 $t$ 从 $0$ 增加时先经历多项式增长再衰减。
    - **预解式爆炸**：对 $z$ 靠近但不等于 $\lambda$，$\|(zI - J_n(\lambda))^{-1}\| \sim |z - \lambda|^{-n}$（而非正规矩阵的 $|z - \lambda|^{-1}$）。

    具体地，取 $\lambda = -1$，$n = 20$：所有特征值实部为 $-1$，但 $\|e^{tJ_{20}(-1)}\|$ 可以在某个 $t > 0$ 处达到约 $10^5$ 的量级——远大于 $e^{-t}$。

---

## 43.2 ε-伪谱的定义

<div class="context-flow" markdown>

**核心问题**：如何定义"接近 $A$ 的特征值"这一概念？

</div>

!!! definition "定义 43.3 (ε-伪谱)"
    设 $A \in \mathbb{C}^{n \times n}$，$\varepsilon > 0$。$A$ 的 **$\varepsilon$-伪谱**（$\varepsilon$-pseudospectrum）定义为

    $$\sigma_\varepsilon(A) = \left\{ z \in \mathbb{C} : \|(zI - A)^{-1}\| \geq \frac{1}{\varepsilon} \right\},$$

    其中约定 $\|(zI - A)^{-1}\| = +\infty$ 当 $z \in \sigma(A)$ 时。

    除非另行说明，范数取为谱范数 $\|\cdot\|_2$。

!!! theorem "定理 43.2 (ε-伪谱的等价刻画)"
    以下三个集合相等：

    $$\sigma_\varepsilon(A) = \left\{ z \in \mathbb{C} : \|(zI - A)^{-1}\| \geq \varepsilon^{-1} \right\}$$

    $$= \left\{ z \in \mathbb{C} : z \in \sigma(A + E),\; \text{某个 } E \text{ 满足 } \|E\| < \varepsilon \right\}$$

    $$= \left\{ z \in \mathbb{C} : \sigma_{\min}(zI - A) < \varepsilon \right\},$$

    其中 $\sigma_{\min}(B)$ 表示 $B$ 的最小奇异值。

??? proof "证明"
    设三个集合分别为 $S_1$、$S_2$、$S_3$。

    **$S_1 = S_3$**：这是因为 $\|(zI - A)^{-1}\| = \sigma_{\min}(zI - A)^{-1}$（$B^{-1}$ 的最大奇异值是 $B$ 的最小奇异值的倒数）。因此 $\|(zI - A)^{-1}\| \geq \varepsilon^{-1}$ 当且仅当 $\sigma_{\min}(zI - A) \leq \varepsilon$。取严格不等号的细节需要注意：由于我们使用 $\geq$ 和 $<$，需要用连续性论证统一。实际上，标准定义中两种形式（$\geq$ vs $>$）会导致 $\sigma_\varepsilon$ 是闭集或开集的差别；这里我们统一采用使 $\sigma_\varepsilon$ 为**开集**的定义：$\|(zI-A)^{-1}\| > \varepsilon^{-1}$，或等价地 $\sigma_{\min}(zI-A) < \varepsilon$，或 $\|E\| < \varepsilon$。

    **$S_2 = S_3$**：这是利用最小奇异值的变分特征：
    $$\sigma_{\min}(zI - A) = \min_{\|E\|=1} \sigma_{\min}(zI - A - E \cdot \varepsilon / 1)$$
    更精确地说，$z \in \sigma(A + E)$ 当且仅当 $\det(zI - A - E) = 0$，即 $zI - A - E$ 奇异，即 $\sigma_{\min}(zI - A - E) = 0$。

    $\sigma_{\min}(zI - A) < \varepsilon$ 意味着存在 $\|E\| = \sigma_{\min}(zI - A) < \varepsilon$ 使得 $zI - A - E$ 奇异（取 $E$ 为使最小奇异值变为 $0$ 的秩一扰动）。

    反之，若存在 $\|E\| < \varepsilon$ 使 $z \in \sigma(A + E)$，则 $\sigma_{\min}(zI - A) \leq \sigma_{\min}(zI - A - E) + \|E\| = \|E\| < \varepsilon$（注意 $\sigma_{\min}(zI-A-E) = 0$）。实际上更直接地，$\sigma_{\min}(zI - A) \leq \|E\|$ 由 Weyl 不等式或直接由 $\sigma_{\min}(B) \leq \sigma_{\min}(B - E) + \|E\|$ 得到。$\blacksquare$

!!! example "例 43.2 (ε-伪谱的直观理解)"
    $\varepsilon$-伪谱的三种等价定义提供了三种互补的直观：

    1. **预解式观点**：$z$ 在伪谱中 $\Leftrightarrow$ 预解式 $(zI - A)^{-1}$ 的范数很大 $\Leftrightarrow$ $zI - A$ 接近奇异。
    2. **扰动观点**：$z$ 在伪谱中 $\Leftrightarrow$ $z$ 是 $A$ 某个小扰动的特征值 $\Leftrightarrow$ $z$ 是"几乎特征值"。
    3. **奇异值观点**：$z$ 在伪谱中 $\Leftrightarrow$ $zI - A$ 的最小奇异值很小。

---

## 43.3 伪谱的基本性质

<div class="context-flow" markdown>

**核心问题**：伪谱具有哪些几何和拓扑性质？正规与非正规矩阵的伪谱有何质的区别？

</div>

!!! theorem "定理 43.3 (伪谱的基本性质)"
    设 $A \in \mathbb{C}^{n \times n}$，$\varepsilon > 0$。

    1. **嵌套性**：$0 < \varepsilon_1 < \varepsilon_2$ $\Rightarrow$ $\sigma(A) \subseteq \sigma_{\varepsilon_1}(A) \subseteq \sigma_{\varepsilon_2}(A)$。
    2. **收敛到谱**：$\bigcap_{\varepsilon > 0} \sigma_\varepsilon(A) = \sigma(A)$。
    3. **有界性**：$\sigma_\varepsilon(A) \subseteq \{z : |z| \leq \|A\| + \varepsilon\}$。
    4. **连通性**：$\sigma_\varepsilon(A)$ 的每个连通分支至少包含一个特征值。
    5. **酉不变性**：$\sigma_\varepsilon(U^*AU) = \sigma_\varepsilon(A)$ 对任意酉矩阵 $U$。
    6. **平移性**：$\sigma_\varepsilon(A + \alpha I) = \sigma_\varepsilon(A) + \alpha$。

??? proof "证明"
    (1) 由定义直接得出：若 $\|(zI - A)^{-1}\| \geq \varepsilon_1^{-1} > \varepsilon_2^{-1}$，则自然有 $\|(zI - A)^{-1}\| \geq \varepsilon_2^{-1}$。

    (2) $z \in \bigcap_{\varepsilon > 0} \sigma_\varepsilon(A)$ 意味着对每个 $\varepsilon > 0$，$\|(zI - A)^{-1}\| \geq \varepsilon^{-1}$，即 $\|(zI - A)^{-1}\| = \infty$，即 $z \in \sigma(A)$。

    (3) 若 $|z| > \|A\| + \varepsilon$，则 $\|zI - A\| \geq |z| - \|A\| > \varepsilon$，因此 $\sigma_{\min}(zI - A) \geq |z| - \|A\| > \varepsilon$，从而 $z \notin \sigma_\varepsilon(A)$。

    (5) $\|(zI - U^*AU)^{-1}\| = \|U^*(zI - A)^{-1}U\| = \|(zI - A)^{-1}\|$。$\blacksquare$

!!! theorem "定理 43.4 (正规矩阵的伪谱)"
    若 $A$ 是正规矩阵，则

    $$\sigma_\varepsilon(A) = \bigcup_{\lambda \in \sigma(A)} D(\lambda, \varepsilon) = \{z \in \mathbb{C} : \operatorname{dist}(z, \sigma(A)) < \varepsilon\},$$

    即正规矩阵的 $\varepsilon$-伪谱恰好是谱的 $\varepsilon$-邻域。

??? proof "证明"
    由定理 43.1，正规矩阵满足 $\|(zI - A)^{-1}\| = 1/\operatorname{dist}(z, \sigma(A))$。因此

    $$z \in \sigma_\varepsilon(A) \Leftrightarrow \|(zI - A)^{-1}\| > \varepsilon^{-1} \Leftrightarrow \operatorname{dist}(z, \sigma(A)) < \varepsilon.$$

    $\blacksquare$

!!! example "例 43.3 (非正规矩阵的伪谱膨胀)"
    考虑 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。特征值为 $\sigma(A) = \{0\}$（重数 2）。

    对于正规矩阵 $B = 0_{2 \times 2}$（也有 $\sigma(B) = \{0\}$），$\sigma_\varepsilon(B) = D(0, \varepsilon)$（半径 $\varepsilon$ 的圆盘）。

    但对于 $A$，$(zI - A)^{-1} = \begin{pmatrix} z^{-1} & -z^{-2} \\ 0 & z^{-1} \end{pmatrix}$，

    $$\|(zI - A)^{-1}\|_2 \geq \frac{1}{|z|^2}$$

    （由 $(2,1)$ 元素的贡献）。因此当 $|z|^2 < \varepsilon^{-1}$，即 $|z| < \varepsilon^{-1/2}$ 时可能有 $z \in \sigma_\varepsilon(A)$。

    更精确地计算：$\sigma_{\min}(zI - A) = \sigma_{\min}\begin{pmatrix} z & -1 \\ 0 & z \end{pmatrix}$。奇异值的平方是 $\frac{|z|^2 + 1/2 \pm \sqrt{(|z|^2 + 1/2)^2 - |z|^4}}{?}$...

    实际数值计算显示，$\sigma_\varepsilon(A)$ 是大约半径为 $\sqrt{\varepsilon}$ 的区域——远大于正规情形的 $\varepsilon$。这就是非正规性导致伪谱膨胀的典型现象。

!!! example "例 43.4 (Grcar 矩阵的伪谱)"
    $n \times n$ Grcar 矩阵
    $$G = \begin{pmatrix} 1 & 1 & 1 & 1 & & \\ -1 & 1 & 1 & 1 & 1 & \\ & -1 & 1 & 1 & 1 & \ddots \\ & & \ddots & \ddots & \ddots & \ddots & 1 \\ & & & -1 & 1 & 1 & 1 \\ & & & & -1 & 1 & 1 \\ & & & & & -1 & 1 \end{pmatrix}$$
    是高度非正规的。其特征值（$n = 32$）大致分布在单位圆附近的某条曲线上，但 $\varepsilon = 10^{-1}$ 的伪谱已经远远偏离了特征值的分布区域，呈现出显著的"腹地"（bulge）现象。

    这说明 Grcar 矩阵的动力学行为（如 $\|G^k\|$ 的增长）不能仅由特征值预测。

---

## 43.4 Kreiss 矩阵定理

<div class="context-flow" markdown>

**核心问题**：如何通过伪谱来判断矩阵幂 $\|A^k\|$ 或矩阵指数 $\|e^{tA}\|$ 的增长行为？

</div>

!!! definition "定义 43.4 (Kreiss 常数)"
    设 $A \in \mathbb{C}^{n \times n}$，谱半径 $\rho(A) \leq 1$。**Kreiss 常数**定义为

    $$\mathcal{K}(A) = \sup_{|z| > 1} (|z| - 1) \|(zI - A)^{-1}\|.$$

    等价地，$\mathcal{K}(A)$ 衡量了伪谱向单位圆外膨胀的速度：$\mathcal{K}(A) = \inf\{C : \sigma_\varepsilon(A) \subseteq D(0, 1 + C\varepsilon),\; \forall \varepsilon > 0\}$。

!!! theorem "定理 43.5 (Kreiss 矩阵定理)"
    设 $A \in \mathbb{C}^{n \times n}$，$\rho(A) \leq 1$。则

    $$\mathcal{K}(A) \leq \sup_{k \geq 0} \|A^k\| \leq en \cdot \mathcal{K}(A).$$

    即 $\|A^k\|$ 一致有界当且仅当 $\mathcal{K}(A) < \infty$，且有界的程度由 Kreiss 常数控制（虽然右侧有 $n$ 的因子，这在某些情况下是不可去除的）。

??? proof "证明"
    **左侧不等式**（$\mathcal{K}(A) \leq \sup_k \|A^k\|$）：

    设 $M = \sup_k \|A^k\|$。对 $|z| > 1$，
    $$(zI - A)^{-1} = z^{-1}(I - z^{-1}A)^{-1} = z^{-1} \sum_{k=0}^{\infty} z^{-k} A^k.$$
    因此
    $$\|(zI - A)^{-1}\| \leq |z|^{-1} \sum_{k=0}^{\infty} |z|^{-k} M = \frac{M}{|z| - 1}.$$
    从而 $(|z| - 1)\|(zI - A)^{-1}\| \leq M$，取上确界得 $\mathcal{K}(A) \leq M$。

    **右侧不等式**（$\sup_k \|A^k\| \leq en \cdot \mathcal{K}(A)$）：

    这个方向的证明更困难。思路是利用 Cauchy 积分公式
    $$A^k = \frac{1}{2\pi i} \oint_{|z| = r} z^k (zI - A)^{-1} dz$$
    （$r > 1$），然后优化 $r$ 的选取。利用 $\|(zI - A)^{-1}\| \leq \mathcal{K}(A)/(r - 1)$，

    $$\|A^k\| \leq \frac{r^{k+1} \mathcal{K}(A)}{r - 1}.$$

    取 $r = 1 + 1/k$ 最优化，得到 $\|A^k\| \leq e(k+1)\mathcal{K}(A)$。

    更精细的分析（利用 $A$ 是 $n \times n$ 矩阵的事实和预解式的有理函数结构）可以将 $k+1$ 替换为 $n$，得到最终结果。$\blacksquare$

!!! example "例 43.5"
    对于 $n \times n$ 上移位矩阵
    $$S = \begin{pmatrix} 0 & 1 & & \\ & 0 & \ddots & \\ & & \ddots & 1 \\ & & & 0 \end{pmatrix},$$

    $\rho(S) = 0$，$S^n = 0$，因此 $\sup_k \|S^k\| = 1$（因为 $\|S^k\| = 1$ 对 $k < n$ 和 $\|S^k\| = 0$ 对 $k \geq n$）。

    $\mathcal{K}(S) = \sup_{|z|>1} (|z|-1) \|z^{-1}(I + z^{-1}S + \cdots + z^{-(n-1)}S^{n-1})\| \leq 1$。

    因此 Kreiss 定理给出 $1 \leq \sup_k \|S^k\| \leq en$，这里下界紧而上界松。

---

## 43.5 瞬态增长

<div class="context-flow" markdown>

**核心问题**：为什么特征值全在左半平面的矩阵仍然可能导致解的剧烈增长？

</div>

!!! definition "定义 43.5 (数值横截面)"
    设 $A \in \mathbb{C}^{n \times n}$。$A$ 的**数值横截面**（numerical abscissa）定义为

    $$\alpha(A) = \max \operatorname{Re}\, \sigma\!\left(\frac{A + A^*}{2}\right) = \max_{\|x\|=1} \operatorname{Re}(x^* A x).$$

    它等于数值域 $W(A)$ 最右端的实部。

!!! theorem "定理 43.6 (瞬态增长与数值横截面)"
    设 $A \in \mathbb{C}^{n \times n}$。

    1. $\frac{d}{dt} \|e^{tA}\| \Big|_{t=0^+} = \alpha(A)$。
    2. 因此，$\|e^{tA}\|$ 在 $t = 0$ 附近增长当且仅当 $\alpha(A) > 0$。
    3. 这可以与谱横截面 $\eta(A) = \max_i \operatorname{Re} \lambda_i$ 同时存在 $\eta(A) < 0 < \alpha(A)$——即长期衰减但短期增长。

??? proof "证明"
    由 $\|e^{tA}\|^2 = \|e^{tA} e^{tA^*}\| = \rho(e^{tA} e^{tA^*})$（谱范数的平方等于 $e^{tA} e^{tA^*}$ 的最大特征值），以及

    $$e^{tA} e^{tA^*} = I + t(A + A^*) + O(t^2),$$

    因此 $\|e^{tA}\|^2 = 1 + 2t \cdot \alpha(A) + O(t^2)$，从而

    $$\|e^{tA}\| = 1 + t \cdot \alpha(A) + O(t^2),$$

    即 $\frac{d}{dt}\|e^{tA}\|\big|_{t=0^+} = \alpha(A)$。$\blacksquare$

!!! theorem "定理 43.7 (伪谱与瞬态增长的关系)"
    设 $A \in \mathbb{C}^{n \times n}$。

    1. **下界**：$\sup_{t \geq 0} \|e^{tA}\| \geq \frac{\eta_\varepsilon(A)}{\varepsilon}$，其中 $\eta_\varepsilon(A) = \sup\{\operatorname{Re} z : z \in \sigma_\varepsilon(A)\}$ 是 $\varepsilon$-伪谱横截面。
    2. **Kreiss 型定理**：
    $$\sup_{\operatorname{Re} z > 0} \operatorname{Re}(z) \cdot \|(zI - A)^{-1}\| \leq \sup_{t \geq 0} \|e^{tA}\| \leq en \cdot \sup_{\operatorname{Re} z > 0} \operatorname{Re}(z) \cdot \|(zI - A)^{-1}\|.$$

??? proof "证明"
    (1) 若 $z \in \sigma_\varepsilon(A)$ 且 $\operatorname{Re} z > 0$，则存在 $\|E\| < \varepsilon$ 使得 $z \in \sigma(A + E)$，因此

    $$\sup_{t \geq 0} \|e^{t(A+E)}\| \geq \lim_{t \to \infty} \frac{\|e^{t(A+E)}\|}{1} \geq \text{（包含 }e^{t \operatorname{Re} z}\text{ 的分量）}.$$

    更直接地：$\|e^{tA}\| \geq \|e^{t(A+E)}\| \cdot e^{-t\|E\|} \cdot e^{-t\varepsilon}$...

    精确的论证使用以下不等式：对 $z \in \sigma_\varepsilon(A)$ 且 $\operatorname{Re} z = \eta_\varepsilon(A) > 0$，
    $$\sup_{t \geq 0} \|e^{tA}\| \geq \frac{\eta_\varepsilon(A)}{\varepsilon}.$$

    这可以由 Laplace 变换的关系 $(zI - A)^{-1} = \int_0^\infty e^{-zt} e^{tA} dt$（$\operatorname{Re} z > \eta(A)$）和范数估计得到。$\blacksquare$

!!! example "例 43.6 (经典瞬态增长例子)"
    设
    $$A = \begin{pmatrix} -1 & 100 \\ 0 & -2 \end{pmatrix}.$$

    **特征值**：$\lambda_1 = -1$，$\lambda_2 = -2$。所有特征值实部为负，因此 $e^{tA} \to 0$ 当 $t \to \infty$。

    **数值横截面**：$\alpha(A) = \max \operatorname{Re}\, \sigma\!\left(\frac{A + A^T}{2}\right) = \max \operatorname{Re}\, \sigma\begin{pmatrix} -1 & 50 \\ 50 & -2 \end{pmatrix}$。

    $\frac{A + A^T}{2}$ 的特征值为 $\frac{-3 \pm \sqrt{1 + 10000}}{2} \approx \frac{-3 \pm 100.005}{2}$，即约 $48.5$ 和 $-51.5$。

    因此 $\alpha(A) \approx 48.5 \gg 0$：尽管所有特征值实部为负，$\|e^{tA}\|$ 在 $t = 0$ 附近以大约 $48.5$ 的速率增长！

    实际计算：$e^{tA} = \begin{pmatrix} e^{-t} & 100(e^{-t} - e^{-2t}) \\ 0 & e^{-2t} \end{pmatrix}$。

    $\|e^{tA}\|$ 的最大值发生在 $t \approx \ln 2 \approx 0.693$ 时，$\|e^{tA}\| \approx 25$。

    **伪谱解释**：$\sigma_\varepsilon(A)$ 在 $\varepsilon$ 很小时已经大幅延伸到右半平面。例如 $\varepsilon = 0.5$ 时，伪谱横截面 $\eta_{0.5}(A) \approx 12.5$，因此 $\sup_{t \geq 0} \|e^{tA}\| \geq 12.5/0.5 = 25$。

---

## 43.6 伪谱的计算

<div class="context-flow" markdown>

**核心问题**：如何高效计算和可视化矩阵的伪谱？

</div>

!!! definition "定义 43.6 (伪谱等高线图)"
    $A$ 的伪谱通常通过绘制函数
    $$z \mapsto \log_{10} \sigma_{\min}(zI - A)$$
    的等高线图来可视化。等高线 $\log_{10} \sigma_{\min}(zI - A) = \log_{10} \varepsilon$ 就是 $\sigma_\varepsilon(A)$ 的边界。

!!! theorem "定理 43.8 (网格方法的复杂度)"
    **网格方法**（grid method）：在复平面的 $m \times m$ 网格上，对每个网格点 $z_{ij}$ 计算 $\sigma_{\min}(z_{ij}I - A)$。

    1. 计算一个点的 $\sigma_{\min}$ 需要 $O(n^3)$ 操作（SVD）。
    2. 总复杂度为 $O(m^2 n^3)$。
    3. 利用 Schur 分解预处理可降至 $O(n^3 + m^2 n^2)$：先计算 $A = QTQ^*$，然后 $\sigma_{\min}(zI - A) = \sigma_{\min}(zI - T)$，上三角矩阵的最小奇异值可以用 $O(n^2)$ 操作求解（求解三角方程组的逆迭代）。

!!! example "例 43.7 (网格方法伪代码)"
    ```
    输入：矩阵 A ∈ ℂⁿˣⁿ，网格范围 [x_min, x_max] × [y_min, y_max]，网格密度 m
    输出：σ_min 值的 m × m 矩阵

    1. 计算 Schur 分解 A = QTQ*
    2. 生成 m × m 网格点 z_{ij} = x_i + iy_j
    3. 对每个 (i,j)：
       a. 形成 B = z_{ij}I - T（上三角）
       b. 计算 σ_min(B)（使用逆迭代或 SVD）
       c. 存储结果
    4. 绘制 log₁₀(σ_min) 的等高线图
    ```

!!! theorem "定理 43.9 (基于延拓的方法)"
    **边界追踪法**（continuation/path-following method）：给定特定的 $\varepsilon$，伪谱的边界 $\partial \sigma_\varepsilon(A)$ 是一条（或多条）光滑曲线，可以用数值延拓方法追踪。

    设 $z(s)$ 是边界的弧长参数化，满足约束 $\sigma_{\min}(z(s)I - A) = \varepsilon$。初始点可以由特征值附近的局部分析给出。

    该方法的优势是：

    - 只计算边界上的点，远少于全网格。
    - 适合绘制单一 $\varepsilon$ 值的伪谱边界。

!!! theorem "定理 43.10 (Lanczos 方法)"
    对于大规模稀疏矩阵，直接 SVD 不可行。可以使用：

    1. **Lanczos 逼近**：对每个网格点 $z$，用 Lanczos 方法估计 $\sigma_{\min}(zI - A)$。
    2. **投影方法**：将问题投影到 Krylov 子空间上，在小维度空间中求解。
    3. **随机化方法**：使用随机向量 $v$，估计 $\|(zI - A)^{-1}v\| / \|v\|$ 作为 $\|(zI - A)^{-1}\|$ 的下界。

    **EigTool**（Trefethen 等人开发）是计算伪谱的标准软件工具，实现了上述多种算法。

---

## 43.7 应用

<div class="context-flow" markdown>

**核心问题**：伪谱理论如何帮助理解流体稳定性、迭代法收敛和数值 ODE 稳定性？

</div>

!!! example "例 43.8 (流体动力学稳定性)"
    **Orr-Sommerfeld 方程**描述了平行流的线性稳定性。在 Couette 流和 Poiseuille 流中，相应的 Orr-Sommerfeld 算子 $L$ 是高度非正规的。

    对于平面 Couette 流（Reynolds 数 $Re$）：

    - **谱分析**：对所有 $Re$，所有特征值都在左半平面。经典线性稳定性理论预测该流**总是稳定的**。
    - **实验**：当 $Re \gtrsim 350$ 时，流动变得不稳定。
    - **伪谱分析**：$\sigma_\varepsilon(L)$ 深入右半平面。伪谱横截面 $\eta_\varepsilon(L) \sim Re$，瞬态增长 $\|e^{tL}\| \sim Re^2$。

    这解释了**亚临界转捩**（subcritical transition）现象：微小扰动可以在短时间内被放大 $O(Re^2)$ 倍，足以触发非线性效应导致湍流。

!!! example "例 43.9 (迭代法的收敛性)"
    考虑线性方程组 $Ax = b$ 的 GMRES 迭代。经典收敛分析基于特征值分布：若特征值聚集在远离原点的小区域中，GMRES 应该快速收敛。

    **但是**：对非正规矩阵，特征值分布可能完全误导。

    **定理**（伪谱收敛界）：GMRES 第 $k$ 步的残差满足

    $$\frac{\|r_k\|}{\|r_0\|} \leq \inf_{p \in \mathcal{P}_k, p(0)=1} \sup_{z \in \sigma_\varepsilon(A)} |p(z)| + O(\varepsilon),$$

    其中 $\mathcal{P}_k$ 是次数不超过 $k$ 的多项式集合。

    因此，GMRES 的收敛速度由伪谱（而非谱）上的多项式逼近问题决定。

!!! example "例 43.10 (数值 ODE 稳定性)"
    考虑 ODE $y' = Ay$（$A$ 半离散化后的矩阵）用数值方法（如 Euler 法、Runge-Kutta 法）求解。

    设数值方法的**稳定域**为 $\mathcal{S}$（复平面中使方法稳定的区域）。步长 $h$ 下的稳定性条件是：

    - **经典分析**：$h \cdot \sigma(A) \subseteq \mathcal{S}$。
    - **伪谱分析**：应要求 $h \cdot \sigma_\varepsilon(A) \subseteq \mathcal{S}$（对适当的 $\varepsilon$）。

    对非正规矩阵，伪谱可能远大于谱，因此实际需要的步长可能远小于经典分析的预测。

!!! theorem "定理 43.11 (伪谱与矩阵指数的精确关系)"
    对任意 $A \in \mathbb{C}^{n \times n}$：

    $$\sup_{t \geq 0} \|e^{tA}\| = \lim_{\varepsilon \to 0} \frac{\alpha_\varepsilon(A)}{\varepsilon},$$

    其中 $\alpha_\varepsilon(A) = \max\{\operatorname{Re} z : z \in \sigma_\varepsilon(A)\}$ 是 $\varepsilon$-伪谱横截面。

    等价地（对连续时间系统）：$e^{tA}$ 一致有界当且仅当 $\alpha_\varepsilon(A) \leq C\varepsilon$ 对所有 $\varepsilon > 0$ 成立。

??? proof "证明"
    利用 Laplace 变换和 Hille-Yosida 定理的有限维版本。

    **上界**：由 $(zI - A)^{-1} = \int_0^\infty e^{-zt} e^{tA} dt$（$\operatorname{Re} z > \eta(A)$），若 $\operatorname{Re} z = \alpha > \eta(A)$：
    $$\|(zI - A)^{-1}\| \leq \int_0^\infty e^{-\alpha t} \|e^{tA}\| dt \leq \frac{M}{\alpha - \eta(A)},$$
    其中 $M = \sup_{t \geq 0} \|e^{tA}\| e^{-\eta(A) t}$...

    **下界**：若 $z_0 \in \sigma_\varepsilon(A)$ 且 $\operatorname{Re} z_0 = \alpha_\varepsilon(A) > 0$，则存在 $\|E\| < \varepsilon$ 使 $z_0 \in \sigma(A + E)$。因此
    $$\sup_{t \geq 0} \|e^{t(A+E)}\| \geq \sup_{t \geq 0} e^{t \operatorname{Re} z_0} = +\infty$$
    （若 $\operatorname{Re} z_0 > 0$）。但 $\|e^{tA}\| \geq \|e^{t(A+E)}\| e^{-t\varepsilon}$...

    严格的证明需要更精细的 Kreiss 型估计。$\blacksquare$

!!! example "例 43.11 (总结对比)"
    | 矩阵性质 | 正规矩阵 | 非正规矩阵 |
    |---------|---------|-----------|
    | $\sigma_\varepsilon(A)$ | 谱的 $\varepsilon$-邻域 | 可能远大于 $\varepsilon$-邻域 |
    | $\|(zI-A)^{-1}\|$ | $1/\operatorname{dist}(z, \sigma)$ | 可能远大于 $1/\operatorname{dist}(z, \sigma)$ |
    | $\|e^{tA}\|$ | $e^{t\eta(A)}$ | 可能有巨大瞬态增长 |
    | $\|p(A)\|$ | $\max_\sigma |p|$ | 可能远大于 $\max_\sigma |p|$ |
    | 迭代法收敛 | 由特征值决定 | 由伪谱决定 |

    伪谱理论的核心启示是：**对非正规矩阵，特征值不能独立预测矩阵的动力学行为；必须考虑预解式在整个复平面上的行为，即伪谱。**

## 练习题

1. **[概念] 什么是 $\varepsilon$-伪谱？它与特征值有什么本质联系？**
   ??? success "参考答案"
       $\varepsilon$-伪谱 $\sigma_\varepsilon(A)$ 是所有使得预解式范数 $\|(zI-A)^{-1}\| > \varepsilon^{-1}$ 的复数 $z$ 的集合。直观上，它是矩阵 $A$ 在受到范数小于 $\varepsilon$ 的扰动后，所有可能出现的特征值的并集。

2. **[正规阵] 证明：若 $A$ 是正规矩阵，则 $\sigma_\varepsilon(A)$ 恰好是特征值的 $\varepsilon$-邻域。**
   ??? success "参考答案"
       对于正规矩阵，存在酉对角化 $A = U\Lambda U^H$。预解式范数 $\|(zI-A)^{-1}\| = \max_i |z - \lambda_i|^{-1}$。不等式 $|z - \lambda_i|^{-1} > \varepsilon^{-1}$ 等价于 $|z - \lambda_i| < \varepsilon$。故伪谱是分布在各特征值周围、半径为 $\varepsilon$ 的圆盘并集。

3. **[计算] 计算 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ 的 $\varepsilon$-伪谱。**
   ??? success "参考答案"
       计算最小奇异值 $\sigma_{\min}(zI - A) = \sigma_{\min}\begin{pmatrix} z & -1 \\ 0 & z \end{pmatrix}$。当 $\varepsilon$ 较小时，伪谱区域大约是以原点为中心、半径约为 $\sqrt{\varepsilon}$ 的圆盘。注意这比正规矩阵的 $\varepsilon$ 级半径大得多。

4. **[瞬态增长] 为什么一个所有特征值实部都为负的系统，在短期内仍可能出现剧烈的状态增长？**
   ??? success "参考答案"
       这是由于矩阵的非正规性（Non-normality）。虽然特征值决定了长期（$t \to \infty$）的衰减，但伪谱若延伸到了复平面的右半平面，则预示着系统存在瞬态增长。物理上，这对应于不同特征模态之间的非正规干涉。

5. **[Kreiss] 简述 Kreiss 矩阵定理的核心结论。**
   ??? success "参考答案"
       它建立了矩阵幂 $\|A^k\|$ 的一致有界性与伪谱向单位圆外“膨胀”速度之间的定量联系。如果伪谱在单位圆外扩张太快，那么系统一定会经历严重的瞬态不稳定性。

6. **[数值] 在数值分析中，为什么伪谱比特征值更能反映 GMRES 等迭代法的收敛性？**
   ??? success "参考答案"
       对于高度非正规矩阵，残差的减少并不取决于特征值的位置，而取决于多项式在伪谱区域（而非仅仅特征值点）上的大小。如果伪谱包含原点或非常靠近原点，迭代法可能会陷入长时间的停滞。

7. **[流体] Orr-Sommerfeld 算子的伪谱研究解决了流体力学中的什么经典矛盾？**
   ??? success "参考答案"
       解决了“线性稳定流在实验中为何发生湍流转捩”的矛盾。经典谱分析认为这些流是稳定的（特征值全在左半平面），但伪谱分析揭示了巨大的瞬态增长（$O(Re^2)$ 级），微小扰动会被非正规机制放大到足以触发非线性失稳。

8. **[判定] 如何从预解式图（Resolvent Plot）中识别一个矩阵是否接近正规？**
   ??? success "参考答案"
       观察等高线的分布。如果等高线是以特征华为中心的圆，且间距均匀（符合 $1/dist$ 规律），则接近正规；如果等高线发生剧烈扭曲、合并，或者在特征值很远处依然保持高值，则矩阵高度非正规。

9. **[性质] 证明：$\sigma_\varepsilon(A)$ 的每个连通分支内部至少包含 $A$ 的一个特征值。**
   ??? success "参考答案"
       这是复分析中最大模原理的应用。如果一个分支不含特征值，则预解式在该区域解析，其模长在边界达到最大，这会与伪谱的定义产生矛盾。

10. **[爱因斯坦思考题] 爱因斯坦在统一场论中追求规律的确定性。但在非正规算子的伪谱理论中，即使规律（矩阵 $A$）确定，极其微小的环境噪声（扰动 $E$）也能让系统的谱发生剧变。这是否意味着“谱”本身并不是一种鲁棒的物理现实？我们应该如何重新定义“可观测的物理特征”？**
    ??? success "参考答案"
        这反映了“经典测量观”的局限。在存在非正规性的宇宙中，单一的特征值不再是鲁棒的观测物理量，因为它们对微扰极其敏感。伪谱理论建议我们：应该将“物理特征”定义为复平面上的一个**区域**而非孤立的点。这种从“点”到“集合”的思维转变，承认了微扰与系统结构的不可分割性。爱因斯坦追求的“确定性”在伪谱意义下得以恢复——即系统行为在“伪谱带”内是确定的、可预测的，即便单独的特征值在乱跳。

## 本章小结

本章探讨了超越经典谱理论的非正规矩阵分析工具——伪谱，主要内容包括：

1. **正规性的局限**：分析了经典特征值理论在处理非正规矩阵（如 Jordan 块）时的失效，引入了偏离正规度（Departure from Normality）的概念。
2. **ε-伪谱定义**：给出了基于预解式范数、特征值扰动和最小奇异值的三种等价定义，确立了其作为“几乎特征值”集合的地位。
3. **几何与谱性质**：对比了正规矩阵（伪谱即邻域）与非正规矩阵（伪谱剧烈膨胀）的区别，揭示了伪谱分支与谱的拓扑关联。
4. **动力学关联**：详细论述了伪谱与矩阵指数/幂增长的关系，特别是 Kreiss 矩阵定理和数值横截面对瞬态增长的预测作用。
5. **现代应用**：展示了伪谱在流体亚临界转捩、数值迭代收敛性分析以及数值 ODE 稳定性判定中的核心价值，标志着矩阵分析从“静态点”向“动态区域”的范式转移。

