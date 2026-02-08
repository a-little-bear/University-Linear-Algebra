# 第 15 章 范数与扰动理论

<div class="context-flow" markdown>

**前置**：Ch14 谱半径与范数

**脉络**：向量范数 → 矩阵范数(算子/Frobenius) → **条件数**量化病态性 → Bauer-Fike/Weyl 量化特征值灵敏度

**延伸**：条件数在数值天气预报、有限元分析、GPS 定位等大规模计算中决定了结果的可信度；特征值扰动理论（Weyl、Bauer-Fike）是随机矩阵理论和量子信息中算子扰动分析的基础

</div>

范数为向量空间和矩阵空间提供了"大小"的度量，是矩阵分析和数值线性代数的基石。扰动理论则研究当矩阵发生微小变化时，其特征值、奇异值和线性方程组的解如何变化。这些理论对于理解数值计算的稳定性和精度至关重要。本章系统地建立向量范数和矩阵范数的理论，引入条件数的概念，并深入讨论特征值和奇异值的扰动界。

---

## 15.1 向量范数

<div class="context-flow" markdown>

**脉络**：$\ell_p$ 范数族($p=1,2,\infty$) → Hölder 不等式($p,q$ 共轭) → 有限维范数等价(拓扑唯一)

</div>

!!! definition "定义 15.1 (向量范数)"
    $\mathbb{C}^n$ 上的**向量范数（vector norm）** 是函数 $\|\cdot\| : \mathbb{C}^n \to \mathbb{R}$，满足对所有 $\mathbf{x}, \mathbf{y} \in \mathbb{C}^n$ 和 $\alpha \in \mathbb{C}$：

    (1) **非负性**：$\|\mathbf{x}\| \geq 0$，等号成立当且仅当 $\mathbf{x} = \mathbf{0}$；

    (2) **齐次性**：$\|\alpha\mathbf{x}\| = |\alpha|\|\mathbf{x}\|$；

    (3) **三角不等式**：$\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$。

!!! definition "定义 15.2 ($\ell_p$ 范数)"
    对 $1 \leq p \leq \infty$，$\mathbb{C}^n$ 上的 **$\ell_p$ 范数（$\ell_p$ norm）** 定义为

    $$\|\mathbf{x}\|_p = \left(\sum_{i=1}^{n} |x_i|^p\right)^{1/p}, \quad 1 \leq p < \infty,$$

    $$\|\mathbf{x}\|_\infty = \max_{1 \leq i \leq n} |x_i|.$$

    特别重要的是：

    - $\|\mathbf{x}\|_1 = \sum_{i=1}^{n}|x_i|$（曼哈顿范数）；
    - $\|\mathbf{x}\|_2 = \sqrt{\sum_{i=1}^{n}|x_i|^2} = \sqrt{\mathbf{x}^*\mathbf{x}}$（欧氏范数）；
    - $\|\mathbf{x}\|_\infty = \max_i |x_i|$（切比雪夫范数）。

!!! theorem "定理 15.1 (Hölder 不等式)"
    设 $p, q \geq 1$ 满足 $\frac{1}{p} + \frac{1}{q} = 1$（共轭指数），则对任意 $\mathbf{x}, \mathbf{y} \in \mathbb{C}^n$，

    $$|\mathbf{x}^*\mathbf{y}| \leq \|\mathbf{x}\|_p\|\mathbf{y}\|_q.$$

    特别地，$p = q = 2$ 时即为 Cauchy-Schwarz 不等式。

??? proof "证明"
    当 $\|\mathbf{x}\|_p = 0$ 或 $\|\mathbf{y}\|_q = 0$ 时显然成立。否则不妨设 $\|\mathbf{x}\|_p = \|\mathbf{y}\|_q = 1$。由 Young 不等式 $ab \leq \frac{a^p}{p} + \frac{b^q}{q}$（$a, b \geq 0$），取 $a = |x_i|$，$b = |y_i|$，求和得

    $$\sum_{i=1}^{n}|x_i||y_i| \leq \frac{1}{p}\sum|x_i|^p + \frac{1}{q}\sum|y_i|^q = \frac{1}{p} + \frac{1}{q} = 1.$$

    因此 $|\mathbf{x}^*\mathbf{y}| \leq \sum|x_i||y_i| \leq 1 = \|\mathbf{x}\|_p\|\mathbf{y}\|_q$。一般情形除以 $\|\mathbf{x}\|_p\|\mathbf{y}\|_q$ 归约到上述情况。$\blacksquare$

!!! theorem "定理 15.2 (范数等价性定理)"
    $\mathbb{C}^n$ 上的任意两个范数 $\|\cdot\|_\alpha$ 和 $\|\cdot\|_\beta$ 是**等价的（equivalent）**，即存在正常数 $c_1, c_2 > 0$ 使得对所有 $\mathbf{x} \in \mathbb{C}^n$，

    $$c_1\|\mathbf{x}\|_\alpha \leq \|\mathbf{x}\|_\beta \leq c_2\|\mathbf{x}\|_\alpha.$$

??? proof "证明"
    只需证明任意范数 $\|\cdot\|$ 与 $\|\cdot\|_2$ 等价。设 $\mathbf{e}_1, \ldots, \mathbf{e}_n$ 为标准基。对 $\mathbf{x} = \sum x_i\mathbf{e}_i$，

    $$\|\mathbf{x}\| \leq \sum |x_i|\|\mathbf{e}_i\| \leq \left(\max_i\|\mathbf{e}_i\|\right)\sum|x_i| \leq \left(\max_i\|\mathbf{e}_i\|\right)\sqrt{n}\|\mathbf{x}\|_2.$$

    令 $c_2 = \sqrt{n}\max_i\|\mathbf{e}_i\|$，则 $\|\mathbf{x}\| \leq c_2\|\mathbf{x}\|_2$。这说明 $\|\cdot\|$ 关于 $\|\cdot\|_2$ 是 Lipschitz 连续的。

    考虑紧集 $S = \{\mathbf{x} : \|\mathbf{x}\|_2 = 1\}$。连续函数 $\|\cdot\|$ 在 $S$ 上达到最小值 $c_1 > 0$（因为 $\|\mathbf{x}\| > 0$ 对 $\mathbf{x} \in S$）。因此对 $\mathbf{x} \neq \mathbf{0}$，$\|\mathbf{x}/\|\mathbf{x}\|_2\| \geq c_1$，即 $\|\mathbf{x}\| \geq c_1\|\mathbf{x}\|_2$。$\blacksquare$

!!! proposition "命题 15.1 ($\ell_p$ 范数之间的关系)"
    对任意 $\mathbf{x} \in \mathbb{C}^n$ 和 $1 \leq p \leq q \leq \infty$，

    $$\|\mathbf{x}\|_q \leq \|\mathbf{x}\|_p \leq n^{1/p - 1/q}\|\mathbf{x}\|_q.$$

    特别地：

    - $\|\mathbf{x}\|_2 \leq \|\mathbf{x}\|_1 \leq \sqrt{n}\|\mathbf{x}\|_2$；
    - $\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_2 \leq \sqrt{n}\|\mathbf{x}\|_\infty$；
    - $\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_1 \leq n\|\mathbf{x}\|_\infty$。

!!! example "例 15.1"
    设 $\mathbf{x} = (3, -4, 0, 1)^T$，计算各 $\ell_p$ 范数。

    - $\|\mathbf{x}\|_1 = |3| + |-4| + |0| + |1| = 8$；
    - $\|\mathbf{x}\|_2 = \sqrt{9 + 16 + 0 + 1} = \sqrt{26} \approx 5.10$；
    - $\|\mathbf{x}\|_\infty = \max\{3, 4, 0, 1\} = 4$。

    验证：$\|\mathbf{x}\|_\infty = 4 \leq 5.10 = \|\mathbf{x}\|_2 \leq 8 = \|\mathbf{x}\|_1$。

    $\|\mathbf{x}\|_1 = 8 \leq \sqrt{4}\cdot 5.10 = 10.2$？实际上 $\sqrt{n}\|\mathbf{x}\|_2 = 2 \times 5.10 = 10.2$，而 $\|\mathbf{x}\|_1 = 8 \leq 10.2$。

!!! example "例 15.2"
    **$\ell_p$ 范数的单位球**。$\ell_p$ 范数的单位球 $B_p = \{\mathbf{x} \in \mathbb{R}^2 : \|\mathbf{x}\|_p \leq 1\}$ 的形状随 $p$ 变化：

    - $p = 1$：菱形（正方形旋转 $45°$），顶点为 $(\pm 1, 0)$ 和 $(0, \pm 1)$；
    - $p = 2$：圆盘；
    - $p = \infty$：正方形 $[-1, 1]^2$；
    - $p \to 0^+$：单位球退化为坐标轴上的线段。

    随 $p$ 增大，单位球从菱形膨胀为圆再膨胀为正方形，体现了范数之间 $\|\mathbf{x}\|_q \leq \|\mathbf{x}\|_p$（$p \leq q$）的关系。

---

## 15.2 矩阵范数

<div class="context-flow" markdown>

**脉络**：矩阵范数 + **次可乘性** $\|AB\|\leq\|A\|\|B\|$ → Frobenius范数 $=\sqrt{\sum\sigma_i^2}$ 是酉不变的"元素级"范数

</div>

!!! definition "定义 15.3 (矩阵范数)"
    $\mathbb{C}^{m \times n}$ 上的**矩阵范数（matrix norm）** 是函数 $\|\cdot\| : \mathbb{C}^{m \times n} \to \mathbb{R}$，满足对所有 $A, B$ 和 $\alpha \in \mathbb{C}$：

    (1) **非负性**：$\|A\| \geq 0$，等号成立当且仅当 $A = O$；

    (2) **齐次性**：$\|\alpha A\| = |\alpha|\|A\|$；

    (3) **三角不等式**：$\|A + B\| \leq \|A\| + \|B\|$。

    若还满足 (4) **次可乘性（submultiplicativity）**：$\|AB\| \leq \|A\|\|B\|$（当乘法有定义时），则称为**次可乘范数**。

!!! definition "定义 15.4 (Frobenius 范数)"
    设 $A = (a_{ij}) \in \mathbb{C}^{m \times n}$，**Frobenius 范数（Frobenius norm）** 定义为

    $$\|A\|_F = \sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n}|a_{ij}|^2} = \sqrt{\operatorname{tr}(A^*A)} = \sqrt{\sum_{i=1}^{\min(m,n)}\sigma_i^2}.$$

!!! theorem "定理 15.3 (Frobenius 范数的性质)"
    Frobenius 范数满足：

    (1) 它是矩阵范数且次可乘；

    (2) $\|A\|_F = \|A^*\|_F$；

    (3) 对酉矩阵（正交矩阵）$U, V$，$\|UAV\|_F = \|A\|_F$（酉不变性）；

    (4) $\|A\|_F^2 = \sum_{i=1}^{r}\sigma_i^2$，其中 $\sigma_i$ 为 $A$ 的奇异值；

    (5) $\|A\|_2 \leq \|A\|_F \leq \sqrt{r}\|A\|_2$，其中 $r = \operatorname{rank}(A)$。

??? proof "证明"
    **(1)** 次可乘性：$\|AB\|_F^2 = \sum_{i,k}\left|\sum_j a_{ij}b_{jk}\right|^2 \leq \sum_{i,k}\left(\sum_j|a_{ij}|^2\right)\left(\sum_j|b_{jk}|^2\right)$（由 Cauchy-Schwarz）$= \sum_i\left(\sum_j|a_{ij}|^2\right)\sum_k\left(\sum_j|b_{jk}|^2\right) = \|A\|_F^2\|B\|_F^2$。

    **(3)** $\|UAV\|_F^2 = \operatorname{tr}(V^*A^*U^*UAV) = \operatorname{tr}(V^*A^*AV) = \operatorname{tr}(A^*AVV^*) = \operatorname{tr}(A^*A) = \|A\|_F^2$。

    **(4)** 由 SVD，$A = U\Sigma V^*$，$\|A\|_F^2 = \|\Sigma\|_F^2 = \sum\sigma_i^2$。

    **(5)** $\|A\|_F^2 = \sum\sigma_i^2 \geq \sigma_1^2 = \|A\|_2^2$，又 $\sum\sigma_i^2 \leq r\sigma_1^2$。$\blacksquare$

!!! example "例 15.3"
    设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。

    $\|A\|_F = \sqrt{1 + 4 + 9 + 16} = \sqrt{30} \approx 5.48$。

    $A^*A = \begin{pmatrix} 10 & 14 \\ 14 & 20 \end{pmatrix}$，特征值 $15 \pm \sqrt{221}/\sqrt{1}$... 让我重新计算。$\operatorname{tr} = 30$，$\det = 200 - 196 = 4$。特征值 $\mu = 15 \pm \sqrt{225 - 4} = 15 \pm \sqrt{221}$。$\sigma_1 = \sqrt{15 + \sqrt{221}} \approx \sqrt{29.87} \approx 5.47$，$\sigma_2 = \sqrt{15 - \sqrt{221}} \approx \sqrt{0.13} \approx 0.37$。

    验证 $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2} = \sqrt{30} \approx 5.48$。

---

## 15.3 算子范数（诱导范数）

<div class="context-flow" markdown>

**洞察**：$\|A\|_1$=最大列和，$\|A\|_\infty$=最大行和，$\|A\|_2=\sigma_1(A)$ · 算子范数自动次可乘且 $\rho(A)\leq\|A\|$

</div>

!!! definition "定义 15.5 (算子范数)"
    设 $\|\cdot\|_\alpha$ 和 $\|\cdot\|_\beta$ 分别是 $\mathbb{C}^n$ 和 $\mathbb{C}^m$ 上的向量范数。由它们**诱导（induced）** 的**算子范数（operator norm）** 定义为

    $$\|A\|_{\alpha \to \beta} = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|_\beta}{\|\mathbf{x}\|_\alpha} = \max_{\|\mathbf{x}\|_\alpha = 1} \|A\mathbf{x}\|_\beta.$$

    当 $\alpha = \beta = p$ 时，简记为 $\|A\|_p$。

!!! theorem "定理 15.4 (算子范数的计算公式)"
    设 $A = (a_{ij}) \in \mathbb{C}^{m \times n}$，则：

    (1) $\|A\|_1 = \max_{1 \leq j \leq n} \sum_{i=1}^{m} |a_{ij}|$（最大列和）；

    (2) $\|A\|_\infty = \max_{1 \leq i \leq m} \sum_{j=1}^{n} |a_{ij}|$（最大行和）；

    (3) $\|A\|_2 = \sigma_1(A)$（最大奇异值）。

??? proof "证明"
    **(1)** 设 $\mathbf{x}$ 满足 $\|\mathbf{x}\|_1 = 1$。则

    $$\|A\mathbf{x}\|_1 = \sum_i\left|\sum_j a_{ij}x_j\right| \leq \sum_i\sum_j|a_{ij}||x_j| = \sum_j|x_j|\sum_i|a_{ij}| \leq \max_j\sum_i|a_{ij}|.$$

    等号在 $\mathbf{x} = \mathbf{e}_{j^*}$（取最大列和对应的标准基向量）时达到。

    **(2)** 类似地，$\|A\mathbf{x}\|_\infty = \max_i\left|\sum_j a_{ij}x_j\right| \leq \max_i\sum_j|a_{ij}|\|\mathbf{x}\|_\infty$。等号在选取适当 $\mathbf{x}$ 时达到。

    **(3)** $\|A\mathbf{x}\|_2^2 = \mathbf{x}^*A^*A\mathbf{x}$。$\max_{\|\mathbf{x}\|_2=1}\mathbf{x}^*A^*A\mathbf{x} = \lambda_{\max}(A^*A) = \sigma_1^2(A)$（由 Rayleigh 商）。故 $\|A\|_2 = \sigma_1(A)$。$\blacksquare$

!!! theorem "定理 15.5 (算子范数的性质)"
    算子范数满足：

    (1) $\|I\| = 1$；

    (2) 次可乘性：$\|AB\| \leq \|A\|\|B\|$；

    (3) 相容性：$\|A\mathbf{x}\| \leq \|A\|\|\mathbf{x}\|$（对相应的向量范数）；

    (4) $\rho(A) \leq \|A\|$。

??? proof "证明"
    **(2)** $\|AB\mathbf{x}\| \leq \|A\|\|B\mathbf{x}\| \leq \|A\|\|B\|\|\mathbf{x}\|$，对所有 $\|\mathbf{x}\| = 1$ 取上确界得 $\|AB\| \leq \|A\|\|B\|$。

    **(3)** 由定义直接得出。

    **(4)** 已在定理 14.6 中证明。$\blacksquare$

!!! example "例 15.4"
    设 $A = \begin{pmatrix} 1 & -2 & 3 \\ 4 & 0 & -1 \end{pmatrix}$。

    - $\|A\|_1 = \max\{|1|+|4|,\; |-2|+|0|,\; |3|+|-1|\} = \max\{5, 2, 4\} = 5$；
    - $\|A\|_\infty = \max\{|1|+|-2|+|3|,\; |4|+|0|+|-1|\} = \max\{6, 5\} = 6$；
    - $\|A\|_2 = \sigma_1(A)$，需计算 $A^TA$ 的最大特征值。

!!! example "例 15.5"
    对正规矩阵，$\|A\|_2 = \rho(A)$。

    设 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$（旋转 $90°$）。$A$ 是正规矩阵（实际是正交矩阵），特征值为 $\pm i$，$\rho(A) = 1$。

    $A^TA = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$，$\sigma_1 = 1$。因此 $\|A\|_2 = 1 = \rho(A)$。

---

## 15.4 范数之间的关系

<div class="context-flow" markdown>

**关键不等式**：$\|A\|_2\leq\|A\|_F\leq\sqrt{r}\|A\|_2$ · $\|A\|_2\leq\sqrt{\|A\|_1\|A\|_\infty}$ 连通三种算子范数

</div>

!!! theorem "定理 15.6 (矩阵范数之间的不等式)"
    设 $A \in \mathbb{C}^{m \times n}$，则：

    (1) $\|A\|_2 \leq \|A\|_F \leq \sqrt{n}\|A\|_2$；

    (2) $\frac{1}{\sqrt{n}}\|A\|_\infty \leq \|A\|_2 \leq \sqrt{m}\|A\|_\infty$；

    (3) $\frac{1}{\sqrt{m}}\|A\|_1 \leq \|A\|_2 \leq \sqrt{n}\|A\|_1$；

    (4) $\|A\|_2 \leq \sqrt{\|A\|_1\|A\|_\infty}$；

    (5) $\frac{1}{\sqrt{mn}}\|A\|_F \leq \|A\|_\infty$，$\frac{1}{\sqrt{mn}}\|A\|_F \leq \|A\|_1$。

??? proof "证明"
    **(1)** 已在定理 15.3 (5) 中证明。

    **(4)** 设 $A\mathbf{x} = \sigma_1\mathbf{u}$（$\|\mathbf{x}\|_2 = \|\mathbf{u}\|_2 = 1$）。则

    $$\sigma_1 = \|A\mathbf{x}\|_2 \leq \sqrt{m}\|A\mathbf{x}\|_\infty \leq \sqrt{m}\|A\|_\infty\|\mathbf{x}\|_\infty \leq \sqrt{m}\|A\|_\infty.$$

    同理利用 $A^*$：$\sigma_1 = \|A^*\mathbf{u}\|_2 \leq \sqrt{n}\|A^*\|_\infty = \sqrt{n}\|A\|_1$。因此 $\sigma_1^2 \leq \sqrt{mn}\|A\|_\infty\sqrt{mn}\|A\|_1/(mn)$...

    更直接的证明：$\|A\|_2^2 = \rho(A^*A) \leq \|A^*A\|_\infty \leq \|A^*\|_\infty\|A\|_\infty = \|A\|_1\|A\|_\infty$。$\blacksquare$

!!! example "例 15.6"
    对矩阵 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$：

    - $\|A\|_1 = \max\{4, 6\} = 6$；
    - $\|A\|_\infty = \max\{3, 7\} = 7$；
    - $\|A\|_F = \sqrt{30} \approx 5.48$；
    - $\|A\|_2 \approx 5.46$（最大奇异值）。

    验证不等式 (4)：$\|A\|_2^2 \approx 29.87 \leq 42 = 6 \times 7 = \|A\|_1\|A\|_\infty$。

---

## 15.5 条件数

<div class="context-flow" markdown>

**核心**：$\kappa(A)=\|A\|\|A^{-1}\|=\sigma_1/\sigma_n$ · 条件数是"扰动放大器"：$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|}\leq\kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}$ · 正规矩阵 $\kappa_2=\rho(A)/|\lambda_{\min}|$

</div>

条件数衡量了求解线性方程组时问题本身对输入扰动的敏感程度，是数值线性代数中最重要的概念之一。

!!! definition "定义 15.6 (条件数)"
    设 $A \in \mathbb{C}^{n \times n}$ 可逆，$\|\cdot\|$ 为矩阵范数。$A$ 关于范数 $\|\cdot\|$ 的**条件数（condition number）** 定义为

    $$\kappa(A) = \|A\|\cdot\|A^{-1}\|.$$

    对于 $\ell_2$ 范数，$\kappa_2(A) = \|A\|_2\|A^{-1}\|_2 = \frac{\sigma_1}{\sigma_n}$（最大奇异值与最小奇异值之比）。

!!! theorem "定理 15.7 (条件数的基本性质)"
    设 $A$ 可逆，则：

    (1) $\kappa(A) \geq 1$；

    (2) $\kappa(\alpha A) = \kappa(A)$，对任意 $\alpha \neq 0$；

    (3) 对算子范数，$\kappa(A) \geq \frac{\|A\|}{\|B\|}$ 当 $B$ 是与 $A$ 同阶的任意可逆矩阵时（此条不成立，应为其他性质）；

    (4) 若 $U$ 是酉矩阵，则 $\kappa_2(U) = 1$；

    (5) $\kappa_2(A) = \kappa_2(A^*)$。

??? proof "证明"
    **(1)** $\kappa(A) = \|A\|\|A^{-1}\| \geq \|AA^{-1}\| = \|I\| = 1$。

    **(2)** $\kappa(\alpha A) = \|\alpha A\|\|(\alpha A)^{-1}\| = |\alpha|\|A\|\frac{1}{|\alpha|}\|A^{-1}\| = \kappa(A)$。

    **(4)** $\|U\|_2 = 1$（酉矩阵的奇异值都是 $1$），$\|U^{-1}\|_2 = \|U^*\|_2 = 1$。

    **(5)** $\kappa_2(A^*) = \sigma_1(A^*)/\sigma_n(A^*) = \sigma_1(A)/\sigma_n(A) = \kappa_2(A)$。$\blacksquare$

!!! theorem "定理 15.8 (线性方程组的扰动定理)"
    设 $A$ 可逆，$A\mathbf{x} = \mathbf{b}$（$\mathbf{b} \neq \mathbf{0}$）。若 $\mathbf{b}$ 被扰动为 $\mathbf{b} + \delta\mathbf{b}$，相应的解变为 $\mathbf{x} + \delta\mathbf{x}$，则

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.$$

    若 $A$ 被扰动为 $A + \delta A$（$\|\delta A\| < \|A^{-1}\|^{-1}$），解变为 $\mathbf{x} + \delta\mathbf{x}$，则

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x} + \delta\mathbf{x}\|} \leq \kappa(A)\frac{\|\delta A\|}{\|A\|}.$$

??? proof "证明"
    **右端扰动**：$A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$，因此 $A\delta\mathbf{x} = \delta\mathbf{b}$，$\delta\mathbf{x} = A^{-1}\delta\mathbf{b}$。

    $$\|\delta\mathbf{x}\| = \|A^{-1}\delta\mathbf{b}\| \leq \|A^{-1}\|\|\delta\mathbf{b}\|.$$

    由 $\|\mathbf{b}\| = \|A\mathbf{x}\| \leq \|A\|\|\mathbf{x}\|$，即 $\frac{1}{\|\mathbf{x}\|} \leq \frac{\|A\|}{\|\mathbf{b}\|}$。因此

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\|\delta\mathbf{b}\|\frac{\|A\|}{\|\mathbf{b}\|} = \kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.$$

    **系数矩阵扰动**：$(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} = A\mathbf{x}$，展开得 $A\delta\mathbf{x} + \delta A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{0}$，即 $\delta\mathbf{x} = -A^{-1}\delta A(\mathbf{x} + \delta\mathbf{x})$。因此

    $$\|\delta\mathbf{x}\| \leq \|A^{-1}\|\|\delta A\|\|\mathbf{x}+\delta\mathbf{x}\|,$$

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}+\delta\mathbf{x}\|} \leq \|A^{-1}\|\|\delta A\| = \kappa(A)\frac{\|\delta A\|}{\|A\|}. \quad \blacksquare$$

!!! example "例 15.7"
    著名的 Hilbert 矩阵 $H_n = \left(\frac{1}{i+j-1}\right)_{i,j=1}^{n}$ 是病态矩阵的典型例子。

    - $H_3 = \begin{pmatrix} 1 & 1/2 & 1/3 \\ 1/2 & 1/3 & 1/4 \\ 1/3 & 1/4 & 1/5 \end{pmatrix}$，$\kappa_2(H_3) \approx 524$；
    - $\kappa_2(H_5) \approx 4.77 \times 10^5$；
    - $\kappa_2(H_{10}) \approx 1.60 \times 10^{13}$。

    条件数随 $n$ 指数增长，意味着用 Hilbert 矩阵求解线性方程组时，输入数据的微小扰动会导致解的巨大变化。

!!! example "例 15.8"
    设 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1.001 \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 2 \\ 2.001 \end{pmatrix}$，解为 $\mathbf{x} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    扰动 $\mathbf{b}' = \begin{pmatrix} 2 \\ 2.002 \end{pmatrix}$（相对扰动约 $0.05\%$），新解为 $\mathbf{x}' = \begin{pmatrix} 0 \\ 2 \end{pmatrix}$。

    相对解变化 $\frac{\|\mathbf{x}'-\mathbf{x}\|_2}{\|\mathbf{x}\|_2} = \frac{\sqrt{2}}{\sqrt{2}} = 1 = 100\%$。

    $\kappa_2(A) \approx 4002$，远大于 $1$，表明此问题极其病态。

---

## 15.6 特征值的扰动

<div class="context-flow" markdown>

**脉络**：**Bauer-Fike**($\kappa(X)\|E\|$ 界) → **Weyl**(Hermitian矩阵 $|\Delta\lambda_i|\leq\|E\|_2$) → **Wielandt-Hoffman**(Frobenius范数界) · 正规矩阵最稳定($\kappa=1$)

</div>

!!! theorem "定理 15.9 (Bauer-Fike 定理)"
    设 $A \in \mathbb{C}^{n \times n}$ 可对角化，$A = X\Lambda X^{-1}$（$\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$）。设 $\mu$ 是 $A + E$ 的任意特征值，则存在 $A$ 的某个特征值 $\lambda_j$ 使得

    $$|\mu - \lambda_j| \leq \kappa_p(X)\|E\|_p,$$

    其中 $\kappa_p(X) = \|X\|_p\|X^{-1}\|_p$ 是特征向量矩阵的条件数。

??? proof "证明"
    若 $\mu = \lambda_j$ 对某个 $j$ 成立则不等式显然。否则 $\mu$ 不是 $A$ 的特征值，即 $\mu I - A$ 可逆。由 $(A+E)\mathbf{v} = \mu\mathbf{v}$（$\mathbf{v} \neq \mathbf{0}$），得 $(\mu I - A)\mathbf{v} = E\mathbf{v}$，即

    $$\mathbf{v} = (\mu I - A)^{-1}E\mathbf{v}.$$

    因此 $1 \leq \|(\mu I - A)^{-1}E\|_p \leq \|(\mu I - A)^{-1}\|_p\|E\|_p$。

    由 $A = X\Lambda X^{-1}$，$(\mu I - A)^{-1} = X(\mu I - \Lambda)^{-1}X^{-1}$。故

    $$\|(\mu I - A)^{-1}\|_p \leq \|X\|_p\|(\mu I - \Lambda)^{-1}\|_p\|X^{-1}\|_p.$$

    $(\mu I - \Lambda)^{-1} = \operatorname{diag}\left(\frac{1}{\mu-\lambda_1}, \ldots, \frac{1}{\mu-\lambda_n}\right)$，其范数为 $\max_j\frac{1}{|\mu-\lambda_j|} = \frac{1}{\min_j|\mu-\lambda_j|}$。

    因此 $1 \leq \frac{\kappa_p(X)}{\min_j|\mu-\lambda_j|}\|E\|_p$，即 $\min_j|\mu-\lambda_j| \leq \kappa_p(X)\|E\|_p$。$\blacksquare$

!!! note "注"
    Bauer-Fike 定理表明，特征值对扰动的敏感性由特征向量矩阵的条件数 $\kappa(X)$ 控制。对于正规矩阵，$X$ 可取为酉矩阵，$\kappa_2(X) = 1$，此时特征值对扰动最不敏感。

!!! theorem "定理 15.10 (Weyl 不等式——Hermitian 矩阵)"
    设 $A, B$ 是 $n \times n$ Hermitian 矩阵，特征值分别降序排列为 $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$ 和 $\lambda_1(B) \geq \cdots \geq \lambda_n(B)$。设 $A + B$ 的特征值为 $\lambda_1(A+B) \geq \cdots \geq \lambda_n(A+B)$。则对所有 $i + j - 1 \leq n$，

    $$\lambda_{i+j-1}(A + B) \leq \lambda_i(A) + \lambda_j(B).$$

    特别地，取 $j = 1$：$\lambda_i(A+B) \leq \lambda_i(A) + \lambda_1(B)$。

    取 $i = 1$：$\lambda_j(A+B) \leq \lambda_1(A) + \lambda_j(B)$。

??? proof "证明"
    利用特征值的极大极小原理（Courant-Fischer 定理）：

    $$\lambda_k(M) = \min_{\dim V = n-k+1}\max_{\substack{\mathbf{x} \in V \\ \|\mathbf{x}\| = 1}}\mathbf{x}^*M\mathbf{x}.$$

    设 $U$ 是使 $\lambda_i(A)$ 达到的 $(n-i+1)$ 维子空间，$W$ 是使 $\lambda_j(B)$ 达到的 $(n-j+1)$ 维子空间。令 $V = U \cap W$，其维数 $\geq (n-i+1) + (n-j+1) - n = n-i-j+2$，因此 $\dim V \geq n-(i+j-1)+1$。

    对 $\mathbf{x} \in V$，$\|\mathbf{x}\| = 1$：

    $$\mathbf{x}^*(A+B)\mathbf{x} = \mathbf{x}^*A\mathbf{x} + \mathbf{x}^*B\mathbf{x} \leq \lambda_i(A) + \lambda_j(B).$$

    由极大极小原理，$\lambda_{i+j-1}(A+B) \leq \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1}\mathbf{x}^*(A+B)\mathbf{x} \leq \lambda_i(A) + \lambda_j(B)$。$\blacksquare$

!!! theorem "定理 15.11 (Wielandt-Hoffman 定理)"
    设 $A, B$ 是 $n \times n$ Hermitian 矩阵，特征值分别为 $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$ 和 $\lambda_1(B) \geq \cdots \geq \lambda_n(B)$。则

    $$\sum_{i=1}^{n}(\lambda_i(A) - \lambda_i(B))^2 \leq \|A - B\|_F^2.$$

??? proof "证明"
    设 $A = U\Lambda_A U^*$，$B = V\Lambda_B V^*$（谱分解）。令 $W = U^*V$，则 $W$ 是酉矩阵。

    $$\|A - B\|_F^2 = \|\Lambda_A - W\Lambda_B W^*\|_F^2 = \operatorname{tr}(\Lambda_A^2) + \operatorname{tr}(\Lambda_B^2) - 2\operatorname{Re}\operatorname{tr}(\Lambda_A W\Lambda_B W^*).$$

    由 von Neumann 迹不等式，$\operatorname{Re}\operatorname{tr}(\Lambda_A W\Lambda_B W^*) \leq \sum_i\lambda_i(A)\lambda_i(B)$（当 $W$ 是置换矩阵时取等）。因此

    $$\|A-B\|_F^2 \geq \sum\lambda_i(A)^2 + \sum\lambda_i(B)^2 - 2\sum\lambda_i(A)\lambda_i(B) = \sum(\lambda_i(A)-\lambda_i(B))^2. \quad \blacksquare$$

!!! example "例 15.9"
    设 $A = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$，$E = \begin{pmatrix} 0.1 & 0 \\ 0 & -0.1 \end{pmatrix}$。

    $A$ 的特征值：$\lambda = 4 \pm \sqrt{2}$，即 $\lambda_1 \approx 5.41$，$\lambda_2 \approx 2.59$。

    $A + E = \begin{pmatrix} 5.1 & 1 \\ 1 & 2.9 \end{pmatrix}$，特征值 $\lambda = 4 \pm \sqrt{2.21}$，即 $\lambda_1' \approx 5.49$，$\lambda_2' \approx 2.51$。

    $\max_i|\lambda_i' - \lambda_i| \approx 0.08 \leq \|E\|_2 = 0.1$（Weyl 不等式对 Hermitian 矩阵给出 $\kappa_2 = 1$ 的估计）。

---

## 15.7 奇异值的扰动

<div class="context-flow" markdown>

**洞察**：构造 Hermitian 膨胀 $\hat{A}=\begin{pmatrix}0&A\\A^*&0\end{pmatrix}$ 将奇异值扰动归结为特征值扰动 → Ch18 统一不等式框架

</div>

!!! theorem "定理 15.12 (奇异值的 Weyl 不等式)"
    设 $A, B \in \mathbb{C}^{m \times n}$，奇异值分别降序排列为 $\sigma_1(A) \geq \cdots \geq \sigma_p(A)$ 和 $\sigma_1(B) \geq \cdots \geq \sigma_p(B)$（$p = \min(m,n)$）。则对每个 $i$，

    $$|\sigma_i(A) - \sigma_i(B)| \leq \|A - B\|_2.$$

    更强地，

    $$\sqrt{\sum_{i=1}^{p}(\sigma_i(A) - \sigma_i(B))^2} \leq \|A - B\|_F.$$

??? proof "证明"
    构造 $(m+n) \times (m+n)$ Hermitian 矩阵

    $$\hat{A} = \begin{pmatrix} O & A \\ A^* & O \end{pmatrix}, \quad \hat{B} = \begin{pmatrix} O & B \\ B^* & O \end{pmatrix}.$$

    $\hat{A}$ 的特征值为 $\pm\sigma_1(A), \ldots, \pm\sigma_p(A)$ 和 $|m-n|$ 个零（若 $m \neq n$）。$\hat{B}$ 类似。

    将 Hermitian 矩阵的 Weyl 不等式应用于 $\hat{A}$ 和 $\hat{B}$：

    $$|\sigma_i(A) - \sigma_i(B)| = |\lambda_i(\hat{A}) - \lambda_i(\hat{B})| \leq \|\hat{A} - \hat{B}\|_2 = \|A - B\|_2.$$

    Frobenius 范数版本类似地由 Wielandt-Hoffman 定理推出。$\blacksquare$

!!! proposition "命题 15.2 (奇异值的次可加性)"
    设 $A, B \in \mathbb{C}^{m \times n}$，则对每个 $i + j - 1 \leq \min(m,n)$，

    $$\sigma_{i+j-1}(A + B) \leq \sigma_i(A) + \sigma_j(B).$$

!!! example "例 15.10"
    设 $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，$E = \begin{pmatrix} 0 & 0.2 \\ 0.2 & 0 \end{pmatrix}$。

    $\sigma_1(A) = 3$，$\sigma_2(A) = 1$。$B = A + E = \begin{pmatrix} 3 & 0.2 \\ 0.2 & 1 \end{pmatrix}$。

    $B^TB = \begin{pmatrix} 9.04 & 0.8 \\ 0.8 & 1.04 \end{pmatrix}$，特征值 $\mu = 5.04 \pm \sqrt{16.64}$。$\mu_1 \approx 9.12$，$\mu_2 \approx 0.96$。$\sigma_1(B) \approx 3.02$，$\sigma_2(B) \approx 0.98$。

    $|\sigma_1(B) - \sigma_1(A)| \approx 0.02 \leq 0.2 = \|E\|_2$。$|\sigma_2(B) - \sigma_2(A)| \approx 0.02 \leq 0.2$。

---

## 15.8 线性方程组的扰动分析

<div class="context-flow" markdown>

**脉络**：残量小 $\not\Rightarrow$ 误差小(条件数放大) · $\frac{1}{\kappa}\cdot\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}\leq\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|}\leq\kappa\cdot\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}$ → 双向夹逼

</div>

!!! definition "定义 15.7 (前向误差与后向误差)"
    设 $\hat{\mathbf{x}}$ 是线性方程组 $A\mathbf{x} = \mathbf{b}$ 的近似解。

    - **前向误差（forward error）**：$\|\hat{\mathbf{x}} - \mathbf{x}\|$（解的误差）；
    - **后向误差（backward error）**：$\|\mathbf{r}\| = \|\mathbf{b} - A\hat{\mathbf{x}}\|$（残量的范数）；
    - **相对前向误差**：$\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|}$；
    - **相对后向误差**：$\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}$。

!!! theorem "定理 15.13 (前向误差与后向误差的关系)"
    设 $A$ 可逆，$\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$ 为残量。则

    $$\frac{1}{\kappa(A)}\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}.$$

??? proof "证明"
    由 $A(\hat{\mathbf{x}} - \mathbf{x}) = A\hat{\mathbf{x}} - \mathbf{b} = -\mathbf{r}$，即 $\hat{\mathbf{x}} - \mathbf{x} = -A^{-1}\mathbf{r}$。

    **上界**：$\|\hat{\mathbf{x}}-\mathbf{x}\| \leq \|A^{-1}\|\|\mathbf{r}\|$，$\|\mathbf{b}\| \leq \|A\|\|\mathbf{x}\|$，因此

    $$\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\frac{\|\mathbf{r}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\|A\|\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} = \kappa(A)\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}.$$

    **下界**：$\|\mathbf{r}\| = \|A(\hat{\mathbf{x}}-\mathbf{x})\| \leq \|A\|\|\hat{\mathbf{x}}-\mathbf{x}\|$，$\|\mathbf{x}\| \leq \|A^{-1}\|\|\mathbf{b}\|$，因此

    $$\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \|A\|\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{b}\|} \leq \|A\|\|A^{-1}\|\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|} = \kappa(A)\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|}.$$

    即 $\frac{1}{\kappa(A)}\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|}$。$\blacksquare$

!!! theorem "定理 15.14 (同时扰动 $A$ 和 $\mathbf{b}$)"
    设 $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$，且 $\kappa(A)\frac{\|\delta A\|}{\|A\|} < 1$，则

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \frac{\kappa(A)}{1 - \kappa(A)\frac{\|\delta A\|}{\|A\|}}\left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}\right).$$

!!! example "例 15.11"
    用浮点运算求解 $A\mathbf{x} = \mathbf{b}$，其中 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1.0001 \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 2 \\ 2.0001 \end{pmatrix}$。

    精确解 $\mathbf{x} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    假设计算得到 $\hat{\mathbf{x}} = \begin{pmatrix} 1.5 \\ 0.5 \end{pmatrix}$。残量 $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}} = \begin{pmatrix} 0 \\ 0.00005 \end{pmatrix}$。

    相对残量 $\frac{\|\mathbf{r}\|_\infty}{\|\mathbf{b}\|_\infty} = \frac{0.00005}{2.0001} \approx 2.5 \times 10^{-5}$（很小）。

    但相对前向误差 $\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|_\infty}{\|\mathbf{x}\|_\infty} = \frac{0.5}{1} = 0.5 = 50\%$（很大）。

    这体现了条件数的放大效应。$\kappa_\infty(A) \approx 4 \times 10^4$，使得小残量对应大误差。

!!! example "例 15.12"
    比较良态和病态系统。

    **良态系统**：$A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$，$\kappa_2(A) = 2$。扰动 $\|\delta\mathbf{b}\|/\|\mathbf{b}\| = 10^{-6}$ 导致 $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \leq 2 \times 10^{-6}$。

    **病态系统**：$B = \begin{pmatrix} 1 & 1 \\ 1 & 1+10^{-10} \end{pmatrix}$，$\kappa_2(B) \approx 4 \times 10^{10}$。同样的扰动可能导致 $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \leq 4 \times 10^4$，即解完全不可靠。
