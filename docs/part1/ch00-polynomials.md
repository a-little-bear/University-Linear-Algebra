# 第 0 章 多项式代数

<div class="context-flow" markdown>

**前置**：实数域、复数域、有理数域的基本概念

**本章脉络**：多项式环 → 整除与 GCD → 互素 → 不可约多项式 → 唯一分解 → 根 → 重根与判别式 → $\mathbb{R}$ 和 $\mathbb{C}$ 上的不可约多项式 → $\mathbb{Q}$ 上的不可约判定 → 多项式插值

**延伸**：多项式代数延伸到形式幂级数环 $\mathbb{F}[[x]]$；代数几何中多项式理想与仿射簇的对应（Hilbert 零点定理）；编码理论中 $\mathbb{F}_q[x]$ 上多项式用于 Reed-Solomon 码的构造

</div>

多项式是代数学中最基本的对象之一。在线性代数中，多项式扮演着核心角色：特征多项式决定特征值，最小多项式刻画线性变换的本质结构，而 Hamilton-Cayley 定理将矩阵与多项式紧密联系。本章系统建立一元多项式的代数理论，为后续章节奠定坚实基础。

---

## 0.1 一元多项式环

<div class="context-flow" markdown>

**核心问题**：如何严格定义多项式？ → 区分多项式（形式表达式）与多项式函数 → 多项式全体构成环

</div>

### 多项式的定义

设 $\mathbb{F}$ 为一个域（本课程中通常取 $\mathbb{F} = \mathbb{Q}, \mathbb{R}$ 或 $\mathbb{C}$）。

!!! definition "定义 0.1 (一元多项式)"
    设 $\mathbb{F}$ 为域，$x$ 为一个不定元（indeterminate）。形如

    $$
    f(x) = a_n x^n + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0
    $$

    的表达式称为 $\mathbb{F}$ 上的一个**一元多项式**（polynomial in one variable），其中 $a_0, a_1, \ldots, a_n \in \mathbb{F}$ 称为**系数**（coefficients）。若 $a_n \neq 0$，则 $n$ 称为 $f(x)$ 的**次数**（degree），记作 $\deg f = n$，$a_n$ 称为**首项系数**（leading coefficient）。若 $a_n = 1$，则称 $f(x)$ 为**首一多项式**（monic polynomial）。零多项式 $f(x) = 0$ 的次数约定为 $-\infty$。

!!! note "注"
    多项式是**形式表达式**，不是函数。两个多项式相等当且仅当对应系数全部相等。在无限域上（如 $\mathbb{Q}, \mathbb{R}, \mathbb{C}$），多项式与多项式函数可以等同看待；但在有限域上（如 $\mathbb{F}_p$），多项式 $x^p - x$ 作为函数恒为零，但作为多项式不是零多项式。

!!! definition "定义 0.2 (多项式环)"
    $\mathbb{F}$ 上所有一元多项式的集合，配上通常的加法和乘法，记作 $\mathbb{F}[x]$，称为 $\mathbb{F}$ 上的**一元多项式环**（polynomial ring）。

### 多项式的运算

设 $f(x) = \sum_{i=0}^n a_i x^i$，$g(x) = \sum_{j=0}^m b_j x^j$。

**加法**：$f(x) + g(x) = \sum_{k=0}^{\max(n,m)} (a_k + b_k) x^k$（不足的高次系数补零）。

**乘法**：$f(x) \cdot g(x) = \sum_{k=0}^{n+m} c_k x^k$，其中 $c_k = \sum_{i+j=k} a_i b_j$。

!!! theorem "定理 0.1 (次数公式)"
    设 $f(x), g(x) \in \mathbb{F}[x]$ 均非零多项式，则

    $$
    \deg(f \cdot g) = \deg f + \deg g.
    $$

    特别地，$\mathbb{F}[x]$ 是**整环**（integral domain）：若 $f(x) \cdot g(x) = 0$，则 $f(x) = 0$ 或 $g(x) = 0$。

??? proof "证明"
    设 $\deg f = n$，$\deg g = m$，首项系数分别为 $a_n \neq 0$ 和 $b_m \neq 0$。乘积 $f \cdot g$ 的最高次项为 $a_n b_m x^{n+m}$。由于 $\mathbb{F}$ 是域（从而是整环），$a_n b_m \neq 0$，故 $\deg(f \cdot g) = n + m$。

    若 $f \cdot g = 0$，则 $\deg(f \cdot g) = -\infty$。若 $f, g$ 均非零，则 $\deg(f \cdot g) = \deg f + \deg g \geq 0$，矛盾。$\blacksquare$

!!! theorem "定理 0.2 (带余除法)"
    设 $f(x), g(x) \in \mathbb{F}[x]$，$g(x) \neq 0$。则存在唯一的 $q(x), r(x) \in \mathbb{F}[x]$，使得

    $$
    f(x) = q(x) \cdot g(x) + r(x),
    $$

    其中 $\deg r < \deg g$（或 $r = 0$）。$q(x)$ 称为**商**（quotient），$r(x)$ 称为**余式**（remainder）。

??? proof "证明"
    **存在性**：对 $\deg f$ 进行归纳。若 $\deg f < \deg g$，取 $q = 0$，$r = f$。否则设 $\deg f = n \geq m = \deg g$，$f$ 的首项系数为 $a_n$，$g$ 的首项系数为 $b_m$。令

    $$
    f_1(x) = f(x) - \frac{a_n}{b_m} x^{n-m} g(x).
    $$

    则 $\deg f_1 < n$。由归纳假设，$f_1 = q_1 g + r$，其中 $\deg r < \deg g$。于是

    $$
    f = \left(\frac{a_n}{b_m} x^{n-m} + q_1\right) g + r.
    $$

    **唯一性**：设 $f = q_1 g + r_1 = q_2 g + r_2$，其中 $\deg r_1, \deg r_2 < \deg g$。则 $(q_1 - q_2)g = r_2 - r_1$。若 $q_1 \neq q_2$，则左边次数 $\geq \deg g$，而右边次数 $< \deg g$，矛盾。故 $q_1 = q_2$，从而 $r_1 = r_2$。$\blacksquare$

!!! example "例 0.1"
    在 $\mathbb{R}[x]$ 中，设 $f(x) = 2x^4 + 3x^3 - x + 5$，$g(x) = x^2 + 1$。求 $f$ 除以 $g$ 的商和余式。

    进行多项式长除法：

    $$
    2x^4 + 3x^3 - x + 5 = (2x^2 + 3x - 2)(x^2 + 1) + (-4x + 7).
    $$

    因此 $q(x) = 2x^2 + 3x - 2$，$r(x) = -4x + 7$。验证：$\deg r = 1 < 2 = \deg g$。

---

## 0.2 整除与最大公因式

<div class="context-flow" markdown>

**核心问题**：多项式之间的整除关系 → 最大公因式的存在与计算 → 辗转相除法（Euclidean 算法）

</div>

!!! definition "定义 0.3 (整除)"
    设 $f(x), g(x) \in \mathbb{F}[x]$。若存在 $q(x) \in \mathbb{F}[x]$ 使得 $f(x) = q(x) \cdot g(x)$，则称 $g(x)$ **整除** $f(x)$，记作 $g(x) \mid f(x)$。此时称 $g(x)$ 是 $f(x)$ 的一个**因式**（factor/divisor），$f(x)$ 是 $g(x)$ 的一个**倍式**（multiple）。

!!! definition "定义 0.4 (相伴)"
    若 $f(x) \mid g(x)$ 且 $g(x) \mid f(x)$，则称 $f(x)$ 与 $g(x)$ **相伴**（associates），记作 $f(x) \sim g(x)$。在 $\mathbb{F}[x]$ 中，$f \sim g$ 当且仅当 $f = cg$，其中 $c \in \mathbb{F}^*$（$c$ 是非零常数）。

!!! definition "定义 0.5 (最大公因式)"
    设 $f(x), g(x) \in \mathbb{F}[x]$，不全为零。多项式 $d(x) \in \mathbb{F}[x]$ 称为 $f(x)$ 与 $g(x)$ 的一个**最大公因式**（greatest common divisor, GCD），如果：

    1. $d(x) \mid f(x)$ 且 $d(x) \mid g(x)$（$d$ 是公因式）；
    2. 若 $c(x) \mid f(x)$ 且 $c(x) \mid g(x)$，则 $c(x) \mid d(x)$（$d$ 是最大的）。

    $f$ 与 $g$ 的首一最大公因式记作 $\gcd(f, g)$ 或 $(f, g)$。

!!! theorem "定理 0.3 (GCD 的存在性与 Bézout 等式)"
    设 $f(x), g(x) \in \mathbb{F}[x]$ 不全为零。则：

    1. $\gcd(f, g)$ 存在且唯一（在相伴意义下）。
    2. 存在 $u(x), v(x) \in \mathbb{F}[x]$，使得

    $$
    \gcd(f, g) = u(x) f(x) + v(x) g(x).
    $$

    此等式称为 **Bézout 等式**（Bézout's identity）。

??? proof "证明"
    考虑集合 $I = \{u(x)f(x) + v(x)g(x) : u, v \in \mathbb{F}[x]\}$。$I$ 是 $\mathbb{F}[x]$ 的一个理想。取 $I$ 中次数最低的非零多项式 $d(x)$（归一化为首一）。

    对任意 $h(x) \in I$，由带余除法 $h = qd + r$，其中 $\deg r < \deg d$。由于 $h, d \in I$，有 $r = h - qd \in I$。由 $d$ 的最低次性，$r = 0$，故 $d \mid h$。

    特别地，$f, g \in I$，所以 $d \mid f$ 且 $d \mid g$。若 $c \mid f$ 且 $c \mid g$，由 $d = uf + vg$ 知 $c \mid d$。故 $d = \gcd(f, g)$。$\blacksquare$

### 辗转相除法（Euclidean 算法）

!!! theorem "定理 0.4 (辗转相除法)"
    设 $f(x), g(x) \in \mathbb{F}[x]$，$g \neq 0$。反复应用带余除法：

    $$
    \begin{aligned}
    f &= q_1 g + r_1, & \deg r_1 &< \deg g, \\
    g &= q_2 r_1 + r_2, & \deg r_2 &< \deg r_1, \\
    r_1 &= q_3 r_2 + r_3, & \deg r_3 &< \deg r_2, \\
    &\;\;\vdots \\
    r_{k-2} &= q_k r_{k-1} + r_k, & \deg r_k &< \deg r_{k-1}, \\
    r_{k-1} &= q_{k+1} r_k.
    \end{aligned}
    $$

    则 $\gcd(f, g) \sim r_k$（最后一个非零余式）。

??? proof "证明"
    由带余除法 $f = q_1 g + r_1$ 可知，$d \mid f$ 且 $d \mid g$ 当且仅当 $d \mid g$ 且 $d \mid r_1$。因此 $\gcd(f, g) = \gcd(g, r_1) = \gcd(r_1, r_2) = \cdots = \gcd(r_{k-1}, r_k) = r_k$（归一化为首一）。

    由于 $\deg g > \deg r_1 > \deg r_2 > \cdots$ 是严格递减的非负整数序列，算法必在有限步终止。$\blacksquare$

!!! example "例 0.2"
    求 $\gcd(x^4 - 1, x^3 - 1)$（在 $\mathbb{Q}[x]$ 中）。

    $$
    x^4 - 1 = x \cdot (x^3 - 1) + (x - 1),
    $$

    $$
    x^3 - 1 = (x^2 + x + 1)(x - 1) + 0.
    $$

    因此 $\gcd(x^4 - 1, x^3 - 1) = x - 1$。

!!! example "例 0.3"
    求 $\gcd(f, g)$ 并求 Bézout 系数，其中 $f = x^3 + x + 1$，$g = x^2 + x$（在 $\mathbb{Q}[x]$ 中）。

    $$
    x^3 + x + 1 = (x - 1)(x^2 + x) + (2x + 1),
    $$

    $$
    x^2 + x = \left(\frac{1}{2}x + \frac{1}{4}\right)(2x + 1) + \frac{3}{4}.
    $$

    余式 $\frac{3}{4}$ 非零且为常数，故 $\gcd(f, g) = 1$。回代得 Bézout 系数：

    $$
    1 = \frac{4}{3}\left[(x^2 + x) - \left(\frac{1}{2}x + \frac{1}{4}\right)(2x+1)\right] = \frac{4}{3}(x^2+x) - \frac{4}{3}\left(\frac{1}{2}x+\frac{1}{4}\right)\left[(x^3+x+1) - (x-1)(x^2+x)\right].
    $$

    整理后可得 $u(x) f(x) + v(x) g(x) = 1$。

---

## 0.3 互素多项式

<div class="context-flow" markdown>

**核心问题**：$\gcd(f, g) = 1$ 意味着什么？ → 互素的等价条件 → 互素多项式的关键性质

</div>

!!! definition "定义 0.6 (互素)"
    若 $\gcd(f(x), g(x)) = 1$，则称 $f(x)$ 与 $g(x)$ **互素**（coprime / relatively prime）。

!!! theorem "定理 0.5 (互素的等价条件)"
    $f(x)$ 与 $g(x)$ 互素的充要条件是：存在 $u(x), v(x) \in \mathbb{F}[x]$，使得

    $$
    u(x) f(x) + v(x) g(x) = 1.
    $$

??? proof "证明"
    **必要性**：若 $\gcd(f, g) = 1$，由 Bézout 等式直接得到。

    **充分性**：若 $uf + vg = 1$，设 $d = \gcd(f, g)$。由 $d \mid f$ 和 $d \mid g$，知 $d \mid (uf + vg) = 1$，故 $d$ 是常数，$\gcd(f, g) = 1$。$\blacksquare$

!!! theorem "定理 0.6 (互素的乘法性质)"
    设 $f(x), g(x), h(x) \in \mathbb{F}[x]$。

    1. 若 $f \mid gh$ 且 $\gcd(f, g) = 1$，则 $f \mid h$。
    2. 若 $f \mid h$ 且 $g \mid h$，且 $\gcd(f, g) = 1$，则 $fg \mid h$。

??? proof "证明"
    1. 由 $\gcd(f, g) = 1$，存在 $u, v$ 使得 $uf + vg = 1$。两边乘以 $h$：$ufh + vgh = h$。由 $f \mid ufh$ 和 $f \mid gh$（故 $f \mid vgh$），得 $f \mid h$。

    2. 由 $f \mid h$，设 $h = fk$。由 $g \mid h = fk$ 且 $\gcd(f, g) = 1$，由 (1) 知 $g \mid k$，设 $k = gl$，则 $h = fgl$，故 $fg \mid h$。$\blacksquare$

!!! theorem "定理 0.7 (多个多项式互素)"
    设 $f_1, f_2, \ldots, f_k \in \mathbb{F}[x]$ 两两互素，且每个 $f_i \mid h$。则 $f_1 f_2 \cdots f_k \mid h$。

??? proof "证明"
    对 $k$ 归纳。$k = 2$ 即定理 0.6 (2)。设 $k - 1$ 时成立，则 $f_1 \cdots f_{k-1} \mid h$。由 $\gcd(f_i, f_k) = 1$ 对 $i = 1, \ldots, k-1$，可以证明 $\gcd(f_1 \cdots f_{k-1}, f_k) = 1$。再由定理 0.6 (2) 得 $f_1 \cdots f_k \mid h$。$\blacksquare$

!!! example "例 0.4"
    设 $f(x) = x^2 - 1 = (x-1)(x+1)$，$g(x) = x^2 + x = x(x+1)$。则

    $$
    \gcd(f, g) = x + 1 \neq 1,
    $$

    所以 $f$ 与 $g$ 不互素。但 $x - 1$ 与 $x$ 互素（它们的 GCD 为 $1$）。

---

## 0.4 不可约多项式

<div class="context-flow" markdown>

**核心问题**：哪些多项式不能进一步分解？ → 不可约多项式——多项式环中的"素元" → 不可约多项式的基本性质

</div>

!!! definition "定义 0.7 (不可约多项式)"
    设 $p(x) \in \mathbb{F}[x]$ 的次数 $\geq 1$。若 $p(x)$ 不能表示为两个次数都小于 $\deg p$ 的多项式之积，则称 $p(x)$ 是 $\mathbb{F}$ 上的**不可约多项式**（irreducible polynomial）。否则称 $p(x)$ 是**可约的**（reducible）。

    等价地，$p(x)$ 不可约当且仅当：$p(x) = f(x) g(x)$ 蕴含 $f(x)$ 或 $g(x)$ 是常数。

!!! note "注"
    不可约性依赖于底域 $\mathbb{F}$。例如 $x^2 + 1$ 在 $\mathbb{R}$ 上不可约，但在 $\mathbb{C}$ 上可约，因为 $x^2 + 1 = (x + i)(x - i)$。

!!! theorem "定理 0.8 (不可约多项式的性质)"
    设 $p(x) \in \mathbb{F}[x]$ 不可约，$f(x) \in \mathbb{F}[x]$。则要么 $p \mid f$，要么 $\gcd(p, f) = 1$。

??? proof "证明"
    设 $d = \gcd(p, f)$。则 $d \mid p$，即 $p = d \cdot q$。由 $p$ 不可约，$d$ 或 $q$ 为常数。若 $d$ 为常数，则 $\gcd(p, f) = 1$；若 $q$ 为常数，则 $d \sim p$，从而 $p \mid f$。$\blacksquare$

!!! theorem "定理 0.9 (不可约多项式的素性)"
    设 $p(x) \in \mathbb{F}[x]$ 不可约。若 $p \mid f_1 f_2 \cdots f_k$，则存在某个 $i$ 使得 $p \mid f_i$。

??? proof "证明"
    对 $k$ 归纳。$k = 1$ 显然。设 $k \geq 2$ 且 $p \mid f_1 \cdots f_k$。若 $p \mid f_k$，结论成立。否则 $\gcd(p, f_k) = 1$（定理 0.8），由定理 0.6 (1) 得 $p \mid f_1 \cdots f_{k-1}$，由归纳假设即得。$\blacksquare$

!!! example "例 0.5"
    在 $\mathbb{R}[x]$ 中，$x^2 + 1$ 是不可约多项式。若 $(x^2+1) \mid f(x)g(x)$，则 $(x^2+1) \mid f(x)$ 或 $(x^2+1) \mid g(x)$。

    在 $\mathbb{Q}[x]$ 中，$x^2 - 2$ 是不可约的（因为若可约则有有理根 $\pm\sqrt{2}$，但 $\sqrt{2} \notin \mathbb{Q}$）。

---

## 0.5 多项式的因式分解

<div class="context-flow" markdown>

**核心问题**：多项式能否唯一地分解为不可约因式之积？ → 类比整数的算术基本定理 → 唯一分解定理

</div>

!!! theorem "定理 0.10 (唯一分解定理)"
    $\mathbb{F}[x]$ 中每个次数 $\geq 1$ 的多项式 $f(x)$ 都可以分解为

    $$
    f(x) = c \cdot p_1(x)^{e_1} p_2(x)^{e_2} \cdots p_s(x)^{e_s},
    $$

    其中 $c \in \mathbb{F}^*$，$p_1, p_2, \ldots, p_s$ 是两两不相伴的首一不可约多项式，$e_1, e_2, \ldots, e_s$ 为正整数。此分解在因式的排列顺序意义下是唯一的。

??? proof "证明"
    **存在性**：对 $\deg f$ 归纳。$\deg f = 1$ 时 $f$ 本身不可约。设 $\deg f \geq 2$。若 $f$ 不可约，分解即为 $f$ 本身。若 $f$ 可约，则 $f = g \cdot h$，其中 $1 \leq \deg g, \deg h < \deg f$。由归纳假设，$g, h$ 各自有不可约分解，合并即得 $f$ 的不可约分解。

    **唯一性**：设 $f = c \cdot p_1^{e_1} \cdots p_s^{e_s} = c' \cdot q_1^{d_1} \cdots q_t^{d_t}$。由 $p_1 \mid f = c' \cdot q_1^{d_1} \cdots q_t^{d_t}$ 和 $p_1$ 不可约的素性（定理 0.9），$p_1 \mid q_j$ 对某个 $j$。由 $q_j$ 不可约知 $p_1 \sim q_j$。将 $p_1$ 和 $q_j$ 配对后消去，对剩余因式归纳即得 $s = t$，且适当排列后 $p_i \sim q_i$，$e_i = d_i$。$\blacksquare$

!!! example "例 0.6"
    在 $\mathbb{Q}[x]$ 中分解 $f(x) = x^4 - 1$。

    $$
    x^4 - 1 = (x^2 - 1)(x^2 + 1) = (x - 1)(x + 1)(x^2 + 1).
    $$

    在 $\mathbb{Q}[x]$ 中，$x - 1$、$x + 1$、$x^2 + 1$ 均不可约，这是 $f$ 的唯一分解。

    在 $\mathbb{C}[x]$ 中，$x^2 + 1 = (x - i)(x + i)$，所以

    $$
    x^4 - 1 = (x - 1)(x + 1)(x - i)(x + i).
    $$

!!! example "例 0.7"
    在 $\mathbb{R}[x]$ 中分解 $f(x) = x^6 - 1$。

    $$
    x^6 - 1 = (x^3 - 1)(x^3 + 1) = (x-1)(x^2+x+1)(x+1)(x^2-x+1).
    $$

    其中 $x^2 + x + 1$ 和 $x^2 - x + 1$ 的判别式分别为 $-3$ 和 $-3$，均 $< 0$，故在 $\mathbb{R}[x]$ 上不可约。

---

## 0.6 多项式的根

<div class="context-flow" markdown>

**核心问题**：多项式何时有根？ → 余数定理与因式定理 → 根的个数上界 → 根与因式分解

</div>

!!! definition "定义 0.8 (根)"
    设 $f(x) \in \mathbb{F}[x]$，$\alpha \in \mathbb{F}$。若 $f(\alpha) = 0$，则称 $\alpha$ 是 $f(x)$ 的一个**根**（root / zero）。

!!! theorem "定理 0.11 (余数定理)"
    设 $f(x) \in \mathbb{F}[x]$，$\alpha \in \mathbb{F}$。则 $f(x)$ 除以 $(x - \alpha)$ 的余式为 $f(\alpha)$。

??? proof "证明"
    由带余除法，$f(x) = q(x)(x - \alpha) + r$，其中 $r \in \mathbb{F}$（因为 $\deg r < \deg(x - \alpha) = 1$，故 $r$ 是常数）。令 $x = \alpha$：$f(\alpha) = q(\alpha) \cdot 0 + r = r$。$\blacksquare$

!!! theorem "定理 0.12 (因式定理)"
    $\alpha$ 是 $f(x)$ 的根当且仅当 $(x - \alpha) \mid f(x)$。

??? proof "证明"
    由余数定理，$f(x) = q(x)(x - \alpha) + f(\alpha)$。$f(\alpha) = 0$ 当且仅当 $f(x) = q(x)(x-\alpha)$，即 $(x - \alpha) \mid f(x)$。$\blacksquare$

!!! definition "定义 0.9 (根的重数)"
    设 $\alpha$ 是 $f(x)$ 的根。若 $(x - \alpha)^k \mid f(x)$ 但 $(x - \alpha)^{k+1} \nmid f(x)$，则称 $\alpha$ 是 $f(x)$ 的 $k$ **重根**（root of multiplicity $k$）。$k = 1$ 时称为**单根**（simple root），$k \geq 2$ 时称为**重根**（multiple root）。

!!! theorem "定理 0.13 (根的个数)"
    $\mathbb{F}[x]$ 中 $n$ 次多项式至多有 $n$ 个根（计重数）。

??? proof "证明"
    对 $n$ 归纳。$n = 0$ 时 $f$ 为非零常数，无根。设 $n \geq 1$ 且 $f$ 有根 $\alpha_1$（重数 $k_1$），则 $f(x) = (x - \alpha_1)^{k_1} g(x)$，其中 $g(\alpha_1) \neq 0$ 且 $\deg g = n - k_1$。$f$ 的其余根都是 $g$ 的根（由 $(x - \alpha_1)^{k_1}$ 在 $\alpha \neq \alpha_1$ 处非零）。由归纳假设，$g$ 至多有 $n - k_1$ 个根（计重数），故 $f$ 至多有 $k_1 + (n - k_1) = n$ 个根。$\blacksquare$

!!! example "例 0.8"
    $f(x) = x^3 - 3x + 2 = (x-1)^2(x+2)$。

    根为 $x = 1$（2 重根）和 $x = -2$（单根），共 $3$ 个根（计重数），等于 $\deg f = 3$。

!!! example "例 0.9"
    利用因式定理证明 Vandermonde 行列式公式的一个特例。

    设 $f(x) = \det \begin{pmatrix} 1 & 1 & 1 \\ a & b & x \\ a^2 & b^2 & x^2 \end{pmatrix}$。将 $f(x)$ 按第三列展开，$f(x)$ 是 $x$ 的二次多项式。

    $f(a) = 0$（第一列与第三列相同），$f(b) = 0$（第二列与第三列相同）。故 $(x - a)(x - b) \mid f(x)$，但 $\deg f = 2$，所以 $f(x) = c(x-a)(x-b)$。比较 $x^2$ 的系数得 $c = a - b$（计算略），从而

    $$
    f(x) = (a - b)(x - a)(x - b).
    $$

---

## 0.7 重根与判别式

<div class="context-flow" markdown>

**核心问题**：如何判断多项式是否有重根？ → 形式导数 → 重根判别准则 → 判别式

</div>

!!! definition "定义 0.10 (形式导数)"
    设 $f(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_1 x + a_0 \in \mathbb{F}[x]$。$f(x)$ 的**形式导数**（formal derivative）定义为

    $$
    f'(x) = n a_n x^{n-1} + (n-1)a_{n-1}x^{n-2} + \cdots + a_1.
    $$

    这是一个纯代数定义，不依赖极限概念。形式导数满足通常的求导法则：$(f + g)' = f' + g'$，$(fg)' = f'g + fg'$。

!!! theorem "定理 0.14 (重根判别准则)"
    设 $f(x) \in \mathbb{F}[x]$，$\deg f \geq 1$，$\alpha \in \overline{\mathbb{F}}$（$\mathbb{F}$ 的代数闭包）。则 $\alpha$ 是 $f(x)$ 的 $k$ 重根（$k \geq 2$）当且仅当 $f(\alpha) = 0$ 且 $f'(\alpha) = 0$。

    更精确地，$\alpha$ 是 $f(x)$ 的 $k$ 重根当且仅当 $f(\alpha) = f'(\alpha) = \cdots = f^{(k-1)}(\alpha) = 0$ 且 $f^{(k)}(\alpha) \neq 0$。

??? proof "证明"
    设 $f(x) = (x - \alpha)^k g(x)$，$g(\alpha) \neq 0$，$k \geq 1$。则

    $$
    f'(x) = k(x - \alpha)^{k-1} g(x) + (x - \alpha)^k g'(x) = (x - \alpha)^{k-1}[k g(x) + (x - \alpha)g'(x)].
    $$

    若 $k \geq 2$：$f'(\alpha) = 0$，即 $(x - \alpha) \mid f'(x)$。

    若 $k = 1$：$f'(\alpha) = 1 \cdot g(\alpha) \neq 0$。

    反之，若 $f(\alpha) = 0$ 且 $f'(\alpha) = 0$，则 $\alpha$ 是 $f$ 的根（重数 $\geq 1$），由上述计算知重数 $\geq 2$。$\blacksquare$

!!! theorem "定理 0.15 (无重根的充要条件)"
    $f(x)$ 没有重根（在 $\overline{\mathbb{F}}$ 中）的充要条件是 $\gcd(f, f') = 1$。

??? proof "证明"
    $f$ 有重根 $\alpha$（$k \geq 2$） $\Leftrightarrow$ $(x - \alpha) \mid f$ 且 $(x - \alpha) \mid f'$ $\Leftrightarrow$ $(x - \alpha) \mid \gcd(f, f')$ $\Leftrightarrow$ $\gcd(f, f') \neq 1$。$\blacksquare$

!!! definition "定义 0.11 (判别式)"
    设 $f(x) = a_n x^n + \cdots + a_0 \in \mathbb{F}[x]$，根（在代数闭包中）为 $\alpha_1, \ldots, \alpha_n$。$f(x)$ 的**判别式**（discriminant）定义为

    $$
    \Delta(f) = a_n^{2n-2} \prod_{i < j} (\alpha_i - \alpha_j)^2.
    $$

    $f$ 有重根当且仅当 $\Delta(f) = 0$。

对于二次多项式 $f(x) = ax^2 + bx + c$：

$$
\Delta = b^2 - 4ac.
$$

!!! example "例 0.10"
    判断 $f(x) = x^3 - 3x + 2$ 是否有重根。

    $f'(x) = 3x^2 - 3 = 3(x-1)(x+1)$。

    $\gcd(f, f') = \gcd(x^3 - 3x + 2, 3x^2 - 3)$。由辗转相除法：

    $$
    x^3 - 3x + 2 = \frac{1}{3}x \cdot (3x^2 - 3) + (-2x + 2),
    $$

    $$
    3x^2 - 3 = \left(-\frac{3}{2}x - \frac{3}{2}\right)(-2x + 2) + 0.
    $$

    $\gcd(f, f') = x - 1 \neq 1$，所以 $f$ 有重根 $x = 1$。验证：$f(x) = (x-1)^2(x+2)$。

!!! example "例 0.11"
    对三次多项式 $f(x) = x^3 + px + q$，其判别式为

    $$
    \Delta = -4p^3 - 27q^2.
    $$

    - $\Delta > 0$：三个不同实根。
    - $\Delta = 0$：有重根。
    - $\Delta < 0$：一个实根，两个共轭复根。

---

## 0.8 实数域与复数域上的不可约多项式

<div class="context-flow" markdown>

**核心问题**：$\mathbb{C}[x]$ 和 $\mathbb{R}[x]$ 中哪些多项式不可约？ → 代数基本定理 → $\mathbb{C}$ 上只有一次不可约 → $\mathbb{R}$ 上只有一次和某些二次不可约

</div>

!!! theorem "定理 0.16 (代数基本定理)"
    $\mathbb{C}[x]$ 中每个次数 $\geq 1$ 的多项式都至少有一个复数根。

!!! note "注"
    代数基本定理的证明需要分析学工具（如 Liouville 定理、复分析中的最大模原理，或拓扑学中的基本群论证），超出了纯代数的范畴。我们在此承认其结论。

!!! theorem "定理 0.17 ($\\mathbb{C}$ 上的不可约多项式)"
    $\mathbb{C}[x]$ 中的不可约多项式恰好是一次多项式。因此每个 $f(x) \in \mathbb{C}[x]$（$\deg f = n \geq 1$）都可以分解为

    $$
    f(x) = a_n(x - \alpha_1)(x - \alpha_2) \cdots (x - \alpha_n),
    $$

    其中 $\alpha_1, \ldots, \alpha_n \in \mathbb{C}$ 是 $f$ 的全部根（计重数），$a_n$ 是首项系数。

??? proof "证明"
    设 $p(x) \in \mathbb{C}[x]$ 不可约且 $\deg p \geq 2$。由代数基本定理，$p$ 有根 $\alpha \in \mathbb{C}$，故 $(x - \alpha) \mid p$，这与 $p$ 不可约矛盾。因此不可约多项式的次数只能为 $1$。

    分解式由唯一分解定理和每个不可约因式为一次式推得。$\blacksquare$

!!! theorem "定理 0.18 (共轭根定理)"
    设 $f(x) \in \mathbb{R}[x]$。若 $\alpha = a + bi \in \mathbb{C}$（$b \neq 0$）是 $f$ 的根，则 $\overline{\alpha} = a - bi$ 也是 $f$ 的根，且 $\alpha$ 与 $\overline{\alpha}$ 的重数相同。

??? proof "证明"
    由 $f(\alpha) = 0$ 取共轭：$\overline{f(\alpha)} = f(\overline{\alpha}) = 0$（利用 $f$ 的系数为实数，$\overline{a_k \alpha^k} = a_k \overline{\alpha}^k$）。

    对于重数，设 $f(x) = (x - \alpha)^k g(x)$，$g(\alpha) \neq 0$。取共轭得 $f(x) = (x - \overline{\alpha})^k \overline{g}(x)$（其中 $\overline{g}$ 表示将 $g$ 的系数取共轭）。由 $f$ 系数为实数，$\overline{g}(x)$ 的系数也为实数，且 $\overline{g}(\overline{\alpha}) = \overline{g(\alpha)} \neq 0$。$\blacksquare$

!!! theorem "定理 0.19 ($\\mathbb{R}$ 上的不可约多项式)"
    $\mathbb{R}[x]$ 中的不可约多项式恰好是：

    1. 一次多项式 $ax + b$（$a \neq 0$）；
    2. 判别式 $< 0$ 的二次多项式 $ax^2 + bx + c$（$a \neq 0$，$b^2 - 4ac < 0$）。

??? proof "证明"
    一次多项式显然不可约。判别式 $< 0$ 的二次多项式在 $\mathbb{R}$ 中无根，故不能有一次因式，因而不可约。

    设 $p(x) \in \mathbb{R}[x]$ 不可约，$\deg p \geq 3$。$p$ 在 $\mathbb{C}$ 中有根 $\alpha$。若 $\alpha \in \mathbb{R}$，则 $(x - \alpha) \mid p$，与不可约矛盾。若 $\alpha = a + bi$（$b \neq 0$），则 $\overline{\alpha}$ 也是 $p$ 的根，$(x - \alpha)(x - \overline{\alpha}) = x^2 - 2ax + a^2 + b^2 \in \mathbb{R}[x]$，且此二次式整除 $p$。由 $\deg p \geq 3$，$p$ 有真因式分解，矛盾。故 $\deg p \leq 2$。

    $\deg p = 2$ 时，若判别式 $\geq 0$ 则有实根而可约。$\blacksquare$

!!! example "例 0.12"
    在 $\mathbb{R}[x]$ 中分解 $f(x) = x^4 + 4$。

    $f$ 在 $\mathbb{R}$ 中无根（$x^4 + 4 > 0$），但 $f$ 是可约的：

    $$
    x^4 + 4 = x^4 + 4x^2 + 4 - 4x^2 = (x^2 + 2)^2 - (2x)^2 = (x^2 + 2x + 2)(x^2 - 2x + 2).
    $$

    $x^2 + 2x + 2$ 和 $x^2 - 2x + 2$ 的判别式分别为 $4 - 8 = -4 < 0$，故均不可约。

---

## 0.9 有理数域上的不可约多项式

<div class="context-flow" markdown>

**核心问题**：$\mathbb{Q}[x]$ 中的不可约判定比 $\mathbb{R}$ 和 $\mathbb{C}$ 困难得多 → 有理根判别法 → Gauss 引理 → Eisenstein 判别法

</div>

!!! theorem "定理 0.20 (有理根定理)"
    设 $f(x) = a_n x^n + \cdots + a_1 x + a_0 \in \mathbb{Z}[x]$，$a_n \neq 0$。若 $\frac{p}{q}$（$p, q$ 互素整数，$q > 0$）是 $f$ 的有理根，则 $p \mid a_0$ 且 $q \mid a_n$。

??? proof "证明"
    由 $f(p/q) = 0$，两边乘以 $q^n$：

    $$
    a_n p^n + a_{n-1} p^{n-1} q + \cdots + a_1 p q^{n-1} + a_0 q^n = 0.
    $$

    $a_n p^n = -q(a_{n-1}p^{n-1} + \cdots + a_0 q^{n-1})$，故 $q \mid a_n p^n$。由 $\gcd(p, q) = 1$ 知 $q \mid a_n$。类似地 $p \mid a_0$。$\blacksquare$

!!! definition "定义 0.12 (本原多项式)"
    整系数多项式 $f(x) = a_n x^n + \cdots + a_0 \in \mathbb{Z}[x]$ 称为**本原的**（primitive），若 $\gcd(a_n, a_{n-1}, \ldots, a_0) = 1$。

!!! theorem "定理 0.21 (Gauss 引理)"
    两个本原多项式的乘积仍是本原的。

??? proof "证明"
    设 $f = \sum a_i x^i$ 和 $g = \sum b_j x^j$ 本原，$h = fg = \sum c_k x^k$。反设 $h$ 不本原，则存在素数 $p$ 使得 $p \mid c_k$ 对所有 $k$。

    由 $f$ 本原，设 $a_r$ 是第一个不被 $p$ 整除的系数（$p \mid a_0, \ldots, a_{r-1}$，$p \nmid a_r$）。类似地设 $b_s$ 是 $g$ 中第一个不被 $p$ 整除的系数。考虑 $c_{r+s}$：

    $$
    c_{r+s} = \sum_{i+j=r+s} a_i b_j = a_r b_s + \sum_{\substack{i+j=r+s \\ i \neq r}} a_i b_j.
    $$

    当 $i < r$ 时 $p \mid a_i$，当 $i > r$ 时 $j = r+s-i < s$ 故 $p \mid b_j$。因此 $p \mid (c_{r+s} - a_r b_s)$，又 $p \mid c_{r+s}$，得 $p \mid a_r b_s$。由 $p$ 为素数，$p \mid a_r$ 或 $p \mid b_s$，矛盾。$\blacksquare$

!!! theorem "定理 0.22 (Gauss 定理：$\\mathbb{Q}$ 上可约蕴含 $\\mathbb{Z}$ 上可约)"
    设 $f(x) \in \mathbb{Z}[x]$ 为本原多项式。若 $f$ 在 $\mathbb{Q}[x]$ 中可约，则 $f$ 可以写成两个次数较低的整系数多项式之积。

??? proof "证明"
    设 $f = gh$，$g, h \in \mathbb{Q}[x]$，$\deg g, \deg h \geq 1$。提取公分母后设 $\frac{a}{b} f = g_0 h_0$，其中 $g_0, h_0 \in \mathbb{Z}[x]$。设 $g_0 = d_1 g_1$，$h_0 = d_2 h_1$，$g_1, h_1$ 本原。则 $\frac{a}{b} f = d_1 d_2 g_1 h_1$。由 Gauss 引理 $g_1 h_1$ 本原，又 $f$ 本原，比较两边内容（所有系数的 GCD）得 $\frac{a}{b} = \pm d_1 d_2$，从而 $f = \pm g_1 h_1$。$\blacksquare$

!!! theorem "定理 0.23 (Eisenstein 判别法)"
    设 $f(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_0 \in \mathbb{Z}[x]$，$n \geq 1$。若存在素数 $p$ 使得：

    1. $p \nmid a_n$；
    2. $p \mid a_{n-1}, p \mid a_{n-2}, \ldots, p \mid a_0$；
    3. $p^2 \nmid a_0$。

    则 $f(x)$ 在 $\mathbb{Q}[x]$ 中不可约。

??? proof "证明"
    反设 $f$ 在 $\mathbb{Q}[x]$ 中可约。由 Gauss 定理，$f = gh$，$g, h \in \mathbb{Z}[x]$，$\deg g = r \geq 1$，$\deg h = s \geq 1$，$r + s = n$。设 $g = b_r x^r + \cdots + b_0$，$h = c_s x^s + \cdots + c_0$。

    由 $a_0 = b_0 c_0$ 和条件 (2)(3)，$p \mid b_0 c_0$ 但 $p^2 \nmid b_0 c_0$，所以 $p$ 恰整除 $b_0, c_0$ 中的一个。不妨设 $p \mid b_0$，$p \nmid c_0$。

    由 $a_n = b_r c_s$ 和条件 (1)，$p \nmid b_r$。设 $b_k$ 是 $g$ 中第一个不被 $p$ 整除的系数（$1 \leq k \leq r < n$）。考虑 $a_k = b_k c_0 + b_{k-1}c_1 + \cdots + b_0 c_k$。由 $p \mid a_k$（条件 (2)，因 $k < n$）和 $p \mid b_0, \ldots, b_{k-1}$，得 $p \mid b_k c_0$，从而 $p \mid b_k$（因 $p \nmid c_0$），矛盾。$\blacksquare$

!!! example "例 0.13"
    证明 $f(x) = x^4 + 3x^3 + 9x + 3$ 在 $\mathbb{Q}$ 上不可约。

    取 $p = 3$：$3 \nmid 1$（首项系数），$3 \mid 3, 0, 9, 3$（其余系数），$9 \nmid 3$（$p^2 = 9 \nmid a_0 = 3$）。由 Eisenstein 判别法，$f$ 在 $\mathbb{Q}$ 上不可约。

!!! example "例 0.14"
    证明 $p$ 次分圆多项式 $\Phi_p(x) = x^{p-1} + x^{p-2} + \cdots + x + 1$（$p$ 为素数）在 $\mathbb{Q}$ 上不可约。

    令 $y = x - 1$，即 $x = y + 1$。

    $$
    \Phi_p(x) = \frac{x^p - 1}{x - 1}, \quad \Phi_p(y + 1) = \frac{(y+1)^p - 1}{y} = \sum_{k=1}^{p} \binom{p}{k} y^{k-1} = y^{p-1} + \binom{p}{1} y^{p-2} + \cdots + \binom{p}{p-1}.
    $$

    由 $p \mid \binom{p}{k}$（$1 \leq k \leq p-1$）且 $p^2 \nmid \binom{p}{1} = p$，Eisenstein 判别法（取素数 $p$）给出 $\Phi_p(y+1)$ 在 $\mathbb{Q}$ 上不可约，从而 $\Phi_p(x)$ 也不可约。

---

## 0.10 多项式插值

<div class="context-flow" markdown>

**核心问题**：给定 $n + 1$ 个点，能否唯一确定一个 $n$ 次多项式通过它们？ → Lagrange 插值 → Newton 插值 → 应用

</div>

!!! theorem "定理 0.24 (多项式插值的存在唯一性)"
    设 $x_0, x_1, \ldots, x_n \in \mathbb{F}$ 两两不同，$y_0, y_1, \ldots, y_n \in \mathbb{F}$ 为给定值。则存在唯一的 $n$ 次以下多项式 $p(x) \in \mathbb{F}[x]$（$\deg p \leq n$）使得

    $$
    p(x_i) = y_i, \quad i = 0, 1, \ldots, n.
    $$

??? proof "证明"
    **唯一性**：设 $p, q$ 都满足条件，$\deg p, \deg q \leq n$。则 $p - q$ 的次数 $\leq n$ 且有 $n + 1$ 个根 $x_0, \ldots, x_n$。由定理 0.13，$p - q = 0$。

    **存在性**：由下面的 Lagrange 插值公式给出。$\blacksquare$

!!! definition "定义 0.13 (Lagrange 基多项式)"
    对给定的互不相同的节点 $x_0, x_1, \ldots, x_n \in \mathbb{F}$，**Lagrange 基多项式**（Lagrange basis polynomial）定义为

    $$
    \ell_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}, \quad i = 0, 1, \ldots, n.
    $$

    其中 $\ell_i(x_j) = \delta_{ij}$（Kronecker delta）。

!!! theorem "定理 0.25 (Lagrange 插值公式)"
    满足插值条件 $p(x_i) = y_i$（$i = 0, 1, \ldots, n$）的唯一多项式为

    $$
    p(x) = \sum_{i=0}^{n} y_i \ell_i(x) = \sum_{i=0}^{n} y_i \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}.
    $$

??? proof "证明"
    由 $\ell_i(x_k) = \delta_{ik}$，直接验证 $p(x_k) = \sum_i y_i \delta_{ik} = y_k$。$\deg p \leq n$ 也显然成立。$\blacksquare$

!!! example "例 0.15"
    求通过点 $(0, 1)$、$(1, 3)$、$(2, 7)$ 的二次插值多项式。

    节点 $x_0 = 0, x_1 = 1, x_2 = 2$，值 $y_0 = 1, y_1 = 3, y_2 = 7$。

    $$
    \ell_0(x) = \frac{(x-1)(x-2)}{(0-1)(0-2)} = \frac{(x-1)(x-2)}{2},
    $$

    $$
    \ell_1(x) = \frac{(x-0)(x-2)}{(1-0)(1-2)} = \frac{x(x-2)}{-1} = -x(x-2),
    $$

    $$
    \ell_2(x) = \frac{(x-0)(x-1)}{(2-0)(2-1)} = \frac{x(x-1)}{2}.
    $$

    $$
    p(x) = 1 \cdot \frac{(x-1)(x-2)}{2} + 3 \cdot [-x(x-2)] + 7 \cdot \frac{x(x-1)}{2}.
    $$

    展开化简：

    $$
    p(x) = \frac{x^2 - 3x + 2}{2} - 3x^2 + 6x + \frac{7x^2 - 7x}{2} = \frac{x^2 - 3x + 2 - 6x^2 + 12x + 7x^2 - 7x}{2} = \frac{2x^2 + 2x + 2}{2} = x^2 + x + 1.
    $$

    验证：$p(0) = 1$，$p(1) = 3$，$p(2) = 7$。

### Newton 插值

!!! definition "定义 0.14 (差商)"
    设 $x_0, x_1, \ldots, x_n$ 两两不同。**差商**（divided difference）递归定义为：

    - 零阶差商：$f[x_i] = f(x_i)$。
    - 一阶差商：$f[x_i, x_{i+1}] = \dfrac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}$。
    - $k$ 阶差商：$f[x_i, x_{i+1}, \ldots, x_{i+k}] = \dfrac{f[x_{i+1}, \ldots, x_{i+k}] - f[x_i, \ldots, x_{i+k-1}]}{x_{i+k} - x_i}$。

!!! theorem "定理 0.26 (Newton 插值公式)"
    插值多项式可以表示为 Newton 形式：

    $$
    p(x) = f[x_0] + f[x_0, x_1](x - x_0) + f[x_0, x_1, x_2](x - x_0)(x - x_1) + \cdots + f[x_0, \ldots, x_n]\prod_{i=0}^{n-1}(x - x_i).
    $$

??? proof "证明"
    令 $p_k(x)$ 为通过前 $k + 1$ 个点的插值多项式。则 $p_0(x) = f[x_0]$，且

    $$
    p_k(x) = p_{k-1}(x) + c_k \prod_{i=0}^{k-1}(x - x_i),
    $$

    其中 $c_k$ 由条件 $p_k(x_k) = y_k$ 确定。归纳可证 $c_k = f[x_0, x_1, \ldots, x_k]$。$\blacksquare$

!!! example "例 0.16"
    用 Newton 插值重做例 0.15（节点 $(0,1), (1,3), (2,7)$）。

    差商表：

    | $x_i$ | $f[x_i]$ | $f[x_i, x_{i+1}]$ | $f[x_i, x_{i+1}, x_{i+2}]$ |
    |---|---|---|---|
    | $0$ | $1$ | | |
    | | | $\frac{3-1}{1-0} = 2$ | |
    | $1$ | $3$ | | $\frac{4-2}{2-0} = 1$ |
    | | | $\frac{7-3}{2-1} = 4$ | |
    | $2$ | $7$ | | |

    $$
    p(x) = 1 + 2(x - 0) + 1 \cdot (x - 0)(x - 1) = 1 + 2x + x^2 - x = x^2 + x + 1.
    $$

    与 Lagrange 方法结果一致。

!!! note "注"
    Newton 插值的优势在于**递推性**：当增加一个新节点时，只需在原多项式基础上增加一项，无需重新计算全部。此外差商具有对称性：$f[x_0, x_1, \ldots, x_k]$ 与节点的排列顺序无关。
