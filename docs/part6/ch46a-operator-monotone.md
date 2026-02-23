# 第 46A 章 算子单调函数与算子凸函数

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 矩阵函数(Ch13) · 矩阵不等式(Ch18) · Löwner 偏序

**本章脉络**：Löwner 偏序回顾 → 算子单调函数定义 → 经典例子（$t^r$, $\log t$, $t/(1+t)$）→ Löwner-Heinz 不等式（完整证明）→ Furuta 不等式 → Löwner 定理（积分表示与 Pick 函数）→ Löwner 矩阵（差商矩阵）→ $n$-单调函数 → 算子凸/凹函数 → Hansen-Pedersen 刻画 → Jensen 算子不等式 → 算子凸与单调的关系 → Choi 定理 → Lieb 凹性定理 → Lieb-Ruskai 定理

**延伸**：算子单调函数是量子信息中量子相对熵和保真度的数学基础；Lieb 凹性定理在量子统计力学中证明了量子熵的强次可加性——这是量子信息论最深刻的结果之一

</div>

当我们将标量函数推广到矩阵函数时，许多在标量世界中显而易见的性质可能不再成立。最典型的例子是：对正实数，$a \geq b > 0$ 意味着 $a^2 \geq b^2$，但对正定矩阵，$A \succeq B \succ 0$ **不**意味着 $A^2 \succeq B^2$。这个看似简单的观察引出了一个深刻的问题：**哪些函数保持矩阵偏序？**

Löwner 在 1934 年的开创性工作中完全回答了这个问题：一个定义在正实数上的连续函数保持所有大小矩阵的偏序，当且仅当它具有特定的积分表示——即它是一个 Pick 函数（在上半平面取正虚部值的解析函数）。这一定理将矩阵分析、复分析和测度论深刻地联系在一起。

本章从算子单调函数出发，经过 Löwner-Heinz 不等式和 Furuta 不等式，到 Löwner 定理的完整陈述，再到算子凸/凹函数理论及其在量子信息中的应用（Lieb 凹性定理）。

---

## 46A.1 Löwner 偏序回顾

<div class="context-flow" markdown>

**核心问题**：正定矩阵上的偏序关系具有哪些基本性质？它与标量序有何本质区别？

</div>

!!! definition "定义 46A.1 (Löwner 偏序)"
    设 $A, B \in \mathbb{C}^{n \times n}$ 为 Hermite 矩阵。定义 **Löwner 偏序**：

    $$A \succeq B \quad \Leftrightarrow \quad A - B \text{ 是半正定的} \quad \Leftrightarrow \quad x^*(A - B)x \geq 0 \;\; \forall x.$$

    严格版本：$A \succ B \Leftrightarrow A - B$ 正定。

    特别地，$A \succeq 0$ 表示 $A$ 半正定，$A \succ 0$ 表示 $A$ 正定。

!!! theorem "定理 46A.1 (Löwner 偏序的基本性质)"
    Löwner 偏序是 Hermite 矩阵集合上的偏序关系，但**不是**全序。具体性质：

    1. **自反性**：$A \succeq A$。
    2. **反对称性**：$A \succeq B$ 且 $B \succeq A$ $\Rightarrow$ $A = B$。
    3. **传递性**：$A \succeq B$ 且 $B \succeq C$ $\Rightarrow$ $A \succeq C$。
    4. **不是全序**：存在 $A, B$ 使得 $A \not\succeq B$ 且 $B \not\succeq A$。
    5. **合同不变**：$A \succeq B$ $\Rightarrow$ $C^*AC \succeq C^*BC$ 对任意 $C$。
    6. **加法保持**：$A \succeq B$ 且 $C \succeq D$ $\Rightarrow$ $A + C \succeq B + D$。
    7. **乘法不保持**：$A \succeq B \succeq 0$ **不**意味着 $A^2 \succeq B^2$。

??? proof "证明"
    (1)-(6) 是半正定锥的标准性质，直接由定义验证。

    **(7)** 反例：$A = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $A - B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} \succeq 0$（特征值 $0, 2$），因此 $A \succeq B \succeq 0$。

    $A^2 = \begin{pmatrix} 5 & 3 \\ 3 & 2 \end{pmatrix}$，$B^2 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $A^2 - B^2 = \begin{pmatrix} 4 & 3 \\ 3 & 2 \end{pmatrix}$，$\det(A^2 - B^2) = 8 - 9 = -1 < 0$。

    因此 $A^2 - B^2$ 不是半正定的，即 $A^2 \not\succeq B^2$。 $\blacksquare$

!!! example "例 46A.1 (Löwner 偏序的反直觉)"
    以下标量性质在矩阵中的命运：

    | 标量性质 | 矩阵版本 | 成立？ |
    |---------|---------|--------|
    | $a \geq b > 0 \Rightarrow a^2 \geq b^2$ | $A \succeq B \succ 0 \Rightarrow A^2 \succeq B^2$ | 否 |
    | $a \geq b > 0 \Rightarrow 1/b \geq 1/a$ | $A \succeq B \succ 0 \Rightarrow B^{-1} \succeq A^{-1}$ | **是** |
    | $a \geq b > 0 \Rightarrow \sqrt{a} \geq \sqrt{b}$ | $A \succeq B \succ 0 \Rightarrow A^{1/2} \succeq B^{1/2}$ | **是** |
    | $a \geq b > 0 \Rightarrow e^a \geq e^b$ | $A \succeq B \Rightarrow e^A \succeq e^B$ | 否 |
    | $a \geq b > 0 \Rightarrow \log a \geq \log b$ | $A \succeq B \succ 0 \Rightarrow \log A \succeq \log B$ | **是** |

    保持偏序的函数（如 $t^{1/2}$, $t^{-1}$, $\log t$）和不保持的（如 $t^2$, $e^t$）之间的精确分界，正是算子单调性理论要回答的问题。

---

## 46A.2 算子单调函数

<div class="context-flow" markdown>

**核心问题**：哪些函数将矩阵偏序保持为矩阵偏序？

</div>

!!! definition "定义 46A.2 (算子单调函数)"
    设 $I \subseteq \mathbb{R}$ 是区间。连续函数 $f: I \to \mathbb{R}$ 称为**算子单调函数**（operator monotone function），如果对任意正整数 $n$ 和谱在 $I$ 中的 Hermite 矩阵 $A, B \in \mathbb{C}^{n \times n}$：

    $$A \succeq B \quad \Rightarrow \quad f(A) \succeq f(B).$$

    若只对固定的 $n$ 要求，则称 $f$ 是 **$n$-单调**的（$n$-monotone）。

    显然算子单调 $\Rightarrow$ 对每个 $n$ 都是 $n$-单调 $\Rightarrow$ $1$-单调（即标量单调递增）。

!!! theorem "定理 46A.2 (算子单调函数的例子)"
    1. **$f(t) = t^r$（$0 < r \leq 1$）** 在 $[0, \infty)$ 上是算子单调的（**Löwner-Heinz 不等式**）。
    2. **$f(t) = \log t$** 在 $(0, \infty)$ 上是算子单调的。
    3. **$f(t) = \frac{t}{1 + t}$** 在 $[0, \infty)$ 上是算子单调的。
    4. **$f(t) = \frac{t - 1}{t + 1}$** 在 $(0, \infty)$ 上是算子单调的。
    5. **$f(t) = t^r$（$r > 1$）** 在 $[0, \infty)$ 上**不是**算子单调的。
    6. **$f(t) = e^t$** 在 $\mathbb{R}$ 上**不是**算子单调的。

---

## 46A.3 Löwner-Heinz 不等式

<div class="context-flow" markdown>

**核心问题**：如何证明 $A \succeq B \succeq 0 \Rightarrow A^r \succeq B^r$（$0 < r \leq 1$）？

</div>

!!! theorem "定理 46A.3 (Löwner-Heinz 不等式)"
    设 $A, B$ 是 $n \times n$ 正半定 Hermite 矩阵。若 $A \succeq B \succeq 0$，则对 $0 < r \leq 1$：
    $$A^r \succeq B^r.$$

    等价地，$f(t) = t^r$（$0 < r \leq 1$）是 $[0, \infty)$ 上的算子单调函数。

    特别地，$r = 1/2$：$A \succeq B \succeq 0 \Rightarrow A^{1/2} \succeq B^{1/2}$。

??? proof "证明"
    **完整证明**分为两个步骤。

    **步骤 1**：积分表示。对 $0 < r < 1$，$t > 0$，有
    $$t^r = \frac{\sin(r\pi)}{\pi} \int_0^\infty \frac{t}{t + s} s^{r-1} \, ds.$$

    这是 Beta 函数的一个经典结果。验证：令 $u = s/t$，
    $$\frac{\sin(r\pi)}{\pi} \int_0^\infty \frac{t}{t + s} s^{r-1} ds = \frac{\sin(r\pi)}{\pi} t^r \int_0^\infty \frac{u^{r-1}}{1 + u} du = t^r \cdot \frac{\sin(r\pi)}{\pi} \cdot \frac{\pi}{\sin(r\pi)} = t^r.$$

    **步骤 2**：矩阵版本的积分。对正定矩阵 $A \succ 0$，由谱映射定理（spectral mapping theorem），
    $$A^r = \frac{\sin(r\pi)}{\pi} \int_0^\infty A(A + sI)^{-1} s^{r-1} \, ds.$$

    因此，只需证明每个"原子"$g_s(A) = A(A + sI)^{-1}$ 是算子单调的。

    **关键引理**：若 $A \succeq B \succeq 0$，则对每个 $s > 0$，
    $$A(A + sI)^{-1} \succeq B(B + sI)^{-1}.$$

    **引理的证明**：注意 $A(A + sI)^{-1} = I - s(A + sI)^{-1}$。因此
    $$A(A + sI)^{-1} - B(B + sI)^{-1} = s[(B + sI)^{-1} - (A + sI)^{-1}].$$

    由 $A \succeq B$，$A + sI \succeq B + sI \succ 0$。逆函数 $t \mapsto t^{-1}$ 是算子反单调的（$T \succeq S \succ 0 \Rightarrow S^{-1} \succeq T^{-1}$，这可以直接验证：$S^{-1} - T^{-1} = S^{-1}(T - S)T^{-1} \succeq 0$），因此 $(B + sI)^{-1} \succeq (A + sI)^{-1}$。

    从而 $A(A + sI)^{-1} - B(B + sI)^{-1} = s[(B + sI)^{-1} - (A + sI)^{-1}] \succeq 0$。

    **步骤 3**：逐积分。
    $$A^r - B^r = \frac{\sin(r\pi)}{\pi} \int_0^\infty [A(A + sI)^{-1} - B(B + sI)^{-1}] s^{r-1} ds \succeq 0,$$
    因为被积函数在每个 $s > 0$ 处都是半正定的。

    **$B$ 半正定但可能奇异的情形**：取 $B_\epsilon = B + \epsilon I \succ 0$，由上述证明 $A^r \succeq B_\epsilon^r$，令 $\epsilon \to 0$，由连续性得 $A^r \succeq B^r$。

    **$r = 1$ 的情形**：$A^1 = A \succeq B = B^1$，这是假设本身，平凡成立。 $\blacksquare$

!!! example "例 46A.2 (Löwner-Heinz 不等式的指数限制)"
    $r > 1$ 时 Löwner-Heinz 不等式**不成立**。取 $A = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $A \succeq B \succeq 0$，但 $A^2 \not\succeq B^2$（见定理 46A.1 的证明）。

    这表明 $t^r$（$r > 1$）不是算子单调的，Löwner-Heinz 不等式中 $r \leq 1$ 的条件是**最优的**。

---

## 46A.4 Furuta 不等式

<div class="context-flow" markdown>

**核心问题**：Löwner-Heinz 不等式能否进一步推广？

</div>

!!! theorem "定理 46A.4 (Furuta 不等式)"
    设 $A \succeq B \succeq 0$。则对 $r \geq 0$，$p \geq 0$，$q \geq 1$，且 $(1 + 2r)q \geq p + 2r$，有

    $$(A^r B^p A^r)^{1/q} \preceq A^{(p+2r)/q},$$
    $$(B^r A^p B^r)^{1/q} \succeq B^{(p+2r)/q}.$$

    特别地：
    - 取 $r = 0$：$B^{p/q} \preceq A^{p/q}$，即 Löwner-Heinz 不等式（当 $p/q \leq 1$ 时）。
    - 取 $p = 1$：$(A^r B A^r)^{1/q} \preceq A^{(1+2r)/q}$。

    Furuta 不等式是 Löwner-Heinz 不等式的重大推广，覆盖了一个二维参数区域。

??? proof "证明"
    **证明概要**（Furuta, 1987）。

    **关键情形**：先证 $q = (p + 2r)/(1 + 2r)$ 的边界情形，然后由 Löwner-Heinz 不等式推导一般情形。

    **步骤 1**（$p = 1$ 的情形）：设 $A \succeq B \succeq 0$。需证
    $$(A^r B A^r)^{(1+2r)^{-1}} \preceq A.$$

    设 $T = A^{-1/2} B A^{-1/2}$，则 $0 \preceq T \preceq I$（因为 $B \preceq A$）。需证
    $$A^{1/2}(A^{r-1/2} B A^{r-1/2})^{(1+2r)^{-1}} A^{1/2} \preceq A,$$
    即 $(A^{r-1/2} B A^{r-1/2})^{(1+2r)^{-1}} \preceq I$。

    令 $S = A^{(2r-1)/2} B A^{(2r-1)/2}$。由 $B \preceq A$ 和 Löwner-Heinz 不等式的已知结果，可以逐步建立 $S^{1/(1+2r)} \preceq I$。

    **步骤 2**（一般 $p$）：由步骤 1 和 Löwner-Heinz 不等式迭代。具体地，对 $A \succeq B \succeq 0$，Löwner-Heinz 给出 $A^s \succeq B^s$（$0 \leq s \leq 1$），将步骤 1 应用于 $A^s$ 和 $B^s$，取适当的 $s$ 得到一般结果。

    **步骤 3**（从边界情形到一般 $q$）：若 $(1+2r)q \geq p + 2r$，则 $(p+2r)/q \leq 1 + 2r$。由边界情形 $(A^r B^p A^r)^{1/q_0} \preceq A^{(p+2r)/q_0}$（$q_0 = (p+2r)/(1+2r)$），对两侧应用 Löwner-Heinz 不等式（指数 $q_0/q \leq 1$）即得。 $\blacksquare$

!!! example "例 46A.3 (Furuta 不等式的参数区域)"
    Furuta 不等式的参数 $(p, q, r)$ 满足 $r \geq 0$, $p \geq 0$, $q \geq 1$, $(1+2r)q \geq p + 2r$。

    在 $(p, q)$ 平面上（固定 $r$），可行区域是 $q \geq 1$ 且 $q \geq (p + 2r)/(1 + 2r)$ 的区域。边界线 $q = (p + 2r)/(1 + 2r)$ 通过点 $(1, 1)$，斜率为 $1/(1+2r)$。

    - $r = 0$：边界线 $q = p$，即 Löwner-Heinz（$p/q \leq 1$）。
    - $r \to \infty$：边界线趋于 $q = 1$（水平线），即对足够大的 $r$，几乎任何 $p$ 都可取 $q = 1$。

---

## 46A.5 Löwner 定理

<div class="context-flow" markdown>

**核心问题**：算子单调函数的完整刻画是什么？

</div>

!!! definition "定义 46A.3 (Pick 函数)"
    函数 $f$ 称为 **Pick 函数**（也叫 Nevanlinna 函数、Herglotz 函数），如果：

    1. $f$ 在上半平面 $\mathbb{C}^+ = \{z : \operatorname{Im} z > 0\}$ 上解析。
    2. $\operatorname{Im} f(z) \geq 0$ 当 $\operatorname{Im} z > 0$（即 $f$ 将上半平面映入上半平面的闭包）。

!!! theorem "定理 46A.5 (Löwner 定理)"
    设 $f: (0, \infty) \to \mathbb{R}$ 是连续函数。则以下三个条件等价：

    **(i)** $f$ 是**算子单调**的：$A \succeq B \succ 0 \Rightarrow f(A) \succeq f(B)$（对所有阶数 $n$）。

    **(ii)** $f$ 具有**积分表示**：
    $$f(t) = a + bt + \int_0^\infty \left(\frac{t}{t + s} - \frac{s}{1 + s}\right) d\mu(s),$$
    其中 $a \in \mathbb{R}$，$b \geq 0$，$\mu$ 是 $(0, \infty)$ 上的正 Borel 测度，满足 $\int_0^\infty \frac{d\mu(s)}{1 + s} < \infty$。

    **(iii)** $f$ 可以解析延拓到上半平面 $\mathbb{C}^+$，成为 **Pick 函数**：$\operatorname{Im} f(z) \geq 0$ 当 $\operatorname{Im} z > 0$。

??? proof "证明"
    **(ii) $\Rightarrow$ (i)**（积分表示 $\Rightarrow$ 算子单调）：

    由积分表示，$f(t)$ 是函数 $t \mapsto a + bt$ 和 $t \mapsto \frac{t}{t+s}$ 的积分叠加。

    - $g_0(t) = a + bt$ 显然是算子单调的（$A \succeq B \Rightarrow a I + b A \succeq a I + b B$）。
    - $g_s(t) = \frac{t}{t+s}$（$s > 0$）的算子单调性：$g_s(A) = A(A + sI)^{-1}$，其算子单调性已在 Löwner-Heinz 不等式的证明中验证（定理 46A.3 步骤 2 的关键引理）。

    由算子单调函数在积分下的稳定性（非负测度的积分保持半正定偏序），$f$ 算子单调。

    **(i) $\Rightarrow$ (iii)**（算子单调 $\Rightarrow$ Pick 函数）：

    $f$ 算子单调意味着 $f$ 对每个 $n$ 都是 $n$-单调的。由 Löwner 矩阵的半正定性（定理 46A.6），$f$ 的差商核在每个有限集上都是正定的。这意味着 $f$ 必须是实解析的。

    对 $z = x + iy$（$y > 0$），考虑 $2 \times 2$ 矩阵 $Z = \begin{pmatrix} x & y \\ y & x \end{pmatrix}$（Hermite 矩阵，特征值 $x + y, x - y$）。$f$ 的 $2$-单调性应用于 $Z$ 的扰动，可以推导出 $\operatorname{Im} f(z) \geq 0$。

    更严格的论证：$f$ 的 $n$-单调性对所有 $n$ 意味着 Löwner 矩阵对所有 $n$ 半正定（定理 46A.6），这是一个关于 $f$ 的差商的正定核条件。由正定核的 Herglotz-Nevanlinna 表示定理，$f$ 延拓为 Pick 函数。

    **(iii) $\Rightarrow$ (ii)**（Pick 函数 $\Rightarrow$ 积分表示）：

    这是经典的 **Nevanlinna 表示定理**。Pick 函数 $f: \mathbb{C}^+ \to \overline{\mathbb{C}^+}$ 有唯一的表示
    $$f(z) = a + bz + \int_0^\infty \left(\frac{1}{s - z} - \frac{s}{1 + s^2}\right) d\mu(s),$$
    其中 $a \in \mathbb{R}$，$b \geq 0$，$\mu$ 是正测度。限制到正实轴 $z = t > 0$，经变量替换得到所述的积分表示。 $\blacksquare$

!!! example "例 46A.4 (验证 Pick 函数)"
    1. **$f(t) = t^r$（$0 < r < 1$）**：$f(z) = z^r = e^{r \log z}$（取上半平面的主值分支）。$\operatorname{Im}(z^r) = |z|^r \sin(r \arg z)$。由 $0 < \arg z < \pi$ 且 $0 < r < 1$，有 $0 < r \arg z < \pi$，因此 $\sin(r \arg z) > 0$。$f$ 是 Pick 函数。

    2. **$f(t) = \log t$**：$f(z) = \log z = \log|z| + i \arg z$。$\operatorname{Im}(\log z) = \arg z \in (0, \pi)$ 当 $\operatorname{Im} z > 0$。$f$ 是 Pick 函数。

    3. **$f(t) = t^2$**：取 $z = -1 + i$，$z^2 = 1 - 1 - 2i = -2i$，$\operatorname{Im}(z^2) = -2 < 0$。因此 $t^2$ 不是 Pick 函数，与其不是算子单调的一致。

    4. **$f(t) = \frac{t}{1+t}$**：$f(z) = \frac{z}{1+z}$。$\operatorname{Im}\left(\frac{z}{1+z}\right) = \frac{\operatorname{Im} z}{|1+z|^2} > 0$ 当 $\operatorname{Im} z > 0$。$f$ 是 Pick 函数。

---

## 46A.6 Löwner 矩阵

<div class="context-flow" markdown>

**核心问题**：如何用有限维矩阵条件刻画 $n$-单调性？

</div>

!!! definition "定义 46A.4 (Löwner 矩阵)"
    设 $f: I \to \mathbb{R}$ 是定义在区间 $I$ 上的函数，$t_1, t_2, \ldots, t_n \in I$ 是 $n$ 个不同的点。$f$ 关于 $t_1, \ldots, t_n$ 的 **Löwner 矩阵**（又称差商矩阵，divided difference matrix）定义为 $n \times n$ 矩阵

    $$L_n(f; t_1, \ldots, t_n) = \left[\frac{f(t_i) - f(t_j)}{t_i - t_j}\right]_{i,j=1}^{n},$$

    其中对角元素定义为 $f'(t_i)$（当 $f$ 可微时），或更一般地通过极限 $\lim_{t_j \to t_i} \frac{f(t_i) - f(t_j)}{t_i - t_j}$。

    $L_n$ 的第 $(i,j)$ 元素 $[f(t_i) - f(t_j)]/(t_i - t_j)$ 称为 $f$ 在 $t_i, t_j$ 处的**一阶差商**。

!!! theorem "定理 46A.6 (Löwner 判据)"
    设 $f: I \to \mathbb{R}$ 连续。则以下等价：

    1. $f$ 是 **$n$-单调**的。
    2. 对 $I$ 中的任意 $n$ 个不同点 $t_1 < t_2 < \cdots < t_n$，Löwner 矩阵 $L_n(f; t_1, \ldots, t_n) \succeq 0$（半正定）。

    因此，$f$ 是**算子单调**的当且仅当 $L_n \succeq 0$ 对**所有** $n$ 和所有不同点成立。

??? proof "证明"
    **必要性**（$n$-单调 $\Rightarrow$ Löwner 矩阵半正定）：

    设 $f$ 是 $n$-单调的，$t_1 < \cdots < t_n$ 是 $I$ 中的不同点。构造对角矩阵 $D = \operatorname{diag}(t_1, \ldots, t_n)$ 和对任意 $\epsilon > 0$，矩阵 $D + \epsilon E$（$E$ 是适当的扰动）。

    取 $A = \operatorname{diag}(t_1, \ldots, t_n)$，$B = A - \epsilon vv^*$（$v \in \mathbb{C}^n$），使得 $A \succeq B$。

    由 $n$-单调性，$f(A) \succeq f(B)$，即 $f(A) - f(B) \succeq 0$。

    $f(A) = \operatorname{diag}(f(t_1), \ldots, f(t_n))$。$f(B)$ 的计算利用一阶扰动展开：
    $$f(B) \approx f(A) - \epsilon \, L_n \circ (vv^*) + O(\epsilon^2),$$
    其中 $L_n \circ (vv^*)$ 表示 Löwner 矩阵与 $vv^*$ 的 Hadamard 乘积。

    更精确地，由矩阵函数的一阶 Frechet 导数公式，
    $$f(A) - f(A - \epsilon vv^*) = \epsilon \, L_n \circ (vv^*) + O(\epsilon^2).$$

    由 $f(A) - f(B) \succeq 0$，得 $L_n \circ (vv^*) \succeq 0$ 对所有 $v$。由 Schur 乘积定理（半正定矩阵的 Hadamard 乘积仍半正定），这等价于 $L_n \succeq 0$。

    **充分性**的证明更为技术性，需要利用 Löwner 矩阵的半正定性来重构 $f$ 的解析延拓。 $\blacksquare$

!!! example "例 46A.5 (Löwner 矩阵的计算)"
    对 $f(t) = \sqrt{t}$ 和点 $t_1 = 1, t_2 = 4, t_3 = 9$：

    $$L_3 = \begin{pmatrix}
    \frac{1}{2\sqrt{1}} & \frac{\sqrt{1} - \sqrt{4}}{1 - 4} & \frac{\sqrt{1} - \sqrt{9}}{1 - 9} \\
    \frac{\sqrt{4} - \sqrt{1}}{4 - 1} & \frac{1}{2\sqrt{4}} & \frac{\sqrt{4} - \sqrt{9}}{4 - 9} \\
    \frac{\sqrt{9} - \sqrt{1}}{9 - 1} & \frac{\sqrt{9} - \sqrt{4}}{9 - 4} & \frac{1}{2\sqrt{9}}
    \end{pmatrix}
    = \begin{pmatrix}
    1/2 & 1/3 & 1/4 \\
    1/3 & 1/4 & 1/5 \\
    1/4 & 1/5 & 1/6
    \end{pmatrix}.$$

    这是一个 Hilbert 型矩阵，可以验证其正定性（所有特征值为正），与 $\sqrt{t}$ 是算子单调的一致。

---

## 46A.7 $n$-单调函数

<div class="context-flow" markdown>

**核心问题**：$n$-单调性是否随 $n$ 严格变强？

</div>

!!! theorem "定理 46A.7 ($n$-单调 $\neq$ $(n+1)$-单调)"
    对每个 $n \geq 1$，存在函数 $f$ 是 $n$-单调但不是 $(n+1)$-单调的。即 $n$-单调函数类随 $n$ **严格递减**：

    $$\mathcal{P}_1 \supsetneq \mathcal{P}_2 \supsetneq \mathcal{P}_3 \supsetneq \cdots \supsetneq \mathcal{P}_\infty,$$

    其中 $\mathcal{P}_n$ 是 $n$-单调函数的集合，$\mathcal{P}_\infty$ 是算子单调函数的集合。

!!! example "例 46A.6 ($n$-单调但非 $(n+1)$-单调的例子)"
    1. **$f(t) = t^2$**：$f$ 是 1-单调的（标量单调递增，$t > 0$），但不是 2-单调的（反例见定理 46A.1 证明 (7)）。

    2. **$f(t) = \frac{t^{n}}{(1+t)^{n-1}}$** 在 $[0, \infty)$ 上是 $n$-单调但不是 $(n+1)$-单调的。

    3. 对 $[0, 1]$ 上的函数，Dobsch (1937) 构造了显式的 $n$-单调非 $(n+1)$-单调的例子：利用 Löwner 矩阵恰好在某个 $(n+1)$ 点配置上退化的函数。

    4. **有理函数的例子**：$f(t) = \frac{t}{1 + t/n}$ 对固定 $n$ 是 $n$-单调的（可以验证 Löwner 矩阵对 $n$ 个点半正定），但 $f(t) = t/(1 + t/n) \to t$（$n \to \infty$），极限 $t$ 是算子单调的。然而对有限 $n$，每个 $f_n$ 只是 $n$-单调的。

---

## 46A.8 算子凸与算子凹函数

<div class="context-flow" markdown>

**核心问题**：除了保持偏序，哪些函数保持矩阵的"凸性结构"？

</div>

!!! definition "定义 46A.5 (算子凸与算子凹函数)"
    连续函数 $f: I \to \mathbb{R}$ 称为**算子凸**（operator convex），如果对谱在 $I$ 中的 Hermite 矩阵 $A, B$ 和 $\lambda \in [0, 1]$：

    $$f(\lambda A + (1 - \lambda)B) \preceq \lambda f(A) + (1 - \lambda) f(B).$$

    $f$ 称为**算子凹**（operator concave），如果 $-f$ 算子凸，即

    $$f(\lambda A + (1 - \lambda)B) \succeq \lambda f(A) + (1 - \lambda) f(B).$$

!!! theorem "定理 46A.8 (算子凸/凹函数的例子)"
    1. **$f(t) = t^r$（$1 \leq r \leq 2$）** 在 $[0, \infty)$ 上是**算子凸**的。
    2. **$f(t) = t^r$（$0 \leq r \leq 1$）** 在 $[0, \infty)$ 上是**算子凹**的。
    3. **$f(t) = \log t$** 在 $(0, \infty)$ 上是**算子凹**的。
    4. **$f(t) = t^{-1}$** 在 $(0, \infty)$ 上是**算子凸**的。
    5. **$f(t) = t \log t$** 在 $(0, \infty)$ 上是**算子凸**的。
    6. **$f(t) = e^t$** 在 $\mathbb{R}$ 上**既不是**算子凸**也不是**算子凹的。

---

## 46A.9 Hansen-Pedersen 刻画

!!! theorem "定理 46A.9 (Hansen-Pedersen 定理)"
    连续函数 $f: [0, \infty) \to \mathbb{R}$ 算子凸当且仅当对所有 $n$ 和所有 $A \in \mathbb{C}^{n \times n}$（$A$ Hermite，谱在 $[0, \infty)$ 中），以下 $2n \times 2n$ **分块矩阵条件**成立：

    $$\begin{pmatrix} f(A) & f(A)^{1/2} \\ f(A)^{1/2} & f(A) \end{pmatrix} \text{ 的适当推广} \succeq 0.$$

    更精确地，$f$ 算子凸当且仅当对任意 $2 \times 2$ 自伴分块矩阵 $\begin{pmatrix} A & C \\ C^* & B \end{pmatrix} \succeq 0$（$A, B$ 谱在定义域中），有

    $$f\left(\begin{pmatrix} A & C \\ C^* & B \end{pmatrix}\right) \succeq \begin{pmatrix} f(A) & 0 \\ 0 & f(B) \end{pmatrix}.$$

    这意味着算子凸性可以通过 **$2 \times 2$ 分块测试**来验证——只需在 $2 \times 2$ 分块矩阵上验证凸性不等式。

??? proof "证明"
    **必要性**（算子凸 $\Rightarrow$ 分块条件）：

    设 $f$ 算子凸，$M = \begin{pmatrix} A & C \\ C^* & B \end{pmatrix} \succeq 0$。取等距算子 $V_1 = \begin{pmatrix} I \\ 0 \end{pmatrix}$ 和 $V_2 = \begin{pmatrix} 0 \\ I \end{pmatrix}$。

    由 $M \succeq 0$ 且 $f$ 算子凸，Jensen 不等式（定理 46A.10）给出
    $$V_i^* f(M) V_i \succeq f(V_i^* M V_i).$$

    $V_1^* M V_1 = A$，$V_2^* M V_2 = B$。因此 $f(M)$ 的对角块满足
    $$[f(M)]_{11} \succeq f(A), \quad [f(M)]_{22} \succeq f(B).$$

    这给出了分块条件的一个方向。完整的等式需要更精细的论证。

    **充分性**的证明利用了 $2 \times 2$ 分块条件递归推导一般凸性不等式。 $\blacksquare$

---

## 46A.10 Jensen 算子不等式

!!! theorem "定理 46A.10 (Jensen 算子不等式)"
    设 $f$ 是算子凸函数，$C_1, \ldots, C_k$ 是满足 $\sum_{i=1}^{k} C_i^* C_i = I$ 的算子。则

    $$f\left(\sum_{i=1}^{k} C_i^* A_i C_i\right) \preceq \sum_{i=1}^{k} C_i^* f(A_i) C_i.$$

    特别地，取 $k = 1$，$C_1 = V$（等距算子，$V^*V = I$）：
    $$f(V^*AV) \preceq V^*f(A)V.$$

    若 $f$ 算子凹，不等式方向反转。

??? proof "证明"
    **$k = 2$ 的情形**：取 $C_1 = \sqrt{\lambda} I$，$C_2 = \sqrt{1-\lambda} I$，$A_1 = A$，$A_2 = B$：
    $$f(\lambda A + (1-\lambda)B) \preceq \lambda f(A) + (1-\lambda)f(B).$$
    这就是算子凸的定义。

    **一般 $k$**：由归纳法。设结论对 $k - 1$ 成立。令 $\lambda = \|C_k\|^2$（假设 $0 < \lambda < 1$），定义 $D_i = C_i / \sqrt{1 - \lambda}$（$i < k$），则 $\sum_{i=1}^{k-1} D_i^* D_i = I$。
    $$\sum_{i=1}^k C_i^* A_i C_i = (1-\lambda) \sum_{i=1}^{k-1} D_i^* A_i D_i + C_k^* A_k C_k.$$

    由算子凸性和归纳假设完成证明。

    **等距算子情形**的另一种直接证明：$V^*V = I$，$P = VV^*$ 是投影。
    $$V^*f(A)V = V^* f(A) V \succeq f(V^*AV),$$
    这可以通过将 $A$ 写成 $V(V^*AV)V^* + (I - P)A(I-P) + $ 交叉项，利用 $f$ 的算子凸性和投影的性质推导。 $\blacksquare$

---

## 46A.11 算子凸与算子单调的关系

!!! theorem "定理 46A.11 (算子凸与算子单调的关系)"
    1. 若 $f$ 在 $(0, \infty)$ 上算子凸且 $f(0) = \lim_{t \to 0^+} f(t) \leq 0$，则 $f$ 是算子单调**递增**的。

    2. 若 $f$ 在 $(0, \infty)$ 上算子单调递增且 $f(0) \leq 0$，则 $f$ 是**算子凹**的。

    3. $f$ 在 $(0, \infty)$ 上算子凸当且仅当 $f'$ 存在且是算子单调**递增**的。

    综合来说，在 $(0, \infty)$ 上有链：
    $$f \text{ 算子凸} \quad \Rightarrow \quad f' \text{ 算子单调} \quad \Rightarrow \quad f' \text{ 算子凹（当 } f'(0) \leq 0\text{）}.$$

!!! example "例 46A.7 (关系的例子)"
    - $f(t) = t^2$：算子凸（$1 \leq 2 \leq 2$），$f'(t) = 2t$ 是算子单调的。
    - $f(t) = t \log t$：算子凸，$f'(t) = 1 + \log t$ 是算子单调的（$\log t$ 算子单调）。
    - $f(t) = -\log t$：算子凸（因为 $\log t$ 算子凹），$f'(t) = -1/t$ 是算子单调递减的。
    - $f(t) = t^{1/2}$：算子单调且 $f(0) = 0$，因此算子凹。

---

## 46A.12 Choi 定理

!!! theorem "定理 46A.12 (Choi 定理)"
    连续函数 $f: [0, \infty) \to [0, \infty)$ 算子凸当且仅当对所有 $n$，所有正半定矩阵 $A \in \mathbb{C}^{n \times n}$，和所有**完全正映射** $\Phi: \mathbb{C}^{n \times n} \to \mathbb{C}^{m \times m}$（$\Phi$ 保持正性且其放大 $\Phi \otimes \operatorname{id}_k$ 也保持正性），有

    $$f(\Phi(A)) \preceq \Phi(f(A)).$$

    反之，若 $f$ 不是算子凸的，则存在完全正映射 $\Phi$ 使得上述不等式不成立。

    Choi 定理将算子凸性与量子信息中的**量子信道**（完全正映射）联系起来：算子凸函数恰好是与所有量子信道"相容"的函数。

---

## 46A.13 Lieb 凹性定理

<div class="context-flow" markdown>

**核心问题**：矩阵幂函数在两个变量上的联合凹性如何？

</div>

!!! theorem "定理 46A.13 (Lieb 凹性定理)"
    设 $K \in \mathbb{C}^{n \times m}$ 是固定矩阵，$0 < p < 1$。则映射

    $$(A, B) \mapsto \operatorname{tr}(K^* A^p K B^{1-p})$$

    在 $(A, B) \in \mathcal{P}_n \times \mathcal{P}_m$（正定矩阵对）上是**联合凹**的。即对 $\lambda \in [0, 1]$：

    $$\operatorname{tr}(K^* (\lambda A_1 + (1-\lambda)A_2)^p K (\lambda B_1 + (1-\lambda)B_2)^{1-p})$$
    $$\geq \lambda \operatorname{tr}(K^* A_1^p K B_1^{1-p}) + (1-\lambda)\operatorname{tr}(K^* A_2^p K B_2^{1-p}).$$

??? proof "证明"
    **证明概要**（Lieb, 1973；简化版本 Ando, 1979）。

    **步骤 1**（利用 Epstein 方法）：定义函数
    $$F(A, B) = \operatorname{tr}(K^* A^p K B^{1-p}).$$

    **步骤 2**（利用积分表示）：由 $A^p = \frac{\sin(p\pi)}{\pi} \int_0^\infty A(A + tI)^{-1} t^{p-1} dt$（Löwner-Heinz 证明中的积分表示），

    $$F(A, B) = \frac{\sin(p\pi)}{\pi} \int_0^\infty \operatorname{tr}(K^* A(A+tI)^{-1} K B^{1-p}) t^{p-1} dt.$$

    因此只需证明对每个 $t > 0$，
    $$G_t(A, B) = \operatorname{tr}(K^* A(A+tI)^{-1} K B^{1-p})$$
    关于 $(A, B)$ 联合凹。

    **步骤 3**（简化）：$A(A+tI)^{-1} = I - t(A + tI)^{-1}$ 关于 $A$ 是算子凹的（因为 $A \mapsto (A+tI)^{-1}$ 算子凸）。$B^{1-p}$ 关于 $B$ 算子凹（$0 < 1-p < 1$，Löwner-Heinz）。

    **步骤 4**（关键技术）：利用 Lieb 的变分公式
    $$\operatorname{tr}(K^* X K Y) = \min_{Z \succ 0} \left[\operatorname{tr}(K^* X K Z) + \operatorname{tr}(K^* Z^{-1} K^{-*} Y)\right] \cdot (\text{某个表达式})$$
    或更直接地，利用 Epstein 的矩阵 Herglotz 函数方法，将联合凹性归结为矩阵值 Pick 函数的性质。

    完整的证明需要矩阵值解析函数论的工具，超出了本章的范围，但核心思想是：$A^p$ 的积分表示将问题约化为分式线性函数 $A(A+tI)^{-1}$ 的凹性，后者可以直接验证。 $\blacksquare$

!!! example "例 46A.8 (Lieb 凹性定理的特殊情形)"
    1. **$K = I$，$n = m$**：$\operatorname{tr}(A^p B^{1-p})$ 在 $(A, B) \succ 0$ 上联合凹。这蕴含了 **Araki-Lieb-Thirring 不等式**的一个形式。

    2. **$p \to 0$**：$\operatorname{tr}(K^* K B) = \operatorname{tr}(K^* K) \cdot \operatorname{tr}(B)$... 不，这不对。$A^p \to I$ 当 $p \to 0$，因此 $\operatorname{tr}(K^* K B^1) = \operatorname{tr}(K^* K B)$，关于 $B$ 线性（当然凹）。

    3. **$B = I$**：$\operatorname{tr}(K^* A^p K)$ 关于 $A$ 凹（$0 < p < 1$）。这就是说 $A \mapsto K^* A^p K$ 的迹是凹函数，即 $A^p$ 的"压缩"的迹是凹的。

---

## 46A.14 Lieb-Ruskai 定理

<div class="context-flow" markdown>

**核心问题**：Lieb 凹性定理在量子信息论中有什么应用？

</div>

!!! theorem "定理 46A.14 (Lieb-Ruskai 定理 — 量子熵的强次可加性)"
    设 $\rho_{ABC}$ 是三部分量子系统 $\mathcal{H}_A \otimes \mathcal{H}_B \otimes \mathcal{H}_C$ 上的密度矩阵。von Neumann 熵 $S(\rho) = -\operatorname{tr}(\rho \log \rho)$ 满足**强次可加性**（strong subadditivity, SSA）：

    $$S(\rho_{ABC}) + S(\rho_B) \leq S(\rho_{AB}) + S(\rho_{BC}),$$

    其中 $\rho_{AB} = \operatorname{tr}_C(\rho_{ABC})$，$\rho_{BC} = \operatorname{tr}_A(\rho_{ABC})$，$\rho_B = \operatorname{tr}_{AC}(\rho_{ABC})$ 是约化密度矩阵。

    等价地，**条件量子互信息**非负：
    $$I(A:C|B) = S(\rho_{AB}) + S(\rho_{BC}) - S(\rho_{ABC}) - S(\rho_B) \geq 0.$$

??? proof "证明"
    **从 Lieb 凹性定理推导**（Lieb-Ruskai, 1973）：

    **步骤 1**：定义**量子相对熵** $D(\rho \| \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma)$（$\rho, \sigma$ 正半定，$\operatorname{supp}(\rho) \subseteq \operatorname{supp}(\sigma)$）。

    **步骤 2**：Lieb 凹性定理（令 $p \to 0$）蕴含**相对熵的联合凸性**：
    $$D(\lambda \rho_1 + (1-\lambda)\rho_2 \| \lambda \sigma_1 + (1-\lambda)\sigma_2) \leq \lambda D(\rho_1 \| \sigma_1) + (1-\lambda) D(\rho_2 \| \sigma_2).$$

    具体地，$D(\rho \| \sigma) = \lim_{p \to 0^+} \frac{1}{p}[\operatorname{tr}(\rho) - \operatorname{tr}(\rho^{1-p} \sigma^p)]$，而 Lieb 凹性定理保证 $\operatorname{tr}(\rho^{1-p} \sigma^p)$ 联合凹，因此 $D(\rho \| \sigma)$ 联合凸。

    **步骤 3**：相对熵的联合凸性蕴含**单调性**（数据处理不等式）：对量子信道 $\Phi$（完全正保迹映射），
    $$D(\Phi(\rho) \| \Phi(\sigma)) \leq D(\rho \| \sigma).$$

    **步骤 4**：强次可加性等价于相对熵在部分迹下的单调性。具体地，取 $\rho = \rho_{ABC}$，$\sigma = I_A \otimes \rho_{BC} / d_A$，$\Phi = \operatorname{tr}_C$，由单调性推导出 SSA。 $\blacksquare$

!!! example "例 46A.9 (强次可加性的信息论意义)"
    强次可加性 $S(\rho_{ABC}) + S(\rho_B) \leq S(\rho_{AB}) + S(\rho_{BC})$ 在量子信息论中意味着：

    1. **量子信息不可能被"免费复制"**：观察系统 $C$ 不会增加 $A$ 与 $B$ 之间的关联。

    2. **量子信道的容量公式**：许多量子通信容量的证明都以 SSA 为关键步骤。

    3. **量子纠错**：量子纠错码的存在性条件与 SSA 的等号条件（量子 Markov 链）密切相关。

    Lieb 和 Ruskai 在 1973 年利用 Lieb 凹性定理（定理 46A.13）证明了 SSA，解决了量子统计力学中一个长期悬而未决的猜想。

---

## 46A.15 习题

!!! example "例 46A.10"
    证明逆函数 $f(t) = t^{-1}$ 在 $(0, \infty)$ 上是算子反单调的：$A \succeq B \succ 0 \Rightarrow B^{-1} \succeq A^{-1}$。

    **提示**：直接计算 $B^{-1} - A^{-1} = B^{-1}(A - B)A^{-1}$，然后利用 Schur 补或直接验证半正定性。

!!! example "例 46A.11"
    验证 $f(t) = \frac{t}{1+t}$ 是 Pick 函数，并写出其积分表示。

!!! example "例 46A.12"
    计算 $f(t) = \log t$ 在点 $t_1 = 1, t_2 = e, t_3 = e^2$ 处的 Löwner 矩阵，并验证其半正定性。

!!! example "例 46A.13"
    证明：若 $f$ 和 $g$ 都是算子单调的，则 $f + g$ 也是算子单调的。$f \cdot g$ 是否一定是算子单调的？给出证明或反例。

!!! example "例 46A.14"
    设 $A = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。验证 $A \succeq B \succ 0$，然后计算 $A^{1/2}$ 和 $B^{1/2}$，验证 $A^{1/2} \succeq B^{1/2}$（Löwner-Heinz 不等式）。

!!! example "例 46A.15"
    证明 $f(t) = t^2$ 不是算子凹的，即存在 $A, B \succ 0$ 和 $\lambda \in (0, 1)$ 使得
    $$(\lambda A + (1-\lambda)B)^2 \not\succeq \lambda A^2 + (1-\lambda)B^2.$$
    （注意：$t^2$ 是算子凸的，因此反向不等式成立。）

!!! example "例 46A.16"
    利用 Jensen 算子不等式（定理 46A.10），证明：对等距算子 $V$（$V^*V = I$）和算子凹函数 $f$，有
    $$f(V^*AV) \succeq V^*f(A)V.$$

!!! example "例 46A.17"
    证明 Furuta 不等式在 $r = 0$ 时退化为 Löwner-Heinz 不等式。

---

## 练习题

****

??? success "参考答案"
    
       注意到 $B^{-1} - A^{-1} = B^{-1}(A - B)A^{-1}$。令 $C = A-B \succeq 0$。
       由 $B \preceq A$ 可得 $A^{-1/2} B A^{-1/2} \preceq I$，即最大特征值 $\le 1$。
       其逆矩阵 $A^{1/2} B^{-1} A^{1/2}$ 的特征值均 $\ge 1$，故 $A^{1/2} B^{-1} A^{1/2} \succeq I$。
       两边同时乘 $A^{-1/2}$ 即得 $B^{-1} \succeq A^{-1}$。

****

??? success "参考答案"
    
       取 $z = -1 + i$（在上半平面），则 $z^2 = (-1+i)^2 = 1 - 1 - 2i = -2i$。
       其虚部为 $-2 < 0$，故 $t^2$ 不是 Pick 函数。

****

??? success "参考答案"
    
       $L_{22} = f'(4) = \frac{1}{2\sqrt{4}} = 0.25$。
       $L_{12} = L_{21} = \frac{\sqrt{4}-\sqrt{1}}{4-1} = \frac{2-1}{3} = 1/3$。
       故 $L_2 = \begin{pmatrix} 0.5 & 1/3 \\ 1/3 & 0.25 \end{pmatrix}$。计算行列式 $\det L_2 = 0.125 - 1/9 > 0$，是正定的。

****

??? success "参考答案"
    
       展开左式：$(\lambda A + (1-\lambda)B)^2 = \lambda^2 A^2 + (1-\lambda)^2 B^2 + \lambda(1-\lambda)(AB + BA)$。
       不等式等价于检查 $\lambda(1-\lambda)(A-B)^2 \succeq 0$，由于对于任何 Hermite 矩阵，平方项 $(A-B)^2$ 总是半正定的，故成立。

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

本章探讨了函数映射下矩阵结构的保持性质：

1. **序的保持**：定义了算子单调性，并识别出 $t^r$ 在 $r > 1$ 时的阶梯壁垒。
2. **解析解析**：利用 Löwner 定理将算子单调性与复分析中的 Pick 函数理论联系起来。
3. **凸性微积分**：发展了算子凸性和针对矩阵算子的 Jensen 型不等式。
4. **信息界限**：应用这些结论确立了迹函数的联合凹性，从而证明了量子力学中的强次可加性。

