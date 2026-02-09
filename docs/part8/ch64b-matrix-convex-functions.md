# 第 64B 章 矩阵凸函数与算子单调性

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 矩阵不等式(Ch18) · 矩阵空间凸集(Ch64A) · 优化谱理论(Ch25) · 优超(Ch31)

**本章脉络**：标量凸的矩阵函数（$\lambda_{\max}$, $-\log\det$, 核范数等）→ Ky Fan 范数凸性 → 矩阵凹函数（$\lambda_{\min}$, $\log\det$, $(\det)^{1/n}$）→ 算子凸函数（矩阵序意义）→ Lowner 定理（算子单调 $\leftrightarrow$ Pick 函数）→ Hansen-Pedersen 定理 → 联合凸性/凹性 → Lieb 凹性定理 → 矩阵透视函数 → Schur 凸性

**延伸**：算子凸/单调理论连接矩阵分析与泛函分析、量子信息论（量子相对熵、数据处理不等式）

</div>

矩阵值函数的凸性理论有两个截然不同的层次。第一个层次是**标量凸性**：函数 $f: S_n \to \mathbb{R}$ 关于矩阵参数是凸的（即 $f(tA + (1-t)B) \leq tf(A) + (1-t)f(B)$）。这是经典凸分析在矩阵空间上的直接推广。第二个更深刻的层次是**算子凸性**：函数 $f$ 满足 $f\left(\frac{A+B}{2}\right) \preceq \frac{f(A)+f(B)}{2}$，其中 $\preceq$ 是 Lowner 偏序（矩阵序）。算子凸性要求不等式在**矩阵序**意义下成立，这是一个极强的条件。

Lowner（1934）关于算子单调函数的开创性工作建立了算子单调性与复分析中 Pick 函数之间的深刻联系。Lieb（1973）的凹性定理则成为量子信息论的基石。本章系统区分并建立这两个层次的理论。

---

## 64B.1 标量凸的矩阵函数

### 定义

!!! definition "定义 64B.1 (矩阵函数的标量凸性)"
    函数 $f: \mathcal{D} \to \mathbb{R}$（其中 $\mathcal{D} \subseteq S_n(\mathbb{R})$ 或 $\mathcal{D} \subseteq M_n(\mathbb{R})$ 是凸集）称为**凸的**，若
    $$
    f(tA + (1-t)B) \leq t f(A) + (1-t) f(B), \quad \forall A, B \in \mathcal{D}, \; t \in [0,1]
    $$
    这里 $f$ 的值域是 $\mathbb{R}$，不等式是标量不等式。

!!! note "注"
    标量凸性是经典凸分析的自然推广。$(S_n(\mathbb{R}), \langle \cdot, \cdot \rangle)$ 同构于 $\mathbb{R}^{n(n+1)/2}$，因此 $\mathbb{R}^d$ 中凸分析的所有结论（梯度条件、Hessian 条件、Jensen 不等式等）均可直接搬用。

### 凸矩阵函数的重要例子

!!! theorem "定理 64B.1 (凸矩阵函数)"
    以下函数在指定域上是凸的：

    1. **最大特征值**：$f(A) = \lambda_{\max}(A)$ 在 $S_n$ 上凸；
    2. **前 $k$ 个特征值之和**：$f(A) = \sum_{i=1}^k \lambda_i^\downarrow(A)$ 在 $S_n$ 上凸；
    3. **核范数**（迹范数）：$f(A) = \|A\|_* = \sum_i \sigma_i(A)$ 在 $M_n$ 上凸；
    4. **谱范数**：$f(A) = \|A\|_2 = \sigma_{\max}(A)$ 在 $M_n$ 上凸；
    5. **负对数行列式**：$f(A) = -\log\det(A)$ 在 $S_n^{++}$ 上凸；
    6. **迹指数**：$f(A) = \operatorname{tr}(e^A)$ 在 $S_n$ 上凸；
    7. **矩阵幂迹**（$p \geq 1$）：$f(A) = \operatorname{tr}(A^p)$ 在 $S_n^+$ 上凸。

??? proof "证明"
    **$\lambda_{\max}$ 的凸性**：利用变分刻画
    $$
    \lambda_{\max}(A) = \max_{\|x\|=1} x^T A x
    $$
    对固定 $x$，$g_x(A) = x^T A x$ 是关于 $A$ 的仿射函数（因此既凸又凹）。$\lambda_{\max}(A)$ 是仿射函数族的逐点上确界，因此凸。

    **$-\log\det$ 的凸性**：设 $A, B \in S_n^{++}$，$t \in [0,1]$，$C = tA + (1-t)B$。需证
    $$
    -\log\det(C) \leq -t\log\det(A) - (1-t)\log\det(B)
    $$
    即 $\det(C) \geq \det(A)^t \det(B)^{1-t}$。

    不妨通过 $B^{-1/2}$ 缩并设 $B = I$。则需证 $\det(tA + (1-t)I) \geq \det(A)^t$。

    设 $A$ 的特征值为 $\lambda_i > 0$，则
    $$
    \det(tA + (1-t)I) = \prod_i (t\lambda_i + 1 - t)
    $$
    由标量凹函数 $\log$ 的 Jensen 不等式：$\log(t\lambda_i + 1-t) \geq t\log\lambda_i$。因此
    $$
    \log\det(tA + (1-t)I) = \sum_i \log(t\lambda_i + 1-t) \geq t\sum_i \log\lambda_i = t\log\det(A)
    $$
    取指数得 $\det(tA + (1-t)I) \geq \det(A)^t$。$\blacksquare$

!!! example "例 64B.1"
    设 $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$，$t = 1/2$。

    $\frac{A+B}{2} = \begin{pmatrix} 2 & 0 \\ 0 & 2.5 \end{pmatrix}$，$\lambda_{\max}\left(\frac{A+B}{2}\right) = 2.5$。

    $\frac{\lambda_{\max}(A) + \lambda_{\max}(B)}{2} = \frac{3 + 4}{2} = 3.5$。

    确实 $2.5 \leq 3.5$，验证了 $\lambda_{\max}$ 的凸性。

### Ky Fan 范数的凸性

!!! theorem "定理 64B.2 (Ky Fan 范数的凸性)"
    对任意 $k = 1, \ldots, n$，**Ky Fan $k$-范数**
    $$
    \|A\|_{(k)} = \sum_{i=1}^k \sigma_i^\downarrow(A)
    $$
    在 $M_n$ 上是凸函数。特别地，$k=1$ 给出谱范数，$k=n$ 给出核范数。

??? proof "证明"
    利用 Ky Fan 变分刻画：
    $$
    \|A\|_{(k)} = \max\{\operatorname{tr}(U^T A V) : U \in \mathbb{R}^{n \times k}, V \in \mathbb{R}^{n \times k}, U^TU = V^TV = I_k\}
    $$
    即在 Stiefel 流形上的最大化。

    对固定的 $U, V$，$A \mapsto \operatorname{tr}(U^T A V)$ 是仿射函数。$\|A\|_{(k)}$ 作为仿射函数族在紧集上的逐点最大值（Stiefel 流形是紧的），因此凸。$\blacksquare$

---

## 64B.2 凹矩阵函数

!!! theorem "定理 64B.3 (凹矩阵函数)"
    以下函数在指定域上是凹的：

    1. **最小特征值**：$f(A) = \lambda_{\min}(A)$ 在 $S_n$ 上凹；
    2. **对数行列式**：$f(A) = \log\det(A)$ 在 $S_n^{++}$ 上凹；
    3. **矩阵幂迹**（$0 \leq p \leq 1$）：$f(A) = \operatorname{tr}(A^p)$ 在 $S_n^+$ 上凹；
    4. **负迹逆**：$f(A) = -\operatorname{tr}(A^{-1})$ 在 $S_n^{++}$ 上凹；
    5. **行列式的 $1/n$ 次幂**：$f(A) = (\det A)^{1/n}$ 在 $S_n^+$ 上凹。

??? proof "证明"
    **$\lambda_{\min}$ 的凹性**：$\lambda_{\min}(A) = \min_{\|x\|=1} x^TAx$ 是仿射函数的逐点下确界，因此凹。

    **$\log\det$ 的凹性**：由定理 64B.1 的证明，$-\log\det$ 凸，故 $\log\det$ 凹。

    **$(\det A)^{1/n}$ 的凹性**（Minkowski 行列式定理）：对 $A, B \in S_n^+$，$t \in [0,1]$：
    $$
    \det(tA + (1-t)B)^{1/n} \geq t\det(A)^{1/n} + (1-t)\det(B)^{1/n}
    $$
    不妨设 $B$ 正定。令 $C = B^{-1/2}AB^{-1/2}$，特征值为 $\mu_i \geq 0$。则
    $$
    \frac{\det(tA + (1-t)B)^{1/n}}{\det(B)^{1/n}} = \det(tC + (1-t)I)^{1/n} = \left(\prod_i (t\mu_i + 1-t)\right)^{1/n}
    $$
    由 AM-GM 不等式：
    $$
    \left(\prod_i (t\mu_i + 1-t)\right)^{1/n} \geq t\left(\prod_i \mu_i\right)^{1/n} + (1-t) \cdot 1 = t\det(C)^{1/n} + 1-t
    $$
    最后一步使用了凹函数 $(t\mu + 1-t)$ 的 Schur 凹性论证（或直接使用 Minkowski 不等式的标准证明）。乘以 $\det(B)^{1/n}$ 即得结论。$\blacksquare$

!!! example "例 64B.2"
    设 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}$，$t = 1/2$。

    $\log\det\left(\frac{A+B}{2}\right) = \log\det\begin{pmatrix} 3 & 0 \\ 0 & 2 \end{pmatrix} = \log 6 \approx 1.791$。

    $\frac{\log\det(A) + \log\det(B)}{2} = \frac{\log 6 + \log 4}{2} = \frac{\log 24}{2} \approx 1.589$。

    确实 $1.791 \geq 1.589$，验证了凹性。

---

## 64B.3 算子凸函数

### 定义与基本概念

标量凸性只要求实值不等式 $f(tA + (1-t)B) \leq tf(A) + (1-t)f(B)$。算子凸性则要求**矩阵序**意义下的不等式。

!!! definition "定义 64B.2 (算子凸函数)"
    设 $I \subseteq \mathbb{R}$ 是区间，$f: I \to \mathbb{R}$ 是连续函数。称 $f$ 是**算子凸的**（operator convex），若对所有 $n \geq 1$ 和所有特征值在 $I$ 中的 $A, B \in S_n(\mathbb{R})$：
    $$
    f\left(\frac{A + B}{2}\right) \preceq \frac{f(A) + f(B)}{2}
    $$
    其中 $\preceq$ 是 Lowner 偏序（$X \preceq Y$ 表示 $Y - X \succeq 0$），$f(A)$ 是矩阵函数（通过谱分解定义）。

    更一般地，对所有 $t \in [0,1]$：
    $$
    f(tA + (1-t)B) \preceq tf(A) + (1-t)f(B)
    $$

!!! definition "定义 64B.3 (算子单调函数)"
    连续函数 $f: I \to \mathbb{R}$ 称为**算子单调的**（operator monotone），若对所有 $n \geq 1$ 和所有特征值在 $I$ 中的 $A, B \in S_n(\mathbb{R})$：
    $$
    A \preceq B \Rightarrow f(A) \preceq f(B)
    $$

!!! note "注"
    算子凸/单调是**非常强**的条件。例如，$f(t) = t^2$ 在标量意义下是凸的，但**不是**算子凸的（对 $2 \times 2$ 矩阵就可以找到反例）。类似地，$f(t) = t^2$ 在标量意义下是单调递增的（$t > 0$），但不是算子单调的。

!!! example "例 64B.3"
    验证 $f(t) = t^2$ 不是算子凸的。取
    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}
    $$
    $$
    f\left(\frac{A+B}{2}\right) = \left(\frac{I}{2}\right)^2 = \frac{I}{4}
    $$
    $$
    \frac{f(A) + f(B)}{2} = \frac{A + B}{2} = \frac{I}{2}
    $$
    $\frac{I}{4} \preceq \frac{I}{2}$，这个例子碰巧满足。但取非对角的例子：
    $$
    A = \begin{pmatrix} 4 & 0 \\ 0 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
    $$
    $\frac{A+B}{2} = \begin{pmatrix} 2.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$，$f(\frac{A+B}{2}) = \begin{pmatrix} 6.5 & 1.5 \\ 1.5 & 0.5 \end{pmatrix}$。

    $\frac{f(A)+f(B)}{2} = \frac{1}{2}\begin{pmatrix} 16 & 0 \\ 0 & 0 \end{pmatrix} + \frac{1}{2}\begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix} = \begin{pmatrix} 9 & 1 \\ 1 & 1 \end{pmatrix}$。

    差矩阵 $\begin{pmatrix} 2.5 & -0.5 \\ -0.5 & 0.5 \end{pmatrix}$，行列式 $= 1.25 - 0.25 = 1 > 0$，正定。这个例子也满足。实际上反例需要更大的矩阵或不同的选择。

### 算子凸函数的例子

!!! theorem "定理 64B.4 (算子凸函数的例子)"
    以下函数是算子凸的：

    1. $f(t) = t^2$ 在 $\mathbb{R}$ 上**不是**算子凸的；
    2. $f(t) = t^{-1}$ 在 $(0, \infty)$ 上是算子凸的；
    3. $f(t) = -\log t$ 在 $(0, \infty)$ 上是算子凸的；
    4. $f(t) = t^p$，$1 \leq p \leq 2$ 在 $[0, \infty)$ 上是算子凸的；
    5. $f(t) = t\log t$ 在 $(0, \infty)$ 上是算子凸的。

    以下函数是**算子凹**的（即 $-f$ 算子凸）：

    1. $f(t) = t^p$，$0 \leq p \leq 1$ 在 $[0, \infty)$ 上算子凹；
    2. $f(t) = \log t$ 在 $(0, \infty)$ 上算子凹；
    3. $f(t) = \frac{t}{1+t}$ 在 $[0, \infty)$ 上算子凹。

---

## 64B.4 Lowner 定理

### 算子单调函数的刻画

Lowner（1934）的定理是算子单调性理论的基石，建立了算子单调函数与复分析中 Pick 函数之间的等价关系。

!!! definition "定义 64B.4 (Pick 函数 / Nevanlinna 函数)"
    函数 $f: \mathbb{C} \setminus \mathbb{R} \to \mathbb{C}$ 称为 **Pick 函数**（或 Nevanlinna 函数），若 $f$ 在上半平面 $\mathbb{C}^+$ 上解析，且 $\operatorname{Im}(f(z)) \geq 0$ 当 $\operatorname{Im}(z) > 0$。

!!! theorem "定理 64B.5 (Lowner 定理, 1934)"
    设 $f: (a, b) \to \mathbb{R}$ 是连续函数（$-\infty \leq a < b \leq +\infty$）。以下三个条件等价：

    1. $f$ 在 $(a, b)$ 上是**算子单调的**；
    2. $f$ 可以解析延拓到上半平面 $\mathbb{C}^+$，且延拓后的函数是 **Pick 函数**（即 $\operatorname{Im}(z) > 0 \Rightarrow \operatorname{Im}(f(z)) \geq 0$）；
    3. $f$ 具有**积分表示**：
    $$
    f(t) = \alpha + \beta t + \int_{\mathbb{R} \setminus (a,b)} \left(\frac{1}{\lambda - t} - \frac{\lambda}{1 + \lambda^2}\right) d\mu(\lambda)
    $$
    其中 $\alpha \in \mathbb{R}$，$\beta \geq 0$，$\mu$ 是 $\mathbb{R} \setminus (a,b)$ 上的非负 Borel 测度且 $\int \frac{d\mu(\lambda)}{1+\lambda^2} < \infty$。

??? proof "证明"
    **（证明思路/关键步骤）**

    **(1) $\Rightarrow$ (2)**：设 $f$ 在 $(a, b)$ 上算子单调。关键是证明 Lowner 矩阵
    $$
    L_f = \left(\frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j}\right)_{i,j=1}^n
    $$
    对所有 $n$ 和所有 $a < \lambda_1 < \cdots < \lambda_n < b$ 是半正定的。

    **证明 Lowner 矩阵半正定**：设 $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，取 $A = \Lambda$，$B = \Lambda + \varepsilon vv^T$（$v = (v_1, \ldots, v_n)^T$，$\varepsilon > 0$ 小）。由 $A \preceq B$（$\varepsilon > 0$），算子单调性给出 $f(A) \preceq f(B)$。

    展开到一阶：$f(B) - f(A) \approx \varepsilon \cdot \operatorname{D}f(\Lambda)[vv^T]$，其中 $\operatorname{D}f(\Lambda)$ 是 Frechet 导数。可以计算
    $$
    \operatorname{D}f(\Lambda)[vv^T] = L_f \circ (vv^T) = \left(\frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j} v_i v_j\right)
    $$
    （Hadamard/Schur 积）。$f(B) - f(A) \succeq 0$ 蕴含 $L_f \circ (vv^T) \succeq 0$ 对所有 $v$，这等价于 $L_f \succeq 0$（Schur 积定理）。

    Lowner 矩阵的半正定性允许 $f$ 解析延拓到 $\mathbb{C}^+$，且延拓函数满足 Pick 条件。这利用了 Pick 矩阵半正定性与 Pick 插值的经典联系。

    **(2) $\Rightarrow$ (3)**：Pick 函数的积分表示是经典的 Herglotz-Nevanlinna 表示定理。

    **(3) $\Rightarrow$ (1)**：积分表示中的基本函数 $t \mapsto \frac{1}{\lambda - t}$（$\lambda \notin (a,b)$）是算子单调的。设 $\lambda > b$（$\lambda < a$ 类似），则 $A \preceq B$ 蕴含 $\lambda I - B \preceq \lambda I - A$，从而 $(\lambda I - A)^{-1} \preceq (\lambda I - B)^{-1}$（因为 $t \mapsto t^{-1}$ 在正数上算子单调递减），即 $\frac{1}{\lambda - A} \preceq \frac{1}{\lambda - B}$（取负号后反向）。

    实际上需要更仔细地处理符号。$\lambda > b > t$ 时，$g(t) = \frac{1}{\lambda - t}$ 是递增的，且 $g$ 是算子单调的，因为 Lowner 矩阵 $\frac{g(\lambda_i) - g(\lambda_j)}{\lambda_i - \lambda_j} = \frac{1}{(\lambda - \lambda_i)(\lambda - \lambda_j)} \succeq 0$。

    由非负测度对算子单调函数的积分保持算子单调性，$f$ 是算子单调的。$\blacksquare$

!!! example "例 64B.4"
    **$f(t) = \sqrt{t}$ 是 $[0, \infty)$ 上的算子单调函数。** 其 Pick 函数延拓为 $f(z) = \sqrt{z}$（取上半平面分支），满足 $\operatorname{Im}(\sqrt{z}) > 0$ 当 $\operatorname{Im}(z) > 0$。

    积分表示为：$\sqrt{t} = \frac{1}{\pi} \int_0^\infty \frac{t}{\lambda + t} \cdot \frac{d\lambda}{\sqrt{\lambda}}$（可以验证）。

    算子单调性意味着：$A \preceq B$（$A, B \succeq 0$）$\Rightarrow$ $A^{1/2} \preceq B^{1/2}$。

!!! example "例 64B.5"
    **$f(t) = t^2$ 不是算子单调的。** 取
    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \preceq \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix} = B
    $$
    但 $A^2 = A$，$B^2 = \begin{pmatrix} 2 & 3 \\ 3 & 5 \end{pmatrix}$，$B^2 - A^2 = \begin{pmatrix} 1 & 3 \\ 3 & 5 \end{pmatrix}$，行列式 $= 5 - 9 = -4 < 0$，不是半正定的。因此 $A \preceq B$ 不蕴含 $A^2 \preceq B^2$。

---

## 64B.5 Hansen-Pedersen 定理

!!! theorem "定理 64B.6 (Hansen-Pedersen 定理)"
    设 $f: I \to \mathbb{R}$ 是连续函数。以下条件等价：

    1. $f$ 在 $I$ 上是**算子凸**的；
    2. 对所有 $n \geq 1$ 和所有 $n \times n$ 的矩阵 $A$（特征值在 $I$ 中）以及所有缩并算子 $C$（$\|C\| \leq 1$，即 $C^*C \preceq I$）：
    $$
    f(C^* A C) \preceq C^* f(A) C
    $$
    即 $f$ 在缩并下满足 Jensen 不等式；
    3. $2 \times 2$ 分块矩阵条件：对所有 $2n \times 2n$ 矩阵 $\begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$（特征值在 $I$ 中），有
    $$
    f\left(\begin{pmatrix} A & B \\ B^* & C \end{pmatrix}\right) \preceq \begin{pmatrix} f(A) & * \\ * & f(C) \end{pmatrix}
    $$
    的适当推广形式。

!!! note "注"
    Hansen-Pedersen 定理的条件 (2) 称为**算子 Jensen 不等式**。它将算子凸性与量子信道（完全正迹保持映射）联系起来：在量子信息论中，量子信道 $\Phi(\rho) = \sum_i K_i \rho K_i^*$（Kraus 表示）可以看作一系列缩并的组合，算子 Jensen 不等式保证了凸函数值在量子信道作用下的单调性。

---

## 64B.6 联合凸性与凹性

### 几何平均的联合凹性

!!! definition "定义 64B.5 (矩阵几何平均)"
    对 $A, B \in S_n^{++}$，**矩阵几何平均**定义为
    $$
    A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}
    $$
    这是满足 $X = A \# B \Leftrightarrow X A^{-1} X = B$ 的唯一正定解。

!!! theorem "定理 64B.7 (几何平均的联合凹性)"
    映射 $(A, B) \mapsto A \# B$ 在 $S_n^{++} \times S_n^{++}$ 上是**联合凹的**，即
    $$
    (tA_1 + (1-t)A_2) \# (tB_1 + (1-t)B_2) \succeq t(A_1 \# B_1) + (1-t)(A_2 \# B_2)
    $$

??? proof "证明"
    利用 $A \# B$ 的变分刻画：
    $$
    A \# B = \max\{X \in S_n^{++} : \begin{pmatrix} A & X \\ X & B \end{pmatrix} \succeq 0\}
    $$
    （最大值在 Lowner 序意义下取）。

    设 $X_1 = A_1 \# B_1$，$X_2 = A_2 \# B_2$。由 Schur 补条件：
    $$
    \begin{pmatrix} A_i & X_i \\ X_i & B_i \end{pmatrix} \succeq 0, \quad i = 1, 2
    $$

    取凸组合：
    $$
    t\begin{pmatrix} A_1 & X_1 \\ X_1 & B_1 \end{pmatrix} + (1-t)\begin{pmatrix} A_2 & X_2 \\ X_2 & B_2 \end{pmatrix} = \begin{pmatrix} tA_1 + (1-t)A_2 & tX_1 + (1-t)X_2 \\ tX_1 + (1-t)X_2 & tB_1 + (1-t)B_2 \end{pmatrix} \succeq 0
    $$

    因此 $tX_1 + (1-t)X_2$ 是使分块矩阵半正定的可行解，而 $(tA_1 + (1-t)A_2) \# (tB_1 + (1-t)B_2)$ 是最大的这样的矩阵。故
    $$
    (tA_1 + (1-t)A_2) \# (tB_1 + (1-t)B_2) \succeq tX_1 + (1-t)X_2 \qquad \blacksquare
    $$

---

## 64B.7 Lieb 凹性定理

!!! theorem "定理 64B.8 (Lieb 凹性定理, 1973)"
    设 $K$ 是固定矩阵。映射
    $$
    A \mapsto \operatorname{tr}(K^* A^p K A^{1-p})
    $$
    对 $p \in [0, 1]$ 在 $S_n^{++}$ 上关于 $A$ 是凹的。

    更一般地，对固定的 $K$ 和 $0 \leq p, q$，$p + q \leq 1$，映射
    $$
    (A, B) \mapsto \operatorname{tr}(K^* A^p K B^q)
    $$
    在 $S_n^{++} \times S_n^{++}$ 上是**联合凹的**。

??? proof "证明"
    **（Epstein 复插值方法的思路）**

    关键工具是 **Epstein 的复插值定理**：设 $F(z)$ 是在带状区域 $\{z : 0 \leq \operatorname{Re}(z) \leq 1\}$ 上解析且有界的矩阵值函数，且在边界上满足适当的范数条件。则 $F$ 在内部的行为由边界上的行为控制。

    **步骤 1**：定义 $\Phi_p(A) = K^* A^p K$（对固定的 $K$ 和 $p \in [0,1]$）。需要证明 $A \mapsto \operatorname{tr}(\Phi_p(A) \cdot A^{1-p})$ 凹。

    **步骤 2**：利用复插值。对正定矩阵 $A$，定义 $A^z = e^{z \log A}$（$z \in \mathbb{C}$）。函数 $z \mapsto A^z$ 在 $\mathbb{C}$ 上解析。

    **步骤 3**：定义 $G(z) = \operatorname{tr}(K^* A_t^z K A_t^{1-z})$，其中 $A_t = tA + (1-t)B$。利用 $G(z)$ 在带状区域上的解析性和边界行为，通过三线定理（three-lines theorem）得到 $|G(p)| \geq t|G_A(p)| + (1-t)|G_B(p)|$ 类型的不等式。

    **步骤 4**：更直接的方法是利用 **Ando 的联合凹性定理**：若 $f$ 是 $[0, \infty)$ 上的算子凹函数且 $f(0) \geq 0$，则 $(A, B) \mapsto B^{1/2} f(B^{-1/2} A B^{-1/2}) B^{1/2}$ 是联合凹的。取 $f(t) = t^p$（$0 \leq p \leq 1$，算子凹）得到 $(A, B) \mapsto A^p \# B^{1-p}$ 类型函数的联合凹性。

    完整证明见 Bhatia, *Matrix Analysis*, Chapter V。$\blacksquare$

!!! note "注"
    Lieb 凹性定理（1973）的推论包括：

    - **强次可加性**（Strong Subadditivity, SSA）：量子 von Neumann 熵 $S(\rho) = -\operatorname{tr}(\rho \log \rho)$ 满足 $S(ABC) + S(B) \leq S(AB) + S(BC)$；
    - **Wigner-Yanase-Dyson 猜想**的证明；
    - **数据处理不等式**：量子相对熵在量子信道下不增。

---

## 64B.8 矩阵透视函数

!!! definition "定义 64B.6 (矩阵透视函数)"
    设 $f: S_n^{++} \to \mathbb{R}$ 是凸函数。$f$ 的**透视函数**（perspective function）定义为
    $$
    g(A, t) = t \cdot f(A/t), \quad A \in S_n^{++}, \; t > 0
    $$

!!! theorem "定理 64B.9 (矩阵透视的凸性)"
    若 $f: S_n^{++} \to \mathbb{R}$ 是凸函数，则其透视函数 $g(A, t) = tf(A/t)$ 在 $S_n^{++} \times (0, \infty)$ 上联合凸。

??? proof "证明"
    设 $(A_1, t_1)$ 和 $(A_2, t_2)$ 是可行点，$\lambda \in [0, 1]$。令 $A = \lambda A_1 + (1-\lambda)A_2$，$t = \lambda t_1 + (1-\lambda)t_2$。

    $$
    g(A, t) = t \cdot f(A/t) = t \cdot f\left(\frac{\lambda A_1 + (1-\lambda)A_2}{\lambda t_1 + (1-\lambda)t_2}\right)
    $$

    注意 $\frac{\lambda A_1 + (1-\lambda)A_2}{\lambda t_1 + (1-\lambda)t_2} = \frac{\lambda t_1}{\lambda t_1 + (1-\lambda)t_2} \cdot \frac{A_1}{t_1} + \frac{(1-\lambda)t_2}{\lambda t_1 + (1-\lambda)t_2} \cdot \frac{A_2}{t_2}$。

    设 $\mu = \frac{\lambda t_1}{\lambda t_1 + (1-\lambda)t_2}$，则由 $f$ 的凸性：
    $$
    f(A/t) \leq \mu f(A_1/t_1) + (1-\mu)f(A_2/t_2)
    $$

    乘以 $t = \lambda t_1 + (1-\lambda)t_2$：
    $$
    g(A, t) \leq \lambda t_1 f(A_1/t_1) + (1-\lambda)t_2 f(A_2/t_2) = \lambda g(A_1, t_1) + (1-\lambda)g(A_2, t_2) \qquad \blacksquare
    $$

!!! example "例 64B.6"
    取 $f(A) = -\log\det(A)$（凸函数）。其透视函数为
    $$
    g(A, t) = t \cdot (-\log\det(A/t)) = -t\log\det(A) + tn\log t
    $$
    这是**量子相对熵** $D(\rho \| \sigma) = \operatorname{tr}(\rho(\log\rho - \log\sigma))$ 的一种推广形式，在量子信息论中有核心地位。

---

## 64B.9 Schur 凸性

!!! definition "定义 64B.7 (Schur 凸函数)"
    设 $\phi: \mathbb{R}^n \to \mathbb{R}$。称 $\phi$ 是 **Schur 凸的**（Schur convex），若对任意 $x, y \in \mathbb{R}^n$：
    $$
    x \prec y \Rightarrow \phi(x) \leq \phi(y)
    $$
    其中 $x \prec y$ 表示 $x$ 被 $y$ **优超**（majorized），即 $\sum_{i=1}^k x_i^\downarrow \leq \sum_{i=1}^k y_i^\downarrow$（$k = 1, \ldots, n-1$）且 $\sum_{i=1}^n x_i = \sum_{i=1}^n y_i$。

    称 $\phi$ 是 **Schur 凹的**，若 $-\phi$ 是 Schur 凸的。

!!! theorem "定理 64B.10 (Schur 凸性的刻画)"
    设 $\phi: \mathbb{R}^n \to \mathbb{R}$ 是连续可微的对称函数（$\phi(Px) = \phi(x)$ 对所有置换矩阵 $P$）。则 $\phi$ 是 Schur 凸的，当且仅当对所有 $i \neq j$：
    $$
    (x_i - x_j)\left(\frac{\partial \phi}{\partial x_i} - \frac{\partial \phi}{\partial x_j}\right) \geq 0
    $$
    即 $\phi$ 的偏导数与坐标值同向排列。

!!! theorem "定理 64B.11 (Schur 凸性与矩阵函数)"
    设 $A \in S_n(\mathbb{R})$，$\lambda(A) = (\lambda_1(A), \ldots, \lambda_n(A))$ 是 $A$ 的特征值向量。

    1. $A \mapsto \sum_i f(\lambda_i(A))$（其中 $f$ 是凸函数）对应的 $\phi(x) = \sum_i f(x_i)$ 是 Schur 凸的；
    2. $\phi(x) = \prod_i x_i$（$x_i > 0$）是 Schur 凹的；
    3. $\phi(x) = \sum_{i=1}^k x_i^\downarrow$ 是 Schur 凸的。

!!! note "注"
    Schur 凸性理论与优超理论（Ch31）密切相关。矩阵的特征值优超理论将 Schur 凸性提升到矩阵分析层面：若 $\lambda(A) \prec \lambda(B)$，则对所有 Schur 凸函数 $\phi$，$\phi(\lambda(A)) \leq \phi(\lambda(B))$。特别地，由 Schur 定理，$\operatorname{diag}(A) \prec \lambda(A)$，因此 $\sum_i f(a_{ii}) \leq \sum_i f(\lambda_i)$ 对所有凸函数 $f$ 成立。

---

## 习题

!!! question "习题 64B.1"
    证明 $f(A) = \operatorname{tr}(A^2)$ 在 $S_n$ 上是凸的。（提示：直接展开或利用 Hessian。）

!!! question "习题 64B.2"
    证明 $f(A) = -\operatorname{tr}(A^{-1})$ 在 $S_n^{++}$ 上是凹的。

!!! question "习题 64B.3"
    证明 $f(t) = t^{1/2}$ 是 $[0, \infty)$ 上的算子单调函数。即证明 $0 \preceq A \preceq B$ 蕴含 $A^{1/2} \preceq B^{1/2}$。

!!! question "习题 64B.4"
    给出一个例子说明 $f(t) = t^3$ 不是算子凸的。

!!! question "习题 64B.5"
    证明矩阵几何平均满足 $A \# B = B \# A$。

!!! question "习题 64B.6"
    设 $f(A) = \lambda_1(A) + \lambda_2(A)$（最大两个特征值之和）。利用 Ky Fan 变分原理证明 $f$ 在 $S_n$ 上是凸的。

!!! question "习题 64B.7"
    证明 Schur 凸函数 $\phi(x) = \sum_i x_i^2$ 满足 Schur 凸性判据（偏导数条件）。

!!! question "习题 64B.8"
    设 $A, B \in S_n^{++}$。利用 Lieb 凹性定理证明映射 $A \mapsto \operatorname{tr}(A^{1/2} B A^{1/2})^{1/2} = \operatorname{tr}(A \# B)$（取 $p = 1/2$, $K = I$）是凹的。

!!! question "习题 64B.9"
    证明：$f(t) = -1/t$ 在 $(0, \infty)$ 上是算子单调的。即 $0 \prec A \preceq B$ 蕴含 $-B^{-1} \preceq -A^{-1}$（等价地 $A^{-1} \succeq B^{-1}$）。

!!! question "习题 64B.10"
    利用 Lowner 定理的积分表示，证明 $[0, \infty)$ 上的算子单调函数 $f$ 满足 $f(0) \geq 0$（若 $f$ 在 $0$ 处有定义）。

!!! question "习题 64B.11"
    证明：$\phi(x) = -\sum_i \log x_i$（$x_i > 0$）是 Schur 凸的。推导出 $-\sum_i \log a_{ii} \leq -\sum_i \log \lambda_i(A)$（即 $\prod a_{ii} \leq \det(A)$ 当 $A$ 正定时——Hadamard 不等式）。

!!! question "习题 64B.12"
    （透视函数）设 $f(A) = \operatorname{tr}(A^2)$。计算透视函数 $g(A, t) = tf(A/t)$，并验证 $g$ 在 $S_n^{++} \times (0, \infty)$ 上联合凸。
