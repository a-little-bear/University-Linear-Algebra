# 第 46B 章 矩阵均值与几何均值

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 算子单调函数与算子凸函数(Ch46A) · Löwner 偏序

**本章脉络**：Kubo-Ando 公理 → Kubo-Ando 定理（均值与算子单调函数的一一对应）→ 算子连接 → 经典均值（算术、调和、几何）→ 几何均值的等价定义（Riccati 方程、变分刻画）→ 几何均值的 10+ 性质 → AM-GM-HM 不等式 → 矩阵幂均值（Lim-Palfia）→ Ando-Li-Mathias 多变量几何均值 → Karcher 均值 → 正定流形的几何（Thompson 度量、Riemannian 结构、Cartan-Hadamard 定理）→ 正半定矩阵的推广

**延伸**：矩阵几何均值在医学影像（扩散张量 MRI 的张量平均）、雷达信号处理（协方差矩阵的 Riemannian 均值）和量子信息（量子态的保真度）中有直接应用

</div>

两个正实数 $a, b > 0$ 的几何均值 $\sqrt{ab}$ 是最自然的"中点"概念之一。但对正定矩阵 $A, B \succ 0$，$\sqrt{AB}$ 甚至不一定是 Hermite 的（因为 $AB$ 不一定 Hermite），更不必说正定了。如何为正定矩阵定义一个"正确的"几何均值，使其继承标量几何均值的所有好性质？

答案是 $A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$——这个看似复杂的公式是唯一满足 Kubo-Ando 公理化框架中"几何均值"角色的定义。Kubo 和 Ando 在 1980 年证明了一个深刻的定理：正定矩阵上的所有"合理的"均值运算与 $[0, \infty)$ 上归一化的算子单调函数之间存在一一对应。

本章从 Kubo-Ando 公理出发，建立矩阵均值的公理化理论，深入探讨几何均值的丰富性质，然后推广到多变量几何均值（Ando-Li-Mathias 均值和 Karcher 均值），最终进入正定矩阵流形的 Riemannian 几何。

---

## 46B.1 Kubo-Ando 矩阵均值公理

<div class="context-flow" markdown>

**核心问题**：如何公理化地定义两个正定矩阵的"均值"？什么样的运算才称得上是"均值"？

</div>

!!! definition "定义 46B.1 (Kubo-Ando 矩阵均值)"
    **矩阵均值**（matrix mean）是正定矩阵对 $(A, B)$ 上的二元运算 $\sigma$，$A \sigma B$ 仍是正定矩阵，满足以下五个公理：

    **(M1) 单调性**：$A \preceq A'$，$B \preceq B'$ $\Rightarrow$ $A \sigma B \preceq A' \sigma B'$。

    **(M2) 变换不等式（传递性）**：$C^*(A \sigma B)C \preceq (C^*AC) \sigma (C^*BC)$，对任意矩阵 $C$。

    **(M3) 从上连续**：若 $A_n \downarrow A$（单调递减趋于 $A$），$B_n \downarrow B$，则 $A_n \sigma B_n \downarrow A \sigma B$。

    **(M4) 归一化**：$I \sigma I = I$。

    **(M5) 对称性（可选）**：$A \sigma B = B \sigma A$。

    公理 (M1)-(M4) 定义**矩阵均值**；加上 (M5) 定义**对称矩阵均值**。

!!! example "例 46B.1 (公理的直觉)"
    - **(M1)**：输入增大，均值增大——这是"均值"最基本的要求。
    - **(M2)**：合同变换（坐标变换）的相容性。对可逆 $C$，等号成立：$C^*(A \sigma B)C = (C^*AC) \sigma (C^*BC)$。
    - **(M3)**：连续性保证极限存在。
    - **(M4)**：两个相同的单位矩阵的均值是单位矩阵。
    - **(M5)**：均值不依赖于输入的顺序。算术均值和几何均值满足 (M5)，但加权版本（$t \neq 1/2$）不满足。

---

## 46B.2 Kubo-Ando 定理

<div class="context-flow" markdown>

**核心问题**：满足 Kubo-Ando 公理的矩阵均值有多少？如何分类？

</div>

!!! theorem "定理 46B.1 (Kubo-Ando 定理)"
    矩阵均值 $\sigma$（满足公理 (M1)-(M4)）与 $[0, \infty)$ 上的算子单调函数 $f$（满足 $f(1) = 1$）之间存在**一一对应**：

    $$A \sigma B = A^{1/2} f(A^{-1/2} B A^{-1/2}) A^{1/2}.$$

    反过来，$f(t) = I \sigma (tI) = 1 \sigma t$（标量情形）。

    因此，矩阵均值的分类**完全等价于**归一化算子单调函数的分类。

??? proof "证明"
    **从均值到函数**：定义 $f(t) = 1 \sigma t$（$t > 0$ 上的标量函数）。

    由公理 (M4)，$f(1) = 1 \sigma 1 = 1$。

    由公理 (M1)，$t_1 \leq t_2 \Rightarrow t_1 I \preceq t_2 I \Rightarrow I \sigma (t_1 I) \preceq I \sigma (t_2 I)$，即 $f(t_1) \leq f(t_2)$。因此 $f$ 单调递增。

    由公理 (M3)，$f$ 从上连续。

    **证明 $f$ 算子单调**：设 $A \succeq B \succ 0$（$n \times n$ Hermite 矩阵）。需证 $f(A) \succeq f(B)$。

    由公理 (M2) 对可逆 $C$ 取等号：$C^*(I \sigma T)C = (C^*C) \sigma (C^*TC)$。取 $C = A^{1/2}$：
    $$A^{1/2}(I \sigma T)A^{1/2} = A \sigma (A^{1/2}TA^{1/2}).$$

    取 $T = A^{-1/2}BA^{-1/2}$：
    $$A^{1/2} f(A^{-1/2}BA^{-1/2}) A^{1/2} = A \sigma B.$$

    由公理 (M1)，$A \succeq B$ 和 $A \succeq A$ 蕴含... 不，我们需要固定一个参数。更精确地：

    取 $A_1 = A_2 = I$，$B_1 = A^{-1/2}BA^{-1/2} \preceq I$（因为 $B \preceq A$），$B_2 = I$。由 (M1)，
    $$I \sigma B_1 \preceq I \sigma I = I,$$
    即 $f(A^{-1/2}BA^{-1/2}) \preceq I = f(I)$。

    这只证明了 $f$ 保持 $\preceq I$ 的关系。完整的算子单调性证明需要更精细地利用 (M2) 的变换性质和 (M3) 的连续性。

    **从函数到均值**：给定算子单调函数 $f$（$f(1) = 1$），定义
    $$A \sigma B = A^{1/2} f(A^{-1/2}BA^{-1/2}) A^{1/2}.$$

    验证五个公理：

    - **(M1) 单调性**：设 $A \preceq A'$，$B \preceq B'$。由 $f$ 的算子单调性和合同变换的性质。具体地，$A^{-1/2}BA^{-1/2} \preceq (A')^{-1/2}B'(A')^{-1/2}$ 的推导需要一些代数（不是直接的），但可以通过 $A \sigma B$ 作为 Riccati 方程解的单调性来证明。

    - **(M2) 变换不等式**：对可逆 $C$，
    $$C^*(A \sigma B)C = C^* A^{1/2} f(A^{-1/2}BA^{-1/2}) A^{1/2} C.$$
    需要证明这等于 $(C^*AC) \sigma (C^*BC)$。设 $D = A^{1/2}C$，则 $C^*AC = D^*D$，
    $$(C^*AC) \sigma (C^*BC) = D^*D \cdot [\text{某个表达式}] = C^* A^{1/2} f(A^{-1/2}BA^{-1/2}) A^{1/2} C.$$
    这需要验证 $f(D^{-1}(C^*BC)(D^*)^{-1}) = f(A^{-1/2}BA^{-1/2})$... 实际上 $D^{-1}C^*BC(D^*)^{-1} = C^{-1}A^{-1/2} \cdot B \cdot A^{-1/2}(C^*)^{-1} \neq A^{-1/2}BA^{-1/2}$。

    正确的验证：$(C^*AC)^{1/2} = |A^{1/2}C|$（极分解的模），计算较为技术性。对可逆 $C$ 的情形，等号可以直接代数验证。

    - **(M3) 连续性**：由 $f$ 的连续性和矩阵运算的连续性。

    - **(M4) 归一化**：$I \sigma I = I^{1/2} f(I^{-1/2}II^{-1/2}) I^{1/2} = f(I) = f(1) \cdot I = I$。

    **一一对应**：不同的算子单调函数 $f_1 \neq f_2$ 给出不同的均值（取 $A = I$ 即得 $f_1(B) \neq f_2(B)$）。反过来，每个均值唯一确定 $f(t) = 1 \sigma t$。 $\blacksquare$

---

## 46B.3 算子连接

!!! definition "定义 46B.2 (算子连接)"
    **算子连接**（operator connection）是正定矩阵对上的二元运算 $\sigma$，满足公理 (M1)-(M3)，但不要求 (M4)（归一化）。

    算子连接与 $[0, \infty)$ 上的**非负算子单调函数**（不要求 $f(1) = 1$）之间存在一一对应：
    $$A \sigma B = A^{1/2} f(A^{-1/2}BA^{-1/2}) A^{1/2}, \quad f(t) = 1 \sigma t.$$

    矩阵均值是算子连接的特殊情形（加上归一化条件 $f(1) = 1$）。

!!! example "例 46B.2 (算子连接的例子)"
    1. **平行和**（parallel sum）：$A : B = (A^{-1} + B^{-1})^{-1}$，对应 $f(t) = t/(1+t)$。这不是均值（$f(1) = 1/2 \neq 1$），但是算子连接。

    2. **缩放均值**：$A \sigma_\alpha B = \alpha \cdot (A \# B)$，对应 $f(t) = \alpha \sqrt{t}$。$\alpha \neq 1$ 时不是均值。

    3. **左投影**：$A \sigma B = A$，对应 $f(t) = 1$。这是平凡的算子连接。

---

## 46B.4 经典矩阵均值

<div class="context-flow" markdown>

**核心问题**：算术、调和、几何均值在矩阵情形如何定义？对应哪些算子单调函数？

</div>

!!! definition "定义 46B.3 (矩阵算术均值)"
    $$A \nabla B = \frac{A + B}{2}.$$

    带权重 $t \in [0, 1]$ 的算术均值：
    $$A \nabla_t B = (1 - t)A + tB.$$

    对应的算子单调函数：$f(t) = (1 + t)/2$（对 $t = 1/2$）；$f_t(s) = 1 - t + ts$。

!!! definition "定义 46B.4 (矩阵调和均值)"
    $$A \,!\, B = 2(A^{-1} + B^{-1})^{-1}.$$

    等价形式：$A \,!\, B = A(A + B)^{-1}B + B(A + B)^{-1}A$。

    带权重版本：$A \,!_t\, B = \left[(1-t)A^{-1} + tB^{-1}\right]^{-1}$。

    对应的算子单调函数：$f(t) = 2t/(1 + t)$。

!!! definition "定义 46B.5 (矩阵几何均值)"
    对正定矩阵 $A, B \succ 0$，定义**矩阵几何均值**为

    $$A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}.$$

    带权重版本（$t \in [0, 1]$）：
    $$A \#_t B = A^{1/2}(A^{-1/2}BA^{-1/2})^t A^{1/2}.$$

    注意 $A \#_0 B = A$，$A \#_1 B = B$，$A \#_{1/2} B = A \# B$。

    对应的算子单调函数：$f(t) = t^{1/2}$（无权重）；$f_t(s) = s^t$（带权重）。

!!! definition "定义 46B.6 (矩阵对数均值)"
    $$A \,\sharp_{\log}\, B = A^{1/2} \frac{A^{-1/2}BA^{-1/2} - I}{\log(A^{-1/2}BA^{-1/2})} A^{1/2}.$$

    对应的算子单调函数：$f(t) = (t-1)/\log t$（$t > 0$, $t \neq 1$；$f(1) = 1$）。

!!! example "例 46B.3 (经典均值的对应表)"
    | 均值名称 | 运算 $A \sigma B$ | 对应函数 $f(t)$ | 算子单调？ |
    |---------|-----------------|----------------|----------|
    | 算术均值 $\nabla$ | $(A + B)/2$ | $(1 + t)/2$ | 是 |
    | 调和均值 $!$ | $2(A^{-1} + B^{-1})^{-1}$ | $2t/(1 + t)$ | 是 |
    | 几何均值 $\#$ | $A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$ | $t^{1/2}$ | 是 |
    | 对数均值 | 见上 | $(t-1)/\log t$ | 是 |

    序关系：$f_!(t) = \frac{2t}{1+t} \leq f_{\#}(t) = \sqrt{t} \leq f_{\log}(t) = \frac{t-1}{\log t} \leq f_{\nabla}(t) = \frac{1+t}{2}$。

---

## 46B.5 几何均值的等价定义

<div class="context-flow" markdown>

**核心问题**：矩阵几何均值除了 $A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$ 之外，还有哪些等价的定义方式？

</div>

!!! theorem "定理 46B.2 (几何均值的等价定义)"
    对正定矩阵 $A, B \succ 0$，以下定义等价：

    1. **标准定义**：$A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$。

    2. **Riccati 方程**：$A \# B$ 是 $XA^{-1}X = B$ 的唯一正定解 $X$。

    3. **变分刻画**：$A \# B = \max\{X \succ 0 : \begin{pmatrix} A & X \\ X & B \end{pmatrix} \succeq 0\}$，最大值在 Löwner 偏序下取。

    4. **中点性质**：$A \# B$ 是正定矩阵流形上 $A$ 到 $B$ 的**测地线中点**（关于 Riemannian 度量 $ds^2 = \operatorname{tr}(X^{-1}dX X^{-1}dX)$）。

    5. **极分解**：$A \# B = |A^{1/2}B^{1/2}|$（矩阵 $A^{1/2}B^{1/2}$ 的极分解中的正部分），其中 $|M| = (M^*M)^{1/2}$。

??? proof "证明"
    **(1) $\Leftrightarrow$ (2)**：设 $X = A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$。则

    \begin{align}
    XA^{-1}X &= A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2} \underbrace{A^{1/2} \cdot A^{-1} \cdot A^{1/2}}_{= I} (A^{-1/2}BA^{-1/2})^{1/2}A^{1/2} \\
    &= A^{1/2}(A^{-1/2}BA^{-1/2})A^{1/2} = B.
    \end{align}

    唯一性：设 $X, Y \succ 0$ 都满足 $XA^{-1}X = YA^{-1}Y = B$。令 $U = A^{-1/2}X A^{-1/2}$，$V = A^{-1/2}Y A^{-1/2}$，则 $U^2 = V^2 = A^{-1/2}BA^{-1/2}$。由 $U, V \succ 0$ 和正定矩阵正平方根的唯一性，$U = V$，即 $X = Y$。

    **(1) $\Leftrightarrow$ (3)**：Schur 补条件：$\begin{pmatrix} A & X \\ X & B \end{pmatrix} \succeq 0$ 当且仅当 $B - XA^{-1}X \succeq 0$，即 $XA^{-1}X \preceq B$。

    最大化 $X$ 在 Löwner 偏序下，使得 $XA^{-1}X \preceq B$。最大值在 $XA^{-1}X = B$ 时取到，即 $X = A \# B$。

    **(1) $\Leftrightarrow$ (5)**：$A^{1/2}B^{1/2}$ 的极分解为 $A^{1/2}B^{1/2} = U|A^{1/2}B^{1/2}|$（$U$ 酉）。

    $|A^{1/2}B^{1/2}| = (B^{1/2}A B^{1/2})^{1/2} = B^{1/2}(B^{-1/2}AB^{-1/2})^{1/2}B^{1/2} = B \# A = A \# B$（利用对称性）。 $\blacksquare$

---

## 46B.6 几何均值的性质

<div class="context-flow" markdown>

**核心问题**：矩阵几何均值满足哪些代数和序性质？

</div>

!!! theorem "定理 46B.3 (矩阵几何均值的性质)"
    设 $A, B \succ 0$。矩阵几何均值 $A \# B$ 满足以下性质：

    1. **对称性**：$A \# B = B \# A$。
    2. **自反性**：$A \# A = A$。
    3. **逆的兼容**：$(A \# B)^{-1} = A^{-1} \# B^{-1}$。
    4. **合同不变**：对可逆 $C$，$C^*(A \# B)C = (C^*AC) \# (C^*BC)$。
    5. **行列式公式**：$\det(A \# B) = [\det(A) \cdot \det(B)]^{1/2}$。
    6. **单调性**：$A \preceq A'$，$B \preceq B'$ $\Rightarrow$ $A \# B \preceq A' \# B'$。
    7. **联合凹性**：映射 $(A, B) \mapsto A \# B$ 在 $\mathcal{P}_n \times \mathcal{P}_n$ 上是联合凹的：
    $$\lambda(A_1 \# B_1) + (1-\lambda)(A_2 \# B_2) \preceq (\lambda A_1 + (1-\lambda)A_2) \# (\lambda B_1 + (1-\lambda)B_2).$$
    8. **自对偶**：$(A \# B)^{-1} = A^{-1} \# B^{-1}$（与性质 3 相同，称为自对偶性）。
    9. **算术-几何均值不等式**：$A \# B \preceq A \nabla B = (A + B)/2$。
    10. **幂等性**：$A \# B = A$ 当且仅当 $A = B$。
    11. **迹不等式**：$\operatorname{tr}(A \# B) \leq [\operatorname{tr}(A) \cdot \operatorname{tr}(B)]^{1/2}$。
    12. **与标量运算的相容**：若 $A$ 和 $B$ 可交换（$AB = BA$），则 $A \# B = (AB)^{1/2} = A^{1/2}B^{1/2}$。

??? proof "证明"
    **(1) 对称性**：设 $X = A \# B$。由 Riccati 方程，$XA^{-1}X = B$。则
    $$XB^{-1}X = X(XA^{-1}X)^{-1}X = X \cdot X^{-1}AX^{-1} \cdot X = A.$$
    因此 $X$ 也满足 $XB^{-1}X = A$，即 $X = B \# A$。

    **(3) 逆的兼容**：
    $(A \# B)^{-1} = [A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}]^{-1} = A^{-1/2}(A^{-1/2}BA^{-1/2})^{-1/2}A^{-1/2}$。

    注意 $(A^{-1/2}BA^{-1/2})^{-1/2} = (A^{1/2}B^{-1}A^{1/2})^{1/2}$，因此
    $(A \# B)^{-1} = A^{-1/2}(A^{1/2}B^{-1}A^{1/2})^{1/2}A^{-1/2}$。

    而 $A^{-1} \# B^{-1} = A^{-1/2}((A^{-1})^{-1/2}B^{-1}(A^{-1})^{-1/2})^{1/2}A^{-1/2} = A^{-1/2}(A^{1/2}B^{-1}A^{1/2})^{1/2}A^{-1/2}$。

    因此 $(A \# B)^{-1} = A^{-1} \# B^{-1}$。

    **(4) 合同不变**：设 $Y = C^*(A \# B)C$。验证 Riccati 方程：
    $$Y(C^*AC)^{-1}Y = C^*(A\#B)C \cdot C^{-1}A^{-1}(C^*)^{-1} \cdot C^*(A\#B)C = C^*(A\#B)A^{-1}(A\#B)C = C^*BC.$$
    因此 $Y = (C^*AC) \# (C^*BC)$。

    **(5) 行列式**：
    \begin{align}
    \det(A \# B) &= \det(A^{1/2}) \cdot \det((A^{-1/2}BA^{-1/2})^{1/2}) \cdot \det(A^{1/2}) \\
    &= \det(A) \cdot [\det(A^{-1/2}BA^{-1/2})]^{1/2} \\
    &= \det(A) \cdot [\det(A)^{-1} \det(B)]^{1/2} \\
    &= [\det(A)]^{1/2} [\det(B)]^{1/2}.
    \end{align}

    **(9) AM-GM**：设 $T = A^{-1/2}BA^{-1/2} \succ 0$。$A \# B \preceq A \nabla B$ 等价于 $T^{1/2} \preceq (I + T)/2$，即 $(I - T^{1/2})^2 \succeq 0$。这显然成立。 $\blacksquare$

---

## 46B.7 AM-GM-HM 不等式

!!! theorem "定理 46B.4 (矩阵算术-几何-调和均值不等式)"
    对正定矩阵 $A, B \succ 0$：

    $$A \,!\, B \preceq A \# B \preceq A \nabla B,$$

    即

    $$2(A^{-1} + B^{-1})^{-1} \preceq A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2} \preceq \frac{A + B}{2}.$$

    这是标量不等式 $\frac{2ab}{a+b} \leq \sqrt{ab} \leq \frac{a+b}{2}$（$a, b > 0$）的矩阵推广。

    更一般地，对带权重 $t \in [0, 1]$：
    $$A \,!_t\, B \preceq A \#_t B \preceq A \nabla_t B.$$

??? proof "证明"
    **右侧**（$A \# B \preceq A \nabla B$）：已在定理 46B.3 性质 (9) 中证明。

    **左侧**（$A \,!\, B \preceq A \# B$）：由逆的兼容性（定理 46B.3 性质 (3)）和右侧不等式：

    $(A \# B)^{-1} = A^{-1} \# B^{-1} \preceq A^{-1} \nabla B^{-1} = \frac{A^{-1} + B^{-1}}{2} = (A \,!\, B)^{-1}.$

    由逆的反单调性（$X^{-1} \preceq Y^{-1} \Leftrightarrow Y \preceq X$，对 $X, Y \succ 0$），

    $A \,!\, B \preceq A \# B.$ $\blacksquare$

!!! example "例 46B.4 (AM-GM-HM 的数值验证)"
    取 $A = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix}$（对角矩阵，因此可交换）。

    - 算术均值：$A \nabla B = \begin{pmatrix} 5/2 & 0 \\ 0 & 5/2 \end{pmatrix}$。
    - 几何均值：$A \# B = (AB)^{1/2} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}$。
    - 调和均值：$A \,!\, B = 2(A^{-1} + B^{-1})^{-1} = \begin{pmatrix} 8/5 & 0 \\ 0 & 8/5 \end{pmatrix}$。

    验证：$8/5 < 2 < 5/2$，即 $A \,!\, B \prec A \# B \prec A \nabla B$。

---

## 46B.8 矩阵幂均值

<div class="context-flow" markdown>

**核心问题**：如何定义在算术均值和几何均值之间插值的矩阵均值？

</div>

!!! definition "定义 46B.7 (矩阵幂均值)"
    对正定矩阵 $A_1, \ldots, A_k \succ 0$，权重 $w = (w_1, \ldots, w_k)$（$w_i > 0$，$\sum w_i = 1$），和参数 $s \in [-1, 1]$，**矩阵幂均值**（matrix power mean）$P_s(w; A_1, \ldots, A_k)$ 定义为方程

    $$X = \sum_{i=1}^k w_i (X \#_s A_i)$$

    的唯一正定解 $X$，其中 $X \#_s A_i = X^{1/2}(X^{-1/2}A_i X^{-1/2})^s X^{1/2}$ 是带参数 $s$ 的加权几何均值。

!!! theorem "定理 46B.5 (Lim-Palfia 定理)"
    Lim 和 Palfia (2012) 证明了以下结果：

    1. 矩阵幂均值方程的正定解 $P_s$ **存在且唯一**。
    2. $P_s$ 关于参数 $s$ **单调递增**：$s_1 \leq s_2 \Rightarrow P_{s_1} \preceq P_{s_2}$。
    3. **极限行为**：
        - $P_1 = \sum w_i A_i$（算术均值）。
        - $P_{-1} = (\sum w_i A_i^{-1})^{-1}$（加权调和均值）。
        - $\lim_{s \to 0} P_s = G(w; A_1, \ldots, A_k)$（加权 Karcher 均值/几何均值）。
    4. $P_s$ 满足 Kubo-Ando 型单调性和合同不变性。

    因此矩阵幂均值提供了从调和均值（$s = -1$）经过几何均值（$s = 0$）到算术均值（$s = 1$）的**连续插值**。

---

## 46B.9 Ando-Li-Mathias 多变量几何均值

<div class="context-flow" markdown>

**核心问题**：如何将两变量几何均值推广到多个正定矩阵？

</div>

!!! theorem "定理 46B.6 (Ando-Li-Mathias 均值)"
    对 $k$ 个正定矩阵 $A_1, \ldots, A_k \succ 0$，**Ando-Li-Mathias (ALM) 多变量几何均值** $G(A_1, \ldots, A_k)$ 通过以下迭代定义：

    设 $A_i^{(0)} = A_i$。递归定义
    $$A_i^{(r+1)} = G_2\left(A_1^{(r)}, \ldots, \widehat{A_i^{(r)}}, \ldots, A_k^{(r)}\right) \quad (i = 1, \ldots, k),$$

    其中 $\widehat{A_i^{(r)}}$ 表示去掉第 $i$ 个，$G_2$ 是 $(k-1)$ 变量的几何均值（递归定义）。

    对 $k = 2$ 的基本情形：$G_2(A, B) = A \# B$。

    **收敛性**：$A_i^{(r)} \to G$（$r \to \infty$）对所有 $i$ 收敛到同一极限 $G = G(A_1, \ldots, A_k)$。

    ALM 均值满足以下性质：
    1. **对称性**：$G(A_{\pi(1)}, \ldots, A_{\pi(k)}) = G(A_1, \ldots, A_k)$ 对任意置换 $\pi$。
    2. **合同不变**：$G(C^*A_1 C, \ldots, C^*A_k C) = C^* G(A_1, \ldots, A_k) C$。
    3. **单调性**：$A_i \preceq B_i$ $\Rightarrow$ $G(A_1, \ldots, A_k) \preceq G(B_1, \ldots, B_k)$。
    4. **行列式**：$\det G(A_1, \ldots, A_k) = [\det(A_1) \cdots \det(A_k)]^{1/k}$。
    5. **自对偶**：$G(A_1^{-1}, \ldots, A_k^{-1})^{-1} = G(A_1, \ldots, A_k)$。

---

## 46B.10 Karcher 均值

!!! definition "定义 46B.8 (Karcher 均值)"
    对正定矩阵 $A_1, \ldots, A_k \succ 0$ 和权重 $w_1, \ldots, w_k > 0$（$\sum w_i = 1$），**Karcher 均值**（也称 Riemannian 重心、Frechet 均值）定义为

    $$G = \arg\min_{X \succ 0} \sum_{i=1}^{k} w_i \, d^2(X, A_i),$$

    其中 $d(X, Y) = \|\log(X^{-1/2}YX^{-1/2})\|_F$ 是正定矩阵流形上的 **Riemannian 距离**（也称 **Thompson 度量**的 Frobenius 版本）。

!!! theorem "定理 46B.7 (Karcher 均值的存在唯一性)"
    Karcher 均值存在且唯一。这一结果是正定矩阵流形 $\mathcal{P}_n = \{X \in \mathbb{C}^{n \times n} : X \succ 0\}$ 上 **Cartan-Hadamard 定理**的推论。

    Karcher 均值等价地由以下**矩阵方程**刻画：
    $$\sum_{i=1}^{k} w_i \log(G^{-1/2} A_i G^{-1/2}) = 0.$$

    这个方程的物理意义是：$G$ 是使得到各 $A_i$ 的"对数偏移"的加权和为零的点——即"重心"条件。

??? proof "证明"
    **存在性**：目标函数 $\phi(X) = \sum w_i d^2(X, A_i)$ 在 $\mathcal{P}_n$ 上连续，且当 $X$ 趋于 $\partial \mathcal{P}_n$（奇异矩阵的边界）或 $\|X\| \to \infty$ 时 $\phi(X) \to \infty$（因为 $d(X, A_i) \to \infty$）。因此 $\phi$ 在 $\mathcal{P}_n$ 上达到最小值。

    **唯一性**：$\mathcal{P}_n$ 配备 Riemannian 度量 $g_X(H, K) = \operatorname{tr}(X^{-1}HX^{-1}K)$ 是一个**Cartan-Hadamard 流形**（完备、单连通、截面曲率 $\leq 0$）。在 Cartan-Hadamard 流形上，距离的平方 $d^2(\cdot, p)$ 是严格凸函数。因此 $\phi = \sum w_i d^2(\cdot, A_i)$ 是严格凸的，最小值唯一。

    **矩阵方程**：$\phi(X)$ 的 Riemannian 梯度为
    $$\operatorname{grad} \phi(X) = -2\sum_{i=1}^k w_i \, X \log(X^{-1/2}A_i X^{-1/2}).$$
    令梯度为零（$X = G$），得到 $\sum w_i \log(G^{-1/2}A_i G^{-1/2}) = 0$。 $\blacksquare$

!!! example "例 46B.5 (Karcher 均值的计算)"
    对 $k = 2$（等权重），Karcher 均值方程为 $\log(G^{-1/2}AG^{-1/2}) + \log(G^{-1/2}BG^{-1/2}) = 0$，即 $G^{-1/2}AG^{-1/2} = (G^{-1/2}BG^{-1/2})^{-1}$，即 $G^{-1/2}AG^{-1/2} \cdot G^{-1/2}BG^{-1/2} = I$，即 $G^{-1/2}ABG^{-1/2} = I$... 这不太对（因为 $AB$ 不一定 Hermite）。

    正确的推导：$\log(G^{-1/2}AG^{-1/2}) = -\log(G^{-1/2}BG^{-1/2})$，即 $G^{-1/2}AG^{-1/2} = (G^{-1/2}BG^{-1/2})^{-1} = G^{1/2}B^{-1}G^{1/2}$。因此 $A = G B^{-1} G$，即 $GA^{-1}G = B$，即 $G = A \# B$。

    这验证了两变量 Karcher 均值就是几何均值 $A \# B$。

---

## 46B.11 正定流形的几何

<div class="context-flow" markdown>

**核心问题**：正定矩阵集合的自然几何结构是什么？

</div>

!!! definition "定义 46B.9 (正定矩阵流形的 Riemannian 结构)"
    正定矩阵集合 $\mathcal{P}_n = \{X \in \mathbb{C}^{n \times n} : X = X^*, X \succ 0\}$ 是一个光滑流形（作为 Hermite 矩阵空间的开子集）。在 $X \in \mathcal{P}_n$ 处的切空间是全体 Hermite 矩阵 $T_X\mathcal{P}_n = \mathcal{H}_n$。

    **Riemannian 度量**：在 $X$ 处的内积定义为
    $$\langle H, K \rangle_X = \operatorname{tr}(X^{-1}HX^{-1}K),$$
    其中 $H, K \in T_X\mathcal{P}_n = \mathcal{H}_n$。

!!! definition "定义 46B.10 (Thompson 度量)"
    正定矩阵上的 **Thompson 度量**定义为
    $$d_T(A, B) = \|\log(A^{-1/2}BA^{-1/2})\|_{\text{op}} = \max\{|\log \lambda_i(A^{-1/2}BA^{-1/2})|\},$$
    其中 $\lambda_i$ 表示特征值，$\|\cdot\|_{\text{op}}$ 是算子范数。

    Riemannian 距离（Frobenius 版本）为
    $$d_R(A, B) = \|\log(A^{-1/2}BA^{-1/2})\|_F = \left(\sum_i [\log \lambda_i(A^{-1/2}BA^{-1/2})]^2\right)^{1/2}.$$

    两种距离都满足：

    1. **合同不变**：$d(C^*AC, C^*BC) = d(A, B)$（可逆 $C$）。
    2. **逆不变**：$d(A^{-1}, B^{-1}) = d(A, B)$。

!!! theorem "定理 46B.8 (Cartan-Hadamard 性质)"
    $(\mathcal{P}_n, \langle \cdot, \cdot \rangle)$ 是一个**完备、单连通、非正曲率**的 Riemannian 流形（Cartan-Hadamard 流形）。具体来说：

    1. **完备性**：每条测地线都可以无限延伸。
    2. **非正截面曲率**：对所有截面，$K \leq 0$。
    3. **测地线**：从 $A$ 到 $B$ 的测地线为
    $$\gamma(t) = A^{1/2}(A^{-1/2}BA^{-1/2})^t A^{1/2} = A \#_t B, \quad t \in [0, 1].$$
    4. **指数映射**：$\exp_A(H) = A^{1/2} \exp(A^{-1/2}HA^{-1/2}) A^{1/2}$。
    5. **Cartan-Hadamard 定理的推论**：任意两点之间存在唯一的最短测地线；$\exp_A$ 是微分同胚。

??? proof "证明"
    **测地线公式**：曲线 $\gamma(t) = A \#_t B = A^{1/2}(A^{-1/2}BA^{-1/2})^t A^{1/2}$ 是测地线，因为它满足测地线方程
    $$\nabla_{\dot\gamma}\dot\gamma = 0,$$
    其中 $\nabla$ 是 Riemannian 联络。

    验证：在 $A = I$ 处，测地线为 $\gamma(t) = B^t = \exp(t \log B)$，$\dot\gamma(t) = B^t \log B$。Riemannian 联络在 $I$ 处简化为普通导数，$\nabla_{\dot\gamma}\dot\gamma = 0$ 可以直接验证。一般情形由合同不变性 $\gamma(t) \mapsto C^*\gamma(t)C$ 得到。

    **非正曲率**：截面曲率的计算利用 $GL_n$ 的李代数结构。正定矩阵流形 $\mathcal{P}_n$ 同构于对称空间 $GL_n(\mathbb{C})/U_n$，其截面曲率可以精确计算：
    $$K(H, K) = -\frac{3}{4}\frac{\|[X^{-1}H, X^{-1}K]\|_X^2}{\|H\|_X^2 \|K\|_X^2 - \langle H, K\rangle_X^2} \leq 0,$$
    其中 $[\cdot, \cdot]$ 是矩阵交换子。 $\blacksquare$

!!! example "例 46B.6 (测地线的几何意义)"
    从 $A$ 到 $B$ 的测地线 $\gamma(t) = A \#_t B$ 有以下性质：

    - $\gamma(0) = A$，$\gamma(1) = B$，$\gamma(1/2) = A \# B$（几何均值 = 测地线中点）。
    - $d_R(A, \gamma(t)) = t \cdot d_R(A, B)$（匀速参数化）。
    - 在对数坐标 $\log A^{-1/2}\gamma(t)A^{-1/2} = t \log(A^{-1/2}BA^{-1/2})$ 下，测地线是直线。

    这解释了为什么几何均值 $A \# B$ 是正定矩阵之间最自然的"中点"——它是曲面上最短路径的中点。

---

## 46B.12 正半定矩阵的推广

<div class="context-flow" markdown>

**核心问题**：几何均值能否推广到正半定（奇异）矩阵？

</div>

!!! theorem "定理 46B.9 (正半定矩阵的几何均值)"
    矩阵几何均值可以通过极限从正定矩阵推广到正半定矩阵：

    对 $A, B \succeq 0$（半正定，可能奇异），定义
    $$A \# B = \lim_{\epsilon \to 0^+} (A + \epsilon I) \# (B + \epsilon I).$$

    该极限存在且满足：

    1. $A \# B \succeq 0$（半正定）。
    2. 当 $A, B \succ 0$ 时，与原定义一致。
    3. **连续性**：$A_n \to A$，$B_n \to B$ $\Rightarrow$ $A_n \# B_n \to A \# B$。
    4. **行列式**：$\det(A \# B) = [\det A \cdot \det B]^{1/2}$（特别地，若 $A$ 或 $B$ 奇异，则 $A \# B$ 奇异）。
    5. **秩**：$\operatorname{rank}(A \# B) \leq \min(\operatorname{rank} A, \operatorname{rank} B)$。

??? proof "证明"
    **极限存在**：$A_\epsilon = A + \epsilon I \succ 0$，$B_\epsilon = B + \epsilon I \succ 0$。$A_\epsilon \# B_\epsilon$ 关于 $\epsilon$ 单调递减（因为 $A_\epsilon \succeq A_{\epsilon'}$ 当 $\epsilon \leq \epsilon'$，由单调性）。又 $A_\epsilon \# B_\epsilon \succeq 0$（有下界），因此极限存在。

    **秩性质**：$A \# B = A^{1/2}(A^{-1/2}BA^{-1/2})^{1/2}A^{1/2}$（通过广义逆的推广）。若 $\operatorname{rank}(A) = r < n$，则 $A^{1/2}$ 的值域是 $r$ 维的，$A \# B$ 的值域包含在其中，因此 $\operatorname{rank}(A \# B) \leq r$。 $\blacksquare$

!!! example "例 46B.7 (正半定情形的计算)"
    取 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$。

    $A_\epsilon = \begin{pmatrix} 1 + \epsilon & 0 \\ 0 & \epsilon \end{pmatrix}$，$B_\epsilon = \begin{pmatrix} \epsilon & 0 \\ 0 & 1 + \epsilon \end{pmatrix}$。

    由于对角矩阵可交换，$A_\epsilon \# B_\epsilon = (A_\epsilon B_\epsilon)^{1/2} = \begin{pmatrix} \sqrt{\epsilon(1+\epsilon)} & 0 \\ 0 & \sqrt{\epsilon(1+\epsilon)} \end{pmatrix}$。

    $\lim_{\epsilon \to 0} A_\epsilon \# B_\epsilon = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} = 0$。

    因此 $A \# B = 0$——当 $A$ 和 $B$ 的值域正交时，几何均值为零。

---

## 46B.13 习题

!!! example "例 46B.8"
    验证 Kubo-Ando 定理：对算术均值 $A \nabla B = (A + B)/2$，对应的算子单调函数 $f(t) = 1 \nabla t = (1 + t)/2$ 确实是算子单调的。

!!! example "例 46B.9"
    证明几何均值的对称性 $A \# B = B \# A$ 的另一种方法：利用等价定义 $A \# B = |A^{1/2}B^{1/2}|$（极分解的模），以及 $|M| = |M^*|$ 的性质。

!!! example "例 46B.10"
    对可交换的正定矩阵 $A, B$（$AB = BA$），证明 $A \# B = (AB)^{1/2} = A^{1/2}B^{1/2}$。

!!! example "例 46B.11"
    计算以下矩阵对的几何均值：
    $$A = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 0 \\ 0 & 9 \end{pmatrix}.$$

!!! example "例 46B.12"
    证明矩阵 AM-GM 不等式 $A \# B \preceq A \nabla B$ 的迹版本：
    $$\operatorname{tr}(A \# B) \leq \frac{\operatorname{tr}(A) + \operatorname{tr}(B)}{2}.$$

!!! example "例 46B.13"
    设 $A, B, C \succ 0$。证明几何均值的**中值性质**：
    $$\lambda_{\min}(A^{-1}B) \cdot A \preceq A \# B \preceq \lambda_{\max}(A^{-1}B) \cdot A,$$
    其中 $\lambda_{\min}, \lambda_{\max}$ 表示最小和最大特征值。

!!! example "例 46B.14"
    对三个正定矩阵 $A, B, C$，利用 Ando-Li-Mathias 迭代的第一步计算：
    $$A^{(1)} = B \# C, \quad B^{(1)} = A \# C, \quad C^{(1)} = A \# B.$$
    验证 $A^{(1)}, B^{(1)}, C^{(1)}$ 比 $A, B, C$ "更接近"。

!!! example "例 46B.15"
    证明 Thompson 度量满足合同不变性：$d_T(C^*AC, C^*BC) = d_T(A, B)$ 对可逆 $C$。

    **提示**：利用 $\log$ 的谱映射性质和合同变换下特征值的不变性。

!!! example "例 46B.16"
    设 $\gamma(t) = A \#_t B$（$t \in [0, 1]$）是正定矩阵流形上的测地线。证明
    $$d_R(A, \gamma(t)) = t \cdot d_R(A, B),$$
    即测地线是匀速参数化的。

!!! example "例 46B.17"
    （Karcher 均值的计算）对三个 $2 \times 2$ 正定矩阵 $A_1 = I$，$A_2 = 2I$，$A_3 = 4I$（等权重 $w_i = 1/3$），计算 Karcher 均值。

    **提示**：由于矩阵都是标量矩阵 $\alpha I$，Karcher 均值也是标量矩阵，等于标量几何均值的矩阵版本：$(1 \cdot 2 \cdot 4)^{1/3} I = 2I$。

!!! example "例 46B.18"
    证明正定矩阵流形 $\mathcal{P}_n$ 上的截面曲率满足 $K \leq 0$。

    **提示**：利用截面曲率公式 $K(H, K) = -\frac{3}{4}\|[X^{-1}H, X^{-1}K]\|^2 / (\ldots)$，其中交换子的范数总是非负的。

!!! example "例 46B.19"
    设 $A \succeq 0$ 是奇异正半定矩阵（$\det A = 0$），$B \succ 0$。计算 $A \# B$ 并证明 $A \# B$ 也是奇异的。

!!! example "例 46B.20"
    （矩阵均值在扩散张量成像中的应用）在 DTI 中，扩散张量 $D_1 = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$，$D_2 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    计算 $D_1 \nabla D_2$，$D_1 \# D_2$，$D_1 \,!\, D_2$，并比较三者的行列式。验证只有几何均值保持 $\det(D_1 \# D_2) = [\det D_1 \cdot \det D_2]^{1/2}$。

---

## 练习题

1. **[基础] 计算 $A = \operatorname{diag}(4, 1)$ 与 $B = \operatorname{diag}(1, 9)$ 的矩阵几何均值 $A \# B$。**
   ??? success "参考答案"
       由于 $A, B$ 都是对角矩阵，它们是可交换的。因此 $A \# B = (AB)^{1/2} = \operatorname{diag}(\sqrt{4 \cdot 1}, \sqrt{1 \cdot 9}) = \operatorname{diag}(2, 3)$。

2. **[Riccati] 验证 $X = A \# B$ 满足 Riccati 方程 $X A^{-1} X = B$。**
   ??? success "参考答案"
       代入定义：$A^{1/2}(A^{-1/2} B A^{-1/2})^{1/2} A^{1/2} \cdot A^{-1} \cdot A^{1/2} (A^{-1/2} B A^{-1/2})^{1/2} A^{1/2} = A^{1/2} (A^{-1/2} B A^{-1/2}) A^{1/2} = B$。

3. **[AM-GM-HM] 在 Löwner 偏序下，将三种经典矩阵均值按从小到大排序。**
   ??? success "参考答案"
       调和均值 $\preceq$ 几何均值 $\preceq$ 算术均值。
       即 $2(A^{-1} + B^{-1})^{-1} \preceq A \# B \preceq \frac{A+B}{2}$。

4. **[行列式] 证明 $\det(A \# B) = \sqrt{\det A \cdot \det B}$。**
   ??? success "参考答案"
       $\det(A \# B) = \det(A^{1/2}) \det((A^{-1/2} B A^{-1/2})^{1/2}) \det(A^{1/2}) = \det(A) \sqrt{\det(A^{-1}B)} = \det(A) \sqrt{\det(A)^{-1} \det B} = \sqrt{\det A \det B}$。

5. **[变换性] 证明：若 $C$ 可逆，则 $C^*(A \# B)C = (C^*AC) \# (C^*BC)$。**
   ??? success "参考答案"
       这是 Kubo-Ando 公理中的变换等式。设 $X = A \# B$，则有 $X A^{-1} X = B$。两边左乘 $C^*$ 右乘 $C$，并在中间插入 $C^{-1}(C^*)^{-1}$ 得：$(C^* X C) (C^* A C)^{-1} (C^* X C) = C^* B C$。由此可见 $C^* X C$ 是 $C^* A C$ 与 $C^* B C$ 的几何均值。

6. **[多变量] 定义 $k$ 个矩阵 $A_1, \dots, A_k$ 的 Karcher 均值。**
   ??? success "参考答案"
       Karcher 均值是正定流形上到各点 Riemannian 距离平方和最小的点：$G = \arg\min_X \sum w_i \|\log(X^{-1/2} A_i X^{-1/2})\|_F^2$。

7. **[度量] 正定矩阵流形上的 Riemannian 度量是如何定义的？**
   ??? success "参考答案"
       在点 $X$ 处，两个切向量（Hermite 矩阵）$H, K$ 的内积定义为 $\langle H, K \rangle_X = \operatorname{tr}(X^{-1} H X^{-1} K)$。

8. **[测地线] 连接 $A$ 到 $B$ 的测地线路径是什么？**
   ??? success "参考答案"
       测地线是加权几何均值路径：$\gamma(t) = A \#_t B = A^{1/2}(A^{-1/2} B A^{-1/2})^t A^{1/2}$。

9. **[交换性] 若 $A$ 与 $B$ 交换，几何均值简化为什么形式？**
   ??? success "参考答案"
       $A \# B = A^{1/2} B^{1/2} = (AB)^{1/2}$。

10. **[奇异情形] 计算 $A = \operatorname{diag}(1, 0)$ 与 $B = \operatorname{diag}(0, 1)$ 的几何均值。**
    ??? success "参考答案"
        使用极限 $\epsilon \to 0$：$(A+\epsilon I) \# (B+\epsilon I) = \operatorname{diag}(\sqrt{\epsilon(1+\epsilon)}, \sqrt{\epsilon(1+\epsilon)})$。当 $\epsilon \to 0$ 时，均值趋于零矩阵。

## 本章小结

本章探讨了代数均值与微分几何的交叉融合：

1. **公理化基础**：确立了 Kubo-Ando 框架作为定义矩阵均值的严谨路径。
2. **谱分析关联**：将矩阵几何均值与 Riccati 方程及算子单调函数理论联系起来。
3. **几何合成**：发展了正定矩阵流形的 Riemannian 几何，识别出几何均值即为测地线的中点。
4. **多变量推广**：通过 Karcher 变分法将二变量均值推广至 $k$ 变量，为现代信号处理提供了数学工具。

