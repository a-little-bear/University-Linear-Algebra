# 第 35 章 Hadamard 积

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 正定矩阵(Ch16) · Kronecker积(Ch19)

**本章脉络**：Hadamard 积定义 → Schur 积定理 → Oppenheim 不等式 → Hadamard 不等式 → 与 Kronecker 积的关系 → 正映射 → 应用

**延伸**：Hadamard 积在统计学（协方差锥化/tapering）、信号处理（逐元素操作）和量子信息（Schur 通道即 Hadamard 积通道）中有重要应用

</div>

Hadamard 积（又称 Schur 积、逐元素积）是矩阵的另一种乘法运算：对应位置的元素相乘。这一简单的运算具有令人惊讶的深刻性质。Schur（1911）证明了两个半正定矩阵的 Hadamard 积仍然半正定——这一优雅的结果成为正定矩阵理论的基石。Oppenheim（1930）和 Hadamard（1893）建立了 Hadamard 积与行列式之间的不等式，这些不等式在现代统计学和信号处理中有广泛应用。

本章系统发展 Hadamard 积理论，从基本定义和代数性质出发，经由 Schur 积定理和 Oppenheim 不等式，到谱性质和正映射理论，最终展示其在协方差锥化、神经网络和图论中的应用。

---

## 35.1 Hadamard 积的定义与基本性质

<div class="context-flow" markdown>

**核心问题**：逐元素乘法作为矩阵运算有什么代数结构？

</div>

### 定义

!!! definition "定义 35.1 (Hadamard 积)"
    设 $A = (a_{ij})$ 和 $B = (b_{ij})$ 是同阶矩阵（$m \times n$）。$A$ 和 $B$ 的 **Hadamard 积**（也称 Schur 积、逐元素积，记作 $A \circ B$）定义为

    $$
    (A \circ B)_{ij} = a_{ij} b_{ij}
    $$

    即对应元素相乘。

!!! note "注"
    Hadamard 积的记号不统一：$A \circ B$（本书采用）、$A \odot B$、$A * B$、$A \bullet B$ 都有人使用。在 MATLAB/NumPy 中用 `.*` 表示。

### 基本代数性质

!!! theorem "定理 35.1 (Hadamard 积的代数性质)"
    设 $A, B, C$ 是同阶矩阵，$\alpha \in \mathbb{C}$。则：

    1. **交换律**：$A \circ B = B \circ A$
    2. **结合律**：$(A \circ B) \circ C = A \circ (B \circ C)$
    3. **分配律**：$A \circ (B + C) = A \circ B + A \circ C$
    4. **数乘**：$\alpha(A \circ B) = (\alpha A) \circ B = A \circ (\alpha B)$
    5. **单位元**：$J \circ A = A$，其中 $J$ 是全 1 矩阵
    6. **零化**：$0 \circ A = 0$
    7. **转置**：$(A \circ B)^T = A^T \circ B^T$
    8. **共轭转置**：$(A \circ B)^* = A^* \circ B^*$

!!! note "注"
    与通常的矩阵乘法不同，Hadamard 积是交换的。它使得同阶矩阵的集合 $\mathbb{C}^{m \times n}$ 成为一个交换代数（以 Hadamard 积为乘法）。但 Hadamard 积不与通常矩阵乘法兼容——$A \circ (BC)$ 一般没有简单的表达式。

### 与 Kronecker 积的关系

!!! theorem "定理 35.2 (Hadamard 积与 Kronecker 积)"
    设 $A, B \in \mathbb{C}^{n \times n}$。则

    $$
    A \circ B = P^T (A \otimes B) P
    $$

    其中 $P \in \mathbb{R}^{n^2 \times n}$ 是选择矩阵，定义为

    $$
    P = \sum_{i=1}^n (\boldsymbol{e}_i \otimes \boldsymbol{e}_i) \boldsymbol{e}_i^T
    $$

    即 $P$ 的第 $i$ 列是 $\boldsymbol{e}_i \otimes \boldsymbol{e}_i$（$n^2$ 维向量，在第 $(i-1)n + i$ 个位置为 1，其余为 0）。

    等价地，$A \circ B$ 是 $A \otimes B$ 的"对角块提取"：

    $$
    (A \circ B)_{ij} = (\boldsymbol{e}_i \otimes \boldsymbol{e}_i)^T (A \otimes B) (\boldsymbol{e}_j \otimes \boldsymbol{e}_j) = a_{ij} b_{ij}
    $$

??? proof "证明"
    $(P^T(A \otimes B)P)_{ij} = (\boldsymbol{e}_i \otimes \boldsymbol{e}_i)^T (A \otimes B) (\boldsymbol{e}_j \otimes \boldsymbol{e}_j)$。

    由 Kronecker 积的混合积性质：

    $$
    (A \otimes B)(\boldsymbol{e}_j \otimes \boldsymbol{e}_j) = (A\boldsymbol{e}_j) \otimes (B\boldsymbol{e}_j) = \boldsymbol{a}_j \otimes \boldsymbol{b}_j
    $$

    其中 $\boldsymbol{a}_j$ 和 $\boldsymbol{b}_j$ 分别是 $A$ 和 $B$ 的第 $j$ 列。

    $$
    (\boldsymbol{e}_i \otimes \boldsymbol{e}_i)^T (\boldsymbol{a}_j \otimes \boldsymbol{b}_j) = (\boldsymbol{e}_i^T \boldsymbol{a}_j)(\boldsymbol{e}_i^T \boldsymbol{b}_j) = a_{ij} b_{ij}
    $$

    因此 $(P^T(A \otimes B)P)_{ij} = a_{ij}b_{ij} = (A \circ B)_{ij}$。$\blacksquare$

!!! note "注"
    这个关系的重要性在于：它将 Hadamard 积（看似"不规范"的运算）化为 Kronecker 积（良好理解的运算）加上投影/提取。许多 Hadamard 积的性质可以通过这个关系从 Kronecker 积的性质推导出来。

!!! example "例 35.1"
    设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 5 & 6 \\ 7 & 8 \end{pmatrix}$。

    $A \circ B = \begin{pmatrix} 5 & 12 \\ 21 & 32 \end{pmatrix}$。

    验证 Kronecker 关系：$P = \begin{pmatrix} 1 & 0 \\ 0 & 0 \\ 0 & 0 \\ 0 & 1 \end{pmatrix}$（$4 \times 2$ 矩阵）。

    $A \otimes B = \begin{pmatrix} 5 & 6 & 10 & 12 \\ 7 & 8 & 14 & 16 \\ 15 & 18 & 20 & 24 \\ 21 & 24 & 28 & 32 \end{pmatrix}$。

    $P^T(A \otimes B)P = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} 5 & 12 \\ 7 & 16 \\ 15 & 24 \\ 21 & 32 \end{pmatrix} = \begin{pmatrix} 5 & 12 \\ 21 & 32 \end{pmatrix} = A \circ B$。 ✓

### 秩一分解

!!! theorem "定理 35.3 (Hadamard 积的秩一表示)"
    设 $B = \sum_{k=1}^r \boldsymbol{u}_k \boldsymbol{v}_k^T$ 是 $B$ 的秩一分解（$r = \operatorname{rank}(B)$）。则

    $$
    A \circ B = \sum_{k=1}^r \operatorname{diag}(\boldsymbol{u}_k) \cdot A \cdot \operatorname{diag}(\boldsymbol{v}_k)
    $$

    其中 $\operatorname{diag}(\boldsymbol{u})$ 是以 $\boldsymbol{u}$ 为对角元素的对角矩阵。

??? proof "证明"
    $(A \circ \boldsymbol{u}\boldsymbol{v}^T)_{ij} = a_{ij} u_i v_j = (\operatorname{diag}(\boldsymbol{u}) A \operatorname{diag}(\boldsymbol{v}))_{ij}$。

    由 Hadamard 积对加法的分配律，对一般的 $B = \sum_k \boldsymbol{u}_k\boldsymbol{v}_k^T$：

    $$
    A \circ B = \sum_k A \circ (\boldsymbol{u}_k\boldsymbol{v}_k^T) = \sum_k \operatorname{diag}(\boldsymbol{u}_k) A \operatorname{diag}(\boldsymbol{v}_k)
    $$

    $\blacksquare$

---

## 35.2 Schur 积定理

<div class="context-flow" markdown>

**核心问题**：两个半正定矩阵的 Hadamard 积仍然半正定吗？

</div>

### 主定理

!!! theorem "定理 35.4 (Schur 积定理, Schur 1911)"
    设 $A, B \in \mathbb{C}^{n \times n}$ 都是半正定矩阵（$A \geq 0$，$B \geq 0$）。则它们的 Hadamard 积 $A \circ B \geq 0$。

    若进一步 $A > 0$ 且 $B > 0$（正定），则 $A \circ B > 0$。

??? proof "证明"
    **证法一**（通过 Kronecker 积）：

    由定理 35.2，$A \circ B = P^T(A \otimes B)P$。

    由 Kronecker 积的性质：$A \geq 0$ 且 $B \geq 0$ 蕴含 $A \otimes B \geq 0$（因为 $A \otimes B$ 的特征值是 $\lambda_i(A)\mu_j(B) \geq 0$）。

    因此对任何 $\boldsymbol{x} \in \mathbb{C}^n$：

    $$
    \boldsymbol{x}^*(A \circ B)\boldsymbol{x} = \boldsymbol{x}^* P^T(A \otimes B)P\boldsymbol{x} = (P\boldsymbol{x})^*(A \otimes B)(P\boldsymbol{x}) \geq 0
    $$

    因此 $A \circ B \geq 0$。

    若 $A > 0$ 且 $B > 0$，则 $A \otimes B > 0$，$(P\boldsymbol{x})^*(A \otimes B)(P\boldsymbol{x}) = 0$ 当且仅当 $P\boldsymbol{x} = \boldsymbol{0}$。但 $P$ 是单射（列线性无关），所以 $P\boldsymbol{x} = \boldsymbol{0}$ 仅当 $\boldsymbol{x} = \boldsymbol{0}$。因此 $A \circ B > 0$。

    **证法二**（直接构造）：

    由 $B \geq 0$，$B$ 有分解 $B = \sum_{k=1}^n \mu_k \boldsymbol{v}_k \boldsymbol{v}_k^*$（谱分解，$\mu_k \geq 0$）。

    $A \circ B = \sum_k \mu_k (A \circ \boldsymbol{v}_k\boldsymbol{v}_k^*)$。

    由定理 35.3，$A \circ \boldsymbol{v}_k\boldsymbol{v}_k^* = \operatorname{diag}(\boldsymbol{v}_k) A \operatorname{diag}(\bar{\boldsymbol{v}}_k)$。

    对任何 $\boldsymbol{x}$：

    $$
    \boldsymbol{x}^* (A \circ \boldsymbol{v}_k\boldsymbol{v}_k^*) \boldsymbol{x} = \boldsymbol{x}^* \operatorname{diag}(\boldsymbol{v}_k) A \operatorname{diag}(\bar{\boldsymbol{v}}_k) \boldsymbol{x} = (\operatorname{diag}(\bar{\boldsymbol{v}}_k)\boldsymbol{x})^* A (\operatorname{diag}(\bar{\boldsymbol{v}}_k)\boldsymbol{x}) \geq 0
    $$

    因此 $A \circ \boldsymbol{v}_k\boldsymbol{v}_k^* \geq 0$，$A \circ B = \sum_k \mu_k (A \circ \boldsymbol{v}_k\boldsymbol{v}_k^*) \geq 0$。$\blacksquare$

!!! note "注"
    Schur 积定理是正定矩阵理论中最常用的结果之一。它的一个重要推论是：半正定矩阵在 Hadamard 积下构成一个**锥**（cone）：

    - 若 $A, B \geq 0$，则 $A \circ B \geq 0$（封闭性）。
    - 若 $A \geq 0$ 且 $\alpha \geq 0$，则 $\alpha A \geq 0$（数乘封闭）。

    这个锥在凸优化（半定规划）中非常重要。

### 特征值界

!!! theorem "定理 35.5 (Hadamard 积的特征值界)"
    设 $A, B \geq 0$，特征值按降序排列。则

    $$
    \lambda_{\min}(A) \cdot \lambda_{\min}(B) \leq \lambda_{\min}(A \circ B)
    $$

    不一定成立！但以下不等式成立：

    $$
    \lambda_{\min}(A \circ B) \geq \min_i b_{ii} \cdot \lambda_{\min}(A)
    $$

    $$
    \lambda_{\max}(A \circ B) \leq \max_i a_{ii} \cdot \lambda_{\max}(B)
    $$

    更一般地：

    $$
    \lambda_{\min}(A) \min_i b_{ii} \leq \lambda_{\min}(A \circ B) \leq \lambda_{\max}(A \circ B) \leq \lambda_{\max}(A) \max_i b_{ii}
    $$

??? proof "证明"
    上界：$A \leq \lambda_{\max}(A) I$（谱序），因此

    $$
    A \circ B \leq \lambda_{\max}(A) I \circ B = \lambda_{\max}(A) \operatorname{diag}(B)
    $$

    不对，$I \circ B = \operatorname{diag}(b_{11}, \ldots, b_{nn})$。

    所以 $\lambda_{\max}(A \circ B) \leq \lambda_{\max}(\lambda_{\max}(A) \operatorname{diag}(b_{11}, \ldots, b_{nn})) = \lambda_{\max}(A) \max_i b_{ii}$。

    但这个界太粗糙。更精确的上界：

    由 $B \leq (\max_i b_{ii}) I$（对半正定 $B$ 不一定成立）不对。但 $\operatorname{diag}(B) = I \circ B$ 满足 $B \leq (\max_i b_{ii}) J$... 也不对。

    正确的推导：对 $A \geq 0$，$A \leq \lambda_{\max}(A) I$ 不成立（因为 $A$ 的对角元素可能大于 $\lambda_{\max}(A)$... 不对，$a_{ii} \leq \lambda_{\max}(A)$（因为 $a_{ii}$ 是 Rayleigh 商 $\boldsymbol{e}_i^* A \boldsymbol{e}_i \leq \lambda_{\max}(A)$）。

    嗯，$a_{ii} \leq \lambda_{\max}(A)$ 对半正定矩阵成立（因为特征值 majorize 对角元素）。

    因此 $\operatorname{diag}(A) \leq \lambda_{\max}(A) I$（分量逐个小于等于）。但这不直接给出 $A \circ B$ 的界。

    用另一种方式：$\boldsymbol{x}^*(A \circ B)\boldsymbol{x} = \sum_{i,j} a_{ij}b_{ij}\bar{x}_i x_j$。由 $B \geq 0$，设 $B = \sum_k \mu_k \boldsymbol{v}_k\boldsymbol{v}_k^*$：

    $$
    \boldsymbol{x}^*(A \circ B)\boldsymbol{x} = \sum_k \mu_k (\boldsymbol{v}_k \circ \boldsymbol{x})^* A (\boldsymbol{v}_k \circ \boldsymbol{x})
    $$

    其中 $\boldsymbol{v}_k \circ \boldsymbol{x} = \operatorname{diag}(\bar{\boldsymbol{v}}_k)\boldsymbol{x}$ 是逐元素积。

    $\leq \sum_k \mu_k \lambda_{\max}(A) \|\boldsymbol{v}_k \circ \boldsymbol{x}\|^2 = \lambda_{\max}(A) \sum_k \mu_k \sum_i |v_{ki}|^2 |x_i|^2$

    $= \lambda_{\max}(A) \sum_i |x_i|^2 \sum_k \mu_k |v_{ki}|^2 = \lambda_{\max}(A) \sum_i |x_i|^2 b_{ii}$

    $\leq \lambda_{\max}(A) \max_i b_{ii} \|\boldsymbol{x}\|^2$。

    因此 $\lambda_{\max}(A \circ B) \leq \lambda_{\max}(A) \max_i b_{ii}$。$\blacksquare$

!!! example "例 35.2"
    设 $A = B = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    $A \circ B = \begin{pmatrix} 4 & 1 \\ 1 & 4 \end{pmatrix}$。

    $A$ 的特征值：$3, 1$。$A \circ B$ 的特征值：$5, 3$。

    界：$\lambda_{\max}(A \circ B) = 5 \leq 3 \cdot 2 = 6 = \lambda_{\max}(A) \max_i b_{ii}$。✓

    $\lambda_{\min}(A \circ B) = 3 \geq 1 \cdot 2 = 2 = \lambda_{\min}(A) \min_i b_{ii}$。✓

    注意 $\lambda_{\min}(A)\lambda_{\min}(B) = 1 \cdot 1 = 1 < 3 = \lambda_{\min}(A \circ B)$。

---

## 35.3 Oppenheim 不等式

<div class="context-flow" markdown>

**核心问题**：Hadamard 积的行列式与原矩阵有什么关系？

</div>

### 主定理

!!! theorem "定理 35.6 (Oppenheim 不等式, 1930)"
    设 $A, B \in \mathbb{C}^{n \times n}$，$A \geq 0$，$B \geq 0$。则

    $$
    \det(A \circ B) \geq \det(A) \cdot \prod_{i=1}^n b_{ii} = \det(A) \cdot \prod_{i=1}^n (B)_{ii}
    $$

    等号成立当且仅当 $A$ 是对角矩阵或 $B$ 是对角矩阵（或 $\det(A) = 0$）。

    由对称性（Hadamard 积的交换律）：

    $$
    \det(A \circ B) \geq \det(B) \cdot \prod_{i=1}^n a_{ii}
    $$

??? proof "证明"
    对 $n$ 用归纳法。

    **基底**：$n = 1$ 时，$\det(A \circ B) = a_{11}b_{11} = \det(A) \cdot b_{11}$。✓

    **归纳步骤**：设命题对 $n-1$ 成立。将 $A$ 和 $B$ 分块：

    $$
    A = \begin{pmatrix} A_{11} & \boldsymbol{a} \\ \boldsymbol{a}^* & a_{nn} \end{pmatrix}, \quad B = \begin{pmatrix} B_{11} & \boldsymbol{b} \\ \boldsymbol{b}^* & b_{nn} \end{pmatrix}
    $$

    其中 $A_{11}, B_{11} \in \mathbb{C}^{(n-1) \times (n-1)}$。假设 $A_{11} > 0$（否则取极限）。

    $A$ 的 Schur 补：$A/A_{11} = a_{nn} - \boldsymbol{a}^*A_{11}^{-1}\boldsymbol{a} \geq 0$（因为 $A \geq 0$）。

    $\det(A) = \det(A_{11}) \cdot (A/A_{11})$。

    对 $A \circ B$ 做类似分块：

    $$
    A \circ B = \begin{pmatrix} A_{11} \circ B_{11} & \boldsymbol{a} \circ \boldsymbol{b} \\ (\boldsymbol{a} \circ \boldsymbol{b})^* & a_{nn}b_{nn} \end{pmatrix}
    $$

    $(A \circ B)$ 的 Schur 补：

    $$
    (A \circ B)/(A_{11} \circ B_{11}) = a_{nn}b_{nn} - (\boldsymbol{a} \circ \boldsymbol{b})^*(A_{11} \circ B_{11})^{-1}(\boldsymbol{a} \circ \boldsymbol{b})
    $$

    $\det(A \circ B) = \det(A_{11} \circ B_{11}) \cdot [(A \circ B)/(A_{11} \circ B_{11})]$。

    由归纳假设：$\det(A_{11} \circ B_{11}) \geq \det(A_{11}) \prod_{i=1}^{n-1} (B_{11})_{ii}$。

    需要证明 Schur 补部分的界。由 Schur 积定理，$A_{11} \circ B_{11} \geq 0$，且

    $$
    (\boldsymbol{a} \circ \boldsymbol{b})^*(A_{11} \circ B_{11})^{-1}(\boldsymbol{a} \circ \boldsymbol{b}) \leq \boldsymbol{a}^*A_{11}^{-1}\boldsymbol{a} \cdot b_{nn}
    $$

    这个关键不等式需要更精细的分析。利用 $B_{11}$ 的对角元素界和 Cauchy-Schwarz 不等式的推广，可以完成证明（具体细节较长，这里省略中间步骤）。

    最终得到：

    $$
    (A \circ B)/(A_{11} \circ B_{11}) \geq (A/A_{11}) \cdot b_{nn}
    $$

    因此

    $$
    \det(A \circ B) \geq \det(A_{11}) \prod_{i=1}^{n-1}(B_{11})_{ii} \cdot (A/A_{11}) \cdot b_{nn} = \det(A) \prod_{i=1}^n b_{ii}
    $$

    $\blacksquare$

!!! example "例 35.3"
    设 $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    $A \circ B = \begin{pmatrix} 6 & 1 \\ 1 & 6 \end{pmatrix}$。

    $\det(A \circ B) = 35$。

    $\det(A) \prod b_{ii} = 8 \times 4 = 32 \leq 35$。✓

    $\det(B) \prod a_{ii} = 3 \times 9 = 27 \leq 35$。✓

---

## 35.4 Hadamard 不等式

<div class="context-flow" markdown>

**核心问题**：正定矩阵的行列式与对角元素有什么关系？

</div>

### 经典 Hadamard 不等式

!!! theorem "定理 35.7 (Hadamard 不等式, 1893)"
    设 $A \in \mathbb{C}^{n \times n}$，$A \geq 0$（半正定）。则

    $$
    \det(A) \leq \prod_{i=1}^n a_{ii}
    $$

    等号成立当且仅当 $A$ 是对角矩阵（或 $A$ 的某个对角元素为 0 且对应行列全为 0）。

??? proof "证明"
    **证法一**（作为 Oppenheim 不等式的推论）：

    在 Oppenheim 不等式中取 $B = I$：

    $$
    \det(A \circ I) \geq \det(A) \cdot \prod_{i=1}^n (I)_{ii} = \det(A)
    $$

    但 $A \circ I = \operatorname{diag}(a_{11}, \ldots, a_{nn})$，故 $\det(A \circ I) = \prod a_{ii}$。

    因此 $\prod a_{ii} \geq \det(A)$。

    **证法二**（直接证明，用 Schur 补递归）：

    对 $n$ 用归纳法。$n = 1$ 显然。设对 $n-1$ 成立。

    若 $a_{nn} = 0$，由 $A \geq 0$，$A$ 的最后一行和最后一列全为 0，$\det(A) = 0 \leq \prod a_{ii}$。

    若 $a_{nn} > 0$，分块 $A = \begin{pmatrix} A_{11} & \boldsymbol{a} \\ \boldsymbol{a}^* & a_{nn} \end{pmatrix}$：

    $\det(A) = a_{nn} \det(A_{11} - \boldsymbol{a}\boldsymbol{a}^*/a_{nn})$。

    $A_{11} - \boldsymbol{a}\boldsymbol{a}^*/a_{nn} \geq 0$（Schur 补半正定）且其对角元素为 $a_{ii} - |a_i|^2/a_{nn} \leq a_{ii}$（$i = 1, \ldots, n-1$，$a_i$ 是 $\boldsymbol{a}$ 的第 $i$ 个分量）。

    由归纳假设：

    $$
    \det(A_{11} - \boldsymbol{a}\boldsymbol{a}^*/a_{nn}) \leq \prod_{i=1}^{n-1}(a_{ii} - |a_i|^2/a_{nn}) \leq \prod_{i=1}^{n-1} a_{ii}
    $$

    因此 $\det(A) \leq a_{nn} \prod_{i=1}^{n-1} a_{ii} = \prod_{i=1}^n a_{ii}$。$\blacksquare$

!!! note "注"
    Hadamard 不等式有等价的几何解释：平行多面体的体积不超过各棱长的乘积。设 $A = (\boldsymbol{a}_1, \ldots, \boldsymbol{a}_n)^T$（行向量），则 $\det(A^*A) = |\det(A)|^2$ 是由行向量 $\boldsymbol{a}_i$ 张成的平行多面体体积的平方。$(A^*A)_{ii} = \|\boldsymbol{a}_i\|^2$。因此

    $$
    |\det(A)|^2 \leq \prod \|\boldsymbol{a}_i\|^2
    $$

    即体积不超过各边长的乘积，等号当且仅当行向量正交。

### Fischer 不等式

!!! theorem "定理 35.8 (Fischer 不等式)"
    设 $A \geq 0$，分块为 $A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}$。则

    $$
    \det(A) \leq \det(A_{11}) \cdot \det(A_{22})
    $$

    这推广了 Hadamard 不等式（后者是将 $A_{11}$ 和 $A_{22}$ 都取为 $1 \times 1$ 块的极端情形）。

??? proof "证明"
    若 $A_{22}$ 奇异，由 $A \geq 0$，$\det(A) = 0 \leq \det(A_{11})\det(A_{22})$。

    若 $A_{22} > 0$：$\det(A) = \det(A_{22}) \det(A/A_{22})$，其中 $A/A_{22} = A_{11} - A_{12}A_{22}^{-1}A_{21} \geq 0$。

    由 $A/A_{22} \leq A_{11}$（因为 $A_{12}A_{22}^{-1}A_{21} \geq 0$），得 $\det(A/A_{22}) \leq \det(A_{11})$。

    因此 $\det(A) \leq \det(A_{22})\det(A_{11})$。$\blacksquare$

!!! example "例 35.4"
    设 $A = \begin{pmatrix} 4 & 1 & 2 \\ 1 & 3 & 1 \\ 2 & 1 & 5 \end{pmatrix}$（可验证 $A > 0$）。

    Hadamard 不等式：$\det(A) \leq 4 \times 3 \times 5 = 60$。

    Fischer 不等式（$A_{11} = \begin{pmatrix} 4 & 1 \\ 1 & 3 \end{pmatrix}$，$A_{22} = (5)$）：

    $\det(A) \leq \det(A_{11}) \cdot \det(A_{22}) = 11 \times 5 = 55$。

    实际 $\det(A) = 4(15-1) - 1(5-2) + 2(1-6) = 56 - 3 - 10 = 43$。

    $43 \leq 55 \leq 60$。Fischer 不等式给出了更紧的界。

---

## 35.5 Hadamard 积的谱性质

<div class="context-flow" markdown>

**核心问题**：$A \circ B$ 的特征值和奇异值如何被 $A$ 和 $B$ 的相应量控制？

</div>

### Horn 不等式

!!! theorem "定理 35.9 (Horn 不等式)"
    设 $A, B \in \mathbb{C}^{n \times n}$，$A \geq 0$，$B \geq 0$。设特征值按降序排列。则对所有 $k = 1, \ldots, n$：

    $$
    \lambda_k(A \circ B) \geq \lambda_k(A) \cdot \lambda_n(B)
    $$

    更一般地，对任何 $i + j - 1 \leq k \leq n$：

    $$
    \lambda_k(A \circ B) \geq \lambda_i(A) \cdot \lambda_j(B)
    $$

    不一定成立（需要更精细的条件）。

!!! theorem "定理 35.10 (奇异值不等式)"
    设 $A, B \in \mathbb{C}^{m \times n}$。则对所有 $k$：

    $$
    \sigma_k(A \circ B) \leq \sigma_k(A) \cdot \max_j \|\boldsymbol{b}_j\|_\infty
    $$

    其中 $\|\boldsymbol{b}_j\|_\infty = \max_i |b_{ij}|$ 是 $B$ 的第 $j$ 列的 $\ell^\infty$ 范数。

    更精确的界：

    $$
    \sigma_k(A \circ B) \leq \sigma_k(A) \cdot \max_{i,j} |b_{ij}| = \sigma_k(A) \cdot \|B\|_{\max}
    $$

### 谱半径

!!! theorem "定理 35.11"
    设 $A, B \geq 0$。则

    $$
    \rho(A \circ B) \leq \rho(A) \cdot \max_i b_{ii}
    $$

    若 $B$ 是相关矩阵（$b_{ii} = 1$），则 $\rho(A \circ B) \leq \rho(A)$。

!!! example "例 35.5"
    设 $A = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$（相关矩阵），$B = \begin{pmatrix} 1 & 0.8 \\ 0.8 & 1 \end{pmatrix}$（相关矩阵）。

    $A \circ B = \begin{pmatrix} 1 & 0.4 \\ 0.4 & 1 \end{pmatrix}$。

    $\rho(A) = 1.5$，$\rho(B) = 1.8$，$\rho(A \circ B) = 1.4$。

    $1.4 \leq 1.5 = \rho(A) \cdot \max b_{ii} = 1.5 \cdot 1$。✓

    Hadamard 积"减弱"了相关性——这正是协方差锥化（tapering）的数学基础。

---

## 35.6 正映射与完全正映射

<div class="context-flow" markdown>

**核心问题**：什么样的矩阵 $M$ 使得映射 $A \mapsto M \circ A$ 保持半正定性？

</div>

### Schur 乘子

!!! definition "定义 35.2 (Schur 乘子)"
    设 $M \in \mathbb{C}^{n \times n}$。线性映射 $\Phi_M: \mathbb{C}^{n \times n} \to \mathbb{C}^{n \times n}$ 定义为

    $$
    \Phi_M(A) = M \circ A
    $$

    称为以 $M$ 为**核**的 **Schur 乘子**（Schur multiplier）。

!!! theorem "定理 35.12 (Schur 乘子保正性的刻画)"
    以下条件等价：

    1. $\Phi_M$ 是正映射（即 $A \geq 0 \Rightarrow M \circ A \geq 0$）。
    2. $M \geq 0$（$M$ 本身半正定）。

??? proof "证明"
    **(2) $\Rightarrow$ (1)**：这就是 Schur 积定理。

    **(1) $\Rightarrow$ (2)**：取 $A = \boldsymbol{x}\boldsymbol{x}^*$（秩一半正定矩阵），则 $M \circ \boldsymbol{x}\boldsymbol{x}^* \geq 0$ 对所有 $\boldsymbol{x}$。

    取 $\boldsymbol{x} = \boldsymbol{1} = (1, \ldots, 1)^T$：$M \circ \boldsymbol{1}\boldsymbol{1}^T = M$。因此 $M \geq 0$。$\blacksquare$

### 完全正映射

!!! definition "定义 35.3 (完全正映射)"
    线性映射 $\Phi: \mathbb{C}^{n \times n} \to \mathbb{C}^{n \times n}$ 称为**完全正**的（completely positive），若对所有 $k \geq 1$，$\Phi \otimes \operatorname{id}_k: \mathbb{C}^{nk \times nk} \to \mathbb{C}^{nk \times nk}$ 也是正映射。

    即 $\Phi$ 是完全正的当且仅当对所有 $k$ 和所有 $A \geq 0 \in \mathbb{C}^{nk \times nk}$，

    $$
    (\Phi \otimes \operatorname{id}_k)(A) \geq 0
    $$

!!! theorem "定理 35.13 (Schur 乘子是完全正的)"
    若 $M \geq 0$，则 Schur 乘子 $\Phi_M$ 不仅是正映射，而且是**完全正映射**。

??? proof "证明"
    需要证明：对所有 $k$ 和 $A \geq 0 \in \mathbb{C}^{nk \times nk}$，$(M \otimes I_k) \circ A \geq 0$。

    这里 $(\Phi_M \otimes \operatorname{id}_k)(A) = (M \otimes I_k) \circ A$... 需要更仔细地定义。

    实际上，$(\Phi_M \otimes \operatorname{id}_k)$ 作用在 $\mathbb{C}^{n \times n} \otimes \mathbb{C}^{k \times k} \cong \mathbb{C}^{nk \times nk}$ 上。若将 $A$ 写成 $n \times n$ 的块矩阵 $A = (A_{ij})$，$A_{ij} \in \mathbb{C}^{k \times k}$，则

    $$
    (\Phi_M \otimes \operatorname{id}_k)(A) = (m_{ij} A_{ij})_{i,j=1}^n
    $$

    这等于 $M \circ_{\text{block}} A$（块 Hadamard 积，每个 $k \times k$ 块乘以对应的标量 $m_{ij}$）。

    由 $M \geq 0$，$M = \sum_l \mu_l \boldsymbol{v}_l \boldsymbol{v}_l^*$。

    $$
    (m_{ij}A_{ij}) = \sum_l \mu_l (v_{li}\bar{v}_{lj} A_{ij}) = \sum_l \mu_l (\operatorname{diag}(\boldsymbol{v}_l) \otimes I_k) A (\operatorname{diag}(\bar{\boldsymbol{v}}_l) \otimes I_k)
    $$

    每一项 $(\operatorname{diag}(\boldsymbol{v}_l) \otimes I_k) A (\operatorname{diag}(\bar{\boldsymbol{v}}_l) \otimes I_k) = C_l A C_l^*$，当 $A \geq 0$ 时 $C_l A C_l^* \geq 0$。

    因此 $\sum_l \mu_l C_l A C_l^* \geq 0$。$\blacksquare$

!!! note "注"
    在量子信息论中，完全正映射是量子通道的数学模型。Schur 乘子 $\Phi_M$（$M \geq 0$）对应的量子通道称为 **Schur 通道**或**相位退相干通道**（phase damping channel）。它保留密度矩阵的对角元素（概率分布），但衰减非对角元素（量子相干性）。

    Kraus 表示：$\Phi_M(A) = \sum_l \mu_l D_l A D_l^*$，其中 $D_l = \operatorname{diag}(\boldsymbol{v}_l)$。

---

## 35.7 应用

<div class="context-flow" markdown>

**核心问题**：Hadamard 积在实际问题中如何使用？

</div>

### 协方差锥化（Tapering）

!!! example "例 35.6 (空间统计中的协方差锥化)"
    在空间统计中，估计大型协方差矩阵 $\Sigma$ 时，样本协方差矩阵 $\hat{\Sigma}$ 往往条件数很大（不稳定）。**协方差锥化**通过 Hadamard 积"压缩"远离对角线的元素：

    $$
    \tilde{\Sigma} = T \circ \hat{\Sigma}
    $$

    其中 $T = (t_{ij})$ 是**锥化矩阵**，满足：

    - $t_{ii} = 1$（保留方差）
    - $|t_{ij}|$ 随 $|i-j|$ 增大而衰减到 0（衰减远距离相关）
    - $T \geq 0$（由 Schur 积定理保证 $\tilde{\Sigma} \geq 0$）

    常用的锥化函数：$t_{ij} = \max(1 - |i-j|/h, 0)$（线性锥化），$t_{ij} = e^{-|i-j|^2/h^2}$（高斯锥化）。

    Schur 积定理保证了锥化后的矩阵仍然半正定，这是锥化方法的理论基石。

!!! note "注"
    Furrer 和 Bengtsson（2007）证明了：在适当的锥化参数 $h$ 下，$\tilde{\Sigma}$ 在算子范数下比 $\hat{\Sigma}$ 更接近真实的 $\Sigma$。这说明 Hadamard 积不仅保持了半正定性，还能改善估计的统计性质。

### 神经网络中的逐元素操作

!!! example "例 35.7 (Neural Tangent Kernel)"
    在深度学习理论中，无限宽神经网络的核函数（Neural Tangent Kernel, NTK）可以通过 Hadamard 积递归定义。

    设 $K^{(0)}$ 是输入层的核矩阵。每经过一层（激活函数 $\sigma$），核矩阵更新为

    $$
    K^{(l+1)} = \dot{\Sigma}^{(l)} \circ K^{(l)} + \Sigma^{(l+1)}
    $$

    其中 $\Sigma^{(l)}$ 和 $\dot{\Sigma}^{(l)}$ 是由激活函数导出的矩阵。

    Schur 积定理保证：若 $\dot{\Sigma}^{(l)} \geq 0$ 和 $K^{(l)} \geq 0$，则乘积项 $\dot{\Sigma}^{(l)} \circ K^{(l)} \geq 0$。这为 NTK 的正定性（从而核方法的适定性）提供了理论保障。

### 图论应用

!!! example "例 35.8 (图 Laplacian 的修改)"
    设 $G$ 是无向图，Laplacian 矩阵 $L = D - A$（$D$ 是度矩阵，$A$ 是邻接矩阵）。

    通过 Hadamard 积修改边权重：设 $W = (w_{ij})$ 是权重矩阵（$w_{ij} \geq 0$），则修改后的邻接矩阵为

    $$
    A' = A \circ W
    $$

    若 $W$ 的结构使得 $W \geq 0$（作为矩阵半正定），则修改后的图具有特殊的谱性质。

    例如，热核 $W = (e^{-|i-j|^2/t})$ 对应于在图上进行热扩散后的权重修改，它保留了图的局部结构同时衰减了长距离连接。

### Hadamard 积与矩阵方程

!!! theorem "定理 35.14 (Hadamard 积与 Lyapunov 方程)"
    矩阵方程 $AX + XA^* = C$ 的解可以表示为

    $$
    \operatorname{vec}(X) = (I \otimes A + \bar{A} \otimes I)^{-1} \operatorname{vec}(C)
    $$

    当 $A = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$ 时，解简化为

    $$
    X = L \circ C, \quad l_{ij} = \frac{1}{\lambda_i + \bar{\lambda}_j}
    $$

    即 Lyapunov 方程的解是 Hadamard 积形式。矩阵 $L = (1/(\lambda_i + \bar{\lambda}_j))$ 称为 **Cauchy 矩阵**（或 **Loewner 矩阵**）。

!!! example "例 35.9"
    设 $A = \operatorname{diag}(1, 2)$，$C = \begin{pmatrix} 4 & 3 \\ 3 & 8 \end{pmatrix}$。求 $AX + XA = C$。

    $L = \begin{pmatrix} 1/(1+1) & 1/(1+2) \\ 1/(2+1) & 1/(2+2) \end{pmatrix} = \begin{pmatrix} 1/2 & 1/3 \\ 1/3 & 1/4 \end{pmatrix}$。

    $X = L \circ C = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。

    验证：$AX + XA = \begin{pmatrix} 2 & 1 \\ 2 & 4 \end{pmatrix} + \begin{pmatrix} 2 & 2 \\ 1 & 4 \end{pmatrix} = \begin{pmatrix} 4 & 3 \\ 3 & 8 \end{pmatrix} = C$。✓

!!! note "注"
    本章建立的 Hadamard 积理论展示了一个重要主题：看似简单的逐元素运算具有深刻的矩阵分析内涵。Schur 积定理（半正定的保持）是核心，Oppenheim 和 Hadamard 不等式是行列式层面的体现，而与 Kronecker 积的关系提供了理论理解的框架。

    读者应特别注意 Hadamard 积与通常矩阵乘法的本质区别：Hadamard 积保持半正定性（Schur 积定理），但通常矩阵乘法不保持——$A \geq 0$ 且 $B \geq 0$ 不蕴含 $AB \geq 0$（除非 $AB$ 是 Hermite 的）。这一对比揭示了两种矩阵乘法在正定锥上的不同几何行为。
