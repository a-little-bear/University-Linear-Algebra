# 第 16 章 正定矩阵

<div class="context-flow" markdown>

**前置**：Hermitian矩阵谱分解(Ch6)

**脉络**：$A\succ 0$ $\Leftrightarrow$ 特征值全正 $\Leftrightarrow$ Cholesky可分解 $\Leftrightarrow$ 二次型正 → **Schur补**是分块正定的钥匙 → Ch18不等式

**延伸**：正定矩阵在统计学（协方差矩阵）、机器学习（核矩阵、Fisher 信息矩阵）、优化（Newton 法的 Hessian）、弹性力学（应力-应变关系）中无处不在；Löwner 偏序在矩阵单调函数和量子信息论中有深刻应用

</div>

正定矩阵是线性代数中最重要的矩阵类之一，它在优化、统计、微分方程、信号处理和机器学习等领域无处不在。正定矩阵的理论联系了二次型、特征值、行列式、矩阵分解等多个核心概念，是线性代数应用的基石。本章系统地建立正定矩阵和半正定矩阵的理论，讨论 Schur 补、Löwner 偏序、正定矩阵的运算性质，以及若干重要的行列式不等式。

---

## 16.1 正定矩阵的定义与等价条件

<div class="context-flow" markdown>

**洞察**：8个等价条件的核心逻辑——二次型正 $\leftrightarrow$ 谱全正 $\leftrightarrow$ $A=C^*C$($C$可逆) $\leftrightarrow$ Cholesky $\leftrightarrow$ 顺序主子式全正

</div>

!!! definition "定义 16.1 (正定矩阵)"
    设 $A \in \mathbb{C}^{n \times n}$ 是 Hermitian 矩阵（$A = A^*$）。若对所有非零向量 $\mathbf{x} \in \mathbb{C}^n$，有

    $$\mathbf{x}^*A\mathbf{x} > 0,$$

    则称 $A$ 为**正定矩阵（positive definite matrix）**，记作 $A \succ 0$。

!!! theorem "定理 16.1 (正定矩阵的等价条件)"
    设 $A \in \mathbb{C}^{n \times n}$ 是 Hermitian 矩阵。以下条件等价：

    (1) $A$ 正定，即 $\mathbf{x}^*A\mathbf{x} > 0$ 对所有 $\mathbf{x} \neq \mathbf{0}$；

    (2) $A$ 的所有特征值为正：$\lambda_1, \lambda_2, \ldots, \lambda_n > 0$；

    (3) $A$ 的所有顺序主子式为正：$\det(A_k) > 0$，$k = 1, 2, \ldots, n$（其中 $A_k$ 是 $A$ 的左上角 $k \times k$ 子矩阵）；

    (4) 存在可逆矩阵 $C$ 使得 $A = C^*C$；

    (5) 存在可逆的下三角矩阵 $L$（对角元为正）使得 $A = LL^*$（Cholesky 分解）；

    (6) 存在可逆的上三角矩阵 $R$ 使得 $A = R^*R$；

    (7) $A$ 的所有主子式（不仅是顺序主子式）为正；

    (8) 存在正定矩阵 $B$ 使得 $A = B^2$（正定平方根）。

??? proof "证明"
    **(1) $\Rightarrow$ (2)**：设 $A\mathbf{v} = \lambda\mathbf{v}$（$\mathbf{v} \neq \mathbf{0}$），则 $\lambda = \frac{\mathbf{v}^*A\mathbf{v}}{\mathbf{v}^*\mathbf{v}} > 0$。

    **(2) $\Rightarrow$ (1)**：$A$ 是 Hermitian 矩阵，设谱分解 $A = U\Lambda U^*$，其中 $\Lambda = \operatorname{diag}(\lambda_1,\ldots,\lambda_n)$，$\lambda_i > 0$。对 $\mathbf{x} \neq \mathbf{0}$，令 $\mathbf{y} = U^*\mathbf{x} \neq \mathbf{0}$，则 $\mathbf{x}^*A\mathbf{x} = \mathbf{y}^*\Lambda\mathbf{y} = \sum_i\lambda_i|y_i|^2 > 0$。

    **(1) $\Rightarrow$ (3)**：$A$ 正定 $\Rightarrow$ 每个顺序主子矩阵 $A_k$ 也正定（限制到子空间）$\Rightarrow$ $A_k$ 的特征值全正 $\Rightarrow$ $\det(A_k) > 0$。

    **(3) $\Rightarrow$ (5)**：用归纳法和 Gauss 消元。对 $n = 1$：$A = (a_{11})$，$a_{11} > 0$，取 $L = (\sqrt{a_{11}})$。假设对 $n-1$ 阶成立。将 $A$ 分块：

    $$A = \begin{pmatrix} A_{n-1} & \mathbf{a} \\ \mathbf{a}^* & a_{nn} \end{pmatrix}.$$

    由归纳假设，$A_{n-1} = L_{n-1}L_{n-1}^*$。令 $\mathbf{l} = L_{n-1}^{-1}\mathbf{a}$，$\ell_{nn} = \sqrt{a_{nn} - \mathbf{l}^*\mathbf{l}}$（由 $\det A > 0$ 和 $\det A_{n-1} > 0$ 可以验证根号下为正）。则

    $$L = \begin{pmatrix} L_{n-1} & \mathbf{0} \\ \mathbf{l}^* & \ell_{nn} \end{pmatrix}, \quad A = LL^*.$$

    **(5) $\Rightarrow$ (4)**：取 $C = L^*$。

    **(4) $\Rightarrow$ (1)**：$\mathbf{x}^*A\mathbf{x} = \mathbf{x}^*C^*C\mathbf{x} = \|C\mathbf{x}\|^2 > 0$（因 $C$ 可逆，$C\mathbf{x} \neq \mathbf{0}$）。

    **(2) $\Rightarrow$ (8)**：设 $A = U\Lambda U^*$，取 $B = U\Lambda^{1/2}U^*$（其中 $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$），则 $B$ 正定且 $B^2 = A$。

    **(8) $\Rightarrow$ (1)**：$\mathbf{x}^*A\mathbf{x} = \mathbf{x}^*B^2\mathbf{x} = \|B\mathbf{x}\|^2 > 0$（$B$ 可逆）。

    **(2) $\Rightarrow$ (7)**：正定矩阵的任意主子矩阵（对应于选取某些行和相同列号）也是正定的，其行列式为正。$\blacksquare$

!!! definition "定义 16.2 (Cholesky 分解)"
    设 $A$ 是正定 Hermitian 矩阵。$A$ 的 **Cholesky 分解（Cholesky decomposition）** 是指唯一的下三角矩阵 $L$（对角元为正）使得 $A = LL^*$。

!!! example "例 16.1"
    对 $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$ 进行 Cholesky 分解。

    $l_{11} = \sqrt{4} = 2$。

    $l_{21} = a_{21}/l_{11} = 2/2 = 1$，$l_{31} = a_{31}/l_{11} = 1/2$。

    $l_{22} = \sqrt{a_{22} - l_{21}^2} = \sqrt{5 - 1} = 2$。

    $l_{32} = (a_{32} - l_{31}l_{21})/l_{22} = (3 - 1/2)/2 = 5/4$。

    $l_{33} = \sqrt{a_{33} - l_{31}^2 - l_{32}^2} = \sqrt{6 - 1/4 - 25/16} = \sqrt{6 - 29/16} = \sqrt{67/16} = \frac{\sqrt{67}}{4}$。

    $$L = \begin{pmatrix} 2 & 0 & 0 \\ 1 & 2 & 0 \\ 1/2 & 5/4 & \sqrt{67}/4 \end{pmatrix}.$$

    可以验证 $LL^T = A$。

!!! example "例 16.2"
    验证正定性的各等价条件。设 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$。

    - **条件 (2)**：特征值 $\lambda = 2 \pm 1 = 3, 1$，全正。
    - **条件 (3)**：$\det(A_1) = 2 > 0$，$\det(A) = 4 - 1 = 3 > 0$。
    - **条件 (4)**：$A = C^TC$，其中 $C = \begin{pmatrix} \sqrt{2} & -1/\sqrt{2} \\ 0 & \sqrt{3/2} \end{pmatrix}$。
    - **条件 (1)**：$\mathbf{x}^TA\mathbf{x} = 2x_1^2 - 2x_1x_2 + 2x_2^2 = (x_1 - x_2)^2 + x_1^2 + x_2^2 > 0$（对 $\mathbf{x} \neq \mathbf{0}$）。

---

## 16.2 半正定矩阵

<div class="context-flow" markdown>

**脉络**：$A\succeq 0$ $\Leftrightarrow$ $\lambda_i\geq 0$ $\Leftrightarrow$ $A=C^*C$ · 半正定锥 $\mathbb{S}_+^n$ 是凸锥 → Ch18 Löwner偏序与优化

</div>

!!! definition "定义 16.3 (半正定矩阵)"
    设 $A \in \mathbb{C}^{n \times n}$ 是 Hermitian 矩阵。若对所有 $\mathbf{x} \in \mathbb{C}^n$，有

    $$\mathbf{x}^*A\mathbf{x} \geq 0,$$

    则称 $A$ 为**半正定矩阵（positive semidefinite matrix）**，记作 $A \succeq 0$。

!!! theorem "定理 16.2 (半正定的等价条件)"
    设 $A$ 是 $n \times n$ Hermitian 矩阵。以下条件等价：

    (1) $A$ 半正定；

    (2) $A$ 的所有特征值非负：$\lambda_i \geq 0$；

    (3) 存在矩阵 $C$（不必方阵或可逆）使得 $A = C^*C$；

    (4) $A$ 的所有主子式非负；

    (5) 存在半正定矩阵 $B$ 使得 $A = B^2$。

??? proof "证明"
    证明与正定情形类似，将严格不等号改为非严格不等号。

    **(1) $\Rightarrow$ (2)**：若 $A\mathbf{v} = \lambda\mathbf{v}$，$\mathbf{v} \neq \mathbf{0}$，则 $\lambda = \frac{\mathbf{v}^*A\mathbf{v}}{\|\mathbf{v}\|^2} \geq 0$。

    **(2) $\Rightarrow$ (3)**：设 $A = U\Lambda U^*$，取 $C = \Lambda^{1/2}U^*$，则 $C^*C = U\Lambda U^* = A$。

    **(3) $\Rightarrow$ (1)**：$\mathbf{x}^*A\mathbf{x} = \|C\mathbf{x}\|^2 \geq 0$。$\blacksquare$

!!! definition "定义 16.4 (半正定锥)"
    全体 $n \times n$ 半正定矩阵构成的集合记为 $\mathbb{S}_+^n$，称为**半正定锥（positive semidefinite cone）**。它是 $\mathbb{S}^n$（$n \times n$ 实对称矩阵空间）中的一个凸锥：

    - 若 $A, B \in \mathbb{S}_+^n$，则 $A + B \in \mathbb{S}_+^n$；
    - 若 $A \in \mathbb{S}_+^n$，$\alpha \geq 0$，则 $\alpha A \in \mathbb{S}_+^n$。

!!! example "例 16.3"
    矩阵 $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$ 是半正定的。

    特征值：$\lambda^2 - 5\lambda = 0$，$\lambda_1 = 5$，$\lambda_2 = 0$。全非负。

    $\mathbf{x}^TA\mathbf{x} = x_1^2 + 4x_1x_2 + 4x_2^2 = (x_1 + 2x_2)^2 \geq 0$。

    $A$ 不正定因为取 $\mathbf{x} = (-2, 1)^T$ 时 $\mathbf{x}^TA\mathbf{x} = 0$。

    $A = C^TC$，其中 $C = (1, 2)$（$1 \times 2$ 矩阵）。

---

## 16.3 Schur 补

<div class="context-flow" markdown>

**核心工具**：$M/A = C-B^*A^{-1}B$ · 分块LDL分解 → $M\succ 0 \Leftrightarrow A\succ 0$ 且 $M/A\succ 0$ · 行列式公式 $\det M=\det A\cdot\det(M/A)$ → Fischer不等式

</div>

Schur 补是分析分块矩阵正定性的核心工具，在控制理论、统计学和优化中有广泛应用。

!!! definition "定义 16.5 (Schur 补)"
    设分块矩阵 $M = \begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$，其中 $A$ 可逆。$A$ 在 $M$ 中的 **Schur 补（Schur complement）** 定义为

    $$M/A = C - B^*A^{-1}B.$$

    类似地，若 $C$ 可逆，$C$ 在 $M$ 中的 Schur 补为 $M/C = A - BC^{-1}B^*$。

!!! theorem "定理 16.3 (分块矩阵的正定性判别)"
    设 $M = \begin{pmatrix} A & B \\ B^* & C \end{pmatrix}$ 是 Hermitian 矩阵。

    (1) 若 $A \succ 0$，则 $M \succ 0 \Leftrightarrow M/A = C - B^*A^{-1}B \succ 0$；

    (2) 若 $C \succ 0$，则 $M \succ 0 \Leftrightarrow M/C = A - BC^{-1}B^* \succ 0$；

    (3) $M \succeq 0$ 且 $A \succ 0$ $\Leftrightarrow$ $A \succ 0$ 且 $M/A \succeq 0$。

??? proof "证明"
    **(1)** 进行分块 LDL 分解。由于 $A \succ 0$，可以写

    $$M = \begin{pmatrix} I & O \\ B^*A^{-1} & I \end{pmatrix}\begin{pmatrix} A & O \\ O & C-B^*A^{-1}B \end{pmatrix}\begin{pmatrix} I & A^{-1}B \\ O & I \end{pmatrix}.$$

    可直接验证右端乘积等于 $M$。

    设 $L = \begin{pmatrix} I & O \\ B^*A^{-1} & I \end{pmatrix}$，$D = \begin{pmatrix} A & O \\ O & M/A \end{pmatrix}$。则 $M = LDL^*$。

    由于 $L$ 可逆，$M$ 正定 $\Leftrightarrow$ $D$ 正定 $\Leftrightarrow$ $A \succ 0$ 且 $M/A \succ 0$。已知 $A \succ 0$，故 $M \succ 0 \Leftrightarrow M/A \succ 0$。

    **(2)** 类似，对 $C$ 做分块消元。$\blacksquare$

!!! theorem "定理 16.4 (Schur 补的行列式公式)"
    设 $A$ 可逆，则

    $$\det\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \det(A)\cdot\det(D - CA^{-1}B).$$

??? proof "证明"
    $$\begin{pmatrix} A & B \\ C & D \end{pmatrix} = \begin{pmatrix} I & O \\ CA^{-1} & I \end{pmatrix}\begin{pmatrix} A & B \\ O & D - CA^{-1}B \end{pmatrix}.$$

    取行列式：$\det(M) = 1 \cdot \det(A)\det(D - CA^{-1}B)$。$\blacksquare$

!!! example "例 16.4"
    设 $M = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 3 & 1 \\ 1 & 1 & 2 \end{pmatrix}$。判断 $M$ 是否正定。

    取 $A = \begin{pmatrix} 4 & 2 \\ 2 & 3 \end{pmatrix}$（正定，因 $4 > 0$，$\det = 8 > 0$），$B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$，$C = (2)$。

    $M/A = C - B^TA^{-1}B = 2 - (1, 1)\frac{1}{8}\begin{pmatrix} 3 & -2 \\ -2 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 1 \end{pmatrix} = 2 - \frac{1}{8}(1, 1)\begin{pmatrix} 1 \\ 2 \end{pmatrix} = 2 - \frac{3}{8} = \frac{13}{8} > 0$。

    因此 $M \succ 0$。

!!! example "例 16.5"
    **线性回归中的 Schur 补**。在最小二乘问题中，正规方程组的系数矩阵 $X^TX$ 是半正定的。当对模型进行分块时（如分为截距项和其他变量），Schur 补自然出现。

    设 $X = (\mathbf{1}, X_2)$，则

    $$(X^TX)/ ((\mathbf{1}^T\mathbf{1})) = X_2^TX_2 - X_2^T\mathbf{1}(\mathbf{1}^T\mathbf{1})^{-1}\mathbf{1}^TX_2 = X_2^T(I - \frac{1}{n}\mathbf{1}\mathbf{1}^T)X_2,$$

    即中心化后的 $X_2^TX_2$，它是样本协方差矩阵（不含因子 $\frac{1}{n-1}$）。

---

## 16.4 Löwner 偏序

<div class="context-flow" markdown>

**脉络**：$A\succeq B \Leftrightarrow A-B\succeq 0$ 定义矩阵上的偏序 · 合同不变性 + **逆序性**($A\succeq B\succ 0 \Rightarrow B^{-1}\succeq A^{-1}$) → Ch18 矩阵单调函数

</div>

!!! definition "定义 16.6 (Löwner 偏序)"
    对 Hermitian 矩阵 $A$ 和 $B$，定义 **Löwner 偏序（Löwner partial order）**：

    $$A \succeq B \quad \Longleftrightarrow \quad A - B \succeq 0 \quad (\text{即 } A - B \text{ 半正定}),$$

    $$A \succ B \quad \Longleftrightarrow \quad A - B \succ 0 \quad (\text{即 } A - B \text{ 正定}).$$

!!! theorem "定理 16.5 (Löwner 偏序的性质)"
    设 $A, B, C$ 是同阶 Hermitian 矩阵。Löwner 偏序满足：

    (1) **自反性**：$A \succeq A$；

    (2) **反对称性**：若 $A \succeq B$ 且 $B \succeq A$，则 $A = B$；

    (3) **传递性**：若 $A \succeq B$ 且 $B \succeq C$，则 $A \succeq C$；

    (4) **加法保持**：若 $A \succeq B$，则 $A + C \succeq B + C$；

    (5) **合同不变性**：若 $A \succeq B$，则对任意矩阵 $P$，$P^*AP \succeq P^*BP$；

    (6) **逆序性**：若 $A \succeq B \succ 0$，则 $B^{-1} \succeq A^{-1}$。

??? proof "证明"
    (1)-(4) 直接由半正定锥的性质得出。

    **(5)** $P^*AP - P^*BP = P^*(A-B)P$。由 $A - B \succeq 0$，对任意 $\mathbf{x}$，$\mathbf{x}^*P^*(A-B)P\mathbf{x} = (P\mathbf{x})^*(A-B)(P\mathbf{x}) \geq 0$。

    **(6)** 由 $A \succeq B \succ 0$，对任意 $\mathbf{x} \neq \mathbf{0}$：

    $$\mathbf{x}^*B^{-1}\mathbf{x} = \mathbf{x}^*B^{-1/2}B^{1/2}B^{-1}B^{1/2}B^{-1/2}\mathbf{x} = (B^{-1/2}\mathbf{x})^*I(B^{-1/2}\mathbf{x}).$$

    更直接地：$A \succeq B \Rightarrow B^{-1/2}AB^{-1/2} \succeq I \Rightarrow (B^{-1/2}AB^{-1/2})^{-1} \preceq I$（因为 $f(t) = 1/t$ 在正半轴上是算子单调递减的）$\Rightarrow B^{1/2}A^{-1}B^{1/2} \preceq I \Rightarrow A^{-1} \preceq B^{-1}$。$\blacksquare$

!!! example "例 16.6"
    设 $A = \begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$。

    $A - B = \begin{pmatrix} 1 & 1 \\ 1 & 2 \end{pmatrix}$，特征值为 $\frac{3 \pm \sqrt{5}}{2} > 0$，因此 $A \succ B$。

    由逆序性，$B^{-1} \succeq A^{-1}$。验证：$B^{-1} = \begin{pmatrix} 1/2 & 0 \\ 0 & 1 \end{pmatrix}$，$A^{-1} = \frac{1}{8}\begin{pmatrix} 3 & -1 \\ -1 & 3 \end{pmatrix}$。

    $B^{-1} - A^{-1} = \begin{pmatrix} 1/2 - 3/8 & 1/8 \\ 1/8 & 1 - 3/8 \end{pmatrix} = \begin{pmatrix} 1/8 & 1/8 \\ 1/8 & 5/8 \end{pmatrix}$，行列式 $= 5/64 - 1/64 = 4/64 > 0$，且 $1/8 > 0$，故半正定（正定）。

---

## 16.5 正定矩阵的运算性质

<div class="context-flow" markdown>

**脉络**：和/数乘保正定 · $AB$特征值全正(但未必Hermitian)

**Schur积定理**($A\circ B\succeq 0$)

**Kronecker积**保正定 → Ch19

</div>

!!! theorem "定理 16.6 (和与数乘)"
    (1) 若 $A \succ 0$ 且 $B \succ 0$，则 $A + B \succ 0$；

    (2) 若 $A \succ 0$ 且 $\alpha > 0$，则 $\alpha A \succ 0$。

!!! theorem "定理 16.7 (积的正定性)"
    设 $A, B$ 是正定 Hermitian 矩阵。则 $AB$ 的特征值全为正实数，但 $AB$ 不一定是 Hermitian 矩阵（除非 $AB = BA$）。

??? proof "证明"
    $AB$ 与 $A^{1/2}BA^{1/2}$ 相似（因为 $AB = A^{1/2}(A^{1/2}B)$ 与 $A^{1/2}BA^{1/2}$ 有相同特征值）。$A^{1/2}BA^{1/2}$ 是正定 Hermitian 矩阵（由合同不变性），其特征值全正。因此 $AB$ 的特征值全正。$\blacksquare$

!!! definition "定义 16.7 (Hadamard 积)"
    设 $A = (a_{ij})$，$B = (b_{ij}) \in \mathbb{C}^{m \times n}$。它们的 **Hadamard 积（Hadamard product，也称 Schur 积）** 定义为

    $$A \circ B = (a_{ij}b_{ij}),$$

    即对应元素相乘。

!!! theorem "定理 16.8 (Schur 积定理)"
    若 $A \succeq 0$ 且 $B \succeq 0$，则 $A \circ B \succeq 0$。若进一步 $A \succ 0$ 且 $B \succ 0$，则 $A \circ B \succ 0$。

??? proof "证明"
    由 $B \succeq 0$，设 $B = \sum_{k=1}^{r} \mathbf{b}_k\mathbf{b}_k^*$（秩一分解）。则

    $$A \circ B = A \circ \left(\sum_k \mathbf{b}_k\mathbf{b}_k^*\right) = \sum_k A \circ (\mathbf{b}_k\mathbf{b}_k^*) = \sum_k D_k A D_k^*,$$

    其中 $D_k = \operatorname{diag}(\mathbf{b}_k)$（以 $\mathbf{b}_k$ 的分量为对角元的对角矩阵）。

    这是因为 $(A \circ \mathbf{b}\mathbf{b}^*)_{ij} = a_{ij}b_ib_j^* = (DAD^*)_{ij}$，其中 $D = \operatorname{diag}(b_1, \ldots, b_n)$。

    每个 $D_kAD_k^*$ 半正定（由合同不变性），它们的和也半正定。

    若 $A \succ 0$ 且 $B \succ 0$，由 $B$ 的对角元全正（$b_{ii} > 0$），至少有一个 $D_k$ 可逆，使得 $D_kAD_k^*$ 正定。$\blacksquare$

!!! theorem "定理 16.9 (Kronecker 积保正定性)"
    若 $A \succ 0$ 且 $B \succ 0$，则 $A \otimes B \succ 0$。

??? proof "证明"
    设 $A$ 的特征值为 $\alpha_1, \ldots, \alpha_m > 0$，$B$ 的特征值为 $\beta_1, \ldots, \beta_n > 0$。$A \otimes B$ 的特征值为 $\alpha_i\beta_j > 0$。或者：$A \otimes B$ 是 Hermitian 矩阵，且 $(\mathbf{x} \otimes \mathbf{y})^*(A \otimes B)(\mathbf{x} \otimes \mathbf{y}) = (\mathbf{x}^*A\mathbf{x})(\mathbf{y}^*B\mathbf{y}) > 0$。一般向量由这类张量积线性组合而成，可用类似推理完成。$\blacksquare$

!!! example "例 16.7"
    设 $A = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} \succ 0$，$B = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 2 \end{pmatrix} \succ 0$。

    **Hadamard 积**：$A \circ B = \begin{pmatrix} 2 & 0.5 \\ 0.5 & 6 \end{pmatrix}$。

    验证正定性：$2 > 0$，$\det = 12 - 0.25 = 11.75 > 0$，故 $A \circ B \succ 0$。

    **普通矩阵积**：$AB = \begin{pmatrix} 2.5 & 3 \\ 2.5 & 6.5 \end{pmatrix}$，不是对称矩阵。特征值为 $\frac{9 \pm \sqrt{49-2.5}}{2}$... 但由定理 16.7，其特征值全正。

---

## 16.6 Hadamard 不等式

<div class="context-flow" markdown>

**洞察**：$\det A\leq\prod a_{ii}$，等号 $\Leftrightarrow$ 对角矩阵 · 两种证法：Cholesky分解 / AM-GM不等式 → Fischer不等式的 $1\times 1$ 特例

</div>

!!! theorem "定理 16.10 (Hadamard 不等式)"
    设 $A = (a_{ij}) \in \mathbb{C}^{n \times n}$ 是正定 Hermitian 矩阵。则

    $$\det(A) \leq \prod_{i=1}^{n} a_{ii},$$

    等号成立当且仅当 $A$ 是对角矩阵。

??? proof "证明"
    **证法一（利用 Cholesky 分解）**：设 $A = LL^*$，$L = (\ell_{ij})$。则 $a_{ii} = \sum_{k=1}^{i}|\ell_{ik}|^2 \geq |\ell_{ii}|^2 = \ell_{ii}^2$。

    $\det(A) = (\det L)^2 = \left(\prod_{i=1}^{n}\ell_{ii}\right)^2 = \prod_{i=1}^{n}\ell_{ii}^2 \leq \prod_{i=1}^{n}a_{ii}$。

    等号成立 $\Leftrightarrow$ 每个 $a_{ii} = \ell_{ii}^2$ $\Leftrightarrow$ $L$ 是对角矩阵 $\Leftrightarrow$ $A$ 是对角矩阵。

    **证法二（利用 AM-GM 不等式）**：设 $D = \operatorname{diag}(a_{11}, \ldots, a_{nn})$，则 $D^{-1/2}AD^{-1/2}$ 是正定矩阵且对角元全为 $1$。设其特征值为 $\mu_1, \ldots, \mu_n > 0$。由 $\operatorname{tr}(D^{-1/2}AD^{-1/2}) = n$，即 $\sum\mu_i = n$。由 AM-GM 不等式：

    $$\frac{\det(A)}{\prod a_{ii}} = \det(D^{-1/2}AD^{-1/2}) = \prod\mu_i \leq \left(\frac{\sum\mu_i}{n}\right)^n = 1.$$

    等号当且仅当 $\mu_1 = \cdots = \mu_n = 1$，即 $D^{-1/2}AD^{-1/2} = I$，即 $A = D$。$\blacksquare$

!!! example "例 16.8"
    设 $A = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$（已验证正定）。

    $\det(A) = 4(30 - 9) - 2(12 - 3) + 1(6 - 5) = 84 - 18 + 1 = 67$。

    $\prod a_{ii} = 4 \times 5 \times 6 = 120$。

    $67 \leq 120$，Hadamard 不等式成立。比值 $67/120 \approx 0.558$，反映了 $A$ 偏离对角矩阵的程度。

---

## 16.7 Fischer 不等式

<div class="context-flow" markdown>

**脉络**：Hadamard不等式的分块推广 · $\det A\leq\det A_{11}\cdot\det A_{22}$，Schur补提供证明 · 等号 $\Leftrightarrow$ 非对角块为零 → Ch18行列式不等式

</div>

Fischer 不等式是 Hadamard 不等式的推广，将对角元替换为分块对角子矩阵的行列式。

!!! theorem "定理 16.11 (Fischer 不等式)"
    设正定 Hermitian 矩阵 $A$ 分块为

    $$A = \begin{pmatrix} A_{11} & A_{12} \\ A_{12}^* & A_{22} \end{pmatrix},$$

    其中 $A_{11}$ 为 $k \times k$，$A_{22}$ 为 $(n-k) \times (n-k)$。则

    $$\det(A) \leq \det(A_{11})\cdot\det(A_{22}),$$

    等号成立当且仅当 $A_{12} = O$（即 $A$ 是分块对角矩阵）。

??? proof "证明"
    由 Schur 补的行列式公式：

    $$\det(A) = \det(A_{11})\cdot\det(A_{22} - A_{12}^*A_{11}^{-1}A_{12}).$$

    由正定性，Schur 补 $S = A_{22} - A_{12}^*A_{11}^{-1}A_{12} \succ 0$。又 $A_{22} - S = A_{12}^*A_{11}^{-1}A_{12} \succeq 0$，即 $A_{22} \succeq S$。

    因此 $A_{22} - S \succeq 0$，其特征值非负。由 $\det(A_{22}) = \det(S)\det(S^{-1}A_{22})$，其中 $S^{-1}A_{22}$ 的特征值都 $\geq 1$（因为 $A_{22} \succeq S \succ 0$ 意味着 $S^{-1/2}A_{22}S^{-1/2} \succeq I$）。故 $\det(S^{-1/2}A_{22}S^{-1/2}) \geq 1$，即 $\det(A_{22}) \geq \det(S)$。

    因此 $\det(A) = \det(A_{11})\det(S) \leq \det(A_{11})\det(A_{22})$。

    等号当且仅当 $S = A_{22}$，即 $A_{12}^*A_{11}^{-1}A_{12} = O$。由 $A_{11} \succ 0$，这等价于 $A_{12} = O$。$\blacksquare$

!!! corollary "推论 16.1"
    对正定矩阵 $A$，更一般地可以将 $A$ 分成多个对角块 $A_{11}, A_{22}, \ldots, A_{kk}$，则

    $$\det(A) \leq \prod_{i=1}^{k}\det(A_{ii}).$$

    Hadamard 不等式是 $k = n$（每块为 $1 \times 1$）的特殊情形。

!!! example "例 16.9"
    设 $A = \begin{pmatrix} 2 & 1 & 0.5 & 0.1 \\ 1 & 3 & 0.3 & 0.2 \\ 0.5 & 0.3 & 4 & 1 \\ 0.1 & 0.2 & 1 & 5 \end{pmatrix}$，分块 $A_{11} = \begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix}$，$A_{22} = \begin{pmatrix} 4 & 1 \\ 1 & 5 \end{pmatrix}$。

    $\det(A_{11}) = 5$，$\det(A_{22}) = 19$。

    Fischer 不等式给出 $\det(A) \leq 5 \times 19 = 95$。

    Hadamard 不等式给出 $\det(A) \leq 2 \times 3 \times 4 \times 5 = 120$。

    Fischer 不等式给出了更紧的界。

---

## 16.8 正定矩阵与最优化

<div class="context-flow" markdown>

**脉络**：Hessian $A\succ 0$ $\Leftrightarrow$ 严格凸 $\Leftrightarrow$ 唯一极小值 · 水平集 = 椭球(半轴 $1/\sqrt{\lambda_i}$) · 正定性判定 $\to$ 凸优化全局最优保证

</div>

正定矩阵在最优化理论中扮演核心角色。二次函数的凸性由 Hessian 矩阵的正定性决定，正定矩阵的水平集是椭球。

!!! definition "定义 16.8 (二次函数)"
    设 $A$ 是 $n \times n$ 实对称矩阵，$\mathbf{b} \in \mathbb{R}^n$，$c \in \mathbb{R}$。**二次函数（quadratic function）** 定义为

    $$f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - \mathbf{b}^T\mathbf{x} + c.$$

!!! theorem "定理 16.12 (二次函数的极值)"
    设 $f(\mathbf{x}) = \frac{1}{2}\mathbf{x}^TA\mathbf{x} - \mathbf{b}^T\mathbf{x} + c$，其中 $A$ 对称。

    (1) $f$ 有唯一全局极小值当且仅当 $A \succ 0$。此时极小值点为 $\mathbf{x}^* = A^{-1}\mathbf{b}$，极小值为 $f(\mathbf{x}^*) = c - \frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b}$。

    (2) $f$ 是凸函数当且仅当 $A \succeq 0$。

??? proof "证明"
    $\nabla f(\mathbf{x}) = A\mathbf{x} - \mathbf{b}$，Hessian 矩阵 $\nabla^2 f = A$。

    **(1)** 若 $A \succ 0$，令 $\nabla f = \mathbf{0}$ 得 $\mathbf{x}^* = A^{-1}\mathbf{b}$。由 $\nabla^2 f = A \succ 0$，$\mathbf{x}^*$ 是严格局部极小值。由凸性（$A \succ 0$ 意味着严格凸），它也是唯一全局极小值。

    $f(\mathbf{x}) = f(\mathbf{x}^*) + \frac{1}{2}(\mathbf{x}-\mathbf{x}^*)^TA(\mathbf{x}-\mathbf{x}^*)$（配方），其中

    $f(\mathbf{x}^*) = \frac{1}{2}\mathbf{b}^TA^{-1}AA^{-1}\mathbf{b} - \mathbf{b}^TA^{-1}\mathbf{b} + c = -\frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b} + c$。

    **(2)** $f$ 凸 $\Leftrightarrow$ 对所有 $\mathbf{x}$，$\nabla^2f(\mathbf{x}) = A \succeq 0$。$\blacksquare$

!!! definition "定义 16.9 (椭球)"
    设 $A \succ 0$，$\mathbf{c} \in \mathbb{R}^n$。以 $\mathbf{c}$ 为中心、$A$ 为"形状矩阵"的**椭球（ellipsoid）** 定义为

    $$\mathcal{E}(A, \mathbf{c}) = \{\mathbf{x} \in \mathbb{R}^n : (\mathbf{x} - \mathbf{c})^TA(\mathbf{x} - \mathbf{c}) \leq 1\}.$$

    椭球的半轴长度为 $A$ 的特征值的倒数平方根 $1/\sqrt{\lambda_i}$，方向为对应的特征向量。体积为 $\operatorname{vol}(\mathcal{E}) = \frac{\omega_n}{\sqrt{\det A}}$，其中 $\omega_n$ 是 $n$ 维单位球的体积。

!!! theorem "定理 16.13 (正定矩阵与凸优化)"
    设约束优化问题 $\min_{\mathbf{x}} f(\mathbf{x})$，$g_i(\mathbf{x}) \leq 0$（$i = 1, \ldots, m$）。若 $f$ 和 $g_i$ 都是二次函数，且 $f$ 的 Hessian 矩阵 $A$ 正定，$g_i$ 的 Hessian 矩阵半正定，则这是一个凸二次规划问题，任意局部最优解也是全局最优解。

!!! example "例 16.10"
    求 $f(\mathbf{x}) = 3x_1^2 + 2x_1x_2 + 2x_2^2 - 4x_1 - 6x_2$ 的最小值。

    写成矩阵形式：$f(\mathbf{x}) = \mathbf{x}^T\begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}\mathbf{x} - \begin{pmatrix} 4 \\ 6 \end{pmatrix}^T\mathbf{x}$。

    $A = \begin{pmatrix} 3 & 1 \\ 1 & 2 \end{pmatrix}$，特征值 $\frac{5 \pm \sqrt{5}}{2} > 0$，正定。

    $\mathbf{x}^* = A^{-1}\mathbf{b} = \frac{1}{5}\begin{pmatrix} 2 & -1 \\ -1 & 3 \end{pmatrix}\begin{pmatrix} 4 \\ 6 \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 2 \\ 14 \end{pmatrix} = \begin{pmatrix} 2/5 \\ 14/5 \end{pmatrix}$。

    $f(\mathbf{x}^*) = -\frac{1}{2}\mathbf{b}^TA^{-1}\mathbf{b} = -\frac{1}{2}\cdot\frac{1}{5}(4,6)\begin{pmatrix} 2 \\ 14 \end{pmatrix} = -\frac{1}{10}(8+84) = -\frac{92}{10} = -9.2$。

!!! example "例 16.11"
    **椭球的几何描述**。设 $A = \begin{pmatrix} 4 & 0 \\ 0 & 1 \end{pmatrix}$，中心 $\mathbf{c} = \mathbf{0}$。

    椭球方程：$4x_1^2 + x_2^2 \leq 1$，即 $\frac{x_1^2}{(1/2)^2} + \frac{x_2^2}{1^2} \leq 1$。

    半轴长度：$1/\sqrt{4} = 1/2$（沿 $x_1$ 轴）和 $1/\sqrt{1} = 1$（沿 $x_2$ 轴）。

    面积 $= \pi \times \frac{1}{2} \times 1 = \frac{\pi}{2} = \frac{\pi}{\sqrt{\det A}} = \frac{\pi}{\sqrt{4}} = \frac{\pi}{2}$。

## 练习题

1. **[定义] 一个实矩阵满足 $\mathbf{x}^TA\mathbf{x} > 0$ 对所有非零 $\mathbf{x}$ 成立，它必定是正定对称矩阵吗？**
   ??? success "参考答案"
       不一定。例如 $A = \begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$，$\mathbf{x}^TA\mathbf{x} = x_1^2 + x_2^2 > 0$，但它不对称。在现代矩阵分析中，“正定矩阵”通常隐含了必须是 Hermite（对称）矩阵的定义前提。

2. **[等价条件] 正定矩阵的行列式有什么特点？**
   ??? success "参考答案"
       必定是大于零的实数，因为它是所有（严格为正的）特征值的乘积。

3. **[分解] 若 $A$ 是正定矩阵，证明存在唯一的可逆下三角矩阵 $L$ 且对角元全正，使得 $A = LL^T$。**
   ??? success "参考答案"
       这就是 Cholesky 分解。其存在性通过顺序主子式全正来保证，其唯一性则通过高斯消元和对角矩阵元素的非负开方唯一性来确立。

4. **[半正定] 两个半正定矩阵之和必定是半正定矩阵吗？**
   ??? success "参考答案"
       是的。对任意 $\mathbf{x}$，$\mathbf{x}^T(A+B)\mathbf{x} = \mathbf{x}^TA\mathbf{x} + \mathbf{x}^TB\mathbf{x} \ge 0 + 0 = 0$。

5. **[运算] 证明如果 $A, B$ 是半正定矩阵，则它们的迹 $\operatorname{tr}(AB) \ge 0$。**
   ??? success "参考答案"
       因为 $A = C^TC$，则 $\operatorname{tr}(AB) = \operatorname{tr}(C^TCB) = \operatorname{tr}(CBC^T)$。由于 $B \succeq 0$，对角元素非负，因此迹非负。

6. **[偏序] Löwner 偏序 $A \succeq B$ 和矩阵元素序 $a_{ij} \ge b_{ij}$ 有何联系？**
   ??? success "参考答案"
       毫无联系！例如 $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ 元素皆为正，但它不是半正定的（特征值为 $3, -1$），而 $B = I$ 元素不全为正，却是半正定的。但对角元满足 $a_{ii} \ge b_{ii}$（取标准基向量）。

7. **[偏序] 证明若 $A \succeq B \succ 0$，则 $B^{-1} \succeq A^{-1}$。**
   ??? success "参考答案"
       这是著名的“逆序性”。通过全等变换和算子单调性 $f(t)=1/t$ 可严格证明。

8. **[Schur补] 在计算多元正态分布的条件方差时，条件协方差矩阵的形式为什么恰好是 Schur 补？**
   ??? success "参考答案"
       因为高斯分布的指数项是一个正定二次型，对部分变量求配方并固定后，其残余部分的二次型系数矩阵正是分块消元后的剩余块，即 Schur 补。

9. **[优化] 为什么凸优化问题的 Hessian 矩阵必须处处半正定？**
   ??? success "参考答案"
       因为泰勒展开的二阶项控制了函数值的局部曲率。半正定保证了函数在任何方向上都不存在“向下弯曲”的趋势，从而保证任何局部极小值都是全局极小值。

10. **[爱因斯坦思考题] 在广义相对论中，时空的“度规张量” $g_{\mu\nu}$ 扮演了极为重要的角色。为什么在一个局部的随动参考系（自由落体电梯）中，只涉及空间部分的子度规张量 $g_{ij}$ 必须是一个正定矩阵？**
    ??? success "参考答案"
        因为物理空间的“长度平方” $ds^2 = g_{ij}dx^i dx^j$ 必须是一个绝对正数，否则空间将退化或允许出现“距离为虚数”的幽灵维度。正定性是时空因果律在纯粹空间切片上保持欧几里得几何特性的终极数学保证，它排除了空间内部的“光锥”。

## 本章小结

本章系统地论述了正定矩阵和半正定矩阵，主要内容包括：

1. **正定性判据**：通过特征值、二次型、Cholesky 分解、矩阵合同和顺序主子式五大维度，给出了刻画矩阵正定性的等价条件。
2. **Schur 补与分块**：证明了分块矩阵正定性的递归下降法则，引入了在统计学（条件分布）和控制论中极其重要的 Schur 补。
3. **Löwner 偏序**：基于半正定锥定义了矩阵上的偏序关系，并探讨了它在加法、合同变换和求逆下的良好性质（逆序性）。
4. **Hadamard 与 Fischer 不等式**：利用正定性导出了行列式的几个著名上限估计，证明了对角矩阵在保持行列式上的极端地位。
5. **最优化几何**：将正定矩阵与二次凸函数及其水平集（椭球）联系起来，为后续的高级优化理论奠定了基础。
