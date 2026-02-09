# 第 34 章 Schur 补

<div class="context-flow" markdown>

**前置**：分块矩阵(Ch2) · 行列式(Ch3) · 正定矩阵(Ch16)

**本章脉络**：Schur 补定义 → 块消元 → 行列式公式 → 正定性判定 → Sherman-Morrison-Woodbury → 商性质与单调性 → 统计应用 → 优化应用

**延伸**：Schur 补在半定规划（线性矩阵不等式 LMI）、统计学（条件分布的协方差 = Schur 补）、数值方法（区域分解法、预条件子构造）中是最核心的矩阵工具之一

</div>

Schur 补是分块矩阵理论的核心概念。当我们对分块矩阵进行块高斯消元时，Schur 补自然出现。它将大矩阵的性质（如行列式、特征值、正定性）化为较小矩阵的性质，是"分而治之"思想在线性代数中的体现。

Schur 补以 Issai Schur（1917）命名，他首先在行列式理论中使用了这一概念。此后 Schur 补在统计学（条件分布）、控制论（Riccati 方程）、优化（半定规划）和数值计算（区域分解）等众多领域中成为不可或缺的工具。

---

## 34.1 Schur 补的定义

<div class="context-flow" markdown>

**核心问题**：如何将分块矩阵的性质化简为较小块的性质？

</div>

### 基本定义

!!! definition "定义 34.1 (Schur 补)"
    设 $M \in \mathbb{C}^{(p+q) \times (p+q)}$ 是分块矩阵

    $$
    M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}
    $$

    其中 $A \in \mathbb{C}^{p \times p}$，$B \in \mathbb{C}^{p \times q}$，$C \in \mathbb{C}^{q \times p}$，$D \in \mathbb{C}^{q \times q}$。

    1. 若 $A$ 可逆，则 $D$ 关于 $A$ 在 $M$ 中的 **Schur 补**定义为
    $$
    M/A = D - CA^{-1}B
    $$

    2. 若 $D$ 可逆，则 $A$ 关于 $D$ 在 $M$ 中的 **Schur 补**定义为
    $$
    M/D = A - BD^{-1}C
    $$

!!! note "注"
    Schur 补的名称和记号 $M/A$ 由 Crabtree 和 Haynsworth（1969）引入。读作"$M$ 对 $A$ 的 Schur 补"或"$A$ 在 $M$ 中的 Schur 补"。注意 $M/A$ 实际上是消去 $A$ 后剩余的部分。

!!! example "例 34.1"
    设 $M = \begin{pmatrix} 2 & 1 & 3 \\ 0 & 4 & 1 \\ 1 & 2 & 5 \end{pmatrix}$，分块为

    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 4 \end{pmatrix}, \quad B = \begin{pmatrix} 3 \\ 1 \end{pmatrix}, \quad C = \begin{pmatrix} 1 & 2 \end{pmatrix}, \quad D = (5)
    $$

    $A^{-1} = \frac{1}{8}\begin{pmatrix} 4 & -1 \\ 0 & 2 \end{pmatrix}$。

    $$
    M/A = D - CA^{-1}B = 5 - \begin{pmatrix} 1 & 2 \end{pmatrix} \frac{1}{8}\begin{pmatrix} 4 & -1 \\ 0 & 2 \end{pmatrix} \begin{pmatrix} 3 \\ 1 \end{pmatrix}
    $$

    $$
    = 5 - \begin{pmatrix} 1 & 2 \end{pmatrix} \frac{1}{8} \begin{pmatrix} 11 \\ 2 \end{pmatrix} = 5 - \frac{15}{8} = \frac{25}{8}
    $$

### 广义 Schur 补

!!! definition "定义 34.2 (广义 Schur 补)"
    当 $A$ 不可逆时，可以用广义逆替代 $A^{-1}$。定义**广义 Schur 补**为

    $$
    M/A = D - CA^\dagger B
    $$

    其中 $A^\dagger$ 是 $A$ 的 Moore-Penrose 逆。更一般地，对任何 $\{1\}$-逆 $G \in A\{1\}$：

    $$
    M/_G A = D - CGA B
    $$

---

## 34.2 块高斯消元与分解

<div class="context-flow" markdown>

**核心问题**：Schur 补如何自然地从块消元中产生？

</div>

### 块 LDU 分解

!!! theorem "定理 34.1 (块 LDU 分解)"
    设 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$，$A$ 可逆。则

    $$
    M = \begin{pmatrix} I & 0 \\ CA^{-1} & I \end{pmatrix} \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} \begin{pmatrix} I & A^{-1}B \\ 0 & I \end{pmatrix}
    $$

    这是块高斯消元的矩阵形式。

??? proof "证明"
    直接展开右侧：

    $$
    \begin{pmatrix} I & 0 \\ CA^{-1} & I \end{pmatrix} \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} = \begin{pmatrix} A & 0 \\ C & M/A \end{pmatrix}
    $$

    $$
    \begin{pmatrix} A & 0 \\ C & M/A \end{pmatrix} \begin{pmatrix} I & A^{-1}B \\ 0 & I \end{pmatrix} = \begin{pmatrix} A & B \\ C & CA^{-1}B + M/A \end{pmatrix}
    $$

    由 $M/A = D - CA^{-1}B$，得 $CA^{-1}B + M/A = D$。因此右侧 $= M$。$\blacksquare$

!!! note "注"
    类似地，若 $D$ 可逆：

    $$
    M = \begin{pmatrix} I & BD^{-1} \\ 0 & I \end{pmatrix} \begin{pmatrix} M/D & 0 \\ 0 & D \end{pmatrix} \begin{pmatrix} I & 0 \\ D^{-1}C & I \end{pmatrix}
    $$

### 块求逆公式

!!! theorem "定理 34.2 (块矩阵求逆)"
    设 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$，$A$ 和 $M/A$ 均可逆。则

    $$
    M^{-1} = \begin{pmatrix} A^{-1} + A^{-1}B(M/A)^{-1}CA^{-1} & -A^{-1}B(M/A)^{-1} \\ -(M/A)^{-1}CA^{-1} & (M/A)^{-1} \end{pmatrix}
    $$

??? proof "证明"
    由块 LDU 分解：

    $$
    M^{-1} = \begin{pmatrix} I & -A^{-1}B \\ 0 & I \end{pmatrix} \begin{pmatrix} A^{-1} & 0 \\ 0 & (M/A)^{-1} \end{pmatrix} \begin{pmatrix} I & 0 \\ -CA^{-1} & I \end{pmatrix}
    $$

    展开计算即得结果。$\blacksquare$

!!! example "例 34.2"
    求 $M = \begin{pmatrix} 2 & 1 \\ 3 & 4 \end{pmatrix}$ 的逆。

    取 $A = (2)$，$B = (1)$，$C = (3)$，$D = (4)$。

    $M/A = D - CA^{-1}B = 4 - 3 \cdot \frac{1}{2} \cdot 1 = \frac{5}{2}$。

    $$
    M^{-1} = \begin{pmatrix} \frac{1}{2} + \frac{1}{2} \cdot \frac{2}{5} \cdot 3 \cdot \frac{1}{2} & -\frac{1}{2} \cdot \frac{2}{5} \\ -\frac{2}{5} \cdot 3 \cdot \frac{1}{2} & \frac{2}{5} \end{pmatrix} = \begin{pmatrix} \frac{4}{5} & -\frac{1}{5} \\ -\frac{3}{5} & \frac{2}{5} \end{pmatrix} = \frac{1}{5}\begin{pmatrix} 4 & -1 \\ -3 & 2 \end{pmatrix}
    $$

    验证：$\det(M) = 8 - 3 = 5$，$M^{-1} = \frac{1}{5}\begin{pmatrix} 4 & -1 \\ -3 & 2 \end{pmatrix}$。✓

---

## 34.3 行列式公式

<div class="context-flow" markdown>

**核心问题**：分块矩阵的行列式如何用 Schur 补来表示？

</div>

### Schur 行列式公式

!!! theorem "定理 34.3 (Schur 行列式公式)"
    设 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$。

    1. 若 $A$ 可逆，则 $\det(M) = \det(A) \cdot \det(M/A) = \det(A) \cdot \det(D - CA^{-1}B)$。
    2. 若 $D$ 可逆，则 $\det(M) = \det(D) \cdot \det(M/D) = \det(D) \cdot \det(A - BD^{-1}C)$。

??? proof "证明"
    **(1)** 由块 LDU 分解（定理 34.1）：

    $$
    \det(M) = \det\begin{pmatrix} I & 0 \\ CA^{-1} & I \end{pmatrix} \cdot \det\begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} \cdot \det\begin{pmatrix} I & A^{-1}B \\ 0 & I \end{pmatrix}
    $$

    第一个和第三个因子的行列式都是 1（下三角和上三角矩阵，对角线全为 1）。第二个因子的行列式是 $\det(A) \cdot \det(M/A)$。

    **(2)** 完全类似，用 $D$ 的块 LDU 分解。$\blacksquare$

!!! example "例 34.3"
    计算 $\det \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 0 \end{pmatrix}$。

    取 $A = \begin{pmatrix} 1 & 2 \\ 4 & 5 \end{pmatrix}$，$\det(A) = 5 - 8 = -3$。

    $A^{-1} = -\frac{1}{3}\begin{pmatrix} 5 & -2 \\ -4 & 1 \end{pmatrix}$。

    $M/A = 0 - (7, 8) \cdot \left(-\frac{1}{3}\right)\begin{pmatrix} 5 & -2 \\ -4 & 1 \end{pmatrix}\begin{pmatrix} 3 \\ 6 \end{pmatrix}$

    $= -(7, 8) \cdot \left(-\frac{1}{3}\right)\begin{pmatrix} 3 \\ -6 \end{pmatrix} = -(7, 8) \cdot \begin{pmatrix} -1 \\ 2 \end{pmatrix} = -(-7 + 16) = -9$

    $\det(M) = (-3)(-9) = 27$。

    验证（直接展开）：$1(0-48) - 2(0-42) + 3(32-35) = -48 + 84 - 9 = 27$。✓

### 行列式的乘法性推广

!!! theorem "定理 34.4"
    设 $A \in \mathbb{C}^{p \times q}$，$B \in \mathbb{C}^{q \times p}$。则

    $$
    \det(I_p + AB) = \det(I_q + BA)
    $$

??? proof "证明"
    考虑分块矩阵

    $$
    M = \begin{pmatrix} I_p & A \\ -B & I_q \end{pmatrix}
    $$

    由定理 34.3：$\det(M) = \det(I_p) \cdot \det(I_q - (-B)I_p^{-1}A) = \det(I_q + BA)$。

    也可以用 $D = I_q$：$\det(M) = \det(I_q) \cdot \det(I_p - A I_q^{-1}(-B)) = \det(I_p + AB)$。

    因此 $\det(I_p + AB) = \det(I_q + BA)$。$\blacksquare$

!!! note "注"
    $\det(I + AB) = \det(I + BA)$ 的推论：$AB$ 和 $BA$ 具有相同的非零特征值（含重数）。这是因为 $\det(\lambda I - AB) = \lambda^{p-q} \det(\lambda I - BA)$（假设 $p \geq q$），而上述恒等式是 $\lambda = -1$ 时的特例（适当调整后）。

---

## 34.4 Schur 补与正定性

<div class="context-flow" markdown>

**核心问题**：分块 Hermite 矩阵的正定性如何用 Schur 补来判定？

</div>

### 正定性判据

!!! theorem "定理 34.5 (Schur 补正定性判据)"
    设 $M = \begin{pmatrix} A & B \\ B^* & D \end{pmatrix}$ 是 Hermite 矩阵（$A \in \mathbb{C}^{p \times p}$，$D \in \mathbb{C}^{q \times q}$ 都是 Hermite 的）。则：

    1. $M > 0$（正定）当且仅当 $A > 0$ 且 $M/A = D - B^*A^{-1}B > 0$。
    2. $M \geq 0$（半正定）当且仅当 $A \geq 0$，$\operatorname{col}(B) \subseteq \operatorname{col}(A)$，且 $M/A = D - B^*A^\dagger B \geq 0$。
    3. 等价地：$M > 0 \Leftrightarrow D > 0$ 且 $M/D = A - BD^{-1}B^* > 0$。

??? proof "证明"
    **(1)** **"$\Rightarrow$"**：$M > 0$ 意味着对所有非零 $\boldsymbol{x} \in \mathbb{C}^p$，

    $$
    \boldsymbol{x}^*A\boldsymbol{x} = \begin{pmatrix} \boldsymbol{x} \\ \boldsymbol{0} \end{pmatrix}^* M \begin{pmatrix} \boldsymbol{x} \\ \boldsymbol{0} \end{pmatrix} > 0
    $$

    因此 $A > 0$。

    现在由块 LDU 分解：

    $$
    M = L \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} L^*
    $$

    其中 $L = \begin{pmatrix} I & 0 \\ B^*A^{-1} & I \end{pmatrix}$。由于 $L$ 可逆，$M > 0$ 等价于 $\begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} > 0$，等价于 $A > 0$ 且 $M/A > 0$。

    **"$\Leftarrow$"**：若 $A > 0$ 且 $M/A > 0$，由上述分解 $M = L \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} L^*$，$L$ 可逆，中间矩阵正定，故 $M$ 正定。

    **(2)** 半正定情形更精细。$M \geq 0$ 意味着 $A \geq 0$（同理取 $\boldsymbol{y} = \boldsymbol{0}$）。

    对任何 $\boldsymbol{y} \in \mathbb{C}^q$，取 $\boldsymbol{x} = -A^\dagger B \boldsymbol{y}$：

    $$
    \begin{pmatrix} \boldsymbol{x} \\ \boldsymbol{y} \end{pmatrix}^* M \begin{pmatrix} \boldsymbol{x} \\ \boldsymbol{y} \end{pmatrix} = \boldsymbol{x}^*A\boldsymbol{x} + 2\operatorname{Re}(\boldsymbol{x}^*B\boldsymbol{y}) + \boldsymbol{y}^*D\boldsymbol{y}
    $$

    当 $\operatorname{col}(B) \subseteq \operatorname{col}(A)$ 时（即 $B = AA^\dagger B$），这等于

    $$
    \boldsymbol{y}^*(D - B^*A^\dagger B)\boldsymbol{y} = \boldsymbol{y}^*(M/A)\boldsymbol{y} \geq 0
    $$

    因此 $M/A \geq 0$。$\operatorname{col}(B) \subseteq \operatorname{col}(A)$ 的必要性来自：若 $B\boldsymbol{y} \notin \operatorname{col}(A)$，可以找到 $\boldsymbol{x}$ 使二次型为负。$\blacksquare$

!!! example "例 34.4"
    判断 $M = \begin{pmatrix} 4 & 2 & 1 \\ 2 & 5 & 3 \\ 1 & 3 & 6 \end{pmatrix}$ 的正定性。

    $A = \begin{pmatrix} 4 & 2 \\ 2 & 5 \end{pmatrix}$，$\det(A) = 16 > 0$，$a_{11} = 4 > 0$，故 $A > 0$。

    $A^{-1} = \frac{1}{16}\begin{pmatrix} 5 & -2 \\ -2 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$。

    $M/A = 6 - (1, 3) \frac{1}{16}\begin{pmatrix} 5 & -2 \\ -2 & 4 \end{pmatrix}\begin{pmatrix} 1 \\ 3 \end{pmatrix} = 6 - (1, 3)\frac{1}{16}\begin{pmatrix} -1 \\ 10 \end{pmatrix} = 6 - \frac{29}{16} = \frac{67}{16} > 0$。

    因此 $M > 0$。

### Haynsworth 惯性公式

!!! theorem "定理 34.6 (Haynsworth 惯性公式, 1968)"
    设 $M = \begin{pmatrix} A & B \\ B^* & D \end{pmatrix}$ 是 Hermite 矩阵，$A$ 可逆。则

    $$
    \operatorname{In}(M) = \operatorname{In}(A) + \operatorname{In}(M/A)
    $$

    其中 $\operatorname{In}(X) = (n_+, n_-, n_0)$（正、负、零特征值的个数），加法是分量逐个相加。

??? proof "证明"
    由块 LDU 分解 $M = L \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} L^*$，$L$ 可逆。由 Sylvester 惯性定律，$M$ 与 $\begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix}$ 有相同的惯性。而块对角矩阵的惯性是各块惯性之和。$\blacksquare$

!!! note "注"
    Haynsworth 惯性公式给出了一种递归计算矩阵惯性的方法：先确定 $A$ 的惯性，再计算 Schur 补 $M/A$ 的惯性。这等价于带号高斯消元。

---

## 34.5 Sherman-Morrison-Woodbury 公式

<div class="context-flow" markdown>

**核心问题**：低秩扰动后的逆矩阵如何高效更新？

</div>

### 一般公式

!!! theorem "定理 34.7 (Sherman-Morrison-Woodbury 公式)"
    设 $A \in \mathbb{C}^{n \times n}$ 可逆，$U \in \mathbb{C}^{n \times k}$，$C \in \mathbb{C}^{k \times k}$ 可逆，$V \in \mathbb{C}^{k \times n}$。若 $A + UCV$ 可逆，则

    $$
    (A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^{-1}U)^{-1}VA^{-1}
    $$

??? proof "证明"
    令 $S = C^{-1} + VA^{-1}U$，设 $S$ 可逆。令 $R = A^{-1} - A^{-1}US^{-1}VA^{-1}$。我们直接验证 $(A + UCV)R = I$。

    展开乘积：

    $$
    (A + UCV)R = (A + UCV)(A^{-1} - A^{-1}US^{-1}VA^{-1})
    $$

    $$
    = I - US^{-1}VA^{-1} + UCVA^{-1} - UCVA^{-1}US^{-1}VA^{-1}
    $$

    将后三项合并，提取公因子 $U$ 和 $VA^{-1}$：

    $$
    = I + U\bigl[-S^{-1} + C - C(VA^{-1}U)S^{-1}\bigr]VA^{-1}
    $$

    其中我们使用了 $UCVA^{-1} = U \cdot C \cdot VA^{-1}$ 以及 $UCVA^{-1}US^{-1}VA^{-1} = UC(VA^{-1}U)S^{-1}VA^{-1}$。

    现在分析中括号内的表达式。由 $S = C^{-1} + VA^{-1}U$，有 $CS = I + C(VA^{-1}U)$，即

    $$
    I + C(VA^{-1}U) = CS
    $$

    因此

    $$
    S^{-1} + C(VA^{-1}U)S^{-1} = (I + C(VA^{-1}U))S^{-1} = CS \cdot S^{-1} = C
    $$

    代入可得

    $$
    -S^{-1} + C - C(VA^{-1}U)S^{-1} = C - \bigl[S^{-1} + C(VA^{-1}U)S^{-1}\bigr] = C - C = 0
    $$

    因此 $(A + UCV)R = I$。类似地可以验证 $R(A + UCV) = I$，故 $R = (A + UCV)^{-1}$。$\blacksquare$

!!! note "注"
    Sherman-Morrison-Woodbury 公式也可通过 Schur 补推导。考虑分块矩阵

    $$
    N = \begin{pmatrix} A & U \\ V & -C^{-1} \end{pmatrix}
    $$

    其关于 $-C^{-1}$ 的 Schur 补为 $N/(-C^{-1}) = A - U(-C)V = A + UCV$，关于 $A$ 的 Schur 补为 $N/A = -C^{-1} - VA^{-1}U = -S$。由块矩阵求逆公式（定理 34.2），$N^{-1}$ 的左上块为 $A^{-1} + A^{-1}U(N/A)^{-1}VA^{-1} = A^{-1} - A^{-1}US^{-1}VA^{-1}$，而这恰好等于 $(N/(-C^{-1}))^{-1} = (A+UCV)^{-1}$，从而得到 Woodbury 公式。

### Sherman-Morrison 公式（秩-1 特例）

!!! theorem "定理 34.8 (Sherman-Morrison 公式)"
    设 $A$ 可逆，$\boldsymbol{u}, \boldsymbol{v} \in \mathbb{C}^n$。若 $1 + \boldsymbol{v}^*A^{-1}\boldsymbol{u} \neq 0$，则

    $$
    (A + \boldsymbol{u}\boldsymbol{v}^*)^{-1} = A^{-1} - \frac{A^{-1}\boldsymbol{u}\boldsymbol{v}^*A^{-1}}{1 + \boldsymbol{v}^*A^{-1}\boldsymbol{u}}
    $$

!!! note "注"
    Sherman-Morrison 公式是 Woodbury 公式中 $k = 1$，$C = 1$ 的特例。它的计算复杂度为 $O(n^2)$（两次矩阵-向量乘法加外积），远低于重新计算逆矩阵的 $O(n^3)$。

!!! example "例 34.5"
    设 $A = I_3$，$\boldsymbol{u} = (1, 0, 0)^T$，$\boldsymbol{v} = (0, 1, 0)^T$。

    $A + \boldsymbol{u}\boldsymbol{v}^* = \begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    $1 + \boldsymbol{v}^*A^{-1}\boldsymbol{u} = 1 + 0 = 1$。

    $(A + \boldsymbol{u}\boldsymbol{v}^*)^{-1} = I - \boldsymbol{u}\boldsymbol{v}^* = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。

    直接求逆验证：$\begin{pmatrix} 1 & 1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}^{-1} = \begin{pmatrix} 1 & -1 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix}$。✓

---

## 34.6 Schur 补的商性质与单调性

<div class="context-flow" markdown>

**核心问题**：Schur 补操作在多级分块下如何复合？正定序下 Schur 补如何变化？

</div>

### 商性质

!!! theorem "定理 34.9 (Schur 补的商性质)"
    设 $M \in \mathbb{C}^{n \times n}$ 分块为

    $$
    M = \begin{pmatrix} A & B & E \\ C & D & F \\ G & H & J \end{pmatrix}
    $$

    记 $\tilde{A} = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$（左上 $(p+q) \times (p+q)$ 块），$\tilde{D} = \begin{pmatrix} D & F \\ H & J \end{pmatrix}$（右下 $(q+r) \times (q+r)$ 块），其中 $A \in \mathbb{C}^{p \times p}$，$D \in \mathbb{C}^{q \times q}$，$J \in \mathbb{C}^{r \times r}$。若 $A$ 和 $\tilde{A}$ 均可逆，则

    $$
    (M/A) / (\tilde{A}/A) = M/\tilde{A}
    $$

    即"先对 $A$ 取 Schur 补，再对 $\tilde{A}/A$ 取 Schur 补"等价于"直接对 $\tilde{A}$ 取 Schur 补"。

??? proof "证明"
    由块 LDU 分解，$M$ 关于 $A$ 的 Schur 补为

    $$
    M/A = \begin{pmatrix} D - CA^{-1}B & F - CA^{-1}E \\ H - GA^{-1}B & J - GA^{-1}E \end{pmatrix} = \begin{pmatrix} \tilde{A}/A & F' \\ H' & J' \end{pmatrix}
    $$

    其中 $\tilde{A}/A = D - CA^{-1}B$。由 Schur 行列式公式（对两种分块应用）：

    $$
    \det(M) = \det(A) \cdot \det(M/A) = \det(A) \cdot \det(\tilde{A}/A) \cdot \det((M/A)/(\tilde{A}/A))
    $$

    另一方面

    $$
    \det(M) = \det(\tilde{A}) \cdot \det(M/\tilde{A})
    $$

    由 $\det(\tilde{A}) = \det(A) \cdot \det(\tilde{A}/A)$，比较两式即得 $\det((M/A)/(\tilde{A}/A)) = \det(M/\tilde{A})$。

    为证明矩阵相等（而非仅行列式相等），可直接展开。$M/\tilde{A}$ 是 $J$ 减去 $(G, H)\tilde{A}^{-1}\begin{pmatrix} E \\ F \end{pmatrix}$。利用 $\tilde{A}$ 的块求逆公式将 $\tilde{A}^{-1}$ 展开，与 $(M/A)/(\tilde{A}/A)$ 的定义逐项比对即可验证两者相等。$\blacksquare$

### Schur 补单调性

!!! theorem "定理 34.10 (Schur 补的单调性)"
    设

    $$
    M_1 = \begin{pmatrix} A_1 & B \\ B^* & D \end{pmatrix} \geq 0, \quad M_2 = \begin{pmatrix} A_2 & B \\ B^* & D \end{pmatrix} \geq 0
    $$

    是 Hermite 半正定矩阵，且 $A_1, A_2$ 均可逆。若 $A_1 \geq A_2 > 0$（即 $A_1 - A_2 \geq 0$ 且 $A_2 > 0$），则

    $$
    M_1/A_1 \leq M_2/A_2
    $$

    即 $D - B^*A_1^{-1}B \leq D - B^*A_2^{-1}B$，等价于 $A_1^{-1} \leq A_2^{-1}$（此为逆矩阵的单调性）。

??? proof "证明"
    由 $A_1 \geq A_2 > 0$ 和逆矩阵的单调反转性质（若 $X \geq Y > 0$ 则 $X^{-1} \leq Y^{-1}$），得 $A_1^{-1} \leq A_2^{-1}$。因此

    $$
    B^*A_1^{-1}B \leq B^*A_2^{-1}B
    $$

    （因为对任何 $\boldsymbol{x}$，$\boldsymbol{x}^*B^*A_1^{-1}B\boldsymbol{x} = (B\boldsymbol{x})^*A_1^{-1}(B\boldsymbol{x}) \leq (B\boldsymbol{x})^*A_2^{-1}(B\boldsymbol{x}) = \boldsymbol{x}^*B^*A_2^{-1}B\boldsymbol{x}$）。

    从而

    $$
    M_1/A_1 = D - B^*A_1^{-1}B \geq D - B^*A_2^{-1}B = M_2/A_2
    $$

    即 $M_2/A_2 \leq M_1/A_1$。

    注意结论的方向：$A$ 块增大时，Schur 补减小，这反映了"消去更多"会留下更小的余量。$\blacksquare$

!!! note "注"
    Schur 补的单调性在优化和统计中有重要应用。例如在最优实验设计中，信息矩阵的 Schur 补给出了参数的条件信息，单调性保证了增加观测信息不会降低参数估计的精度。

---

## 34.7 Schur 补在统计学中的应用

<div class="context-flow" markdown>

**核心问题**：条件分布的协方差矩阵与 Schur 补有什么关系？

</div>

### 多元正态分布的条件分布

!!! theorem "定理 34.11 (条件正态分布)"
    设 $\boldsymbol{X} = \begin{pmatrix} \boldsymbol{X}_1 \\ \boldsymbol{X}_2 \end{pmatrix} \sim \mathcal{N}\left(\begin{pmatrix} \boldsymbol{\mu}_1 \\ \boldsymbol{\mu}_2 \end{pmatrix}, \begin{pmatrix} \Sigma_{11} & \Sigma_{12} \\ \Sigma_{21} & \Sigma_{22} \end{pmatrix}\right)$

    是多元正态分布，其中 $\Sigma_{22} > 0$。则 $\boldsymbol{X}_1 | \boldsymbol{X}_2 = \boldsymbol{x}_2$ 服从正态分布：

    $$
    \boldsymbol{X}_1 | \boldsymbol{X}_2 = \boldsymbol{x}_2 \sim \mathcal{N}(\boldsymbol{\mu}_{1|2}, \Sigma_{1|2})
    $$

    其中

    $$
    \boldsymbol{\mu}_{1|2} = \boldsymbol{\mu}_1 + \Sigma_{12}\Sigma_{22}^{-1}(\boldsymbol{x}_2 - \boldsymbol{\mu}_2)
    $$

    $$
    \Sigma_{1|2} = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21} = \Sigma / \Sigma_{22}
    $$

    即**条件协方差矩阵恰好是 Schur 补**。

??? proof "证明"
    设 $\Sigma > 0$。联合密度函数为

    $$
    f(\boldsymbol{x}) \propto \exp\left(-\frac{1}{2}(\boldsymbol{x} - \boldsymbol{\mu})^T \Sigma^{-1} (\boldsymbol{x} - \boldsymbol{\mu})\right)
    $$

    由块矩阵求逆公式（定理 34.2），$\Sigma^{-1}$ 的左上块为 $(\Sigma/\Sigma_{22})^{-1} = \Sigma_{1|2}^{-1}$。

    将指数中的二次型展开并配方。设 $\boldsymbol{z}_1 = \boldsymbol{x}_1 - \boldsymbol{\mu}_1$，$\boldsymbol{z}_2 = \boldsymbol{x}_2 - \boldsymbol{\mu}_2$。

    $$
    Q = (\boldsymbol{x} - \boldsymbol{\mu})^T \Sigma^{-1} (\boldsymbol{x} - \boldsymbol{\mu})
    $$

    利用块 LDU 分解 $\Sigma^{-1}$：

    $$
    Q = (\boldsymbol{z}_1 - \Sigma_{12}\Sigma_{22}^{-1}\boldsymbol{z}_2)^T \Sigma_{1|2}^{-1} (\boldsymbol{z}_1 - \Sigma_{12}\Sigma_{22}^{-1}\boldsymbol{z}_2) + \boldsymbol{z}_2^T \Sigma_{22}^{-1} \boldsymbol{z}_2
    $$

    给定 $\boldsymbol{X}_2 = \boldsymbol{x}_2$ 后，第二项是常数，第一项给出条件分布：

    $$
    \boldsymbol{X}_1 | \boldsymbol{X}_2 \sim \mathcal{N}(\boldsymbol{\mu}_1 + \Sigma_{12}\Sigma_{22}^{-1}\boldsymbol{z}_2, \Sigma_{1|2})
    $$

    $\blacksquare$

!!! note "注"
    注意条件协方差 $\Sigma_{1|2} = \Sigma_{11} - \Sigma_{12}\Sigma_{22}^{-1}\Sigma_{21}$ **不依赖于** $\boldsymbol{x}_2$ 的值。这是正态分布的特殊性质：给定条件后，方差不变，只有均值移动。

### 偏相关系数

!!! definition "定义 34.3 (偏相关系数)"
    变量 $X_i$ 和 $X_j$ 在给定其余变量条件下的**偏相关系数**为

    $$
    \rho_{ij \cdot \text{rest}} = -\frac{(\Sigma^{-1})_{ij}}{\sqrt{(\Sigma^{-1})_{ii} \cdot (\Sigma^{-1})_{jj}}}
    $$

    这与 Schur 补的关系：精度矩阵 $\Sigma^{-1}$ 的元素可以通过递归 Schur 补来计算。

!!! example "例 34.6"
    设 $\boldsymbol{X} = (X_1, X_2, X_3)^T$ 的协方差矩阵为

    $$
    \Sigma = \begin{pmatrix} 1 & 0.5 & 0.3 \\ 0.5 & 1 & 0.4 \\ 0.3 & 0.4 & 1 \end{pmatrix}
    $$

    $X_1 | X_3 = x_3$ 的条件方差：

    取 $\Sigma_{11} = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix}$，$\Sigma_{22} = (1)$，$\Sigma_{12} = \begin{pmatrix} 0.3 \\ 0.4 \end{pmatrix}$。

    $$
    \Sigma_{1|2} = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix} - \begin{pmatrix} 0.3 \\ 0.4 \end{pmatrix} (1)^{-1} \begin{pmatrix} 0.3 & 0.4 \end{pmatrix} = \begin{pmatrix} 0.91 & 0.38 \\ 0.38 & 0.84 \end{pmatrix}
    $$

    因此 $\operatorname{Var}(X_1 | X_3) = 0.91$，$\operatorname{Var}(X_2 | X_3) = 0.84$，$\operatorname{Cov}(X_1, X_2 | X_3) = 0.38$。

---

## 34.8 Schur 补在优化中的应用

<div class="context-flow" markdown>

**核心问题**：如何利用 Schur 补将非线性约束转化为线性矩阵不等式？

</div>

### Schur 补引理

!!! theorem "定理 34.12 (Schur 补引理 / S-过程)"
    以下非线性矩阵不等式

    $$
    A - BD^{-1}B^* \geq 0, \quad D > 0
    $$

    等价于线性矩阵不等式（LMI）

    $$
    \begin{pmatrix} A & B \\ B^* & D \end{pmatrix} \geq 0
    $$

    这将非线性约束（含 $D^{-1}$）转化为线性约束（矩阵半正定），是半定规划（SDP）的核心工具。

!!! example "例 34.7"
    **最小方差无偏估计**：在信号处理中，约束

    $$
    \min_{\boldsymbol{w}} \boldsymbol{w}^* R \boldsymbol{w} \quad \text{s.t.} \quad \boldsymbol{a}^* \boldsymbol{w} = 1
    $$

    （$R$ 是协方差矩阵，$\boldsymbol{a}$ 是导向向量）的对偶问题可以用 Schur 补写为 SDP：

    $$
    \max_{t} \quad t \quad \text{s.t.} \quad \begin{pmatrix} R & \boldsymbol{a} \\ \boldsymbol{a}^* & t^{-1} \end{pmatrix} \geq 0
    $$

    进一步通过变量替换 $s = t^{-1}$ 线性化。

### 控制论中的应用

!!! theorem "定理 34.13 (有界实引理 / Bounded Real Lemma)"
    线性系统 $\dot{\boldsymbol{x}} = A\boldsymbol{x} + B\boldsymbol{u}$，$\boldsymbol{y} = C\boldsymbol{x} + D\boldsymbol{u}$ 的 $H_\infty$ 范数不超过 $\gamma$（即 $\|G\|_\infty < \gamma$）当且仅当存在 $P > 0$ 使得

    $$
    \begin{pmatrix} A^TP + PA + C^TC & PB + C^TD \\ B^TP + D^TC & D^TD - \gamma^2 I \end{pmatrix} < 0
    $$

    这是一个 LMI，可以高效地用 SDP 求解。Schur 补将这个 LMI 与 Riccati 不等式联系起来。

!!! note "注"
    Schur 补在现代优化和控制论中的核心地位源于一个简单但深刻的事实：它将**矩阵的逆运算**（非线性操作）转化为**矩阵不等式**（凸约束）。这使得原本难以处理的非凸优化问题变成了凸优化问题。

    本章的 Schur 补理论将分块矩阵的代数结构与正定性、行列式、条件分布等概念统一在一个框架中。读者应特别记住三个核心结果：行列式公式（定理 34.3）、正定性判据（定理 34.5）和 Sherman-Morrison-Woodbury 公式（定理 34.7）。它们是矩阵分析中使用频率最高的工具。
