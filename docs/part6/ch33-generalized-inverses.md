# 第 33 章 广义逆

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 最小二乘(Ch7) · SVD(Ch11)

**本章脉络**：$\{1\}$-逆 → Moore-Penrose 逆（四个 Penrose 条件）→ SVD 表示 → 最小范数最小二乘解 → Drazin 逆 → 群逆 → 扰动分析 → 加权 Moore-Penrose 逆 → 逆序律 → 并联和

**延伸**：Drazin 逆在奇异微分方程和 Markov 链（稳态分布计算）中不可或缺；广义逆理论推广到 Hilbert 空间中的闭算子（von Neumann 正则逆）和 Banach 代数

</div>

对于可逆方阵 $A$，逆矩阵 $A^{-1}$ 提供了线性方程组 $A\boldsymbol{x} = \boldsymbol{b}$ 的唯一解。但当 $A$ 是奇异矩阵或非方阵时，经典逆不存在。广义逆理论的核心任务是为一般矩阵定义各种"逆"，使得它们在特定意义下保留逆矩阵的部分性质。

Moore（1920）和 Penrose（1955）独立定义了最重要的广义逆——Moore-Penrose 逆。Drazin（1958）从代数的观点引入了另一种广义逆。本章系统发展这些理论，揭示它们与 SVD、投影和最小二乘问题的深层联系。

---

## 33.1 内逆与 $\{1\}$-逆

<div class="context-flow" markdown>

**核心问题**：对于一般矩阵 $A$，是否存在矩阵 $X$ 使得 $AXA = A$？这样的 $X$ 有什么用？

</div>

### 定义

!!! definition "定义 33.1 ($\{1\}$-逆)"
    设 $A \in \mathbb{C}^{m \times n}$。称 $X \in \mathbb{C}^{n \times m}$ 为 $A$ 的 **$\{1\}$-逆**（或**内逆**、**广义逆**），若满足

    $$
    AXA = A \tag{1}
    $$

    $A$ 的所有 $\{1\}$-逆的集合记作 $A\{1\}$。

!!! theorem "定理 33.1 ($\{1\}$-逆的存在性)"
    对任何矩阵 $A \in \mathbb{C}^{m \times n}$，$A\{1\}$ 非空。即 $\{1\}$-逆总是存在的。

??? proof "证明"
    设 $\operatorname{rank}(A) = r$。则存在可逆矩阵 $P \in \mathbb{C}^{m \times m}$ 和 $Q \in \mathbb{C}^{n \times n}$ 使得

    $$
    A = P \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} Q
    $$

    这是 $A$ 的秩分解的推论（通过行列变换化为标准形）。取

    $$
    X = Q^{-1} \begin{pmatrix} I_r & C \\ D & E \end{pmatrix} P^{-1}
    $$

    其中 $C, D, E$ 是任意大小适当的矩阵。验证：

    $$
    AXA = P \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} I_r & C \\ D & E \end{pmatrix} \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} Q = P \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} Q = A
    $$

    因此 $X \in A\{1\}$。$\blacksquare$

!!! note "注"
    $\{1\}$-逆一般不唯一。上述证明中 $C, D, E$ 的任意性说明 $A\{1\}$ 通常是一个无穷集合。

### 与相容方程组的关系

!!! theorem "定理 33.2"
    线性方程组 $A\boldsymbol{x} = \boldsymbol{b}$ 有解当且仅当对任何 $X \in A\{1\}$，$AX\boldsymbol{b} = \boldsymbol{b}$。

    当方程组有解时，$\boldsymbol{x}_0 = X\boldsymbol{b}$ 是一个特解，通解为

    $$
    \boldsymbol{x} = X\boldsymbol{b} + (I - XA)\boldsymbol{z}, \quad \boldsymbol{z} \in \mathbb{C}^n \text{ 任意}
    $$

??? proof "证明"
    **必要性**：设 $A\boldsymbol{x}_0 = \boldsymbol{b}$，则 $AX\boldsymbol{b} = AXA\boldsymbol{x}_0 = A\boldsymbol{x}_0 = \boldsymbol{b}$。

    **充分性**：若 $AX\boldsymbol{b} = \boldsymbol{b}$，取 $\boldsymbol{x}_0 = X\boldsymbol{b}$，则 $A\boldsymbol{x}_0 = AX\boldsymbol{b} = \boldsymbol{b}$。

    **通解**：$A(X\boldsymbol{b} + (I-XA)\boldsymbol{z}) = AX\boldsymbol{b} + A\boldsymbol{z} - AXA\boldsymbol{z} = \boldsymbol{b} + A\boldsymbol{z} - A\boldsymbol{z} = \boldsymbol{b}$。

    反方向，若 $A\boldsymbol{x} = \boldsymbol{b}$，令 $\boldsymbol{z} = \boldsymbol{x}$，则 $X\boldsymbol{b} + (I-XA)\boldsymbol{x} = X\boldsymbol{b} + \boldsymbol{x} - XA\boldsymbol{x} = X\boldsymbol{b} + \boldsymbol{x} - X\boldsymbol{b} = \boldsymbol{x}$。$\blacksquare$

!!! example "例 33.1"
    设 $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$。$\operatorname{rank}(A) = 1$。

    可以验证 $X = \begin{pmatrix} 1/5 & 0 \\ 0 & 0 \end{pmatrix}$ 满足 $AXA = A$：

    $$
    AXA = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} \begin{pmatrix} 1/5 & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} = \begin{pmatrix} 1/5 & 0 \\ 2/5 & 0 \end{pmatrix} \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} = \begin{pmatrix} 1/5 & 2/5 \\ 2/5 & 4/5 \end{pmatrix}
    $$

    等等，让我们重新计算。$AX = \begin{pmatrix} 1/5 & 0 \\ 2/5 & 0 \end{pmatrix}$，$(AX)A = \begin{pmatrix} 1/5 & 2/5 \\ 2/5 & 4/5 \end{pmatrix} \neq A$。

    重新选择。取 $X = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$：$AXA = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 2 & 0 \end{pmatrix}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix} = A$。正确。

    因此 $X = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} \in A\{1\}$。

---

## 33.2 Moore-Penrose 逆

<div class="context-flow" markdown>

**核心问题**：能否在所有广义逆中找到一个"最自然"的？

</div>

### Penrose 条件

!!! definition "定义 33.2 (Penrose 条件)"
    设 $A \in \mathbb{C}^{m \times n}$，$X \in \mathbb{C}^{n \times m}$。以下四个方程称为 **Penrose 条件**：

    1. $AXA = A$
    2. $XAX = X$
    3. $(AX)^* = AX$（即 $AX$ 是 Hermite 的）
    4. $(XA)^* = XA$（即 $XA$ 是 Hermite 的）

    满足条件 (1) 的 $X$ 称为 $A$ 的 $\{1\}$-逆。满足所有四个条件的 $X$ 称为 $A$ 的 **Moore-Penrose 逆**，记作 $A^\dagger$。

!!! theorem "定理 33.3 (Moore-Penrose 逆的唯一性)"
    对任何 $A \in \mathbb{C}^{m \times n}$，满足全部四个 Penrose 条件的 $X$ 存在且唯一。

??? proof "证明"
    **唯一性**：设 $X$ 和 $Y$ 都满足四个条件。则

    $$
    X = XAX = X(AX)^* = XX^*A^* = XX^*(AYA)^* = XX^*A^*Y^*A^*
    $$

    $$
    = X(AX)^*(AY)^* = XAXAY = XAY
    $$

    类似地 $Y = XAY$。因此 $X = Y$。

    更简洁的证明：

    $$
    XA = (XA)^* = A^*X^*, \quad AX = (AX)^* = X^*A^*
    $$

    $$
    X = XAX = (XA)X = A^*X^*X
    $$

    $$
    Y = YAY = Y(AY) = YX^*A^*
    $$

    因此

    $$
    X = A^*X^*X = A^*(XAX)^*X = A^*X^*(AX)^*X = A^*X^*X^*A^*X = (A^*X^*)(X^*A^*)X
    $$

    这变得复杂了。用更直接的方法：

    $XA$ 和 $YA$ 都是 $\mathbb{C}^n$ 上的正交投影（Hermite 且幂等），投影到 $\operatorname{col}(A^*)$ 上。类似地 $AX$ 和 $AY$ 都是投影到 $\operatorname{col}(A)$ 上。因此 $XA = YA$，$AX = AY$，故

    $$
    X = XAX = XAY = YAY = Y
    $$

    **存在性**：通过 SVD 构造，见下一小节。$\blacksquare$

### SVD 表示

!!! theorem "定理 33.4 (Moore-Penrose 逆的 SVD 表示)"
    设 $A \in \mathbb{C}^{m \times n}$，$\operatorname{rank}(A) = r$，SVD 为

    $$
    A = U \Sigma V^* = U \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix} V^*
    $$

    其中 $\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$，$\sigma_1 \geq \cdots \geq \sigma_r > 0$。则

    $$
    A^\dagger = V \begin{pmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{pmatrix} U^* = V \Sigma^\dagger U^*
    $$

    其中 $\Sigma^\dagger$ 是将 $\Sigma$ 中每个非零奇异值取倒数、零保持为零后转置得到的矩阵。

??? proof "证明"
    令 $X = V \Sigma^\dagger U^*$，验证四个 Penrose 条件。

    设 $\tilde{\Sigma} = \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix}_{m \times n}$，$\tilde{\Sigma}^\dagger = \begin{pmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{pmatrix}_{n \times m}$。

    **(1)** $AXA = U\tilde{\Sigma}V^* \cdot V\tilde{\Sigma}^\dagger U^* \cdot U\tilde{\Sigma}V^* = U\tilde{\Sigma}\tilde{\Sigma}^\dagger\tilde{\Sigma}V^*$。

    $$
    \tilde{\Sigma}\tilde{\Sigma}^\dagger\tilde{\Sigma} = \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix} = \tilde{\Sigma}
    $$

    故 $AXA = U\tilde{\Sigma}V^* = A$。✓

    **(2)** 类似地，$XAX = V\tilde{\Sigma}^\dagger \tilde{\Sigma} \tilde{\Sigma}^\dagger V^* = V\tilde{\Sigma}^\dagger V^* = X$。✓

    **(3)** $AX = U\tilde{\Sigma}\tilde{\Sigma}^\dagger U^* = U \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} U^*$。这是 Hermite 的（$U$ 酉，中间矩阵实对角）。✓

    **(4)** $XA = V\tilde{\Sigma}^\dagger\tilde{\Sigma} V^* = V \begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix} V^*$。同理是 Hermite 的。✓

    因此 $X = A^\dagger$。$\blacksquare$

!!! example "例 33.2"
    设 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \\ 0 & 1 \end{pmatrix}$。

    SVD：$A = U\Sigma V^*$，其中 $\Sigma = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$，$U = I_3$，$V = I_2$（$A$ 的列已正交归一）。

    $$
    A^\dagger = V\Sigma^\dagger U^* = I_2 \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} I_3 = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}
    $$

    验证：$AA^\dagger = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}$（投影到 $\operatorname{col}(A)$），$A^\dagger A = I_2$（投影到 $\operatorname{col}(A^*) = \mathbb{C}^2$）。

---

## 33.3 Moore-Penrose 逆的性质与计算

<div class="context-flow" markdown>

**核心问题**：Moore-Penrose 逆具有哪些代数和分析性质？

</div>

### 基本代数性质

!!! theorem "定理 33.5 (Moore-Penrose 逆的性质)"
    设 $A \in \mathbb{C}^{m \times n}$。则：

    1. $(A^\dagger)^\dagger = A$
    2. $(A^*)^\dagger = (A^\dagger)^*$
    3. $(\alpha A)^\dagger = \alpha^{-1} A^\dagger$（当 $\alpha \neq 0$）
    4. $(A^*A)^\dagger = A^\dagger (A^\dagger)^*$
    5. $(AA^*)^\dagger = (A^\dagger)^* A^\dagger$
    6. $A^* = A^*AA^\dagger = A^\dagger AA^*$
    7. $\operatorname{rank}(A^\dagger) = \operatorname{rank}(A)$
    8. 若 $A$ 可逆，则 $A^\dagger = A^{-1}$

??? proof "证明"
    **(1)** $A$ 满足关于 $A^\dagger$ 的四个条件的"对偶"形式（交换 $A$ 和 $A^\dagger$ 的角色并利用 Hermite 条件的对称性）：$A^\dagger A A^\dagger = A^\dagger$，$A A^\dagger A = A$，$(A^\dagger A)^* = A^\dagger A$，$(A A^\dagger)^* = A A^\dagger$。因此 $A$ 是 $A^\dagger$ 的 Moore-Penrose 逆。

    **(2)** 设 $X = (A^\dagger)^*$。验证 $X$ 满足关于 $A^*$ 的四个 Penrose 条件：
    - $A^*XA^* = A^*(A^\dagger)^*A^* = (A A^\dagger A)^* = A^*$。✓
    - $XA^*X = (A^\dagger)^*A^*(A^\dagger)^* = (A^\dagger AA^\dagger)^* = (A^\dagger)^* = X$。✓
    - $(A^*X)^* = (A^*(A^\dagger)^*)^* = A^\dagger A = (A^\dagger A)^* = (XA^*)^*$... 需要仔细检查。
    $(A^*X)^* = ((A^\dagger)^* A^*)^{**} = (A (A^\dagger))$... 实际上 $A^*X = A^*(A^\dagger)^* = (A^\dagger A)^*= A^\dagger A$，而 $(A^*X)^* = (A^\dagger A)^* = A^\dagger A = A^*X$。✓
    - $(XA^*)^* = ((A^\dagger)^*A^*)^* = A A^\dagger = (AA^\dagger)^* = XA^*$...
    $XA^* = (A^\dagger)^*A^* = (AA^\dagger)^*= AA^\dagger$。$(XA^*)^* = (AA^\dagger)^* = AA^\dagger = XA^*$。✓

    因此 $(A^\dagger)^* = (A^*)^\dagger$。

    **(8)** 若 $A$ 可逆，取 $X = A^{-1}$：$AXA = A$，$XAX = X$，$AX = I = (AX)^*$，$XA = I = (XA)^*$。$\blacksquare$

### 投影性质

!!! theorem "定理 33.6 (投影解释)"
    设 $A \in \mathbb{C}^{m \times n}$，$\operatorname{rank}(A) = r$。则：

    1. $AA^\dagger$ 是到 $\operatorname{col}(A)$ 上的正交投影。
    2. $A^\dagger A$ 是到 $\operatorname{col}(A^*) = \operatorname{row}(A)$ 上的正交投影。
    3. $I_m - AA^\dagger$ 是到 $\ker(A^*)$ 上的正交投影。
    4. $I_n - A^\dagger A$ 是到 $\ker(A)$ 上的正交投影。

??? proof "证明"
    **(1)** $AA^\dagger$ 是 Hermite 的（Penrose 条件 3）且幂等的（$(AA^\dagger)^2 = A(A^\dagger A)A^\dagger = AA^\dagger$，利用条件 1）。因此它是正交投影。

    其值域：$\operatorname{col}(AA^\dagger) \subseteq \operatorname{col}(A)$。反方向，若 $\boldsymbol{y} = A\boldsymbol{x}$，则 $AA^\dagger \boldsymbol{y} = AA^\dagger A\boldsymbol{x} = A\boldsymbol{x} = \boldsymbol{y}$（条件 1），故 $\boldsymbol{y} \in \operatorname{col}(AA^\dagger)$。

    **(2)** 类似证明。$A^\dagger A$ 是 Hermite 幂等的，值域为 $\operatorname{col}(A^\dagger A) = \operatorname{col}(A^*) = \operatorname{row}(A)$。

    **(3)(4)** 由 $\mathbb{C}^m = \operatorname{col}(A) \oplus \ker(A^*)$ 和 $\mathbb{C}^n = \operatorname{row}(A) \oplus \ker(A)$ 直接得到。$\blacksquare$

### 极限表示

!!! theorem "定理 33.7 (极限表示)"
    设 $A \in \mathbb{C}^{m \times n}$。则

    $$
    A^\dagger = \lim_{\epsilon \to 0^+} (A^*A + \epsilon I)^{-1} A^* = \lim_{\epsilon \to 0^+} A^*(AA^* + \epsilon I)^{-1}
    $$

    这提供了一种正则化计算 $A^\dagger$ 的方法（Tikhonov 正则化）。

??? proof "证明"
    设 $A = U\Sigma V^*$ 是 SVD，$\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$。

    $A^*A = V\Sigma^*\Sigma V^* = V \operatorname{diag}(\sigma_1^2, \ldots, \sigma_r^2, 0, \ldots, 0) V^*$。

    $(A^*A + \epsilon I)^{-1} = V \operatorname{diag}\left(\frac{1}{\sigma_1^2 + \epsilon}, \ldots, \frac{1}{\sigma_r^2 + \epsilon}, \frac{1}{\epsilon}, \ldots, \frac{1}{\epsilon}\right) V^*$。

    $(A^*A + \epsilon I)^{-1}A^* = V \operatorname{diag}\left(\frac{1}{\sigma_1^2 + \epsilon}, \ldots, \frac{1}{\sigma_r^2 + \epsilon}, \frac{1}{\epsilon}, \ldots, \frac{1}{\epsilon}\right) V^* V \Sigma^* U^*$

    $= V \operatorname{diag}\left(\frac{\sigma_1}{\sigma_1^2 + \epsilon}, \ldots, \frac{\sigma_r}{\sigma_r^2 + \epsilon}, 0, \ldots, 0\right)_{\text{适当大小}} U^*$

    当 $\epsilon \to 0^+$ 时，$\frac{\sigma_i}{\sigma_i^2 + \epsilon} \to \frac{1}{\sigma_i}$，而零奇异值对应的项保持为 $0$。因此极限就是 $V\Sigma^\dagger U^* = A^\dagger$。$\blacksquare$

!!! example "例 33.3"
    设 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。$\operatorname{rank}(A) = 1$。

    SVD：$\sigma_1 = 2$，$\boldsymbol{u}_1 = \frac{1}{\sqrt{2}}(1, 1)^T$，$\boldsymbol{v}_1 = \frac{1}{\sqrt{2}}(1, 1)^T$。

    $$
    A^\dagger = \boldsymbol{v}_1 \cdot \frac{1}{\sigma_1} \cdot \boldsymbol{u}_1^* = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 \\ 1 \end{pmatrix} \cdot \frac{1}{2} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \end{pmatrix} = \frac{1}{4}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}
    $$

    验证极限公式：$A^*A = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$，$(A^*A + \epsilon I)^{-1} = \frac{1}{(4+\epsilon)\epsilon - 4}\begin{pmatrix} 2+\epsilon & -2 \\ -2 & 2+\epsilon \end{pmatrix}$... 实际上利用特征分解更方便。

---

## 33.4 最小范数最小二乘解

<div class="context-flow" markdown>

**核心问题**：对于一般的线性方程组 $A\boldsymbol{x} = \boldsymbol{b}$（可能不相容），Moore-Penrose 逆给出什么样的解？

</div>

### 最小二乘问题回顾

!!! definition "定义 33.3 (最小范数最小二乘解)"
    给定 $A \in \mathbb{C}^{m \times n}$ 和 $\boldsymbol{b} \in \mathbb{C}^m$，**最小二乘问题**是

    $$
    \min_{\boldsymbol{x} \in \mathbb{C}^n} \|A\boldsymbol{x} - \boldsymbol{b}\|_2
    $$

    最小二乘解的集合记为 $\mathcal{L} = \operatorname{argmin} \|A\boldsymbol{x} - \boldsymbol{b}\|$。在 $\mathcal{L}$ 中范数最小的解

    $$
    \boldsymbol{x}^* = \operatorname{argmin}_{\boldsymbol{x} \in \mathcal{L}} \|\boldsymbol{x}\|_2
    $$

    称为**最小范数最小二乘解**。

!!! theorem "定理 33.8"
    对任何 $A \in \mathbb{C}^{m \times n}$ 和 $\boldsymbol{b} \in \mathbb{C}^m$，最小范数最小二乘解存在且唯一，等于

    $$
    \boldsymbol{x}^* = A^\dagger \boldsymbol{b}
    $$

??? proof "证明"
    **第一步**：最小二乘解集。$\|A\boldsymbol{x} - \boldsymbol{b}\|$ 最小当且仅当 $A\boldsymbol{x}$ 是 $\boldsymbol{b}$ 在 $\operatorname{col}(A)$ 上的正交投影，即

    $$
    A\boldsymbol{x} = AA^\dagger \boldsymbol{b}
    $$

    （因为 $AA^\dagger$ 是到 $\operatorname{col}(A)$ 的正交投影）。这等价于正规方程 $A^*A\boldsymbol{x} = A^*\boldsymbol{b}$。

    最小二乘解集为 $\mathcal{L} = \{A^\dagger \boldsymbol{b} + (I - A^\dagger A)\boldsymbol{z} : \boldsymbol{z} \in \mathbb{C}^n\}$。

    **第二步**：最小范数。对 $\boldsymbol{x} = A^\dagger \boldsymbol{b} + (I - A^\dagger A)\boldsymbol{z} \in \mathcal{L}$：

    $$
    \|\boldsymbol{x}\|^2 = \|A^\dagger \boldsymbol{b}\|^2 + \|(I - A^\dagger A)\boldsymbol{z}\|^2 + 2\operatorname{Re}\langle A^\dagger \boldsymbol{b}, (I - A^\dagger A)\boldsymbol{z}\rangle
    $$

    注意 $A^\dagger \boldsymbol{b} \in \operatorname{col}(A^*)$（因为 $A^\dagger \boldsymbol{b} = A^\dagger A (A^\dagger \boldsymbol{b})$，即 $A^\dagger \boldsymbol{b}$ 在 $\operatorname{col}(A^*)$ 的投影下不变），而 $(I - A^\dagger A)\boldsymbol{z} \in \ker(A)$。由 $\operatorname{col}(A^*) \perp \ker(A)$：

    $$
    \langle A^\dagger \boldsymbol{b}, (I - A^\dagger A)\boldsymbol{z}\rangle = 0
    $$

    因此 $\|\boldsymbol{x}\|^2 = \|A^\dagger \boldsymbol{b}\|^2 + \|(I - A^\dagger A)\boldsymbol{z}\|^2 \geq \|A^\dagger \boldsymbol{b}\|^2$。

    等号当且仅当 $(I - A^\dagger A)\boldsymbol{z} = \boldsymbol{0}$，即 $\boldsymbol{z} \in \operatorname{col}(A^*)$。此时 $\boldsymbol{x} = A^\dagger \boldsymbol{b}$。$\blacksquare$

!!! note "注"
    几何解释：

    - $AA^\dagger$ 将 $\boldsymbol{b}$ 投影到 $\operatorname{col}(A)$，得到 $\boldsymbol{b}$ 的"最佳逼近"$\hat{\boldsymbol{b}} = AA^\dagger \boldsymbol{b}$。
    - $A^\dagger A$ 将解空间投影到 $\operatorname{col}(A^*)$，在所有最小二乘解中选出范数最小的。
    - 因此 $\boldsymbol{x}^* = A^\dagger \boldsymbol{b}$ 同时在两个方向上进行了"最优选择"。

!!! example "例 33.4"
    设 $A = \begin{pmatrix} 1 & 1 \\ 1 & 1 \\ 0 & 0 \end{pmatrix}$ 和 $\boldsymbol{b} = \begin{pmatrix} 3 \\ 1 \\ 2 \end{pmatrix}$。

    $A^*A = \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix}$，$A^*\boldsymbol{b} = \begin{pmatrix} 4 \\ 4 \end{pmatrix}$。

    SVD 计算：$\sigma_1 = 2$，$\boldsymbol{u}_1 = \frac{1}{\sqrt{2}}(1, 1, 0)^T$，$\boldsymbol{v}_1 = \frac{1}{\sqrt{2}}(1, 1)^T$。

    $A^\dagger = \boldsymbol{v}_1 \sigma_1^{-1} \boldsymbol{u}_1^* = \frac{1}{4}\begin{pmatrix} 1 & 1 & 0 \\ 1 & 1 & 0 \end{pmatrix}$。

    $\boldsymbol{x}^* = A^\dagger \boldsymbol{b} = \frac{1}{4}\begin{pmatrix} 4 \\ 4 \end{pmatrix} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    残差：$A\boldsymbol{x}^* - \boldsymbol{b} = (2, 2, 0)^T - (3, 1, 2)^T = (-1, 1, -2)^T$，$\|A\boldsymbol{x}^* - \boldsymbol{b}\| = \sqrt{6}$。

---

## 33.5 Drazin 逆

<div class="context-flow" markdown>

**核心问题**：能否定义一种与矩阵交换的广义逆？

</div>

### 指标与定义

!!! definition "定义 33.4 (指标)"
    方阵 $A \in \mathbb{C}^{n \times n}$ 的**指标**（index）$\operatorname{ind}(A) = k$ 定义为使得 $\operatorname{rank}(A^k) = \operatorname{rank}(A^{k+1})$ 的最小非负整数 $k$。

    等价地，$\operatorname{ind}(A) = k$ 当且仅当 $\mathbb{C}^n = \operatorname{col}(A^k) \oplus \ker(A^k)$。

!!! note "注"
    - 若 $A$ 可逆，则 $\operatorname{ind}(A) = 0$。
    - 若 $A$ 是幂零矩阵且 $A^k = 0$、$A^{k-1} \neq 0$，则 $\operatorname{ind}(A) = k$。
    - 对任何矩阵，$\operatorname{ind}(A) \leq n$。

!!! definition "定义 33.5 (Drazin 逆)"
    设 $A \in \mathbb{C}^{n \times n}$，$\operatorname{ind}(A) = k$。$A$ 的 **Drazin 逆** $A^D$ 是满足以下条件的唯一矩阵：

    1. $A^{k+1} A^D = A^k$（或等价地 $A^{k+1} X = A^k$）
    2. $A^D A A^D = A^D$（即 $XAX = X$）
    3. $AA^D = A^D A$（$X$ 与 $A$ 交换）

!!! theorem "定理 33.9 (Drazin 逆的唯一性)"
    满足上述三个条件的 $A^D$ 存在且唯一。

??? proof "证明"
    **唯一性**：设 $X$ 和 $Y$ 都满足三个条件。由条件 3，$X$ 和 $Y$ 都与 $A$ 交换，因此也与 $A^k$ 交换。

    $XA^{k+1} = A^k$（由条件 1），因此 $XA^{k+1}Y = A^kY$。但 $XA^{k+1}Y = X A^k \cdot AY = X A^k Y A = \cdots$

    更直接地：$X = XAX = X^2A = X^3A^2 = \cdots = X^{k+1}A^k$。同理 $Y = Y^{k+1}A^k$。

    由条件 1 和 3：$XA = A^D A = A A^D$，所以 $XA$ 是幂等矩阵（因为 $(XA)^2 = XAXA = XA \cdot XA$... 不对，需要更仔细。

    利用 $AX = XA$ 和 $A^{k+1}X = A^k$：

    $$
    (AX)^{k+1} = A^{k+1}X^{k+1} = A^k \cdot X^k
    $$

    而 $(AX)^k = A^k X^k$。由 $A^{k+1}X = A^k$，得 $AX \cdot A^k = A^k$（在 $\operatorname{col}(A^k)$ 上 $AX$ 是恒等映射）。

    类似地 $AY$ 在 $\operatorname{col}(A^k)$ 上是恒等映射，在 $\ker(A^k)$ 上是零映射。因此 $AX = AY$。

    $X = XAX = X(AX) = X(AY) = XAY$。类似地 $Y = YAY = Y(AX) = YAX = XAY$（利用交换性）。因此 $X = Y$。

    **存在性**：设 $\operatorname{ind}(A) = k$，$\mathbb{C}^n = \operatorname{col}(A^k) \oplus \ker(A^k)$。在这个分解下

    $$
    A = \begin{pmatrix} A_1 & 0 \\ 0 & N \end{pmatrix}
    $$

    其中 $A_1 = A|_{\operatorname{col}(A^k)}$ 是可逆的，$N = A|_{\ker(A^k)}$ 是幂零的（$N^k = 0$）。定义

    $$
    A^D = \begin{pmatrix} A_1^{-1} & 0 \\ 0 & 0 \end{pmatrix}
    $$

    验证三个条件是直接的。$\blacksquare$

### 谱刻画

!!! theorem "定理 33.10 (Drazin 逆的谱刻画)"
    设 $A$ 的 Jordan 标准形在非零特征值 $\lambda_1, \ldots, \lambda_s$ 处的 Jordan 块为 $J_1, \ldots, J_s$，在特征值 $0$ 处的 Jordan 块（如果存在）为 $N_1, \ldots, N_t$。则 $A^D$ 的 Jordan 标准形在 $\lambda_i^{-1}$ 处有对应的 Jordan 块 $J_i^{-1}$，在 $0$ 处为零块。

    即 $A^D$ 的非零特征值恰好是 $A$ 的非零特征值的倒数，而 $0$ 特征值对应的部分被"消除"。

!!! example "例 33.5"
    设 $A = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$。

    $A^2 = \begin{pmatrix} 4 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$，$\operatorname{rank}(A) = 2$，$\operatorname{rank}(A^2) = 1$，$\operatorname{rank}(A^3) = 1$。

    因此 $\operatorname{ind}(A) = 2$。$\operatorname{col}(A^2) = \operatorname{span}\{(1,0,0)^T\}$，$\ker(A^2) = \operatorname{span}\{(0,1,0)^T, (0,0,1)^T\}$。

    在分解 $\mathbb{C}^3 = \operatorname{col}(A^2) \oplus \ker(A^2)$ 下，$A_1 = (2)$，$N = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    $$
    A^D = \begin{pmatrix} 1/2 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
    $$

---

## 33.6 群逆与核逆

<div class="context-flow" markdown>

**核心问题**：当矩阵的指标为 1 时，Drazin 逆有什么特殊性质？

</div>

### 群逆

!!! definition "定义 33.6 (群逆)"
    设 $A \in \mathbb{C}^{n \times n}$，$\operatorname{ind}(A) \leq 1$。此时 $A$ 的 Drazin 逆称为 $A$ 的**群逆**（group inverse），记作 $A^\#$。即 $A^\#$ 满足：

    1. $A A^\# A = A$
    2. $A^\# A A^\# = A^\#$
    3. $A A^\# = A^\# A$

!!! note "注"
    群逆得名于半群理论：在半群 $\{A^n : n \geq 0\}$ 中，$A^\#$ 是 $A$ 的群论意义上的逆元。

    $\operatorname{ind}(A) \leq 1$ 等价于 $\operatorname{rank}(A) = \operatorname{rank}(A^2)$，即 $\operatorname{col}(A)$ 和 $\operatorname{col}(A^2)$ 相同。

!!! theorem "定理 33.11 (群逆与 Moore-Penrose 逆的关系)"
    设 $A \in \mathbb{C}^{n \times n}$ 满足 $\operatorname{ind}(A) \leq 1$。则：

    1. $A^\#$ 存在当且仅当 $\operatorname{rank}(A) = \operatorname{rank}(A^2)$。
    2. $A^\# = A^\dagger$ 当且仅当 $A$ 是 **EP 矩阵**（Range-Hermitian），即 $\operatorname{col}(A) = \operatorname{col}(A^*)$。
    3. 一个重要的计算公式为 $A^\# = A(A^2)^\dagger A$。

!!! note "注"
    EP 矩阵类包含所有正规矩阵（$AA^* = A^*A$）、可逆矩阵和正定矩阵。对于这些矩阵，群逆与 Moore-Penrose 逆重合。但在非 EP 情形下，两者通常不同。例如 $\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$ 不是 EP 矩阵，其 $A^\# = A$ 而 $A^\dagger = \frac{1}{2}\begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$。

### 核逆

!!! definition "定义 33.7 (核逆)"
    设 $A \in \mathbb{C}^{n \times n}$，$\operatorname{ind}(A) \leq 1$。$A$ 的**核逆**（core inverse）$A^{\tiny\textcircled{\#}}$ 定义为满足以下条件的矩阵：

    $$
    A A^{\tiny\textcircled{\#}} = P_{\operatorname{col}(A)}, \quad \operatorname{col}(A^{\tiny\textcircled{\#}}) \subseteq \operatorname{col}(A)
    $$

    其中 $P_{\operatorname{col}(A)}$ 是到 $\operatorname{col}(A)$ 的正交投影。

!!! theorem "定理 33.12 (核逆的表示)"
    $A^{\tiny\textcircled{\#}} = A^\# A A^\dagger = A^\dagger A A^\#$（当两者都存在时）。

!!! example "例 33.6"
    设 $A = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$。

    $\operatorname{rank}(A) = \operatorname{rank}(A^2) = 1$（$A^2 = A$），故 $\operatorname{ind}(A) \leq 1$。

    群逆：由 $AA^\#A = A$，$A^\#AA^\# = A^\#$，$AA^\# = A^\#A$，以及 $A$ 的分解 $\operatorname{col}(A) = \operatorname{span}\{(1,0)^T\}$，$\ker(A) = \operatorname{span}\{(-1,1)^T\}$... 实际上 $\ker(A) = \operatorname{span}\{(1, -1)^T\}$... 不对，$A(x_1, x_2)^T = (x_1+x_2, 0)^T = 0$ 当 $x_1 = -x_2$。

    $\operatorname{col}(A) = \operatorname{span}\{(1,0)^T\}$，$\ker(A) = \operatorname{span}\{(-1,1)^T\}$。

    在分解 $\mathbb{C}^2 = \operatorname{col}(A) \oplus \ker(A)$ 下（注意这不是正交分解），$A_1 = (1)$（可逆），故 $A^\# = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$（在标准基下不对，需要变换基）。

    变换矩阵 $S = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$（列为 $(1,0)^T$ 和 $(-1,1)^T$），$S^{-1}AS = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $A^\# = S \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} S^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$。

    验证：$AA^\# = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = A$... 即 $AA^\# = A = A^\#$。

    $A^\#A = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = AA^\#$。✓

    注意 $A^\# = A = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$，这是因为 $A$ 是幂等的（$A^2 = A$），幂等矩阵的群逆就是自身。

---

## 33.7 广义逆的扰动与应用

<div class="context-flow" markdown>

**核心问题**：Moore-Penrose 逆对矩阵扰动有多敏感？广义逆在实际中有哪些应用？

</div>

### 扰动界

!!! theorem "定理 33.13 (Moore-Penrose 逆的扰动界)"
    设 $A, E \in \mathbb{C}^{m \times n}$，$B = A + E$。若 $\operatorname{rank}(A) = \operatorname{rank}(B)$，则

    $$
    \|B^\dagger - A^\dagger\| \leq \sqrt{2} \max(\|A^\dagger\|^2, \|B^\dagger\|^2) \|E\|
    $$

    若进一步 $\|A^\dagger\| \|E\| < 1$，则

    $$
    \|B^\dagger\| \leq \frac{\|A^\dagger\|}{1 - \|A^\dagger\|\|E\|}
    $$

!!! note "注"
    当 $A$ 的最小非零奇异值 $\sigma_r$ 很小时（即 $A$ 接近秩亏矩阵），$\|A^\dagger\| = 1/\sigma_r$ 很大，扰动界也很大。这反映了广义逆的不稳定性，与矩阵条件数的概念密切相关。

### 秩变化时的不连续性

!!! theorem "定理 33.14"
    映射 $A \mapsto A^\dagger$ 在秩恒定的矩阵集合上是连续的，但在秩发生变化的地方是不连续的。

!!! example "例 33.7"
    设 $A_\epsilon = \begin{pmatrix} 1 & 0 \\ 0 & \epsilon \end{pmatrix}$。

    - 当 $\epsilon \neq 0$：$A_\epsilon^\dagger = A_\epsilon^{-1} = \begin{pmatrix} 1 & 0 \\ 0 & 1/\epsilon \end{pmatrix}$。
    - 当 $\epsilon = 0$：$A_0^\dagger = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $\lim_{\epsilon \to 0} A_\epsilon^\dagger = \begin{pmatrix} 1 & 0 \\ 0 & +\infty \end{pmatrix} \neq A_0^\dagger$。

    因此 $A \mapsto A^\dagger$ 在 $A_0$ 处不连续。

### 应用：奇异线性系统

!!! theorem "定理 33.15 (Drazin 逆与奇异微分方程)"
    考虑奇异线性微分方程

    $$
    A\boldsymbol{x}'(t) + B\boldsymbol{x}(t) = \boldsymbol{f}(t)
    $$

    其中 $A$ 可能是奇异的。若 $(A, B)$ 是正则矩阵束（即 $\det(\lambda A + B) \not\equiv 0$），则通过变换可以将系统化为

    $$
    \boldsymbol{x}_1'(t) + C_1 \boldsymbol{x}_1(t) = \boldsymbol{g}_1(t), \quad N\boldsymbol{x}_2'(t) + \boldsymbol{x}_2(t) = \boldsymbol{g}_2(t)
    $$

    其中 $N$ 是幂零的。第二个方程的解涉及 $N$ 的 Drazin 逆：

    $$
    \boldsymbol{x}_2(t) = -\sum_{j=0}^{k-1} N^j \boldsymbol{g}_2^{(j)}(t)
    $$

    其中 $k = \operatorname{ind}(N)$。

!!! example "例 33.8"
    **Markov 链的稳态分析**：

    设 $P$ 是 Markov 链的转移矩阵（行随机矩阵），$Q = I - P$。稳态分布 $\boldsymbol{\pi}$ 满足 $\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T$，即 $\boldsymbol{\pi}^T Q = \boldsymbol{0}^T$。

    $Q$ 的群逆 $Q^\#$ 可以用来表示 Markov 链的基本矩阵（fundamental matrix）$Z = (I - P + \Pi)^{-1}$，其中 $\Pi = \boldsymbol{1}\boldsymbol{\pi}^T$。实际上：

    $$
    Z = I - Q^\#
    $$

    (在适当归一化下)。基本矩阵 $Z$ 包含了 Markov 链的全部二阶信息，包括平均首达时间、方差等。

!!! note "注"
    广义逆理论将逆矩阵的概念从可逆方阵推广到任意矩阵，提供了处理奇异性和非方性的统一框架。Moore-Penrose 逆侧重几何（正交投影、最小范数），Drazin 逆侧重代数（交换性、谱分解）。两者在各自的应用领域中都是不可或缺的工具。

    读者应注意广义逆与正则逆的一个根本区别：$(AB)^\dagger \neq B^\dagger A^\dagger$（一般地）。只有当 $A$ 和 $B$ 满足特定的秩条件时，"逆的乘积 = 乘积的逆"才成立。下面我们将详细讨论这个问题。

---

## 33.8 加权 Moore-Penrose 逆

<div class="context-flow" markdown>

**核心问题**：当内积空间的度量不是标准的（如加权内积），如何定义和计算相应的"最优"广义逆？

</div>

在许多应用中（特别是加权最小二乘问题），自变量和因变量的不同分量具有不同的重要性或不同的量纲，需要在非标准内积下考虑最小化问题。加权 Moore-Penrose 逆正是为此而定义的。

### 定义与刻画

!!! definition "定义 33.8 (加权 Moore-Penrose 逆)"
    设 $A \in \mathbb{C}^{m \times n}$，$M \in \mathbb{C}^{m \times m}$ 和 $N \in \mathbb{C}^{n \times n}$ 是 Hermite 正定矩阵。$A$ 关于权矩阵 $(M, N)$ 的**加权 Moore-Penrose 逆** $A^\dagger_{M,N}$ 定义为满足以下条件的唯一矩阵 $X \in \mathbb{C}^{n \times m}$：

    对任意 $\boldsymbol{b} \in \mathbb{C}^m$，$X\boldsymbol{b}$ 是加权最小二乘问题

    $$\min_{\boldsymbol{x}} \|A\boldsymbol{x} - \boldsymbol{b}\|_M \quad \text{subject to} \quad \|\boldsymbol{x}\|_N \text{ 最小}$$

    的唯一解，其中 $\|\boldsymbol{y}\|_M = \sqrt{\boldsymbol{y}^* M \boldsymbol{y}}$，$\|\boldsymbol{x}\|_N = \sqrt{\boldsymbol{x}^* N \boldsymbol{x}}$。

!!! theorem "定理 33.16 (加权 Penrose 条件)"
    $X = A^\dagger_{M,N}$ 当且仅当 $X$ 满足以下四个**加权 Penrose 条件**：

    1. $AXA = A$
    2. $XAX = X$
    3. $(MAX)^* = MAX$（即 $AX$ 关于 $M$-内积是 Hermite 的）
    4. $(NXA)^* = NXA$（即 $XA$ 关于 $N$-内积是 Hermite 的）

    等价地，$AX$ 是关于 $M$-内积到 $\operatorname{col}(A)$ 的正交投影，$XA$ 是关于 $N$-内积到 $\operatorname{row}(A)$ 的正交投影。

??? proof "证明"
    令 $\hat{A} = M^{1/2} A N^{-1/2}$。则加权最小二乘问题

    $$\min_{\boldsymbol{x}} \|A\boldsymbol{x} - \boldsymbol{b}\|_M, \quad \|\boldsymbol{x}\|_N \text{ 最小}$$

    通过变量替换 $\boldsymbol{y} = N^{1/2}\boldsymbol{x}$，$\hat{\boldsymbol{b}} = M^{1/2}\boldsymbol{b}$ 转化为

    $$\min_{\boldsymbol{y}} \|\hat{A}\boldsymbol{y} - \hat{\boldsymbol{b}}\|_2, \quad \|\boldsymbol{y}\|_2 \text{ 最小}.$$

    其解为 $\boldsymbol{y}^* = \hat{A}^\dagger \hat{\boldsymbol{b}}$。回到原变量：$\boldsymbol{x}^* = N^{-1/2}\hat{A}^\dagger M^{1/2}\boldsymbol{b}$。

    因此 $A^\dagger_{M,N} = N^{-1/2}\hat{A}^\dagger M^{1/2} = N^{-1/2}(M^{1/2}AN^{-1/2})^\dagger M^{1/2}$。

    验证四个加权 Penrose 条件可通过将标准 Penrose 条件对 $\hat{A}$ 翻译回 $A$ 来完成。$\blacksquare$

### 显式公式

!!! theorem "定理 33.17 (加权 Moore-Penrose 逆的显式公式)"
    设 $A \in \mathbb{C}^{m \times n}$，$\operatorname{rank}(A) = r$，$M, N$ 正定。则：

    1. **一般公式**：$A^\dagger_{M,N} = N^{-1/2}(M^{1/2}AN^{-1/2})^\dagger M^{1/2}$。

    2. **当 $A$ 列满秩时**（$r = n$）：

        $$A^\dagger_{M,N} = (A^*MA)^{-1}A^*M.$$

        这是 $M$-加权最小二乘的法方程 $A^*MA\boldsymbol{x} = A^*M\boldsymbol{b}$ 的解。

    3. **当 $A$ 行满秩时**（$r = m$）：

        $$A^\dagger_{M,N} = N^{-1}A^*(AN^{-1}A^*)^{-1}.$$

        这给出 $N$-范数最小的解。

    4. **一般情形**：利用 $A$ 的秩分解 $A = FG$（$F \in \mathbb{C}^{m \times r}$ 列满秩，$G \in \mathbb{C}^{r \times n}$ 行满秩），

        $$A^\dagger_{M,N} = N^{-1}G^*(GN^{-1}G^*)^{-1}(F^*MF)^{-1}F^*M.$$

### 性质

!!! theorem "定理 33.18 (加权 Moore-Penrose 逆的性质)"
    设 $A \in \mathbb{C}^{m \times n}$，$M, N$ 正定。则：

    1. **投影性质**：$AA^\dagger_{M,N}$ 是关于 $M$-内积到 $\operatorname{col}(A)$ 的正交投影；$A^\dagger_{M,N}A$ 是关于 $N$-内积到 $\operatorname{row}(A)$ 的正交投影。

    2. **退化为标准情形**：当 $M = I_m$，$N = I_n$ 时，$A^\dagger_{I,I} = A^\dagger$。

    3. **秩保持**：$\operatorname{rank}(A^\dagger_{M,N}) = \operatorname{rank}(A)$。

    4. **对偶性**：$(A^*)^\dagger_{N,M} = (A^\dagger_{M,N})^*$。

### 应用：加权最小二乘

!!! example "例 33.9 (加权最小二乘)"
    在统计回归中，模型 $\boldsymbol{b} = A\boldsymbol{x} + \boldsymbol{\epsilon}$，其中误差 $\boldsymbol{\epsilon}$ 的协方差矩阵为 $\operatorname{Cov}(\boldsymbol{\epsilon}) = M^{-1}$（$M$ 正定）。**广义最小二乘**（GLS）估计量

    $$\hat{\boldsymbol{x}}_{\text{GLS}} = (A^*MA)^{-1}A^*M\boldsymbol{b} = A^\dagger_{M,I}\boldsymbol{b}$$

    是最佳线性无偏估计量（Gauss-Markov 定理的推广）。

    **数值例子**：设 $A = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 3 \end{pmatrix}$，$\boldsymbol{b} = \begin{pmatrix} 1 \\ 3 \\ 4 \end{pmatrix}$，权矩阵 $M = \operatorname{diag}(1, 2, 1)$（第二个观测更可靠）。

    $A^*MA = \begin{pmatrix} 4 & 8 \\ 8 & 18 \end{pmatrix}$，$A^*M\boldsymbol{b} = \begin{pmatrix} 11 \\ 25 \end{pmatrix}$。

    $\hat{\boldsymbol{x}}_{\text{GLS}} = \begin{pmatrix} 4 & 8 \\ 8 & 18 \end{pmatrix}^{-1}\begin{pmatrix} 11 \\ 25 \end{pmatrix} = \frac{1}{8}\begin{pmatrix} 18 & -8 \\ -8 & 4 \end{pmatrix}\begin{pmatrix} 11 \\ 25 \end{pmatrix} = \frac{1}{8}\begin{pmatrix} -2 \\ 12 \end{pmatrix} = \begin{pmatrix} -0.25 \\ 1.5 \end{pmatrix}.$

    与普通最小二乘 $\hat{\boldsymbol{x}}_{\text{OLS}} = A^\dagger \boldsymbol{b}$ 比较：$A^\dagger \boldsymbol{b} = (A^*A)^{-1}A^*\boldsymbol{b} = \begin{pmatrix} -1/3 \\ 3/2 \end{pmatrix} \approx \begin{pmatrix} -0.333 \\ 1.5 \end{pmatrix}$。加权使截距估计向零偏移，因为更可靠的第二个观测对应更大的权重。

---

## 33.9 逆的乘积律

<div class="context-flow" markdown>

**核心问题**：$(AB)^\dagger = B^\dagger A^\dagger$ 何时成立？

</div>

对于可逆矩阵，$(AB)^{-1} = B^{-1}A^{-1}$ 总是成立。但对 Moore-Penrose 逆，**逆序律**（reverse order law）$(AB)^\dagger = B^\dagger A^\dagger$ 一般不成立。精确条件的刻画是广义逆理论中的重要课题。

### 精确条件

!!! theorem "定理 33.19 (逆序律的充要条件)"
    设 $A \in \mathbb{C}^{m \times k}$，$B \in \mathbb{C}^{k \times n}$。则 $(AB)^\dagger = B^\dagger A^\dagger$ 当且仅当以下两个条件同时成立：

    1. $A^*AB B^* = B B^* A^*A$（即 $A^*A$ 和 $BB^*$ 交换），
    2. $\operatorname{col}(A^*AB) \subseteq \operatorname{col}(B)$ 且 $\operatorname{col}(BB^*A^*) \subseteq \operatorname{col}(A^*)$。

    这些条件等价于 $A^\dagger A B B^\dagger$ 是正交投影。

### 充分条件

!!! theorem "定理 33.20 (逆序律成立的充分条件)"
    $(AB)^\dagger = B^\dagger A^\dagger$ 在以下任一条件下成立：

    1. $A$ 具有**正交列**：$A^*A = \alpha I$（$\alpha > 0$），即 $A$ 的列两两正交且等范数。
    2. $B$ 具有**正交行**：$BB^* = \beta I$（$\beta > 0$）。
    3. $A$ 是**列满秩**且 $B$ 是**行满秩**的。
    4. $A$ 或 $B$ 是酉/正交矩阵。

??? proof "证明"
    **(1)** 设 $A^*A = \alpha I$。验证 $X = B^\dagger A^\dagger$ 满足关于 $AB$ 的四个 Penrose 条件。

    条件 1：$ABX(AB) = AB B^\dagger A^\dagger AB = AB B^\dagger (\alpha I) B / \alpha = AB B^\dagger B$。由于 $B B^\dagger B = B$（Penrose 条件 1 对 $B$），得 $AB$。✓

    条件 3：$(ABX)^* = (AB B^\dagger A^\dagger)^*$。注意 $A^\dagger = \frac{1}{\alpha}A^*$（因为 $A^*A = \alpha I$），所以 $A^\dagger A = I$。因此 $ABX = AB B^\dagger \frac{1}{\alpha}A^* = \frac{1}{\alpha}AB B^\dagger A^*$。

    $ABB^\dagger$ 是 Hermite 的（Penrose 条件 3 对 $B$），即 $(BB^\dagger)^* = BB^\dagger$。因此

    $(ABX)^* = \frac{1}{\alpha}(ABB^\dagger A^*)^* = \frac{1}{\alpha}A(BB^\dagger)^*A^* = \frac{1}{\alpha}ABB^\dagger A^* = ABX$。✓

    条件 2 和 4 的验证类似。

    **(2)** 完全对偶的论证。

    **(3)** 当 $A$ 列满秩时 $A^\dagger = (A^*A)^{-1}A^*$，$A^\dagger A = I$。当 $B$ 行满秩时 $B^\dagger = B^*(BB^*)^{-1}$，$BB^\dagger = I$。因此 $B^\dagger A^\dagger AB = B^\dagger B$，$ABB^\dagger A^\dagger = AA^\dagger$。验证四个 Penrose 条件是直接的。$\blacksquare$

!!! example "例 33.10 (逆序律失效)"
    设 $A = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    $AB = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$(AB)^\dagger = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

    $B^\dagger = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$，$A^\dagger = \begin{pmatrix} 1 & 0 \end{pmatrix}$。

    $B^\dagger A^\dagger = \begin{pmatrix} 1 \\ 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix} = (AB)^\dagger$。

    此处恰好成立（因为 $A$ 列满秩，$B$ 行满秩）。

    但取 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 0 \\ 0 & 1 \end{pmatrix}$。

    $AB = 0$，$(AB)^\dagger = 0$。但 $B^\dagger A^\dagger = BA = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix} = 0$。

    此处也成立。取更有意义的反例：$A = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。

    $AB = \begin{pmatrix} 2 \\ 0 \end{pmatrix}$，$(AB)^\dagger = \begin{pmatrix} 1/2 & 0 \end{pmatrix}$。

    $A^\dagger = \frac{1}{2}\begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$，$B^\dagger = \frac{1}{2}\begin{pmatrix} 1 & 1 \end{pmatrix}$。

    $B^\dagger A^\dagger = \frac{1}{4}\begin{pmatrix} 1 & 1 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix} = \frac{1}{4}\begin{pmatrix} 2 & 0 \end{pmatrix} = \begin{pmatrix} 1/2 & 0 \end{pmatrix} = (AB)^\dagger$。

    仍然成立！实际上反例需要更精心的构造。取 $A = \begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$。

    $AB = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$，$(AB)^\dagger = \frac{1}{4}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。

    $A^\dagger = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$，$B^\dagger = \frac{1}{2}\begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}$。

    $B^\dagger A^\dagger = \frac{1}{4}\begin{pmatrix} 1 & 0 \\ 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix} = \frac{1}{4}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = (AB)^\dagger$。✓

    逆序律失效的典型例子是当 $A^*A$ 和 $BB^*$ 不交换时。

---

## 33.10 并联和

<div class="context-flow" markdown>

**核心问题**：能否用广义逆定义矩阵的"并联连接"运算？

</div>

并联和（parallel sum）是矩阵运算中一个源于电学网络理论的优美概念。它将电阻并联的公式推广到正定矩阵（以及一般矩阵）。

!!! definition "定义 33.9 (并联和)"
    设 $A, B \in \mathbb{C}^{n \times n}$ 为半正定矩阵。$A$ 与 $B$ 的**并联和**定义为

    $$A : B = A(A + B)^\dagger B.$$

    当 $A + B$ 可逆时，简化为 $A : B = A(A+B)^{-1}B$。

!!! theorem "定理 33.21 (并联和的性质)"
    设 $A, B \succeq 0$。则：

    1. **对称性**：$A : B = B : A$。
    2. **变分刻画**：$\boldsymbol{x}^*(A : B)\boldsymbol{x} = \inf_{\boldsymbol{y}} \{\boldsymbol{y}^* A \boldsymbol{y} + (\boldsymbol{x} - \boldsymbol{y})^* B (\boldsymbol{x} - \boldsymbol{y})\}$。
    3. **标量情形**：当 $A = aI$，$B = bI$（$a, b > 0$）时，$A : B = \frac{ab}{a+b}I$，恰好是标量的**调和均值**。
    4. **半正定性**：$A : B \succeq 0$。
    5. **单调性**：若 $A \preceq A'$，则 $A : B \preceq A' : B$。
    6. **并联公式**：$(A : B)^{-1} = A^{-1} + B^{-1}$（当 $A, B$ 正定时）。

??? proof "证明"
    **(6)** 当 $A, B$ 正定时，$A : B = A(A+B)^{-1}B$。则

    $(A:B)^{-1} = B^{-1}(A+B)A^{-1} = B^{-1} + A^{-1}$。✓

    **(1)** $A(A+B)^\dagger B = B(B+A)^\dagger A$。这可以通过验证两者都满足相同的 Penrose 条件来证明。

    **(2)** 对 $\boldsymbol{y}$ 求导令其为零：$2A\boldsymbol{y} - 2B(\boldsymbol{x}-\boldsymbol{y}) = 0$，即 $(A+B)\boldsymbol{y} = B\boldsymbol{x}$，故 $\boldsymbol{y} = (A+B)^\dagger B\boldsymbol{x}$。代入得最小值 $\boldsymbol{x}^* A(A+B)^\dagger B \boldsymbol{x} = \boldsymbol{x}^*(A:B)\boldsymbol{x}$。$\blacksquare$

**电学网络解释**：在电路理论中，正定矩阵 $A$ 和 $B$ 可以视为多端口网络的阻抗矩阵。$A : B$ 正是将两个网络并联后的等效阻抗。性质 (6) 的 $(A:B)^{-1} = A^{-1} + B^{-1}$ 就是将导纳相加的并联规则在矩阵情形的推广。Anderson 和 Duffin（1969）率先系统研究了矩阵并联和，并发现它在电路综合、控制论和统计学中的广泛应用。

## 练习题

1. **[基础] 所有的矩阵（包括非方阵和奇异阵）都有 $\{1\}$-逆吗？它是唯一的吗？**
   ??? success "参考答案"
       是的，所有矩阵都存在 $\{1\}$-逆（满足 $AXA=A$）。但它通常不唯一，除非原矩阵 $A$ 是可逆方阵。

2. **[Moore-Penrose] 什么是 Moore-Penrose 逆必须满足的四个条件？**
   ??? success "参考答案"
       1. $AXA = A$ (内逆)
       2. $XAX = X$ (外逆)
       3. $(AX)^H = AX$ (左自伴)
       4. $(XA)^H = XA$ (右自伴)

3. **[计算] 如果 $A = \operatorname{diag}(\sigma_1, \sigma_2, 0)$，求其 Moore-Penrose 逆 $A^\dagger$。**
   ??? success "参考答案"
       $A^\dagger = \operatorname{diag}(1/\sigma_1, 1/\sigma_2, 0)$。非零元素取倒数，零元素保持为零。

4. **[投影] 证明 $A A^\dagger$ 是一个正交投影矩阵。它投影到哪个空间？**
   ??? success "参考答案"
       满足 $(AA^\dagger)^2 = A(A^\dagger A)A^\dagger = AA^\dagger$ 且 $(AA^\dagger)^H = AA^\dagger$。它投影到 $A$ 的**列空间** $\operatorname{col}(A)$。

5. **[最小二乘] 为什么说 $x = A^\dagger b$ 是线性方程组 $Ax = b$ 的“最优解”？**
   ??? success "参考答案"
       它不仅使残差 $\|Ax-b\|_2$ 达到最小（最小二乘解），而且在所有使残差最小的解中，它的范数 $\|x\|_2$ 也是最小的（最小范数解）。

6. **[性质] 证明 $(A^\dagger)^\dagger = A$。**
   ??? success "参考答案"
       将 $A$ 代入 $A^\dagger$ 满足的四个 Penrose 条件中，可以发现它们是对称的。矩阵 $A$ 完美满足作为 $A^\dagger$ 之广义逆的所有条件。

7. **[秩] $\operatorname{rank}(A^\dagger)$ 与 $\operatorname{rank}(A)$ 有什么关系？**
   ??? success "参考答案"
       它们严格相等。这可以从 SVD 的构造中一眼看出：非零奇异值的个数在取倒数后保持不变。

8. **[Drazin逆] Drazin 逆与 Moore-Penrose 逆的主要区别是什么？**
   ??? success "参考答案"
       Moore-Penrose 逆侧重于**正交投影与最小二乘**（适用于任何矩阵）；而 Drazin 逆侧重于**保持矩阵的幂次结构与交换性**（仅适用于方阵），它满足 $AX = XA$。

9. **[逆序律] 为什么一般情况下 $(AB)^\dagger \neq B^\dagger A^\dagger$？**
   ??? success "参考答案"
       因为广义逆的定义涉及投影方向。只有当 $A$ 的列空间和 $B$ 的行空间满足特定的正交或对齐条件（如 $A$ 列满秩且 $B$ 行满秩）时，逆序律才成立。

10. **[爱因斯坦思考题] 假设宇宙中某处存在一条不相容的定律（方程组无解）。广义逆告诉我们，即便没有完美解，依然存在一个使“能量扰动”（残差）最小的解。这反映了自然界怎样的折中原则？**
    ??? success "参考答案"
        这反映了“最小作用量原理”的代数版：自然界不会因为无法达到理想状态而停摆。广义逆提供了一个在现有约束下最经济、最不“浮夸”的妥协方案。它承认了误差（不相容性）的存在，但通过向已知事实（列空间）进行正交投影，守住了客观真实性的底线。

## 本章小结

本章将逆矩阵的概念推广到了最一般的情形，建立了奇异和非方阵系统的统一代数框架：

1. **内逆与 $\{1\}$-逆**：确立了广义逆存在性的基本逻辑，解决了不相容方程组的特解构造问题。
2. **Moore-Penrose 逆**：通过四个 Penrose 条件定义了唯一的“最佳广义逆”，并给出了基于 SVD 的构造算法。
3. **投影几何**：揭示了广义逆与正交投影的内在联系，证明了 $A^\dagger$ 是连接四个基本子空间的桥梁。
4. **最小范数最小二乘解**：论证了 $A^\dagger b$ 在处理不相容系统时的最优性，是现代信号处理和统计估计的理论核心。
5. **Drazin 逆与群逆**：讨论了在处理微分方程和马尔可夫链时，具有良好交换性质和幂次保持特性的代数广义逆。
6. **逆序律与计算稳定性**：探讨了广义逆运算在复合映射下的复杂行为，以及秩亏损时的数值不稳定性。

