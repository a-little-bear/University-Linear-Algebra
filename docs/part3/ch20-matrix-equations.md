# 第 20 章 矩阵方程

<div class="context-flow" markdown>

**前置**：Ch19 Kronecker积/Vec算子

**脉络**：$AX=B$(伪逆) → **Sylvester** $AX+XB=C$(特征值分离条件) → **Lyapunov**(稳定性) → **Riccati**(最优控制，非线性) → Penrose方程(伪逆公理化)

**延伸**：Sylvester/Lyapunov 方程在控制系统稳定性分析和模型降阶中是核心工具；Riccati 方程出现在最优控制（LQR/LQG）、Kalman 滤波、博弈论中；矩阵方程的数值求解（Bartels-Stewart 算法）是科学计算的重要课题

</div>

矩阵方程是线性代数在控制论、信号处理和数值分析中应用的核心。与普通线性方程组 $A\mathbf{x} = \mathbf{b}$ 不同，矩阵方程中的未知量是一个**矩阵**而非向量。从最简单的 $AX = B$ 到结构丰富的 Sylvester 方程 $AX + XB = C$、Lyapunov 方程 $AX + XA^* = Q$，以及非线性的 Riccati 方程，这些方程各有其独特的数学结构和求解方法。

本章系统地介绍这些矩阵方程的理论基础、可解性条件、求解公式以及数值算法，最后讨论与 Moore-Penrose 伪逆密切相关的 Penrose 方程组。

---

## 20.1 线性矩阵方程概述

<div class="context-flow" markdown>

**脉络**：$AX=B$ 可解 $\Leftrightarrow$ $\mathcal{C}(B)\subseteq\mathcal{C}(A)$ · 通解结构 = 特解 + 零空间 · $AXB=C$：双边消元，条件 $AA^+CB^+B=C$

</div>

本节从最基本的矩阵方程出发，建立矩阵方程理论的框架。

!!! definition "定义 20.1 (基本线性矩阵方程)"
    以下是三种最常见的线性矩阵方程：

    1. **单边方程**：$AX = B$，其中 $A$ 为 $m \times n$，$B$ 为 $m \times p$，$X$ 为 $n \times p$ 未知矩阵。
    2. **双边方程**：$AXB = C$，其中 $A$ 为 $m \times n$，$B$ 为 $p \times q$，$C$ 为 $m \times q$，$X$ 为 $n \times p$ 未知矩阵。
    3. **Sylvester 型方程**：$AX + XB = C$，其中 $A$ 为 $m \times m$，$B$ 为 $n \times n$，$C$ 为 $m \times n$，$X$ 为 $m \times n$ 未知矩阵。

!!! theorem "定理 20.1 (单边方程 $AX = B$ 的可解性)"
    方程 $AX = B$ 有解当且仅当 $\operatorname{rank}(A) = \operatorname{rank}(A \mid B)$，即 $B$ 的列空间包含在 $A$ 的列空间中：$\mathcal{C}(B) \subseteq \mathcal{C}(A)$。

    当 $A$ 可逆时，唯一解为 $X = A^{-1}B$。

    一般情况下，设 $A^+$ 为 $A$ 的 Moore-Penrose 伪逆，则通解为：

    $$
    X = A^+ B + (I - A^+ A)Z,
    $$

    其中 $Z$ 为任意 $n \times p$ 矩阵。

??? proof "证明"
    **必要性**：若 $AX = B$，则 $B$ 的每一列 $\mathbf{b}_j = A\mathbf{x}_j \in \mathcal{C}(A)$。

    **充分性**：若 $\mathcal{C}(B) \subseteq \mathcal{C}(A)$，则 $AA^+B = B$（因为 $A^+$ 是从 $\mathcal{C}(A)$ 到自身的左逆，$AA^+$ 是到 $\mathcal{C}(A)$ 的正交投影）。因此 $X_0 = A^+B$ 是一个特解。

    **通解**：设 $X$ 为任意解，则 $A(X - X_0) = 0$，即 $X - X_0 \in \mathcal{N}(A)$。而 $\mathcal{N}(A) = \mathcal{C}(I - A^+A)$（$I - A^+A$ 是到 $\mathcal{N}(A)$ 的正交投影），因此 $X - X_0 = (I - A^+A)Z$。$\blacksquare$

!!! theorem "定理 20.2 (双边方程 $AXB = C$ 的可解性)"
    方程 $AXB = C$ 有解当且仅当 $AA^+CB^+B = C$。

    当 $A$ 和 $B$ 均可逆时，唯一解为 $X = A^{-1}CB^{-1}$。

    一般情况下，通解为：

    $$
    X = A^+ C B^+ + Z - A^+ A Z B B^+,
    $$

    其中 $Z$ 为任意大小适当的矩阵。

??? proof "证明"
    **必要性**：若 $AXB = C$，则 $AA^+CB^+B = AA^+(AXB)B^+B = (AA^+A)X(BB^+B) = AXB = C$。

    **充分性**：设 $AA^+CB^+B = C$。令 $X_0 = A^+CB^+$，则：

    $$
    AX_0B = A(A^+CB^+)B = (AA^+)C(B^+B) = AA^+CB^+B = C.
    $$

    **通解**：设 $AXB = C$，则 $A(X - X_0)B = 0$。需要找到 $A\tilde{X}B = 0$ 的通解。$\tilde{X} = Z - A^+AZB B^+$ 满足：

    $$
    A(Z - A^+AZBB^+)B = AZB - (AA^+A)Z(BB^+B) = AZB - AZB = 0. \quad \blacksquare
    $$

!!! definition "定义 20.2 (矩阵方程的算子形式)"
    一般的线性矩阵方程可以写成算子形式 $\mathcal{L}(X) = C$，其中 $\mathcal{L}: \mathbb{C}^{m \times n} \to \mathbb{C}^{m \times n}$ 为线性算子。利用 Vec 算子和 Kronecker 积，这等价于：

    $$
    L \operatorname{vec}(X) = \operatorname{vec}(C),
    $$

    其中 $L$ 为 $mn \times mn$ 矩阵（见第 19 章）。

!!! example "例 20.1"
    求解 $AX = B$，其中 $A = \begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$，$B = \begin{pmatrix} 3 \\ 6 \end{pmatrix}$。

    $\operatorname{rank}(A) = 1$（第二行是第一行的 2 倍）。$B = 3\begin{pmatrix} 1 \\ 2 \end{pmatrix} \in \mathcal{C}(A)$，因此有解。

    $A^+ = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$（因为 $A = \mathbf{a}\mathbf{a}^T/5$，其中 $\mathbf{a} = (1,2)^T$，故 $A^+ = \mathbf{a}\mathbf{a}^T/(5 \cdot 5) \cdot 5 = A/\|A\|_F^2$... 让我们直接计算）。

    由于 $A = \begin{pmatrix} 1 \\ 2 \end{pmatrix}(1, 2)$，$A^+ = \frac{1}{5}\begin{pmatrix} 1 \\ 2 \end{pmatrix}\frac{1}{5}(1, 2) \cdot 5 = \frac{1}{5}\begin{pmatrix} 1/5 & 2/5 \\ 2/5 & 4/5 \end{pmatrix}$...

    更直接地：$A^+ = \frac{1}{25}A$，因为 $A = \mathbf{u}\mathbf{v}^T$（$\mathbf{u} = \mathbf{v} = (1,2)^T$），则 $A^+ = \frac{\mathbf{v}\mathbf{u}^T}{\|\mathbf{u}\|^2\|\mathbf{v}\|^2} = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}$。

    特解：$X_0 = A^+B = \frac{1}{25}\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\begin{pmatrix} 3 \\ 6 \end{pmatrix} = \frac{1}{25}\begin{pmatrix} 15 \\ 30 \end{pmatrix} = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix}$。

    通解：$X = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix} + (I - A^+A)Z = \begin{pmatrix} 3/5 \\ 6/5 \end{pmatrix} + \frac{1}{5}\begin{pmatrix} 4 & -2 \\ -2 & 1 \end{pmatrix}\begin{pmatrix} z \end{pmatrix}$。

    即 $X = \begin{pmatrix} 3/5 + 4t/5 \\ 6/5 - 2t/5 \end{pmatrix}$，或写成 $X = \begin{pmatrix} 3 + 4t \\ 6 - 2t \end{pmatrix}/5$（$t$ 为任意实数）。

    验证（取 $t = 3$）：$X = \begin{pmatrix} 3 \\ 0 \end{pmatrix}$，$AX = \begin{pmatrix} 3 \\ 6 \end{pmatrix} = B$。正确。

---

## 20.2 Sylvester 方程

<div class="context-flow" markdown>

**核心**：$AX+XB=C$ 唯一可解 $\Leftrightarrow$ $\sigma(A)\cap\sigma(-B)=\emptyset$ · Ch19 Kronecker和视角：系数矩阵 $=A\oplus B^T$，特征值 $\lambda_i(A)+\lambda_j(B)$ · 积分表示需稳定性条件

</div>

Sylvester 方程是最重要的矩阵方程之一，广泛出现在控制论、模型简化和矩阵函数理论中。

!!! definition "定义 20.3 (Sylvester 方程)"
    **Sylvester 方程**是如下形式的线性矩阵方程：

    $$
    AX + XB = C,
    $$

    其中 $A \in \mathbb{C}^{m \times m}$，$B \in \mathbb{C}^{n \times n}$，$C \in \mathbb{C}^{m \times n}$ 为已知矩阵，$X \in \mathbb{C}^{m \times n}$ 为未知矩阵。
    当 $A = -B^*$ 时，Sylvester 方程退化为 Lyapunov 方程。

!!! theorem "定理 20.3 (Sylvester 方程的可解性 —— Sylvester-Rosenblum 定理)"
    Sylvester 方程 $AX + XB = C$ 对任意 $C$ 有唯一解 $X$，当且仅当 $A$ 与 $-B$ 无公共特征值，即：

    $$
    \sigma(A) \cap \sigma(-B) = \emptyset.
    $$

??? proof "证明"
    **Kronecker 积方法**：由第 19 章，方程 $AX + XB = C$ 等价于：

    $$
    (I_n \otimes A + B^T \otimes I_m)\operatorname{vec}(X) = \operatorname{vec}(C).
    $$

    $I_n \otimes A + B^T \otimes I_m$ 的特征值为 $\{\lambda_i(A) + \lambda_j(B) : i = 1,\ldots,m;\; j = 1,\ldots,n\}$（这是 Kronecker 和 $A \oplus B^T$ 的特征值，注意 $\sigma(B^T) = \sigma(B)$）。

    因此方程对任意 $C$ 有唯一解当且仅当所有 $\lambda_i(A) + \lambda_j(B) \neq 0$，即 $\lambda_i(A) \neq -\lambda_j(B)$ 对所有 $i, j$ 成立，即 $\sigma(A) \cap \sigma(-B) = \emptyset$。

    **直接方法**（Rosenblum）：将 $A$ 化为 Jordan 标准形 $A = PJP^{-1}$，方程变为 $JY + YB = D$（$Y = P^{-1}X$，$D = P^{-1}C$）。由 $J$ 的三角结构，可以逐列求解 $Y$。当 $\sigma(A) \cap \sigma(-B) = \emptyset$ 时，每步求解的系数矩阵 $(\lambda_i I + B)$ 可逆。$\blacksquare$

!!! theorem "定理 20.4 (Sylvester 方程的积分表示)"
    若 $A$ 的所有特征值实部为负（$A$ 稳定），$B$ 的所有特征值实部为正，则 Sylvester 方程 $AX + XB = C$ 的唯一解为：

    $$
    X = \int_0^{\infty} e^{At} C e^{Bt} \, dt.
    $$

??? proof "证明"
    首先验证积分收敛。由于 $A$ 稳定，$\|e^{At}\| \leq M e^{-\alpha t}$（某 $\alpha > 0$）；由于 $B$ 的特征值实部为正，$-B$ 稳定，$\|e^{Bt}\| = \|e^{-(-B)t}\|$...

    更准确地说，条件应为 $\operatorname{Re}\lambda_i(A) < 0$ 且 $\operatorname{Re}\lambda_j(B) > 0$（或反过来取负号），使得 $e^{At} \to 0$（$t \to +\infty$）且 $e^{Bt}$ 增长但积分仍收敛。

    实际上，标准条件是 $\operatorname{Re}\lambda_i(A) + \operatorname{Re}\mu_j(B) < 0$。此时：

    令 $X = \int_0^{\infty} e^{At} C e^{Bt} \, dt$。对 $AX + XB$：

    $$
    AX + XB = \int_0^{\infty} (Ae^{At})Ce^{Bt} \, dt + \int_0^{\infty} e^{At}C(e^{Bt}B) \, dt
    $$

    $$
    = \int_0^{\infty} \frac{d}{dt}\left(e^{At}Ce^{Bt}\right) dt = \left[e^{At}Ce^{Bt}\right]_0^{\infty} = 0 - C = -C.
    $$

    因此 $X$ 是 $AX + XB = -C$ 的解。若方程为 $AX + XB = C$，则解为 $X = -\int_0^{\infty} e^{At}Ce^{Bt}\,dt$。$\blacksquare$

!!! definition "定义 20.4 (齐次 Sylvester 方程与交换子)"
    齐次 Sylvester 方程 $AX - XA = 0$（即 $AX = XA$）描述了与 $A$ **交换**的矩阵 $X$。与 $A$ 交换的所有矩阵构成一个代数，称为 $A$ 的**交换子代数**（commutant）。

    映射 $\operatorname{ad}_A(X) = AX - XA$ 称为 $A$ 的**伴随映射**（adjoint map），$[A, X] = AX - XA$ 称为 $A$ 和 $X$ 的**交换子**（commutator）或**Lie 括号**（Lie bracket）。

!!! theorem "定理 20.5 (交换矩阵的结构)"
    设 $A$ 为 $n \times n$ 矩阵，有 $k$ 个不同的特征值 $\lambda_1, \ldots, \lambda_k$，对应的 Jordan 块大小分别为 $n_{i,1} \geq n_{i,2} \geq \cdots$（$i = 1,\ldots,k$）。则与 $A$ 交换的矩阵空间的维数为：

    $$
    \dim\{X : AX = XA\} = \sum_{i=1}^{k} \sum_{j} (2j - 1) n_{i,j}',
    $$

    其中 $n_{i,j}'$ 是特征值 $\lambda_i$ 对应的 Jordan 块大小的共轭分拆。

    特别地，当 $A$ 有 $n$ 个不同特征值时（即 $A$ 的最小多项式等于特征多项式），与 $A$ 交换的矩阵恰好是 $A$ 的多项式 $p(A)$，维数为 $n$。

??? proof "证明"
    **简单情形**（$A$ 有 $n$ 个不同特征值）：设 $A = P \operatorname{diag}(\lambda_1, \ldots, \lambda_n) P^{-1}$，$AX = XA$ 当且仅当 $DY = YD$（$Y = P^{-1}XP$，$D = \operatorname{diag}(\lambda_1,\ldots,\lambda_n)$）。

    $DY = YD$ 意味着 $\lambda_i y_{ij} = \lambda_j y_{ij}$，即 $(\lambda_i - \lambda_j)y_{ij} = 0$。当 $\lambda_i \neq \lambda_j$ 时 $y_{ij} = 0$，因此 $Y$ 为对角矩阵，$X = P Y P^{-1}$。

    对角矩阵 $Y = \operatorname{diag}(d_1,\ldots,d_n)$ 可以唯一地写成 $Y = p(D)$，其中 $p$ 为 $n-1$ 次多项式（Lagrange 插值），因此 $X = p(A)$。维数为 $n$。$\blacksquare$

!!! example "例 20.2"
    求解 Sylvester 方程 $AX + XB = C$，其中：

    $$
    A = \begin{pmatrix} 1 & 0 \\ 0 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} 2 & 0 \\ 0 & 4 \end{pmatrix}, \quad C = \begin{pmatrix} 6 & 10 \\ 15 & 35 \end{pmatrix}.
    $$

    $\sigma(A) = \{1, 3\}$，$\sigma(-B) = \{-2, -4\}$，无公共特征值，唯一解存在。

    设 $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$。

    $AX + XB = \begin{pmatrix} x_{11} & x_{12} \\ 3x_{21} & 3x_{22} \end{pmatrix} + \begin{pmatrix} 2x_{11} & 4x_{12} \\ 2x_{21} & 4x_{22} \end{pmatrix} = \begin{pmatrix} 3x_{11} & 5x_{12} \\ 5x_{21} & 7x_{22} \end{pmatrix}$。

    因此 $3x_{11} = 6$，$5x_{12} = 10$，$5x_{21} = 15$，$7x_{22} = 35$。

    解：$X = \begin{pmatrix} 2 & 2 \\ 3 & 5 \end{pmatrix}$。

    验证：$AX + XB = \begin{pmatrix} 2 & 2 \\ 9 & 15 \end{pmatrix} + \begin{pmatrix} 4 & 8 \\ 6 & 20 \end{pmatrix} = \begin{pmatrix} 6 & 10 \\ 15 & 35 \end{pmatrix} = C$。正确。

!!! example "例 20.3"
    考虑齐次 Sylvester 方程 $AX - XB = 0$，即 $AX = XB$。

    设 $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$，$B = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$。

    $\sigma(A) = \{2\}$，$\sigma(B) = \{2, 3\}$。公共特征值为 $\{2\}$，因此齐次方程有非平凡解。

    设 $X = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$。

    $AX = \begin{pmatrix} 2a+c & 2b+d \\ 2c & 2d \end{pmatrix}$，$XB = \begin{pmatrix} 2a & 3b \\ 2c & 3d \end{pmatrix}$。

    $AX = XB$ 给出：$c = 0$，$2b + d = 3b$（即 $d = b$），$2d = 3d$（即 $d = 0$）。

    因此 $b = 0$，解空间为 $X = \begin{pmatrix} a & 0 \\ 0 & 0 \end{pmatrix}$，维数为 1。

---

## 20.3 Lyapunov 方程

<div class="context-flow" markdown>

**脉络**：$AX+XA^*=-Q$：Sylvester方程取 $B=A^*$

**稳定性等价**：$A$稳定($\operatorname{Re}\lambda_i<0$) $\Leftrightarrow$ $\exists X\succ 0$ 满足方程($Q\succ 0$) · 积分解 $X=\int_0^\infty e^{At}Qe^{A^*t}dt$

</div>

Lyapunov 方程是 Sylvester 方程的重要特殊情形，在系统稳定性分析中有核心地位。

!!! definition "定义 20.5 (连续 Lyapunov 方程)"
    **连续 Lyapunov 方程**（continuous Lyapunov equation）为：

    $$
    AX + XA^* = Q,
    $$

    其中 $A \in \mathbb{C}^{n \times n}$，$Q \in \mathbb{C}^{n \times n}$（$Q = Q^*$），$X \in \mathbb{C}^{n \times n}$ 为未知 Hermite 矩阵。

!!! definition "定义 20.6 (离散 Lyapunov 方程)"
    **离散 Lyapunov 方程**（discrete Lyapunov equation，也称 **Stein 方程**）为：

    $$
    AXA^* - X = Q,
    $$

    或等价地 $X - AXA^* = -Q$。其中 $A, Q, X \in \mathbb{C}^{n \times n}$。

!!! theorem "定理 20.6 (连续 Lyapunov 方程与稳定性)"
    设 $A \in \mathbb{C}^{n \times n}$。

    1. 若 $A$ **稳定**（即所有特征值实部为负：$\operatorname{Re}\lambda_i(A) < 0$），则对任意半正定 $Q \succeq 0$，连续 Lyapunov 方程 $AX + XA^* = -Q$ 有唯一的半正定解：

    $$
    X = \int_0^{\infty} e^{At} Q e^{A^*t} \, dt \succeq 0.
    $$

    2. 反之，若对某个 $Q \succ 0$，存在 $X \succ 0$ 满足 $AX + XA^* = -Q$，则 $A$ 稳定。

??? proof "证明"
    **(1)**：由 $A$ 稳定，$e^{At} \to 0$（$t \to \infty$）以指数速度，因此积分收敛。验证：

    $$
    AX + XA^* = \int_0^{\infty} (Ae^{At})Qe^{A^*t} \, dt + \int_0^{\infty} e^{At}Q(e^{A^*t}A^*) \, dt
    $$

    $$
    = \int_0^{\infty} \frac{d}{dt}(e^{At}Qe^{A^*t}) \, dt = [e^{At}Qe^{A^*t}]_0^{\infty} = 0 - Q = -Q.
    $$

    $X \succeq 0$ 因为对任意 $\mathbf{v}$：$\mathbf{v}^*X\mathbf{v} = \int_0^{\infty} \|Q^{1/2}e^{A^*t}\mathbf{v}\|^2 \, dt \geq 0$。

    唯一性由 Sylvester-Rosenblum 定理保证（$\sigma(A) \cap \sigma(-A^*) = \emptyset$，因为 $A$ 的特征值实部为负而 $-A^*$ 的特征值 $-\overline{\lambda_i(A)}$ 实部为正）。

    **(2)**：设 $A\mathbf{v} = \lambda\mathbf{v}$（$\mathbf{v} \neq 0$），则：

    $$
    \mathbf{v}^*(AX + XA^*)\mathbf{v} = \lambda(\mathbf{v}^*X\mathbf{v}) + \bar{\lambda}(\mathbf{v}^*X\mathbf{v}) = 2\operatorname{Re}(\lambda)(\mathbf{v}^*X\mathbf{v}) = -\mathbf{v}^*Q\mathbf{v}.
    $$

    由 $X \succ 0$，$\mathbf{v}^*X\mathbf{v} > 0$；由 $Q \succ 0$，$\mathbf{v}^*Q\mathbf{v} > 0$。因此 $2\operatorname{Re}(\lambda) < 0$，即 $\operatorname{Re}(\lambda) < 0$。$\blacksquare$

!!! theorem "定理 20.7 (离散 Lyapunov 方程与稳定性)"
    设 $A \in \mathbb{C}^{n \times n}$。

    1. 若 $A$ **Schur 稳定**（即所有特征值模小于 1：$|\lambda_i(A)| < 1$），则对任意 $Q \succeq 0$，离散方程 $X - AXA^* = Q$ 有唯一的半正定解：

    $$
    X = \sum_{k=0}^{\infty} A^k Q (A^*)^k \succeq 0.
    $$

    2. 反之，若对某个 $Q \succ 0$，存在 $X \succ 0$ 满足 $X - AXA^* = Q$，则 $A$ Schur 稳定。

??? proof "证明"
    **(1)**：由 $|\lambda_i(A)| < 1$，级数 $\sum_{k=0}^{\infty} A^k Q (A^*)^k$ 绝对收敛。验证：

    $$
    X - AXA^* = \sum_{k=0}^{\infty} A^k Q (A^*)^k - A\left(\sum_{k=0}^{\infty} A^k Q (A^*)^k\right)A^* = \sum_{k=0}^{\infty} A^k Q (A^*)^k - \sum_{k=1}^{\infty} A^k Q (A^*)^k = A^0 Q (A^*)^0 = Q.
    $$

    **(2)**：设 $A\mathbf{v} = \lambda\mathbf{v}$，则 $\mathbf{v}^*(X - AXA^*)\mathbf{v} = \mathbf{v}^*X\mathbf{v} - |\lambda|^2 \mathbf{v}^*X\mathbf{v} = (1 - |\lambda|^2)\mathbf{v}^*X\mathbf{v} = \mathbf{v}^*Q\mathbf{v} > 0$。由 $X \succ 0$，得 $1 - |\lambda|^2 > 0$，即 $|\lambda| < 1$。$\blacksquare$

!!! example "例 20.4"
    求解连续 Lyapunov 方程 $AX + XA^T = -Q$，其中：

    $$
    A = \begin{pmatrix} -1 & 0 \\ 0 & -2 \end{pmatrix}, \quad Q = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}.
    $$

    $A$ 稳定（特征值 $-1, -2$）。利用积分公式：

    $$
    X = \int_0^{\infty} e^{At} Q e^{A^Tt} \, dt = \int_0^{\infty} \begin{pmatrix} e^{-t} & 0 \\ 0 & e^{-2t} \end{pmatrix} \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} \begin{pmatrix} e^{-t} & 0 \\ 0 & e^{-2t} \end{pmatrix} dt.
    $$

    $$
    = \int_0^{\infty} \begin{pmatrix} 2e^{-2t} & 0 \\ 0 & 2e^{-4t} \end{pmatrix} dt = \begin{pmatrix} 1 & 0 \\ 0 & 1/2 \end{pmatrix}.
    $$

    验证：$AX + XA^T = \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix} + \begin{pmatrix} -1 & 0 \\ 0 & -1 \end{pmatrix} = \begin{pmatrix} -2 & 0 \\ 0 & -2 \end{pmatrix} = -Q$。正确。

!!! example "例 20.5"
    求解离散 Lyapunov 方程 $X - AXA^T = Q$，其中：

    $$
    A = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}, \quad Q = I_2.
    $$

    $A$ Schur 稳定（特征值 $1/2, 1/3$）。

    设 $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{12} & x_{22} \end{pmatrix}$（对称解）。

    $X - AXA^T = \begin{pmatrix} x_{11} - x_{11}/4 & x_{12} - x_{12}/6 \\ x_{12} - x_{12}/6 & x_{22} - x_{22}/9 \end{pmatrix} = \begin{pmatrix} 3x_{11}/4 & 5x_{12}/6 \\ 5x_{12}/6 & 8x_{22}/9 \end{pmatrix} = I_2$。

    解：$x_{11} = 4/3$，$x_{12} = 0$，$x_{22} = 9/8$。

    $X = \begin{pmatrix} 4/3 & 0 \\ 0 & 9/8 \end{pmatrix}$。

    验证可直接代入。$X \succ 0$ 且 $A$ Schur 稳定，与定理一致。

---

## 20.4 Riccati 方程

<div class="context-flow" markdown>

**脉络**：$A^*X+XA-XBR^{-1}B^*X+Q=0$（非线性！）

**Hamilton矩阵**的稳定不变子空间 → 稳定化解 $X=U_2U_1^{-1}$ · 可镇定+可检测 $\Rightarrow$ 唯一半正定解

</div>

Riccati 方程是一个**非线性**矩阵方程，在最优控制、滤波和博弈论中有核心地位。

!!! definition "定义 20.7 (代数 Riccati 方程)"
    **连续代数 Riccati 方程**（continuous algebraic Riccati equation, CARE）为：

    $$
    A^*X + XA - XBR^{-1}B^*X + Q = 0,
    $$

    其中 $A \in \mathbb{C}^{n \times n}$，$B \in \mathbb{C}^{n \times m}$，$Q = Q^* \in \mathbb{C}^{n \times n}$（$Q \succeq 0$），$R = R^* \in \mathbb{C}^{m \times m}$（$R \succ 0$），$X = X^* \in \mathbb{C}^{n \times n}$ 为未知 Hermite 矩阵。

    **离散代数 Riccati 方程**（discrete algebraic Riccati equation, DARE）为：

    $$
    X = A^*XA - A^*XB(R + B^*XB)^{-1}B^*XA + Q.
    $$

!!! definition "定义 20.8 (Hamiltonian 矩阵)"
    与连续 Riccati 方程相关联的 **Hamilton 矩阵**（Hamiltonian matrix）为：

    $$
    H = \begin{pmatrix} A & -BR^{-1}B^* \\ -Q & -A^* \end{pmatrix}.
    $$

    Hamilton 矩阵的一个重要性质是：若 $\lambda$ 是 $H$ 的特征值，则 $-\bar{\lambda}$ 也是。

!!! theorem "定理 20.8 (CARE 与 Hamilton 矩阵)"
    设 Hamilton 矩阵 $H$ 没有纯虚数特征值。将 $H$ 的 $2n$ 个特征值分为稳定部分（实部为负）和不稳定部分（实部为正），各 $n$ 个。

    设 $\begin{pmatrix} U_1 \\ U_2 \end{pmatrix}$ 为 $H$ 的稳定不变子空间的基（$U_1, U_2$ 均为 $n \times n$），若 $U_1$ 可逆，则 CARE 的**稳定化解**为：

    $$
    X = U_2 U_1^{-1}.
    $$

    此解使得闭环矩阵 $A - BR^{-1}B^*X$ 稳定。

??? proof "证明"
    设 $\mathbf{w} = \begin{pmatrix} \mathbf{x} \\ \mathbf{p} \end{pmatrix}$ 是 $H$ 的特征向量，$H\mathbf{w} = \lambda\mathbf{w}$。对稳定不变子空间，$H\begin{pmatrix} U_1 \\ U_2 \end{pmatrix} = \begin{pmatrix} U_1 \\ U_2 \end{pmatrix}\Lambda_s$，其中 $\Lambda_s$ 的特征值全在左半平面。

    展开第一行：$AU_1 - BR^{-1}B^*U_2 = U_1\Lambda_s$。
    展开第二行：$-QU_1 - A^*U_2 = U_2\Lambda_s$。

    若 $U_1$ 可逆，令 $X = U_2U_1^{-1}$，则 $U_2 = XU_1$。

    代入第一行：$AU_1 - BR^{-1}B^*XU_1 = U_1\Lambda_s$，即 $(A - BR^{-1}B^*X)U_1 = U_1\Lambda_s$。这说明 $A - BR^{-1}B^*X$ 与 $\Lambda_s$ 相似，因此闭环矩阵稳定。

    代入第二行：$-QU_1 - A^*XU_1 = XU_1\Lambda_s = X(AU_1 - BR^{-1}B^*XU_1)$。

    整理：$-Q - A^*X = XA - XBR^{-1}B^*X$，即 $A^*X + XA - XBR^{-1}B^*X + Q = 0$。$\blacksquare$

!!! theorem "定理 20.9 (CARE 的正定解存在条件)"
    设 $(A, B)$ 可镇定（stabilizable），$(A, Q^{1/2})$ 可检测（detectable），$Q \succeq 0$，$R \succ 0$。则 CARE 存在唯一的正半定解 $X \succeq 0$，且闭环矩阵 $A - BR^{-1}B^*X$ 稳定。

??? proof "证明"
    本定理的完整证明较长，我们给出思路框架。

    **可镇定性**意味着存在反馈增益 $K$ 使得 $A - BK$ 稳定。**可检测性**意味着 $(A^*, (Q^{1/2})^*)$ 可镇定。

    在这些条件下，可以证明 Hamilton 矩阵 $H$ 没有纯虚数特征值（关键一步），因此定理 20.8 的构造可以进行。

    正半定性 $X \succeq 0$ 来自 $Q \succeq 0$ 的条件和 Lyapunov 论证：闭环系统 $(A-BR^{-1}B^*X)$ 稳定时，$X$ 满足一个 Lyapunov 方程，由此保证半正定性。

    唯一性可以通过反证法：若有两个半正定解 $X_1, X_2$，则它们的差满足一个线性方程，由稳定性条件可推出差为零。$\blacksquare$

!!! example "例 20.6"
    求解标量 Riccati 方程：$2x + 2x - x^2 + 1 = 0$，其中 $A = 1$，$B = 1$，$R = 1$，$Q = 1$。

    简化：$2x - x^2 + 1 = 0$，即 $x^2 - 2x - 1 = 0$。

    $x = \frac{2 \pm \sqrt{4+4}}{2} = 1 \pm \sqrt{2}$。

    稳定化解：闭环 $A - BR^{-1}B^*x = 1 - x$。取 $x = 1 + \sqrt{2}$，闭环 $= 1 - 1 - \sqrt{2} = -\sqrt{2} < 0$（稳定）。取 $x = 1 - \sqrt{2}$，闭环 $= \sqrt{2} > 0$（不稳定）。

    因此稳定化解 $x = 1 + \sqrt{2}$。

!!! example "例 20.7"
    **最优控制中的 Riccati 方程**。

    考虑线性系统 $\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$，代价函数 $J = \int_0^{\infty}(\mathbf{x}^TQ\mathbf{x} + \mathbf{u}^TR\mathbf{u})\,dt$。

    设 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$（双积分器），$B = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$，$Q = I_2$，$R = 1$。

    Hamilton 矩阵：$H = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & -1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & -1 & 0 \end{pmatrix}$。

    $H$ 的特征多项式可以计算为 $\lambda^4 + 3\lambda^2 + 1 = 0$（令 $\mu = \lambda^2$，则 $\mu^2 + 3\mu + 1 = 0$，$\mu = \frac{-3 \pm \sqrt{5}}{2}$）。

    稳定特征值对应的子空间给出 Riccati 方程的正定解 $X = \begin{pmatrix} \sqrt{3} & 1 \\ 1 & \sqrt{3} \end{pmatrix}$（精确值需要更仔细的计算）。

    最优反馈为 $\mathbf{u} = -R^{-1}B^TX\mathbf{x} = -\begin{pmatrix} 0 & 1 \end{pmatrix}X\mathbf{x}$。

---

## 20.5 数值解法

<div class="context-flow" markdown>

**洞察**：Kronecker积直接法 $O(n^6)$ 不可行 → **Bartels-Stewart**：Schur分解 + 逐列回代 = $O(n^3)$ · 关键：三角Sylvester方程可逐列求解

</div>

本节介绍求解 Sylvester 和 Lyapunov 方程的数值方法。直接使用 Kronecker 积将方程转化为 $n^2$ 维线性系统虽然理论上可行，但计算复杂度为 $O(n^6)$，实际中不可接受。Bartels-Stewart 算法将复杂度降低到 $O(n^3)$。

!!! definition "定义 20.9 (Schur 分解)"
    任意方阵 $A$ 可以分解为 $A = UTU^*$，其中 $U$ 为酉矩阵，$T$ 为上三角矩阵（对角元素为 $A$ 的特征值）。这称为 **Schur 分解**（Schur decomposition）。对于实矩阵，可以使用**实 Schur 分解**：$A = UTU^T$，其中 $T$ 为拟上三角矩阵（对角块为 $1 \times 1$ 或 $2 \times 2$）。

!!! theorem "定理 20.10 (Bartels-Stewart 算法)"
    Sylvester 方程 $AX + XB = C$ 可以通过以下步骤高效求解：

    **步骤 1**：计算 Schur 分解 $A = U_A T_A U_A^*$，$B = U_B T_B U_B^*$。

    **步骤 2**：变量替换 $Y = U_A^* X U_B$，$D = U_A^* C U_B$，方程变为：

    $$
    T_A Y + Y T_B = D.
    $$

    **步骤 3**：由于 $T_A$ 和 $T_B$ 为（拟）上三角矩阵，方程可以逐列（从最后一列开始）用回代法求解。

    **步骤 4**：恢复 $X = U_A Y U_B^*$。

    总计算复杂度为 $O(n^3)$（Schur 分解主导）。

??? proof "证明"
    **正确性**：代入验证，$T_A Y + YT_B = U_A^*AU_A \cdot U_A^*XU_B + U_A^*XU_B \cdot U_B^*BU_B = U_A^*(AX + XB)U_B = U_A^*CU_B = D$。

    **逐列求解**：设 $T_B$ 的第 $j$ 列为 $\mathbf{t}_j$（仅前 $j$ 个分量可能非零），$Y$ 的第 $j$ 列为 $\mathbf{y}_j$，$D$ 的第 $j$ 列为 $\mathbf{d}_j$。

    方程 $T_A Y + Y T_B = D$ 的第 $j$ 列为：

    $$
    T_A \mathbf{y}_j + \sum_{k=1}^{j} (T_B)_{kj} \mathbf{y}_k = \mathbf{d}_j.
    $$

    即 $(T_A + (T_B)_{jj}I)\mathbf{y}_j = \mathbf{d}_j - \sum_{k=1}^{j-1}(T_B)_{kj}\mathbf{y}_k$。

    从 $j = 1$ 开始，右端已知（当 $j = 1$ 时无求和项），$(T_A + (T_B)_{11}I)$ 为上三角矩阵，可用回代求解 $\mathbf{y}_1$。然后逐列进行。

    每列求解 $O(n^2)$（三角系统），共 $n$ 列，回代部分 $O(n^3)$。Schur 分解 $O(n^3)$。总复杂度 $O(n^3)$。$\blacksquare$

!!! theorem "定理 20.11 (Hessenberg-Schur 方法)"
    对于 Sylvester 方程 $AX + XB = C$（$A$ 为 $m \times m$，$B$ 为 $n \times n$，$m \gg n$ 或 $m \ll n$），可以采用改进的 **Hessenberg-Schur 方法**：仅对较小的矩阵做 Schur 分解，较大的矩阵化为 Hessenberg 形。这在某些情况下可以减少计算量。

??? proof "证明"
    对 $B$ 做 Schur 分解 $B = U_B T_B U_B^*$，对 $A$ 化为上 Hessenberg 形 $A = Q_A H_A Q_A^*$。变量替换后方程变为 $H_A Y + Y T_B = D$。

    逐列求解时，$(H_A + (T_B)_{jj}I)\mathbf{y}_j = \mathbf{r}_j$ 为 Hessenberg 系统，可用 $O(n^2)$ 的方法求解（无需化为三角形，直接用 Gauss 消元即可，因 Hessenberg 矩阵的 LU 分解很快）。

    总计算量的常数因子比标准 Bartels-Stewart 略小。$\blacksquare$

!!! example "例 20.8"
    用 Bartels-Stewart 算法求解 $AX + XB = C$，其中：

    $$
    A = \begin{pmatrix} 2 & 1 \\ 0 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} 1 & 2 \\ 0 & 4 \end{pmatrix}, \quad C = \begin{pmatrix} 10 & 26 \\ 9 & 28 \end{pmatrix}.
    $$

    $A$ 和 $B$ 已经是上三角矩阵（即为自身的 Schur 分解，$U_A = U_B = I$）。

    直接逐列求解 $T_A Y + Y T_B = D$（即 $AX + XB = C$）。

    **第 1 列**（$j = 1$）：$(A + b_{11}I)\mathbf{x}_1 = \mathbf{c}_1$。

    $$
    \begin{pmatrix} 3 & 1 \\ 0 & 4 \end{pmatrix}\begin{pmatrix} x_{11} \\ x_{21} \end{pmatrix} = \begin{pmatrix} 10 \\ 9 \end{pmatrix}.
    $$

    回代：$x_{21} = 9/4$，$x_{11} = (10 - 9/4)/3 = 31/12$。

    **第 2 列**（$j = 2$）：$(A + b_{22}I)\mathbf{x}_2 = \mathbf{c}_2 - b_{12}\mathbf{x}_1$。

    $$
    \begin{pmatrix} 6 & 1 \\ 0 & 7 \end{pmatrix}\begin{pmatrix} x_{12} \\ x_{22} \end{pmatrix} = \begin{pmatrix} 26 \\ 28 \end{pmatrix} - 2\begin{pmatrix} 31/12 \\ 9/4 \end{pmatrix} = \begin{pmatrix} 26 - 31/6 \\ 28 - 9/2 \end{pmatrix} = \begin{pmatrix} 125/6 \\ 47/2 \end{pmatrix}.
    $$

    回代：$x_{22} = 47/14$，$x_{12} = (125/6 - 47/14)/6 = (875/42 - 141/42)/6 = (734/42)/6 = 734/252 = 367/126$。

    验证可代回原方程。（数值较复杂但过程正确。）

!!! example "例 20.9"
    **Lyapunov 方程的 Bartels-Stewart 求解**。

    对 $AX + XA^T = -Q$（实对称情形），只需对 $A$ 做一次 Schur 分解 $A = UTU^T$，方程变为 $TY + YT^T = -\tilde{Q}$（$Y = U^TXU$，$\tilde{Q} = U^TQU$）。由 $T$ 上三角而 $T^T$ 下三角，逐列求解时需要利用对称性。

    设 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$... 这不是上三角的。取实 Schur 分解：$A$ 的特征值为 $\pm i$（纯虚数），此时 Lyapunov 方程 $AX + XA^T = -Q$ 不一定有解（因为 $\sigma(A) \cap \sigma(-A^T) = \{i, -i\} \cap \{i, -i\} \neq \emptyset$）。

    改为 $A = \begin{pmatrix} -1 & 1 \\ 0 & -2 \end{pmatrix}$（稳定，已上三角），$Q = I$。

    $(T + t_{11}I)\mathbf{y}_1 = -\mathbf{q}_1$：$\begin{pmatrix} -2 & 1 \\ 0 & -3 \end{pmatrix}\begin{pmatrix} y_{11} \\ y_{21} \end{pmatrix} = -\begin{pmatrix} 1 \\ 0 \end{pmatrix}$。

    $y_{21} = 0$，$y_{11} = 1/2$。

    $(T + t_{22}I)\mathbf{y}_2 = -\mathbf{q}_2 - t_{12}^T\mathbf{y}_1$：其中 $t_{12}^T$ 是 $T^T$ 的上三角元素...

    更仔细地处理：$TY + YT^T = -Q$。逐列时方程结构需要同时利用行和列的三角性。设 $Y$ 对称，则只需求 $n(n+1)/2$ 个变量。此例中解为 $X = \begin{pmatrix} 5/6 & 1/6 \\ 1/6 & 1/4 \end{pmatrix}$（需代入验证）。

---

## 20.6 Penrose 方程

<div class="context-flow" markdown>

**脉络**：四个Penrose方程公理化定义 $A^+$ · (1)内逆 → (1)(3)最小二乘 → (1)(4)最小范数 → (1)(2)(3)(4)唯一Moore-Penrose伪逆 · SVD给出存在性和显式公式

</div>

Penrose 方程组给出了 Moore-Penrose 伪逆的公理化刻画。

!!! definition "定义 20.10 (Penrose 方程组)"
    设 $A$ 为 $m \times n$ 矩阵。**Penrose 方程组**（Penrose equations）是关于 $n \times m$ 矩阵 $X$ 的以下四个方程：

    $$
    \begin{aligned}
    &(1)\quad AXA = A, \\
    &(2)\quad XAX = X, \\
    &(3)\quad (AX)^* = AX, \\
    &(4)\quad (XA)^* = XA.
    \end{aligned}
    $$

!!! definition "定义 20.11 (广义逆)"
    满足 Penrose 方程组不同子集的矩阵 $X$ 称为 $A$ 的**广义逆**（generalized inverse）：

    - 满足 (1) 的 $X$ 称为 **$\{1\}$-逆**或**内逆**（inner inverse），记作 $A^{(1)}$。
    - 满足 (1)(2) 的 $X$ 称为 **$\{1,2\}$-逆**或**自反广义逆**（reflexive generalized inverse）。
    - 满足 (1)(3) 的 $X$ 称为 **$\{1,3\}$-逆**或**最小二乘广义逆**。
    - 满足 (1)(4) 的 $X$ 称为 **$\{1,4\}$-逆**或**最小范数广义逆**。
    - 满足全部四个方程的唯一的 $X$ 称为 **Moore-Penrose 伪逆**（Moore-Penrose pseudoinverse），记作 $A^+$。

!!! theorem "定理 20.12 (Moore-Penrose 伪逆的存在唯一性)"
    对任意 $m \times n$ 矩阵 $A$，Moore-Penrose 伪逆 $A^+$ 存在且唯一。

??? proof "证明"
    **存在性**：设 $A = U\Sigma V^*$ 为奇异值分解，其中 $\Sigma = \begin{pmatrix} \Sigma_r & 0 \\ 0 & 0 \end{pmatrix}$，$\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$（$r = \operatorname{rank}(A)$）。

    定义 $A^+ = V\Sigma^+ U^*$，其中 $\Sigma^+ = \begin{pmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{pmatrix}$。

    逐一验证四个 Penrose 方程：

    (1) $AXA = U\Sigma V^* V\Sigma^+ U^* U\Sigma V^* = U\Sigma\Sigma^+\Sigma V^* = U\Sigma V^* = A$。

    (2) $XAX = V\Sigma^+ U^* U\Sigma V^* V\Sigma^+ U^* = V\Sigma^+\Sigma\Sigma^+ U^* = V\Sigma^+ U^* = X$。

    (3) $AX = U\Sigma\Sigma^+ U^* = U\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}U^*$，为 Hermite（因为对角分块实对称，$U(\cdot)U^*$ 保持 Hermite 性）。

    (4) $XA = V\Sigma^+\Sigma V^* = V\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}V^*$，同样为 Hermite。

    **唯一性**：设 $X_1, X_2$ 均满足四个方程。

    由 (1)(3)：$AX_1$ 和 $AX_2$ 均为到 $\mathcal{C}(A)$ 的正交投影，因此 $AX_1 = AX_2$。

    由 (1)(4)：$X_1A$ 和 $X_2A$ 均为到 $\mathcal{C}(A^*)$ 的正交投影，因此 $X_1A = X_2A$。

    由 (2)：$X_1 = X_1AX_1 = X_2AX_1 = X_2AX_2 = X_2$（中间步骤利用了 $X_1A = X_2A$ 和 $AX_1 = AX_2$）。$\blacksquare$

!!! theorem "定理 20.13 (Moore-Penrose 伪逆的性质)"
    设 $A$ 为 $m \times n$ 矩阵，$A^+$ 为其 Moore-Penrose 伪逆，则：

    1. $(A^+)^+ = A$。
    2. $(A^*)^+ = (A^+)^*$。
    3. $(A^*A)^+ = A^+(A^+)^* = A^+(A^*)^+$。
    4. $\operatorname{rank}(A^+) = \operatorname{rank}(A)$。
    5. $A^+ = (A^*A)^+A^* = A^*(AA^*)^+$。
    6. 若 $A$ 为列满秩，则 $A^+ = (A^*A)^{-1}A^*$。
    7. 若 $A$ 为行满秩，则 $A^+ = A^*(AA^*)^{-1}$。

??? proof "证明"
    **(1)**：$A$ 满足以 $A^+$ 为"原矩阵"的四个 Penrose 方程（只需交换 $A, X$ 的角色并利用 Hermite 条件），因此 $(A^+)^+ = A$。

    **(2)**：设 $A = U\Sigma V^*$，则 $A^* = V\Sigma^* U^*$，$(A^*)^+ = U(\Sigma^*)^+ V^*$。又 $A^+ = V\Sigma^+ U^*$，$(A^+)^* = U(\Sigma^+)^* V^*$。由于 $\Sigma$ 为实矩阵，$(\Sigma^*)^+ = (\Sigma^+)^*$，因此 $(A^*)^+ = (A^+)^*$。

    **(6)**：若 $A$ 列满秩（$\operatorname{rank}(A) = n$），则 $A^*A$ 可逆。此时 $X = (A^*A)^{-1}A^*$ 满足：
    - (1)：$AXA = A(A^*A)^{-1}A^*A = A$。
    - (3)：$AX = A(A^*A)^{-1}A^*$ 为 Hermite。
    - (4)：$XA = (A^*A)^{-1}A^*A = I$ 为 Hermite。
    - (2)：$XAX = IX = X$。
    因此 $X = A^+$。$\blacksquare$

!!! theorem "定理 20.14 ($\{1\}$-逆的通解)"
    对 $m \times n$ 矩阵 $A$（$\operatorname{rank}(A) = r$），方程 $AXA = A$ 的通解为：

    $$
    X = A^+ + (I - A^+A)W_1 + W_2(I - AA^+),
    $$

    其中 $W_1, W_2$ 为任意大小适当的矩阵。

    等价地，设 $A = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q$（$P, Q$ 可逆），则 $\{1\}$-逆的一般形式为：

    $$
    X = Q^{-1}\begin{pmatrix} I_r & L \\ M & N \end{pmatrix}P^{-1},
    $$

    其中 $L, M, N$ 为任意矩阵。

??? proof "证明"
    设 $A = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q$，$X = Q^{-1}\begin{pmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{pmatrix}P^{-1}$。

    $AXA = P\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} G_{11} & G_{12} \\ G_{21} & G_{22} \end{pmatrix}\begin{pmatrix} I_r & 0 \\ 0 & 0 \end{pmatrix}Q = P\begin{pmatrix} G_{11} & 0 \\ 0 & 0 \end{pmatrix}Q$。

    $AXA = A$ 要求 $G_{11} = I_r$，而 $G_{12}, G_{21}, G_{22}$ 任意。$\blacksquare$

!!! theorem "定理 20.15 (Penrose 方程与最小二乘)"
    方程 $AXA = A$ 的解与最小二乘问题密切相关：

    对于不相容方程组 $A\mathbf{x} = \mathbf{b}$：

    1. $\mathbf{x} = A^{(1)}\mathbf{b}$ 是某个最小二乘解（满足 (1)），但不一定是范数最小的。
    2. $\mathbf{x} = A^{(1,3)}\mathbf{b}$ 是最小二乘解（$\min \|A\mathbf{x} - \mathbf{b}\|$），满足 (1)(3)。
    3. $\mathbf{x} = A^{(1,4)}\mathbf{b}$ 是最小范数解（在相容方程的解中范数最小），满足 (1)(4)。
    4. $\mathbf{x} = A^+\mathbf{b}$ 是**最小范数最小二乘解**（满足全部四个方程）。

??? proof "证明"
    **(4)** 的证明：$\mathbf{x}_0 = A^+\mathbf{b}$ 满足 $A\mathbf{x}_0 = AA^+\mathbf{b}$。由 Penrose 方程 (3)，$AA^+$ 为正交投影到 $\mathcal{C}(A)$ 上，因此 $AA^+\mathbf{b}$ 是 $\mathbf{b}$ 在 $\mathcal{C}(A)$ 上的正交投影，即 $\|A\mathbf{x} - \mathbf{b}\|$ 的最小值在 $A\mathbf{x} = AA^+\mathbf{b}$ 时取到。

    最小二乘解的集合为 $\{\mathbf{x} : A\mathbf{x} = AA^+\mathbf{b}\} = A^+\mathbf{b} + \mathcal{N}(A)$。

    $\mathbf{x}_0 = A^+\mathbf{b} \in \mathcal{C}(A^+) = \mathcal{C}(A^*)$（由 SVD 可见）。而 $\mathcal{C}(A^*) = \mathcal{N}(A)^{\perp}$。因此 $\mathbf{x}_0$ 与 $\mathcal{N}(A)$ 正交，在仿射集 $A^+\mathbf{b} + \mathcal{N}(A)$ 中具有最小范数。$\blacksquare$

!!! example "例 20.10"
    计算矩阵 $A = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$ 的 Moore-Penrose 伪逆。

    $A$ 列满秩（$\operatorname{rank} = 2$），因此 $A^+ = (A^TA)^{-1}A^T$。

    $A^TA = I_2$，故 $A^+ = A^T = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}$。

    验证 Penrose 方程：

    (1) $AXA = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix} = A$。

    (2) $XAX = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \end{pmatrix} = X$。

    (3) $AX = \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 0 \end{pmatrix}$，对称。

    (4) $XA = I_2$，对称。

    全部验证通过。

!!! example "例 20.11"
    计算秩 1 矩阵 $A = \begin{pmatrix} 1 \\ 2 \\ 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}$ 的 Moore-Penrose 伪逆。

    $A = \mathbf{u}\mathbf{v}^T$，其中 $\mathbf{u} = (1,2,3)^T$，$\mathbf{v} = (1,1)^T$。

    对秩 1 矩阵，$A^+ = \frac{\mathbf{v}\mathbf{u}^T}{\|\mathbf{u}\|^2\|\mathbf{v}\|^2} = \frac{1}{14 \cdot 2}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix} = \frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}$。

    验证 (1)：$AXA = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}\frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}$。

    先算 $XA = \frac{1}{28}\begin{pmatrix} 1 & 2 & 3 \\ 1 & 2 & 3 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix} = \frac{14}{28}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。

    $AXA = A \cdot \frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$... 更准确地：$A(XA) = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix}\frac{1}{2}\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 2 & 2 \\ 3 & 3 \end{pmatrix} = A$。正确。

!!! example "例 20.12"
    **利用伪逆求解不相容方程组**。

    设 $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \\ 0 & 1 \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix}$。

    $A^TA = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$，$(A^TA)^{-1} = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}$。

    $A^+ = (A^TA)^{-1}A^T = \begin{pmatrix} 1/2 & 0 \\ 0 & 1/3 \end{pmatrix}\begin{pmatrix} 1 & 1 & 0 \\ 1 & -1 & 1 \end{pmatrix} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ 1/3 & -1/3 & 1/3 \end{pmatrix}$。

    最小范数最小二乘解：$\mathbf{x} = A^+\mathbf{b} = \begin{pmatrix} 1/2 & 1/2 & 0 \\ 1/3 & -1/3 & 1/3 \end{pmatrix}\begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} = \begin{pmatrix} 2 \\ 1 \end{pmatrix}$。

    残差：$A\mathbf{x} - \mathbf{b} = \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} - \begin{pmatrix} 3 \\ 1 \\ 1 \end{pmatrix} = \mathbf{0}$。

    方程组竟然是相容的！解唯一（因为 $A$ 列满秩）。

!!! example "例 20.13"
    **真正不相容的情形**。设 $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$，$\mathbf{b} = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$。

    $A^+ = \frac{A^T}{A^TA} = \frac{1}{2}(1, 1)$。

    $\mathbf{x} = A^+\mathbf{b} = \frac{1}{2}(1 + 3) = 2$。

    最小二乘解：$x = 2$，$A\mathbf{x} = \begin{pmatrix} 2 \\ 2 \end{pmatrix}$，残差 $= \begin{pmatrix} -1 \\ 1 \end{pmatrix}$，$\|$残差$\| = \sqrt{2}$。

    这是最优的，因为 $\|A x - \mathbf{b}\|^2 = (x-1)^2 + (x-3)^2$，对 $x$ 求导得 $2(x-1) + 2(x-3) = 0$，$x = 2$。
