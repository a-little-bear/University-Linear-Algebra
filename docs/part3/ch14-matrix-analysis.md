# 第 14 章 矩阵分析

矩阵分析是线性代数与分析学的交汇领域，它将极限、连续性、微分等分析工具引入矩阵空间，为研究矩阵的渐近行为、函数性质和谱结构提供了强大的理论框架。本章从矩阵空间的拓扑结构出发，依次讨论矩阵级数的收敛性、谱半径的刻画、特征值的定位（Gershgorin 圆盘定理）、矩阵微积分以及特征值与奇异值之间的深层联系。这些内容构成了数值线性代数、控制理论和优化理论的数学基础。

---

## 14.1 矩阵空间的拓扑

矩阵空间 $\mathbb{C}^{m \times n}$（或 $\mathbb{R}^{m \times n}$）可以自然地视为有限维赋范线性空间。通过在该空间上引入范数，我们可以讨论矩阵序列的收敛、矩阵函数的连续性等分析概念。

!!! definition "定义 14.1 (矩阵序列的收敛)"
    设 $\{A_k\}_{k=1}^{\infty}$ 是 $\mathbb{C}^{m \times n}$ 中的矩阵序列，$A \in \mathbb{C}^{m \times n}$。称序列 $\{A_k\}$ **收敛（convergent）** 于 $A$，记作 $\lim_{k \to \infty} A_k = A$，若对 $\mathbb{C}^{m \times n}$ 上的某个（从而任意一个）矩阵范数 $\|\cdot\|$，有

    $$\lim_{k \to \infty} \|A_k - A\| = 0.$$

!!! definition "定义 14.2 (逐元素收敛)"
    矩阵序列 $\{A_k\}$ **逐元素收敛（entrywise convergence）** 于 $A$，是指对每个 $i, j$，有 $\lim_{k \to \infty} (A_k)_{ij} = A_{ij}$。

!!! theorem "定理 14.1 (收敛等价性)"
    设 $\{A_k\} \subset \mathbb{C}^{m \times n}$，$A \in \mathbb{C}^{m \times n}$。则以下条件等价：

    (1) $\{A_k\}$ 在某个矩阵范数下收敛于 $A$；

    (2) $\{A_k\}$ 在任意矩阵范数下收敛于 $A$；

    (3) $\{A_k\}$ 逐元素收敛于 $A$。

??? proof "证明"
    **(1) $\Rightarrow$ (2)**：由于 $\mathbb{C}^{m \times n}$ 是有限维空间，其上所有范数等价。设 $\|\cdot\|_\alpha$ 和 $\|\cdot\|_\beta$ 是两个范数，则存在常数 $c_1, c_2 > 0$ 使得 $c_1\|X\|_\alpha \leq \|X\|_\beta \leq c_2\|X\|_\alpha$ 对所有 $X$ 成立。若 $\|A_k - A\|_\alpha \to 0$，则 $\|A_k - A\|_\beta \leq c_2\|A_k - A\|_\alpha \to 0$。

    **(2) $\Rightarrow$ (3)**：取 $\|\cdot\|$ 为 Frobenius 范数，则 $|(A_k)_{ij} - A_{ij}| \leq \|A_k - A\|_F \to 0$。

    **(3) $\Rightarrow$ (1)**：若逐元素收敛成立，则取 $\|\cdot\|_{\max} = \max_{i,j}|a_{ij}|$（这是一个范数），有 $\|A_k - A\|_{\max} = \max_{i,j}|(A_k)_{ij} - A_{ij}| \to 0$。由范数等价性，在任意范数下也收敛。$\blacksquare$

!!! definition "定义 14.3 (矩阵空间的赋范空间结构)"
    $\mathbb{C}^{m \times n}$ 配备矩阵范数 $\|\cdot\|$ 后构成一个**赋范线性空间（normed linear space）**。由于 $\mathbb{C}^{m \times n}$ 同构于 $\mathbb{C}^{mn}$（作为向量空间），它是有限维的，因此是**完备的（complete）**，即构成**Banach 空间（Banach space）**。

!!! theorem "定理 14.2 (矩阵空间的紧性)"
    $\mathbb{C}^{m \times n}$ 中的有界闭集是紧集。特别地，每个有界矩阵序列都有收敛子列。

??? proof "证明"
    $\mathbb{C}^{m \times n}$ 同构于 $\mathbb{C}^{mn}$，后者是有限维赋范空间。由 Heine-Borel 定理，有限维赋范空间中的有界闭集是紧集。有界序列位于某个闭球中，由紧性可提取收敛子列。$\blacksquare$

!!! example "例 14.1"
    考虑矩阵序列 $A_k = \begin{pmatrix} 1/k & 2/k \\ 0 & 1/k^2 \end{pmatrix}$。

    逐元素地，$(A_k)_{11} = 1/k \to 0$，$(A_k)_{12} = 2/k \to 0$，$(A_k)_{21} = 0 \to 0$，$(A_k)_{22} = 1/k^2 \to 0$。

    因此 $\lim_{k \to \infty} A_k = O$（零矩阵）。

    在 Frobenius 范数下验证：$\|A_k\|_F = \sqrt{1/k^2 + 4/k^2 + 0 + 1/k^4} = \sqrt{5/k^2 + 1/k^4} \to 0$。

!!! example "例 14.2"
    设 $A = \begin{pmatrix} 0.5 & 0.1 \\ 0 & 0.3 \end{pmatrix}$，考虑序列 $\{A^k\}$。

    由于 $A$ 是上三角矩阵，其特征值为对角元素 $\lambda_1 = 0.5$，$\lambda_2 = 0.3$，均满足 $|\lambda_i| < 1$。

    直接计算 $A^2 = \begin{pmatrix} 0.25 & 0.08 \\ 0 & 0.09 \end{pmatrix}$，$A^3 = \begin{pmatrix} 0.125 & 0.049 \\ 0 & 0.027 \end{pmatrix}$。

    可以观察到各元素逐渐趋于 $0$，因此 $A^k \to O$。这与后面 14.7 节的理论一致。

---

## 14.2 矩阵级数

类比数值级数，我们可以定义矩阵级数及其收敛性。最重要的矩阵级数之一是 Neumann 级数，它给出了 $(I - A)^{-1}$ 的级数表示。

!!! definition "定义 14.4 (矩阵级数的收敛)"
    设 $\{A_k\}_{k=0}^{\infty} \subset \mathbb{C}^{n \times n}$。矩阵级数 $\sum_{k=0}^{\infty} A_k$ **收敛（convergent）**，是指其部分和序列 $S_N = \sum_{k=0}^{N} A_k$ 收敛。若 $\sum_{k=0}^{\infty} \|A_k\|$ 收敛（其中 $\|\cdot\|$ 为某个矩阵范数），则称级数**绝对收敛（absolutely convergent）**。

!!! theorem "定理 14.3 (绝对收敛蕴含收敛)"
    在 $\mathbb{C}^{n \times n}$ 中，绝对收敛的矩阵级数必收敛。

??? proof "证明"
    设 $\sum_{k=0}^{\infty}\|A_k\|$ 收敛。对 $N > M$，有

    $$\left\|\sum_{k=M+1}^{N} A_k\right\| \leq \sum_{k=M+1}^{N}\|A_k\| \to 0 \quad (M, N \to \infty).$$

    因此部分和 $\{S_N\}$ 是 Cauchy 序列。由 $\mathbb{C}^{n \times n}$ 的完备性，该序列收敛。$\blacksquare$

!!! definition "定义 14.5 (Neumann 级数)"
    设 $A \in \mathbb{C}^{n \times n}$，级数 $\sum_{k=0}^{\infty} A^k$ 称为 **Neumann 级数（Neumann series）**。

!!! theorem "定理 14.4 (Neumann 级数收敛定理)"
    设 $A \in \mathbb{C}^{n \times n}$。以下条件等价：

    (1) Neumann 级数 $\sum_{k=0}^{\infty} A^k$ 收敛；

    (2) $\rho(A) < 1$（谱半径小于 $1$）；

    (3) $A^k \to O$（$k \to \infty$）。

    当上述条件成立时，$I - A$ 可逆，且

    $$\sum_{k=0}^{\infty} A^k = (I - A)^{-1}.$$

??? proof "证明"
    **(2) $\Rightarrow$ (3)**：见定理 14.10。

    **(3) $\Rightarrow$ (1)**：设 $S_N = \sum_{k=0}^{N} A^k$。注意到

    $$(I - A)S_N = S_N(I - A) = I - A^{N+1}.$$

    由 $A^{N+1} \to O$，右端趋于 $I$。这说明当 $N$ 充分大时 $I - A^{N+1}$ 可逆，从而 $I - A$ 可逆（因为 $I - A$ 不依赖 $N$，它要么可逆要么不可逆）。

    更精确地，由 $(I-A)S_N = I - A^{N+1}$，若 $I - A$ 可逆则 $S_N = (I-A)^{-1}(I - A^{N+1}) \to (I-A)^{-1}$。

    下证 $I - A$ 确实可逆：若 $I - A$ 不可逆，则 $1$ 是 $A$ 的特征值，设 $A\mathbf{v} = \mathbf{v}$（$\mathbf{v} \neq \mathbf{0}$），则 $A^k \mathbf{v} = \mathbf{v}$ 对所有 $k$ 成立，矛盾于 $A^k \to O$。

    **(1) $\Rightarrow$ (2)**：若 $\sum A^k$ 收敛，则 $A^k \to O$（级数收敛的必要条件）。若存在特征值 $\lambda$ 满足 $|\lambda| \geq 1$，设 $A\mathbf{v} = \lambda\mathbf{v}$，则 $A^k\mathbf{v} = \lambda^k \mathbf{v}$，$\|A^k\mathbf{v}\| = |\lambda|^k\|\mathbf{v}\| \geq \|\mathbf{v}\| > 0$，矛盾于 $A^k \to O$。故 $\rho(A) < 1$。$\blacksquare$

!!! example "例 14.3"
    设 $A = \begin{pmatrix} 0 & 1/2 \\ 1/3 & 0 \end{pmatrix}$。

    特征多项式为 $\lambda^2 - 1/6 = 0$，特征值 $\lambda = \pm 1/\sqrt{6}$。因此 $\rho(A) = 1/\sqrt{6} < 1$，Neumann 级数收敛。

    $$(I - A)^{-1} = \begin{pmatrix} 1 & -1/2 \\ -1/3 & 1 \end{pmatrix}^{-1} = \frac{1}{1-1/6}\begin{pmatrix} 1 & 1/2 \\ 1/3 & 1 \end{pmatrix} = \frac{6}{5}\begin{pmatrix} 1 & 1/2 \\ 1/3 & 1 \end{pmatrix} = \begin{pmatrix} 6/5 & 3/5 \\ 2/5 & 6/5 \end{pmatrix}.$$

!!! theorem "定理 14.5 (Neumann 级数的范数估计)"
    设 $\|\cdot\|$ 为 $\mathbb{C}^{n \times n}$ 上的次可乘范数，且 $\|A\| < 1$。则 $I - A$ 可逆，且

    $$\|(I - A)^{-1}\| \leq \frac{1}{1 - \|A\|}, \qquad \|(I-A)^{-1} - I\| \leq \frac{\|A\|}{1 - \|A\|}.$$

??? proof "证明"
    由 $\|A\| < 1$ 和次可乘性，$\|A^k\| \leq \|A\|^k$，因此 $\sum_{k=0}^{\infty}\|A^k\| \leq \sum_{k=0}^{\infty}\|A\|^k = \frac{1}{1-\|A\|}$，Neumann 级数绝对收敛。于是

    $$\|(I-A)^{-1}\| = \left\|\sum_{k=0}^{\infty}A^k\right\| \leq \sum_{k=0}^{\infty}\|A^k\| \leq \frac{1}{1-\|A\|}.$$

    对于第二个不等式：

    $$\|(I-A)^{-1} - I\| = \left\|\sum_{k=1}^{\infty}A^k\right\| \leq \sum_{k=1}^{\infty}\|A\|^k = \frac{\|A\|}{1-\|A\|}. \quad \blacksquare$$

!!! example "例 14.4"
    利用 Neumann 级数近似求解线性方程组。设 $A = I - E$，其中 $E$ 是"小"扰动矩阵。若 $\|E\| < 1$，则

    $$A^{-1} = (I - E)^{-1} = I + E + E^2 + \cdots$$

    取前几项作为近似：$A^{-1} \approx I + E + E^2$。

    例如 $E = \begin{pmatrix} 0.1 & 0.05 \\ 0.02 & 0.1 \end{pmatrix}$，$\|E\|_\infty = 0.15 < 1$。

    则 $A^{-1} \approx I + E + E^2 = \begin{pmatrix} 1.111 & 0.06 \\ 0.024 & 1.111 \end{pmatrix}$（保留三位小数）。

---

## 14.3 谱半径

谱半径是矩阵分析中最核心的概念之一，它刻画了矩阵特征值的"大小"，决定了矩阵幂的渐近行为。

!!! definition "定义 14.6 (谱半径)"
    设 $A \in \mathbb{C}^{n \times n}$ 的特征值为 $\lambda_1, \lambda_2, \ldots, \lambda_n$（含重数）。$A$ 的**谱半径（spectral radius）** 定义为

    $$\rho(A) = \max_{1 \leq i \leq n} |\lambda_i|.$$

!!! theorem "定理 14.6 (谱半径与范数的关系)"
    设 $\|\cdot\|$ 是 $\mathbb{C}^{n \times n}$ 上的任意次可乘范数，则

    $$\rho(A) \leq \|A\|.$$

??? proof "证明"
    设 $\lambda$ 是 $A$ 的特征值，$\mathbf{v}$ 是对应的单位特征向量（在某个向量范数下 $\|\mathbf{v}\| = 1$）。考虑秩一矩阵 $B = \mathbf{v}\mathbf{w}^*$，其中 $\mathbf{w}$ 选取使得 $\mathbf{w}^*\mathbf{v} = 1$。

    另一种更直接的证明：设 $\lambda$ 是模最大的特征值，$A\mathbf{v} = \lambda\mathbf{v}$。构造矩阵 $X = \mathbf{v}\mathbf{e}_1^T$（其中 $\mathbf{e}_1$ 是标准基向量），则 $AX = \lambda\mathbf{v}\mathbf{e}_1^T = \lambda X$。由次可乘性 $\|AX\| \leq \|A\|\|X\|$，因此 $|\lambda|\|X\| \leq \|A\|\|X\|$。由于 $X \neq O$，得 $|\lambda| \leq \|A\|$。对所有特征值取最大值即得 $\rho(A) \leq \|A\|$。$\blacksquare$

!!! theorem "定理 14.7 (Gelfand 公式)"
    设 $A \in \mathbb{C}^{n \times n}$，$\|\cdot\|$ 是任意次可乘范数，则

    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}.$$

    此极限与范数的选取无关。

??? proof "证明"
    **下界**：由 $\rho(A^k) = \rho(A)^k$（因为 $A$ 的特征值 $\lambda_i$ 对应 $A^k$ 的特征值 $\lambda_i^k$）和定理 14.6，有

    $$\rho(A)^k = \rho(A^k) \leq \|A^k\|,$$

    因此 $\rho(A) \leq \|A^k\|^{1/k}$。故 $\rho(A) \leq \liminf_{k\to\infty}\|A^k\|^{1/k}$。

    **上界**：设 $A$ 的 Jordan 标准形为 $A = PJP^{-1}$，其中 $J = \operatorname{diag}(J_1, J_2, \ldots, J_s)$。对任意 $\varepsilon > 0$，令 $D_\varepsilon = \operatorname{diag}(1, \varepsilon, \varepsilon^2, \ldots, \varepsilon^{n-1})$，则 $D_\varepsilon^{-1} J D_\varepsilon$ 的超对角线元素被乘以 $\varepsilon$。取 $\varepsilon$ 充分小，可以使 $\|D_\varepsilon^{-1}JD_\varepsilon\|$ 在某个范数下任意接近 $\rho(A)$。

    更精确地，设 $B_\varepsilon = (PD_\varepsilon)^{-1}A(PD_\varepsilon)$，则 $\|B_\varepsilon\| \leq \rho(A) + \varepsilon$（在适当范数下）。因此

    $$\|A^k\|^{1/k} = \|(PD_\varepsilon)B_\varepsilon^k(PD_\varepsilon)^{-1}\|^{1/k} \leq \left(\|PD_\varepsilon\|\cdot\|(PD_\varepsilon)^{-1}\|\right)^{1/k}(\rho(A)+\varepsilon).$$

    当 $k \to \infty$ 时，$\left(\|PD_\varepsilon\|\cdot\|(PD_\varepsilon)^{-1}\|\right)^{1/k} \to 1$。故 $\limsup_{k\to\infty}\|A^k\|^{1/k} \leq \rho(A) + \varepsilon$。由 $\varepsilon$ 的任意性即得 $\limsup_{k\to\infty}\|A^k\|^{1/k} \leq \rho(A)$。

    综合上下界，极限存在且等于 $\rho(A)$。$\blacksquare$

!!! proposition "命题 14.1 (谱半径的基本性质)"
    设 $A \in \mathbb{C}^{n \times n}$，则：

    (1) $\rho(A^T) = \rho(A)$，$\rho(\bar{A}) = \rho(A)$，$\rho(A^*) = \rho(A)$；

    (2) $\rho(\alpha A) = |\alpha|\rho(A)$，对任意 $\alpha \in \mathbb{C}$；

    (3) $\rho(A^k) = \rho(A)^k$，对任意正整数 $k$；

    (4) 若 $A$ 是正规矩阵（$A^*A = AA^*$），则 $\rho(A) = \|A\|_2$（算子 2-范数）；

    (5) 对任意 $\varepsilon > 0$，存在次可乘范数 $\|\cdot\|$ 使得 $\|A\| \leq \rho(A) + \varepsilon$。

!!! example "例 14.5"
    计算矩阵 $A = \begin{pmatrix} 2 & 1 \\ 0 & -3 \end{pmatrix}$ 的谱半径。

    特征值为 $\lambda_1 = 2$，$\lambda_2 = -3$，因此 $\rho(A) = \max\{|2|, |-3|\} = 3$。

    用 Gelfand 公式验证：$A^2 = \begin{pmatrix} 4 & -1 \\ 0 & 9 \end{pmatrix}$，$\|A^2\|_\infty = \max\{5, 9\} = 9$，$\|A^2\|_\infty^{1/2} = 3$。

    $A^4 = \begin{pmatrix} 16 & -10 \\ 0 & 81 \end{pmatrix}$，$\|A^4\|_\infty^{1/4} = 81^{1/4} = 3$。序列趋于 $\rho(A) = 3$。

---

## 14.4 Gershgorin 圆盘定理

Gershgorin 圆盘定理是特征值定位理论中最优雅、最实用的结果之一。它仅利用矩阵元素就给出了特征值的位置估计。

!!! definition "定义 14.7 (Gershgorin 圆盘)"
    设 $A = (a_{ij}) \in \mathbb{C}^{n \times n}$。第 $i$ 个 **Gershgorin 圆盘（Gershgorin disc）** 定义为

    $$D_i = \left\{z \in \mathbb{C} : |z - a_{ii}| \leq R_i\right\}, \quad R_i = \sum_{j \neq i} |a_{ij}|,$$

    其中 $R_i$ 称为第 $i$ 行的**删去行和（deleted row sum）**。

!!! theorem "定理 14.8 (Gershgorin 圆盘定理)"
    设 $A = (a_{ij}) \in \mathbb{C}^{n \times n}$，则 $A$ 的每个特征值至少属于一个 Gershgorin 圆盘，即

    $$\sigma(A) \subseteq \bigcup_{i=1}^{n} D_i.$$

??? proof "证明"
    设 $\lambda$ 是 $A$ 的特征值，$\mathbf{x} = (x_1, x_2, \ldots, x_n)^T$ 是对应的特征向量。选取 $p$ 使得 $|x_p| = \max_{1 \leq i \leq n}|x_i|$。由 $\mathbf{x} \neq \mathbf{0}$ 知 $|x_p| > 0$。

    由 $A\mathbf{x} = \lambda\mathbf{x}$ 的第 $p$ 个分量：

    $$\sum_{j=1}^{n} a_{pj}x_j = \lambda x_p.$$

    因此

    $$(\lambda - a_{pp})x_p = \sum_{j \neq p} a_{pj}x_j.$$

    取模并利用三角不等式：

    $$|\lambda - a_{pp}| \cdot |x_p| \leq \sum_{j \neq p} |a_{pj}|\cdot|x_j| \leq |x_p| \sum_{j \neq p} |a_{pj}| = |x_p| R_p.$$

    除以 $|x_p| > 0$ 得 $|\lambda - a_{pp}| \leq R_p$，即 $\lambda \in D_p$。$\blacksquare$

!!! theorem "定理 14.9 (Gershgorin 圆盘的连通分量)"
    若 $n$ 个 Gershgorin 圆盘的并集可以分成 $k$ 个互不相交的连通区域，每个连通区域分别由 $m_1, m_2, \ldots, m_k$ 个圆盘的并构成（$m_1 + m_2 + \cdots + m_k = n$），则第 $j$ 个连通区域恰好包含 $A$ 的 $m_j$ 个特征值（计重数）。

??? proof "证明"
    考虑矩阵族 $A(t) = D + t(A - D)$，其中 $D = \operatorname{diag}(a_{11}, \ldots, a_{nn})$，$t \in [0, 1]$。当 $t = 0$ 时 $A(0) = D$，特征值为 $a_{11}, \ldots, a_{nn}$；当 $t = 1$ 时 $A(1) = A$。

    $A(t)$ 的第 $i$ 个 Gershgorin 圆盘的半径为 $tR_i$。当 $t$ 从 $0$ 连续增大到 $1$ 时，圆盘连续膨胀。特征值是特征多项式系数的连续函数，因此特征值关于 $t$ 连续变化。

    在 $t = 0$ 时，$m_j$ 个圆盘（退化为点）对应 $m_j$ 个特征值。由于特征值连续变化，且不同连通区域之间无交集，特征值不能从一个连通区域跳到另一个。因此每个连通区域恰含 $m_j$ 个特征值。$\blacksquare$

!!! corollary "推论 14.1 (严格对角占优矩阵可逆)"
    若 $A$ 是**严格对角占优矩阵（strictly diagonally dominant matrix）**，即对每个 $i$ 都有 $|a_{ii}| > R_i = \sum_{j \neq i}|a_{ij}|$，则 $A$ 可逆。

??? proof "证明"
    若 $A$ 不可逆，则 $0$ 是 $A$ 的特征值。由 Gershgorin 定理，$0 \in D_i$ 对某个 $i$ 成立，即 $|a_{ii}| \leq R_i$，矛盾于严格对角占优条件。$\blacksquare$

!!! example "例 14.6"
    设 $A = \begin{pmatrix} 4 & -1 & 0 \\ 1 & 5 & -1 \\ 0 & -1 & 3 \end{pmatrix}$。

    三个 Gershgorin 圆盘为：

    - $D_1$：圆心 $4$，半径 $1$，即 $\{z : |z - 4| \leq 1\} = [3, 5]$；
    - $D_2$：圆心 $5$，半径 $2$，即 $\{z : |z - 5| \leq 2\} = [3, 7]$；
    - $D_3$：圆心 $3$，半径 $1$，即 $\{z : |z - 3| \leq 1\} = [2, 4]$。

    三个圆盘的并集为 $[2, 7]$（一个连通区域），因此三个特征值均在 $[2, 7]$ 中。

    实际计算特征值约为 $\lambda \approx 2.38, 4.18, 5.44$，确实都在 $[2, 7]$ 中。

!!! example "例 14.7"
    设 $A = \begin{pmatrix} 10 & 0.1 & 0.2 \\ 0.1 & 20 & 0.3 \\ 0.2 & 0.1 & 30 \end{pmatrix}$。

    三个 Gershgorin 圆盘为：

    - $D_1$：$|z - 10| \leq 0.3$，即 $[9.7, 10.3]$；
    - $D_2$：$|z - 20| \leq 0.4$，即 $[19.6, 20.4]$；
    - $D_3$：$|z - 30| \leq 0.3$，即 $[29.7, 30.3]$。

    三个圆盘互不相交，因此由定理 14.9，每个圆盘恰含一个特征值。特别地，$A$ 有三个彼此相距较远的实特征值，分别在 $10, 20, 30$ 附近。

---

## 14.5 矩阵微积分

将微积分的概念推广到矩阵值函数和关于矩阵变量的函数，是矩阵分析的重要内容，在优化、控制理论和统计学中有广泛应用。

!!! definition "定义 14.8 (矩阵值函数的导数)"
    设 $A(t) = (a_{ij}(t))$ 是定义在区间 $I \subseteq \mathbb{R}$ 上的矩阵值函数，其中每个 $a_{ij}(t)$ 是 $t$ 的实值函数。$A(t)$ 关于 $t$ 的**导数**定义为

    $$\frac{dA}{dt} = \left(\frac{da_{ij}}{dt}\right),$$

    即逐元素求导。

!!! theorem "定理 14.10 (矩阵值函数的求导法则)"
    设 $A(t)$，$B(t)$ 是可微的矩阵值函数，$c(t)$ 是可微标量函数，则：

    (1) $\frac{d}{dt}(A + B) = \frac{dA}{dt} + \frac{dB}{dt}$；

    (2) $\frac{d}{dt}(cA) = \frac{dc}{dt}A + c\frac{dA}{dt}$；

    (3) $\frac{d}{dt}(AB) = \frac{dA}{dt}B + A\frac{dB}{dt}$（注意顺序）；

    (4) 若 $A(t)$ 可逆，则 $\frac{d}{dt}A^{-1} = -A^{-1}\frac{dA}{dt}A^{-1}$。

??? proof "证明"
    (1)—(3) 由逐元素求导直接验证。

    **(4)** 由 $A(t)A^{-1}(t) = I$，两端对 $t$ 求导：

    $$\frac{dA}{dt}A^{-1} + A\frac{dA^{-1}}{dt} = O.$$

    解出 $\frac{dA^{-1}}{dt} = -A^{-1}\frac{dA}{dt}A^{-1}$。$\blacksquare$

!!! definition "定义 14.9 (标量函数关于矩阵的导数)"
    设 $f : \mathbb{R}^{m \times n} \to \mathbb{R}$ 是标量值函数。$f$ 关于矩阵 $X = (x_{ij})$ 的**导数（derivative）**（或**梯度**）定义为

    $$\frac{\partial f}{\partial X} = \left(\frac{\partial f}{\partial x_{ij}}\right) \in \mathbb{R}^{m \times n}.$$

!!! proposition "命题 14.2 (常用矩阵导数公式)"
    设 $A$ 为常数矩阵，$X$ 为矩阵变量，$\mathbf{x}$ 为向量变量，则：

    (1) $\frac{\partial}{\partial \mathbf{x}}(\mathbf{a}^T\mathbf{x}) = \mathbf{a}$；

    (2) $\frac{\partial}{\partial \mathbf{x}}(\mathbf{x}^T A\mathbf{x}) = (A + A^T)\mathbf{x}$；若 $A$ 对称则等于 $2A\mathbf{x}$；

    (3) $\frac{\partial}{\partial X}\operatorname{tr}(AX) = A^T$；

    (4) $\frac{\partial}{\partial X}\operatorname{tr}(X^TAX) = (A + A^T)X$；

    (5) $\frac{\partial}{\partial X}\ln\det(X) = X^{-T}$（当 $X$ 可逆时）。

!!! example "例 14.8"
    设 $f(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} + \mathbf{b}^T\mathbf{x} + c$，其中 $A$ 是 $n \times n$ 对称矩阵，$\mathbf{b} \in \mathbb{R}^n$，$c \in \mathbb{R}$。

    则 $\frac{\partial f}{\partial \mathbf{x}} = 2A\mathbf{x} + \mathbf{b}$。

    令 $\frac{\partial f}{\partial \mathbf{x}} = \mathbf{0}$，得驻点 $\mathbf{x}^* = -\frac{1}{2}A^{-1}\mathbf{b}$（当 $A$ 可逆时）。

    二阶导数（Hessian 矩阵）为 $\frac{\partial^2 f}{\partial \mathbf{x}\partial \mathbf{x}^T} = 2A$。若 $A$ 正定则 $\mathbf{x}^*$ 是严格极小值点。

!!! example "例 14.9"
    **矩阵指数的导数**。设 $A$ 为常数矩阵，则矩阵指数 $e^{tA} = \sum_{k=0}^{\infty}\frac{(tA)^k}{k!}$ 满足

    $$\frac{d}{dt}e^{tA} = Ae^{tA} = e^{tA}A.$$

    这可以通过逐项求导得到：

    $$\frac{d}{dt}\sum_{k=0}^{\infty}\frac{t^k A^k}{k!} = \sum_{k=1}^{\infty}\frac{kt^{k-1}A^k}{k!} = A\sum_{k=1}^{\infty}\frac{t^{k-1}A^{k-1}}{(k-1)!} = Ae^{tA}.$$

    这是线性常微分方程组 $\frac{d\mathbf{x}}{dt} = A\mathbf{x}$ 的基解矩阵。

---

## 14.6 矩阵的特征值与奇异值的关系

矩阵的特征值和奇异值是两组不同但密切相关的谱信息。奇异值由 $A^*A$ 的特征值决定，而特征值直接来自 $A$ 本身。本节探讨它们之间的深层联系。

!!! definition "定义 14.10 (奇异值)"
    设 $A \in \mathbb{C}^{m \times n}$，$A$ 的**奇异值（singular values）** $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_{\min(m,n)} \geq 0$ 定义为 $A^*A$ 的特征值的非负平方根，即 $\sigma_i = \sqrt{\lambda_i(A^*A)}$。

!!! theorem "定理 14.11 (特征值与奇异值的模不等式)"
    设 $A \in \mathbb{C}^{n \times n}$ 的特征值按模降序排列为 $|\lambda_1| \geq |\lambda_2| \geq \cdots \geq |\lambda_n|$，奇异值降序排列为 $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_n$。则对每个 $k = 1, 2, \ldots, n$，

    $$\prod_{i=1}^{k} |\lambda_i| \leq \prod_{i=1}^{k} \sigma_i.$$

    特别地，$|\lambda_1| \leq \sigma_1$，且 $\prod_{i=1}^{n}|\lambda_i| = \prod_{i=1}^{n}\sigma_i = |\det A|$。

??? proof "证明"
    由 Schur 分解 $A = UTU^*$（$U$ 酉矩阵，$T$ 上三角），有 $A^*A = UT^*T U^*$。

    对 $k = 1$ 的情形：$|\lambda_1| \leq \|A\|_2 = \sigma_1$，因为 $\|A\|_2$ 是算子 2-范数。

    对 $k = n$ 的情形：$\prod_{i=1}^{n}|\lambda_i| = |\det A| = \prod_{i=1}^{n}\sigma_i$，因为 $|\det A|^2 = \det(A^*A) = \prod \sigma_i^2$。

    一般情形需要用到复合矩阵（compound matrix）的理论。设 $C_k(A)$ 为 $A$ 的第 $k$ 阶复合矩阵，其特征值为 $A$ 的所有 $k$ 阶特征值之积 $\lambda_{i_1}\cdots\lambda_{i_k}$。则

    $$\prod_{i=1}^{k}|\lambda_i| \leq \|C_k(A)\|_2 = \sigma_1(C_k(A)) = \prod_{i=1}^{k}\sigma_i(A). \quad \blacksquare$$

!!! theorem "定理 14.12 (Weyl 不等式——奇异值版本)"
    设 $A, B \in \mathbb{C}^{m \times n}$，奇异值分别降序排列为 $\sigma_1(A) \geq \cdots$ 和 $\sigma_1(B) \geq \cdots$，则对所有 $i$，

    $$|\sigma_i(A) - \sigma_i(B)| \leq \|A - B\|_2.$$

    特别地，$|\sigma_1(A) - \sigma_1(B)| \leq \|A - B\|_2$。

??? proof "证明"
    将在第 15 章中给出完整证明。此处给出思路：利用奇异值的极大极小刻画

    $$\sigma_i(A) = \min_{\dim V = n-i+1}\max_{\substack{\mathbf{x} \in V \\ \|\mathbf{x}\|=1}}\|A\mathbf{x}\|,$$

    以及三角不等式 $\|A\mathbf{x}\| \leq \|B\mathbf{x}\| + \|(A-B)\mathbf{x}\| \leq \|B\mathbf{x}\| + \|A-B\|_2$，可以推出结论。$\blacksquare$

!!! example "例 14.10"
    设 $A = \begin{pmatrix} 3 & 1 \\ 0 & 2 \end{pmatrix}$。

    **特征值**：$\lambda_1 = 3$，$\lambda_2 = 2$；$|\lambda_1| = 3$，$|\lambda_2| = 2$。

    **奇异值**：$A^*A = \begin{pmatrix} 9 & 3 \\ 3 & 5 \end{pmatrix}$，特征值为 $\frac{14 \pm \sqrt{36}}{2} = \frac{14 \pm 6}{2}$，即 $10$ 和 $4$。奇异值 $\sigma_1 = \sqrt{10} \approx 3.16$，$\sigma_2 = 2$。

    验证：$|\lambda_1| = 3 \leq 3.16 = \sigma_1$。$|\lambda_1||\lambda_2| = 6 = \sigma_1\sigma_2 = \sqrt{10}\cdot 2 = 2\sqrt{10} \approx 6.32$？

    等一下，让我重新验证。$|\det A| = |3 \cdot 2 - 1 \cdot 0| = 6$。$\sigma_1\sigma_2 = \sqrt{\det(A^*A)} = \sqrt{45 - 9} = \sqrt{36} = 6$。所以 $\prod|\lambda_i| = 6 = \prod\sigma_i = 6$。

    重新计算：$A^*A = \begin{pmatrix} 9 & 3 \\ 3 & 5 \end{pmatrix}$，$\det(A^*A) = 45 - 9 = 36$，$\operatorname{tr}(A^*A) = 14$。特征值满足 $\mu^2 - 14\mu + 36 = 0$，$\mu = 7 \pm \sqrt{13}$。所以 $\sigma_1 = \sqrt{7+\sqrt{13}} \approx 3.21$，$\sigma_2 = \sqrt{7-\sqrt{13}} \approx 1.87$。$\sigma_1\sigma_2 = \sqrt{36} = 6 = |\lambda_1\lambda_2|$。

    $|\lambda_1| = 3 \leq 3.21 = \sigma_1$。$|\lambda_1||\lambda_2| = 6 = \sigma_1\sigma_2$。不等式成立。

---

## 14.7 矩阵的极限与收敛

矩阵幂 $A^k$ 的收敛行为完全由谱半径决定，这一结果在迭代法、Markov 链和动力系统理论中有基本的重要性。

!!! theorem "定理 14.13 (矩阵幂收敛的充要条件)"
    设 $A \in \mathbb{C}^{n \times n}$，则

    $$\lim_{k \to \infty} A^k = O \quad \Longleftrightarrow \quad \rho(A) < 1.$$

??? proof "证明"
    **($\Leftarrow$)**：设 $\rho(A) < 1$。由 Gelfand 公式，对任意满足 $\rho(A) < r < 1$ 的 $r$，存在正整数 $N$ 使得当 $k \geq N$ 时 $\|A^k\|^{1/k} < r$，即 $\|A^k\| < r^k$。因此 $\|A^k\| \to 0$，即 $A^k \to O$。

    **($\Rightarrow$)**：设 $A^k \to O$。若 $\rho(A) \geq 1$，则存在特征值 $\lambda$ 满足 $|\lambda| \geq 1$。设 $A\mathbf{v} = \lambda\mathbf{v}$（$\mathbf{v} \neq \mathbf{0}$），则 $A^k\mathbf{v} = \lambda^k\mathbf{v}$，$\|A^k\mathbf{v}\| = |\lambda|^k\|\mathbf{v}\| \geq \|\mathbf{v}\| > 0$。但 $A^k \to O$ 意味着 $A^k\mathbf{v} \to \mathbf{0}$，矛盾。$\blacksquare$

!!! theorem "定理 14.14 (矩阵幂的有界性)"
    设 $A \in \mathbb{C}^{n \times n}$，则序列 $\{A^k\}_{k=0}^{\infty}$ 有界（即 $\sup_k \|A^k\| < \infty$）当且仅当 $\rho(A) \leq 1$，且模为 $1$ 的特征值都是半单的（即其代数重数等于几何重数，或等价地，在 Jordan 标准形中对应 $1 \times 1$ 的 Jordan 块）。

??? proof "证明"
    设 $A = PJP^{-1}$，其中 $J$ 是 Jordan 标准形。$A^k = PJ^kP^{-1}$。

    **充分性**：若 $\rho(A) \leq 1$ 且模为 $1$ 的特征值半单，则每个 Jordan 块 $J_i$ 要么对应 $|\lambda_i| < 1$（此时 $J_i^k \to O$），要么是 $1 \times 1$ 块 $(\lambda_i)$，其中 $|\lambda_i| = 1$（此时 $|J_i^k| = |\lambda_i|^k = 1$）。无论哪种情况 $J_i^k$ 有界，因此 $J^k$ 有界，从而 $A^k = PJ^kP^{-1}$ 有界。

    **必要性**：若 $\rho(A) > 1$，由 Gelfand 公式 $\|A^k\|^{1/k} \to \rho(A) > 1$，$\|A^k\| \to \infty$。若存在模为 $1$ 的特征值 $\lambda$ 对应 $s \times s$（$s \geq 2$）的 Jordan 块，则 $J_i^k$ 的 $(1,2)$ 元素为 $k\lambda^{k-1}$，模为 $k \to \infty$，不有界。$\blacksquare$

!!! proposition "命题 14.3 (矩阵幂收敛的速率)"
    若 $\rho(A) < 1$，则对任意满足 $\rho(A) < r < 1$ 的 $r$，存在常数 $C > 0$ 使得

    $$\|A^k\| \leq Cr^k, \quad \forall k \geq 0.$$

    即 $A^k$ 以指数速率趋于零，速率由 $\rho(A)$ 决定。

!!! example "例 14.11"
    判断以下矩阵的幂是否收敛于零矩阵。

    (a) $A = \begin{pmatrix} 0.5 & 0.3 \\ 0.1 & 0.4 \end{pmatrix}$。

    特征多项式：$\lambda^2 - 0.9\lambda + 0.17 = 0$，$\lambda = \frac{0.9 \pm \sqrt{0.81 - 0.68}}{2} = \frac{0.9 \pm \sqrt{0.13}}{2}$。

    $\sqrt{0.13} \approx 0.361$，故 $\lambda_1 \approx 0.63$，$\lambda_2 \approx 0.27$。$\rho(A) \approx 0.63 < 1$，因此 $A^k \to O$。

    (b) $B = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$。

    特征值为 $\pm i$，$\rho(B) = 1$。$B$ 是正规矩阵，特征值半单，因此 $\{B^k\}$ 有界但不趋于零。实际上 $B^2 = -I$，$B^4 = I$，序列以周期 $4$ 振荡。

!!! example "例 14.12"
    **Jacobi 迭代法的收敛性**。设线性方程组 $A\mathbf{x} = \mathbf{b}$，分裂 $A = D - L - U$（$D$ 为对角部分，$-L$ 为严格下三角，$-U$ 为严格上三角）。Jacobi 迭代矩阵为 $T_J = D^{-1}(L + U)$。

    Jacobi 迭代 $\mathbf{x}^{(k+1)} = T_J\mathbf{x}^{(k)} + D^{-1}\mathbf{b}$ 收敛的充要条件是 $\rho(T_J) < 1$。

    例如对矩阵 $A = \begin{pmatrix} 4 & 1 \\ 1 & 3 \end{pmatrix}$，$T_J = \begin{pmatrix} 0 & -1/4 \\ -1/3 & 0 \end{pmatrix}$，特征值为 $\pm\frac{1}{2\sqrt{3}}$，$\rho(T_J) = \frac{1}{2\sqrt{3}} \approx 0.289 < 1$，因此 Jacobi 迭代收敛。
