# 第 17 章 非负矩阵与 Perron-Frobenius 理论

非负矩阵理论是线性代数中一个美丽而深刻的分支，其核心是 Perron-Frobenius 定理——它揭示了非负矩阵的谱结构具有令人惊叹的规律性。该理论最初源于 Perron（1907）对正矩阵的研究和 Frobenius（1912）对不可约非负矩阵的推广，后来成为概率论（Markov 链）、经济学（投入产出模型）、网络分析（PageRank 算法）和人口动力学等领域的数学基础。

---

## 17.1 非负矩阵的定义

!!! definition "定义 17.1 (非负矩阵与正矩阵)"
    设 $A = (a_{ij}) \in \mathbb{R}^{m \times n}$。

    - 若 $a_{ij} \geq 0$ 对所有 $i, j$ 成立，则称 $A$ 为**非负矩阵（nonnegative matrix）**，记作 $A \geq O$。
    - 若 $a_{ij} > 0$ 对所有 $i, j$ 成立，则称 $A$ 为**正矩阵（positive matrix）**，记作 $A > O$。
    - 类似地定义非负向量 $\mathbf{x} \geq \mathbf{0}$ 和正向量 $\mathbf{x} > \mathbf{0}$。

!!! definition "定义 17.2 (矩阵的元素序)"
    对同型矩阵 $A, B$，定义 $A \geq B$ 表示 $A - B \geq O$（即 $a_{ij} \geq b_{ij}$ 对所有 $i, j$）。

!!! note "注"
    矩阵的元素序 $A \geq B$（逐元素比较）与 Löwner 偏序 $A \succeq B$（$A - B$ 半正定）是完全不同的概念，切勿混淆。

!!! proposition "命题 17.1 (非负矩阵的基本性质)"
    (1) 若 $A \geq O$，$B \geq O$，则 $A + B \geq O$；

    (2) 若 $A \geq O$，$\alpha \geq 0$，则 $\alpha A \geq O$；

    (3) 若 $A \geq O$，$B \geq O$，尺寸相容，则 $AB \geq O$；

    (4) 若 $A \geq O$，则 $A^k \geq O$ 对所有正整数 $k$。

!!! example "例 17.1"
    在实际应用中，非负矩阵大量出现：

    - **转移概率矩阵**：Markov 链中的转移矩阵 $P = (p_{ij})$，$p_{ij} \geq 0$，$\sum_j p_{ij} = 1$；
    - **邻接矩阵**：有向图的邻接矩阵 $A$，$a_{ij} \in \{0, 1\}$；
    - **投入产出矩阵**：Leontief 经济模型中的技术系数矩阵。

---

## 17.2 Perron 定理（正矩阵）

!!! theorem "定理 17.1 (Perron 定理)"
    设 $A > O$ 是 $n \times n$ 正矩阵（$n \geq 2$）。则：

    (1) $\rho(A) > 0$ 且 $\rho(A)$ 是 $A$ 的特征值（称为 **Perron 根（Perron root）**）；

    (2) $\rho(A)$ 是 $A$ 的代数单特征值（重数为 $1$）；

    (3) 存在正向量 $\mathbf{v} > \mathbf{0}$ 使得 $A\mathbf{v} = \rho(A)\mathbf{v}$（**Perron 向量**）；

    (4) $A$ 的任何其他特征值 $\lambda$ 满足 $|\lambda| < \rho(A)$（严格不等号）；

    (5) Perron 向量（在归一化后）是唯一的非负特征向量。

??? proof "证明"
    **存在性**。考虑紧集 $\Delta = \{\mathbf{x} \in \mathbb{R}^n : \mathbf{x} \geq \mathbf{0},\; \sum x_i = 1\}$（单纯形）。定义

    $$r(\mathbf{x}) = \min_{i: x_i > 0} \frac{(A\mathbf{x})_i}{x_i}.$$

    由于 $A > O$，对 $\mathbf{x} \in \Delta$（$\mathbf{x} \neq \mathbf{0}$），$A\mathbf{x} > \mathbf{0}$。在 $\Delta$ 的内部 $\Delta^\circ = \{\mathbf{x} > \mathbf{0}, \sum x_i = 1\}$ 上，$r(\mathbf{x})$ 定义良好且连续。

    取 $\rho^* = \sup_{\mathbf{x} \in \Delta} r(\mathbf{x})$。可以证明上确界在某个 $\mathbf{v}^* \in \Delta^\circ$ 处达到，且 $A\mathbf{v}^* = \rho^*\mathbf{v}^*$。

    实际上，$\rho^* = \rho(A)$。若存在特征值 $\lambda$ 满足 $|\lambda| > \rho^*$，设 $A\mathbf{w} = \lambda\mathbf{w}$（$\mathbf{w}$ 可能为复数），则 $A|\mathbf{w}| \geq |A\mathbf{w}| = |\lambda||\mathbf{w}|$（分量意义），由此可推出 $r(|\mathbf{w}|/\||{\mathbf{w}}\||_1) \geq |\lambda| > \rho^*$，矛盾。

    **唯一性（代数单根）**。假设 $\mathbf{w}$ 是对应于 $\rho(A)$ 的另一个线性无关特征向量。取实部和虚部，可以找到实特征向量。但实特征向量中必存在分量为零或负的，设 $A\mathbf{u} = \rho(A)\mathbf{u}$，$\mathbf{u}$ 有非正分量。由 $A > O$，$\rho(A)\mathbf{u} = A\mathbf{u}$，取绝对值 $\rho(A)|\mathbf{u}| \leq A|\mathbf{u}|$，且在某些分量上严格不等（因为 $\mathbf{u}$ 有正有负分量混合）。这与 $\rho(A)$ 的最大性矛盾。

    **严格支配性**。设 $\lambda$ 是另一特征值，$|\lambda| = \rho(A)$。设 $A\mathbf{w} = \lambda\mathbf{w}$。由三角不等式 $\rho(A)|\mathbf{w}| = |\lambda\mathbf{w}| = |A\mathbf{w}| \leq A|\mathbf{w}|$。若某处严格不等式成立则 $r(|\mathbf{w}|) > \rho(A)$，矛盾。若处处等号成立，需要 $A\mathbf{w}$ 各分量符号一致，这只在 $\mathbf{w}$ 为正向量或负向量时才可能，从而 $\lambda = \rho(A)$（正实数）。$\blacksquare$

!!! example "例 17.2"
    设 $A = \begin{pmatrix} 2 & 1 \\ 3 & 4 \end{pmatrix} > O$。

    特征多项式：$\lambda^2 - 6\lambda + 5 = 0$，$\lambda_1 = 5$，$\lambda_2 = 1$。

    $\rho(A) = 5$，Perron 根。$\lambda_2 = 1 < 5$。

    Perron 向量：$A\mathbf{v} = 5\mathbf{v}$，$(A - 5I)\mathbf{v} = \mathbf{0}$，$\begin{pmatrix} -3 & 1 \\ 3 & -1 \end{pmatrix}\mathbf{v} = \mathbf{0}$，$\mathbf{v} = t\begin{pmatrix} 1 \\ 3 \end{pmatrix}$。

    归一化 Perron 向量 $\mathbf{v} = \begin{pmatrix} 1/4 \\ 3/4 \end{pmatrix}$（使 $\sum v_i = 1$），确实为正向量。

!!! example "例 17.3"
    $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{pmatrix}$ 是正矩阵。

    特征值约为 $\lambda_1 \approx 16.12$，$\lambda_2 \approx -1.12$，$\lambda_3 \approx 0$。

    $\rho(A) \approx 16.12$，是代数单特征值。$|\lambda_2| \approx 1.12 < 16.12$，$|\lambda_3| \approx 0 < 16.12$。

    对应的 Perron 向量（归一化后各分量为正）约为 $\mathbf{v} \approx (0.164, 0.387, 0.610)^T > \mathbf{0}$。

---

## 17.3 不可约矩阵

!!! definition "定义 17.3 (可约与不可约矩阵)"
    设 $A \in \mathbb{R}^{n \times n}$（$n \geq 2$）。若存在置换矩阵 $P$ 使得

    $$P^TAP = \begin{pmatrix} B & C \\ O & D \end{pmatrix},$$

    其中 $B, D$ 是方阵且 $B$ 至少 $1 \times 1$，$D$ 至少 $1 \times 1$，则称 $A$ **可约（reducible）**。否则称 $A$ **不可约（irreducible）**。

    $1 \times 1$ 的非零矩阵定义为不可约。

!!! definition "定义 17.4 (有向图的解释)"
    $n \times n$ 非负矩阵 $A$ 的**关联有向图（associated directed graph）** $G(A)$ 以 $\{1, 2, \ldots, n\}$ 为顶点集，若 $a_{ij} > 0$ 则从 $i$ 到 $j$ 有一条有向边。

    $A$ 不可约当且仅当 $G(A)$ 是**强连通的（strongly connected）**，即对任意两个顶点 $i, j$，存在从 $i$ 到 $j$ 的有向路径。

!!! theorem "定理 17.2 (不可约的等价条件)"
    设 $A \geq O$ 是 $n \times n$ 矩阵（$n \geq 2$）。以下条件等价：

    (1) $A$ 不可约；

    (2) $G(A)$ 强连通；

    (3) $(I + A)^{n-1} > O$（所有元素为正）；

    (4) 不存在指标集 $S \subsetneq \{1, \ldots, n\}$（$S \neq \emptyset$）使得 $a_{ij} = 0$ 对所有 $i \in S, j \notin S$ 成立。

??? proof "证明"
    **(1) $\Leftrightarrow$ (2)**：标准结果，可约矩阵对应的有向图不强连通。

    **(2) $\Leftrightarrow$ (3)**：$(I+A)^{n-1}$ 的 $(i,j)$ 元素为 $\sum_{k=0}^{n-1}\binom{n-1}{k}(A^k)_{ij}$ 的某种组合。$(I+A)^{n-1}_{ij} > 0$ 当且仅当存在从 $i$ 到 $j$ 的长度不超过 $n-1$ 的路径（包括长度为 $0$ 的自环），这等价于强连通。$\blacksquare$

!!! example "例 17.4"
    $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$ 是不可约的。

    关联有向图：$1 \to 2 \to 3 \to 1$，构成一个环，强连通。

    $B = \begin{pmatrix} 1 & 1 \\ 0 & 2 \end{pmatrix}$ 是可约的（已是上三角形式，对应图 $1 \to 2$ 但无 $2 \to 1$）。

---

## 17.4 Perron-Frobenius 定理（不可约非负矩阵）

!!! theorem "定理 17.3 (Perron-Frobenius 定理)"
    设 $A \geq O$ 是 $n \times n$ 不可约非负矩阵。则：

    (1) $\rho(A) > 0$ 且 $\rho(A)$ 是 $A$ 的特征值（Perron 根）；

    (2) 对应于 $\rho(A)$ 的特征向量可取为正向量 $\mathbf{v} > \mathbf{0}$（Perron 向量）；

    (3) $\rho(A)$ 作为特征值的代数重数为 $1$；

    (4) 若 $A$ 有 $h$ 个模等于 $\rho(A)$ 的特征值，则它们恰好是

    $$\rho(A), \rho(A)e^{2\pi i/h}, \rho(A)e^{4\pi i/h}, \ldots, \rho(A)e^{2(h-1)\pi i/h},$$

    即均匀分布在半径为 $\rho(A)$ 的圆上。$h$ 称为 $A$ 的**循环指数（index of imprimitivity）**；

    (5) 若 $\mathbf{x} \geq \mathbf{0}$，$\mathbf{x} \neq \mathbf{0}$，$A\mathbf{x} \leq \lambda\mathbf{x}$（元素意义），则 $\lambda \geq \rho(A)$；若 $A\mathbf{x} \geq \lambda\mathbf{x}$，则 $\lambda \leq \rho(A)$。

??? proof "证明"
    证明的核心思路如下。

    **存在性（Perron 根与正特征向量）**。定义 $\rho^* = \sup\{r \geq 0 : \exists \mathbf{x} \geq \mathbf{0}, \mathbf{x} \neq \mathbf{0}, A\mathbf{x} \geq r\mathbf{x}\}$。可以证明 $\rho^* = \rho(A)$ 且该上确界可达。

    设 $\mathbf{x}^*$ 是使上确界达到的非负向量。由 $A\mathbf{x}^* \geq \rho^*\mathbf{x}^*$ 和不可约性，可以证明 $\mathbf{x}^* > \mathbf{0}$（因为若某分量为零，不可约性保证经过若干次乘以 $A$ 后所有分量变为正）。进而 $A\mathbf{x}^* = \rho^*\mathbf{x}^*$（否则可以增大 $r$）。

    **代数单根**。设 $\rho = \rho(A)$。由正特征向量 $\mathbf{v} > \mathbf{0}$ 的存在，可以构造对角矩阵 $D = \operatorname{diag}(\mathbf{v})$，$B = D^{-1}AD$ 是行和为 $\rho$ 的非负矩阵。设 $\mathbf{w}$ 是 $B$ 对应于 $\rho$ 的另一特征向量，用类似于 Perron 定理的论证可证明 $\mathbf{w}$ 必须是 $\mathbf{v}$ 的标量倍，故代数重数为 $1$。

    **循环结构**。不可约矩阵的关联有向图中，所有环的长度的最大公约数 $h$ 就是循环指数。通过将顶点集按照到某固定顶点的距离模 $h$ 分类，可以将 $A$ 变换为循环分块形式，从而证明特征值具有旋转对称性。$\blacksquare$

!!! example "例 17.5"
    设 $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 2 \\ 3 & 0 & 0 \end{pmatrix}$，不可约（$1\to 2\to 3\to 1$）。

    特征多项式：$\lambda^3 - 6 = 0$，$\lambda = \sqrt[3]{6}, \sqrt[3]{6}\omega, \sqrt[3]{6}\omega^2$，其中 $\omega = e^{2\pi i/3}$。

    $\rho(A) = \sqrt[3]{6} \approx 1.817$。三个特征值的模都等于 $\sqrt[3]{6}$，均匀分布在圆上。循环指数 $h = 3$（图中唯一的环长度为 $3$，$\gcd = 3$）。

    Perron 向量：$A\mathbf{v} = \sqrt[3]{6}\mathbf{v}$，解得 $\mathbf{v} \propto (\sqrt[3]{36}, 2\sqrt[3]{6}, 6)^T / c > \mathbf{0}$。

---

## 17.5 本原矩阵

!!! definition "定义 17.5 (本原矩阵)"
    非负不可约矩阵 $A$ 称为**本原的（primitive）**，若其循环指数 $h = 1$，即 $\rho(A)$ 是唯一的模等于 $\rho(A)$ 的特征值。

    等价地，$A \geq O$ 本原当且仅当存在正整数 $k$ 使得 $A^k > O$（所有元素为正）。

!!! theorem "定理 17.4 (本原矩阵的判别)"
    设 $A \geq O$ 不可约。以下条件等价：

    (1) $A$ 本原；

    (2) 循环指数 $h = 1$；

    (3) 存在 $k$ 使得 $A^k > O$；

    (4) 关联有向图中所有环长度的最大公约数为 $1$。

!!! theorem "定理 17.5 (本原矩阵的收敛性)"
    设 $A \geq O$ 是 $n \times n$ 本原矩阵，Perron 根为 $\rho = \rho(A)$。设 $\mathbf{v}$ 和 $\mathbf{w}^T$ 分别是右 Perron 向量和左 Perron 向量，归一化使得 $\mathbf{w}^T\mathbf{v} = 1$。则

    $$\lim_{k \to \infty} \left(\frac{A}{\rho}\right)^k = \mathbf{v}\mathbf{w}^T.$$

??? proof "证明"
    由 Jordan 分解。设 $A$ 的特征值为 $\rho = \lambda_1 > |\lambda_2| \geq \cdots \geq |\lambda_n|$（本原性保证严格不等号）。$A$ 的 Jordan 分解为

    $$A = \rho\mathbf{v}\mathbf{w}^T + \sum_{i=2}^{n}\lambda_i P_i + N_i,$$

    其中 $P_i$ 是投影，$N_i$ 是幂零部分。因此

    $$\left(\frac{A}{\rho}\right)^k = \mathbf{v}\mathbf{w}^T + \sum_{i=2}^{n}\left(\frac{\lambda_i}{\rho}\right)^k(\cdots).$$

    由 $|\lambda_i/\rho| < 1$，后面的项趋于零（幂零部分贡献多项式增长，但被指数衰减淹没）。$\blacksquare$

!!! example "例 17.6"
    $A = \begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix}$ 是正矩阵，自动本原。

    特征值：$\lambda_1 = 1$，$\lambda_2 = 0.2$。$\rho(A) = 1$。

    右 Perron 向量 $\mathbf{v} = \begin{pmatrix} 3/8 \\ 5/8 \end{pmatrix}$（归一化使分量之和为 $1$）。左 Perron 向量 $\mathbf{w}^T = (1, 1)$（因为行和为 $1$）。$\mathbf{w}^T\mathbf{v} = 1$。

    $A^k \to \mathbf{v}\mathbf{w}^T = \begin{pmatrix} 3/8 & 3/8 \\ 5/8 & 5/8 \end{pmatrix}$。验证：$A^2 = \begin{pmatrix} 0.40 & 0.60 \\ 0.36 & 0.64 \end{pmatrix}$，$A^{10} \approx \begin{pmatrix} 0.375 & 0.375 \\ 0.625 & 0.625 \end{pmatrix}$（误差来自 $0.2^{10} \approx 10^{-7}$）。

    等一下，让我重新计算。$A\mathbf{v} = \mathbf{v}$，$\begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix}\begin{pmatrix} a \\ b \end{pmatrix} = \begin{pmatrix} a \\ b \end{pmatrix}$，$0.5a + 0.5b = a$，$0.3a + 0.7b = b$，均给出 $a = b$... 不对。$-0.5a + 0.5b = 0$ 给出 $a = b$。但 $\mathbf{w}^TA = \mathbf{w}^T$，$\begin{pmatrix} w_1 & w_2 \end{pmatrix}\begin{pmatrix} 0.5 & 0.5 \\ 0.3 & 0.7 \end{pmatrix} = \begin{pmatrix} w_1 & w_2 \end{pmatrix}$，$0.5w_1 + 0.3w_2 = w_1$ 给出 $w_2 = \frac{5}{3}w_1$。取 $\mathbf{w}^T = (3, 5)$，则 $\mathbf{w}^T\mathbf{v} = 3+5 = 8$（若 $\mathbf{v} = (1,1)^T$）。归一化：$\mathbf{v} = (1,1)^T$，$\mathbf{w}^T = (3/8, 5/8)$，$\mathbf{w}^T\mathbf{v} = 1$。

    $\mathbf{v}\mathbf{w}^T = \begin{pmatrix} 3/8 & 5/8 \\ 3/8 & 5/8 \end{pmatrix}$。这就是 $A^k$ 的极限。

---

## 17.6 随机矩阵

!!! definition "定义 17.6 (随机矩阵)"
    设 $P = (p_{ij}) \in \mathbb{R}^{n \times n}$。

    - 若 $P \geq O$ 且每行之和为 $1$（即 $P\mathbf{1} = \mathbf{1}$），则称 $P$ 为**行随机矩阵（row stochastic matrix）**；
    - 若 $P \geq O$ 且每列之和为 $1$（即 $\mathbf{1}^TP = \mathbf{1}^T$），则称 $P$ 为**列随机矩阵（column stochastic matrix）**；
    - 若 $P$ 既是行随机又是列随机矩阵，则称为**双随机矩阵（doubly stochastic matrix）**。

!!! proposition "命题 17.2 (随机矩阵的谱性质)"
    (1) 行随机矩阵的谱半径为 $1$，且 $\mathbf{1} = (1, 1, \ldots, 1)^T$ 是对应于特征值 $1$ 的右特征向量；

    (2) 所有特征值 $\lambda$ 满足 $|\lambda| \leq 1$；

    (3) 双随机矩阵的左 Perron 向量为 $\frac{1}{n}\mathbf{1}^T$。

!!! theorem "定理 17.6 (Birkhoff 定理)"
    $n \times n$ 双随机矩阵的集合是一个凸多面体（称为 **Birkhoff 多面体**），其顶点恰好是全体 $n \times n$ 置换矩阵。即每个双随机矩阵都可以表示为置换矩阵的凸组合。

??? proof "证明"
    **极点是置换矩阵**：设 $P$ 是双随机矩阵，且是 Birkhoff 多面体的极点。若 $P$ 不是置换矩阵，则 $P$ 有至少两个非零非一的元素。利用 König 定理可以找到两个不同的置换矩阵 $Q_1, Q_2$ 使得 $P = \alpha Q_1 + (1-\alpha)Q_2$（$0 < \alpha < 1$），矛盾于 $P$ 是极点。

    **置换矩阵是极点**：若置换矩阵 $Q = \alpha P_1 + (1-\alpha)P_2$，由 $Q$ 的元素为 $0$ 或 $1$ 以及 $P_1, P_2$ 的元素在 $[0,1]$ 中，可推出 $P_1 = P_2 = Q$。

    由 Krein-Milman 定理（或有限维的 Minkowski 定理），凸多面体等于其极点的凸包。$\blacksquare$

!!! example "例 17.7"
    双随机矩阵 $P = \begin{pmatrix} 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \\ 1/3 & 1/3 & 1/3 \end{pmatrix} = \frac{1}{3}(I + Q + Q^2)$，其中 $Q = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$，是三个置换矩阵 $I, Q, Q^2$ 的凸组合（各权重 $1/3$）。

---

## 17.7 Markov 链

!!! definition "定义 17.7 (Markov 链与转移矩阵)"
    设 $\{X_t\}_{t=0,1,2,\ldots}$ 是取值在有限状态空间 $\{1, 2, \ldots, n\}$ 上的随机过程。若

    $$\Pr(X_{t+1} = j \mid X_t = i, X_{t-1}, \ldots, X_0) = \Pr(X_{t+1} = j \mid X_t = i) = p_{ij},$$

    则 $\{X_t\}$ 是**（齐次）Markov 链（Markov chain）**，$P = (p_{ij})$ 是其**转移矩阵（transition matrix）**。$P$ 是行随机矩阵。

!!! definition "定义 17.8 (平稳分布)"
    行向量 $\boldsymbol{\pi}^T = (\pi_1, \ldots, \pi_n)$ 称为 Markov 链的**平稳分布（stationary distribution）**，若

    $$\boldsymbol{\pi}^T P = \boldsymbol{\pi}^T, \quad \boldsymbol{\pi} \geq \mathbf{0}, \quad \sum_i \pi_i = 1.$$

    即 $\boldsymbol{\pi}$ 是 $P^T$ 的对应于特征值 $1$ 的非负归一化左特征向量。

!!! theorem "定理 17.7 (Markov 链的收敛定理)"
    设 $P$ 是不可约且非周期（即本原）的转移矩阵。则：

    (1) 存在唯一的平稳分布 $\boldsymbol{\pi}^T$，且 $\boldsymbol{\pi} > \mathbf{0}$；

    (2) 对任意初始分布 $\mathbf{p}_0^T$，$\mathbf{p}_0^T P^k \to \boldsymbol{\pi}^T$（$k \to \infty$）；

    (3) $P^k \to \mathbf{1}\boldsymbol{\pi}^T$（每一行都收敛到 $\boldsymbol{\pi}^T$）。

??? proof "证明"
    由 $P$ 是本原的行随机矩阵，Perron 根 $\rho(P) = 1$，对应的右 Perron 向量为 $\mathbf{1}$。左 Perron 向量 $\boldsymbol{\pi}^T$（归一化使 $\sum\pi_i = 1$）就是平稳分布。

    由本原矩阵的收敛性定理（定理 17.5），$(P/\rho)^k = P^k \to \mathbf{1}\boldsymbol{\pi}^T$（因为 $\rho = 1$）。

    对任意初始分布 $\mathbf{p}_0^T$，$\mathbf{p}_0^TP^k \to \mathbf{p}_0^T\mathbf{1}\boldsymbol{\pi}^T = 1 \cdot \boldsymbol{\pi}^T = \boldsymbol{\pi}^T$。$\blacksquare$

!!! example "例 17.8"
    **天气模型**。设状态 $1$ = 晴天，$2$ = 雨天。转移矩阵

    $$P = \begin{pmatrix} 0.9 & 0.1 \\ 0.5 & 0.5 \end{pmatrix}.$$

    $P > O$，因此本原。平稳分布满足 $\boldsymbol{\pi}^TP = \boldsymbol{\pi}^T$：

    $0.9\pi_1 + 0.5\pi_2 = \pi_1$，$0.1\pi_1 + 0.5\pi_2 = \pi_2$，$\pi_1 + \pi_2 = 1$。

    第一个方程给出 $\pi_2 = 0.2\pi_1$，代入 $\pi_1 + 0.2\pi_1 = 1$，$\pi_1 = 5/6$，$\pi_2 = 1/6$。

    长期来看，约 $83.3\%$ 的天数是晴天，$16.7\%$ 是雨天。

---

## 17.8 非负矩阵的谱性质

!!! definition "定义 17.9 (循环指数)"
    不可约非负矩阵 $A$ 的**循环指数（index of imprimitivity）**（或**周期**）$h$ 定义为关联有向图 $G(A)$ 中所有回路长度的最大公约数。

    若 $h = 1$，$A$ 是本原的；若 $h > 1$，$A$ 是**周期的（cyclic）** 或**非本原的（imprimitive）**。

!!! theorem "定理 17.8 (谱的旋转对称性)"
    设 $A \geq O$ 是不可约的，循环指数为 $h$。则 $A$ 的谱关于原点具有 $h$ 重旋转对称性，即若 $\lambda$ 是 $A$ 的特征值，则 $\lambda e^{2\pi i/h}$ 也是 $A$ 的特征值。

    等价地，$A$ 的特征多项式仅含 $\lambda^h$ 的幂次，即 $p(\lambda) = \lambda^r q(\lambda^h)$ 对某个多项式 $q$。

??? proof "证明"
    由 Perron-Frobenius 定理，模等于 $\rho(A)$ 的特征值为 $\rho(A)e^{2k\pi i/h}$（$k = 0, 1, \ldots, h-1$）。

    更一般地，设 $\omega = e^{2\pi i/h}$，存在对角矩阵 $D = \operatorname{diag}(\omega^{d_1}, \ldots, \omega^{d_n})$（$d_i$ 取决于顶点 $i$ 在循环分类中的位置）使得

    $$D^{-1}AD = \omega A.$$

    这意味着 $A$ 和 $\omega A$ 相似，因此有相同的谱。即 $\lambda$ 是特征值 $\Rightarrow$ $\omega\lambda$ 也是特征值。$\blacksquare$

!!! proposition "命题 17.3 (循环标准形)"
    设 $A \geq O$ 不可约，循环指数 $h$。则存在置换矩阵 $P$ 使得

    $$P^TAP = \begin{pmatrix} O & A_{12} & O & \cdots & O \\ O & O & A_{23} & \cdots & O \\ \vdots & & & \ddots & \vdots \\ O & O & O & \cdots & A_{h-1,h} \\ A_{h1} & O & O & \cdots & O \end{pmatrix},$$

    即 $A$ 可以置换为循环分块形式，非零块只出现在"超对角"和左下角位置。

!!! example "例 17.9"
    $A = \begin{pmatrix} 0 & 2 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 0 & 0 & 1 \\ 4 & 0 & 0 & 0 \end{pmatrix}$。

    有向图：$1 \to 2 \to 3 \to 4 \to 1$，唯一回路长度为 $4$，循环指数 $h = 4$。

    特征多项式：$\lambda^4 - 24 = 0$，$\lambda = \sqrt[4]{24}\cdot\omega^k$（$k = 0,1,2,3$），其中 $\omega = i$。

    四个特征值 $\sqrt[4]{24}, i\sqrt[4]{24}, -\sqrt[4]{24}, -i\sqrt[4]{24}$ 在复平面上构成正方形，体现了 $4$ 重旋转对称性。

!!! example "例 17.10"
    **PageRank 算法的数学模型**。设 $n$ 个网页的链接关系由有向图描述，邻接矩阵的列归一化得到列随机矩阵 $H$（$H$ 的第 $j$ 列表示从页面 $j$ 出发的链接概率分布）。Google 矩阵为

    $$G = \alpha H + (1-\alpha)\frac{1}{n}\mathbf{1}\mathbf{1}^T,$$

    其中 $\alpha \in (0, 1)$（通常取 $\alpha = 0.85$）。

    $G$ 是列随机正矩阵（$G > O$），因此本原。由 Perron-Frobenius 定理，$\rho(G) = 1$ 是代数单根，对应的正特征向量 $\boldsymbol{\pi}$（归一化后）就是各页面的 PageRank 值。第二大特征值模 $|\lambda_2| \leq \alpha < 1$，保证了幂迭代的快速收敛。
