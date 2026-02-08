# 第 12 章 Jordan 标准形

并非所有矩阵都可以对角化。当一个矩阵的特征值的几何重数小于代数重数时，它不存在由特征向量组成的基，因而无法对角化。然而，每个方阵都可以化为一种"几乎对角"的标准形式——**Jordan 标准形**（Jordan Normal Form / Jordan Canonical Form）。Jordan 标准形是对角化的自然推广，它在矩阵理论、微分方程、控制论等领域具有核心地位。本章将从不变子空间出发，逐步构建 Jordan 标准形理论。

---

## 12.1 不变子空间

不变子空间（invariant subspace）是理解线性变换结构的基本工具，也是 Jordan 分解理论的起点。

!!! definition "定义 12.1 (不变子空间 Invariant Subspace)"
    设 $A$ 为 $n \times n$ 矩阵（或 $T$ 为线性空间 $V$ 上的线性变换），$W$ 为 $\mathbb{R}^n$（或 $V$）的子空间。若对任意 $\mathbf{w} \in W$ 都有 $A\mathbf{w} \in W$（即 $A(W) \subseteq W$），则称 $W$ 为 $A$ 的**不变子空间**。

!!! example "例 12.1"
    设 $A$ 为 $n \times n$ 矩阵，以下都是 $A$ 的不变子空间：

    1. $\{0\}$ 和 $\mathbb{R}^n$（平凡不变子空间）；
    2. $\ker(A) = \{\mathbf{x} : A\mathbf{x} = \mathbf{0}\}$（零空间）；
    3. $\operatorname{Im}(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\}$（像空间 / 列空间）；
    4. 特征空间 $E_\lambda = \ker(A - \lambda I)$；
    5. $\operatorname{span}\{\mathbf{v}\}$，其中 $\mathbf{v}$ 为 $A$ 的特征向量。

    **验证 (2)：** 若 $A\mathbf{x} = \mathbf{0}$，则 $A(A\mathbf{x}) = A\mathbf{0} = \mathbf{0}$，故 $A\mathbf{x} = \mathbf{0} \in \ker(A)$。

!!! definition "定义 12.2 (不变子空间直和分解)"
    设 $W_1, W_2, \ldots, W_k$ 都是 $A$ 的不变子空间，若
    $$
    \mathbb{R}^n = W_1 \oplus W_2 \oplus \cdots \oplus W_k
    $$
    （即每个 $\mathbf{x} \in \mathbb{R}^n$ 可以唯一表示为 $\mathbf{x} = \mathbf{w}_1 + \cdots + \mathbf{w}_k$，$\mathbf{w}_i \in W_i$），则称 $\mathbb{R}^n$ 为 $A$ 的不变子空间的**直和分解**。

!!! theorem "定理 12.1 (不变子空间与分块对角化)"
    设 $\mathbb{R}^n = W_1 \oplus W_2 \oplus \cdots \oplus W_k$ 为 $A$ 的不变子空间直和分解，取 $W_i$ 的基 $\mathcal{B}_i$，合并为 $\mathcal{B} = \mathcal{B}_1 \cup \cdots \cup \mathcal{B}_k$。设 $P$ 为基 $\mathcal{B}$ 的矩阵（列为各基向量），则
    $$
    P^{-1}AP = \begin{pmatrix} A_1 & & \\ & A_2 & \\ & & \ddots & \\ & & & A_k \end{pmatrix},
    $$
    其中 $A_i$ 为 $A$ 限制在 $W_i$ 上的矩阵表示（$\dim W_i \times \dim W_i$）。

??? proof "证明"
    设 $\dim W_i = n_i$，$\mathcal{B}_i = \{\mathbf{w}_{i,1}, \ldots, \mathbf{w}_{i,n_i}\}$。由于 $W_i$ 是 $A$ 的不变子空间，$A\mathbf{w}_{i,j} \in W_i$，因此 $A\mathbf{w}_{i,j}$ 可以表示为 $\mathcal{B}_i$ 中向量的线性组合。

    设 $P = [\mathbf{w}_{1,1}, \ldots, \mathbf{w}_{1,n_1}, \ldots, \mathbf{w}_{k,1}, \ldots, \mathbf{w}_{k,n_k}]$，则
    $$
    AP = P \begin{pmatrix} A_1 & & \\ & \ddots & \\ & & A_k \end{pmatrix},
    $$
    其中 $A_i$ 的第 $j$ 列是 $A\mathbf{w}_{i,j}$ 在基 $\mathcal{B}_i$ 下的坐标。由 $P$ 可逆即得结论。$\blacksquare$

!!! theorem "定理 12.2 (特征空间是不变子空间)"
    设 $\lambda$ 为矩阵 $A$ 的特征值，则特征空间 $E_\lambda = \ker(A - \lambda I)$ 是 $A$ 的不变子空间。更一般地，对任意多项式 $p$，$\ker(p(A))$ 是 $A$ 的不变子空间。

??? proof "证明"
    设 $\mathbf{v} \in E_\lambda$，即 $A\mathbf{v} = \lambda\mathbf{v}$。则 $A(\lambda\mathbf{v}) = \lambda(A\mathbf{v}) = \lambda^2\mathbf{v}$。但我们要证的是 $A\mathbf{v} \in E_\lambda$：$(A - \lambda I)(A\mathbf{v}) = A(A\mathbf{v}) - \lambda(A\mathbf{v}) = A(\lambda\mathbf{v}) - \lambda(\lambda\mathbf{v}) = \lambda^2\mathbf{v} - \lambda^2\mathbf{v} = \mathbf{0}$。因此 $A\mathbf{v} \in E_\lambda$。

    对一般情况：设 $p(A)\mathbf{v} = \mathbf{0}$，由于 $A$ 与 $p(A)$ 可交换，$p(A)(A\mathbf{v}) = A(p(A)\mathbf{v}) = A\mathbf{0} = \mathbf{0}$，故 $A\mathbf{v} \in \ker(p(A))$。$\blacksquare$

---

## 12.2 广义特征向量

当特征值的几何重数小于代数重数时，特征向量不足以构成一组基。广义特征向量（generalized eigenvectors）通过放宽条件来弥补这一不足。

!!! definition "定义 12.3 (广义特征向量 Generalized Eigenvector)"
    设 $\lambda$ 为矩阵 $A$ 的特征值。若非零向量 $\mathbf{v}$ 满足
    $$
    (A - \lambda I)^k \mathbf{v} = \mathbf{0}
    $$
    对某个正整数 $k$ 成立，则称 $\mathbf{v}$ 为 $A$ 对应于 $\lambda$ 的**广义特征向量**（秩为 $k$，如果 $(A-\lambda I)^{k-1}\mathbf{v} \neq \mathbf{0}$）。

!!! definition "定义 12.4 (广义特征空间 Generalized Eigenspace)"
    矩阵 $A$ 对应于特征值 $\lambda$ 的**广义特征空间**为
    $$
    G_\lambda = \ker(A - \lambda I)^n = \{\mathbf{v} \in \mathbb{R}^n : (A - \lambda I)^n \mathbf{v} = \mathbf{0}\},
    $$
    其中 $n$ 为矩阵的阶数。等价地，$G_\lambda = \bigcup_{k=1}^{\infty} \ker(A - \lambda I)^k$。

!!! theorem "定理 12.3 (广义特征空间的维数)"
    设 $\lambda$ 为 $A$ 的特征值，代数重数为 $m$，则
    $$
    \dim G_\lambda = m.
    $$
    即广义特征空间的维数等于特征值的代数重数。

??? proof "证明"
    由 Jordan 标准形定理（定理 12.8），$A$ 相似于 Jordan 矩阵 $J$。设 $A = PJP^{-1}$，则 $(A - \lambda I)^n = P(J - \lambda I)^n P^{-1}$。对 Jordan 矩阵 $J$，$\ker(J - \lambda I)^n$ 恰好由所有对应于 $\lambda$ 的 Jordan 块的列空间张成，其维数为 $\lambda$ 的代数重数 $m$。由相似变换保维数，$\dim G_\lambda = m$。$\blacksquare$

!!! theorem "定理 12.4 (广义特征空间的直和分解)"
    设 $A$ 为 $n \times n$ 矩阵（在复数域上），特征值为 $\lambda_1, \ldots, \lambda_s$（互不相同），则
    $$
    \mathbb{C}^n = G_{\lambda_1} \oplus G_{\lambda_2} \oplus \cdots \oplus G_{\lambda_s}.
    $$

??? proof "证明"
    **维数验证：** 由定理 12.3，$\sum \dim G_{\lambda_i} = \sum m_i = n$（代数重数之和等于 $n$）。

    **直和性：** 需证明若 $\mathbf{v}_1 + \cdots + \mathbf{v}_s = \mathbf{0}$（$\mathbf{v}_i \in G_{\lambda_i}$），则 $\mathbf{v}_i = \mathbf{0}$。设 $(A - \lambda_i I)^{m_i} \mathbf{v}_i = \mathbf{0}$。对 $i \neq j$，令 $p_j(x) = \prod_{i \neq j}(x - \lambda_i)^{m_i}$。则 $p_j(A)\mathbf{v}_i = \mathbf{0}$（$i \neq j$），而 $p_j(A)\mathbf{v}_j = c_j \mathbf{v}_j$（其中 $c_j = \prod_{i \neq j}(\lambda_j - \lambda_i)^{m_i} \neq 0$）。

    对 $\mathbf{v}_1 + \cdots + \mathbf{v}_s = \mathbf{0}$ 两边作用 $p_j(A)$：$c_j \mathbf{v}_j = \mathbf{0}$，故 $\mathbf{v}_j = \mathbf{0}$。$\blacksquare$

!!! example "例 12.2"
    设 $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$，求 $A$ 的广义特征空间。

    **解：** $A$ 的特征值为 $\lambda = 2$（代数重数 2）。

    $(A - 2I) = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$\ker(A - 2I) = \operatorname{span}\left\{\begin{pmatrix}1\\0\end{pmatrix}\right\}$，维数为 1（几何重数）。

    $(A - 2I)^2 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$，$\ker(A - 2I)^2 = \mathbb{R}^2$，维数为 2。

    因此 $G_2 = \mathbb{R}^2$。广义特征向量 $\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$ 满足 $(A-2I)\mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix} \neq \mathbf{0}$，但 $(A-2I)^2\mathbf{v}_2 = \mathbf{0}$，故 $\mathbf{v}_2$ 是秩为 2 的广义特征向量。

---

## 12.3 幂零矩阵

幂零矩阵（nilpotent matrix）是理解 Jordan 标准形的关键。每个 Jordan 块可以分解为一个标量矩阵加一个幂零矩阵。

!!! definition "定义 12.5 (幂零矩阵 Nilpotent Matrix)"
    若存在正整数 $k$ 使得 $N^k = 0$，则称矩阵 $N$ 为**幂零矩阵**。满足 $N^k = 0$ 的最小正整数 $k$ 称为 $N$ 的**幂零指数**（nilpotency index）。

!!! theorem "定理 12.5 (幂零矩阵的性质)"
    设 $N$ 为 $n \times n$ 幂零矩阵，幂零指数为 $k$，则：

    1. $N$ 的所有特征值都为 $0$；
    2. $\operatorname{tr}(N) = 0$，$\det(N) = 0$；
    3. $k \le n$；
    4. $N$ 的特征多项式为 $p(\lambda) = \lambda^n$；
    5. $I - N$ 可逆，且 $(I - N)^{-1} = I + N + N^2 + \cdots + N^{k-1}$。

??? proof "证明"
    **(1)** 设 $\lambda$ 为 $N$ 的特征值，$\mathbf{v}$ 为对应特征向量。则 $N\mathbf{v} = \lambda\mathbf{v}$，$N^k\mathbf{v} = \lambda^k\mathbf{v}$。由 $N^k = 0$，得 $\lambda^k\mathbf{v} = \mathbf{0}$。因为 $\mathbf{v} \neq \mathbf{0}$，所以 $\lambda^k = 0$，即 $\lambda = 0$。

    **(2)** 由 (1)，所有特征值为零，故 $\operatorname{tr}(N) = 0$（特征值之和），$\det(N) = 0$（特征值之积）。

    **(3)** 考虑子空间链 $\{0\} \subseteq \ker(N) \subseteq \ker(N^2) \subseteq \cdots$。若 $\ker(N^j) = \ker(N^{j+1})$，则对所有 $i \ge j$，$\ker(N^i) = \ker(N^j)$。由 $\ker(N^k) = \mathbb{R}^n$（因为 $N^k = 0$），且维数严格递增直到稳定，$k \le n$。

    **(4)** 由 (1)，$N$ 的特征多项式只有零根，故 $p(\lambda) = \lambda^n$。

    **(5)** $(I - N)(I + N + \cdots + N^{k-1}) = I - N^k = I$。$\blacksquare$

!!! example "例 12.3"
    设 $N = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 0 & 0 & 0 \end{pmatrix}$。验证 $N$ 是幂零矩阵并求其幂零指数。

    **解：**
    $$
    N^2 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}, \qquad N^3 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
    $$
    $N^2 \neq 0$ 而 $N^3 = 0$，故幂零指数为 $3$。

!!! theorem "定理 12.6 (幂零矩阵的 Jordan 标准形)"
    设 $N$ 为 $n \times n$ 幂零矩阵，幂零指数为 $k$。则 $N$ 相似于分块对角矩阵
    $$
    J = \operatorname{diag}(J_{n_1}(0), J_{n_2}(0), \ldots, J_{n_t}(0)),
    $$
    其中 $J_{n_i}(0)$ 为 $n_i \times n_i$ 的幂零 Jordan 块（对角线为 0，超对角线为 1），$n_1 \ge n_2 \ge \cdots \ge n_t \ge 1$，$\sum n_i = n$，$n_1 = k$。

??? proof "证明"
    这是 Jordan 标准形定理（定理 12.8）在 $\lambda = 0$ 情形的特例。幂零指数 $k$ 等于最大 Jordan 块的阶数 $n_1$，因为 $J_{n_1}(0)^{n_1} = 0$ 但 $J_{n_1}(0)^{n_1-1} \neq 0$。详细证明见定理 12.8。

---

## 12.4 Jordan 块

Jordan 块（Jordan block）是构成 Jordan 标准形的基本单元。

!!! definition "定义 12.6 (Jordan 块 Jordan Block)"
    $k \times k$ 的 **Jordan 块** $J_k(\lambda)$ 定义为
    $$
    J_k(\lambda) = \begin{pmatrix}
    \lambda & 1 & 0 & \cdots & 0 \\
    0 & \lambda & 1 & \cdots & 0 \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & \lambda & 1 \\
    0 & \cdots & 0 & 0 & \lambda
    \end{pmatrix} = \lambda I_k + N_k,
    $$
    其中 $N_k$ 为 $k \times k$ 的**基本幂零矩阵**（超对角线元素全为 1，其余为 0）。

!!! theorem "定理 12.7 (Jordan 块的性质)"
    $J_k(\lambda)$ 具有以下性质：

    1. 唯一特征值为 $\lambda$，代数重数 $k$，几何重数 $1$；
    2. $(J_k(\lambda) - \lambda I)^k = N_k^k = 0$，但 $(J_k(\lambda) - \lambda I)^{k-1} \neq 0$；
    3. $J_k(\lambda)^m = \sum_{j=0}^{k-1} \binom{m}{j} \lambda^{m-j} N_k^j$（$m \ge k-1$）；
    4. 当 $\lambda \neq 0$ 时，$J_k(\lambda)$ 可逆，$J_k(\lambda)^{-1} = \frac{1}{\lambda}(I + \frac{1}{\lambda}N_k)^{-1}$。

??? proof "证明"
    **(1)** $\det(J_k(\lambda) - \mu I) = (\lambda - \mu)^k$，故唯一特征值为 $\lambda$，代数重数 $k$。$(J_k(\lambda) - \lambda I)\mathbf{v} = N_k\mathbf{v} = \mathbf{0}$ 的解空间为 $\operatorname{span}\{\mathbf{e}_1\}$，几何重数为 1。

    **(2)** $J_k(\lambda) - \lambda I = N_k$，而 $N_k^k = 0$、$N_k^{k-1} \neq 0$（$N_k^{k-1}$ 的 $(1,k)$ 元素为 1）。

    **(3)** 由二项式定理，$J_k(\lambda)^m = (\lambda I + N_k)^m = \sum_{j=0}^{m} \binom{m}{j} \lambda^{m-j} N_k^j$。由于 $N_k^j = 0$（$j \ge k$），求和实际只到 $k-1$。

    **(4)** $J_k(\lambda) = \lambda(I + \frac{1}{\lambda}N_k)$。由定理 12.5(5)，$I + \frac{1}{\lambda}N_k$ 可逆。$\blacksquare$

!!! example "例 12.4"
    计算 $J_3(2)^4$。

    **解：** $J_3(2) = 2I + N_3$，其中 $N_3 = \begin{pmatrix}0&1&0\\0&0&1\\0&0&0\end{pmatrix}$。

    由公式：
    $$
    J_3(2)^4 = \sum_{j=0}^{2} \binom{4}{j} 2^{4-j} N_3^j = \binom{4}{0}2^4 I + \binom{4}{1}2^3 N_3 + \binom{4}{2}2^2 N_3^2
    $$
    $$
    = 16I + 32N_3 + 24N_3^2 = \begin{pmatrix}16&32&24\\0&16&32\\0&0&16\end{pmatrix}.
    $$

---

## 12.5 Jordan 标准形定理

!!! definition "定义 12.7 (Jordan 矩阵 Jordan Matrix)"
    **Jordan 矩阵**是分块对角矩阵
    $$
    J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_s}(\lambda_s)),
    $$
    其中每个 $J_{k_i}(\lambda_i)$ 是 Jordan 块。不同块可以有相同的 $\lambda$ 值。

!!! theorem "定理 12.8 (Jordan 标准形定理)"
    设 $A$ 为 $n \times n$ 复数矩阵。则存在可逆矩阵 $P$ 使得
    $$
    P^{-1}AP = J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s)),
    $$
    其中 $J$ 为 Jordan 矩阵。

    而且，在**不计 Jordan 块的排列顺序**的意义下，$J$ 是唯一的。即 Jordan 块的大小和对应的特征值是由 $A$ 唯一确定的。

??? proof "证明"
    **证明思路：**（完整证明较长，这里给出关键步骤）

    **存在性：**

    第一步（广义特征空间分解）：由定理 12.4，$\mathbb{C}^n = G_{\lambda_1} \oplus \cdots \oplus G_{\lambda_s}$。每个 $G_{\lambda_i}$ 是 $A$ 的不变子空间，$A$ 限制在 $G_{\lambda_i}$ 上的矩阵为 $A|_{G_{\lambda_i}}$。

    第二步（化为幂零情形）：在 $G_{\lambda_i}$ 上，$A|_{G_{\lambda_i}} - \lambda_i I$ 是幂零的（因为 $(A - \lambda_i I)^{m_i}|_{G_{\lambda_i}} = 0$，其中 $m_i$ 为 $\lambda_i$ 的代数重数）。

    第三步（幂零矩阵的 Jordan 形）：设 $N$ 为 $d \times d$ 幂零矩阵。通过构造 Jordan 链（Jordan chains）可以找到基使得 $N$ 在此基下为 Jordan 形。

    具体地，选取 $\mathbf{v}$ 使得 $(A-\lambda I)^{k-1}\mathbf{v} \neq \mathbf{0}$ 但 $(A-\lambda I)^k\mathbf{v} = \mathbf{0}$。定义 Jordan 链：
    $$
    \mathbf{v}, (A-\lambda I)\mathbf{v}, (A-\lambda I)^2\mathbf{v}, \ldots, (A-\lambda I)^{k-1}\mathbf{v}.
    $$
    这 $k$ 个向量线性无关，且 $A$ 在它们张成的子空间上表示为 $J_k(\lambda)$。

    反复应用此过程，可将整个广义特征空间分解为 Jordan 链的直和。

    **唯一性：**

    Jordan 块的结构由以下不变量确定：对每个特征值 $\lambda$ 和每个正整数 $j$，
    $$
    r_j = \operatorname{rank}(A - \lambda I)^j
    $$
    是相似不变量。大小为 $k$ 的 Jordan 块 $J_k(\lambda)$ 的个数为
    $$
    n_k = r_{k-2} - 2r_{k-1} + r_k
    $$
    （约定 $r_{-1} = n$, $r_0 = n - \dim\ker(A-\lambda I)$），这由 $A$ 唯一确定。$\blacksquare$

!!! theorem "定理 12.9 (Jordan 标准形与矩阵性质)"
    设 $A$ 的 Jordan 标准形为 $J = \operatorname{diag}(J_{k_1}(\lambda_1), \ldots, J_{k_s}(\lambda_s))$，则：

    1. $A$ 可对角化当且仅当所有 Jordan 块的大小为 $1$（即 $k_i = 1$）；
    2. 特征值 $\lambda$ 的几何重数等于对应于 $\lambda$ 的 Jordan 块的个数；
    3. 特征值 $\lambda$ 的代数重数等于对应于 $\lambda$ 的所有 Jordan 块大小之和；
    4. $\det(A) = \prod_{i=1}^s \lambda_i^{k_i}$，$\operatorname{tr}(A) = \sum_{i=1}^s k_i \lambda_i$。

??? proof "证明"
    **(1)** $A$ 可对角化 $\Leftrightarrow$ 对每个特征值，几何重数 = 代数重数 $\Leftrightarrow$ 每个 Jordan 块都是 $1 \times 1$ 的。

    **(2)** 特征值 $\lambda$ 对应的 Jordan 块为 $J_{k_1}(\lambda), \ldots, J_{k_t}(\lambda)$。每个 $J_{k_i}(\lambda)$ 的特征空间维数为 1，因此总几何重数为 $t$（块的个数）。

    **(3)** 各 Jordan 块 $J_{k_i}(\lambda)$ 贡献代数重数 $k_i$，总代数重数为 $\sum k_i$。

    **(4)** 相似变换不改变行列式和迹。$\det(J) = \prod \det(J_{k_i}(\lambda_i)) = \prod \lambda_i^{k_i}$。$\operatorname{tr}(J) = \sum k_i\lambda_i$。$\blacksquare$

!!! example "例 12.5"
    一个 $5 \times 5$ 矩阵 $A$ 的特征值为 $\lambda_1 = 2$（代数重数 3，几何重数 2）和 $\lambda_2 = -1$（代数重数 2，几何重数 1）。确定 $A$ 的 Jordan 标准形。

    **解：**

    对 $\lambda_1 = 2$：代数重数 3，几何重数 2，因此有 2 个 Jordan 块，大小之和为 3。可能为 $(2,1)$（一个 $2\times 2$ 块和一个 $1\times 1$ 块）。

    对 $\lambda_2 = -1$：代数重数 2，几何重数 1，因此有 1 个 Jordan 块，大小为 2。

    Jordan 标准形为
    $$
    J = \begin{pmatrix}
    2 & 1 & & & \\
    0 & 2 & & & \\
    & & 2 & & \\
    & & & -1 & 1 \\
    & & & 0 & -1
    \end{pmatrix}.
    $$

---

## 12.6 最小多项式

最小多项式（minimal polynomial）与 Jordan 标准形有着密切联系。

!!! definition "定义 12.8 (最小多项式 Minimal Polynomial)"
    矩阵 $A$ 的**最小多项式** $m_A(\lambda)$ 是满足 $m_A(A) = 0$ 的次数最低的首一多项式。

!!! theorem "定理 12.10 (Cayley-Hamilton 定理与最小多项式)"
    设 $A$ 的特征多项式为 $p_A(\lambda)$，最小多项式为 $m_A(\lambda)$，则：

    1. $m_A(\lambda)$ 整除 $p_A(\lambda)$；
    2. $m_A(\lambda)$ 与 $p_A(\lambda)$ 有相同的根（即相同的特征值集合）；
    3. $A$ 可对角化当且仅当 $m_A(\lambda)$ 无重根。

??? proof "证明"
    **(1)** 由 Cayley-Hamilton 定理，$p_A(A) = 0$。由最小多项式的定义，$m_A$ 是使 $m_A(A) = 0$ 的最低次首一多项式。用多项式除法 $p_A = q \cdot m_A + r$，则 $r(A) = p_A(A) - q(A)m_A(A) = 0$。若 $r \neq 0$，则 $\deg r < \deg m_A$，矛盾。故 $r = 0$，$m_A | p_A$。

    **(2)** 设 $\lambda_0$ 是 $p_A$ 的根，$\mathbf{v}$ 为对应特征向量。$m_A(A)\mathbf{v} = m_A(\lambda_0)\mathbf{v} = \mathbf{0}$。因为 $\mathbf{v} \neq \mathbf{0}$，所以 $m_A(\lambda_0) = 0$。反向由 (1) 显然。

    **(3)** $A$ 可对角化 $\Leftrightarrow$ 所有 Jordan 块为 $1 \times 1$ $\Leftrightarrow$ $m_A(\lambda) = \prod(\lambda - \lambda_i)$（无重根）。

    详细来说，$m_A(\lambda) = \prod_{i=1}^s (\lambda - \lambda_i)^{d_i}$，其中 $d_i$ 是 $\lambda_i$ 对应的最大 Jordan 块的大小。$A$ 可对角化 $\Leftrightarrow$ 所有 $d_i = 1$ $\Leftrightarrow$ $m_A$ 无重根。$\blacksquare$

!!! theorem "定理 12.11 (最小多项式与 Jordan 形的关系)"
    设 $A$ 的 Jordan 标准形中，特征值 $\lambda_i$ 对应的最大 Jordan 块大小为 $d_i$，则
    $$
    m_A(\lambda) = \prod_{i=1}^{s} (\lambda - \lambda_i)^{d_i}.
    $$

??? proof "证明"
    设 $q(\lambda) = \prod_{i=1}^s (\lambda - \lambda_i)^{d_i}$。需证 $q(A) = 0$ 且 $q$ 是满足此条件的最低次首一多项式。

    $q(J) = \operatorname{diag}(q(J_{k_1}(\lambda_1)), \ldots)$。对每个 Jordan 块 $J_k(\lambda_i)$，$q(J_k(\lambda_i))$ 包含因子 $(J_k(\lambda_i) - \lambda_i I)^{d_i} = N_k^{d_i}$。由于 $k \le d_i$（$d_i$ 是最大块大小），$N_k^{d_i} = 0$，故 $q(J_k(\lambda_i)) = 0$。因此 $q(J) = 0$，$q(A) = Pq(J)P^{-1} = 0$。

    若 $\deg q$ 可以更低，则存在某个 $\lambda_i$ 使得 $(lambda-\lambda_i)$ 的幂次小于 $d_i$，但这会导致最大 Jordan 块 $J_{d_i}(\lambda_i)$ 不被零化，矛盾。$\blacksquare$

!!! example "例 12.6"
    设 $A = \begin{pmatrix}3&1&0\\0&3&0\\0&0&5\end{pmatrix}$，求其最小多项式。

    **解：** $A$ 已经是 Jordan 形：$J_2(3) \oplus J_1(5)$。

    特征值 $\lambda_1 = 3$（最大块大小 2），$\lambda_2 = 5$（最大块大小 1）。

    最小多项式：$m_A(\lambda) = (\lambda - 3)^2(\lambda - 5)$。

    验证：$p_A(\lambda) = (\lambda-3)^2(\lambda-5)$。此处最小多项式等于特征多项式。

---

## 12.7 Jordan 标准形的计算

本节通过系统方法和详细例题展示如何计算 Jordan 标准形和过渡矩阵。

!!! definition "定义 12.9 (Jordan 链 Jordan Chain)"
    设 $\lambda$ 为 $A$ 的特征值。如果向量序列 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k$ 满足
    $$
    (A - \lambda I)\mathbf{v}_1 = \mathbf{0}, \quad (A - \lambda I)\mathbf{v}_j = \mathbf{v}_{j-1}, \quad j = 2, \ldots, k,
    $$
    则称 $\{\mathbf{v}_k, \mathbf{v}_{k-1}, \ldots, \mathbf{v}_1\}$ 为一条**Jordan 链**，长度为 $k$。其中 $\mathbf{v}_1$ 是特征向量，$\mathbf{v}_2, \ldots, \mathbf{v}_k$ 是广义特征向量。

!!! definition "定义 12.10 (Jordan 标准形的计算步骤)"
    给定 $n \times n$ 矩阵 $A$，计算其 Jordan 标准形的系统方法：

    1. 求特征多项式，确定特征值 $\lambda_i$ 及其代数重数 $m_i$；
    2. 对每个 $\lambda_i$，逐步计算 $\ker(A-\lambda_i I)^j$（$j = 1, 2, \ldots$）的维数，确定 Jordan 块结构；
    3. 大小为 $k$ 的 Jordan 块个数 = $\dim\ker(A-\lambda I)^k - 2\dim\ker(A-\lambda I)^{k-1} + \dim\ker(A-\lambda I)^{k-2}$；
    4. 构造 Jordan 链以确定过渡矩阵 $P$。

!!! example "例 12.7"
    求矩阵 $A = \begin{pmatrix} 5 & 4 & 2 & 1 \\ 0 & 1 & -1 & -1 \\ -1 & -1 & 3 & 0 \\ 1 & 1 & -1 & 2 \end{pmatrix}$ 的 Jordan 标准形。

    **解：**

    **第一步：** 特征多项式（可通过行列式展开或其他方法计算）：
    $$
    p_A(\lambda) = (\lambda - 1)^2(\lambda - 4)(\lambda - 3).
    $$
    特征值：$\lambda_1 = 1$（代数重数 2），$\lambda_2 = 4$（代数重数 1），$\lambda_3 = 3$（代数重数 1）。

    **第二步：** 对 $\lambda_1 = 1$：
    $$
    A - I = \begin{pmatrix} 4&4&2&1 \\ 0&0&-1&-1 \\ -1&-1&2&0 \\ 1&1&-1&1 \end{pmatrix}.
    $$
    计算 $\operatorname{rank}(A-I)$。通过行化简可得 $\operatorname{rank}(A-I) = 3$，故 $\dim\ker(A-I) = 1$。

    几何重数为 1 < 代数重数 2，因此有一个 $2 \times 2$ Jordan 块 $J_2(1)$。

    对 $\lambda_2 = 4$ 和 $\lambda_3 = 3$：代数重数均为 1，因此各有一个 $1 \times 1$ 块。

    **第三步：** Jordan 标准形为
    $$
    J = \begin{pmatrix}
    1 & 1 & & \\
    0 & 1 & & \\
    & & 4 & \\
    & & & 3
    \end{pmatrix}.
    $$

!!! example "例 12.8"
    求矩阵 $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ 的 Jordan 标准形和过渡矩阵。

    **解：** $A$ 本身已经是 Jordan 形 $J_2(2)$。

    特征值 $\lambda = 2$（代数重数 2，几何重数 1）。

    $(A - 2I) = \begin{pmatrix}0&1\\0&0\end{pmatrix}$，$\ker(A-2I) = \operatorname{span}\left\{\begin{pmatrix}1\\0\end{pmatrix}\right\}$。

    特征向量 $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$。广义特征向量 $\mathbf{v}_2$ 满足 $(A-2I)\mathbf{v}_2 = \mathbf{v}_1$：
    $$
    \begin{pmatrix}0&1\\0&0\end{pmatrix}\mathbf{v}_2 = \begin{pmatrix}1\\0\end{pmatrix}, \quad \Rightarrow \quad \mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}.
    $$

    过渡矩阵 $P = [\mathbf{v}_1, \mathbf{v}_2] = \begin{pmatrix}1&0\\0&1\end{pmatrix} = I$。

    验证：$P^{-1}AP = A = J_2(2)$。

!!! example "例 12.9"
    求矩阵 $A = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 8 & -12 & 6 \end{pmatrix}$ 的 Jordan 标准形。

    **解：** 特征多项式为
    $$
    p_A(\lambda) = -\lambda^3 + 6\lambda^2 - 12\lambda + 8 = -(\lambda - 2)^3.
    $$
    唯一特征值 $\lambda = 2$，代数重数 3。

    $$
    A - 2I = \begin{pmatrix}-2&1&0\\0&-2&1\\8&-12&4\end{pmatrix}.
    $$

    行化简得 $\operatorname{rank}(A - 2I) = 2$，$\dim\ker(A-2I) = 1$，几何重数为 1。

    $(A-2I)^2$ 的秩：计算
    $$
    (A-2I)^2 = \begin{pmatrix}4&-4&1\\-8&8&-2\\-16&16&-4\end{pmatrix} + \cdots
    $$
    实际计算得 $\operatorname{rank}(A-2I)^2 = 1$，$\dim\ker(A-2I)^2 = 2$。

    $(A-2I)^3 = 0$，$\dim\ker(A-2I)^3 = 3$。

    大小为 3 的 Jordan 块个数 = $3 - 2 \times 2 + 1 = 0$...

    重新使用公式：Jordan 块计数由递增序列 $d_j = \dim\ker(A-\lambda I)^j$ 决定。$d_0 = 0, d_1 = 1, d_2 = 2, d_3 = 3$。增量为 $\Delta_1 = 1, \Delta_2 = 1, \Delta_3 = 1$。

    大小 $\ge k$ 的块的个数 = $\Delta_k$。因此大小 $\ge 1$ 的块有 1 个，大小 $\ge 2$ 的块有 1 个，大小 $\ge 3$ 的块有 1 个。故恰好有 1 个大小为 3 的 Jordan 块。

    Jordan 标准形：$J = J_3(2) = \begin{pmatrix}2&1&0\\0&2&1\\0&0&2\end{pmatrix}$。

---

## 12.8 Jordan 形的应用

### 12.8.1 矩阵幂的计算

!!! theorem "定理 12.12 (利用 Jordan 形计算矩阵幂)"
    设 $A = PJP^{-1}$，则 $A^n = PJ^nP^{-1}$。而
    $$
    J^n = \operatorname{diag}(J_{k_1}(\lambda_1)^n, \ldots, J_{k_s}(\lambda_s)^n),
    $$
    其中
    $$
    J_k(\lambda)^n = \begin{pmatrix}
    \lambda^n & \binom{n}{1}\lambda^{n-1} & \binom{n}{2}\lambda^{n-2} & \cdots & \binom{n}{k-1}\lambda^{n-k+1} \\
    0 & \lambda^n & \binom{n}{1}\lambda^{n-1} & \cdots & \binom{n}{k-2}\lambda^{n-k+2} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & \lambda^n & \binom{n}{1}\lambda^{n-1} \\
    0 & \cdots & 0 & 0 & \lambda^n
    \end{pmatrix}.
    $$

??? proof "证明"
    $J_k(\lambda)^n = (\lambda I + N_k)^n = \sum_{j=0}^{k-1}\binom{n}{j}\lambda^{n-j}N_k^j$。而 $N_k^j$ 的 $(p,q)$ 元素为 $\delta_{p+j,q}$（第 $j$ 条超对角线为 1）。因此 $J_k(\lambda)^n$ 的 $(p,q)$ 元素为：
    $$
    [J_k(\lambda)^n]_{pq} = \begin{cases} \binom{n}{q-p}\lambda^{n-q+p} & \text{if } q \ge p, \\ 0 & \text{if } q < p. \end{cases}
    $$
    其中约定 $\binom{n}{j} = 0$ 当 $j > n$ 或 $j < 0$。$\blacksquare$

!!! example "例 12.10"
    设 $A = \begin{pmatrix}3&1\\0&3\end{pmatrix}$，求 $A^{100}$。

    **解：** $A = J_2(3)$，因此
    $$
    A^{100} = J_2(3)^{100} = \begin{pmatrix} 3^{100} & 100 \cdot 3^{99} \\ 0 & 3^{100} \end{pmatrix}.
    $$

### 12.8.2 微分方程组

!!! theorem "定理 12.13 (Jordan 形与线性常微分方程组)"
    线性常系数微分方程组 $\mathbf{x}'(t) = A\mathbf{x}(t)$ 的解为 $\mathbf{x}(t) = e^{At}\mathbf{x}(0)$。利用 Jordan 分解 $A = PJP^{-1}$：
    $$
    e^{At} = Pe^{Jt}P^{-1} = P \operatorname{diag}(e^{J_{k_1}(\lambda_1)t}, \ldots, e^{J_{k_s}(\lambda_s)t}) P^{-1},
    $$
    其中
    $$
    e^{J_k(\lambda)t} = e^{\lambda t}\begin{pmatrix}
    1 & t & \frac{t^2}{2!} & \cdots & \frac{t^{k-1}}{(k-1)!} \\
    0 & 1 & t & \cdots & \frac{t^{k-2}}{(k-2)!} \\
    \vdots & & \ddots & \ddots & \vdots \\
    0 & \cdots & 0 & 1 & t \\
    0 & \cdots & 0 & 0 & 1
    \end{pmatrix}.
    $$

??? proof "证明"
    $e^{J_k(\lambda)t} = e^{(\lambda I + N_k)t} = e^{\lambda t I} e^{N_k t}$（因为 $\lambda I$ 和 $N_k$ 可交换）。

    $e^{\lambda t I} = e^{\lambda t} I$，$e^{N_k t} = \sum_{j=0}^{k-1} \frac{t^j}{j!} N_k^j$（因为 $N_k^k = 0$）。

    因此 $e^{J_k(\lambda)t} = e^{\lambda t} \sum_{j=0}^{k-1} \frac{t^j}{j!} N_k^j$，矩阵元素如上所述。$\blacksquare$

!!! example "例 12.11"
    求微分方程组 $\mathbf{x}' = \begin{pmatrix}2&1\\0&2\end{pmatrix}\mathbf{x}$，初始条件 $\mathbf{x}(0) = \begin{pmatrix}1\\3\end{pmatrix}$ 的解。

    **解：** $A = J_2(2)$，因此
    $$
    e^{At} = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}.
    $$

    $$
    \mathbf{x}(t) = e^{At}\mathbf{x}(0) = e^{2t}\begin{pmatrix}1&t\\0&1\end{pmatrix}\begin{pmatrix}1\\3\end{pmatrix} = e^{2t}\begin{pmatrix}1+3t\\3\end{pmatrix}.
    $$

    即 $x_1(t) = (1+3t)e^{2t}$，$x_2(t) = 3e^{2t}$。

---

## 本章小结

本章系统介绍了 Jordan 标准形理论：

1. **不变子空间**提供了将线性变换分块研究的框架；
2. **广义特征向量**弥补了特征向量不足的问题，$\dim G_\lambda$ 等于代数重数；
3. **Jordan 块** $J_k(\lambda) = \lambda I + N_k$ 是"几乎对角"的基本单元；
4. **Jordan 标准形定理**保证每个方阵（在复数域上）相似于唯一的 Jordan 矩阵；
5. **最小多项式**反映了最大 Jordan 块的大小，$A$ 可对角化当且仅当最小多项式无重根；
6. **应用**包括高效计算矩阵幂 $A^n$ 和求解线性常微分方程组。

Jordan 标准形是矩阵理论中最深刻的结果之一，为后续的矩阵函数理论奠定了基础。
