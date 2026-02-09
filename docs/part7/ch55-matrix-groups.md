# 第 55 章 矩阵群与经典 Lie 群

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch2) · 行列式 (Ch3) · 特征值 (Ch6) · 矩阵指数 (Ch13) · 正交矩阵/酉矩阵 (Ch7-Ch8)

**本章脉络**：矩阵群定义 → $\mathrm{GL}(n)$ → $\mathrm{SL}(n)$ → $\mathrm{O}(n)/\mathrm{SO}(n)$ → $\mathrm{U}(n)/\mathrm{SU}(n)$ → $\mathrm{Sp}(2n)$ → Lie 代数（切空间） → 指数映射 → 伴随表示 → Baker-Campbell-Hausdorff 公式 → $\mathfrak{sl}(2)$ 表示论 → 根系与 Dynkin 图 → Weyl 完全可约性 → Peter-Weyl 定理 → 例外 Lie 群

**延伸**：矩阵 Lie 群是微分几何、粒子物理（规范对称性群 $SU(3) \times SU(2) \times U(1)$）和机器人学（$SE(3)$ 刚体运动群）的核心数学结构

</div>

前几章中我们遇到了许多由特殊性质的矩阵构成的集合——正交矩阵、酉矩阵、辛矩阵等——它们对矩阵乘法构成群。这些"矩阵群"不仅是代数对象（群），还是几何对象（光滑流形），两者的统一就是 **Lie 群**的概念。本章系统梳理经典矩阵群，然后引入 Lie 代数——Lie 群在单位元处的切空间。Lie 代数是 Lie 群的"无穷小版本"，它用线性代数的工具（矩阵加法、矩阵括号 $[X,Y] = XY - YX$）来编码群的局部结构。指数映射 $\exp: \mathfrak{g} \to G$ 是连接 Lie 代数与 Lie 群的桥梁。最后，我们介绍伴随表示和 Baker-Campbell-Hausdorff 公式，它们揭示了群乘法在 Lie 代数层面的反映。

---

## 55.1 矩阵群的定义

<div class="context-flow" markdown>

**核心问题**：什么是矩阵群？如何从一般线性群中用方程"切出"经典群？

</div>

!!! definition "定义 55.1 (矩阵群)"
    一个**矩阵群**（matrix group）或**线性群**（linear group）是一般线性群 $\mathrm{GL}(n, \mathbb{F})$（$\mathbb{F} = \mathbb{R}$ 或 $\mathbb{C}$）的一个子群 $G$，即 $G \subseteq \mathrm{GL}(n, \mathbb{F})$ 满足：

    (i) $I_n \in G$。

    (ii) 若 $A, B \in G$，则 $AB \in G$。

    (iii) 若 $A \in G$，则 $A^{-1} \in G$。

!!! theorem "定理 55.1 (闭子群定理 / Cartan 定理)"
    $\mathrm{GL}(n, \mathbb{R})$（或 $\mathrm{GL}(n, \mathbb{C})$）的任何**闭**子群都是 Lie 群（即具有光滑流形结构，使得群运算是光滑的）。

??? proof "证明"
    这是 Lie 群理论的基础定理之一，完整证明见微分几何教材。这里给出思路。

    设 $G$ 是 $\mathrm{GL}(n)$ 的闭子群。定义
    $$\mathfrak{g} = \{X \in M_n : e^{tX} \in G, \forall t \in \mathbb{R}\}.$$

    **第一步**：证明 $\mathfrak{g}$ 是 $M_n$ 的线性子空间。若 $X, Y \in \mathfrak{g}$，则利用 Lie-Trotter 公式
    $$e^{t(X+Y)} = \lim_{m \to \infty} (e^{tX/m} e^{tY/m})^m,$$
    由 $G$ 的闭性知 $e^{t(X+Y)} \in G$，故 $X + Y \in \mathfrak{g}$。标量乘法的封闭性是直接的。

    **第二步**：证明 $\mathfrak{g}$ 对 Lie 括号封闭。利用
    $$[X, Y] = \lim_{t \to 0} \frac{e^{tX}e^{tY}e^{-tX}e^{-tY} - I}{t^2},$$
    以及更精确的
    $$e^{t^2[X,Y]} = \lim_{m \to \infty} (e^{tX/\sqrt{m}} e^{tY/\sqrt{m}} e^{-tX/\sqrt{m}} e^{-tY/\sqrt{m}})^m,$$
    可知 $e^{t[X,Y]} \in G$，故 $[X,Y] \in \mathfrak{g}$。

    **第三步**：在单位元附近，$\exp$ 将 $\mathfrak{g}$ 的邻域微分同胚到 $G$ 的邻域，由此赋予 $G$ 流形结构。$\blacksquare$

!!! remark "注记"
    闭子群定理的威力在于：要证明一个矩阵群是 Lie 群，只需验证它是 $\mathrm{GL}(n)$ 的闭子群。所有经典矩阵群都是由多项式方程定义的（如 $\det A = 1$，$A^T A = I$ 等），因此自动是闭子群，从而是 Lie 群。

---

## 55.2 一般线性群 $\mathrm{GL}(n)$

<div class="context-flow" markdown>

**核心问题**：$\mathrm{GL}(n)$ 作为所有经典矩阵群的"母群"具有什么基本结构？

</div>

!!! definition "定义 55.2 (一般线性群)"
    **实一般线性群**
    $$\mathrm{GL}(n, \mathbb{R}) = \{A \in M_n(\mathbb{R}) : \det A \neq 0\}.$$
    **复一般线性群**
    $$\mathrm{GL}(n, \mathbb{C}) = \{A \in M_n(\mathbb{C}) : \det A \neq 0\}.$$

!!! theorem "定理 55.2 ($\mathrm{GL}(n)$ 的基本性质)"
    **(a)** $\mathrm{GL}(n, \mathbb{R})$ 是 $M_n(\mathbb{R}) \cong \mathbb{R}^{n^2}$ 中的开集（因为 $\det$ 是连续函数，$\mathrm{GL}(n) = \det^{-1}(\mathbb{R} \setminus \{0\})$）。因此 $\dim \mathrm{GL}(n, \mathbb{R}) = n^2$。

    **(b)** $\mathrm{GL}(n, \mathbb{R})$ 有恰好两个连通分支：
    $$\mathrm{GL}^+(n, \mathbb{R}) = \{A : \det A > 0\}, \quad \mathrm{GL}^-(n, \mathbb{R}) = \{A : \det A < 0\}.$$

    **(c)** $\mathrm{GL}(n, \mathbb{C})$ 是连通的（且是 $M_n(\mathbb{C}) \cong \mathbb{C}^{n^2}$ 中的开集），维数为 $2n^2$（作为实流形）。

    **(d)** $\mathrm{GL}(n)$ 的 Lie 代数是整个矩阵空间：$\mathfrak{gl}(n) = M_n$。Lie 括号为 $[X, Y] = XY - YX$。

??? proof "证明"
    **(b)** $\det: \mathrm{GL}(n, \mathbb{R}) \to \mathbb{R}^*$ 是连续满射。$\mathbb{R}^* = \mathbb{R}_{>0} \sqcup \mathbb{R}_{<0}$ 有两个连通分支，因此 $\mathrm{GL}(n)$ 至少有两个连通分支。要证明恰好两个，需要说明 $\mathrm{GL}^+(n)$ 连通。

    任取 $A \in \mathrm{GL}^+(n)$。利用极分解 $A = OP$（$O \in \mathrm{SO}(n)$，$P$ 正定对称），以及 $\mathrm{SO}(n)$ 连通和正定矩阵空间凸（从而连通），可以连续地将 $A$ 变形为 $I$。

    **(c)** 对复矩阵，$\det: \mathrm{GL}(n, \mathbb{C}) \to \mathbb{C}^*$。$\mathbb{C}^* = \mathbb{C} \setminus \{0\}$ 连通，但这不直接证明 $\mathrm{GL}(n, \mathbb{C})$ 连通。正确的论证：任取 $A \in \mathrm{GL}(n, \mathbb{C})$，其 Jordan 标准形 $J$ 的对角线元素非零，可以连续地将每个对角元素沿 $\mathbb{C}^*$ 中的路径变形到 $1$，将每个上三角元素变形到 $0$，从而将 $A$ 变形到 $I$。

    **(d)** $\mathrm{GL}(n)$ 是 $M_n$ 中的开集，因此其切空间（在任何点）等于 $M_n$ 本身。$\blacksquare$

!!! example "例 55.1"
    $\mathrm{GL}(1, \mathbb{R}) = \mathbb{R}^* = \mathbb{R} \setminus \{0\}$，有两个连通分支 $\mathbb{R}_{>0}$ 和 $\mathbb{R}_{<0}$。$\mathrm{GL}(1, \mathbb{C}) = \mathbb{C}^*$，连通。

---

## 55.3 特殊线性群 $\mathrm{SL}(n)$

<div class="context-flow" markdown>

**核心问题**：行列式为 1 的矩阵构成怎样的群？

</div>

!!! definition "定义 55.3 (特殊线性群)"
    $$\mathrm{SL}(n, \mathbb{F}) = \{A \in \mathrm{GL}(n, \mathbb{F}) : \det A = 1\}.$$

!!! theorem "定理 55.3 ($\mathrm{SL}(n)$ 的性质)"
    **(a)** $\mathrm{SL}(n, \mathbb{R})$ 和 $\mathrm{SL}(n, \mathbb{C})$ 都是连通的 Lie 群。

    **(b)** $\dim \mathrm{SL}(n) = n^2 - 1$。

    **(c)** $\mathrm{SL}(n)$ 的 Lie 代数是
    $$\mathfrak{sl}(n) = \{X \in M_n : \operatorname{tr} X = 0\}.$$

??? proof "证明"
    **(c)** 若 $e^{tX} \in \mathrm{SL}(n)$ 对所有 $t$，则 $\det(e^{tX}) = e^{t \operatorname{tr} X} = 1$，故 $\operatorname{tr} X = 0$。反之，$\operatorname{tr} X = 0$ 蕴含 $\det(e^{tX}) = 1$。

    **(b)** $\dim \mathfrak{sl}(n) = n^2 - 1$（迹零矩阵的空间比全矩阵空间少一维）。

    **(a) 连通性**：对 $\mathrm{SL}(n, \mathbb{R})$，任取 $A \in \mathrm{SL}(n)$，极分解 $A = OP$，$O \in \mathrm{SO}(n) \subset \mathrm{SL}(n)$，$P$ 正定且 $\det P = 1$。$P = e^S$（$S$ 对称，$\operatorname{tr} S = 0$）可以沿 $e^{tS}$ 连续变形到 $I$。$\mathrm{SO}(n)$ 连通，故 $\mathrm{SL}(n)$ 连通。$\blacksquare$

!!! example "例 55.2"
    $\mathrm{SL}(2, \mathbb{R})$ 是 3 维 Lie 群。其 Lie 代数 $\mathfrak{sl}(2, \mathbb{R})$ 有标准基
    $$e = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}, \quad f = \begin{bmatrix} 0 & 0 \\ 1 & 0 \end{bmatrix}, \quad h = \begin{bmatrix} 1 & 0 \\ 0 & -1 \end{bmatrix},$$
    满足 Lie 括号关系 $[h, e] = 2e$，$[h, f] = -2f$，$[e, f] = h$。这是半单 Lie 代数理论的最基本范例。

---

## 55.4 正交群 $\mathrm{O}(n)$ 与旋转群 $\mathrm{SO}(n)$

<div class="context-flow" markdown>

**核心问题**：保持欧氏内积的线性变换群具有什么结构？

</div>

!!! definition "定义 55.4 (正交群和特殊正交群)"
    $$\mathrm{O}(n) = \{Q \in M_n(\mathbb{R}) : Q^T Q = I\},$$
    $$\mathrm{SO}(n) = \{Q \in \mathrm{O}(n) : \det Q = 1\}.$$

!!! theorem "定理 55.4 ($\mathrm{O}(n)$ 和 $\mathrm{SO}(n)$ 的性质)"
    **(a)** $\mathrm{O}(n)$ 是**紧** Lie 群，维数 $\frac{n(n-1)}{2}$。

    **(b)** $\mathrm{O}(n)$ 有恰好两个连通分支：$\mathrm{SO}(n)$（$\det = +1$）和 $\det = -1$ 的分支。

    **(c)** $\mathrm{SO}(n)$ 是连通的。$\mathrm{SO}(1) = \{1\}$，$\mathrm{SO}(2) \cong S^1$（圆周），$\mathrm{SO}(3)$ 同胚于 $\mathbb{RP}^3$（实射影空间）。

    **(d)** Lie 代数为
    $$\mathfrak{o}(n) = \mathfrak{so}(n) = \{X \in M_n(\mathbb{R}) : X^T + X = 0\}$$
    即**反对称矩阵**的空间。$\dim \mathfrak{so}(n) = \frac{n(n-1)}{2}$。

??? proof "证明"
    **紧性**：$\mathrm{O}(n)$ 是有界的（$Q^TQ = I$ 意味着每列是单位向量，故 $\|Q\|_F = \sqrt{n}$），且是闭的（由连续方程 $Q^TQ = I$ 定义），因此紧。

    **Lie 代数**：若 $e^{tX} \in \mathrm{O}(n)$，则 $(e^{tX})^T e^{tX} = I$，即 $e^{tX^T} e^{tX} = I$。在 $t = 0$ 处对 $t$ 求导：$X^T + X = 0$，即 $X$ 反对称。反之，$X^T = -X$ 蕴含 $(e^{tX})^T = e^{tX^T} = e^{-tX} = (e^{tX})^{-1}$。

    **维数**：$n \times n$ 反对称矩阵由上三角部分的 $\frac{n(n-1)}{2}$ 个元素决定。

    **连通分支**：$\det: \mathrm{O}(n) \to \{+1, -1\}$ 是连续满射，故 $\mathrm{O}(n)$ 至少两个连通分支。$\mathrm{SO}(n)$ 连通可由归纳法证明：考虑 $\mathrm{SO}(n)$ 对 $S^{n-1}$ 的传递作用，稳定子群同构于 $\mathrm{SO}(n-1)$，得到纤维丛 $\mathrm{SO}(n-1) \to \mathrm{SO}(n) \to S^{n-1}$。$S^{n-1}$ 在 $n \geq 2$ 时连通，归纳基础 $\mathrm{SO}(2) \cong S^1$ 连通。由纤维丛的长正合序列和归纳假设知 $\mathrm{SO}(n)$ 连通。$\blacksquare$

!!! theorem "定理 55.5 (Cayley 变换)"
    设 $X$ 是 $n \times n$ 反对称实矩阵。则
    $$Q = (I - X)(I + X)^{-1}$$
    是正交矩阵且 $\det Q = +1$（即 $Q \in \mathrm{SO}(n)$），当且仅当 $-1$ 不是 $X$ 的特征值（此时 $I + X$ 可逆）。

    反之，若 $Q \in \mathrm{SO}(n)$ 且 $-1$ 不是 $Q$ 的特征值，则 $X = (I - Q)(I + Q)^{-1}$ 是反对称矩阵。

    Cayley 变换在 $\mathfrak{so}(n)$（去掉奇异集）和 $\mathrm{SO}(n)$（去掉 $-1$ 特征值的矩阵）之间建立了有理映射的双射。

??? proof "证明"
    设 $X^T = -X$，$Q = (I-X)(I+X)^{-1}$。

    $$Q^T = ((I+X)^{-1})^T(I-X)^T = ((I+X)^T)^{-1}(I-X)^T = (I-X)^{-1}(I+X).$$

    $$Q^T Q = (I-X)^{-1}(I+X)(I-X)(I+X)^{-1}.$$

    因为 $X$ 反对称，$(I+X)$ 和 $(I-X)$ 交换：$(I+X)(I-X) = I - X^2 = (I-X)(I+X)$。所以

    $$Q^T Q = (I-X)^{-1}(I-X)(I+X)(I+X)^{-1} = I. \checkmark$$

    对 $\det Q$：$\det Q = \det(I-X)/\det(I+X)$。$\det(I-X) = \det((I+X)^T) = \det(I+X)$（因为 $(I-X)^T = I+X$），故 $\det Q = 1$。$\blacksquare$

!!! example "例 55.3"
    $n = 3$：取 $X = \begin{bmatrix} 0 & -\theta_3 & \theta_2 \\ \theta_3 & 0 & -\theta_1 \\ -\theta_2 & \theta_1 & 0 \end{bmatrix}$（反对称矩阵由向量 $\theta = (\theta_1, \theta_2, \theta_3)$ 参数化）。

    Cayley 变换 $Q = (I - X)(I + X)^{-1}$ 给出 $\mathrm{SO}(3)$ 中的旋转矩阵。这在计算机图形学中作为 **Cayley-Rodrigues 参数化**使用。

---

## 55.5 酉群 $\mathrm{U}(n)$ 与特殊酉群 $\mathrm{SU}(n)$

<div class="context-flow" markdown>

**核心问题**：保持 Hermite 内积的线性变换群有什么结构？

</div>

!!! definition "定义 55.5 (酉群和特殊酉群)"
    $$\mathrm{U}(n) = \{U \in M_n(\mathbb{C}) : U^* U = I\},$$
    $$\mathrm{SU}(n) = \{U \in \mathrm{U}(n) : \det U = 1\}.$$

!!! theorem "定理 55.6 ($\mathrm{U}(n)$ 和 $\mathrm{SU}(n)$ 的性质)"
    **(a)** $\mathrm{U}(n)$ 是紧、连通 Lie 群。$\dim_\mathbb{R} \mathrm{U}(n) = n^2$。

    **(b)** $\mathrm{SU}(n)$ 是紧、连通、单连通 Lie 群。$\dim_\mathbb{R} \mathrm{SU}(n) = n^2 - 1$。

    **(c)** $\mathrm{U}(n)$ 的 Lie 代数是**反 Hermite 矩阵**空间：
    $$\mathfrak{u}(n) = \{X \in M_n(\mathbb{C}) : X^* + X = 0\}.$$

    **(d)** $\mathrm{SU}(n)$ 的 Lie 代数是**迹零反 Hermite 矩阵**空间：
    $$\mathfrak{su}(n) = \{X \in M_n(\mathbb{C}) : X^* + X = 0, \operatorname{tr} X = 0\}.$$

??? proof "证明"
    Lie 代数的推导与正交群完全类似：若 $e^{tX} \in \mathrm{U}(n)$，则 $(e^{tX})^* e^{tX} = I$，即 $e^{tX^*} e^{tX} = I$。在 $t = 0$ 求导得 $X^* + X = 0$。加上 $\det e^{tX} = e^{t \operatorname{tr} X} = 1$ 的条件给出 $\operatorname{tr} X = 0$（注意 $\operatorname{tr} X$ 已经是纯虚数，条件即 $\operatorname{tr} X = 0$）。

    **维数**：$\mathfrak{u}(n)$ 中的矩阵 $X = (x_{jk})$ 满足 $x_{jk} = -\bar{x}_{kj}$。对角元素为纯虚数（$n$ 个实参数），上三角元素 $x_{jk}$（$j < k$）有 $\frac{n(n-1)}{2}$ 个，每个贡献 2 个实参数。总维数 $n + n(n-1) = n^2$。对 $\mathfrak{su}(n)$，加上 $\operatorname{tr} X = 0$（一个实约束，因为迹已经纯虚），维数 $n^2 - 1$。

    **连通性**：$\mathrm{U}(n)$ 传递地作用在 $S^{2n-1}$ 上，稳定子群为 $\mathrm{U}(n-1)$，纤维丛 $\mathrm{U}(n-1) \to \mathrm{U}(n) \to S^{2n-1}$。归纳基础 $\mathrm{U}(1) \cong S^1$ 连通。类似地，$\mathrm{SU}(n)$ 连通且单连通（$\mathrm{SU}(1) = \{1\}$，$\mathrm{SU}(2) \cong S^3$）。$\blacksquare$

!!! example "例 55.4"
    $\mathrm{SU}(2)$ 是 3 维紧单连通 Lie 群，微分同胚于 $S^3$。其元素可以参数化为
    $$U = \begin{bmatrix} \alpha & -\bar{\beta} \\ \beta & \bar{\alpha} \end{bmatrix}, \quad |\alpha|^2 + |\beta|^2 = 1.$$

    其 Lie 代数 $\mathfrak{su}(2)$ 有标准基（Pauli 矩阵的 $i$ 倍）：
    $$e_1 = \begin{bmatrix} i & 0 \\ 0 & -i \end{bmatrix}, \quad e_2 = \begin{bmatrix} 0 & 1 \\ -1 & 0 \end{bmatrix}, \quad e_3 = \begin{bmatrix} 0 & i \\ i & 0 \end{bmatrix}.$$

    括号关系：$[e_1, e_2] = 2e_3$，$[e_2, e_3] = 2e_1$，$[e_3, e_1] = 2e_2$。这与 $\mathfrak{so}(3)$ 同构（$\mathrm{SU}(2)$ 是 $\mathrm{SO}(3)$ 的二重覆盖）。

!!! theorem "定理 55.7 ($\mathrm{SU}(2)$ 与 $\mathrm{SO}(3)$ 的关系)"
    存在满射 Lie 群同态 $\pi: \mathrm{SU}(2) \to \mathrm{SO}(3)$，核为 $\{\pm I\}$。因此 $\mathrm{SO}(3) \cong \mathrm{SU}(2)/\{\pm I\}$，$\mathrm{SU}(2)$ 是 $\mathrm{SO}(3)$ 的**万有覆叠群**。

??? proof "证明"
    考虑 $\mathrm{SU}(2)$ 在 $\mathfrak{su}(2)$（3 维实向量空间，带有 Killing 内积）上的伴随表示 $\mathrm{Ad}: \mathrm{SU}(2) \to \mathrm{GL}(\mathfrak{su}(2))$。伴随表示保持 Killing 内积，故 $\mathrm{Ad}(\mathrm{SU}(2)) \subseteq \mathrm{O}(\mathfrak{su}(2)) \cong \mathrm{O}(3)$。由 $\mathrm{SU}(2)$ 连通，像在 $\mathrm{SO}(3)$ 中。

    $\mathrm{Ad}$ 的核是 $\mathrm{SU}(2)$ 的中心 $Z(\mathrm{SU}(2)) = \{\pm I\}$。

    由维数比较（$\dim \mathrm{SU}(2) = 3 = \dim \mathrm{SO}(3)$），$\mathrm{Ad}$ 是满射。$\blacksquare$

---

## 55.6 辛群 $\mathrm{Sp}(2n)$

<div class="context-flow" markdown>

**核心问题**：从 Lie 代数的视角如何理解辛群？

</div>

!!! definition "定义 55.6 (辛群与辛 Lie 代数)"
    辛群
    $$\mathrm{Sp}(2n, \mathbb{R}) = \{M \in \mathrm{GL}(2n, \mathbb{R}) : M^T J M = J\}$$
    已在第 53 章详细讨论。其 Lie 代数为
    $$\mathfrak{sp}(2n) = \{X \in M_{2n}(\mathbb{R}) : X^T J + J X = 0\}.$$
    等价地，$X \in \mathfrak{sp}(2n)$ 当且仅当 $JX$ 是对称矩阵。

!!! theorem "定理 55.8 ($\mathfrak{sp}(2n)$ 的结构)"
    $\mathfrak{sp}(2n)$ 的元素可以写成分块形式
    $$X = \begin{bmatrix} A & B \\ C & -A^T \end{bmatrix},$$
    其中 $B = B^T$，$C = C^T$ 是对称矩阵，$A$ 任意。

    $\dim \mathfrak{sp}(2n) = n^2 + 2 \cdot \frac{n(n+1)}{2} = n^2 + n(n+1) = n(2n+1)$。

!!! remark "注记"
    注意术语的混淆：在某些文献中，"$\mathrm{Sp}(n)$"指的是**紧辛群**（即保持四元数 Hermite 内积的 $n \times n$ 四元数矩阵群），它是 $\mathrm{Sp}(2n, \mathbb{C}) \cap \mathrm{U}(2n)$，是紧 Lie 群。本章的 $\mathrm{Sp}(2n, \mathbb{R})$ 是非紧的。

---

## 55.7 Lie 代数与切空间

<div class="context-flow" markdown>

**核心问题**：Lie 群在单位元的切空间为什么自然具有 Lie 代数结构？

</div>

!!! definition "定义 55.7 (矩阵 Lie 群的 Lie 代数)"
    设 $G$ 是矩阵 Lie 群（$\mathrm{GL}(n)$ 的闭子群）。$G$ 的 **Lie 代数** $\mathfrak{g}$ 定义为
    $$\mathfrak{g} = T_I G = \{X \in M_n : e^{tX} \in G, \forall t \in \mathbb{R}\},$$
    即 $G$ 在单位元 $I$ 处的切空间。

!!! theorem "定理 55.9 (Lie 代数的基本性质)"
    **(a)** $\mathfrak{g}$ 是 $M_n$ 的实线性子空间。

    **(b)** $\mathfrak{g}$ 对 **Lie 括号** $[X, Y] = XY - YX$ 封闭。

    **(c)** Lie 括号满足：
    - 双线性性：$[\alpha X + \beta Y, Z] = \alpha[X,Z] + \beta[Y,Z]$。
    - 反对称性：$[X, Y] = -[Y, X]$。
    - **Jacobi 恒等式**：$[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$。

??? proof "证明"
    **(a)** 和 **(b)** 已在定理 55.1 的证明中说明。

    **(c)** 双线性性和反对称性由 $[X,Y] = XY - YX$ 直接验证。

    Jacobi 恒等式：展开左边
    $$[X,[Y,Z]] = X(YZ-ZY) - (YZ-ZY)X = XYZ - XZY - YZX + ZYX,$$
    $$[Y,[Z,X]] = YZX - YXZ - ZXY + XZY,$$
    $$[Z,[X,Y]] = ZXY - ZYX - XYZ + YXZ.$$
    三式相加，所有 12 项两两抵消，和为 $0$。$\blacksquare$

!!! definition "定义 55.8 (Lie 代数同态)"
    设 $\mathfrak{g}, \mathfrak{h}$ 是 Lie 代数。线性映射 $\varphi: \mathfrak{g} \to \mathfrak{h}$ 称为 **Lie 代数同态**，若 $\varphi([X,Y]) = [\varphi(X), \varphi(Y)]$ 对所有 $X, Y \in \mathfrak{g}$。

!!! theorem "定理 55.10 (Lie 群同态诱导 Lie 代数同态)"
    设 $\Phi: G \to H$ 是矩阵 Lie 群之间的光滑群同态。则其在单位元的微分
    $$\phi = d\Phi_I: \mathfrak{g} \to \mathfrak{h}$$
    是 Lie 代数同态，且 $\Phi(e^X) = e^{\phi(X)}$ 对所有 $X \in \mathfrak{g}$。

??? proof "证明"
    由 $\Phi(e^{tX}) \in H$ 对所有 $t$，在 $t = 0$ 求导得 $\phi(X) = \frac{d}{dt}\Big|_{t=0} \Phi(e^{tX}) \in \mathfrak{h}$。

    $\Phi(e^{tX})$ 是 $H$ 中过 $I$ 的一参数子群，导数为 $\phi(X)$，故 $\Phi(e^{tX}) = e^{t\phi(X)}$（一参数子群的唯一性）。

    对 Lie 括号的保持：利用 $\frac{d^2}{ds\,dt}\Big|_{s=t=0} e^{sX}e^{tY}e^{-sX} = [X,Y]$ 和 $\Phi$ 的群同态性质推导。$\blacksquare$

下表汇总经典矩阵群及其 Lie 代数。

| 矩阵群 $G$ | 定义条件 | Lie 代数 $\mathfrak{g}$ | $\dim_\mathbb{R}$ |
|:---|:---|:---|:---|
| $\mathrm{GL}(n, \mathbb{R})$ | $\det A \neq 0$ | $\mathfrak{gl}(n) = M_n(\mathbb{R})$ | $n^2$ |
| $\mathrm{SL}(n, \mathbb{R})$ | $\det A = 1$ | $\mathfrak{sl}(n) = \{X : \operatorname{tr} X = 0\}$ | $n^2 - 1$ |
| $\mathrm{O}(n)$ | $Q^T Q = I$ | $\mathfrak{so}(n) = \{X : X^T = -X\}$ | $\frac{n(n-1)}{2}$ |
| $\mathrm{SO}(n)$ | $Q^T Q = I, \det Q = 1$ | $\mathfrak{so}(n) = \{X : X^T = -X\}$ | $\frac{n(n-1)}{2}$ |
| $\mathrm{U}(n)$ | $U^* U = I$ | $\mathfrak{u}(n) = \{X : X^* = -X\}$ | $n^2$ |
| $\mathrm{SU}(n)$ | $U^* U = I, \det U = 1$ | $\mathfrak{su}(n) = \{X : X^* = -X, \operatorname{tr} X = 0\}$ | $n^2 - 1$ |
| $\mathrm{Sp}(2n, \mathbb{R})$ | $M^T J M = J$ | $\mathfrak{sp}(2n) = \{X : X^T J + JX = 0\}$ | $n(2n+1)$ |

---

## 55.8 指数映射

<div class="context-flow" markdown>

**核心问题**：矩阵指数如何连接 Lie 代数与 Lie 群？它何时是满射？

</div>

!!! definition "定义 55.9 (指数映射)"
    矩阵 Lie 群 $G$ 的**指数映射** $\exp: \mathfrak{g} \to G$ 定义为
    $$\exp(X) = e^X = \sum_{k=0}^\infty \frac{X^k}{k!}.$$

!!! theorem "定理 55.11 (指数映射的基本性质)"
    **(a)** $\exp(0) = I$。

    **(b)** $\exp((s+t)X) = \exp(sX)\exp(tX)$ 对所有 $s, t \in \mathbb{R}$。即 $t \mapsto \exp(tX)$ 是 $G$ 中的**一参数子群**。

    **(c)** $\frac{d}{dt}\Big|_{t=0} \exp(tX) = X$。因此 $d(\exp)_0 = \mathrm{id}_\mathfrak{g}$，即 $\exp$ 在原点的微分是恒等映射。

    **(d)** 由反函数定理，$\exp$ 将 $\mathfrak{g}$ 中 $0$ 的某个邻域微分同胚到 $G$ 中 $I$ 的某个邻域。

    **(e)** 若 $XY = YX$，则 $e^{X+Y} = e^X e^Y$。

!!! theorem "定理 55.12 (紧连通群上的满射性)"
    若 $G$ 是**紧连通** Lie 群，则 $\exp: \mathfrak{g} \to G$ 是**满射**。

??? proof "证明"
    $G$ 紧连通意味着 $G$ 上存在双不变 Riemannian 度量。在此度量下，$\exp$ 恰好是 Riemannian 指数映射。由 Hopf-Rinow 定理，紧连通 Riemannian 流形上任意两点可以用测地线连接，因此 $\exp$ 满射。

    更初等的证明（对矩阵群）：$G$ 紧连通，每个 $g \in G$ 可对角化（在酉情形下），特征值在单位圆上，可以取对数。$\blacksquare$

!!! example "例 55.5 (指数映射不满射的例子)"
    $\mathrm{SL}(2, \mathbb{R})$ 是非紧连通 Lie 群。矩阵
    $$A = \begin{bmatrix} -1 & 1 \\ 0 & -1 \end{bmatrix} \in \mathrm{SL}(2, \mathbb{R})$$
    不在 $\exp(\mathfrak{sl}(2, \mathbb{R}))$ 的像中。

    **证明**：假设 $A = e^X$，$X \in \mathfrak{sl}(2)$，$\operatorname{tr} X = 0$。$A$ 的特征值为 $-1$（二重）。$e^X$ 的特征值为 $e^\lambda$（$\lambda$ 为 $X$ 的特征值）。$\operatorname{tr} X = 0$ 意味着 $X$ 的特征值为 $\lambda, -\lambda$。$e^\lambda = -1$ 且 $e^{-\lambda} = -1$ 要求 $\lambda = (2k+1)\pi i$ 且 $-\lambda = (2m+1)\pi i$，故 $\lambda = -(2m+1)\pi i$，加上第一个条件得 $(2k+1) = -(2m+1)$，即 $k + m = -1$。取 $k = 0, m = -1$：$\lambda = \pi i$。则 $X$ 的特征值为 $\pm \pi i$，$X$ 可对角化为 $\begin{bmatrix} \pi i & 0 \\ 0 & -\pi i \end{bmatrix}$，但 $X$ 必须是实矩阵。实矩阵特征值为 $\pm \pi i$ 时，$e^X$ 的特征值为 $e^{\pm \pi i} = -1$，且 $e^X$ 实际上是 $-I$ 的旋转版本，可以验证 $e^X = -I$（因为共轭特征值对应 $2 \times 2$ 旋转块 $R(\pi) = -I$）。但 $A \neq -I$（因为 $A$ 有非零的 $(1,2)$ 元素）。这导致矛盾。

!!! example "例 55.6 ($\mathrm{SO}(3)$ 的指数映射)"
    $\mathfrak{so}(3)$ 中的反对称矩阵 $\hat{\omega}$ 可由向量 $\omega = (\omega_1, \omega_2, \omega_3)^T$ 参数化：
    $$\hat{\omega} = \begin{bmatrix} 0 & -\omega_3 & \omega_2 \\ \omega_3 & 0 & -\omega_1 \\ -\omega_2 & \omega_1 & 0 \end{bmatrix}.$$

    **Rodrigues 公式**：$\theta = \|\omega\|$，$\hat{n} = \hat{\omega}/\theta$，则
    $$e^{\hat{\omega}} = I + \frac{\sin\theta}{\theta}\hat{\omega} + \frac{1 - \cos\theta}{\theta^2}\hat{\omega}^2.$$
    这给出绕轴 $\omega/\|\omega\|$ 旋转角度 $\theta = \|\omega\|$ 的旋转矩阵。

    $\exp: \mathfrak{so}(3) \to \mathrm{SO}(3)$ 是满射（因为 $\mathrm{SO}(3)$ 紧连通），但不是单射：$\hat{\omega}$ 和 $\hat{\omega} + 2\pi \hat{n}$ 给出相同的旋转。

---

## 55.9 伴随表示

<div class="context-flow" markdown>

**核心问题**：Lie 群如何通过共轭作用在自身的 Lie 代数上？

</div>

!!! definition "定义 55.10 (伴随表示 $\mathrm{Ad}$)"
    对矩阵 Lie 群 $G$，**伴随表示** $\mathrm{Ad}: G \to \mathrm{GL}(\mathfrak{g})$ 定义为
    $$\mathrm{Ad}_g(X) = gXg^{-1}, \quad g \in G, X \in \mathfrak{g}.$$
    这是一个群同态（从 $G$ 到 $\mathfrak{g}$ 上的可逆线性变换群）。

!!! definition "定义 55.11 (小伴随表示 $\mathrm{ad}$)"
    $\mathrm{Ad}$ 的微分定义了**小伴随表示** $\mathrm{ad}: \mathfrak{g} \to \mathfrak{gl}(\mathfrak{g})$：
    $$\mathrm{ad}_X(Y) = [X, Y], \quad X, Y \in \mathfrak{g}.$$

!!! theorem "定理 55.13 (伴随表示的性质)"
    **(a)** $\mathrm{Ad}$ 是 Lie 群同态：$\mathrm{Ad}_{gh} = \mathrm{Ad}_g \circ \mathrm{Ad}_h$。

    **(b)** $\mathrm{ad}$ 是 Lie 代数同态：$\mathrm{ad}_{[X,Y]} = [\mathrm{ad}_X, \mathrm{ad}_Y] = \mathrm{ad}_X \circ \mathrm{ad}_Y - \mathrm{ad}_Y \circ \mathrm{ad}_X$。

    **(c)** $\mathrm{Ad}_{e^X} = e^{\mathrm{ad}_X}$，即
    $$e^X Y e^{-X} = Y + [X,Y] + \frac{1}{2}[X,[X,Y]] + \frac{1}{3!}[X,[X,[X,Y]]] + \cdots = \sum_{k=0}^\infty \frac{1}{k!}(\mathrm{ad}_X)^k(Y).$$

??? proof "证明"
    **(a)** $\mathrm{Ad}_{gh}(X) = (gh)X(gh)^{-1} = g(hXh^{-1})g^{-1} = \mathrm{Ad}_g(\mathrm{Ad}_h(X))$。

    **(b)** 这等价于 Jacobi 恒等式：$\mathrm{ad}_{[X,Y]}(Z) = [[X,Y],Z]$，而 $[\mathrm{ad}_X, \mathrm{ad}_Y](Z) = [X,[Y,Z]] - [Y,[X,Z]]$。Jacobi 恒等式 $[[X,Y],Z] = [X,[Y,Z]] - [Y,[X,Z]]$ 保证两者相等。

    **(c)** 定义 $f(t) = e^{tX} Y e^{-tX}$。则 $f(0) = Y$，$f'(t) = e^{tX}(XY - YX)e^{-tX} = e^{tX}[X,Y]e^{-tX} = \mathrm{Ad}_{e^{tX}}([X,Y])$。更一般地，$f^{(k)}(0) = (\mathrm{ad}_X)^k(Y)$。Taylor 展开给出结论。$\blacksquare$

!!! definition "定义 55.12 (Killing 型)"
    Lie 代数 $\mathfrak{g}$ 上的 **Killing 型**（Killing form）是双线性型
    $$B(X, Y) = \operatorname{tr}(\mathrm{ad}_X \circ \mathrm{ad}_Y), \quad X, Y \in \mathfrak{g}.$$

!!! theorem "定理 55.14 (Killing 型的性质)"
    **(a)** $B$ 是对称双线性型：$B(X, Y) = B(Y, X)$。

    **(b)** $B$ 是 $\mathrm{Ad}$-不变的：$B(\mathrm{Ad}_g X, \mathrm{Ad}_g Y) = B(X, Y)$。

    **(c)** $B$ 是 $\mathrm{ad}$-不变的：$B([Z,X], Y) + B(X, [Z,Y]) = 0$。

    **(d)** (Cartan 判据) $\mathfrak{g}$ 是**半单**的（即没有非零可交换理想）当且仅当 Killing 型非退化。

??? proof "证明"
    **(a)** $B(X,Y) = \operatorname{tr}(\mathrm{ad}_X \mathrm{ad}_Y) = \operatorname{tr}(\mathrm{ad}_Y \mathrm{ad}_X) = B(Y,X)$（迹的轮换性）。

    **(c)** $B([Z,X],Y) = \operatorname{tr}(\mathrm{ad}_{[Z,X]} \mathrm{ad}_Y) = \operatorname{tr}([\mathrm{ad}_Z, \mathrm{ad}_X] \mathrm{ad}_Y) = \operatorname{tr}(\mathrm{ad}_Z \mathrm{ad}_X \mathrm{ad}_Y - \mathrm{ad}_X \mathrm{ad}_Z \mathrm{ad}_Y)$。
    $B(X,[Z,Y]) = \operatorname{tr}(\mathrm{ad}_X \mathrm{ad}_{[Z,Y]}) = \operatorname{tr}(\mathrm{ad}_X (\mathrm{ad}_Z \mathrm{ad}_Y - \mathrm{ad}_Y \mathrm{ad}_Z))$。
    两者相加，利用 $\operatorname{tr}(ABC) = \operatorname{tr}(CAB)$，各项抵消。$\blacksquare$

!!! example "例 55.7"
    对 $\mathfrak{sl}(n, \mathbb{C})$，Killing 型为 $B(X,Y) = 2n \operatorname{tr}(XY)$。它是非退化的，反映了 $\mathfrak{sl}(n)$ 是半单 Lie 代数。

    对 $\mathfrak{so}(n)$，Killing 型为 $B(X,Y) = (n-2)\operatorname{tr}(XY)$（$n \geq 3$ 时非退化）。

---

## 55.10 Baker-Campbell-Hausdorff 公式

<div class="context-flow" markdown>

**核心问题**：$e^X e^Y$ 何时等于 $e^Z$？$Z$ 如何用 $X, Y$ 及其 Lie 括号表示？

</div>

!!! theorem "定理 55.15 (Baker-Campbell-Hausdorff 公式)"
    设 $X, Y \in \mathfrak{g}$，$\|X\|$ 和 $\|Y\|$ 足够小。则存在 $Z \in \mathfrak{g}$ 使得 $e^X e^Y = e^Z$，且
    $$Z = \log(e^X e^Y) = X + Y + \frac{1}{2}[X,Y] + \frac{1}{12}\big([X,[X,Y]] + [Y,[Y,X]]\big) + \cdots$$
    其中高阶项完全由 $X, Y$ 的嵌套 Lie 括号决定。

??? proof "证明"
    **第一步（存在性）**：由 $\exp$ 在 $0$ 附近是微分同胚，当 $X, Y$ 足够小时，$e^X e^Y$ 落在 $\exp$ 的像的邻域中，故可以定义 $Z = \log(e^X e^Y)$。

    **第二步（前几项计算）**：利用 $e^X = I + X + \frac{X^2}{2} + \cdots$ 和 $\log(I + W) = W - \frac{W^2}{2} + \cdots$：

    $$e^X e^Y = (I + X + \frac{X^2}{2} + \cdots)(I + Y + \frac{Y^2}{2} + \cdots)$$
    $$= I + (X+Y) + (XY + \frac{X^2}{2} + \frac{Y^2}{2}) + \cdots$$

    令 $W = e^X e^Y - I = (X+Y) + (XY + \frac{X^2}{2} + \frac{Y^2}{2}) + \cdots$

    $$Z = \log(I+W) = W - \frac{W^2}{2} + \frac{W^3}{3} - \cdots$$

    **一阶项**：$Z_1 = X + Y$。

    **二阶项**：$W$ 的二阶部分为 $XY + \frac{X^2}{2} + \frac{Y^2}{2}$。$W^2$ 的二阶部分为 $(X+Y)^2 = X^2 + XY + YX + Y^2$。所以
    $$Z_2 = XY + \frac{X^2}{2} + \frac{Y^2}{2} - \frac{1}{2}(X^2 + XY + YX + Y^2) = \frac{1}{2}(XY - YX) = \frac{1}{2}[X,Y].$$

    **三阶项**：经过类似但更繁琐的计算，
    $$Z_3 = \frac{1}{12}[X,[X,Y]] + \frac{1}{12}[Y,[Y,X]] = \frac{1}{12}[X,[X,Y]] - \frac{1}{12}[Y,[X,Y]].$$

    **关键定性结论**（Dynkin 1947）：$Z$ 的每一项都可以表示为 $X$ 和 $Y$ 的嵌套 Lie 括号的有限线性组合。这意味着 $Z \in \mathfrak{g}$（因为 $\mathfrak{g}$ 对 Lie 括号封闭）。$\blacksquare$

!!! theorem "定理 55.16 (BCH 公式的收敛性)"
    BCH 级数在 $\|X\| + \|Y\| < \log 2$ 时绝对收敛（此处 $\|\cdot\|$ 为某个与 Lie 群相容的矩阵范数）。

!!! example "例 55.8"
    若 $[X, Y] = 0$（$X, Y$ 交换），则 BCH 公式的所有高阶项消失，$Z = X + Y$，即 $e^X e^Y = e^{X+Y}$。

!!! example "例 55.9"
    若 $[X, [X, Y]] = 0$ 且 $[Y, [X, Y]] = 0$（$[X,Y]$ 与 $X, Y$ 都交换），则
    $$e^X e^Y = e^{X + Y + \frac{1}{2}[X,Y]}.$$
    这一特殊情形在量子力学中经常出现（例如 Weyl 关系 $e^{i\alpha \hat{q}} e^{i\beta \hat{p}} = e^{-i\alpha\beta\hbar/2} e^{i(\alpha\hat{q} + \beta\hat{p})}$，其中 $[\hat{q}, \hat{p}] = i\hbar$ 是标量，自动与一切交换）。

!!! theorem "定理 55.17 (Zassenhaus 公式)"
    $e^{X+Y}$ 也可以分解为无穷乘积：
    $$e^{X+Y} = e^X e^Y e^{-\frac{1}{2}[X,Y]} e^{\frac{1}{3}[Y,[X,Y]] + \frac{1}{6}[X,[X,Y]]} \cdots$$
    每个因子都是嵌套括号的指数。这在数值方法（算子分裂法）中有重要应用。

!!! remark "注记"
    BCH 公式的核心意义是：**Lie 群的群乘法完全由 Lie 代数的括号运算决定**（至少在单位元附近）。这是为什么研究 Lie 群可以"下降"到研究 Lie 代数——一个纯线性代数对象——的根本原因。

---

## 55.11 $\mathfrak{sl}(2)$ 的不可约表示分类

<div class="context-flow" markdown>

**核心问题**：Lie 代数 $\mathfrak{sl}(2, \mathbb{C})$ 有哪些不可约表示？如何显式构造？

</div>

$\mathfrak{sl}(2, \mathbb{C})$ 的表示论是整个半单 Lie 代数表示论的基石。所有半单 Lie 代数的有限维表示理论最终都归结为对若干 $\mathfrak{sl}(2)$ 子代数的分析。

!!! definition "定义 55.13 (权空间与权向量)"
    设 $(\rho, V)$ 是 $\mathfrak{sl}(2, \mathbb{C})$ 的有限维表示，即 $\rho: \mathfrak{sl}(2) \to \mathfrak{gl}(V)$ 是 Lie 代数同态。利用例 55.2 中的标准基 $\{e, f, h\}$，定义 $V$ 的**权空间**（weight space）：

    $$V_\lambda = \{v \in V : \rho(h) v = \lambda v\}, \quad \lambda \in \mathbb{C}.$$

    非零的 $V_\lambda$ 中的向量称为权 $\lambda$ 的**权向量**。标量 $\lambda$ 称为**权**（weight）。

!!! theorem "定理 55.18 (权空间的基本性质)"
    设 $(\rho, V)$ 是 $\mathfrak{sl}(2, \mathbb{C})$ 的有限维表示。

    **(a)** $V$ 是权空间的直和：$V = \bigoplus_{\lambda} V_\lambda$（因为 $\rho(h)$ 在有限维空间上可对角化）。

    **(b)** $\rho(e)$ 将 $V_\lambda$ 映入 $V_{\lambda+2}$：若 $v \in V_\lambda$，则 $\rho(e)v \in V_{\lambda+2}$。

    **(c)** $\rho(f)$ 将 $V_\lambda$ 映入 $V_{\lambda-2}$：若 $v \in V_\lambda$，则 $\rho(f)v \in V_{\lambda-2}$。

    因此 $\rho(e)$ 是"升权算子"，$\rho(f)$ 是"降权算子"。

??? proof "证明"
    **(a)** 由 $\rho(h)$ 在有限维复向量空间上的 Jordan 分解。利用 $[h, e] = 2e$ 和 $[h, f] = -2f$，可以证明 $\rho(h)$ 的所有广义特征值的广义特征空间就是特征空间（即 $\rho(h)$ 半单/可对角化）。

    **(b)** 设 $v \in V_\lambda$。利用 $[h, e] = 2e$：
    $$\rho(h)\rho(e)v = \rho(e)\rho(h)v + \rho([h,e])v = \rho(e)(\lambda v) + 2\rho(e)v = (\lambda + 2)\rho(e)v.$$
    故 $\rho(e)v \in V_{\lambda+2}$。**(c)** 的证明完全类似。$\blacksquare$

!!! definition "定义 55.14 (最高权向量)"
    表示 $(\rho, V)$ 中的非零向量 $v_0 \in V_\lambda$ 称为**最高权向量**（highest weight vector），如果 $\rho(e)v_0 = 0$。对应的权 $\lambda$ 称为**最高权**。

!!! theorem "定理 55.19 ($\mathfrak{sl}(2)$ 不可约表示的分类)"
    **(a)** 对每个非负整数 $n \in \mathbb{Z}_{\geq 0}$，存在唯一（至同构）的 $\mathfrak{sl}(2, \mathbb{C})$ 的 $(n+1)$ 维不可约表示 $V_n$。

    **(b)** $V_n$ 有基 $\{v_0, v_1, \ldots, v_n\}$，$\mathfrak{sl}(2)$ 的作用为：

    $$\rho(h) v_k = (n - 2k) v_k, \quad k = 0, 1, \ldots, n,$$

    $$\rho(e) v_k = k(n - k + 1) v_{k-1}, \quad k = 1, \ldots, n, \quad \rho(e) v_0 = 0,$$

    $$\rho(f) v_k = v_{k+1}, \quad k = 0, \ldots, n-1, \quad \rho(f) v_n = 0.$$

    **(c)** $V_n$ 的权为 $n, n-2, n-4, \ldots, -n+2, -n$（每个权的重数为 1）。最高权为 $n$。

    **(d)** 每个 $\mathfrak{sl}(2)$ 的有限维不可约表示都同构于某个 $V_n$。

??? proof "证明"
    **存在性（显式构造）。** 设 $V_n$ 是以 $\{v_0, v_1, \ldots, v_n\}$ 为基的 $(n+1)$ 维向量空间。定义 $\rho(h), \rho(e), \rho(f)$ 的作用如定理所述。需要验证 Lie 括号关系被保持：

    - $[\rho(h), \rho(e)] = 2\rho(e)$：$\rho(h)\rho(e)v_k - \rho(e)\rho(h)v_k = (n-2(k-1))k(n-k+1)v_{k-1} - (n-2k)k(n-k+1)v_{k-1} = 2k(n-k+1)v_{k-1} = 2\rho(e)v_k$。
    - $[\rho(h), \rho(f)] = -2\rho(f)$：类似验证。
    - $[\rho(e), \rho(f)] = \rho(h)$：$\rho(e)\rho(f)v_k - \rho(f)\rho(e)v_k = (k+1)(n-k)v_k - k(n-k+1)v_k = (n-2k)v_k = \rho(h)v_k$。

    **唯一性。** 设 $W$ 是最高权为 $n$ 的不可约表示，$w_0$ 为最高权向量。令 $w_k = \rho(f)^k w_0 / k!$（适当归一化）。由不可约性，$\{w_0, w_1, \ldots\}$ 生成 $W$。利用权的约束，$\rho(f)^{n+1}w_0 = 0$（否则权无下界，矛盾于有限维），故 $\dim W = n + 1$。括号关系唯一确定了 $\rho(e), \rho(f), \rho(h)$ 在此基上的作用，从而 $W \cong V_n$。

    **不可约性。** 设 $W \subset V_n$ 是非零子表示。$W$ 是权空间的直和，取 $W$ 中最高权的权向量 $w$。反复施加 $\rho(f)$ 得到所有基向量（因为 $\rho(f)v_k = v_{k+1}$ 非零直到 $k = n$），故 $W = V_n$。$\blacksquare$

!!! example "例 55.10"
    **低维不可约表示：**

    - $V_0$：1 维平凡表示。$\rho(e) = \rho(f) = \rho(h) = 0$。
    - $V_1$：2 维标准表示（$\mathfrak{sl}(2)$ 在 $\mathbb{C}^2$ 上的自然作用）。基 $v_0, v_1$，权为 $1, -1$。
    - $V_2$：3 维伴随表示（$\mathfrak{sl}(2)$ 对自身的伴随作用 $\mathrm{ad}$）。基 $v_0, v_1, v_2$，权为 $2, 0, -2$。

    也可以将 $V_n$ 实现为 $\mathbb{C}^2$ 上 $n$ 次齐次多项式空间 $\mathbb{C}[x, y]_n$：
    $$\rho(e) = x\frac{\partial}{\partial y}, \quad \rho(f) = y\frac{\partial}{\partial x}, \quad \rho(h) = x\frac{\partial}{\partial x} - y\frac{\partial}{\partial y}.$$
    基向量 $v_k = \frac{x^{n-k}y^k}{k!}$，$\dim \mathbb{C}[x,y]_n = n + 1$。

!!! theorem "定理 55.20 ($\mathfrak{sl}(2)$ 表示的完全可约性)"
    $\mathfrak{sl}(2, \mathbb{C})$ 的每个有限维表示都是不可约表示的直和：

    $$V \cong V_{n_1} \oplus V_{n_2} \oplus \cdots \oplus V_{n_k}.$$

    分解由 $V$ 的权的重数唯一确定。

??? proof "证明"
    这是 Weyl 完全可约性定理（定理 55.23）的特殊情形。但对 $\mathfrak{sl}(2)$ 有更初等的证明。

    **Casimir 算子法**：定义 $C = \rho(h)^2 + 2\rho(e)\rho(f) + 2\rho(f)\rho(e) \in \mathrm{End}(V)$。可以验证 $C$ 与 $\rho(\mathfrak{sl}(2))$ 中的所有算子交换（$C$ 是泛包络代数的中心元素）。由 Schur 引理，$C$ 在每个不可约分量 $V_n$ 上是标量 $n(n+2)$。

    若 $V$ 不是完全可约的，则存在不可约子表示 $W \subset V$，使得 $V/W$ 有不可约子表示 $V_m$，但 $V$ 中对应的子空间不是 $W$ 的直和补。然而 $C$ 在 $W$ 上为标量 $n(n+2)$，在提升到 $V$ 中的 $V_m$ 部分为标量 $m(m+2)$。当 $n \neq m$ 时，$C - n(n+2)I$ 的核给出直和分裂。当 $n = m$ 时，利用 $\rho(e), \rho(f)$ 的具体公式可以直接证明扩张分裂。$\blacksquare$

---

## 55.12 根系初步

<div class="context-flow" markdown>

**核心问题**：如何将 $\mathfrak{sl}(2)$ 的权空间分解推广到一般半单 Lie 代数？

</div>

!!! definition "定义 55.15 (Cartan 子代数)"
    半单 Lie 代数 $\mathfrak{g}$ 的**Cartan 子代数** $\mathfrak{h}$ 是极大交换半单子代数。对矩阵 Lie 代数，$\mathfrak{h}$ 由可同时对角化的元素组成。

    例如，$\mathfrak{sl}(n, \mathbb{C})$ 的 Cartan 子代数为对角迹零矩阵空间，维数 $\ell = n - 1$。$\ell = \dim \mathfrak{h}$ 称为 $\mathfrak{g}$ 的**秩**（rank）。

!!! theorem "定理 55.21 (根分解)"
    设 $\mathfrak{g}$ 是复半单 Lie 代数，$\mathfrak{h}$ 为 Cartan 子代数。则 $\mathfrak{g}$ 有**根分解**（root decomposition）：

    $$\mathfrak{g} = \mathfrak{h} \oplus \bigoplus_{\alpha \in \Phi} \mathfrak{g}_\alpha,$$

    其中 $\Phi \subset \mathfrak{h}^* \setminus \{0\}$ 是**根系**（root system），$\mathfrak{g}_\alpha = \{X \in \mathfrak{g} : [H, X] = \alpha(H) X, \forall H \in \mathfrak{h}\}$ 是**根空间**。

    基本性质：
    - 每个根空间 $\mathfrak{g}_\alpha$ 都是一维的。
    - 若 $\alpha \in \Phi$，则 $-\alpha \in \Phi$。
    - $[\mathfrak{g}_\alpha, \mathfrak{g}_\beta] \subseteq \mathfrak{g}_{\alpha+\beta}$（若 $\alpha + \beta \in \Phi$），否则为 $\{0\}$。
    - 对每个 $\alpha \in \Phi$，$\mathfrak{g}_\alpha \oplus \mathfrak{g}_{-\alpha} \oplus [\mathfrak{g}_\alpha, \mathfrak{g}_{-\alpha}]$ 构成 $\mathfrak{sl}(2)$ 的一个副本。

!!! example "例 55.11"
    **经典 Lie 代数的根系：**

    - **$A_{n-1}$ 型**（$\mathfrak{sl}(n)$）：根为 $\alpha_{ij} = \varepsilon_i - \varepsilon_j$（$i \neq j$），其中 $\varepsilon_i$ 是 $\mathfrak{h}^*$ 的标准基。共有 $n(n-1)$ 个根，秩 $n - 1$。根空间 $\mathfrak{g}_{\varepsilon_i - \varepsilon_j}$ 由矩阵单位 $E_{ij}$ 生成。

    - **$B_n$ 型**（$\mathfrak{so}(2n+1)$）：根为 $\pm \varepsilon_i \pm \varepsilon_j$（$i \neq j$）和 $\pm \varepsilon_i$。共有 $2n^2$ 个根。

    - **$C_n$ 型**（$\mathfrak{sp}(2n)$）：根为 $\pm \varepsilon_i \pm \varepsilon_j$（$i \neq j$）和 $\pm 2\varepsilon_i$。共有 $2n^2$ 个根。

    - **$D_n$ 型**（$\mathfrak{so}(2n)$）：根为 $\pm \varepsilon_i \pm \varepsilon_j$（$i \neq j$）。共有 $2n(n-1)$ 个根。

!!! definition "定义 55.16 (Weyl 群)"
    根系 $\Phi$ 的 **Weyl 群** $W$ 是由根反射 $s_\alpha$（$\alpha \in \Phi$）生成的 $\mathfrak{h}^*$ 上的有限群：

    $$s_\alpha(\lambda) = \lambda - \frac{2\langle \lambda, \alpha \rangle}{\langle \alpha, \alpha \rangle}\alpha, \quad \lambda \in \mathfrak{h}^*.$$

    $W$ 是有限 Coxeter 群。例如，$A_{n-1}$ 的 Weyl 群是对称群 $S_n$（置换 $\varepsilon_1, \ldots, \varepsilon_n$），$B_n$ 和 $C_n$ 的 Weyl 群是超八面体群 $S_n \ltimes (\mathbb{Z}/2)^n$。

!!! definition "定义 55.17 (Dynkin 图)"
    选取根系 $\Phi$ 的一组**单根**（simple roots）$\Delta = \{\alpha_1, \ldots, \alpha_\ell\}$（$\ell$ = 秩），则 **Dynkin 图**是以单根为顶点的图，顶点 $\alpha_i$ 和 $\alpha_j$ 之间的连线数由 Cartan 矩阵 $A_{ij} = \frac{2\langle \alpha_i, \alpha_j \rangle}{\langle \alpha_j, \alpha_j \rangle}$ 决定。

!!! theorem "定理 55.22 (Killing-Cartan 分类)"
    不可约根系（从而单 Lie 代数）的 Dynkin 图恰好有以下类型：

    - **经典型**：$A_n$（$n \geq 1$）、$B_n$（$n \geq 2$）、$C_n$（$n \geq 3$）、$D_n$（$n \geq 4$），分别对应 $\mathfrak{sl}(n+1)$、$\mathfrak{so}(2n+1)$、$\mathfrak{sp}(2n)$、$\mathfrak{so}(2n)$。

    - **例外型**：$G_2$、$F_4$、$E_6$、$E_7$、$E_8$。

    没有其他的。这是数学中最深刻、最优美的分类定理之一。

---

## 55.13 Weyl 完全可约性定理

<div class="context-flow" markdown>

**核心问题**：半单 Lie 代数的有限维表示是否总能分解为不可约表示的直和？

</div>

!!! theorem "定理 55.23 (Weyl 完全可约性定理)"
    设 $\mathfrak{g}$ 是**半单** Lie 代数（即 Killing 型非退化）。则 $\mathfrak{g}$ 的每个有限维表示都是**完全可约**的，即可以分解为不可约子表示的直和。

    等价表述：若 $W \subset V$ 是子表示，则存在子表示 $W'$ 使得 $V = W \oplus W'$。

??? proof "证明（概要）"
    **第一步（Casimir 算子）。** 利用 $\mathfrak{g}$ 上的 Killing 型 $B$ 非退化，取 $\mathfrak{g}$ 的基 $\{X_i\}$ 和对偶基 $\{Y_i\}$（即 $B(X_i, Y_j) = \delta_{ij}$）。定义 **Casimir 算子**

    $$C_\rho = \sum_i \rho(X_i)\rho(Y_i) \in \mathrm{End}(V).$$

    由 Killing 型的 $\mathrm{ad}$-不变性，$C_\rho$ 与 $\rho(\mathfrak{g})$ 中的所有算子交换。由 Schur 引理，$C_\rho$ 在每个不可约子表示上是标量。

    **第二步（一维扩张分裂）。** 设 $0 \to W \to V \to \mathbb{C} \to 0$ 是短正合列（$\mathbb{C}$ 是平凡表示）。$C_\rho$ 在 $\mathbb{C}$ 上为 $0$，在 $W$ 的不可约分量上为非零标量（对非平凡不可约表示可以计算 $C_\rho$ 的值非零）。因此 $C_\rho$ 的核给出分裂。

    **第三步（一般情形归约）。** 设 $W \subset V$ 是子表示。考虑 $\mathrm{Hom}(V/W, V)$ 上的 $\mathfrak{g}$-表示。利用第二步的结论（对 $\mathrm{Hom}$ 空间中的一维子表示），找到 $\mathfrak{g}$-不变的投影 $V \to W$，其核就是 $W$ 的直和补。$\blacksquare$

!!! remark "注记"
    Weyl 完全可约性定理对**可解**或**幂零** Lie 代数一般不成立。例如，上三角矩阵 Lie 代数 $\mathfrak{b} \subset \mathfrak{sl}(2)$ 在 $\mathbb{C}^2$ 上的自然表示不完全可约：子空间 $\mathbb{C} e_1$ 是子表示，但没有 $\mathfrak{b}$-不变的补。

---

## 55.14 紧 Lie 群上的调和分析

<div class="context-flow" markdown>

**核心问题**：紧 Lie 群上的函数如何用表示论来分析？

</div>

!!! theorem "定理 55.24 (Peter-Weyl 定理)"
    设 $G$ 是紧 Lie 群，$\hat{G}$ 是 $G$ 的所有（有限维）不可约酉表示的等价类集合。则 $L^2(G)$（关于 Haar 测度的平方可积函数空间）有正交分解：

    $$L^2(G) = \widehat{\bigoplus_{\pi \in \hat{G}}} \, d_\pi \cdot V_\pi,$$

    其中 $d_\pi = \dim \pi$，求和是 Hilbert 空间的直和。具体地，$\pi$ 的矩阵元素 $\pi_{ij}(g)$（$1 \leq i, j \leq d_\pi$）组成 $L^2(G)$ 的完备正交系（适当归一化后）。

    这是 Fourier 分析从 $S^1 \cong U(1)$ 到一般紧群的推广：$U(1)$ 的不可约表示是 $e^{in\theta}$（$n \in \mathbb{Z}$），Peter-Weyl 定理退化为经典 Fourier 级数。

!!! theorem "定理 55.25 (极大环面定理)"
    设 $G$ 是紧连通 Lie 群。

    **(a)** $G$ 包含**极大环面** $T$（即极大连通交换子群，同构于 $(S^1)^\ell$），$\ell$ 称为 $G$ 的秩。

    **(b)** 每个 $g \in G$ 都共轭于 $T$ 中的某个元素：$G = \bigcup_{x \in G} xTx^{-1}$。

    **(c)** $G$ 的不可约表示由其在 $T$ 上的限制（即**最高权**）唯一确定。

    例如，$U(n)$ 的极大环面是对角酉矩阵 $T = \{diag(e^{i\theta_1}, \ldots, e^{i\theta_n})\}$，秩为 $n$。

---

## 55.15 例外 Lie 群

<div class="context-flow" markdown>

**核心问题**：除了经典矩阵群，还有哪些单 Lie 群？

</div>

Killing-Cartan 分类（定理 55.22）表明，除了经典系列 $A_n, B_n, C_n, D_n$ 之外，还有五个例外单 Lie 代数：

!!! definition "定义 55.18 (例外 Lie 群)"

    | 例外型 | 秩 | $\dim \mathfrak{g}$ | 简要描述 |
    |:------|:---|:-------------------|:---------|
    | $G_2$ | 2 | 14 | 八元数自同构群 $\mathrm{Aut}(\mathbb{O})$；与 7 维叉积相关 |
    | $F_4$ | 4 | 52 | 八元数射影平面 $\mathbb{OP}^2$ 的等距群 |
    | $E_6$ | 6 | 78 | 出现在弦论的 Calabi-Yau 紧化中 |
    | $E_7$ | 7 | 133 | 与超引力理论中的对称性相关 |
    | $E_8$ | 8 | 248 | 最大的例外 Lie 群；其根格是 8 维中最密球堆积 |

    例外 Lie 群虽然不能用经典矩阵群（正交、酉、辛）直接描述，但都可以实现为某种矩阵群（由 Cartan 定理保证）。$E_8$ 的最小忠实表示是 248 维的（恰好等于 $\dim \mathfrak{e}_8$，即伴随表示本身就是最小表示）。

!!! remark "注记"
    例外 Lie 群在纯数学（代数几何、数论）和理论物理（弦论、M 理论）中都有深刻的应用。$E_8 \times E_8$ 是杂化弦理论的规范群。$G_2$ 流形是 M 理论紧化的重要几何结构。

---

## 本章小结

| 矩阵群 | 定义条件 | Lie 代数 | 维数 | 连通/紧 |
|:---|:---|:---|:---|:---|
| $\mathrm{GL}(n, \mathbb{R})$ | $\det \neq 0$ | $M_n(\mathbb{R})$ | $n^2$ | 不连通/不紧 |
| $\mathrm{SL}(n, \mathbb{R})$ | $\det = 1$ | $\operatorname{tr} = 0$ | $n^2-1$ | 连通/不紧 |
| $\mathrm{O}(n)$ | $Q^TQ=I$ | 反对称 | $\frac{n(n-1)}{2}$ | 不连通/紧 |
| $\mathrm{SO}(n)$ | $Q^TQ=I, \det=1$ | 反对称 | $\frac{n(n-1)}{2}$ | 连通/紧 |
| $\mathrm{U}(n)$ | $U^*U=I$ | 反 Hermite | $n^2$ | 连通/紧 |
| $\mathrm{SU}(n)$ | $U^*U=I, \det=1$ | 反 Hermite + $\operatorname{tr}=0$ | $n^2-1$ | 连通/紧/单连通 |
| $\mathrm{Sp}(2n)$ | $M^TJM=J$ | $X^TJ+JX=0$ | $n(2n+1)$ | 连通/不紧 |

关键概念：

- **Lie 代数** $\mathfrak{g}$：切空间 + Lie 括号 $[X,Y] = XY - YX$。
- **指数映射** $\exp: \mathfrak{g} \to G$：在原点局部微分同胚；紧连通群上满射。
- **伴随表示**：$\mathrm{Ad}_g(X) = gXg^{-1}$，$\mathrm{ad}_X(Y) = [X,Y]$。
- **BCH 公式**：$\log(e^X e^Y) = X + Y + \frac{1}{2}[X,Y] + \cdots$（纯括号表达式）。
- **$\mathfrak{sl}(2)$ 表示论**：不可约表示由最高权 $n \in \mathbb{Z}_{\geq 0}$ 分类，维度 $n + 1$。
- **根系与 Dynkin 图**：半单 Lie 代数的结构编码在有限组合数据中。
- **Weyl 定理**：半单 Lie 代数的有限维表示完全可约。
- **Peter-Weyl 定理**：紧群上的 $L^2$ 函数按不可约表示分解。

---

## 习题

!!! exercise "习题 55.1"
    证明 $\mathrm{GL}(n, \mathbb{C})$ 连通。（提示：利用 Jordan 标准形将任意可逆矩阵连续变形到 $I$。）

!!! exercise "习题 55.2"
    计算 $\mathfrak{sl}(2, \mathbb{R})$ 的 Killing 型。验证它是非退化的，从而 $\mathfrak{sl}(2)$ 半单。在基 $\{e, f, h\}$（例 55.2）下写出 Killing 型的矩阵。

!!! exercise "习题 55.3"
    设 $X \in \mathfrak{so}(3)$，$\|X\| = \theta$。利用 $X^3 = -\theta^2 X$，推导 Rodrigues 公式
    $$e^X = I + \frac{\sin\theta}{\theta}X + \frac{1-\cos\theta}{\theta^2}X^2.$$

!!! exercise "习题 55.4"
    证明 Cayley 变换的性质：若 $X^T = -X$（反对称），则 $Q = (I-X)(I+X)^{-1}$ 满足 $Q^TQ = I$ 且 $\det Q = 1$。

!!! exercise "习题 55.5"
    证明 $\mathrm{U}(n)$ 上的指数映射是满射。（提示：酉矩阵可以对角化，对角元素在单位圆上，可以取对数。）

!!! exercise "习题 55.6"
    对 $X = \begin{bmatrix} 0 & -1 \\ 1 & 0 \end{bmatrix}$，$Y = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}$，计算 BCH 公式的前三阶项 $Z = X + Y + \frac{1}{2}[X,Y] + \cdots$，并与 $\log(e^X e^Y)$ 的直接计算比较。

!!! exercise "习题 55.7"
    证明：$\mathrm{O}(n)$ 的两个连通分支 $\mathrm{SO}(n)$ 和 $\{Q \in \mathrm{O}(n) : \det Q = -1\}$ 作为流形微分同胚（但只有前者是子群）。

!!! exercise "习题 55.8"
    设 $G$ 是紧连通矩阵 Lie 群。证明 $\mathfrak{g}$ 上的 Killing 型 $B$ 是半负定的（即 $B(X,X) \leq 0$）。（提示：$G$ 紧意味着 $\mathrm{Ad}$ 保持某个内积，故 $\mathrm{ad}_X$ 是反对称的。）

!!! exercise "习题 55.9"
    验证标准粒子物理模型的规范群 $SU(3) \times SU(2) \times U(1)$ 的总维数为 $8 + 3 + 1 = 12$。

!!! exercise "习题 55.10"
    设 $SE(3)$ 是三维刚体运动群（旋转 + 平移），由 $4 \times 4$ 矩阵
    $$\begin{bmatrix} R & t \\ 0 & 1 \end{bmatrix}, \quad R \in \mathrm{SO}(3), t \in \mathbb{R}^3$$
    构成。

    (a) 验证 $SE(3)$ 是矩阵群。

    (b) 写出其 Lie 代数 $\mathfrak{se}(3)$ 的一般元素形式。

    (c) 计算 $\dim SE(3) = 6$。

!!! exercise "习题 55.11"
    证明 Jacobi 恒等式 $[X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0$ 对矩阵 Lie 括号 $[A,B] = AB - BA$ 成立。
