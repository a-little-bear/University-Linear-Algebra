# 第 31 章 Majorization 与双随机矩阵

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · 凸集与凸优化 (Ch25) · 矩阵不等式 (Ch18)

**本章脉络**：优序 (Majorization) 的几何定义 $\to$ Hardy-Littlewood-Pólya 定理 $\to$ 双随机矩阵定义与性质 $\to$ Birkhoff 定理（置换矩阵的凸包） $\to$ Schur-Horn 定理（谱与对角元的纽带） $\to$ 算子单调与 Schur 凸性 $\to$ Robin 不等式 $\to$ 熵的 Majorization 性质 $\to$ 算子级数中的应用

**延伸**：Majorization 是量子信息论中纯态转换、量子信道容量以及统计学中不平等度量（Gini 系数）的统一数学框架；它量化了“分布的均匀程度”

</div>

在线性代数中，我们经常需要比较两个向量的“分散程度”或“混乱程度”。**优序理论**（Majorization）为此提供了一个强有力的数学框架。它不仅连接了矩阵的特征值与对角元素（Schur-Horn 定理），还通过**双随机矩阵**（Doubly Stochastic Matrices）将这种序关系与凸几何联系起来。本章将揭示这种隐藏在矩阵数值分布背后的深刻规律。

---

## 31.1 优序 (Majorization) 的定义

!!! definition "定义 31.1 (向量优序)"
    设 $x, y \in \mathbb{R}^n$，将其分量按非递增顺序排列为 $x_{[1]} \ge x_{[2]} \ge \cdots \ge x_{[n]}$。
    称 $y$ **优于** $x$（或 $x$ 被 $y$ 优序），记作 $x \prec y$，如果：
    1.  对 $k=1, \ldots, n-1$，有 $\sum_{i=1}^k x_{[i]} \le \sum_{i=1}^k y_{[i]}$。
    2.  总和相等：$\sum_{i=1}^n x_i = \sum_{i=1}^n y_i$。

!!! intuition "直观理解"
    $x \prec y$ 意味着 $x$ 比 $y$ “更均匀”或“更不集中”。例如，$(1/n, \ldots, 1/n) \prec x \prec (1, 0, \ldots, 0)$ 对任何和为 1 的非负向量 $x$ 都成立。

---

## 31.2 双随机矩阵与 Birkhoff 定理

!!! definition "定义 31.2 (双随机矩阵)"
    矩阵 $P \in M_n(\mathbb{R})$ 称为**双随机矩阵**，如果其元素非负且每行、每列之和均等于 1。

!!! theorem "定理 31.1 (Birkhoff-von Neumann 定理)"
    $n$ 阶双随机矩阵的全集 $\Omega_n$ 是 $n$ 阶**置换矩阵**（Permutation Matrices）的凸包。
    这意味着任何双随机矩阵都可以表示为置换矩阵的凸线性组合：$P = \sum \alpha_i P_{\sigma_i}$。

!!! theorem "定理 31.2 (Hardy-Littlewood-Pólya)"
    对于 $x, y \in \mathbb{R}^n$， $x \prec y$ 当且仅当存在双随机矩阵 $P$ 使得 $x = Py$。

---

## 31.3 Schur-Horn 定理

!!! theorem "定理 31.3 (Schur-Horn 定理)"
    设 $A$ 是 $n$ 阶 Hermite 矩阵，$\mathbf{d} = (a_{11}, \ldots, a_{nn})$ 为其对角元向量，$\boldsymbol{\lambda} = (\lambda_1, \ldots, \lambda_n)$ 为其特征值向量。则：
    $$\mathbf{d} \prec \boldsymbol{\lambda}$$
    反之，若 $\mathbf{d} \prec \boldsymbol{\lambda}$，则必存在以 $\mathbf{d}$ 为对角元、$\boldsymbol{\lambda}$ 为特征值的 Hermite 矩阵。

---

## 31.4 Schur 凸性

!!! definition "定义 31.3 (Schur 凸函数)"
    函数 $\phi: \mathbb{R}^n \to \mathbb{R}$ 称为 **Schur 凸**的，如果 $x \prec y \Rightarrow \phi(x) \le \phi(y)$。
    **例子**：$\phi(x) = \sum x_i^2$ 和 $\phi(x) = \sum x_i \log x_i$ 都是 Schur 凸的。

---

## 练习题

**1. [基础] 判定 $(2, 2, 2) \prec (3, 2, 1)$ 是否成立。**

??? success "参考答案"
    **验证步骤：**
    1. **排序**：$x = (2, 2, 2)$ 已降序，$y = (3, 2, 1)$ 已降序。
    2. **检查总和**：$\sum x_i = 2+2+2 = 6$；$\sum y_i = 3+2+1 = 6$。相等。
    3. **检查部分和**：
       - $k=1$：$x_1 = 2 \le y_1 = 3$（满足）。
       - $k=2$：$x_1+x_2 = 4 \le y_1+y_2 = 5$（满足）。
    **结论**：由于所有条件均满足， $(2, 2, 2) \prec (3, 2, 1)$ 成立。这说明均匀分布被非均匀分布优序。

**2. [Birkhoff] 将双随机矩阵 $\begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$ 分解为置换矩阵的组合。**

??? success "参考答案"
    **构造：**
    1. 寻找 2 阶置换矩阵：$P_1 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$（单位阵），$P_2 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。
    2. 尝试线性组合：$0.5 P_1 + 0.5 P_2 = \begin{pmatrix} 0.5 & 0 \\ 0 & 0.5 \end{pmatrix} + \begin{pmatrix} 0 & 0.5 \\ 0.5 & 0 \end{pmatrix} = \begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$。
    **结论**：成功分解。这验证了 Birkhoff 定理，即双随机矩阵位于置换矩阵构成的凸包内。

**3. [谱性质] 证明双随机矩阵的最大特征值为 1。**

??? success "参考答案"
    **证明：**
    1. **存在性**：令 $\mathbf{1} = (1, \ldots, 1)^T$。由于每行之和为 1，$P \mathbf{1} = \mathbf{1}$，说明 1 是特征值。
    2. **模最大性**：双随机矩阵是非负矩阵且其无穷范数 $\|P\|_\infty = \max (\text{行和}) = 1$。
    3. 根据谱半径性质：$\rho(P) \le \|P\|_\infty = 1$。
    **结论**：谱半径恰好为 1。由 Perron-Frobenius 定理，该特征值对应于平稳概率分布。

**4. [Schur-Horn] 若对称阵对角元为 $(5, 5)$，其特征值可能为 $(10, 0)$ 吗？**

??? success "参考答案"
    **判定：**
    1. 根据 Schur-Horn 定理，对角元向量必须被特征值向量优序。
    2. 检查：$(5, 5) \prec (10, 0)$？
       - 总和：$5+5=10 = 10+0$（相等）。
       - 部分和：$5 \le 10$（成立）。
    **结论**：可以。存在这样的矩阵，例如 $A = \begin{pmatrix} 5 & 5 \\ 5 & 5 \end{pmatrix}$，其特征值正是 10 和 0。

**5. [熵] 证明香农熵 $H(p) = -\sum p_i \log p_i$ 是 Schur 凹函数。**

??? success "参考答案"
    **推导：**
    1. 考虑函数 $f(x) = x \log x$。其二阶导数 $f''(x) = 1/x > 0$（当 $x > 0$），故 $f(x)$ 是凸函数。
    2. 根据性质，若 $f$ 是凸的，则 $\phi(p) = \sum f(p_i)$ 是 Schur 凸函数。
    3. 既然 $\sum p_i \log p_i$ 是 Schur 凸的，那么其相反数 $-\sum p_i \log p_i$ 必然是 **Schur 凹**的。
    **意义**：这在物理上解释了为什么概率分布越均匀（优序越小），系统熵越大。

**6. [计算] 计算 $x = \begin{pmatrix} 1/3 & 2/3 \\ 2/3 & 1/3 \end{pmatrix} \begin{pmatrix} 3 \\ 0 \end{pmatrix}$，并验证 $x \prec (3, 0)^T$。**

??? success "参考答案"
    **计算：**
    $x = \begin{pmatrix} (1/3) \cdot 3 + (2/3) \cdot 0 \\ (2/3) \cdot 3 + (1/3) \cdot 0 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}$。
    **验证：**
    1. 排序后向量为 $(2, 1)$。
    2. 检查：$(2, 1) \prec (3, 0)$？
       - $2+1=3 = 3+0$。
       - $2 \le 3$。
    **结论**：成立。这是 Hardy-Littlewood-Pólya 定理的直接应用：双随机矩阵作用于向量会使结果“更均匀”。

**7. [极值] 在所有和为 1 的非负向量中，哪一个在优序关系中最小？**

??? success "参考答案"
    **结论：**
    均匀分布向量 $\mathbf{u} = (1/n, 1/n, \ldots, 1/n)^T$。
    **理由**：对于任何满足和为 1 的降序向量 $x$，由于平均值是 $1/n$，其前 $k$ 项之和 $\sum_{i=1}^k x_i$ 必然大于等于 $k \cdot (1/n)$。这意味着 $\mathbf{u} \prec x$ 对任何 $x$ 都成立。

**8. [迹] 证明：若 $A \ge 0$ 是双随机的，则 $\operatorname{tr}(A) \le n$。**

??? success "参考答案"
    **证明：**
    1. 双随机矩阵的每个元素 $a_{ij}$ 必须满足 $0 \le a_{ij} \le 1$（因为行/列和为 1 且非负）。
    2. 特别地，对角元 $a_{ii} \le 1$。
    3. 迹是所有对角元之和：$\operatorname{tr}(A) = \sum_{i=1}^n a_{ii} \le \sum_{i=1}^n 1 = n$。
    **补充**：仅当 $A$ 是置换矩阵且不含错排时取等号（即 $A=I$）。

**9. [凸性] 判定 $\phi(x) = \max(x_i)$ 是否为 Schur 凸函数。**

??? success "参考答案"
    **判定：**
    1. 若 $x \prec y$，则 $y$ 的分布比 $x$ 更极端。
    2. $x \prec y$ 的定义直接包含了 $x_{[1]} \le y_{[1]}$。
    3. 而 $\phi(x) = x_{[1]}$。
    **结论**：是的。最大值算子是 Schur 凸的，它衡量了分布的最顶端集中程度。

**10. [应用] 在经济学中，Majorization 如何描述财富分配？**

??? success "参考答案"
    **解释：**
    1. 财富向量 $x \prec y$ 意味着 $x$ 的财富分配比 $y$ 更平均。
    2. 洛伦兹曲线（Lorenz Curve）的上方包容关系正是优序关系的几何表现。
    3. 如果收入分布从 $y$ 变为 $x$，说明基尼系数（一个 Schur 凸函数）下降，社会不平等程度降低。
