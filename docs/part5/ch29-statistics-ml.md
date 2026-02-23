# 第 29 章 线性代数在统计与机器学习中的应用

<div class="context-flow" markdown>

**前置**：特征值/SVD(Ch6,11) · 正定矩阵(Ch16) · 最小二乘(Ch25) · 矩阵分解(Ch10)

**脉络**：协方差矩阵(对称半正定) → PCA(特征分解/SVD) → 线性回归(正规方程/QR) → LDA(广义特征值) → 核方法(Gram矩阵) → SVM(二次规划) → 神经网络(矩阵链乘/梯度) → 推荐系统(低秩矩阵分解)

</div>

统计学与机器学习是线性代数最重要的应用领域之一。从多元统计的协方差分析到深度学习的反向传播，线性代数提供了统一的数学语言和高效的计算工具。本章将系统阐述线性代数在统计与机器学习中的核心作用，涵盖多元统计基础、主成分分析、线性回归与正则化、线性判别分析、核方法、支持向量机、神经网络以及推荐系统中的矩阵分解。

---

## 29.1 多元统计基础与协方差矩阵

<div class="context-flow" markdown>

**线性代数视角**：$n$ 个样本、$p$ 个特征 → 数据矩阵 $X \in \mathbb{R}^{n \times p}$ → 协方差矩阵 $\Sigma = \frac{1}{n-1}X_c^T X_c$（半正定, Ch16）→ 特征值 = 各方向方差 → Mahalanobis 距离 = $\Sigma^{-1}$ 加权内积

**链接**：Ch8 内积空间 · Ch16 正定矩阵

</div>

多元统计的基础是将数据组织为矩阵，并通过协方差矩阵刻画变量间的相关结构。

!!! definition "定义 29.1 (数据矩阵与中心化)"
    设有 $n$ 个样本，每个样本有 $p$ 个特征。**数据矩阵**（data matrix）$X \in \mathbb{R}^{n \times p}$ 的第 $i$ 行 $\mathbf{x}_i^T$ 为第 $i$ 个样本的特征向量。样本均值为

    $$
    \bar{\mathbf{x}} = \frac{1}{n}\sum_{i=1}^n \mathbf{x}_i = \frac{1}{n}X^T \mathbf{1}_n.
    $$

    **中心化数据矩阵**（centered data matrix）定义为

    $$
    X_c = X - \mathbf{1}_n \bar{\mathbf{x}}^T = \left(I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T\right)X = H X,
    $$

    其中 $H = I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T$ 为**中心化矩阵**（centering matrix），它是一个正交投影矩阵（$H^2 = H$，$H^T = H$）。

!!! definition "定义 29.2 (样本协方差矩阵)"
    **样本协方差矩阵**（sample covariance matrix）定义为

    $$
    S = \frac{1}{n-1}X_c^T X_c \in \mathbb{R}^{p \times p}.
    $$

    $S$ 是对称半正定矩阵。当 $n > p$ 且数据不在任何超平面上时，$S$ 为正定矩阵。$S$ 的第 $(i,j)$ 元素为第 $i$ 和第 $j$ 个特征的样本协方差：

    $$
    S_{ij} = \frac{1}{n-1}\sum_{k=1}^n (x_{ki} - \bar{x}_i)(x_{kj} - \bar{x}_j).
    $$

!!! theorem "定理 29.1 (协方差矩阵的谱性质)"
    设 $S \in \mathbb{R}^{p \times p}$ 为样本协方差矩阵，其特征分解为

    $$
    S = Q \Lambda Q^T = \sum_{j=1}^p \lambda_j \mathbf{q}_j \mathbf{q}_j^T,
    $$

    其中 $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_p \ge 0$。则：

    1. 数据在方向 $\mathbf{q}_j$ 上的方差为 $\lambda_j$。
    2. **总方差**（total variance）$= \operatorname{tr}(S) = \sum_{j=1}^p \lambda_j$。
    3. **广义方差**（generalized variance）$= \det(S) = \prod_{j=1}^p \lambda_j$。
    4. $\operatorname{rank}(S) = \operatorname{rank}(X_c) \le \min(n-1, p)$。

??? proof "证明"
    对于方向 $\mathbf{v}$（$\|\mathbf{v}\| = 1$），数据在该方向上的投影为 $X_c \mathbf{v}$，其方差为

    $$
    \frac{1}{n-1}\|X_c \mathbf{v}\|^2 = \frac{1}{n-1}\mathbf{v}^T X_c^T X_c \mathbf{v} = \mathbf{v}^T S \mathbf{v}.
    $$

    当 $\mathbf{v} = \mathbf{q}_j$ 时，$\mathbf{q}_j^T S \mathbf{q}_j = \lambda_j$。总方差 $\sum_j \operatorname{Var}(X_c \mathbf{e}_j) = \sum_j S_{jj} = \operatorname{tr}(S) = \sum_j \lambda_j$。由于 $S = \frac{1}{n-1}X_c^T X_c$，$\operatorname{rank}(S) = \operatorname{rank}(X_c)$。中心化去掉一个自由度，故 $\operatorname{rank}(X_c) \le \min(n-1, p)$。$\blacksquare$

!!! definition "定义 29.3 (Mahalanobis 距离)"
    设 $S$ 为正定协方差矩阵。点 $\mathbf{x}$ 与 $\mathbf{y}$ 之间的 **Mahalanobis 距离**为

    $$
    d_M(\mathbf{x}, \mathbf{y}) = \sqrt{(\mathbf{x} - \mathbf{y})^T S^{-1} (\mathbf{x} - \mathbf{y})}.
    $$

    这等价于先对数据做白化变换 $\mathbf{z} = S^{-1/2}\mathbf{x}$，再计算欧氏距离。Mahalanobis 距离考虑了各特征的尺度和相关性。

!!! example "例 29.1"
    **协方差矩阵的计算与谱分析。** 设三个二维样本为 $\mathbf{x}_1 = (1,2)^T$，$\mathbf{x}_2 = (3,4)^T$，$\mathbf{x}_3 = (5,6)^T$。数据矩阵与均值：

    $$
    X = \begin{pmatrix}1&2\\3&4\\5&6\end{pmatrix}, \quad \bar{\mathbf{x}} = (3, 4)^T.
    $$

    中心化数据：

    $$
    X_c = \begin{pmatrix}-2&-2\\0&0\\2&2\end{pmatrix}, \quad S = \frac{1}{2}X_c^T X_c = \frac{1}{2}\begin{pmatrix}8&8\\8&8\end{pmatrix} = \begin{pmatrix}4&4\\4&4\end{pmatrix}.
    $$

    特征值 $\lambda_1 = 8$，$\lambda_2 = 0$，对应特征向量 $\mathbf{q}_1 = \frac{1}{\sqrt{2}}(1,1)^T$。数据完全分布在 $(1,1)^T$ 方向上，与 $\operatorname{rank}(S) = 1$ 一致。

!!! example "例 29.2"
    **Mahalanobis 距离的几何意义。** 设

    $$
    S = \begin{pmatrix}4&2\\2&3\end{pmatrix}, \quad S^{-1} = \frac{1}{8}\begin{pmatrix}3&-2\\-2&4\end{pmatrix}.
    $$

    点 $\mathbf{a} = (0,0)^T$ 和 $\mathbf{b} = (2,1)^T$ 的 Mahalanobis 距离为

    $$
    d_M^2 = (2,1)\frac{1}{8}\begin{pmatrix}3&-2\\-2&4\end{pmatrix}\begin{pmatrix}2\\1\end{pmatrix} = \frac{1}{8}(2,1)\begin{pmatrix}4\\0\end{pmatrix} = 1.
    $$

    而欧氏距离为 $\sqrt{5} \approx 2.24$。Mahalanobis 距离较小，因为 $(2,1)^T$ 方向是数据方差较大的方向。

---

## 29.2 主成分分析（PCA）

<div class="context-flow" markdown>

**核心**：PCA = 对协方差矩阵做特征分解 = 对中心化数据矩阵做 SVD → 前 $k$ 个主成分 = 方差最大的 $k$ 个正交方向 = 最佳秩-$k$ 近似

**链接**：Ch11 SVD · Ch6 特征值 · Eckart-Young 定理

</div>

主成分分析是多元统计中最基本的降维方法，其数学本质是特征值分解或奇异值分解。

!!! definition "定义 29.4 (主成分)"
    设中心化数据矩阵为 $X_c \in \mathbb{R}^{n \times p}$，样本协方差矩阵 $S = \frac{1}{n-1}X_c^T X_c$ 的特征分解为 $S = Q\Lambda Q^T$。第 $j$ 个**主成分**（principal component）为

    $$
    \mathbf{z}_j = X_c \mathbf{q}_j \in \mathbb{R}^n, \quad j = 1, \ldots, p,
    $$

    其中 $\mathbf{q}_j$ 为 $S$ 的第 $j$ 个特征向量（按特征值降序排列）。$\mathbf{q}_j$ 称为第 $j$ 个**主成分方向**（loading vector）。保留前 $k$ 个主成分得到降维表示 $Z_k = X_c Q_k \in \mathbb{R}^{n \times k}$。

!!! theorem "定理 29.2 (PCA 的最优性)"
    前 $k$ 个主成分方向 $\mathbf{q}_1, \ldots, \mathbf{q}_k$ 是如下优化问题的解：

    $$
    \max_{\substack{V \in \mathbb{R}^{p \times k} \\ V^T V = I_k}} \operatorname{tr}(V^T S V).
    $$

    等价地，在所有秩为 $k$ 的矩阵中，$\hat{X}_c = X_c Q_k Q_k^T$ 使得重构误差 $\|X_c - \hat{X}_c\|_F^2$ 最小。前 $k$ 个主成分解释的方差比例为

    $$
    \frac{\sum_{j=1}^k \lambda_j}{\sum_{j=1}^p \lambda_j}.
    $$

??? proof "证明"
    设 $X_c = U\Sigma Q^T$ 为 SVD，其中 $\sigma_j^2 = (n-1)\lambda_j$。对于任意满足 $V^T V = I_k$ 的 $V$，

    $$
    \operatorname{tr}(V^T S V) = \frac{1}{n-1}\operatorname{tr}(V^T X_c^T X_c V) = \frac{1}{n-1}\|X_c V\|_F^2.
    $$

    由 Eckart-Young 定理（Ch11），$X_c Q_k Q_k^T$ 是 $X_c$ 的最佳秩-$k$ 近似。投影到 $Q_k$ 的 $k$ 列张成的子空间上保留最大的 Frobenius 范数，等价于 $\operatorname{tr}(V^T S V)$ 最大化在 $V = Q_k$ 处取到。重构误差

    $$
    \|X_c - X_c Q_k Q_k^T\|_F^2 = \sum_{j=k+1}^p \sigma_j^2 = (n-1)\sum_{j=k+1}^p \lambda_j. \quad \blacksquare
    $$

!!! theorem "定理 29.3 (通过 SVD 计算 PCA)"
    设中心化数据矩阵 $X_c = U\Sigma V^T$ 为 SVD。则：

    1. 主成分方向为 $V$ 的列。
    2. 主成分得分为 $Z = U\Sigma$，即 $Z$ 的第 $j$ 列为 $\sigma_j \mathbf{u}_j$。
    3. 协方差矩阵的特征值 $\lambda_j = \sigma_j^2/(n-1)$。

    当 $n \gg p$ 时，直接对 $S$ 做特征分解（$O(p^3)$）更高效；当 $p \gg n$ 时，先计算 $X_c X_c^T$ 的特征分解（$O(n^3)$）或使用截断 SVD 更高效。

??? proof "证明"
    由 $X_c = U\Sigma V^T$，

    $$
    S = \frac{1}{n-1}X_c^T X_c = \frac{1}{n-1}V\Sigma^T U^T U\Sigma V^T = V\left(\frac{\Sigma^2}{n-1}\right)V^T.
    $$

    这正是 $S$ 的特征分解，特征值 $\lambda_j = \sigma_j^2/(n-1)$，特征向量为 $V$ 的列。主成分得分 $Z = X_c V = U\Sigma V^T V = U\Sigma$。$\blacksquare$

!!! example "例 29.3"
    **PCA 降维。** 设 $4$ 个三维数据点的中心化数据矩阵为

    $$
    X_c = \begin{pmatrix}2&1&0\\-1&0&1\\0&-1&-1\\-1&0&0\end{pmatrix}.
    $$

    协方差矩阵 $S = \frac{1}{3}X_c^T X_c = \frac{1}{3}\begin{pmatrix}6&2&-1\\2&2&-1\\-1&-1&2\end{pmatrix}$。设特征值（近似）为 $\lambda_1 \approx 2.54$，$\lambda_2 \approx 0.98$，$\lambda_3 \approx 0.15$。前两个主成分解释的方差比例约为 $(2.54 + 0.98)/3.67 \approx 96\%$，因此可以用 $k = 2$ 个主成分有效表示原始三维数据。

---

## 29.3 线性回归与正则化

<div class="context-flow" markdown>

**线性代数核心**：$\hat{\boldsymbol{\beta}} = (X^TX)^{-1}X^T\mathbf{y}$（正规方程）= 投影 $\mathbf{y}$ 到 $\operatorname{col}(X)$ → 岭回归 = 谱过滤（$\sigma_i^2 \to \sigma_i^2 + \lambda$） → LASSO = $\ell_1$ 稀疏

**链接**：Ch25 最小二乘 · Ch11 SVD · Ch7 正交投影

</div>

线性回归是最基本的统计建模工具，其求解和分析完全依赖线性代数。

!!! definition "定义 29.5 (线性回归模型)"
    **线性回归模型**（linear regression model）假设响应变量 $\mathbf{y} \in \mathbb{R}^n$ 与设计矩阵 $X \in \mathbb{R}^{n \times p}$ 之间的关系为

    $$
    \mathbf{y} = X\boldsymbol{\beta} + \boldsymbol{\varepsilon},
    $$

    其中 $\boldsymbol{\beta} \in \mathbb{R}^p$ 为回归系数，$\boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, \sigma^2 I_n)$ 为噪声。**最小二乘估计**（ordinary least squares, OLS）为

    $$
    \hat{\boldsymbol{\beta}}_{\text{OLS}} = \arg\min_{\boldsymbol{\beta}} \|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 = (X^T X)^{-1}X^T \mathbf{y}.
    $$

!!! theorem "定理 29.4 (OLS 的几何与谱解释)"
    设 $X = U\Sigma V^T$ 为 SVD，$\operatorname{rank}(X) = r$。则：

    1. **几何解释**：$\hat{\mathbf{y}} = X\hat{\boldsymbol{\beta}}_{\text{OLS}} = H\mathbf{y}$，其中 $H = X(X^TX)^{-1}X^T = U_r U_r^T$ 为到 $\operatorname{col}(X)$ 的正交投影矩阵（帽子矩阵）。
    2. **谱解释**：$\hat{\boldsymbol{\beta}}_{\text{OLS}} = \sum_{j=1}^r \frac{\mathbf{u}_j^T \mathbf{y}}{\sigma_j}\mathbf{v}_j$。
    3. **方差**：$\operatorname{Var}(\hat{\boldsymbol{\beta}}_{\text{OLS}}) = \sigma^2(X^TX)^{-1} = \sigma^2 V \Sigma^{-2} V^T$，故小奇异值 $\sigma_j$ 导致估计方差大。

??? proof "证明"
    正规方程 $X^TX\hat{\boldsymbol{\beta}} = X^T\mathbf{y}$。代入 SVD：$(V\Sigma^2 V^T)\hat{\boldsymbol{\beta}} = V\Sigma U^T\mathbf{y}$，故 $\hat{\boldsymbol{\beta}} = V\Sigma^{-1}U^T\mathbf{y} = \sum_j \frac{\mathbf{u}_j^T\mathbf{y}}{\sigma_j}\mathbf{v}_j$。帽子矩阵 $H = X(X^TX)^{-1}X^T = U\Sigma V^T \cdot V\Sigma^{-2}V^T \cdot V\Sigma U^T = UU^T$（取前 $r$ 列）。方差由 $\hat{\boldsymbol{\beta}} = (X^TX)^{-1}X^T\mathbf{y}$ 和 $\operatorname{Var}(\mathbf{y}) = \sigma^2 I$ 得 $\operatorname{Var}(\hat{\boldsymbol{\beta}}) = \sigma^2(X^TX)^{-1}X^T X (X^TX)^{-1} = \sigma^2(X^TX)^{-1}$。$\blacksquare$

!!! definition "定义 29.6 (岭回归与 LASSO)"
    **岭回归**（ridge regression）在 OLS 基础上加 $\ell_2$ 惩罚：

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = \arg\min_{\boldsymbol{\beta}} \|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 + \lambda\|\boldsymbol{\beta}\|_2^2 = (X^TX + \lambda I)^{-1}X^T\mathbf{y}.
    $$

    **LASSO**（Least Absolute Shrinkage and Selection Operator）加 $\ell_1$ 惩罚：

    $$
    \hat{\boldsymbol{\beta}}_{\text{LASSO}} = \arg\min_{\boldsymbol{\beta}} \frac{1}{2}\|\mathbf{y} - X\boldsymbol{\beta}\|_2^2 + \lambda\|\boldsymbol{\beta}\|_1.
    $$

    LASSO 倾向于产生稀疏解（部分 $\beta_j = 0$），因此同时实现特征选择。

!!! theorem "定理 29.5 (岭回归的谱过滤解释)"
    设 $X = U\Sigma V^T$，则岭回归的解可以写为

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = \sum_{j=1}^r \frac{\sigma_j^2}{\sigma_j^2 + \lambda} \cdot \frac{\mathbf{u}_j^T \mathbf{y}}{\sigma_j} \mathbf{v}_j.
    $$

    因子 $\frac{\sigma_j^2}{\sigma_j^2 + \lambda}$ 是**谱过滤器**：对大奇异值（$\sigma_j^2 \gg \lambda$）几乎不衰减，对小奇异值（$\sigma_j^2 \ll \lambda$）强烈压缩。这是处理多重共线性的关键机制。

??? proof "证明"
    $(X^TX + \lambda I)^{-1} = (V\Sigma^2 V^T + \lambda VV^T)^{-1} = V(\Sigma^2 + \lambda I)^{-1}V^T$（假设 $V$ 是方阵或相应地处理）。于是

    $$
    \hat{\boldsymbol{\beta}}_{\text{ridge}} = V(\Sigma^2 + \lambda I)^{-1}\Sigma U^T\mathbf{y} = \sum_{j=1}^r \frac{\sigma_j}{\sigma_j^2+\lambda}(\mathbf{u}_j^T\mathbf{y})\mathbf{v}_j = \sum_{j=1}^r \frac{\sigma_j^2}{\sigma_j^2+\lambda}\cdot\frac{\mathbf{u}_j^T\mathbf{y}}{\sigma_j}\mathbf{v}_j. \quad \blacksquare
    $$

!!! example "例 29.4"
    **多重共线性与岭回归。** 设设计矩阵

    $$
    X = \begin{pmatrix}1&1\\1&1.01\\1&0.99\end{pmatrix}, \quad \mathbf{y} = \begin{pmatrix}2\\2.5\\1.5\end{pmatrix}.
    $$

    $X$ 的两列几乎线性相关，$\sigma_1 \approx 2.45$，$\sigma_2 \approx 0.01$。OLS 解对噪声极度敏感（方差 $\propto 1/\sigma_2^2 \approx 10^4$）。取 $\lambda = 0.1$，岭回归将小奇异值方向的贡献压缩，系数 $\frac{0.01^2}{0.01^2 + 0.1} \approx 0.001$，从而获得稳定的估计。

---

## 29.4 线性判别分析（LDA）

<div class="context-flow" markdown>

**核心**：Fisher 判别 = 找方向 $\mathbf{w}$ 使类间散布/类内散布之比最大 → 广义特征值问题 $S_B\mathbf{w} = \lambda S_W\mathbf{w}$

**链接**：Ch6 广义特征值 · Ch7 Rayleigh 商

</div>

线性判别分析通过求解广义特征值问题来找到最佳的分类投影方向。

!!! definition "定义 29.7 (类内与类间散布矩阵)"
    设有 $C$ 个类别，第 $c$ 类有 $n_c$ 个样本，类均值 $\boldsymbol{\mu}_c$，总均值 $\boldsymbol{\mu}$。**类内散布矩阵**（within-class scatter matrix）和**类间散布矩阵**（between-class scatter matrix）分别为

    $$
    S_W = \sum_{c=1}^C \sum_{\mathbf{x} \in \text{class } c} (\mathbf{x} - \boldsymbol{\mu}_c)(\mathbf{x} - \boldsymbol{\mu}_c)^T, \quad S_B = \sum_{c=1}^C n_c (\boldsymbol{\mu}_c - \boldsymbol{\mu})(\boldsymbol{\mu}_c - \boldsymbol{\mu})^T.
    $$

    总散布矩阵 $S_T = S_W + S_B$。注意 $\operatorname{rank}(S_B) \le C - 1$。

!!! theorem "定理 29.6 (Fisher 判别准则)"
    Fisher 线性判别寻找投影方向 $\mathbf{w}$ 使 **Fisher 准则函数**最大化：

    $$
    J(\mathbf{w}) = \frac{\mathbf{w}^T S_B \mathbf{w}}{\mathbf{w}^T S_W \mathbf{w}}.
    $$

    最优 $\mathbf{w}$ 是广义特征值问题 $S_B \mathbf{w} = \lambda S_W \mathbf{w}$ 的最大特征值对应的特征向量。当 $S_W$ 可逆时，等价于 $S_W^{-1}S_B$ 的最大特征值对应特征向量。

    对于多类（$C > 2$）问题，前 $C - 1$ 个广义特征向量给出最佳的 $(C-1)$ 维判别子空间。

??? proof "证明"
    $J(\mathbf{w})$ 是广义 Rayleigh 商。对其求导并令梯度为零：

    $$
    \frac{\partial}{\partial \mathbf{w}}\frac{\mathbf{w}^T S_B\mathbf{w}}{\mathbf{w}^T S_W\mathbf{w}} = 0 \implies S_B\mathbf{w} = J(\mathbf{w}) S_W \mathbf{w}.
    $$

    这正是广义特征值问题 $S_B\mathbf{w} = \lambda S_W\mathbf{w}$，$J(\mathbf{w})$ 的最大值为最大广义特征值 $\lambda_1$。由于 $\operatorname{rank}(S_B) \le C - 1$，至多有 $C - 1$ 个非零广义特征值。$\blacksquare$

!!! theorem "定理 29.7 (二分类 LDA 的闭式解)"
    对于二分类问题（$C = 2$），$S_B = n_1 n_2/(n_1+n_2) \cdot (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)^T$（秩 1）。最优判别方向为

    $$
    \mathbf{w}^* \propto S_W^{-1}(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2).
    $$

    分类规则：将新样本 $\mathbf{x}$ 投影到 $\mathbf{w}^*$ 上，根据投影值与阈值 $\frac{1}{2}\mathbf{w}^{*T}(\boldsymbol{\mu}_1 + \boldsymbol{\mu}_2)$ 的比较进行分类。

??? proof "证明"
    由于 $S_B\mathbf{w} \propto (\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)^T\mathbf{w}$，$S_B\mathbf{w}$ 总是 $(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)$ 方向。从 $S_B\mathbf{w} = \lambda S_W\mathbf{w}$ 得 $\mathbf{w} \propto S_W^{-1}(\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2)$。$\blacksquare$

!!! example "例 29.5"
    **二分类 LDA 计算。** 设两类二维数据，类 1 样本 $(0,0)^T$, $(2,2)^T$（均值 $\boldsymbol{\mu}_1 = (1,1)^T$），类 2 样本 $(4,2)^T$, $(6,4)^T$（均值 $\boldsymbol{\mu}_2 = (5,3)^T$）。

    $$
    S_W = \begin{pmatrix}1&1\\1&1\end{pmatrix} + \begin{pmatrix}1&1\\1&1\end{pmatrix} = \begin{pmatrix}2&2\\2&2\end{pmatrix}.
    $$

    $S_W$ 奇异，加正则化 $S_W + 0.01I$。$\boldsymbol{\mu}_1 - \boldsymbol{\mu}_2 = (-4,-2)^T$。最优方向 $\mathbf{w}^* \propto (S_W + 0.01I)^{-1}(-4,-2)^T$。投影后两类的均值分离被最大化。

---

## 29.5 核方法与特征空间

<div class="context-flow" markdown>

**核心思想**：通过非线性映射 $\phi: \mathbb{R}^p \to \mathcal{H}$（高维/无穷维）将数据线性化 → 核技巧：$k(\mathbf{x},\mathbf{y}) = \langle\phi(\mathbf{x}),\phi(\mathbf{y})\rangle$ 避免显式计算 $\phi$ → Gram 矩阵 $K$ 取代 $X^TX$

**链接**：Ch8 内积空间 · Ch16 正定矩阵 · Mercer 定理

</div>

核方法通过隐式的高维映射将线性方法推广到非线性问题，其数学基础是正定核与再生核 Hilbert 空间。

!!! definition "定义 29.8 (正定核与 Gram 矩阵)"
    函数 $k: \mathbb{R}^p \times \mathbb{R}^p \to \mathbb{R}$ 称为**正定核**（positive definite kernel），若对任意 $n$ 个点 $\mathbf{x}_1, \ldots, \mathbf{x}_n$，**Gram 矩阵**

    $$
    K = [k(\mathbf{x}_i, \mathbf{x}_j)]_{n \times n}
    $$

    是对称半正定矩阵。常用核函数包括：

    - **线性核**：$k(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T\mathbf{y}$。
    - **多项式核**：$k(\mathbf{x}, \mathbf{y}) = (\mathbf{x}^T\mathbf{y} + c)^d$，$c \ge 0$，$d \in \mathbb{N}$。
    - **高斯（RBF）核**：$k(\mathbf{x}, \mathbf{y}) = \exp\!\left(-\frac{\|\mathbf{x}-\mathbf{y}\|^2}{2\sigma^2}\right)$。

!!! theorem "定理 29.8 (Mercer 定理，有限维版本)"
    对称函数 $k: \mathbb{R}^p \times \mathbb{R}^p \to \mathbb{R}$ 是正定核当且仅当存在映射 $\phi: \mathbb{R}^p \to \mathcal{H}$（$\mathcal{H}$ 为某 Hilbert 空间）使得

    $$
    k(\mathbf{x}, \mathbf{y}) = \langle \phi(\mathbf{x}), \phi(\mathbf{y}) \rangle_{\mathcal{H}}.
    $$

    等价地，$K$ 半正定 $\Leftrightarrow$ $K$ 可分解为 $K = \Phi\Phi^T$，其中 $\Phi$ 的第 $i$ 行为 $\phi(\mathbf{x}_i)^T$。

??? proof "证明"
    （充分性）若 $k(\mathbf{x},\mathbf{y}) = \langle\phi(\mathbf{x}),\phi(\mathbf{y})\rangle$，则对任意 $\mathbf{c} \in \mathbb{R}^n$，$\mathbf{c}^T K\mathbf{c} = \sum_{i,j}c_ic_j\langle\phi(\mathbf{x}_i),\phi(\mathbf{x}_j)\rangle = \|\sum_i c_i\phi(\mathbf{x}_i)\|^2 \ge 0$，故 $K$ 半正定。

    （必要性）由 $K$ 半正定，存在特征分解 $K = Q\Lambda Q^T$，取 $\Phi = Q\Lambda^{1/2}$，则 $K = \Phi\Phi^T$。更一般地，可以构造再生核 Hilbert 空间 $\mathcal{H}$，令 $\phi(\mathbf{x}) = k(\cdot, \mathbf{x})$。$\blacksquare$

!!! definition "定义 29.9 (核 PCA)"
    **核 PCA**（kernel PCA）在特征空间 $\mathcal{H}$ 中执行 PCA。设中心化 Gram 矩阵为

    $$
    \tilde{K} = HKH, \quad H = I_n - \frac{1}{n}\mathbf{1}_n\mathbf{1}_n^T.
    $$

    对 $\tilde{K}$ 做特征分解 $\tilde{K} = Q\Lambda Q^T$，核主成分为 $\boldsymbol{\alpha}_j = \frac{1}{\sqrt{\lambda_j}}\mathbf{q}_j$。新样本 $\mathbf{x}$ 在第 $j$ 个核主成分上的投影为

    $$
    z_j = \sum_{i=1}^n \alpha_{ji} \tilde{k}(\mathbf{x}_i, \mathbf{x}),
    $$

    其中 $\tilde{k}$ 为中心化核函数。

!!! example "例 29.6"
    **多项式核的特征映射。** 对于二次核 $k(\mathbf{x},\mathbf{y}) = (\mathbf{x}^T\mathbf{y})^2$（$\mathbf{x}, \mathbf{y} \in \mathbb{R}^2$），显式特征映射为

    $$
    \phi(x_1, x_2) = (x_1^2,\; \sqrt{2}\,x_1 x_2,\; x_2^2)^T.
    $$

    验证：$\phi(\mathbf{x})^T\phi(\mathbf{y}) = x_1^2 y_1^2 + 2x_1 x_2 y_1 y_2 + x_2^2 y_2^2 = (x_1 y_1 + x_2 y_2)^2 = (\mathbf{x}^T\mathbf{y})^2$。核技巧的优势在于：即使 $\phi$ 映射到极高维空间（如高斯核对应无穷维），计算仍只需评估核函数 $k(\mathbf{x},\mathbf{y})$。

---

## 29.6 支持向量机

<div class="context-flow" markdown>

**核心**：最大间隔分类器 = 约束优化（二次规划） → 对偶问题仅涉及内积 $\mathbf{x}_i^T\mathbf{x}_j$ → 核化 → 支持向量 = $\alpha_i > 0$ 的样本

**链接**：Ch25 优化 · 核方法(29.5)

</div>

支持向量机是核方法最成功的应用之一，其训练归结为一个凸二次规划问题。

!!! definition "定义 29.10 (硬间隔 SVM)"
    给定线性可分的训练数据 $\{(\mathbf{x}_i, y_i)\}_{i=1}^n$，$y_i \in \{+1, -1\}$，**硬间隔支持向量机**（hard-margin SVM）求解

    $$
    \min_{\mathbf{w}, b} \frac{1}{2}\|\mathbf{w}\|^2, \quad \text{s.t.} \quad y_i(\mathbf{w}^T\mathbf{x}_i + b) \ge 1, \quad i = 1, \ldots, n.
    $$

    几何间隔为 $2/\|\mathbf{w}\|$。**软间隔 SVM** 引入松弛变量 $\xi_i \ge 0$：

    $$
    \min_{\mathbf{w}, b, \boldsymbol{\xi}} \frac{1}{2}\|\mathbf{w}\|^2 + C\sum_{i=1}^n \xi_i, \quad \text{s.t.} \quad y_i(\mathbf{w}^T\mathbf{x}_i + b) \ge 1 - \xi_i.
    $$

!!! theorem "定理 29.9 (SVM 对偶问题)"
    SVM 的 Lagrange 对偶问题为

    $$
    \max_{\boldsymbol{\alpha}} \sum_{i=1}^n \alpha_i - \frac{1}{2}\sum_{i,j=1}^n \alpha_i \alpha_j y_i y_j \mathbf{x}_i^T\mathbf{x}_j, \quad \text{s.t.} \quad 0 \le \alpha_i \le C, \quad \sum_{i=1}^n \alpha_i y_i = 0.
    $$

    用矩阵形式写为

    $$
    \max_{\boldsymbol{\alpha}} \mathbf{1}^T\boldsymbol{\alpha} - \frac{1}{2}\boldsymbol{\alpha}^T Q\boldsymbol{\alpha}, \quad Q_{ij} = y_i y_j \mathbf{x}_i^T\mathbf{x}_j,
    $$

    其中 $Q = \operatorname{diag}(\mathbf{y}) \cdot X X^T \cdot \operatorname{diag}(\mathbf{y})$ 是半正定矩阵。核化 SVM 只需将 $\mathbf{x}_i^T\mathbf{x}_j$ 替换为 $k(\mathbf{x}_i, \mathbf{x}_j)$。

??? proof "证明"
    引入 Lagrange 乘子 $\alpha_i \ge 0$，拉格朗日函数为

    $$
    L = \frac{1}{2}\|\mathbf{w}\|^2 - \sum_i \alpha_i[y_i(\mathbf{w}^T\mathbf{x}_i + b) - 1].
    $$

    对 $\mathbf{w}$ 求导置零：$\mathbf{w} = \sum_i \alpha_i y_i \mathbf{x}_i$。对 $b$ 求导置零：$\sum_i \alpha_i y_i = 0$。代回 $L$：

    $$
    L = \sum_i \alpha_i - \frac{1}{2}\sum_{i,j}\alpha_i\alpha_j y_i y_j \mathbf{x}_i^T\mathbf{x}_j.
    $$

    $Q$ 半正定因为 $Q = D_y X X^T D_y$，其中 $D_y = \operatorname{diag}(\mathbf{y})$，故 $\mathbf{c}^T Q\mathbf{c} = \|X^T D_y\mathbf{c}\|^2 \ge 0$。$\blacksquare$

!!! example "例 29.7"
    **简单 SVM 求解。** 设两类数据：类 $+1$：$(2,1)^T$，$(3,2)^T$；类 $-1$：$(0,0)^T$，$(1,1)^T$。对偶问题的 $Q$ 矩阵：

    $$
    Q_{ij} = y_i y_j \mathbf{x}_i^T\mathbf{x}_j, \quad X X^T = \begin{pmatrix}5&8&0&3\\8&13&0&5\\0&0&0&0\\3&5&0&2\end{pmatrix}.
    $$

    由 KKT 条件，支持向量（$\alpha_i > 0$ 的样本）位于间隔边界上。对于这组数据，最优超平面将两类以最大间隔分开。解对偶 QP 后，$\mathbf{w} = \sum_i \alpha_i y_i \mathbf{x}_i$，决策函数 $f(\mathbf{x}) = \operatorname{sign}(\mathbf{w}^T\mathbf{x} + b)$。

---

## 29.7 神经网络中的线性代数

<div class="context-flow" markdown>

**核心**：前向传播 = 矩阵链乘 + 逐元素非线性 → 反向传播 = 链式法则中的 Jacobian 矩阵乘法 → 权重矩阵的条件数影响训练稳定性

**链接**：Ch2 矩阵乘法 · Ch6 特征值 · Ch15 范数

</div>

深度神经网络的每一层核心运算都是矩阵乘法，理解这一点对于高效训练和分析至关重要。

!!! definition "定义 29.11 (全连接神经网络的矩阵表示)"
    一个 $L$ 层**全连接神经网络**（fully-connected neural network）可表示为复合映射

    $$
    f(\mathbf{x}) = \sigma_L(W_L \sigma_{L-1}(W_{L-1} \cdots \sigma_1(W_1 \mathbf{x} + \mathbf{b}_1) \cdots + \mathbf{b}_{L-1}) + \mathbf{b}_L),
    $$

    其中 $W_\ell \in \mathbb{R}^{d_\ell \times d_{\ell-1}}$ 为第 $\ell$ 层权重矩阵，$\mathbf{b}_\ell \in \mathbb{R}^{d_\ell}$ 为偏置，$\sigma_\ell$ 为逐元素激活函数。对于 mini-batch 输入 $X \in \mathbb{R}^{d_0 \times m}$（$m$ 个样本），第 $\ell$ 层的前向计算为

    $$
    Z_\ell = W_\ell A_{\ell-1} + \mathbf{b}_\ell \mathbf{1}^T, \quad A_\ell = \sigma_\ell(Z_\ell),
    $$

    其中 $A_0 = X$。核心运算是矩阵乘法 $W_\ell A_{\ell-1}$。

!!! theorem "定理 29.10 (反向传播的矩阵形式)"
    设损失函数 $\mathcal{L}$ 对第 $\ell$ 层预激活值的梯度为 $\delta_\ell = \frac{\partial \mathcal{L}}{\partial Z_\ell} \in \mathbb{R}^{d_\ell \times m}$。则反向传播的递推关系为

    $$
    \delta_\ell = (W_{\ell+1}^T \delta_{\ell+1}) \odot \sigma_\ell'(Z_\ell),
    $$

    其中 $\odot$ 为逐元素乘法（Hadamard 积）。权重梯度为

    $$
    \frac{\partial \mathcal{L}}{\partial W_\ell} = \frac{1}{m}\delta_\ell A_{\ell-1}^T.
    $$

    核心运算仍是矩阵乘法：$W_{\ell+1}^T \delta_{\ell+1}$（维度 $d_\ell \times m$）和 $\delta_\ell A_{\ell-1}^T$（维度 $d_\ell \times d_{\ell-1}$）。

??? proof "证明"
    由链式法则，$\frac{\partial \mathcal{L}}{\partial Z_\ell} = \frac{\partial A_\ell}{\partial Z_\ell} \cdot \frac{\partial Z_{\ell+1}}{\partial A_\ell} \cdot \frac{\partial \mathcal{L}}{\partial Z_{\ell+1}}$。由于 $A_\ell = \sigma_\ell(Z_\ell)$（逐元素），$\frac{\partial A_\ell}{\partial Z_\ell}$ 是对角的（逐元素意义下），贡献 $\sigma_\ell'(Z_\ell)$。由 $Z_{\ell+1} = W_{\ell+1}A_\ell + \mathbf{b}_{\ell+1}\mathbf{1}^T$，$\frac{\partial Z_{\ell+1}}{\partial A_\ell} = W_{\ell+1}$，故 $\delta_\ell = (W_{\ell+1}^T\delta_{\ell+1})\odot\sigma_\ell'(Z_\ell)$。对于权重，$\frac{\partial \mathcal{L}}{\partial W_\ell} = \delta_\ell \frac{\partial Z_\ell}{\partial W_\ell}$，由 $Z_\ell = W_\ell A_{\ell-1}$，得 $\frac{\partial \mathcal{L}}{\partial W_\ell} = \frac{1}{m}\delta_\ell A_{\ell-1}^T$。$\blacksquare$

!!! theorem "定理 29.11 (梯度消失/爆炸与奇异值)"
    对于线性网络（$\sigma_\ell = \text{id}$），从第 $L$ 层到第 $\ell$ 层的梯度传播涉及矩阵积 $\prod_{k=\ell+1}^L W_k^T$。设 $W_k$ 的最大奇异值为 $\sigma_{\max}(W_k)$，最小奇异值为 $\sigma_{\min}(W_k)$。则：

    1. 若 $\sigma_{\max}(W_k) > 1$ 对多数 $k$ 成立，梯度可能**指数爆炸**。
    2. 若 $\sigma_{\max}(W_k) < 1$ 对多数 $k$ 成立，梯度可能**指数消失**。
    3. 当 $W_k$ 为正交矩阵（所有奇异值为 1）时，梯度范数保持不变——这是正交初始化和正交约束训练的理论基础。

??? proof "证明"
    对于线性网络，$\delta_\ell = \prod_{k=\ell+1}^L W_k^T \cdot \delta_L$。由奇异值的次可乘性，$\|\delta_\ell\| \le \prod_{k=\ell+1}^L \sigma_{\max}(W_k) \|\delta_L\|$。类似地，$\|\delta_\ell\| \ge \prod_{k=\ell+1}^L \sigma_{\min}(W_k) \|\delta_L\|$。当所有 $W_k$ 正交时，$\sigma_{\max} = \sigma_{\min} = 1$，故 $\|\delta_\ell\| = \|\delta_L\|$。$\blacksquare$

!!! example "例 29.8"
    **梯度传播分析。** 考虑 $5$ 层线性网络，每层权重 $W_k \in \mathbb{R}^{100 \times 100}$。

    **情形 1**：$W_k = 1.1 I$。经过 $4$ 层，梯度放大 $1.1^4 \approx 1.46$ 倍。经过 $50$ 层则放大 $1.1^{50} \approx 117$ 倍——梯度爆炸。

    **情形 2**：$W_k = 0.9 I$。经过 $50$ 层，梯度衰减至 $0.9^{50} \approx 0.005$——梯度消失。

    **情形 3**：$W_k$ 为随机正交矩阵。无论层数多少，$\|\delta_\ell\| = \|\delta_L\|$，梯度完美保持。这解释了为什么正交初始化（如 Saxe et al., 2014）能显著改善深层网络训练。

---

## 29.8 矩阵分解在推荐系统中的应用

<div class="context-flow" markdown>

**核心**：用户-物品评分矩阵 $R$ 稀疏且不完整 → 低秩近似 $R \approx UV^T$ → NMF 给出可解释的非负因子 → 正则化防过拟合

**链接**：Ch11 SVD · Ch10 矩阵分解 · Ch17 非负矩阵

</div>

推荐系统是矩阵分解在工业界最成功的应用之一，其核心思想是利用评分矩阵的低秩结构进行协同过滤。

!!! definition "定义 29.12 (矩阵分解模型)"
    设 $R \in \mathbb{R}^{m \times n}$ 为**用户-物品评分矩阵**，其中 $R_{ij}$ 为用户 $i$ 对物品 $j$ 的评分（大量缺失值）。**矩阵分解模型**（matrix factorization model）寻找低秩近似

    $$
    R \approx U V^T, \quad U \in \mathbb{R}^{m \times k}, \quad V \in \mathbb{R}^{n \times k},
    $$

    其中 $k \ll \min(m,n)$ 为**隐因子**（latent factor）个数。$U$ 的第 $i$ 行 $\mathbf{u}_i^T$ 为用户 $i$ 的隐向量，$V$ 的第 $j$ 行 $\mathbf{v}_j^T$ 为物品 $j$ 的隐向量，预测评分 $\hat{R}_{ij} = \mathbf{u}_i^T\mathbf{v}_j$。

!!! definition "定义 29.13 (正则化矩阵分解)"
    带正则化的矩阵分解求解如下优化问题：

    $$
    \min_{U, V} \sum_{(i,j) \in \Omega} (R_{ij} - \mathbf{u}_i^T\mathbf{v}_j)^2 + \lambda(\|U\|_F^2 + \|V\|_F^2),
    $$

    其中 $\Omega$ 为已观测评分的索引集，$\lambda > 0$ 为正则化参数。常用算法包括**交替最小二乘**（Alternating Least Squares, ALS）：

    - 固定 $V$，对每个用户 $i$ 更新 $\mathbf{u}_i = (V_{\Omega_i}^T V_{\Omega_i} + \lambda I)^{-1}V_{\Omega_i}^T \mathbf{r}_i$。
    - 固定 $U$，对每个物品 $j$ 更新 $\mathbf{v}_j = (U_{\Omega_j}^T U_{\Omega_j} + \lambda I)^{-1}U_{\Omega_j}^T \mathbf{r}_j$。

    其中 $\Omega_i$ 为用户 $i$ 评过分的物品集，$V_{\Omega_i}$ 为对应行构成的子矩阵。

!!! theorem "定理 29.12 (SVD 与矩阵补全)"
    若完整评分矩阵 $R$ 恰好秩为 $k$，且观测集 $\Omega$ 满足一定条件（如随机均匀采样且样本量 $|\Omega| \ge C \cdot (m+n)k \log^2(m+n)$），则 $R$ 可以通过核范数最小化精确恢复：

    $$
    \min_{M} \|M\|_* \quad \text{s.t.} \quad M_{ij} = R_{ij}, \quad (i,j) \in \Omega,
    $$

    其中 $\|M\|_* = \sum_j \sigma_j(M)$ 为核范数（奇异值之和），是秩函数的凸松弛。

??? proof "证明"
    （证明概要）核范数是秩函数在 $\{M : \|M\|_{\text{op}} \le 1\}$ 上的凸包。在 $R$ 满足非相干性条件（incoherence condition）且 $|\Omega|$ 足够大时，核范数最小化的解是唯一的且等于 $R$。这一结果由 Candes 和 Recht（2009）证明，其关键工具是矩阵浓度不等式和对偶证书构造。完整证明需要随机矩阵理论（Ch23）的工具。$\blacksquare$

!!! definition "定义 29.14 (非负矩阵分解 NMF)"
    **非负矩阵分解**（Nonnegative Matrix Factorization, NMF）要求因子矩阵非负：

    $$
    R \approx WH, \quad W \in \mathbb{R}_{\ge 0}^{m \times k}, \quad H \in \mathbb{R}_{\ge 0}^{k \times n}.
    $$

    非负约束使得分解具有"部分-整体"的可解释性：每个隐因子代表一种"主题"或"模式"，用户/物品由这些主题的非负加权组合表示。常用的乘法更新规则为

    $$
    W \leftarrow W \odot \frac{RH^T}{WHH^T}, \quad H \leftarrow H \odot \frac{W^TR}{W^TWH},
    $$

    其中除法为逐元素运算。

!!! theorem "定理 29.13 (NMF 乘法更新的收敛性)"
    Lee-Seung 乘法更新规则使目标函数 $\|R - WH\|_F^2$ 单调不增。具体而言，定义辅助函数

    $$
    G(h, h') = F(h') + F'(h')(h - h') + \frac{(WHH^T)_{ij}}{h'_{ij}}(h - h')^2,
    $$

    则 $G(h, h') \ge F(h)$ 且 $G(h', h') = F(h')$。乘法更新等价于 $h^{(t+1)} = \arg\min_h G(h, h^{(t)})$，由辅助函数性质保证 $F(h^{(t+1)}) \le F(h^{(t)})$。

??? proof "证明"
    对 $\|R - WH\|_F^2$ 关于 $H_{ij}$ 求导：$\frac{\partial F}{\partial H_{ij}} = -2(W^TR)_{ij} + 2(W^TWH)_{ij}$。乘法更新可写为

    $$
    H_{ij} \leftarrow H_{ij} \cdot \frac{(W^TR)_{ij}}{(W^TWH)_{ij}}.
    $$

    定义辅助函数（二次上界），在 $h = h'$ 处相切。最小化辅助函数给出更新规则。由于每步都不增加辅助函数值，而辅助函数是目标函数的上界，故目标函数单调不增。有界递减序列收敛。$\blacksquare$

!!! example "例 29.9"
    **简单矩阵分解示例。** 设 $3$ 个用户对 $4$ 部电影的评分矩阵（? 表示缺失）：

    $$
    R = \begin{pmatrix}5&3&?&1\\4&?&?&1\\1&1&?&5\end{pmatrix}.
    $$

    设隐因子 $k = 2$。初始化 $U \in \mathbb{R}^{3 \times 2}$，$V \in \mathbb{R}^{4 \times 2}$。经过 ALS 迭代，假设收敛到

    $$
    U \approx \begin{pmatrix}2.1&0.3\\1.8&0.2\\0.5&2.0\end{pmatrix}, \quad V \approx \begin{pmatrix}2.3&0.2\\1.4&0.1\\0.8&0.5\\0.1&2.4\end{pmatrix}.
    $$

    预测缺失评分 $\hat{R}_{13} = \mathbf{u}_1^T\mathbf{v}_3 \approx 2.1 \times 0.8 + 0.3 \times 0.5 = 1.83$，$\hat{R}_{22} = 1.8 \times 1.4 + 0.2 \times 0.1 = 2.54$。隐因子可解释为：第一维可能对应"动作片偏好"，第二维对应"文艺片偏好"。

## 练习题

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

本章横跨了古典统计学与现代机器学习，全方位展示了线性代数作为数据科学“操作系统”的基石作用：

1. **多元统计的几何化**：用协方差矩阵和马氏距离将随机变量的联合分布直观化为高维椭球，将统计学里的不确定性转化为矩阵的特征值缩放。
2. **主成分分析 (PCA)**：从方差最大化和最优低秩逼近的双重视角推导了 PCA，证明了它就是中心化数据 SVD 分解的直接推论，并探讨了它在数据降维去噪中的统治地位。
3. **核方法 (Kernel Methods)**：通过 Mercer 定理和半正定核矩阵，将原本只能处理线性问题的内积工具，神奇地跃迁到了无穷维的非线性特征空间。
4. **深度学习的反向传播**：将人工神经网络的层次结构抽象为交替的仿射变换与非线性激活，揭示了梯度下降在本质上就是雅可比矩阵链式乘积的连续回传。
5. **推荐系统的矩阵分解**：分析了协同过滤中如何利用非负矩阵分解 (NMF) 等低秩技术，在极度稀疏的大矩阵中挖掘出隐藏的用户-物品双向潜变量特征。
