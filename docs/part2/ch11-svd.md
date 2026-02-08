# 第 11 章 奇异值分解

奇异值分解（Singular Value Decomposition, SVD）是线性代数中最重要的矩阵分解之一。与特征值分解不同，SVD 适用于**任意**矩阵——不要求方阵，也不要求可对角化。SVD 将一个 $m \times n$ 矩阵分解为两个正交矩阵与一个对角矩阵的乘积，清晰地揭示了线性变换的几何本质：旋转、缩放、再旋转。SVD 在数据科学、信号处理、数值计算等领域有着极其广泛的应用，是现代应用数学的核心工具。

---

## 11.1 奇异值的定义

对于任意 $m \times n$ 实矩阵 $A$，矩阵 $A^T A$ 是一个 $n \times n$ 的半正定对称矩阵，因此其所有特征值都是非负实数。奇异值正是从这些特征值中提取出来的。

!!! definition "定义 11.1 (奇异值 Singular Value)"
    设 $A$ 为 $m \times n$ 实矩阵，$A^T A$ 的特征值为 $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_n \ge 0$。则 $A$ 的**奇异值**（singular values）定义为
    $$
    \sigma_i = \sqrt{\lambda_i}, \quad i = 1, 2, \ldots, n.
    $$
    按降序排列 $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_n \ge 0$。

!!! note "注"
    矩阵 $A^T A$ 半正定的证明：对任意 $\mathbf{x} \in \mathbb{R}^n$，
    $$
    \mathbf{x}^T (A^T A) \mathbf{x} = (A\mathbf{x})^T (A\mathbf{x}) = \|A\mathbf{x}\|^2 \ge 0.
    $$
    因此 $A^T A$ 的特征值均非负，奇异值定义合理。

!!! theorem "定理 11.1 (奇异值与矩阵范数的关系)"
    设 $A$ 为 $m \times n$ 实矩阵，$\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_n \ge 0$ 为其奇异值，则：

    1. $\|A\|_2 = \sigma_1$（谱范数等于最大奇异值）；
    2. $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2 + \cdots + \sigma_n^2}$（Frobenius 范数等于奇异值的平方和的平方根）；
    3. $\operatorname{rank}(A) = $ 非零奇异值的个数。

??? proof "证明"
    **(1)** 谱范数定义为
    $$
    \|A\|_2 = \max_{\|\mathbf{x}\|=1} \|A\mathbf{x}\|.
    $$
    由于 $\|A\mathbf{x}\|^2 = \mathbf{x}^T A^T A \mathbf{x}$，而 $A^T A$ 是对称半正定矩阵，设其特征值为 $\lambda_1 \ge \cdots \ge \lambda_n \ge 0$，对应的单位正交特征向量为 $\mathbf{v}_1, \ldots, \mathbf{v}_n$。任意单位向量 $\mathbf{x} = \sum c_i \mathbf{v}_i$（$\sum c_i^2 = 1$），则
    $$
    \|A\mathbf{x}\|^2 = \sum_{i=1}^n \lambda_i c_i^2 \le \lambda_1 \sum c_i^2 = \lambda_1.
    $$
    取 $\mathbf{x} = \mathbf{v}_1$ 时等号成立。因此 $\|A\|_2 = \sqrt{\lambda_1} = \sigma_1$。

    **(2)** $\|A\|_F^2 = \operatorname{tr}(A^T A) = \sum_{i=1}^n \lambda_i = \sum_{i=1}^n \sigma_i^2$。

    **(3)** $\operatorname{rank}(A) = \operatorname{rank}(A^T A)$（因为 $A\mathbf{x} = \mathbf{0}$ 当且仅当 $A^T A \mathbf{x} = \mathbf{0}$），而对称矩阵的秩等于非零特征值的个数。

!!! definition "定义 11.2 (左奇异向量与右奇异向量)"
    设 $A$ 为 $m \times n$ 实矩阵。

    - $A^T A$ 的单位正交特征向量 $\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n \in \mathbb{R}^n$ 称为 $A$ 的**右奇异向量**（right singular vectors）；
    - $A A^T$ 的单位正交特征向量 $\mathbf{u}_1, \mathbf{u}_2, \ldots, \mathbf{u}_m \in \mathbb{R}^m$ 称为 $A$ 的**左奇异向量**（left singular vectors）。

!!! theorem "定理 11.2 ($A^TA$ 与 $AA^T$ 的非零特征值相同)"
    设 $A$ 为 $m \times n$ 实矩阵，则 $A^T A$ 与 $A A^T$ 具有相同的非零特征值（计重数）。

??? proof "证明"
    设 $\lambda \neq 0$ 是 $A^T A$ 的特征值，$\mathbf{v}$ 为对应的特征向量，则 $A^T A \mathbf{v} = \lambda \mathbf{v}$。令 $\mathbf{u} = A\mathbf{v}$，则 $\mathbf{u} \neq \mathbf{0}$（因为若 $A\mathbf{v} = \mathbf{0}$，则 $\lambda \mathbf{v} = A^T A \mathbf{v} = \mathbf{0}$，与 $\lambda \neq 0$、$\mathbf{v} \neq \mathbf{0}$ 矛盾），且
    $$
    A A^T \mathbf{u} = A A^T (A \mathbf{v}) = A (A^T A \mathbf{v}) = A (\lambda \mathbf{v}) = \lambda (A\mathbf{v}) = \lambda \mathbf{u}.
    $$
    因此 $\lambda$ 也是 $A A^T$ 的特征值。反向同理可证。

!!! example "例 11.1"
    求矩阵 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的奇异值。

    **解：** 计算
    $$
    A^T A = \begin{pmatrix} 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}.
    $$
    特征多项式为
    $$
    \det(A^T A - \lambda I) = (2 - \lambda)^2 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 3)(\lambda - 1).
    $$
    特征值为 $\lambda_1 = 3$，$\lambda_2 = 1$。

    因此奇异值为 $\sigma_1 = \sqrt{3}$，$\sigma_2 = 1$。

---

## 11.2 奇异值分解定理

奇异值分解定理是本章的核心结果，它保证任意矩阵都可以被分解为标准形式。

!!! theorem "定理 11.3 (奇异值分解定理 SVD Theorem)"
    设 $A$ 为 $m \times n$ 实矩阵，$\operatorname{rank}(A) = r$。则存在 $m \times m$ 正交矩阵 $U$、$n \times n$ 正交矩阵 $V$ 以及 $m \times n$ 矩阵 $\Sigma$，使得
    $$
    A = U \Sigma V^T,
    $$
    其中 $\Sigma$ 的形式为
    $$
    \Sigma = \begin{pmatrix} \sigma_1 & & & 0 & \cdots & 0 \\ & \sigma_2 & & 0 & \cdots & 0 \\ & & \ddots & \vdots & & \vdots \\ & & & \sigma_r & \cdots & 0 \\ 0 & \cdots & 0 & 0 & \cdots & 0 \\ \vdots & & \vdots & \vdots & & \vdots \\ 0 & \cdots & 0 & 0 & \cdots & 0 \end{pmatrix}_{m \times n},
    $$
    且 $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$。

??? proof "证明"
    **构造性证明：**

    **第一步：** 由于 $A^T A$ 是 $n \times n$ 对称半正定矩阵，由谱定理，存在正交矩阵 $V = [\mathbf{v}_1, \ldots, \mathbf{v}_n]$ 使得
    $$
    A^T A = V \operatorname{diag}(\lambda_1, \ldots, \lambda_n) V^T,
    $$
    其中 $\lambda_1 \ge \cdots \ge \lambda_r > 0 = \lambda_{r+1} = \cdots = \lambda_n$。令 $\sigma_i = \sqrt{\lambda_i}$（$i = 1, \ldots, r$）。

    **第二步：** 对 $i = 1, \ldots, r$，定义
    $$
    \mathbf{u}_i = \frac{1}{\sigma_i} A \mathbf{v}_i.
    $$
    验证 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ 是正交单位组：
    $$
    \mathbf{u}_i^T \mathbf{u}_j = \frac{1}{\sigma_i \sigma_j} \mathbf{v}_i^T A^T A \mathbf{v}_j = \frac{1}{\sigma_i \sigma_j} \mathbf{v}_i^T (\lambda_j \mathbf{v}_j) = \frac{\lambda_j}{\sigma_i \sigma_j} \delta_{ij} = \delta_{ij}.
    $$

    **第三步：** 将 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ 扩充为 $\mathbb{R}^m$ 的标准正交基 $\{\mathbf{u}_1, \ldots, \mathbf{u}_m\}$，令 $U = [\mathbf{u}_1, \ldots, \mathbf{u}_m]$。

    **第四步：** 验证 $A = U \Sigma V^T$。由构造，$A\mathbf{v}_i = \sigma_i \mathbf{u}_i$（$i \le r$），且对 $i > r$，$A\mathbf{v}_i = \mathbf{0}$（因为 $\|A\mathbf{v}_i\|^2 = \mathbf{v}_i^T A^T A \mathbf{v}_i = \lambda_i = 0$）。因此
    $$
    AV = U\Sigma, \quad \text{即} \quad A = U\Sigma V^T. \qquad \blacksquare
    $$

!!! definition "定义 11.3 (SVD 的各分量)"
    在 SVD 分解 $A = U\Sigma V^T$ 中：

    - $U = [\mathbf{u}_1, \ldots, \mathbf{u}_m]$ 为 $m \times m$ 正交矩阵，其列为**左奇异向量**；
    - $V = [\mathbf{v}_1, \ldots, \mathbf{v}_n]$ 为 $n \times n$ 正交矩阵，其列为**右奇异向量**；
    - $\Sigma$ 为 $m \times n$ 对角矩阵，对角元素为奇异值 $\sigma_1 \ge \cdots \ge \sigma_r > 0$。

!!! theorem "定理 11.4 (SVD 的外积展开形式)"
    设 $A = U\Sigma V^T$ 为 $A$ 的 SVD，$\operatorname{rank}(A) = r$，则
    $$
    A = \sum_{i=1}^{r} \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$
    即 $A$ 可以表示为 $r$ 个秩一矩阵的加权和。

??? proof "证明"
    由 $A = U\Sigma V^T$，展开得
    $$
    A = \sum_{i=1}^{\min(m,n)} \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$
    当 $i > r$ 时 $\sigma_i = 0$，因此求和实际只到 $r$。每个 $\mathbf{u}_i \mathbf{v}_i^T$ 是 $m \times n$ 矩阵，且 $\operatorname{rank}(\mathbf{u}_i \mathbf{v}_i^T) = 1$。

!!! example "例 11.2"
    对矩阵 $A = \begin{pmatrix} 3 & 0 \\ 0 & 2 \end{pmatrix}$ 求 SVD。

    **解：**
    $$
    A^T A = \begin{pmatrix} 9 & 0 \\ 0 & 4 \end{pmatrix},
    $$
    特征值 $\lambda_1 = 9$，$\lambda_2 = 4$；特征向量 $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$，$\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$。

    奇异值 $\sigma_1 = 3$，$\sigma_2 = 2$。

    左奇异向量：$\mathbf{u}_1 = \frac{1}{3}A\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$，$\mathbf{u}_2 = \frac{1}{2}A\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$。

    因此
    $$
    A = \begin{pmatrix}1&0\\0&1\end{pmatrix}\begin{pmatrix}3&0\\0&2\end{pmatrix}\begin{pmatrix}1&0\\0&1\end{pmatrix} = I \cdot \Sigma \cdot I.
    $$
    对角矩阵本身就是其 SVD（$U = V = I$）。

---

## 11.3 SVD 的几何意义

SVD 的几何解释极其直观：任何线性变换都可以分解为**旋转（或反射）**、**缩放**、**再旋转（或反射）**三步。

!!! definition "定义 11.4 (正交变换的几何意义)"
    正交矩阵 $Q$（$Q^TQ = I$）所表示的线性变换是保距的（isometric），即 $\|Q\mathbf{x}\| = \|\mathbf{x}\|$。当 $\det Q = 1$ 时为旋转，当 $\det Q = -1$ 时为旋转加反射。

!!! theorem "定理 11.5 (线性变换的几何分解)"
    设 $A = U\Sigma V^T$ 为 $m \times n$ 矩阵 $A$ 的 SVD。则线性变换 $\mathbf{x} \mapsto A\mathbf{x}$ 可以分解为三步：

    1. $V^T$：在 $\mathbb{R}^n$ 中旋转（将标准基变为右奇异向量方向）；
    2. $\Sigma$：沿坐标轴方向按 $\sigma_i$ 缩放（从 $\mathbb{R}^n$ 到 $\mathbb{R}^m$）；
    3. $U$：在 $\mathbb{R}^m$ 中旋转（将坐标轴方向变为左奇异向量方向）。

??? proof "证明"
    对任意 $\mathbf{x} \in \mathbb{R}^n$，令 $\mathbf{y} = V^T \mathbf{x}$。由于 $V$ 是正交矩阵，$\|\mathbf{y}\| = \|\mathbf{x}\|$，因此第一步是等距变换（旋转或反射）。

    第二步，$\mathbf{z} = \Sigma \mathbf{y}$ 将 $\mathbf{y}$ 的第 $i$ 个分量乘以 $\sigma_i$，这是沿坐标轴的缩放。

    第三步，$A\mathbf{x} = U\mathbf{z}$，$U$ 是正交矩阵，因此也是等距变换。

    综合三步，单位球 $\{\mathbf{x} : \|\mathbf{x}\| = 1\}$ 在 $V^T$ 下仍为单位球，在 $\Sigma$ 下变为以 $\sigma_i$ 为半轴长的超椭球，在 $U$ 下旋转到最终位置。

!!! note "注"
    在二维情形中，单位圆 $\|\mathbf{x}\| = 1$ 经过线性变换 $A$ 后变为一个椭圆，椭圆的半长轴和半短轴长度恰好是 $\sigma_1$ 和 $\sigma_2$，方向分别为 $\mathbf{u}_1$ 和 $\mathbf{u}_2$。

!!! example "例 11.3"
    设 $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$。描述 $A$ 对单位圆的作用。

    **解：** $A$ 已经是对角矩阵，其 SVD 为 $A = I \cdot \begin{pmatrix}2&0\\0&1\end{pmatrix} \cdot I$。

    单位圆 $x^2 + y^2 = 1$ 在 $A$ 作用下变为椭圆 $\frac{x^2}{4} + y^2 = 1$，长半轴为 $\sigma_1 = 2$（沿 $x$ 轴方向），短半轴为 $\sigma_2 = 1$（沿 $y$ 轴方向）。

---

## 11.4 SVD 的计算方法

本节通过详细例题展示 SVD 的计算步骤。

!!! definition "定义 11.5 (SVD 计算的标准流程)"
    给定 $m \times n$ 矩阵 $A$，SVD 的计算分为以下步骤：

    1. 计算 $A^T A$；
    2. 求 $A^T A$ 的特征值 $\lambda_1 \ge \cdots \ge \lambda_n \ge 0$ 和对应的单位正交特征向量 $\mathbf{v}_1, \ldots, \mathbf{v}_n$；
    3. 奇异值 $\sigma_i = \sqrt{\lambda_i}$；
    4. 对每个 $\sigma_i > 0$，计算 $\mathbf{u}_i = \frac{1}{\sigma_i} A \mathbf{v}_i$；
    5. 将 $\{\mathbf{u}_1, \ldots, \mathbf{u}_r\}$ 扩充为 $\mathbb{R}^m$ 的标准正交基；
    6. 组装 $U$、$\Sigma$、$V$。

!!! example "例 11.4"
    求矩阵 $A = \begin{pmatrix} 1 & 1 \\ 1 & -1 \\ 0 & 0 \end{pmatrix}$ 的 SVD。

    **解：**

    **第一步：** 计算 $A^T A$：
    $$
    A^T A = \begin{pmatrix} 1&1&0 \\ 1&-1&0 \end{pmatrix}\begin{pmatrix} 1&1 \\ 1&-1 \\ 0&0 \end{pmatrix} = \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix}.
    $$

    **第二步：** 特征值 $\lambda_1 = \lambda_2 = 2$，取正交特征向量 $\mathbf{v}_1 = \begin{pmatrix}1\\0\end{pmatrix}$，$\mathbf{v}_2 = \begin{pmatrix}0\\1\end{pmatrix}$。

    **第三步：** 奇异值 $\sigma_1 = \sigma_2 = \sqrt{2}$。

    **第四步：** 计算左奇异向量：
    $$
    \mathbf{u}_1 = \frac{1}{\sqrt{2}} A \mathbf{v}_1 = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\1\\0\end{pmatrix}, \quad
    \mathbf{u}_2 = \frac{1}{\sqrt{2}} A \mathbf{v}_2 = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-1\\0\end{pmatrix}.
    $$

    **第五步：** 扩充 $\{\mathbf{u}_1, \mathbf{u}_2\}$，取 $\mathbf{u}_3 = \begin{pmatrix}0\\0\\1\end{pmatrix}$。

    **第六步：** 因此
    $$
    A = \begin{pmatrix} \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 0 \\ \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 0 \\ 0 & 0 & 1 \end{pmatrix} \begin{pmatrix} \sqrt{2} & 0 \\ 0 & \sqrt{2} \\ 0 & 0 \end{pmatrix} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}^T.
    $$

!!! example "例 11.5"
    求矩阵 $A = \begin{pmatrix} 4 & 0 \\ 3 & -5 \end{pmatrix}$ 的 SVD。

    **解：**

    **第一步：**
    $$
    A^T A = \begin{pmatrix} 4&3 \\ 0&-5 \end{pmatrix}\begin{pmatrix} 4&0 \\ 3&-5 \end{pmatrix} = \begin{pmatrix} 25 & -15 \\ -15 & 25 \end{pmatrix}.
    $$

    **第二步：** 特征多项式：
    $$
    \det(A^T A - \lambda I) = (25-\lambda)^2 - 225 = \lambda^2 - 50\lambda + 400 = (\lambda - 40)(\lambda - 10).
    $$
    特征值 $\lambda_1 = 40$，$\lambda_2 = 10$。

    对 $\lambda_1 = 40$：$(A^T A - 40I)\mathbf{v} = \mathbf{0}$ 给出 $\begin{pmatrix}-15&-15\\-15&-15\end{pmatrix}\mathbf{v} = \mathbf{0}$，得 $\mathbf{v}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\-1\end{pmatrix}$。

    对 $\lambda_2 = 10$：$(A^T A - 10I)\mathbf{v} = \mathbf{0}$ 给出 $\begin{pmatrix}15&-15\\-15&15\end{pmatrix}\mathbf{v} = \mathbf{0}$，得 $\mathbf{v}_2 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$。

    **第三步：** $\sigma_1 = \sqrt{40} = 2\sqrt{10}$，$\sigma_2 = \sqrt{10}$。

    **第四步：**
    $$
    \mathbf{u}_1 = \frac{1}{2\sqrt{10}} A \mathbf{v}_1 = \frac{1}{2\sqrt{10}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}4\\8\end{pmatrix} = \frac{1}{\sqrt{5}}\begin{pmatrix}1\\2\end{pmatrix},
    $$
    $$
    \mathbf{u}_2 = \frac{1}{\sqrt{10}} A \mathbf{v}_2 = \frac{1}{\sqrt{10}} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}4\\-2\end{pmatrix} = \frac{1}{\sqrt{5}}\begin{pmatrix}2\\-1\end{pmatrix}.
    $$

    **第五步：** 组装：
    $$
    A = \frac{1}{\sqrt{5}}\begin{pmatrix}1&2\\2&-1\end{pmatrix} \begin{pmatrix}2\sqrt{10}&0\\0&\sqrt{10}\end{pmatrix} \frac{1}{\sqrt{2}}\begin{pmatrix}1&1\\-1&1\end{pmatrix}^T.
    $$

---

## 11.5 紧凑 SVD 与截断 SVD

当矩阵的秩远小于其维度时，完整的 SVD 包含了大量冗余信息。紧凑 SVD 和截断 SVD 提供了更经济的表示。

!!! definition "定义 11.6 (紧凑 SVD / Compact SVD)"
    设 $A = U\Sigma V^T$ 为 $m \times n$ 矩阵 $A$ 的完整 SVD，$\operatorname{rank}(A) = r$。令 $U_r = [\mathbf{u}_1, \ldots, \mathbf{u}_r]$（$m \times r$），$\Sigma_r = \operatorname{diag}(\sigma_1, \ldots, \sigma_r)$（$r \times r$），$V_r = [\mathbf{v}_1, \ldots, \mathbf{v}_r]$（$n \times r$），则
    $$
    A = U_r \Sigma_r V_r^T
    $$
    称为 $A$ 的**紧凑 SVD**（compact SVD）或**经济 SVD**（economy SVD）。

!!! definition "定义 11.7 (截断 SVD / Truncated SVD)"
    设 $\operatorname{rank}(A) = r$，对 $k < r$，令
    $$
    A_k = \sum_{i=1}^{k} \sigma_i \mathbf{u}_i \mathbf{v}_i^T = U_k \Sigma_k V_k^T,
    $$
    其中 $U_k$、$\Sigma_k$、$V_k$ 分别取前 $k$ 列/行。$A_k$ 称为 $A$ 的**秩 $k$ 截断 SVD**（rank-$k$ truncated SVD）。

!!! theorem "定理 11.6 (截断 SVD 的近似误差)"
    设 $A$ 的奇异值为 $\sigma_1 \ge \cdots \ge \sigma_r > 0$，$A_k$ 为秩 $k$ 截断近似，则
    $$
    \|A - A_k\|_2 = \sigma_{k+1}, \qquad \|A - A_k\|_F = \sqrt{\sigma_{k+1}^2 + \cdots + \sigma_r^2}.
    $$

??? proof "证明"
    由 $A = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ 以及 $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$，得
    $$
    A - A_k = \sum_{i=k+1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T.
    $$
    这本身就是 $A - A_k$ 的 SVD（奇异值为 $\sigma_{k+1}, \ldots, \sigma_r$），因此
    $$
    \|A - A_k\|_2 = \sigma_{k+1}, \qquad \|A - A_k\|_F = \sqrt{\sum_{i=k+1}^r \sigma_i^2}. \qquad \blacksquare
    $$

!!! example "例 11.6"
    设 $A$ 的奇异值为 $10, 5, 2, 0.1$。用秩 2 截断近似时，谱范数误差和 Frobenius 范数误差分别是多少？相对 Frobenius 误差是多少？

    **解：**
    $$
    \|A - A_2\|_2 = \sigma_3 = 2, \qquad \|A - A_2\|_F = \sqrt{2^2 + 0.1^2} = \sqrt{4.01} \approx 2.0025.
    $$
    $$
    \|A\|_F = \sqrt{100 + 25 + 4 + 0.01} = \sqrt{129.01} \approx 11.359.
    $$
    相对误差为 $\frac{\|A - A_2\|_F}{\|A\|_F} \approx \frac{2.0025}{11.359} \approx 17.6\%$。

---

## 11.6 低秩近似

SVD 提供了最优低秩近似，这一结论由 Eckart-Young-Mirsky 定理精确表述。

!!! theorem "定理 11.7 (Eckart-Young-Mirsky 定理)"
    设 $A$ 为 $m \times n$ 实矩阵，奇异值为 $\sigma_1 \ge \cdots \ge \sigma_r > 0$。对任意秩不超过 $k$（$k < r$）的矩阵 $B$，有
    $$
    \|A - A_k\|_2 \le \|A - B\|_2, \qquad \|A - A_k\|_F \le \|A - B\|_F,
    $$
    其中 $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ 为 $A$ 的秩 $k$ 截断 SVD。即 $A_k$ 是 $A$ 在谱范数和 Frobenius 范数下的**最佳秩 $k$ 近似**。

??? proof "证明"
    **谱范数情形的证明：**

    设 $B$ 为任意秩不超过 $k$ 的矩阵。$\ker(B)$ 的维数至少为 $n - k$。考虑子空间 $W = \operatorname{span}\{\mathbf{v}_1, \ldots, \mathbf{v}_{k+1}\}$，$\dim W = k + 1$。

    由维数公式，$W \cap \ker(B) \neq \{\mathbf{0}\}$（因为 $(k+1) + (n-k) = n + 1 > n$）。取单位向量 $\mathbf{w} \in W \cap \ker(B)$，则 $B\mathbf{w} = \mathbf{0}$，因此
    $$
    \|A - B\|_2^2 \ge \|(A-B)\mathbf{w}\|^2 = \|A\mathbf{w}\|^2.
    $$

    写 $\mathbf{w} = \sum_{i=1}^{k+1} c_i \mathbf{v}_i$（$\sum c_i^2 = 1$），则
    $$
    \|A\mathbf{w}\|^2 = \left\|\sum_{i=1}^{k+1} c_i \sigma_i \mathbf{u}_i\right\|^2 = \sum_{i=1}^{k+1} c_i^2 \sigma_i^2 \ge \sigma_{k+1}^2 \sum c_i^2 = \sigma_{k+1}^2.
    $$
    因此 $\|A - B\|_2 \ge \sigma_{k+1} = \|A - A_k\|_2$。$\blacksquare$

!!! example "例 11.7"
    设 $A = \begin{pmatrix} 3 & 2 & 2 \\ 2 & 3 & -2 \end{pmatrix}$，求 $A$ 的最佳秩一近似。

    **解：** 计算
    $$
    A A^T = \begin{pmatrix} 17 & 8 \\ 8 & 17 \end{pmatrix}.
    $$
    特征值 $\lambda_1 = 25$，$\lambda_2 = 9$，故 $\sigma_1 = 5$，$\sigma_2 = 3$。

    对 $\lambda_1 = 25$：$\mathbf{u}_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix}$。

    对应右奇异向量：
    $$
    \mathbf{v}_1 = \frac{1}{\sigma_1} A^T \mathbf{u}_1 = \frac{1}{5} \cdot \frac{1}{\sqrt{2}} \begin{pmatrix}5\\5\\0\end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\\0\end{pmatrix}.
    $$

    最佳秩一近似为
    $$
    A_1 = \sigma_1 \mathbf{u}_1 \mathbf{v}_1^T = 5 \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1\\1\end{pmatrix} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1&1&0\end{pmatrix} = \begin{pmatrix} \frac{5}{2} & \frac{5}{2} & 0 \\ \frac{5}{2} & \frac{5}{2} & 0 \end{pmatrix}.
    $$

    近似误差 $\|A - A_1\|_2 = \sigma_2 = 3$。

---

## 11.7 Moore-Penrose 伪逆

对于不可逆或非方阵的矩阵，Moore-Penrose 伪逆（pseudoinverse）提供了"最佳逆"的概念，可通过 SVD 简洁地构造。

!!! definition "定义 11.8 (Moore-Penrose 伪逆)"
    设 $A$ 为 $m \times n$ 实矩阵。$A$ 的 **Moore-Penrose 伪逆**（pseudoinverse）$A^+$ 是满足以下四个条件的唯一 $n \times m$ 矩阵：

    1. $A A^+ A = A$；
    2. $A^+ A A^+ = A^+$；
    3. $(A A^+)^T = A A^+$（即 $A A^+$ 是对称矩阵）；
    4. $(A^+ A)^T = A^+ A$（即 $A^+ A$ 是对称矩阵）。

!!! theorem "定理 11.8 (Moore-Penrose 伪逆的存在唯一性)"
    对任意 $m \times n$ 实矩阵 $A$，满足定义 11.8 中四个条件的矩阵 $A^+$ 存在且唯一。

??? proof "证明"
    **存在性：** 设 $A = U\Sigma V^T$ 为 $A$ 的 SVD，$\operatorname{rank}(A) = r$。定义
    $$
    \Sigma^+ = \begin{pmatrix} \sigma_1^{-1} & & \\ & \ddots & \\ & & \sigma_r^{-1} \\ & \mathbf{0} & \end{pmatrix}_{n \times m},
    $$
    令 $A^+ = V \Sigma^+ U^T$。下面验证四个条件。

    (1) $A A^+ A = (U\Sigma V^T)(V\Sigma^+ U^T)(U\Sigma V^T) = U\Sigma\Sigma^+\Sigma V^T = U\Sigma V^T = A$。

    这里 $\Sigma\Sigma^+\Sigma = \Sigma$ 成立，因为对角元素满足 $\sigma_i \cdot \sigma_i^{-1} \cdot \sigma_i = \sigma_i$（$i \le r$），其余为零。

    (2)-(4) 类似验证。

    **唯一性：** 设 $B_1, B_2$ 都满足四个条件。利用条件 (1)(3) 可证 $AB_1 = AB_2$，再由 (1)(4) 可证 $B_1A = B_2A$，最终由 (2) 得 $B_1 = B_1AB_1 = B_2AB_2 = B_2$。$\blacksquare$

!!! definition "定义 11.9 (伪逆的 SVD 公式)"
    设 $A = U\Sigma V^T$ 为 $A$ 的 SVD，$\operatorname{rank}(A) = r$，则
    $$
    A^+ = V \Sigma^+ U^T = \sum_{i=1}^{r} \frac{1}{\sigma_i} \mathbf{v}_i \mathbf{u}_i^T.
    $$

!!! theorem "定理 11.9 (伪逆与最小二乘)"
    设 $A$ 为 $m \times n$ 矩阵，$\mathbf{b} \in \mathbb{R}^m$。则 $\mathbf{x}^* = A^+ \mathbf{b}$ 是线性方程组 $A\mathbf{x} = \mathbf{b}$ 的**最小范数最小二乘解**，即在所有使 $\|A\mathbf{x} - \mathbf{b}\|$ 最小的 $\mathbf{x}$ 中，$\mathbf{x}^*$ 的范数 $\|\mathbf{x}^*\|$ 最小。

??? proof "证明"
    设 $A = U\Sigma V^T$，$\operatorname{rank}(A) = r$。令 $\mathbf{c} = U^T \mathbf{b}$，$\mathbf{y} = V^T \mathbf{x}$。由于 $U, V$ 是正交矩阵，
    $$
    \|A\mathbf{x} - \mathbf{b}\|^2 = \|\Sigma\mathbf{y} - \mathbf{c}\|^2 = \sum_{i=1}^r (\sigma_i y_i - c_i)^2 + \sum_{i=r+1}^m c_i^2.
    $$
    第二项与 $\mathbf{x}$ 无关。最小化第一项得 $y_i = c_i / \sigma_i$（$i = 1, \ldots, r$）。

    在此约束下，$\|\mathbf{x}\|^2 = \|\mathbf{y}\|^2 = \sum_{i=1}^r y_i^2 + \sum_{i=r+1}^n y_i^2$。最小化 $\|\mathbf{x}\|$ 需取 $y_i = 0$（$i = r+1, \ldots, n$）。

    因此最小范数最小二乘解为 $\mathbf{y}^* = (\frac{c_1}{\sigma_1}, \ldots, \frac{c_r}{\sigma_r}, 0, \ldots, 0)^T$，即 $\mathbf{x}^* = V\mathbf{y}^* = V\Sigma^+ U^T \mathbf{b} = A^+ \mathbf{b}$。$\blacksquare$

!!! example "例 11.8"
    设 $A = \begin{pmatrix}1\\2\end{pmatrix}$，求 $A^+$ 并求 $A\mathbf{x} = \begin{pmatrix}3\\4\end{pmatrix}$ 的最小二乘解。

    **解：** $A^T A = (5)$，$\sigma_1 = \sqrt{5}$。$\mathbf{v}_1 = (1)$，$\mathbf{u}_1 = \frac{1}{\sqrt{5}}\begin{pmatrix}1\\2\end{pmatrix}$。

    $$
    A^+ = \frac{1}{\sigma_1} \mathbf{v}_1 \mathbf{u}_1^T = \frac{1}{\sqrt{5}} \cdot (1) \cdot \frac{1}{\sqrt{5}}\begin{pmatrix}1&2\end{pmatrix} = \frac{1}{5}\begin{pmatrix}1&2\end{pmatrix} = \begin{pmatrix}\frac{1}{5}&\frac{2}{5}\end{pmatrix}.
    $$

    最小二乘解：
    $$
    x^* = A^+ \mathbf{b} = \frac{1}{5}(1 \cdot 3 + 2 \cdot 4) = \frac{11}{5} = 2.2.
    $$

    验证：$A x^* = \begin{pmatrix}2.2\\4.4\end{pmatrix}$，残差 $\mathbf{b} - Ax^* = \begin{pmatrix}0.8\\-0.4\end{pmatrix}$，与 $A$ 的列空间正交：$A^T(\mathbf{b} - Ax^*) = (1)(0.8) + (2)(-0.4) = 0$。

---

## 11.8 SVD 的应用

SVD 在众多领域有着重要应用，本节介绍几个典型场景。

### 11.8.1 数据压缩

!!! definition "定义 11.10 (SVD 图像压缩)"
    将 $m \times n$ 灰度图像视为矩阵 $A$，其秩 $k$ 截断 SVD 近似 $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$ 可用于图像压缩。存储 $A_k$ 需要 $k(m + n + 1)$ 个数，而原始矩阵需要 $mn$ 个数。当 $k \ll \min(m,n)$ 时，压缩比为
    $$
    \rho = \frac{k(m + n + 1)}{mn}.
    $$

!!! example "例 11.9"
    一幅 $1000 \times 800$ 的灰度图像，若取秩 $k = 50$ 的截断 SVD 近似，压缩比是多少？

    **解：**
    $$
    \rho = \frac{50 \times (1000 + 800 + 1)}{1000 \times 800} = \frac{50 \times 1801}{800000} = \frac{90050}{800000} \approx 11.3\%.
    $$
    即只需存储原数据量的约 $11.3\%$。

### 11.8.2 主成分分析简介

!!! theorem "定理 11.10 (SVD 与 PCA 的关系)"
    设数据矩阵 $X$（$n \times p$，已中心化）的 SVD 为 $X = U\Sigma V^T$。则：

    1. $X$ 的协方差矩阵为 $S = \frac{1}{n-1}X^TX = \frac{1}{n-1}V\Sigma^2 V^T$；
    2. $V$ 的列即为主成分方向（主轴）；
    3. 第 $i$ 个主成分的方差为 $\frac{\sigma_i^2}{n-1}$；
    4. 主成分得分矩阵为 $XV = U\Sigma$。

??? proof "证明"
    由 $X = U\Sigma V^T$，得 $X^T X = V\Sigma^T U^T U \Sigma V^T = V\Sigma^T\Sigma V^T = V\Sigma^2 V^T$（其中 $\Sigma^2$ 指 $\Sigma^T\Sigma$，是 $p \times p$ 对角矩阵）。因此
    $$
    S = \frac{1}{n-1}X^T X = V \left(\frac{\Sigma^2}{n-1}\right) V^T.
    $$
    这正是协方差矩阵的谱分解。$V$ 的列是 $S$ 的特征向量（主成分方向），对应特征值 $\frac{\sigma_i^2}{n-1}$（第 $i$ 个主成分的方差）。

    主成分得分 $Z = XV = U\Sigma V^T V = U\Sigma$。$\blacksquare$

### 11.8.3 条件数

!!! definition "定义 11.11 (矩阵条件数 Condition Number)"
    设 $A$ 为 $m \times n$ 矩阵，$\operatorname{rank}(A) = r$。$A$ 的**条件数**（condition number）定义为
    $$
    \kappa(A) = \frac{\sigma_1}{\sigma_r} = \frac{\sigma_{\max}}{\sigma_{\min}},
    $$
    其中 $\sigma_1$ 和 $\sigma_r$ 分别为 $A$ 的最大和最小非零奇异值。

!!! theorem "定理 11.11 (条件数与数值稳定性)"
    对于线性方程组 $A\mathbf{x} = \mathbf{b}$（$A$ 可逆），若 $\mathbf{b}$ 有相对扰动 $\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}$，则解的相对扰动满足
    $$
    \frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \kappa(A) \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.
    $$
    条件数越大，问题越**病态**（ill-conditioned）；条件数接近 1 时，问题**良态**（well-conditioned）。

??? proof "证明"
    设 $A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$，则 $A\delta\mathbf{x} = \delta\mathbf{b}$，故 $\delta\mathbf{x} = A^{-1}\delta\mathbf{b}$。

    $$
    \|\delta\mathbf{x}\| = \|A^{-1}\delta\mathbf{b}\| \le \|A^{-1}\|_2 \|\delta\mathbf{b}\|.
    $$

    又 $\|\mathbf{b}\| = \|A\mathbf{x}\| \le \|A\|_2 \|\mathbf{x}\|$，故 $\frac{1}{\|\mathbf{x}\|} \le \frac{\|A\|_2}{\|\mathbf{b}\|}$。因此
    $$
    \frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \le \|A\|_2 \|A^{-1}\|_2 \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = \frac{\sigma_1}{\sigma_n} \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|} = \kappa(A) \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}. \qquad \blacksquare
    $$

!!! example "例 11.10"
    设 $A = \begin{pmatrix}1&1\\1&1.0001\end{pmatrix}$，求条件数并分析数值稳定性。

    **解：**
    $$
    A^T A = \begin{pmatrix}2&2.0001\\2.0001&2.00020001\end{pmatrix}.
    $$
    特征值约为 $\lambda_1 \approx 4.00020001$ 和 $\lambda_2 \approx 0.000000005$（精确计算可得）。

    粗略估计：$\sigma_1 \approx 2.00005$，$\sigma_2 \approx 0.00005$。

    $$
    \kappa(A) = \frac{\sigma_1}{\sigma_2} \approx \frac{2.00005}{0.00005} \approx 40000.
    $$

    条件数约为 $40000$，说明该矩阵是病态的。右端 $\mathbf{b}$ 的微小扰动可能导致解的巨大变化。

---

## 本章小结

本章系统介绍了奇异值分解（SVD），包括：

1. **奇异值**来源于 $A^T A$ 的特征值的平方根；
2. **SVD 定理**保证任意矩阵 $A = U\Sigma V^T$；
3. **几何意义**：线性变换 = 旋转 + 缩放 + 旋转；
4. **Eckart-Young-Mirsky 定理**：截断 SVD 给出最优低秩近似；
5. **Moore-Penrose 伪逆** $A^+ = V\Sigma^+ U^T$ 给出最小范数最小二乘解；
6. **应用**涵盖数据压缩、PCA、条件数分析等领域。

SVD 将线性代数的理论优美性与实际计算的有效性完美统一，是后续学习数值线性代数、机器学习等课程的基石。
