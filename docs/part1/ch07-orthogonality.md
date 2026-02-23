# 第 07 章 正交性与最小二乘

<div class="context-flow" markdown>

**前置**：向量空间 (Ch04) · 线性变换 (Ch05)

**本章脉络**：内积与范数 $\to$ 正交性定义 $\to$ 正交向量组与标准正交基 $\to$ 正交投影矩阵 $\to$ Gram-Schmidt 正交化过程 $\to$ QR 分解 $\to$ 最小二乘法 (Least Squares) $\to$ 正规方程 $\to$ 正交矩阵的性质 $\to$ 最佳逼近定理

**延伸**：正交性为抽象的向量空间引入了度量（长度和角度），是信号处理（小波、傅里叶）、数值计算稳定性及 SVD (Ch11) 的基石

</div>

如果线性无关是空间的“骨架”，那么正交性就是它的“标尺”。正交性不仅简化了坐标表示，更通过投影机制解决了现实中由于数据噪声导致的“无解”方程组问题。本章将展示如何通过正交化过程将任意一组基转化为最完美的计算基。

---

## 07.1 内积、范数与正交性

!!! definition "定义 07.1 (内积与长度)"
    对于 $\mathbb{R}^n$ 中的向量 $\mathbf{u}, \mathbf{v}$，其**点积**（内积）定义为 $\mathbf{u} \cdot \mathbf{v} = \mathbf{u}^T \mathbf{v}$。向量的**长度**（范数）定义为 $\|\mathbf{v}\| = \sqrt{\mathbf{v} \cdot \mathbf{v}}$。

!!! definition "定义 07.2 (正交性)"
    若 $\mathbf{u} \cdot \mathbf{v} = 0$，则称向量 $\mathbf{u}, \mathbf{v}$ **正交**。
    **勾股定理**：$\mathbf{u}, \mathbf{v}$ 正交 $\iff$ $\|\mathbf{u} + \mathbf{v}\|^2 = \|\mathbf{u}\|^2 + \|\mathbf{v}\|^2$。

---

## 07.2 正交投影与 Gram-Schmidt

!!! theorem "定理 07.1 (最佳逼近定理)"
    设 $W$ 是 $V$ 的子空间。对于 $V$ 中的向量 $\mathbf{y}$，其在 $W$ 上的正交投影 $\hat{\mathbf{y}} = \operatorname{proj}_W \mathbf{y}$ 是 $W$ 中距离 $\mathbf{y}$ 最近的向量。

!!! algorithm "算法 07.1 (Gram-Schmidt 过程)"
    给定一组基 $\{\mathbf{x}_1, \ldots, \mathbf{x}_n\}$，构造标准正交基 $\{\mathbf{q}_1, \ldots, \mathbf{q}_n\}$：
    1.  $\mathbf{v}_1 = \mathbf{x}_1$
    2.  $\mathbf{v}_2 = \mathbf{x}_2 - \frac{\mathbf{x}_2 \cdot \mathbf{v}_1}{\mathbf{v}_1 \cdot \mathbf{v}_1} \mathbf{v}_1$
    3.  以此类推，最后单位化：$\mathbf{q}_i = \mathbf{v}_i / \|\mathbf{v}_i\|$。

---

## 07.3 QR 分解

!!! theorem "定理 07.2 (QR 分解)"
    若 $m \times n$ 矩阵 $A$ 的列线性无关，则 $A$ 可以分解为 $A = QR$。
    - $Q$ 是 $m \times n$ 矩阵，其列向量为 $A$ 的列空间的标准正交基。
    - $R$ 是 $n \times n$ 上三角可逆矩阵。

---

## 07.4 最小二乘法

!!! definition "定义 07.3 (最小二乘解)"
    对于方程 $Ax = \mathbf{b}$，若无解，我们寻找 $\hat{x}$ 使得 $\|\mathbf{b} - A\hat{x}\|$ 达到最小。
    $\hat{x}$ 满足**正规方程**（Normal Equations）：
    $$A^T A \hat{x} = A^T \mathbf{b}$$

---

## 练习题

**1. [内积] 计算 $(1, 2, 3)$ 与 $(1, 0, -1)$ 是否正交。**

??? success "参考答案"
    **计算点积：**
    $(1)(1) + (2)(0) + (3)(-1) = 1 + 0 - 3 = -2$。
    **结论**：由于内积不为 0，这两个向量**不正交**。它们之间的夹角为钝角。

**2. [单位化] 将向量 $(3, 4)$ 单位化。**

??? success "参考答案"
    **步骤 1：计算长度（范数）。**
    $\|\mathbf{v}\| = \sqrt{3^2 + 4^2} = \sqrt{9 + 16} = 5$。
    **步骤 2：缩放向量。**
    单位向量 $\mathbf{u} = \frac{1}{\|\mathbf{v}\|} \mathbf{v} = \frac{1}{5} (3, 4) = (0.6, 0.8)$。
    **结论**：单位化后的向量方向不变，长度为 1。

**3. [投影] 求 $(1, 2)$ 在 $(1, 0)$ 方向上的投影。**

??? success "参考答案"
    **利用投影公式：**
    $\hat{y} = \frac{\mathbf{y} \cdot \mathbf{u}}{\mathbf{u} \cdot \mathbf{u}} \mathbf{u}$。
    1. 分子：$(1, 2) \cdot (1, 0) = 1$。
    2. 分母：$(1, 0) \cdot (1, 0) = 1$。
    **结论**：$\hat{y} = 1 \cdot (1, 0) = (1, 0)$。几何上，这对应于将点垂直投影到 $x$ 轴上。

**4. [正交阵] 证明正交矩阵 $Q$ 保持长度不变，即 $\|Q\mathbf{x}\| = \|\mathbf{x}\|$。**

??? success "参考答案"
    **证明过程：**
    1. 考虑长度的平方：$\|Q\mathbf{x}\|^2 = (Q\mathbf{x})^T(Q\mathbf{x})$。
    2. 利用转置性质：$= \mathbf{x}^T Q^T Q \mathbf{x}$。
    3. 利用正交阵定义 $Q^T Q = I$：$= \mathbf{x}^T I \mathbf{x}$。
    4. $= \mathbf{x}^T \mathbf{x} = \|\mathbf{x}\|^2$。
    **结论**：正交变换是一种**等距变换**（旋转或反射），不会拉伸空间。

**5. [G-S过程] 对 $(1, 1), (0, 1)$ 进行正交化。**

??? success "参考答案"
    **步骤 1：确定第一个基向量。**
    $\mathbf{v}_1 = (1, 1)$。
    **步骤 2：求第二个正交分量。**
    $\mathbf{v}_2 = (0, 1) - \operatorname{proj}_{\mathbf{v}_1}(0, 1) = (0, 1) - \frac{(0,1)\cdot(1,1)}{(1,1)\cdot(1,1)}(1, 1) = (0, 1) - \frac{1}{2}(1, 1) = (-0.5, 0.5)$。
    **步骤 3：单位化。**
    $\mathbf{q}_1 = \frac{1}{\sqrt{2}}(1, 1), \mathbf{q}_2 = \frac{1}{\sqrt{0.5}}( -0.5, 0.5) = \frac{1}{\sqrt{2}}(-1, 1)$。

**6. [最小二乘] 为什么 $A^T A$ 在 $A$ 列满秩时是可逆的？**

??? success "参考答案"
    **证明思路：**
    1. 设 $A^T A \mathbf{x} = 0$。
    2. 两边左乘 $\mathbf{x}^T$：$\mathbf{x}^T A^T A \mathbf{x} = (A\mathbf{x})^T(A\mathbf{x}) = \|A\mathbf{x}\|^2 = 0$。
    3. 由范数性质知 $A\mathbf{x} = \mathbf{0}$。
    4. 由于 $A$ 是列满秩的，其列向量线性无关，故其零空间只有零向量，即 $\mathbf{x} = \mathbf{0}$。
    **结论**：这意味着 $A^T A$ 的核只有零向量，因此是非奇异（可逆）的。

**7. [QR] 若 $A=QR$，则 $A^T A$ 等于什么？**

??? success "参考答案"
    **代数推导：**
    1. $A^T A = (QR)^T (QR)$。
    2. $= R^T Q^T Q R$。
    3. 由于 $Q$ 的列是标准正交的，$Q^T Q = I$。
    **结论**：$A^T A = R^T R$。这是对称正定矩阵的 Cholesky 分解的变体。

**8. [性质] 证明正交矩阵的行列式必为 $\pm 1$。**

??? success "参考答案"
    **证明：**
    1. 已知 $Q^T Q = I$。
    2. 取行列式：$\det(Q^T Q) = \det(I) = 1$。
    3. 利用乘法性质：$\det(Q^T)\det(Q) = 1$。
    4. 由于 $\det(Q^T) = \det(Q)$，得 $[\det(Q)]^2 = 1$。
    **结论**：$\det(Q) = 1$（旋转）或 $-1$（包含反射的旋转）。

**9. [几何] 最小二乘解 $\hat{x}$ 对应的残差向量 $\mathbf{b} - A\hat{x}$ 与 $A$ 的列空间是什么关系？**

??? success "参考答案"
    **结论：**
    残差向量与 $A$ 的列空间 $C(A)$ **正交**。
    **理由**：根据正规方程 $A^T(A\hat{x} - \mathbf{b}) = 0$，这意味着 $A$ 的每一行（即 $A^T$ 的列）与残差向量的点积均为 0。因此残差向量垂直于 $A$ 的所有列向量所张成的空间。

**10. [应用] 为什么在数值计算中 QR 分解比正规方程求解最小二乘更稳定？**

??? success "参考答案"
    **稳定性分析：**
    1. 直接计算 $A^T A$ 会导致**条件数平方倍增加**：$\kappa(A^T A) = \kappa(A)^2$。这会极大放大舍入误差。
    2. QR 分解使用等距变换（Householder 或 Givens），正交矩阵的条件数始终为 1，不会改变问题的病态程度。
    **结论**：QR 分解保持了数据的数值精度，是工业级最小二乘求解器的标准方法。
