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

正交性确立了线性空间的度量几何：

1.  **坐标的纯净化**：Gram-Schmidt 过程展示了如何从杂乱的基中提取出相互独立、标准化的“纯净”维度。
2.  **误差的最优解**：最小二乘法通过正交投影，为现实中因噪声而矛盾的观测数据找到了数学上最合理的“折中”方案。
3.  **稳定性保障**：正交矩阵 $Q$ 因其保持能量（范数）的特性，成为了现代数值线性代数算法（如 Householder 变换、QR 算法）的底层核心。
