# 第 12 章 Jordan 标准形

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · 矩阵分解 (Ch10) · 多项式代数 (Ch00) · 空间分解专题 (Ch13)

**本章脉络**：对角化的局限性（亏损矩阵） $\to$ Jordan 块的定义 $\to$ 广义特征向量与 Jordan 链 $\to$ 空间的根子空间分解 $\to$ Jordan 标准形 (JCF) 的存在性与唯一性 $\to$ 最小多项式与 JCF 的深刻联系 $\to$ 确定 JCF 的步骤（秩法、Weyr 特征） $\to$ 矩阵幂与级数的 Jordan 分析 $\to$ 数值不稳定性

**延伸**：Jordan 标准形是相似变换下的终极标准形；它完美揭示了线性算子在不可对角化时的“准对角”结构，是解决线性微分方程组 (Ch26) 和矩阵函数 (Ch13) 理论推导的必然路径

</div>

并非所有的方阵都能对角化。当一个特征值的几何重数小于代数重数时，矩阵被称为“亏损”的。**Jordan 标准形**（Jordan Canonical Form, JCF）为这类矩阵提供了最接近对角形的结构。它通过引入“1”步进结构（Jordan 块），将空间的退化程度量化为代数链的长度。本章将深入方阵结构的最终判决——JCF。

---

## 12.1 Jordan 块与广义特征向量

!!! definition "定义 12.1 (Jordan 块)"
    一个 $k$ 阶 **Jordan 块** $J_k(\lambda)$ 是一个对角线上全是 $\lambda$，紧邻对角线上方全为 1，其余为 0 的方阵：
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & 0 & \lambda & 1 \\ 0 & 0 & 0 & 0 & \lambda \end{pmatrix}$$

!!! definition "定义 12.2 (广义特征向量与 Jordan 链)"
    若向量 $\mathbf{v}_k$ 满足 $(A - \lambda I)^k \mathbf{v}_k = \mathbf{0}$ 但 $(A - \lambda I)^{k-1} \mathbf{v}_k \neq \mathbf{0}$，则称其为 **$k$ 阶广义特征向量**。
    序列 $\{ \mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k \}$ 构成一条 **Jordan 链**，其中 $\mathbf{v}_1$ 是普通特征向量。

---

## 12.2 Jordan 标准形定理

!!! theorem "定理 12.1 (JCF 存在唯一性)"
    每一个复方阵 $A$ 都相似于一个 Jordan 标准形 $J$，且 $J$ 在不计块的排列顺序下是**唯一**的：
    $$P^{-1} A P = J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_m}(\lambda_m))$$
    - 每个 Jordan 块对应一个线性无关的特征向量（几何重数）。
    - 属于同一特征值的 Jordan 块的总阶数等于其代数重数。

---

## 12.3 最小多项式与 JCF

!!! theorem "定理 12.2 (最小多项式判据)"
    使得 $m(A) = O$ 的最低次首一多项式 $m(\lambda)$ 称为 $A$ 的**最小多项式**。
    - **对角化判定**：$A$ 可对角化 $\iff$ $m(\lambda)$ 没有重根。
    - **块大小判定**：特征值 $\lambda_i$ 在 $m(\lambda)$ 中的重数，等于 $A$ 的 JCF 中属于 $\lambda_i$ 的**最大 Jordan 块的阶数**。

---

## 练习题

**1. [Jordan块] 写出 $J_2(5)$ 的平方及其特征值。**

??? success "参考答案"
    **计算步骤：**
    1. $J_2(5) = \begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix}$。
    2. 矩阵乘法：$\begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix} \begin{pmatrix} 5 & 1 \\ 0 & 5 \end{pmatrix} = \begin{pmatrix} 25 & 10 \\ 0 & 25 \end{pmatrix}$。
    **特征值分析**：
    由于 $J_2(5)^2$ 是上三角阵，特征值直接从对角线读取，为 $25$（代数重数为 2）。注意结果仍然是一个 Jordan 块的变形（非对角阵）。

**2. [对角化] 判定 $A = \begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ 是否可对角化。**

??? success "参考答案"
    **判定逻辑：**
    1. 特征值为 2，代数重数 $\alpha = 2$。
    2. 计算特征空间 $E_2$ 的维数（几何重数）：
       $A - 2I = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。
       其秩为 1，故零空间维数 $\gamma = 2 - 1 = 1$。
    3. 由于 $\gamma < \alpha$（特征向量不足），矩阵**不可对角化**。它本身就是一个 2 阶 Jordan 块。

**3. [JCF判定] 若 $3 \times 3$ 矩阵 $A$ 的特征值均为 0，且 $\operatorname{rank}(A)=1$，求其 JCF。**

??? success "参考答案"
    **分析过程：**
    1. 代数重数总和为 3。
    2. 几何重数 $\gamma = n - \operatorname{rank}(A) = 3 - 1 = 2$。
    3. 几何重数 2 意味着 JCF 中总共有 **2 个 Jordan 块**。
    4. 将数字 3 拆分为 2 个正整数之和：只能是 $2 + 1$。
    **结论**：JCF 为 $\operatorname{diag}(J_2(0), J_1(0)) = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$。

**4. [最小多项式] 已知 $J = \operatorname{diag}(J_3(2), J_2(2))$，求其特征多项式与最小多项式。**

??? success "参考答案"
    **计算：**
    1. **特征多项式** $p(\lambda)$：所有块阶数之和。$(\lambda - 2)^{3+2} = (\lambda - 2)^5$。
    2. **最小多项式** $m(\lambda)$：最大块的阶数。对应特征值 2 的最大块阶数为 3。
    **结论**：$m(\lambda) = (\lambda - 2)^3$。

**5. [性质] 若 $A$ 的特征多项式为 $(\lambda-1)^4$，最小多项式为 $(\lambda-1)^2$，求 $A$ 的 JCF 所有可能形式。**

??? success "参考答案"
    **约束条件：**
    1. 块的总阶数为 4。
    2. 最大块的阶数必须为 2。
    **可能组合：**
    - 方案 A：$2 + 2$。即 $\operatorname{diag}(J_2(1), J_2(1))$。
    - 方案 B：$2 + 1 + 1$。即 $\operatorname{diag}(J_2(1), J_1(1), J_1(1))$。
    **结论**：共有两种非相似的可能结构。

**6. [幂零性] 描述 $J_k(0)^n$ 在 $n \ge k$ 时的结果。**

??? success "参考答案"
    **结论：**
    结果为**零矩阵**。
    **理由**：$J_k(0)$ 是严格上三角阵，每次乘方都会将对角线上方的 1 向右上角推移一格。经过 $k$ 次推移后，所有元素都会移出矩阵边界。这证明了特征值为 0 的 Jordan 块是幂零的。

**7. [秩法] 如何通过秩计算确定特征值 $\lambda$ 对应的阶数为 1 的块的个数？**

??? success "参考答案"
    **公式推导：**
    设 $n_k$ 为 $k$ 阶块的个数。根据秩的性质：
    $n_1 = \operatorname{rank}(A-\lambda I)^2 - 2\operatorname{rank}(A-\lambda I) + \operatorname{rank}(A-\lambda I)^0$
    更一般地，可以通过计算相邻幂次的秩差之差（二阶差分）来确定每一阶块的具体数量。

**8. [特征空间] 广义特征空间与普通特征空间有什么区别？**

??? success "参考答案"
    **对比：**
    - **特征空间**：被算子作用后仅发生缩放的向量集合（$(A-\lambda I)v=0$）。
    - **广义特征空间**：被算子多次作用后最终消失的向量集合（$(A-\lambda I)^k v=0$）。
    在 JCF 理论中，广义特征空间包含了整个 Jordan 链，从而填补了对角化时缺失的维度。

**9. [唯一性] 若 $A$ 与 $B$ 的 JCF 相同，它们是否相似？**

??? success "参考答案"
    **结论：**
    **是的**。
    **理由**：JCF 是相似变换下的全系不变量。如果两个矩阵相似于同一个标准形，那么根据等价关系的传递性，它们彼此相似。这意味着 JCF 提供了方阵相似类的完美分类。

**10. [数值] 为什么在科学计算软件中很少直接求解 JCF？**

??? success "参考答案"
    **稳定性分析：**
    1. JCF 对矩阵元素的扰动极度敏感。
    2. 例如 $\begin{pmatrix} 0 & 1 \\ \epsilon & 0 \end{pmatrix}$ 只要 $\epsilon \neq 0$，它就有两个互异特征值 $\pm\sqrt{\epsilon}$ 且可对角化。
    3. 但当 $\epsilon = 0$ 时，它突然变为一个 2 阶 Jordan 块。
    **结论**：这种非连续性导致 JCF 在带误差的浮点运算下极不稳定，通常被 **Schur 分解**或 **SVD** 取代。

## 本章小结

Jordan 标准形是方阵结构的最终判决：

1.  **亏损的补完**：它通过引入 Jordan 链，完美填补了非对角化矩阵在特征向量个数上的缺失，确立了空间的根子空间结构。
2.  **多项式的深度**：最小多项式与 JCF 块大小的对应关系，揭示了矩阵作为多项式根时的几何深度，是分析算子性质的利刃。
3.  **结构的唯一性**：JCF 确立了方阵相似类的分类标准，是分析矩阵函数、幂序列收敛性最精确的理论框架。
