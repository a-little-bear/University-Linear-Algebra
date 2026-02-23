# 第 12 章 Jordan 标准形

<div class="context-flow" markdown>

**前置**：特征值 (Ch06) · 矩阵分解 (Ch10) · 多项式代数 (Ch00)

**本章脉络**：对角化的局限性（亏损矩阵） $\to$ Jordan 块的定义 $\to$ 广义特征向量与 Jordan 链 $\to$ Jordan 标准形 (JCF) 的存在性与唯一性 $\to$ 最小多项式与特征多项式的关系 $\to$ 确定 JCF 的步骤（秩法、Weyr 特征） $\to$ 矩阵幂与级数的 Jordan 分析 $\to$ 数值不稳定性

**延伸**：Jordan 标准形是相似变换下的终极标准形；它在解决线性微分方程组 (Ch26) 和矩阵函数 (Ch13) 的理论分析中不可或缺

</div>

并非所有的方阵都能对角化。当一个特征值的几何重数小于代数重数时，矩阵被称为“亏损”的。Jordan 标准形（Jordan Canonical Form, JCF）为这类矩阵提供了最接近对角形的结构。它不仅揭示了线性算子的深层结构，也是矩阵分析论中最重要的理论工具。

---

## 12.1 Jordan 块与广义特征向量

!!! definition "定义 12.1 (Jordan 块)"
    一个 $k$ 阶 **Jordan 块** $J_k(\lambda)$ 是一个对角线上全是 $\lambda$，紧邻对角线上方全为 1，其余为 0 的方阵：
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & 0 & \cdots & 0 \\ 0 & \lambda & 1 & \cdots & 0 \\ \vdots & \vdots & \ddots & \ddots & \vdots \\ 0 & 0 & 0 & \lambda & 1 \\ 0 & 0 & 0 & 0 & \lambda \end{pmatrix}$$

!!! definition "定义 12.2 (广义特征向量)"
    若向量 $\mathbf{v}$ 满足 $(A - \lambda I)^k \mathbf{v} = \mathbf{0}$ 但 $(A - \lambda I)^{k-1} \mathbf{v} \neq \mathbf{0}$，则称其为对应于 $\lambda$ 的 **$k$ 阶广义特征向量**。

---

## 12.2 Jordan 标准形定理

!!! theorem "定理 12.1 (Jordan 标准形定理)"
    每一个复方阵 $A$ 都相似于一个 Jordan 标准形 $J$，且 $J$ 在不计块的排列顺序下是唯一的：
    $$P^{-1} A P = J = \operatorname{diag}(J_{k_1}(\lambda_1), J_{k_2}(\lambda_2), \ldots, J_{k_m}(\lambda_m))$$
    - 每个 Jordan 块对应一个特征向量（几何重数）。
    - 属于同一特征值的 Jordan 块的总阶数等于其代数重数。

---

## 12.3 最小多项式

!!! definition "定义 12.3 (最小多项式)"
    使得 $m(A) = O$ 的次数最低的首一多项式 $m(\lambda)$ 称为 $A$ 的**最小多项式**。
    **性质**：
    1.  $m(\lambda)$ 整除特征多项式 $p(\lambda)$。
    2.  $A$ 可对角化 $\iff$ $m(\lambda)$ 没有重根。
    3.  特征值 $\lambda_i$ 在 $m(\lambda)$ 中的重数等于对应的最大 Jordan 块的阶数。

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

Jordan 标准形是方阵结构的最终判决：

1.  **亏损的补完**：它通过引入“1”步进结构（Jordan 块），完美填补了非对角化矩阵在特征向量个数上的缺失。
2.  **多项式的深度**：最小多项式与 JCF 块大小的对应关系，揭示了矩阵作为多项式根时的几何深度。
3.  **结构的唯一性**：JCF 确立了方阵相似类的分类标准，是分析矩阵函数、幂序列收敛性最精确的理论框架。
