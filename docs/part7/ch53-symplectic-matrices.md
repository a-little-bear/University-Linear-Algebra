# 第 53 章 辛矩阵与辛几何

<div class="context-flow" markdown>

**前置**：矩阵代数(Ch2) · 正交性(Ch7) · 特征多项式(Ch6) · 二次型(Ch9)

**本章脉络**：辛形式 $J$ → 辛群 $Sp(2n, \mathbb{R})$ → 辛矩阵 → 特征值性质 → Williamson 定理 → Hamilton 矩阵 → 辛基底 → 在经典力学中的应用

**延伸**：辛矩阵是 Hamiltonian 力学（保持相空间体积）和线性光学系统的数学基础

</div>

辛矩阵保持一个非退化的交替双线性形式（通常由矩阵 $J = \begin{pmatrix} 0 & I \\ -I & 0 \end{pmatrix}$ 表示）。正交矩阵保持欧几里得内积（距离），而辛矩阵保持相空间中的“辛面积”。这种保持性质是经典力学运动定律的基础。

---

## 53.1 辛群与辛形式

!!! definition "定义 53.1 (辛矩阵)"
    若 $2n \times 2n$ 实矩阵 $M$ 满足：
    $$M^T J M = J, \quad \text{其中 } J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$$
    则称 $M$ 为**辛矩阵**。所有此类矩阵构成的集合称为**辛群** $Sp(2n, \mathbb{R})$。

!!! theorem "定理 53.1 (特征值对称性)"
    若 $\lambda$ 是辛矩阵 $M$ 的特征值，则 $1/\lambda$, $\bar{\lambda}$ 和 $1/\bar{\lambda}$ 也是具有相同重数的特征值。

---

## 练习题

1. **[基础] 证明每个辛矩阵 $M$ 的行列式 $\det M = 1$。**
   ??? success "参考答案"
       由 $M^T J M = J$ 取行列式得 $(\det M)^2 \det J = \det J$。由于 $J$ 非奇异，故 $(\det M)^2 = 1$。更深入的论证（如利用 Pfaffian 性质）可以证明 $\det M$ 必须为 $+1$ 而非 $-1$，因为 $Sp(2n, \mathbb{R})$ 是连通群且包含单位阵。

2. **[维数] 计算 Lie 代数 $\mathfrak{sp}(2n, \mathbb{R})$ 的维数。**
   ??? success "参考答案"
       该 Lie 代数由满足 $AJ + JA^T = 0$ 的矩阵 $A$ 组成。这一条件意味着 $JA$ 是对称矩阵。由于 $2n \times 2n$ 对称矩阵空间的维数为 $\frac{2n(2n+1)}{2}$，故 $\mathfrak{sp}(2n, \mathbb{R})$ 的维数为 $n(2n+1)$。

3. **[Hamilton矩阵] 定义 Hamilton 矩阵及其与辛矩阵的关系。**
   ??? success "参考答案"
       若 $JA$ 是对称矩阵，则称 $A$ 为 Hamilton 矩阵。这类矩阵是辛群的无穷小生成元：若 $A$ 是 Hamilton 矩阵，则矩阵指数 $e^{tA}$ 对于所有 $t \in \mathbb{R}$ 都是辛矩阵。

4. **[谱结构] 若 $\lambda = 2+i$ 是一个特征值，分析 $M$ 的整个特征值集合。**
   ??? success "参考答案"
       该谱必须包含四元组 $\{2+i, 2-i, \frac{2-i}{5}, \frac{2+i}{5}\}$。辛矩阵的特征值总是以四元组形式出现（除非它们落在单位圆或实轴上）。这一对称性确保了特征值之积为 1，与 $\det M = 1$ 一致。

5. **[Williamson定理] 阐述 Williamson 定理对正定矩阵的意义。**
   ??? success "参考答案"
       任何 $2n \times 2n$ 的正定矩阵 $V$ 都可以通过辛合同变换 $M^T V M$ 化为对角形 $\operatorname{diag}(d_1, \dots, d_n, d_1, \dots, d_n)$。其中的 $d_i > 0$ 被称为 $V$ 的**辛特征值**。

6. **[复合] 证明两个辛矩阵的乘积仍为辛矩阵。**
   ??? success "参考答案"
       设 $M_1, M_2 \in Sp(2n)$。则 $(M_1 M_2)^T J (M_1 M_2) = M_2^T (M_1^T J M_1) M_2 = M_2^T J M_2 = J$。故辛性质在矩阵乘法下保持。

7. **[正交辛] 证明一个矩阵既是正交的又是辛的，当且仅当它具有形式 $\begin{pmatrix} A & B \\ -B & A \end{pmatrix}$，其中 $A+iB$ 是复酉矩阵。**
   ??? success "参考答案"
       这源于同构关系 $Sp(2n, \mathbb{R}) \cap O(2n) \cong U(n)$。它表明酉变换是唯一既保持欧几里得度规又保持辛结构的线性变换。

8. **[求逆] 证明若 $M$ 是辛矩阵，其逆矩阵 $M^{-1}$ 也是辛矩阵。**
   ??? success "参考答案"
       $M^T J M = J \implies J^{-1} (M^T)^{-1} J M^{-1} = I \implies J (M^{-1})^T (-J) M^{-1} = I$。整理可得 $(M^{-1})^T J M^{-1} = J$，证明 $M^{-1}$ 满足辛矩阵定义。

9. **[相空间] 在经典力学中，为什么 Hamiltonian 流的 Jacobian 矩阵总是辛矩阵？**
   ??? success "参考答案"
       Hamilton 方程 $\dot{z} = J \nabla H(z)$ 描述了一个保持辛形式 $\omega = \sum dp_i \wedge dq_i$ 的流。根据 Liouville 定理，该流的线性化（即 Jacobian）必须保持 $\omega$，这等价于其为辛矩阵。

10. **[光学] 辛矩阵在线性光学系统中有何应用？**
    ??? success "参考答案"
        在射线光学中，ABCD 矩阵描述了位置和角度 $(x, \theta)$ 的变换。光展量（Etendue）守恒要求行列式为 1，而对于高维系统，传递矩阵必须是辛矩阵，以保持光束在相空间中的体积。

## 本章小结

本章探讨了保持非欧几里得相空间几何的矩阵：

1. **定义对称性**：通过保持标准反对称形式 $J$ 定义了辛矩阵。
2. **谱四元组**：分析了 $Sp(2n)$ 中特征值特有的倒数与共轭对称性。
3. **Hamilton 关联**：将辛群与其生成元（物理中使用的 Hamilton 矩阵）联系起来。
4. **辛规范化**：利用 Williamson 定理将谱分析扩展到了辛变换下的正定矩阵。
