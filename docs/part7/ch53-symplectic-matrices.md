# 第 53 章 辛矩阵与 Hamilton 矩阵

<div class="context-flow" markdown>

**前置**：对偶空间 (Ch13A) · 矩阵群 (Ch55) · 正定矩阵 (Ch16)

**本章脉络**：辛形式 (Symplectic Form) 与标准反对称矩阵 $J$ $\to$ 辛矩阵定义与性质 $\to$ 辛群 $Sp(2n, \mathbb{R})$ $\to$ 谱对称性（特征值成对出现） $\to$ Hamilton 矩阵定义与性质 $\to$ 指数映射与 Cayley 变换 $\to$ Williamson 定理（正定阵的辛对角化） $\to$ 应用：Hamilton 力学（相空间体积守恒）、辛几何积分器、量子高斯态

**延伸**：辛矩阵是描述物理系统“哈密顿演化”的代数语言；它保持了相空间的面积（或体积）不变，是研究保守系统、机器人动力学和量子计算中不可或缺的数学工具

</div>

在正交变换保持“欧氏距离”的同时，有一类矩阵保持着一种更微妙的“面积结构”。这类矩阵被称为**辛矩阵**（Symplectic Matrices）。它们与物理学中的哈密顿力学息息相关，揭示了动力系统相空间演化的内在约束。本章将从辛形式的定义出发，确立辛矩阵与 Hamilton 矩阵的对偶关系，并介绍其在现代数值计算中的核心作用。

---

## 53.1 辛形式与标准算子 $J$

!!! definition "定义 53.1 (标准辛算子 $J$)"
    在 $2n$ 维空间中，定义标准反对称矩阵 $J$：
    $$J = \begin{pmatrix} 0 & I_n \\ -I_n & 0 \end{pmatrix}$$
    满足 $J^2 = -I$ 且 $J^T = -J = J^{-1}$。

!!! definition "定义 53.2 (辛矩阵)"
    矩阵 $M \in M_{2n}(\mathbb{R})$ 称为**辛矩阵**，如果满足：
    $$M^T J M = J$$
    这意味着 $M$ 保持辛双线性型 $\omega(u, v) = u^T J v$ 不变。

---

## 53.2 辛群及其谱性质

!!! theorem "定理 53.1 (辛矩阵的性质)"
    1.  **行列式**：所有辛矩阵的行列式均为 $\det(M) = 1$。
    2.  **谱对称性**：若 $\lambda$ 是辛矩阵的特征值，则 $1/\lambda, \bar{\lambda}, 1/\bar{\lambda}$ 也必然是特征值。
    3.  **群结构**：所有 $2n \times 2n$ 辛矩阵在乘法下构成**辛群** $Sp(2n, \mathbb{R})$。

---

## 53.3 Hamilton 矩阵

!!! definition "定义 53.3 (Hamilton 矩阵)"
    矩阵 $H$ 称为 **Hamilton 矩阵**，如果满足：
    $$(JH)^T = JH$$
    这等价于 $H^T J + JH = 0$。
    **关系**：Hamilton 矩阵构成了辛群的李代数 $\mathfrak{sp}(2n)$。

---

## 53.4 Williamson 定理

!!! theorem "定理 53.2 (Williamson 定理)"
    对于任意 $2n \times 2n$ 实正定矩阵 $A \succ 0$，存在辛矩阵 $M$ 使得：
    $$M^T A M = \operatorname{diag}(d_1, \ldots, d_n, d_1, \ldots, d_n)$$
    其中 $d_i > 0$ 称为 $A$ 的**辛特征值**（Symplectic eigenvalues）。

---

## 练习题

1. **[计算] 验证 $J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 本身是否是辛矩阵。**
   ??? success "参考答案"
       $J^T J J = (-J)J^2 = (-J)(-I) = J$。成立。

2. **[行列式] 证明：若 $M$ 是辛矩阵，则 $M$ 是非奇异的。**
   ??? success "参考答案"
       由 $M^T J M = J$，取行列式得 $\det(M)^2 \det(J) = \det(J)$。由于 $\det(J)=1 \neq 0$，故 $\det(M)^2=1 \implies \det(M) = \pm 1$。

3. **[特征值] 解释为什么辛矩阵的特征值不能仅包含一个 2。**
   ??? success "参考答案"
       由于特征值成对出现，若有 2，则必有 $1/2 = 0.5$。

4. **[Hamilton] 判定 $\begin{pmatrix} A & B \\ C & -A^T \end{pmatrix}$ 何时是 Hamilton 矩阵（其中 $B, C$ 对称）。**
   ??? success "参考答案"
       始终是。代入 $H^T J + JH = 0$ 验证，或观察 $JH = \begin{pmatrix} C & -A^T \\ -A & -B \end{pmatrix}$ 的对称性。

5. **[指数映射] 证明：若 $H$ 是 Hamilton 矩阵，则 $e^{H}$ 是辛矩阵。**
   ??? success "参考答案"
       $(e^H)^T J e^H = e^{H^T} J e^H = J e^{-H} e^H = J$。利用了 $H^T J = -JH \implies H^T = -JHJ^{-1}$。

6. **[Williamson] 求 $A = \operatorname{diag}(2, 2)$ 的辛特征值。**
   ??? success "参考答案"
       由于它已经是 Williamson 标准形形式，辛特征值为 $d_1 = 2$。

7. **[性质] 证明辛矩阵的逆也是辛矩阵。**
   ??? success "参考答案"
       对 $M^T J M = J$ 两边取逆并整理得 $(M^{-1})^T J M^{-1} = J$。

8. **[应用] 什么是辛积分器（Symplectic Integrator）？**
    ??? success "参考答案"
        一种数值求解 Hamilton 方程的算法，其每一步迭代矩阵都是辛矩阵，从而保证数值解长时间保持系统的能量和相空间体积不变。

9. **[关系] 辛矩阵与酉矩阵的交集是什么？**
    ??? success "参考答案"
        交集是同构于 $U(n)$ 的子群。

10. **[量子应用] 在量子光学中，辛变换对应于什么物理过程？**
    ??? success "参考答案"
        对应于高斯态的线性正则变换（如压缩、位移和波束分束），保持正则对易关系不变。

## 本章小结

辛矩阵与 Hamilton 矩阵确立了保守系统的代数骨架：

1.  **几何的不变量**：辛矩阵不仅是代数对象，更是保持微分几何中“辛面积”不变的算子，揭示了物理演化中的刚性约束。
2.  **谱的对称美**：特征值的四重对称性反映了 Hamilton 系统中能量守恒与时间反演对称性的深刻联系。
3.  **计算的守护者**：通过 Williamson 定理和辛算法，线性代数为复杂力学系统的数值模拟提供了“结构保持”的终极保证，防止了数值耗散引发的物理错误。
