# 第 55B 章 李代数

<div class="context-flow" markdown>

**前置**：矩阵群 (Ch55) · 矩阵分析 (Ch14) · 矩阵积与交换子 (Ch02)

**本章脉络**：李代数的抽象定义（Lie Bracket） $\to$ 双线性、反对称性与 Jacobi 恒等式 $\to$ 矩阵李代数作为矩阵群的切空间 $\to$ 经典群的李代数（$\mathfrak{gl}, \mathfrak{sl}, \mathfrak{so}, \mathfrak{su}, \mathfrak{sp}$） $\to$ 指数映射及其逆（对数映射） $\to$ 伴随表示 $\operatorname{ad}_X(Y) = [X, Y]$ $\to$ 结构常数 $\to$ 物理应用：量子力学中的对易关系、角动量算子 $\to$ 控制理论：非完整系统的可达性分析

**延伸**：李代数将“群”这种复杂的几何对象局部线性化为“向量空间”；它是研究连续对称性的终极代数工具，支撑着从广义相对论到规范场论的全部数学表述

</div>

如果说矩阵群描述的是“全局变换”，那么**李代数**（Lie Algebras）描述的则是“无穷小变换”。通过研究矩阵群在单位元处的切空间，我们获得了一个具有特殊乘法运算（李括号）的向量空间。李代数的优越性在于，它将复杂的群乘法转化为了线性的括号运算，极大地简化了对称性的分析。

---

## 55B.1 李代数的定义与公理

!!! definition "定义 55B.1 (李代数)"
    一个向量空间 $\mathfrak{g}$ 配备二元运算 $[\cdot, \cdot]: \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}$（称为**李括号**），若满足：
    1.  **双线性**：$[ax+by, z] = a[x, z] + b[y, z]$。
    2.  **反对称性**：$[x, y] = -[y, x]$（意味着 $[x, x] = 0$）。
    3.  **Jacobi 恒等式**：$[x, [y, z]] + [y, [z, x]] + [z, [x, y]] = 0$。

!!! note "矩阵李括号"
    对于方阵，$ [A, B] = AB - BA $（交换子）是满足上述公理的标准实现。

---

## 55B.2 经典矩阵群的李代数

!!! theorem "定理 55B.1 (经典李代数分类)"
    1.  **$\mathfrak{gl}(n)$**：一般线性群的李代数，为所有 $n \times n$ 矩阵。
    2.  **$\mathfrak{sl}(n)$**：特殊线性群的李代数，由**迹为 0** 的矩阵构成（$\operatorname{tr}(X) = 0$）。
    3.  **$\mathfrak{so}(n)$**：正交群的李代数，由**反对称**矩阵构成（$X^T = -X$）。
    4.  **$\mathfrak{su}(n)$**：特殊酉群的李代数，由**迹为 0 的反 Hermite** 矩阵构成（$X^* = -X, \operatorname{tr}(X)=0$）。

---

## 55B.3 伴随表示与结构常数

!!! definition "定义 55B.2 (伴随表示 $\operatorname{ad}$)"
    对于 $X \in \mathfrak{g}$，算子 $\operatorname{ad}_X: \mathfrak{g} \to \mathfrak{g}$ 定义为 $\operatorname{ad}_X(Y) = [X, Y]$。
    这是李代数在其自身上的线性表示，其矩阵形式的迹与 **Killing 型**（判定半单性）密切相关。

!!! technique "结构常数"
    选定基 $\{E_i\}$，则 $[E_i, E_j] = \sum C_{ij}^k E_k$。常数 $C_{ij}^k$ 包含了李代数的所有局时代数信息。

---

## 练习题

1. **[基础] 计算 $\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$ 与 $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的李括号。**

   ??? success "参考答案"
       $AB - BA = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} - \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ -2 & 0 \end{pmatrix}$。

2. **[Jacobi] 验证 $[x, [y, y]] = 0$ 是否符合 Jacobi 恒等式。**

   ??? success "参考答案"
       由于 $[y, y]=0$，左边为 0。Jacobi 恒等式在其中一个项为 0 时显然成立。

3. **[迹性质] 证明：若 $A, B \in \mathfrak{sl}(n)$，则 $[A, B] \in \mathfrak{sl}(n)$。**

   ??? success "参考答案"
       $\operatorname{tr}(AB - BA) = \operatorname{tr}(AB) - \operatorname{tr}(BA) = 0$。满足迹为 0 的条件。

4. **[反对称] 证明 $\mathfrak{so}(3)$ 的维数是 3。**

   ??? success "参考答案"
       $3 \times 3$ 反对称阵的形式为 $\begin{pmatrix} 0 & a & b \\ -a & 0 & c \\ -b & -c & 0 \end{pmatrix}$，由 3 个独立参数 $a, b, c$ 确定。

5. **[量子力学] Pauli 矩阵的对易关系 $[i\sigma_x, i\pi_y]$ 对应哪个李代数？**

   ??? success "参考答案"
       对应 $\mathfrak{su}(2)$（或同构的 $\mathfrak{so}(3)$）。它们代表了角动量的无穷小生成元。

6. **[指数映射] 描述 $\exp: \mathfrak{so}(n) \to SO(n)$ 的几何意义。**

   ??? success "参考答案"
       它将一个旋转速度向量（反对称阵）映射为一个有限旋转矩阵。

7. **[中心] 什么是李代数的中心？**

   ??? success "参考答案"
       与所有元素括号运算均为 0 的元素集合 $Z(\mathfrak{g}) = \{x \in \mathfrak{g} : [x, y]=0, \forall y\}$。

8. **[伴随] 证明 $\operatorname{ad}_{[X, Y]} = [\operatorname{ad}_X, \operatorname{ad}_Y]$。**

   ??? success "参考答案"
       这等价于 Jacobi 恒等式：$[ [X, Y], Z ] = [X, [Y, Z]] - [Y, [X, Z]]$。

9. **[结构常数] 若结构常数全为 0，该李代数是什么类型？**

   ??? success "参考答案"
       **阿贝尔李代数**（Abelian）。此时李括号恒为 0，对应于交换群。

10. **[控制理论] 在机器人路径规划中，为什么需要计算 Lie 括号？**

   ??? success "参考答案"
        两个向量场（控制方向）的 Lie 括号代表了通过切换控制方向可以到达的新维度（非完整约束产生的侧向移动），是判定系统是否完全可控的关键。

## 本章小结

李代数是线性代数在对称性研究中的最高形态：

1.  **无穷小的解析力**：它通过切空间技术，将复杂的非线性群流形降维为线性的向量空间，使得对称性的分类变得可计算。
2.  **运算的统一性**：李括号抽象了矩阵交换子的本质，确立了物理算子（如动量、自旋）之间相互干涉的代数规则。
3.  **结构的映射**：通过指数映射与伴随表示，李代数建立了一套“以动代静”的分析模式，证明了局部的无穷小变化足以决定全局的变换规律。
