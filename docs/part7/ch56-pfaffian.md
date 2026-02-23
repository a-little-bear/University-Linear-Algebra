# 第 56 章 Pfaffian 与反对称阵

<div class="context-flow" markdown>

**Prerequisites**: 行列式 (Ch03) · 矩阵群 (Ch55) · 辛矩阵 (Ch53) · 组合数学基础

**Chapter Outline**: 反对称矩阵的特殊性质 $\to$ Pfaffian 的定义 $\to$ 与行列式的关系 $\operatorname{Pf}(A)^2 = \det(A)$ $\to$ 代数性质（相似变换下的不变性） $\to$ 基于完美匹配的显式公式 $\to$ 柯西-比内公式的 Pfaffian 版本 $\to$ 应用：平面图完美匹配计数（FKT 算法）、统计物理中的 Ising 模型、Majorana 费米子

**Extension**: Pfaffian 是反对称矩阵的“平方根”；它在处理具有奇偶对称性的物理系统和组合计数问题中展现出比行列式更精细的描述能力

</div>

对于对称矩阵，我们有正定性和特征值理论。而对于**反对称矩阵**（Skew-symmetric Matrices），最深刻的标量函数不是行列式，而是 **Pfaffian**。由于反对称阵的行列式总是某个多项式的平方，Pfaffian 恰好提取了那个多项式。本章将揭示这个“行列式的代数平方根”如何连接了线性代数、图论与量子物理。

---

## 56.1 Pfaffian 的定义

!!! definition "定义 56.1 (Pfaffian)"
    设 $A$ 是一个 $2n \times 2n$ 的反对称矩阵。其 **Pfaffian** $\operatorname{Pf}(A)$ 定义为：
    $$\operatorname{Pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    对于奇数阶反对称阵，约定其 Pfaffian 为 0。

!!! example "例 56.1"
    对于 $2 \times 2$ 阵 $A = \begin{pmatrix} 0 & a \\ -a & 0 \end{pmatrix}$，$\operatorname{Pf}(A) = a$。
    注意 $\det(A) = a^2 = \operatorname{Pf}(A)^2$。

---

## 56.2 核心性质

!!! theorem "定理 56.1 (与行列式的关系)"
    对于任何 $2n \times 2n$ 反对称矩阵 $A$：
    $$\det(A) = [\operatorname{Pf}(A)]^2$$

!!! theorem "定理 56.2 (变量替换性质)"
    对于任意 $2n \times 2n$ 矩阵 $M$：
    $$\operatorname{Pf}(M A M^T) = \det(M) \operatorname{Pf}(A)$$
    这说明 Pfaffian 在旋转变换下保持不变（当 $\det(M)=1$ 时）。

---

## 56.3 组合意义：FKT 算法

!!! technique "完美匹配计数"
    设 $G$ 是一个平面图，$A$ 为其关联的 Pfaffian 取向矩阵。则 $G$ 的完美匹配个数恰好等于 $|\operatorname{Pf}(A)|$。这就是著名的 **FKT 算法**，它允许我们在多项式时间内计算平面图的完美匹配总数，而一般图的这一问题是 #P-完全的。

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

Pfaffian 是反对称代数中的精致结构：


****：它在代数上完美解释了反对称阵行列式必须是非负平方项的本质，确立了这种特殊对称性下的基本标量。

****：通过 FKT 算法，Pfaffian 将原本属于指数级难度的计数问题拉回到多项式复杂度的安全区，展现了代数对组合结构的巨大降维能力。

****：从经典热力学到前沿拓扑超导，Pfaffian 作为描述配对（Pairing）现象的天然语言，证明了矩阵论在揭示物质微观秩序方面的深刻性。
