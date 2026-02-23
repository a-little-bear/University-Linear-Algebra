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

1. **[基础] 计算 $\begin{pmatrix} 0 & 3 \\ -3 & 0 \end{pmatrix}$ 的 Pfaffian。**
   ??? success "参考答案"
       $\operatorname{Pf} = 3$。

2. **[行列式] 若 $\det(A) = 16$，且 $A$ 是反对称阵，$\operatorname{Pf}(A)$ 可能的值是多少？**
   ??? success "参考答案"
       $\pm 4$。符号取决于矩阵元素的具体排列。

3. **[性质] 证明：若 $A$ 是 $3 \times 3$ 反对称阵，则 $\operatorname{Pf}(A) = 0$。**
   ??? success "参考答案"
       奇数阶反对称阵的行列式必为 0，故其 Pfaffian（作为平方根）也必为 0。

4. **[计算] 求 $J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} \oplus \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 的 Pfaffian。**
   ??? success "参考答案"
       对于分块对角反对称阵，$\operatorname{Pf}(A \oplus B) = \operatorname{Pf}(A)\operatorname{Pf}(B)$。故 $\operatorname{Pf}(J) = 1 \cdot 1 = 1$。

5. **[Cayley] 证明：对于 $4 \times 4$ 反对称阵，其 Pfaffian 有 3 项。**
   ??? success "参考答案"
       $\operatorname{Pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$。这是通过枚举匹配项得到的。

6. **[物理] 在超导理论中，Pfaffian 态对应什么？**
   ??? success "参考答案"
       对应于具有非阿贝尔统计特性的 Majorana 束缚态，是拓扑量子计算的候选者。

7. **[特征值] 反对称矩阵的特征值有什么特点？**
   ??? success "参考答案"
       特征值成共轭纯虚数对出现：$\pm i\lambda_1, \pm i\lambda_2, \ldots$。Pfaffian 等于这些虚部之积：$\prod \lambda_j$。

8. **[伴随] 证明 $\operatorname{Pf}(k A) = k^n \operatorname{Pf}(A)$。**
   ??? success "参考答案"
       由行列式性质 $\det(kA) = k^{2n} \det(A)$，取平方根即得。

9. **[组合] 为什么一般图的匹配计数比平面图难？**
   ??? success "参考答案"
       因为一般图不存在统一的 Pfaffian 取向（即无法给边赋正负号使得所有回路贡献一致），这反映了平面性在代数结构上的独特性。

10. **[应用] 简述 Pfaffian 在统计力学 Ising 模型中的作用。**
    ??? success "参考答案"
        通过将配分函数表示为格点关联矩阵的 Pfaffian，可以将统计求和转化为矩阵行列式的计算，从而求得精确解。

## 本章小结

Pfaffian 是反对称代数中的精致结构：

1.  **开方之美**：它在代数上完美解释了反对称阵行列式必须是非负平方项的本质，确立了这种特殊对称性下的基本标量。
2.  **组合的捷径**：通过 FKT 算法，Pfaffian 将原本属于指数级难度的计数问题拉回到多项式复杂度的安全区，展现了代数对组合结构的巨大降维能力。
3.  **物理的桥梁**：从经典热力学到前沿拓扑超导，Pfaffian 作为描述配对（Pairing）现象的天然语言，证明了矩阵论在揭示物质微观秩序方面的深刻性。
