# 第 56 章 Pfaffian 与反对称矩阵

<div class="context-flow" markdown>

**前置**：矩阵代数(Ch2) · 行列式(Ch3) · 外代数(Ch49) · 图论(Ch27)

**本章脉络**：Pfaffian 定义 → 与行列式的关系 $\det A = (\operatorname{pf} A)^2$ → 反对称矩阵的性质 → 图的完美匹配 → 平面图的 FKT 算法 → Pfaffian 定向 → 在统计力学中的应用（Ising 模型）

**延伸**：Pfaffian 是计算平面图完美匹配数量的核心工具，也描述了某些物质拓扑相（如 Majorana 费米子）的状态

</div>

对于反对称矩阵 $A$（满足 $A^T = -A$），其行列式总是其条目的某个多项式的平方。这个多项式被称为 **Pfaffian**，记作 $\operatorname{pf}(A)$。行列式涉及所有置换，而 Pfaffian 仅对集合划分为不相交对的情况求和，这在矩阵线性代数与组合匹配理论之间建立了一道深层的联系。

---

## 56.1 定义与基本关系

!!! definition "定义 56.1 (Pfaffian)"
    对于 $2n \times 2n$ 反对称矩阵 $A$，其 Pfaffian 定义为：
    $$\operatorname{pf}(A) = \frac{1}{2^n n!} \sum_{\sigma \in S_{2n}} \operatorname{sgn}(\sigma) \prod_{i=1}^n a_{\sigma(2i-1), \sigma(2i)}$$
    对于奇数维矩阵，Pfaffian 定义为零。

!!! theorem "定理 56.1 (Pfaffian-行列式恒等式)"
    对于任何反对称矩阵 $A$：
    $$\det(A) = [\operatorname{pf}(A)]^2$$

---

## 练习题

1. **[基础计算] 计算 $A = \begin{pmatrix} 0 & a \\ -a & 0 \end{pmatrix}$ 的 Pfaffian。**
   ??? success "参考答案"
       $\operatorname{pf}(A) = a$。注意 $\det(A) = a^2$，符合平方根恒等式。

2. **[展开式] 写出 $4 \times 4$ 反对称矩阵 Pfaffian 的显式展开项。**
   ??? success "参考答案"
       $\operatorname{pf}(A) = a_{12}a_{34} - a_{13}a_{24} + a_{14}a_{23}$。这对应于将 4 个顶点配对成完美匹配的三种可能方式。

3. **[图论] 解释 Pfaffian 如何与图的完美匹配相关联。**
   ??? success "参考答案"
       若 $A$ 是具有 Pfaffian 定向的图的邻接矩阵，则 $\operatorname{pf}(A)$ 恰好等于该图的完美匹配数量。这为组合计数问题提供了线性代数解法。

4. **[不变性] 置换矩阵 $P$ 对 $\operatorname{pf}(P A P^T)$ 有何影响？**
   ??? success "参考答案"
       $\operatorname{pf}(P A P^T) = \det(P) \operatorname{pf}(A)$。这意味着 Pfaffian 在顶点重编号下保持不变（除非置换是奇置换，此时变号）。

5. **[外代数] 利用双向量的楔积表达 Pfaffian。**
   ??? success "参考答案"
       设 $\omega = \sum_{i<j} a_{ij} e_i \wedge e_j$。则 $n$ 阶外幂满足 $\frac{1}{n!} \omega^n = \operatorname{pf}(A) e_1 \wedge \dots \wedge e_{2n}$。这为 Pfaffian 提供了几何基础。

6. **[奇数维] 证明奇数阶反对称矩阵的行列式必为零。**
   ??? success "参考答案"
       $\det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)$。若 $n$ 为奇数，则 $\det(A) = -\det(A)$，故 $\det(A) = 0$。

7. **[特征值] 描述实反对称矩阵的谱。**
   ??? success "参考答案"
       其特征值均为纯虚数，且成共轭对 $\pm i\theta_j$ 出现。由于 $\det A = \prod (i\theta_j)(-i\theta_j) = \prod \theta_j^2$，故其行列式总长非负。

8. **[FKT算法] FKT 算法的核心思想是什么？**
   ??? success "参考答案"
       Fisher-Kasteleyn-Temperley 算法通过构造一种特殊的定向（Pfaffian 定向），使得定向邻接矩阵的项在 Pfaffian 求和中不会因正负号相互抵消，从而在多项式时间内计算平面图的匹配数。

9. **[规范形] 叙述反对称矩阵在合同变换下的规范形。**
   ??? success "参考答案"
       对于每个反对称矩阵 $A$，都存在非奇异阵 $P$ 使得 $P A P^T$ 为分块对角阵，其对角块为 $\begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ 或零。

10. **[物理学] Pfaffian 在统计力学中出现在何处？**
    ??? success "参考答案"
        在 Ising 模型中，平面点阵上的配分函数可以表示为一个 Pfaffian。它还刻画了拓扑超导体（如 Majorana 费米子）的基态。

## 本章小结

本章探讨了反对称矩阵特有的多项式平方根结构：

1. **多项式平方根**：定义了 Pfaffian 作为反对称矩阵行列式的规范“平方根”。
2. **组合匹配**：确立了 Pfaffian 项与图中完美匹配之间的一一对应关系。
3. **拓扑算法**：详述了利用线性代数计算平面图匹配数的 FKT 算法。
4. **不变量理论**：利用外幂与合同变换形式化了 Pfaffian 的代数性质。
