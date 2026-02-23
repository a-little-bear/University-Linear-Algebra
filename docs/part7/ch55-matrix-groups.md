# 第 55 章 矩阵群与经典 Lie 群

<div class="context-flow" markdown>

**前置**：矩阵代数(Ch2) · 群与对称性 · 矩阵指数(Ch13) · 微分学

**本章脉络**：矩阵 Lie 群 ($GL, SL, O, SO, U, SU, Sp$) → 指数映射 $\exp: \mathfrak{g} \to G$ → Lie 代数作为切空间 → 对易子与 Lie 括号 → 伴随表示 → Baker-Campbell-Hausdorff 公式 → 群拓扑与代数结构的关系

**延伸**：矩阵群是物理学（标准模型）和几何学（黎曼流形）中描述对称性的基本语言

</div>

矩阵 Lie 群是同时具有光滑流形结构的矩阵群。每个此类群都有一个关联的 **Lie 代数**，它捕捉了群在单位元附近的无穷小结构。连接群与其代数的纽带是**矩阵指数**，它将代数中的“速度”映射为群中的“旋转”或“位移”。

---

## 55.1 核心群与代数

!!! definition "定义 55.1 (矩阵 Lie 群)"
    若子群 $G \subseteq GL(n, \mathbb{C})$ 是 $GL(n, \mathbb{C})$ 的闭子集，则称 $G$ 为矩阵 Lie 群。

!!! theorem "定理 55.1 (Lie 代数)"
    群 $G$ 的 Lie 代数 $\mathfrak{g}$ 由满足对于所有 $t \in \mathbb{R}$ 都有 $e^{tX} \in G$ 的所有矩阵 $X$ 组成。该代数配备了 **Lie 括号** $[X, Y] = XY - YX$。

---

## 练习题

1. **[特殊线性群] 确定群 $SL(n, \mathbb{C})$ 的 Lie 代数 $\mathfrak{sl}(n, \mathbb{C})$。**
   ??? success "参考答案"
       $A \in SL(n, \mathbb{C}) \implies \det A = 1$。利用恒等式 $\det(e^{tX}) = e^{t \operatorname{tr}(X)}$，条件 $\det(e^{tX}) = 1$ 对所有 $t$ 成立蕴含 $\operatorname{tr}(X) = 0$。因此 $\mathfrak{sl}(n, \mathbb{C})$ 由所有迹为零的矩阵组成。

2. **[正交群] 证明 Lie 代数 $\mathfrak{so}(n, \mathbb{R})$ 由反对称矩阵组成。**
   ??? success "参考答案"
       $R \in SO(n) \implies R^T R = I$。设 $R(t) = e^{tX}$。则 $(e^{tX})^T e^{tX} = e^{tX^T} e^{tX} = I$。在 $t=0$ 处求导得 $X^T + X = 0$，故 $X^T = -X$。

3. **[酉群] 刻画特殊酉群 $SU(n)$ 的 Lie 代数 $\mathfrak{su}(n)$。**
   ??? success "参考答案"
       它由迹为零的反 Hermitian 矩阵组成（即 $X^* = -X$ 且 $\operatorname{tr}(X) = 0$）。该代数的维数为 $n^2 - 1$。

4. **[BCH公式] 写出 $e^X e^Y = e^Z$ 的 Baker-Campbell-Hausdorff 公式的前几项。**
   ??? success "参考答案"
       $Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}([X, [X, Y]] + [Y, [Y, X]]) + \dots$。这表明群乘法在局部上由 Lie 括号唯一确定。

5. **[伴随表示] 定义伴随表示 $\operatorname{Ad}_g(X)$。**
   ??? success "参考答案"
       $\operatorname{Ad}_g(X) = g X g^{-1}$。该算子通过群作用将 Lie 代数的一个元素映射到同代数的另一个元素，且保持 Lie 括号不变。

6. **[对易子] 证明 Lie 括号满足 Jacobi 恒等式：$[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$。**
   ??? success "参考答案"
       利用 $[X, Y] = XY - YX$ 展开并验证 12 个项成对抵消。该恒等式是 Lie 代数的基本公理。

7. **[中心] Lie 代数 $\mathfrak{gl}(n, \mathbb{C})$ 的中心是什么？**
   ??? success "参考答案"
       中心由与代数中所有其他矩阵都对易的矩阵组成。对于 $\mathfrak{gl}(n)$，这恰好是所有标量矩阵 $cI$。

8. **[Spin群] 阐述 $SU(2)$ 与 $SO(3)$ 通过同态建立的联系。**
   ??? success "参考答案"
       存在一个从 $SU(2)$ 到 $SO(3)$ 的 2 对 1 满同态。这意味着 $SU(2)$ 是 $SO(3)$ 的单位覆盖空间，这一事实在量子力学描述电子自旋时至关重要。

9. **[指数映射] 指数映射对于 $GL(n, \mathbb{C})$ 总是满射吗？**
   ??? success "参考答案"
       是的，对于 $GL(n, \mathbb{C})$，每个可逆矩阵都存在矩阵对数。然而对于 $SL(2, \mathbb{R})$，它不是满射的（例如，具有负特征值且含有非平凡 Jordan 块的矩阵没有实对数）。

10. **[物理学] 为什么 Lie 群是粒子物理标准模型的核心？**
    ??? success "参考答案"
        基本相互作用的对称性由规范群 $U(1)$（电磁）、$SU(2)$（弱力）和 $SU(3)$（强力）描述。粒子对应于这些矩阵群的不可约表示。

## 本章小结

本章探讨了矩阵群中代数、拓扑与几何的综合：

1. **无穷小分析**：利用矩阵指数将群线性化为 Lie 代数。
2. **标准对称性**：通过迹与对称性约束刻画了正交、酉及特殊线性群的代数。
3. **代数结构**：定义了 Lie 括号作为捕捉非对易性的核心运算。
4. **全局-局部对偶**：通过 BCH 公式和 Lie 定理展示了局部代数结构如何决定全局群行为。
