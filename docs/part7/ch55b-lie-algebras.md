# 第 55B 章 李代数与无限小生成元

<div class="context-flow" markdown>

**前置**：矩阵群(Ch55) · 矩阵指数(Ch13) · 向量空间(Ch4) · 微分几何

**本章脉络**：抽象李代数定义 → 作为括号的换位子 → 结构常数 → 李代数的表示 → 伴随表示 $\operatorname{ad}$ → 可解与幂零代数 → Engel 定理与 Lie 定理 → Killing 型 → Cartan 判据 → 半单李代数与根系

**延伸**：李代数是李群的“线性化”版本；它们允许利用纯线性代数工具（特征值、根、权）来研究连续对称性

</div>

虽然李群是流形，但**李代数**是向量空间。通过在单位元处线性化群作用，我们将非线性几何问题转化为线性代数问题。括号运算 $[X, Y]$ 替代了群乘法，捕捉了对称群的一阶非对角性。

---

## 55B.1 公理与基本结构

!!! definition "定义 55B.1 (李代数)"
    域 $K$ 上的李代数 $\mathfrak{g}$ 是一个配备了双线性映射 $[\cdot, \cdot]: \mathfrak{g} \times \mathfrak{g} \to \mathfrak{g}$（李括号）的向量空间，满足：
    1. **反对称性**：$[X, Y] = -[Y, X]$
    2. **Jacobi 恒等式**：$[X, [Y, Z]] + [Y, [Z, X]] + [Z, [X, Y]] = 0$

!!! theorem "定理 55B.3 (Engel 定理)"
    一个有限维李代数（在其伴随表示中）由幂零算子组成，当且仅当该代数是**幂零**的。

---

## 练习题

1. **[基础] 验证具有换位子括号 $[A, B] = AB - BA$ 的 $M_n(K)$ 满足 Jacobi 恒等式。**
   ??? success "参考答案"
       展开括号：$[A, [B, C]] = A(BC-CB) - (BC-CB)A = ABC - ACB - BCA + CBA$。对 $(A, B, C)$ 的所有循环置换求和，可以发现 12 个项成对抵消（例如 $ABC$ 与 $[B, [C, A]]$ 中的 $-ABC$ 抵消）。

2. **[结构常数] 若 $[e_i, e_j] = \sum c_{ij}^k e_k$，反对称性如何约束结构常数 $c_{ij}^k$？**
   ??? success "参考答案"
       $c_{ij}^k = -c_{ji}^k$。这减少了定义李代数结构所需的独立参数数量。

3. **[幂零性] 证明严格上三角矩阵代数是幂零的。**
   ??? success "参考答案"
       两个严格上三角矩阵的括号运算会增加零次对角线的数量。经过 $n$ 步嵌套括号后，结果矩阵必为零。这是幂零李代数的原型。

4. **[伴随表示] 定义伴随算子 $\operatorname{ad}_X: \mathfrak{g} \to \mathfrak{g}$。**
   ??? success "参考答案"
       $\operatorname{ad}_X(Y) = [X, Y]$。映射 $X \mapsto \operatorname{ad}_X$ 是李代数到其自身自同态空间 $\mathfrak{gl}(\mathfrak{g})$ 的一个表示，称为**伴随表示**。

5. **[Killing型] 什么是 Killing 型 $B(X, Y)$？**
   ??? success "参考答案"
       $B(X, Y) = \operatorname{tr}(\operatorname{ad}_X \circ \operatorname{ad}_Y)$。这种对称双线性型为李代数提供了几何度量。Cartan 证明了一个代数是半单的当且仅当其 Killing 型是非退化的。

6. **[可解性] 叙述关于 $\mathbb{C}$ 上可解李代数的 Lie 定理。**
   ??? success "参考答案"
       可解李代数的每个有限维表示都有一个公共特征向量。这意味着代数中的所有矩阵都可以在某个基下同时上三角化。

7. **[半单性] 定义半单李代数。**
   ??? success "参考答案"
       若一个代数没有非零的可解理想，则称其为半单的。这类代数（如 $\mathfrak{sl}_n$）是李理论的“原子”，可以通过根系和 Dynkin 图进行完全分类。

8. **[计算] 在 $\mathfrak{sl}_2(\mathbb{C})$ 中计算 $[H, E]$，其中 $H = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}, E = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。**
   ??? success "参考答案"
       $[H, E] = HE - EH = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} - \begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 2 \\ 0 & 0 \end{pmatrix} = 2E$。这表明 $E$ 是 $H$ 的对应于根（特征值）2 的根向量。

9. **[Cartan子代数] 什么是 Cartan 子代数 $\mathfrak{h} \subseteq \mathfrak{g}$？**
   ??? success "参考答案"
       它是一个由可对角化元素构成的极大交换子代数。它充当了代数的“对角”部分，允许将 $\mathfrak{g}$ 分解为各个根空间。

10. **[表示论] 将李代数表示与李群表示联系起来。**
    ??? success "参考答案"
        若 $\pi: G \to GL(V)$ 是一个群表示，则其微分 $d\pi: \mathfrak{g} \to \mathfrak{gl}(V)$ 是一个李代数表示。对于单连通群，这种对应是双射，允许我们利用线性代数来分类群对称性。

## 本章小结

本章探讨了无穷小对称性的线性框架：

1. **代数基础**：通过括号运算和 Jacobi 恒等式定义了李代数，将非交换概念线性化。
2. **结构层次**：根据导列和理想将代数分类为幂零、可解和半单类型。
3. **迹度量**：引入 Killing 型作为判定半单性和全局代数结构的决定性工具。
4. **谱分析**：将表示的研究与根系几何及公共特征向量（Lie 定理）联系起来。
