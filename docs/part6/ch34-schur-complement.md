# 第 34 章 Schur 补矩阵

<div class="context-flow" markdown>

**前置**：矩阵求逆(Ch2) · 行列式(Ch3) · 分块矩阵 · 正定性(Ch16)

**本章脉络**：Schur 补定义 → 分块 LU 分解 → 分块矩阵行列式公式 → 矩阵求逆引理 (Woodbury) → 通过 Schur 补判定正定性 → 惯性定理与 Haynsworth 定律 → 在高斯过程中的应用

**延伸**：Schur 补是部分高斯消元和高维协方差矩阵降维的代数引擎

</div>

**Schur 补**自然地产生于对分块矩阵进行高斯消元的过程中。给定一个分块矩阵，其中一个块的 Schur 补捕捉了在考虑了交叉相关性后，另一个块的“剩余”信息。它是计算大型系统逆矩阵和行列式的基本工具，并为分块矩阵的正定性提供了决定性的判据。

---

## 34.1 定义与分块分解

!!! definition "定义 34.1 (Schur 补)"
    设 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ 为分块矩阵。若 $A$ 可逆，则 $A$ 在 $M$ 中的 Schur 补定义为：
    $$M/A = D - C A^{-1} B$$

!!! theorem "定理 34.1 (行列式公式)"
    分块矩阵 $M$ 的行列式等于子块 $A$ 的行列式与其 Schur 补行列式的乘积：
    $$\det M = \det A \cdot \det(D - C A^{-1} B)$$

---

## 练习题

1. **[基础] 计算 $M = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ 中 $A = (2)$ 的 Schur 补。**
   ??? success "参考答案"
       $M/A = 2 - (1)(2)^{-1}(1) = 2 - 0.5 = 1.5$。注意 $\det M = 2(1.5) = 3$。

2. **[分块求逆] 写出利用 Schur 补表示的 $M^{-1}$ 公式。**
   ??? success "参考答案"
       $M^{-1} = \begin{pmatrix} A^{-1} + A^{-1}B(M/A)^{-1}CA^{-1} & -A^{-1}B(M/A)^{-1} \\ -(M/A)^{-1}CA^{-1} & (M/A)^{-1} \end{pmatrix}$。这允许通过求逆较小的子块来完成大型矩阵的求逆。

3. **[正定性] 证明 $M = \begin{pmatrix} A & B \\ B^* & D \end{pmatrix} \succ 0$ 当且仅当 $A \succ 0$ 且 $M/A \succ 0$。**
   ??? success "参考答案"
       利用合同变换 $M = \begin{pmatrix} I & 0 \\ B^* A^{-1} & I \end{pmatrix} \begin{pmatrix} A & 0 \\ 0 & M/A \end{pmatrix} \begin{pmatrix} I & A^{-1}B \\ 0 & I \end{pmatrix}$。根据 Sylvester 惯性定律，$M$ 与分块对角阵具有相同数量的正特征值。故 $M \succ 0 \iff A \succ 0$ 且 $M/A \succ 0$。

4. **[Woodbury] 简述 Woodbury 矩阵恒等式（矩阵求逆引理）。**
   ??? success "参考答案"
       $(A + UCV)^{-1} = A^{-1} - A^{-1}U(C^{-1} + VA^{-1}U)^{-1}VA^{-1}$。它将受扰动矩阵的逆与其扰动项的 Schur 补联系起来。

5. **[变分性质] 将 Schur 补表达为一个最小化问题的结果。**
   ??? success "参考答案"
       对于 $M \succ 0$，Schur 补 $M/A$ 出现在二次型最小化中：$\min_x \begin{pmatrix} x \\ y \end{pmatrix}^T \begin{pmatrix} A & B \\ B^T & D \end{pmatrix} \begin{pmatrix} x \\ y \end{pmatrix} = y^T (M/A) y$。最优解 $x = -A^{-1}By$ 消除了交叉项。

6. **[惯性定律] $M$ 的惯性与 $A$ 及其 Schur 补的惯性有何关系？**
   ??? success "参考答案"
       $\operatorname{In}(M) = \operatorname{In}(A) + \operatorname{In}(M/A)$。这一关于正、负、零特征值数量的加法法则对稳定性分析至关重要。

7. **[正交投影] 从投影角度解释 Schur 补。**
   ??? success "参考答案"
       Schur 补 $D - C A^{-1} B$ 代表了 $D$ 中与 $B$ 的列空间正交的部分（在加权内积意义下）。

8. **[概率论] 在高斯分布 $\begin{pmatrix} X_1 \\ X_2 \end{pmatrix} \sim N(\mu, \Sigma)$ 中，给定 $X_1$ 后 $X_2$ 的条件协方差是什么？**
   ??? success "参考答案"
       条件协方差恰好是 Schur 补 $\Sigma_{22} - \Sigma_{21} \Sigma_{11}^{-1} \Sigma_{12}$。它捕捉了观测到 $X_1$ 后 $X_2$ 剩余的不确定性。

9. **[单调性] 证明若 $M \succeq 0$，则 $M/A \preceq D$。**
   ??? success "参考答案"
       $M/A = D - C A^{-1} B$。由于 $A \succeq 0 \implies A^{-1} \succeq 0$，项 $C A^{-1} B$（对称情形下为 $B^* A^{-1} B$）是半正定的。从 $D$ 中减去一个半正定阵，结果必然在 Lowner 序下小于等于 $D$。

10. **[秩] $M$ 的秩与 $A$ 及其 Schur 补的秩有何关系？**
    ??? success "参考答案"
        $\operatorname{rank}(M) = \operatorname{rank}(A) + \operatorname{rank}(M/A)$。当 $B$ 的列空间包含在 $A$ 的列空间中且 $C$ 的行空间包含在 $A$ 的行空间中时，该等式成立。

## 本章小结

本章探讨了分块矩阵通过 Schur 补进行的降维与化简：

1. **分解核心**：将 Schur 补定位为分块 LU 分解和高斯消元的核心。
2. **谱判据**：确立了分区矩阵正定性与惯性的决定性条件。
3. **求逆引理**：形式化了 Woodbury 恒等式，作为低秩扰动下更新逆矩阵的强力工具。
4. **统计深度**：展示了其在高斯条件化和概率方差缩减中的基础作用。
