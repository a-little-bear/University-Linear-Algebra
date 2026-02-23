# 第 12 章 Jordan 标准形

<div class="context-flow" markdown>

**前置**：特征值与重数(Ch6) · 相似变换(Ch5) · 极小多项式(Ch13b)

**本章脉络**：对角化的局限性 → 广义特征向量 → Jordan 块 $J_k(\lambda)$ → Jordan 标准形定理 → 计算步骤（链式构造） → 唯一性（初等因子） → 应用（矩阵幂、微分方程解的稳定性）

**延伸**：Jordan 标准形是相似变换下最精简的结构，它揭示了矩阵在无法对角化时的“骨折”补偿机制（次对角线的 1）

</div>

当矩阵的几何重数小于代数重数时，对角化宣告失败。**Jordan 标准形**是这种失败情况下的最优代数补偿。它通过在对角线上方引入 1（广义特征向量的耦合），建立了一个通用的相似变换终点。

---

## 12.1 核心结构与定理

!!! definition "定义 12.1 (Jordan 块)"
    $k$ 阶 Jordan 块 $J_k(\lambda)$ 是一个对角线上均为 $\lambda$，次对角线上全为 1 的方阵：
    $$J_k(\lambda) = \begin{pmatrix} \lambda & 1 & & \\ & \lambda & \ddots & \\ & & \ddots & 1 \\ & & & \lambda \end{pmatrix}$$

!!! theorem "定理 12.1 (Jordan 标准形定理)"
    每一个复方阵 $A$ 都相似于一个 Jordan 标准形 $J$，且 $J$ 在 Jordan 块的排列次序外是唯一的。

---

## 练习题

1. **[基础] 写出 $J_2(3)$ 的具体形式。**
   ??? success "参考答案"
       $J_2(3) = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$。

2. **[特征值判定] 矩阵 $\begin{pmatrix} 2 & 1 \\ 0 & 2 \end{pmatrix}$ 的特征值是多少？它能对角化吗？**
   ??? success "参考答案"
       特征值为 2（重数 2）。它本身就是一个 Jordan 块。由于其特征向量只有 1 个（几何重数 1 < 代数重数 2），故不可对角化。

3. **[Jordan块数量] 矩阵 $A$ 的 Jordan 标准形中，对应特征值 $\lambda$ 的 Jordan 块总数由什么决定？**
   ??? success "参考答案"
       等于 $\lambda$ 的**几何重数**，即零空间 $\ker(\lambda I - A)$ 的维数。每个 Jordan 块对应一条独立的广义特征向量链。

4. **[计算] 求 $A = \begin{pmatrix} 5 & 4 \\ -1 & 1 \end{pmatrix}$ 的 Jordan 标准形。**
   ??? success "参考答案"
       特征多项式 $(\lambda-5)(\lambda-1)+4 = \lambda^2-6\lambda+9 = (\lambda-3)^2$。
       特征值 $\lambda=3$（重数 2）。
       计算几何重数：$3I-A = \begin{pmatrix} -2 & -4 \\ 1 & 2 \end{pmatrix}$，秩为 1，故零空间维数为 1。
       对应的 Jordan 标准形为 $J_2(3) = \begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}$。

5. **[广义特征向量] 定义 $\lambda$ 的 2 阶广义特征向量。**
   ??? success "参考答案"
       满足 $(A-\lambda I)^2 v = 0$ 但 $(A-\lambda I) v \neq 0$ 的向量 $v$。它是链 $[(A-\lambda I)v, v]$ 的末端。

6. **[幂运算] 计算 $J = \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$ 的 $n$ 次幂。**
   ??? success "参考答案"
       $J^n = \begin{pmatrix} \lambda^n & n\lambda^{n-1} \\ 0 & \lambda^n \end{pmatrix}$。注意右上角的项是利用二项式展开或导数性质得到的。

7. **[唯一性] 为什么 Jordan 标准形在相似变换下是唯一的？**
   ??? success "参考答案"
       因为 Jordan 块的大小和数量由各级核空间的维数序列 $\dim \ker(A-\lambda I)^k$ 唯一确定。这些维数是相似不变量。

8. **[迹与行列式] 若 $A$ 的 Jordan 标准形由 $J_2(1)$ 和 $J_1(2)$ 组成，求 $\det(A)$。**
   ??? success "参考答案"
       特征值为 $\{1, 1, 2\}$。$\det(A) = 1 \cdot 1 \cdot 2 = 2$。

9. **[最小多项式初步] 若 $A$ 的 Jordan 标准形中最大的一块 $\lambda$ 的大小是 $k$，则它的最小多项式中 $(\lambda-x)$ 的幂次是多少？**
   ??? success "参考答案"
       幂次恰好为 $k$。最小多项式捕捉了使每个 Jordan 块变为零矩阵所需的最小幂次。

10. **[应用] Jordan 标准形如何解释微分方程 $\dot{x} = Ax$ 中解的 $t e^{\lambda t}$ 项？**
    ??? success "参考答案"
        当 $A$ 不可对角化时，其指数矩阵 $e^{At} = P e^{Jt} P^{-1}$。Jordan 块 $J_k(\lambda)$ 的指数包含 $t^m e^{\lambda t}$ 项（$m < k$），这对应于共振或临界阻尼情况下的解增长。

## 本章小结

Jordan 标准形是相似类理论的终点：

1. **结构完备性**：解决了所有矩阵在相似变换下的分类问题。
2. **广义链逻辑**：通过广义特征向量链，将退化的线性算子理顺。
3. **计算基石**：它是研究矩阵函数（尤其是指数函数）和长期动态行为的最简模型。
