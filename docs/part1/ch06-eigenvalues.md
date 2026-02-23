# 第 06 章 特征值与特征向量

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Determinants (Ch3) · Linear Transformations (Ch5)

**本章脉络**：特征值方程 $Av = \lambda v$ → 特征多项式 → 几何重数与代数重数 → 矩阵对角化条件 → 相似变换 $P^{-1}AP = D$ → 对称矩阵的谱定理 → Cayley-Hamilton 定理

**延伸**：特征值揭示了线性算子在某些特定方向上的纯粹缩放行为，是分析动态系统稳定性的金钥匙

</div>

特征值理论是线性代数中最深刻的部分之一。它将复杂的矩阵运算降维为标量的乘法。通过寻找矩阵的“固有频率”（特征值）和“固有振型”（特征向量），我们可以解构算子的核心结构。

---

## 06.1 核心定义

!!! definition "定义 06.1 (特征值与特征向量)"
    对于 $n \times n$ 矩阵 $A$，若存在非零向量 $v$ 和标量 $\lambda$ 使得：
    $$Av = \lambda v$$
    则称 $\lambda$ 为 $A$ 的**特征值**，$v$ 为对应于 $\lambda$ 的**特征向量**。

!!! theorem "定理 06.3 (特征值与迹、行列式的关系)"
    1. $\sum \lambda_i = \operatorname{tr}(A)$
    2. $\prod \lambda_i = \det(A)$

---

## 练习题

1. **[基础计算] 计算 $A = \begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix}$ 的特征值。**
   ??? success "参考答案"
       特征方程：$\det(\lambda I - A) = (\lambda-4)(\lambda-3) - 2 = \lambda^2 - 7\lambda + 10 = 0$。
       解得 $\lambda_1 = 5, \lambda_2 = 2$。

2. **[特征向量] 求 $A = \begin{pmatrix} 4 & 1 \\ 2 & 3 \end{pmatrix}$ 对应 $\lambda=5$ 的特征向量。**
   ??? success "参考答案"
       解 $(5I-A)v = 0 \implies \begin{pmatrix} 1 & -1 \\ -2 & 2 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = 0$。
       得 $x_1 = x_2$。特征向量为 $c \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ (取 $c \neq 0$)。

3. **[对角化判定] 判定 $A = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ 是否可以对角化。**
   ??? success "参考答案"
       特征值为 $\lambda = 1, 1$（代数重数为 2）。
       求解特征向量：$(I-A)v = 0 \implies \begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix} v = 0 \implies x_2 = 0$。
       特征空间维数（几何重数）为 1。由于几何重数 < 代数重数，矩阵不可对角化。

4. **[迹与行列式] 若 $3 \times 3$ 矩阵的特征值为 1, 2, 3，求 $\operatorname{tr}(A)$ 和 $\det(A)$。**
   ??? success "参考答案"
       $\operatorname{tr}(A) = 1+2+3 = 6$。
       $\det(A) = 1 \cdot 2 \cdot 3 = 6$。

5. **[谱定理初步] 证明：实对称矩阵的特征值均为实数。**
   ??? success "参考答案"
       设 $Av = \lambda v$。则 $\bar{v}^T Av = \lambda \bar{v}^T v$。
       取共轭转置并利用 $A^T=A$，得 $\bar{v}^T A v = \bar{\lambda} \bar{v}^T v$。
       故 $\lambda \|v\|^2 = \bar{\lambda} \|v\|^2$。因 $v \neq 0$，必有 $\lambda = \bar{\lambda}$，即 $\lambda$ 为实数。

6. **[幂运算] 若 $A$ 的特征值为 $\lambda$，证明 $A^k$ 的特征值为 $\lambda^k$。**
   ??? success "参考答案"
       $Av = \lambda v \implies A(Av) = A(\lambda v) = \lambda (Av) = \lambda^2 v$。
       通过数学归纳法可得 $A^k v = \lambda^k v$。

7. **[Cayley-Hamilton] 验证 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$ 满足其特征方程。**
   ??? success "参考答案"
       特征多项式 $p(\lambda) = (\lambda-1)(\lambda-2) = \lambda^2 - 3\lambda + 2$。
       $p(A) = A^2 - 3A + 2I = \begin{pmatrix} 1 & 0 \\ 0 & 4 \end{pmatrix} - \begin{pmatrix} 3 & 0 \\ 0 & 6 \end{pmatrix} + \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$。

8. **[逆矩阵特征值] 若 $A$ 可逆且特征值为 $\lambda$，证明 $A^{-1}$ 的特征值为 $1/\lambda$。**
   ??? success "参考答案"
       $Av = \lambda v \implies v = A^{-1}(\lambda v) = \lambda A^{-1} v$。
       由于 $A$ 可逆，$\lambda \neq 0$。两边除以 $\lambda$ 得 $A^{-1} v = (1/\lambda) v$。

9. **[相似性] 证明相似矩阵具有相同的特征值。**
   ??? success "参考答案"
       $\det(\lambda I - P^{-1}AP) = \det(P^{-1}(\lambda I - A)P) = \det(P^{-1})\det(\lambda I - A)\det(P) = \det(\lambda I - A)$。特征多项式相同，故特征值相同。

10. **[代数重数] 矩阵 $A = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 2 \end{pmatrix}$ 的特征值 2 的代数重数是多少？**
    ??? success "参考答案"
        特征多项式为 $(\lambda-2)^3$。因此特征值 2 的代数重数是 3。

## 本章小结

特征值理论是线性代数的灵魂：

1. **解耦作用**：对角化本质上是将复杂的耦合系统分解为相互独立的线性分量。
2. **重数博弈**：代数重数与几何重数的匹配是矩阵能否简化的关键。
3. **结构守恒**：迹和行列式作为特征值的函数，是不变量理论的基石。
