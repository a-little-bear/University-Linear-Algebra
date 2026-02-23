# 第 34 章 Schur 补

<div class="context-flow" markdown>

**前置**：矩阵基础 (Ch02) · 正定矩阵 (Ch16) · 行列式 (Ch03)

**本章脉络**：分块矩阵的消元动机 $\to$ Schur 补的定义 $\to$ 行列式分解公式 $\to$ 逆矩阵的分块公式（Banachiewicz 公式） $\to$ 正定性判定准则 $\to$ 惯性加性公式（Haynsworth） $\to$ 矩阵方程中的应用 $\to$ 统计学意义（偏相关与条件方差）

**延伸**：Schur 补是分块矩阵计算的“手术刀”；它不仅是大型线性方程组分治求解的核心，也是理解现代概率论中高斯过程与核方法 (Ch29) 的数学钥匙

</div>

在处理大规模系统时，我们通常将矩阵划分为子块。**Schur 补**（Schur Complement）正是通过局部消元得到的关键中间结构。它不仅给出了分块矩阵行列式与逆的显式表达，还深刻揭示了矩阵各部分之间的相关性。

---

## 34.1 定义与行列式公式

!!! definition "定义 34.1 (Schur 补)"
    设分块矩阵 $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$。若 $A$ 可逆，则 $D$ 关于 $A$ 的 **Schur 补** 定义为：
    $$S = D - C A^{-1} B$$

!!! theorem "定理 34.1 (行列式分解公式)"
    分块矩阵 $M$ 的行列式可以分解为：
    $$\det(M) = \det(A) \det(D - C A^{-1} B)$$
    这说明整个系统的“体积”等于左上角子系统的体积乘以其 Schur 补的体积。

---

## 34.2 逆矩阵的分块公式

!!! technique "技术：Banachiewicz 公式"
    若 $A$ 和 Schur 补 $S$ 均可逆，则 $M$ 的逆矩阵为：
    $$M^{-1} = \begin{pmatrix} A^{-1} + A^{-1} B S^{-1} C A^{-1} & -A^{-1} B S^{-1} \\ -S^{-1} C A^{-1} & S^{-1} \end{pmatrix}$$
    这是求解结构化大系统的数值基石。

---

## 34.3 正定性与惯性

!!! theorem "定理 34.2 (正定性判据)"
    设 $M = \begin{pmatrix} A & B \\ B^T & C \end{pmatrix}$ 为对称阵。则：
    $$M \succ 0 \iff A \succ 0 \text{ 且 } C - B^T A^{-1} B \succ 0$$

!!! theorem "定理 34.3 (Haynsworth 惯性公式)"
    矩阵 $M$ 的正、负惯性指数等于 $A$ 的惯性指数与 Schur 补 $S$ 的惯性指数之和：
    $$\operatorname{In}(M) = \operatorname{In}(A) + \operatorname{In}(D - C A^{-1} B)$$

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

Schur 补是分块代数的核心语法：


****：它展示了如何通过局部的可逆信息，将高维系统的复杂性压缩到低维的余项中，是分布式算法的数学前提。

****：Schur 补公式确立了分块矩阵全局稳定性与其子系统局部稳定性及交互项之间的定量联系。

****：在概率论中，Schur 补揭示了高斯分布在条件化过程中的代数本质，是处理变量间“剩余相关性”的终极工具。
