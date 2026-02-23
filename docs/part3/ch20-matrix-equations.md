# 第 20 章 矩阵方程

<div class="context-flow" markdown>

**前置**：Kronecker 积(Ch19) · 特征值(Ch6) · 矩阵分析(Ch14) · 矩阵稳定性(Ch36)

**本章脉络**：Sylvester 方程 $AX - XB = C$ → 解的存在性条件 → Lyapunov 方程 $AX + XA^* = Q$ → 稳定性判定与惯性定理 → 代数 Riccati 方程 (ARE) → 迭代解法 → 应用（现代控制理论中的能控性与能观性）

**延伸**：矩阵方程是将静态代数转化为动态特性的语言，Lyapunov 方程的解直接映射了动力系统的“能量”分布

</div>

矩阵方程研究的是以矩阵为未知量的代数方程。它们通常出现在控制系统设计、平衡态分析和数值分析中。最著名的 Sylvester 方程和 Lyapunov 方程为我们理解算子之间的相互作用提供了数学框架。

---

## 20.1 Sylvester 与 Lyapunov 方程

!!! definition "定义 20.1 (Sylvester 方程)"
    矩阵方程 $AX - XB = C$ 称为 Sylvester 方程。

!!! theorem "定理 20.3 (解的存在唯一性)"
    Sylvester 方程 $AX - XB = C$ 有唯一解，当且仅当 $A$ 与 $B$ 没有公共的特征值：$\sigma(A) \cap \sigma(B) = \emptyset$。

---

## 练习题

1. **[基础判定] 方程 $AX - XA = 0$ 总是存在非零解吗？**
   ??? success "参考答案"
       是的。任何与 $A$ 对易的矩阵都是该方程的解。特别地，$X=I$ 和 $X=A$ 总是解。由于 $\sigma(A) \cap \sigma(A) \neq \emptyset$（除非为空集），由定理 20.3 可知该方程没有唯一解。

2. **[Lyapunov方程] 设 $A$ 是 Hurwitz 稳定的。证明对于任何 $Q$，Lyapunov 方程 $AX + XA^* = Q$ 的解可以表示为积分形式。**
   ??? success "参考答案"
       解为 $X = -\int_0^\infty e^{At} Q e^{A^*t} dt$。由于 $A$ 是稳定的，指数项随时间衰减，保证了积分的收敛性。

3. **[计算] 求解 $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} X - X \begin{pmatrix} 0 & 0 \\ 0 & 3 \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。**
   ??? success "参考答案"
       设 $X = \begin{pmatrix} x_{11} & x_{12} \\ x_{21} & x_{22} \end{pmatrix}$。
       代入得：$\begin{pmatrix} 1x_{11} & 1x_{12} \\ 2x_{21} & 2x_{22} \end{pmatrix} - \begin{pmatrix} 0 & 3x_{12} \\ 0 & 3x_{22} \end{pmatrix} = \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$。
       解方程组：$x_{11}=1, -2x_{12}=1, 2x_{21}=1, -x_{22}=1$。
       故 $X = \begin{pmatrix} 1 & -0.5 \\ 0.5 & -1 \end{pmatrix}$。

4. **[稳定性判定] 证明：若存在 $P \succ 0$ 使得 $A^T P + PA \prec 0$，则 $A$ 是稳定的。**
   ??? success "参考答案"
       这是 Lyapunov 稳定性定理。考虑能量函数 $V(x) = x^T P x$。
       其导数 $\dot{V} = x^T (A^T P + PA) x$。由于 $A^T P + PA$ 是负定的，$\dot{V} < 0$，说明能量随时间减少，系统收敛。

5. **[迹的应用] 在 Lyapunov 方程 $AX + XA^T + BB^T = 0$ 中，$X$ 的迹代表什么？**
   ??? success "参考答案"
       若该方程描述控制系统的能控性 Gramian，$X$ 的迹（ eigenvalues 之和）代表了从各个方向输入控制能量的“平均增益”或总体能控性程度。

6. **[Riccati方程初步] 什么是代数 Riccati 方程 (ARE)？它与 Sylvester 方程的区别是什么？**
   ??? success "参考答案"
       ARE 形式如 $A^T P + PA - PBR^{-1}B^T P + Q = 0$。与线性 Sylvester 方程不同，ARE 是**二次**矩阵方程。它在最优控制（LQR）中起核心作用。

7. **[能控性] 若系统能控性矩阵为 $W_c$，它满足什么矩阵方程？**
   ??? success "参考答案"
       满足 Lyapunov 方程 $A W_c + W_c A^T + BB^T = 0$。解 $W_c$ 的正定性直接对应系统的能控性。

8. **[解的对称性] 若 $A$ 稳定且 $Q$ 是对称阵，Lyapunov 方程 $AX + XA^T = Q$ 的解 $X$ 一定是对称阵吗？**
   ??? success "参考答案"
       是的。对原方程取转置得 $XA^T + AX^T = Q^T = Q$。由于解是唯一的，且 $X$ 和 $X^T$ 满足同样的方程，故 $X = X^T$。

9. **[初等算子] Sylvester 方程可以看作什么线性算子的逆问题？**
   ??? success "参考答案"
       可以看作作用在矩阵空间上的初等算子 $\mathcal{L}(X) = AX - XB$ 的求逆问题。其特征值是 $\lambda_i(A) - \mu_j(B)$。

10. **[应用] 在图像复原中，矩阵方程如何模拟退化过程？**
    ??? success "参考答案"
        模糊过程常建模为 $Y = AXB + N$，其中 $A, B$ 分别代表水平和垂直方向的模糊算子。复原图像即求解该矩阵方程（通常结合正则化项）。

## 本章小结

矩阵方程是线性代数的动态逻辑：

1. **谱的相互作用**：解的存亡取决于两个算子谱集的隔离程度。
2. **能量映射**：Lyapunov 方程在代数与动力学稳定性之间架起了直接的桥梁。
3. **优化核心**：从线性方程到二次方程（Riccati），矩阵方程构成了现代控制理论的计算底座。
