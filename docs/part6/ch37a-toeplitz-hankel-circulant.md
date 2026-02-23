# 第 37A 章 Toeplitz·Hankel·循环矩阵

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 傅里叶变换(Ch30) · 多项式(Ch0)

**本章脉络**：Toeplitz 矩阵（常数对角线） → Hankel 矩阵（常数反对角线） → 循环矩阵 → 循环矩阵的特征值与 DFT → 离散卷积 → 位移结构 → Levinson-Durbin 递归 → Szegő 极限定理

**延伸**：这些矩阵是数字信号处理（滤波）、时间序列分析和求解离散微分方程的计算骨架

</div>

具有结构对称性的矩阵，如 **Toeplitz**、**Hankel** 和**循环矩阵**，产生于组件间相互作用仅取决于其相对位移的系统中。Toeplitz 矩阵具有常数对角线，而循环矩阵是一种特殊情况，其中每一行都是前一行的循环移位。这些结构与**离散傅里叶变换 (DFT)** 和多项式乘法密切相关，允许使用超快算法（$O(n \log n)$）代替标准的 $O(n^2)$ 或 $O(n^3)$。

---

## 37A.1 结构定义与谱性质

!!! definition "定义 37A.1 (Toeplitz 与循环矩阵)"
    若 $T_{i,j} = t_{i-j}$，则矩阵 $T$ 是 **Toeplitz 矩阵**。
    若 $C_{i,j} = c_{(j-i) \pmod n}$，则矩阵 $C$ 是**循环矩阵**。

!!! theorem "定理 37A.1 (循环矩阵的特征值)"
    每个循环矩阵 $C$ 都可以被 DFT 矩阵 $F$ 酉对角化。其特征值是 $C$ 第一行的离散傅里叶变换：
    $$C = F^* \Lambda F, \quad \text{其中 } \Lambda = \operatorname{diag}(F \mathbf{c})$$

---

## 练习题

1. **[基础] 写出由向量 $(1, 2, 3)$ 生成的 $3 \times 3$ 循环矩阵。**
   ??? success "参考答案"
       $C = \begin{pmatrix} 1 & 2 & 3 \\ 3 & 1 & 2 \\ 2 & 3 & 1 \end{pmatrix}$。

2. **[Toeplitz与Hankel] 证明 $H$ 是 Hankel 矩阵当且仅当 $HJ$ 是 Toeplitz 矩阵，其中 $J$ 是交换矩阵（反对角单位阵）。**
   ??? success "参考答案"
       交换矩阵 $J$ 会反转列的顺序。由于 Hankel 矩阵具有常数反对角线，反转列序会将这些常数项对齐到标准对角线上，从而形成 Toeplitz 矩阵。

3. **[DFT映射] 计算 $3 \times 3$ 循环移位矩阵 $S = \begin{pmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{pmatrix}$ 的特征值。**
   ??? success "参考答案"
       这是一个第一行为 $(0, 1, 0)$ 的循环矩阵。其特征值为单位根：$1, \omega, \omega^2$，其中 $\omega = e^{2\pi i/3}$。

4. **[卷积] 解释为什么循环矩阵与向量的乘法等价于离散循环卷积。**
   ??? success "参考答案"
       乘积项 $(Cx)_i = \sum c_{(j-i) \pmod n} x_j$ 恰好符合序列 $c$ 与 $x$ 的循环卷积公式。利用 FFT，该运算可以在 $O(n \log n)$ 时间内完成。

5. **[求逆] 对 $n \times n$ Toeplitz 矩阵求逆的复杂度是多少？**
   ??? success "参考答案"
       普通求逆为 $O(n^3)$，但利用 Levinson-Durbin 或 Bareiss 算法可以在 $O(n^2)$ 内完成。利用位移结构（Displacement structure）甚至可以达到 $O(n \log^2 n)$。

6. **[对易性] 证明所有 $n \times n$ 循环矩阵彼此对易。**
   ??? success "参考答案"
       由于所有循环矩阵都由同一个酉矩阵 $F$ 对角化，故 $C_1 C_2 = (F^* \Lambda_1 F)(F^* \Lambda_2 F) = F^* (\Lambda_1 \Lambda_2) F = F^* (\Lambda_2 \Lambda_1) F = C_2 C_1$。

7. **[平稳过程] 在统计学中，Toeplitz 矩阵与平稳随机过程有何联系？**
   ??? success "参考答案"
       宽平稳（WSS）离散时间过程的协方差矩阵是 Toeplitz 矩阵，因为协方差 $E[X_t X_{t+k}]$ 仅取决于时间间隔 $k$，而与绝对时间 $t$ 无关。

8. **[谱密度] 简述 Szegő 极限定理。**
   ??? success "参考答案"
       当 $n \to \infty$ 时，由函数 $f(\theta)$ 生成的 Toeplitz 矩阵的特征值分布趋于函数 $f$ 的值的分布。这建立了离散谱与连续功率谱密度之间的桥梁。

9. **[Vandermonde] 循环矩阵与 Vandermonde 矩阵有何关系？**
   ??? success "参考答案"
       DFT 矩阵 $F$ 实际上是一种特殊的 Vandermonde 矩阵，其节点取在 $n$ 次单位根上。

10. **[乘积] 两个 Toeplitz 矩阵的乘积是否仍为 Toeplitz 矩阵？**
    ??? success "参考答案"
        不一定。与循环矩阵集不同，Toeplitz 矩阵集对乘法并不封闭。但两个循环矩阵的乘积仍是循环的。

## 本章小结

本章探讨了由平移不变性定义的矩阵结构：

1. **位移对称性**：定义了 Toeplitz、Hankel 和循环矩阵作为不同平移不变程度的线性算子表示。
2. **傅里叶对偶**：确立了循环矩阵作为 DFT 代数对偶的地位，实现了卷积运算的对角化加速。
3. **递归求逆**：概述了 Levinson-Durbin 等算法，展示了如何利用结构冗余实现平方阶复杂度的计算。
4. **渐近分析**：通过 Szegő 理论将大型 Toeplitz 矩阵的谱属性与连续频域函数联系起来。
