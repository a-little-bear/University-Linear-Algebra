# 第 47B 章 Fréchet 导数与矩阵函数分析

<div class="context-flow" markdown>

**前置**：矩阵函数(Ch13) · 矩阵微积分(Ch47A) · 范数(Ch15) · 线性算子

**本章脉络**：Fréchet 导数定义 → 矩阵函数的条件数 → Fréchet 导数的 Kronecker 形式 → Daleckii-Krein 定理（除商差分） → 通过 SVD 估计导数 → 复步法 (Complex Step) → 在稳定性与灵敏度分析中的应用

**延伸**：Fréchet 导数量化了矩阵函数（如 $\exp(A)$ 或 $\log(A)$）对输入矩阵 $A$ 微扰的敏感程度

</div>

虽然矩阵微积分 (Ch47A) 通常处理矩阵的标量函数，但 **Fréchet 导数** 描述了矩阵函数 $f(A)$ 本身在输入矩阵 $A$ 受到微扰时如何变化。对于非线性函数 $f$，Fréchet 导数 $L_f(A, E)$ 是一个线性算子，它将微扰 $E$ 映射到 $f(A)$ 的一阶变化量。这对于计算矩阵函数的**条件数**以及分析矩阵算法的数值稳定性至关重要。

---

## 47B.1 Fréchet 导数

!!! definition "定义 47B.1 (Fréchet 导数)"
    函数 $f$ 在 $A$ 处的 Fréchet 导数是一个线性算子 $L_f(A, \cdot)$，满足：
    $$f(A + E) = f(A) + L_f(A, E) + o(\|E\|)$$
    值 $L_f(A, E)$ 是 $f$ 在 $A$ 处沿方向 $E$ 的方向导数。

!!! theorem "定理 47B.1 (Daleckii-Krein 定理)"
    若 $A = Q \Lambda Q^*$ 可对角化，则在特征向量基下，Fréchet 导数的条目为：
    $$(Q^* L_f(A, E) Q)_{ij} = f[\lambda_i, \lambda_j] (Q^* E Q)_{ij}$$
    其中 $f[\lambda_i, \lambda_j]$ 是**除商差分**：
    $$f[\lambda_i, \lambda_j] = \begin{cases} \frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j} & \text{若 } \lambda_i \neq \lambda_j \\ f'(\lambda_i) & \text{若 } \lambda_i = \lambda_j \end{cases}$$

---

## 练习题

1. **[基础] 计算 $f(A) = A^2$ 的 Fréchet 导数。**
   ??? success "参考答案"
       展开 $(A+E)^2 = A^2 + AE + EA + E^2$。关于 $E$ 的线性部分为 $L_f(A, E) = AE + EA$。该算子通过从两侧乘以 $A$ 来作用于 $E$。

2. **[逆矩阵] 计算 $f(A) = A^{-1}$ 的 Fréchet 导数。**
   ??? success "参考答案"
       利用 $d(A^{-1}) = -A^{-1}(dA)A^{-1}$，Fréchet 导数即为算子 $L_f(A, E) = -A^{-1} E A^{-1}$。它描述了当矩阵微扰时逆矩阵的移动规律。

3. **[条件数] 定义矩阵函数的相对条件数 $\kappa_f(A)$。**
   ??? success "参考答案"
       $\kappa_f(A) = \frac{\|L_f(A)\| \|A\|}{\|f(A)\|}$，其中 $\|L_f(A)\|$ 是诱导算子范数。它衡量了输出函数值对输入矩阵相对误差的敏感度。

4. **[对易性] 在什么条件下 $L_f(A, E) = f'(A) E$？**
   ??? success "参考答案"
       该等式成立当且仅当 $A$ 与 $E$ 对易（$AE = EA$）。非对易性是矩阵导数表现为算子而非简单乘数的主要原因。

5. **[Daleckii-Krein] 利用 Daleckii-Krein 定理求 $f(A) = e^A$ 在 $A = \operatorname{diag}(\lambda_1, \lambda_2)$ 处的 $L_f(A, E)$。**
   ??? success "参考答案"
       $L_f(A, E)$ 的条目为：当 $i \neq j$ 时，$L_{ij} = \frac{e^{\lambda_i} - e^{\lambda_j}}{\lambda_i - \lambda_j} E_{ij}$；当 $i=j$ 时，$L_{ii} = e^{\lambda_i} E_{ii}$。该公式将指数函数的敏感度与特征值的分散程度联系起来。

6. **[Kronecker形式] 将算子 $L_f(A, E) = AE + EA$ 写为其 $n^2 \times n^2$ 的 Kronecker 矩阵形式。**
   ??? success "参考答案"
       其矩阵表示为 $K_f(A) = I \otimes A + A^T \otimes I$。该矩阵作用于 $\operatorname{vec}(E)$ 即产生 $\operatorname{vec}(L_f(A, E))$。

7. **[指数恒等式] 矩阵指数函数在原点 $A=0$ 处的 Fréchet 导数是什么？**
   ??? success "参考答案"
       $L_{\exp}(0, E) = E$。在一阶近似下，指数映射在零点附近表现为单位算子。

8. **[复步法] 描述利用复步法估计 $L_f(A, E)$ 的优势。**
   ??? success "参考答案"
       $L_f(A, E) \approx \operatorname{Im}(f(A + i h E) / h)$。由于不涉及两个接近数值的相减，它免疫了困扰标准有限差分法的数值抵消误差。

9. **[复合函数] 写出矩阵 Fréchet 导数的链式法则。**
   ??? success "参考答案"
       对于 $h(A) = g(f(A))$，导数是线性算子的复合：$L_h(A, E) = L_g(f(A), L_f(A, E))$。

10. **[对称性] 证明若 $A$ 和 $E$ 均是 Hermitian 的，则 $L_f(A, E)$ 也是 Hermitian 的（假设 $f$ 在 $\mathbb{R}$ 上取实值）。**
    ??? success "参考答案"
        由 Daleckii-Krein 公式，$(Q^* L Q)_{ij} = f[\lambda_i, \lambda_j] (Q^* E Q)_{ij}$。由于 $f[\lambda_i, \lambda_j] = f[\lambda_j, \lambda_i]$ 且 $Q^*EQ$ 是 Hermitian 的，其乘积仍为 Hermitian。故 $L$ 是 Hermitian 的。

## 本章小结

本章探讨了矩阵函数的灵敏度与微扰理论：

1. **算子微积分**：定义了 Fréchet 导数作为捕捉矩阵函数对微扰的一阶响应的线性算子。
2. **谱灵敏度**：利用 Daleckii-Krein 定理建立了导数与特征值除商差分之间的联系。
3. **稳定性量化**：确立了条件数作为评估矩阵函数计算数值可靠性的决定性指标。
4. **计算工具**：概述了 Kronecker 表示与复步法，用于实际估计矩阵导数。
