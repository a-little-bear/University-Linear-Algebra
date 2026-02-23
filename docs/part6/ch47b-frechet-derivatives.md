# 第 47B 章 Fréchet 导数与高阶理论

<div class="context-flow" markdown>

**前置**：矩阵微积分基础 (Ch47A) · 矩阵分析 (Ch14) · 泛函分析初步

**本章脉络**：从有限维微积分到算子微积分 $\to$ Gateaux 导数（方向导数）定义 $\to$ Fréchet 导数（全微分）的定义与线性算子表示 $\to$ 矩阵函数的 Fréchet 导数 $L_f(A, E)$ $\to$ 核心恒等式：Kronecker 积形式的导数向量化 $\to$ 链式法则的高阶形式 $\to$ 逆矩阵与矩阵指数的高阶导数 $\to$ 应用：算法的条件数分析、非线性矩阵方程的牛顿法、矩阵流形上的优化算法

**延伸**：Fréchet 导数是研究矩阵函数“灵敏度”的终极语言；它不仅能告诉我们函数值如何随输入变化，还通过其线性算子的性质揭示了算子空间内部的微观曲率，是高级数值稳定性理论的基石

</div>

在初等矩阵微积分中，我们习惯于处理标量函数的梯度。然而，当映射本身是矩阵到矩阵（如 $f(A) = e^A$ 或 $f(A) = A^{-1}$）时，导数不再是一个简单的矩阵，而是一个**线性算子**。**Fréchet 导数** 正是这种无限小线性逼近的严格表达。本章将介绍如何刻画这些作用于矩阵空间的算子，并利用 Kronecker 积技术实现它们的数值化计算。

---

## 47B.1 Fréchet 导数与方向导数

!!! definition "定义 47B.1 (Fréchet 导数)"
    映射 $f: M_n \to M_n$ 在 $A$ 处的 Fréchet 导数是一个线性算子 $L_f(A, \cdot)$，满足：
    $$f(A + E) = f(A) + L_f(A, E) + o(\|E\|)$$
    其中 $E$ 是微小扰动。

!!! definition "定义 47B.2 (Gateaux 导数)"
    若全微分存在，其在方向 $E$ 上的导数可以通过极限定义：
    $$L_f(A, E) = \lim_{h \to 0} \frac{f(A + hE) - f(A)}{h} = \left. \frac{d}{dt} f(A + tE) \right|_{t=0}$$

---

## 47B.2 Kronecker 积表示

!!! theorem "定理 47B.1 (导数的向量化形式)"
    Fréchet 导数作为线性算子，其作用可以利用 Kronecker 积矩阵 $K_f(A)$ 表示：
    $$\operatorname{vec}(L_f(A, E)) = K_f(A) \operatorname{vec}(E)$$
    **意义**：这一公式将抽象的算子作用转化为标准的矩阵-向量乘法，是数值软件计算矩阵导数的基础。

---

## 47B.3 典型函数的 Fréchet 导数

!!! example "例 47B.1 (逆矩阵与指数)"
    1.  **逆矩阵** $f(A) = A^{-1}$：$L_f(A, E) = -A^{-1} E A^{-1}$。
    2.  **矩阵幂** $f(A) = A^2$：$L_f(A, E) = AE + EA$。
    3.  **矩阵指数** $f(A) = e^A$：其导数由积分公式 $\int_0^1 e^{A(1-s)} E e^{As} ds$ 给出。

---

## 练习题

**1. [基础] 利用极限定义求 $f(A) = A^2$ 的 Fréchet 导数。**

??? success "参考答案"
    **计算步骤：**
    1. 计算 $f(A+hE) = (A+hE)(A+hE) = A^2 + h(AE + EA) + h^2 E^2$。
    2. $f(A+hE) - f(A) = h(AE + EA) + h^2 E^2$。
    3. 除以 $h$ 并取极限 $h \to 0$：
    **结论**：$L_f(A, E) = AE + EA$。注意：由于矩阵不交换，结果不能写成 $2AE$。

**2. [向量化] 将上题中的 $L_f(A, E) = AE + EA$ 写成 Kronecker 积矩阵 $K_f(A)$ 的形式。**

??? success "参考答案"
    **推导：**
    1. $\operatorname{vec}(AE) = (I \otimes A) \operatorname{vec}(E)$。
    2. $\operatorname{vec}(EA) = (A^T \otimes I) \operatorname{vec}(E)$。
    3. $\operatorname{vec}(L) = (I \otimes A + A^T \otimes I) \operatorname{vec}(E)$。
    **结论**：$K_f(A) = I \otimes A + A^T \otimes I$。这正是 $A$ 与 $A^T$ 的 Kronecker 和。

**3. [逆矩阵] 证明：若 $f(A) = A^{-1}$，其 Fréchet 导数满足 $L_f(A, E) = -A^{-1} E A^{-1}$。**

??? success "参考答案"
    **证明：**
    1. 考虑恒等式 $(A+E)(A+E)^{-1} = I$。
    2. 展开左侧：$(A+E)(A^{-1} + L_f + \cdots) = I + E A^{-1} + A L_f + \cdots = I$。
    3. 忽略二阶项，令一阶项之和为 0：$E A^{-1} + A L_f = O$。
    4. 解出 $L_f$：$A L_f = -E A^{-1} \implies L_f = -A^{-1} E A^{-1}$。

**4. [二阶导数] 求 $f(A) = A^2$ 的二阶 Fréchet 导数 $D^2 f(A)(E, H)$。**

??? success "参考答案"
    **推导：**
    1. 一阶导数为 $L(A, E) = AE + EA$。
    2. 对 $A$ 沿方向 $H$ 再次求导：
    3. $\delta L = (A+H)E + E(A+H) - (AE + EA) = HE + EH$。
    **结论**：$D^2 f(A)(E, H) = HE + EH$。在此例中，二阶导数与 $A$ 无关（因为原函数是二次的）。

**5. [行列式] 计算 $\det(A)$ 的 Fréchet 导数（假设 $A$ 可逆）。**

??? success "参考答案"
    **利用 Jacobi 公式：**
    $\delta \det(A) = \det(A) \operatorname{tr}(A^{-1} E)$。
    由于结果是一个标量，这与 $47A$ 中的梯度公式 $\nabla \det A = \det(A) (A^{-1})^T$ 是等价的（内积意义下）。

**6. [链式法则] 若 $h(A) = f(g(A))$，写出其 Fréchet 导数的复合形式。**

??? success "参考答案"
    **公式：**
    $L_h(A, E) = L_f(g(A), L_g(A, E))$。
    **意义**：这说明导数算子的复合遵循线性算子的嵌套逻辑，而不是简单的矩阵乘法。

**7. [矩阵指数] 为什么 $e^A$ 的导数不是 $e^A E$？**

??? success "参考答案"
    **理由：**
    只有当扰动 $E$ 与 $A$ 交换时， $e^{A+E} = e^A e^E$ 这种指数律才成立。
    在一般情况下，$A$ 与 $E$ 不交换，导致泰勒级数展开中的各项位置无法任意移动。必须使用积分公式（如例 47B.1）来处理这种不交换性。

**8. [应用] 什么是矩阵函数的“条件数”？它如何通过 Fréchet 导数定义？**

??? success "参考答案"
    **定义：**
    $\kappa_f(A) = \frac{\|L_f(A, \cdot)\| \cdot \|A\|}{\|f(A)\|}$。
    其中 $\|L_f\|$ 是线性算子的算子范数。
    它衡量了矩阵函数计算对输入扰动的灵敏度。

**9. [计算] 求 $f(A) = A^3$ 的一阶导数。**

??? success "参考答案"
    **推导：**
    $\delta(AAA) = (\delta A)AA + A(\delta A)A + AA(\delta A)$。
    **结论**：$L_f(A, E) = EAA + AEA + AAE$。

**10. [数值算法] 简述如何利用复步扰动法（Complex Step）近似计算 Fréchet 导数。**

??? success "参考答案"
    **思想：**
    利用 $f(A + i h E) \approx f(A) + i h L_f(A, E)$。
    取虚部：$L_f(A, E) \approx \frac{\operatorname{Im}(f(A + i h E))}{h}$。
    这种方法避免了传统有限差分法的减法消去误差，具有极高的数值精度。

## 本章小结

Fréchet 导数是算子空间中的全微分工具：

1.  **线性逼近的本质**：它将复杂的非线性矩阵函数局部简化为作用于扰动矩阵的线性算子，确立了算子灵敏度分析的标准语言。
2.  **向量化的桥梁**：Kronecker 积矩阵 $K_f(A)$ 的引入，实现了从抽象算子理论到具体数值矩阵计算的跨越，是开发高性能数学库的基础。
3.  **不交换性的精算**：通过处理 $A^{-1}$ 和 $e^A$ 等核心函数的导数推导，本章揭示了矩阵微积分与标量微积分最深刻的区别——对顺序和交换性的严格尊重。
