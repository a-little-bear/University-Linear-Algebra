# 第 19 章 Kronecker 积与 Vec 算子

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch02) · 矩阵方程初步 (Ch20) · 多线性代数与张量 (Ch21)

**本章脉络**：从普通乘法到张量积 $\to$ Kronecker 积的定义与分块结构 $\to$ 核心代数性质（结合律、转置、逆） $\to$ 混合乘积性质 (Mixed-product property) $\to$ 谱性质：Kronecker 积的特征值与迹 $\to$ Vec 算子及其线性本质 $\to$ 核心恒等式：$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ $\to$ Kronecker 和 ($A \oplus B$) $\to$ 应用：线性矩阵方程（Sylvester, Lyapunov）的向量化求解、多变量系统的建模

**延伸**：Kronecker 积是张量积在矩阵范畴下的具象化；它通过“维度乘法”构建了一个描述多物理量相互耦合的代数空间，是研究量子纠缠 (Ch28) 与高维信号处理的必备工具

</div>

在处理涉及多个矩阵相互作用的复杂方程（如 $AX + XB = C$）时，传统的矩阵算术往往难以直接给出闭式解。**Kronecker 积**与 **Vec 算子**提供了一种强大的工具，能够将“矩阵方程”转化为标准的“向量方程”。这种降维打击的策略，使得我们能利用经典的线性系统理论解决高维的算子相互作用问题。

---

## 19.1 Kronecker 积的定义与性质

!!! definition "定义 19.1 (Kronecker 积)"
    设 $A$ 为 $m \times n$ 矩阵，$B$ 为 $p \times q$ 矩阵。它们的 **Kronecker 积** $A \otimes B$ 是一个 $mp \times nq$ 的分块矩阵：
    $$A \otimes B = \begin{pmatrix} a_{11}B & a_{12}B & \cdots & a_{1n}B \\ a_{21}B & a_{22}B & \cdots & a_{2n}B \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1}B & a_{m2}B & \cdots & a_{mn}B \end{pmatrix}$$

!!! theorem "定理 19.1 (混合乘积性质)"
    若矩阵维数匹配，满足：
    $$(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$$
    **意义**：这一性质允许我们将复杂的复合变换解构为两个低维空间的独立变换。

---

## 19.2 谱性质与迹

!!! theorem "定理 19.2 (特征值与迹)"
    1.  **特征值**：若 $A$ 的特征值为 $\{\lambda_i\}$，$B$ 的特征值为 $\{\mu_j\}$，则 $A \otimes B$ 的特征值为所有可能的乘积 $\{\lambda_i \mu_j\}$。
    2.  **迹**：$\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$。
    3.  **行列式**：$\det(A \otimes B) = (\det A)^p (\det B)^m$（其中 $B$ 为 $p$ 阶阵，$A$ 为 $m$ 阶阵）。

---

## 19.3 Vec 算子与向量化恒等式

!!! definition "定义 19.2 (Vec 算子)"
    对于 $m \times n$ 矩阵 $X$，$\operatorname{vec}(X)$ 是将 $X$ 的列按顺序堆叠成的一个 $mn \times 1$ 的向量。

!!! theorem "定理 19.3 (向量化恒等式)"
    对于适当维数的矩阵 $A, X, B$：
    $$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$$
    **意义**：这是矩阵方程求解的“金钥匙”，它成功地将未知矩阵 $X$ 从夹缝中剥离出来，转化为线性方程组的标准形式。

---

## 练习题

**1. [计算] 计算 $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} \otimes \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。**

??? success "参考答案"
    **计算步骤：**
    根据定义，将左边矩阵的每个元素乘以右边矩阵：
    $$A \otimes B = \begin{pmatrix} 1 \cdot B & 0 \cdot B \\ 0 \cdot B & 2 \cdot B \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 2 \\ 0 & 0 & 2 & 0 \end{pmatrix}$$

**2. [特征值] 若 $A$ 的特征值为 1 和 2，$B$ 的特征值为 3 和 4。求 $A \otimes B$ 的所有特征值。**

??? success "参考答案"
    **定理应用：**
    Kronecker 积的特征值是两个矩阵特征值的所有可能对的乘积。
    1. $1 \cdot 3 = 3$
    2. $1 \cdot 4 = 4$
    3. $2 \cdot 3 = 6$
    4. $2 \cdot 4 = 8$
    **结论**：特征值为 $\{3, 4, 6, 8\}$。

**3. [Kronecker和] 求上题中 $A \oplus B = A \otimes I + I \otimes B$ 的特征值。**

??? success "参考答案"
    **定理应用：**
    Kronecker 和的特征值是特征值对的**和**。
    1. $1 + 3 = 4$
    2. $1 + 4 = 5$
    3. $2 + 3 = 5$
    4. $2 + 4 = 6$
    **结论**：特征值为 $\{4, 5, 5, 6\}$。这常用于求解 $AX + XB = C$ 时的谱分析。

**4. [向量化] 将矩阵方程 $AX = B$ 转化为标准的 $My = f$ 形式（其中 $y = \operatorname{vec}(X)$）。**

??? success "参考答案"
    **推导：**
    1. 方程可以写为 $AXI = B$（其中 $I$ 是适当维数的单位阵）。
    2. 利用公式 $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$。
    3. 令 $B$（公式中的）$= I$，$X$（公式中的）$= X$，$A$（公式中的）$= A$。
    **结论**：$(I \otimes A) \operatorname{vec}(X) = \operatorname{vec}(B)$。

**5. [迹性质] 证明 $\operatorname{tr}(A \otimes B) = \operatorname{tr}(B \otimes A)$。**

??? success "参考答案"
    **证明：**
    1. $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$。
    2. 由于标量乘法满足交换律，$\operatorname{tr}(A)\operatorname{tr}(B) = \operatorname{tr}(B)\operatorname{tr}(A)$。
    3. $\operatorname{tr}(B \otimes A) = \operatorname{tr}(B)\operatorname{tr}(A)$。
    **结论**：虽然 $A \otimes B \neq B \otimes A$，但它们的迹是相等的。

**6. [逆矩阵] 若 $A, B$ 均可逆，求 $(A \otimes B)^{-1}$。**

??? success "参考答案"
    **利用混合乘积性质：**
    1. 设逆矩阵为 $A^{-1} \otimes B^{-1}$。
    2. 计算 $(A \otimes B)(A^{-1} \otimes B^{-1}) = (A A^{-1}) \otimes (B B^{-1})$。
    3. $= I \otimes I = I$。
    **结论**：$(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$。

**7. [秩] 证明 $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$。**

??? success "参考答案"
    **证明思路：**
    1. 利用 SVD：若 $A = U_1 \Sigma_1 V_1^*$ 且 $B = U_2 \Sigma_2 V_2^*$。
    2. 则 $A \otimes B = (U_1 \otimes U_2) (\Sigma_1 \otimes \Sigma_2) (V_1 \otimes V_2)^*$。
    3. 由于 $U_1 \otimes U_2$ 和 $V_1 \otimes V_2$ 仍然是酉矩阵，这构成了 $A \otimes B$ 的 SVD。
    4. 奇异值矩阵为 $\Sigma_1 \otimes \Sigma_2$，其非零元个数显然是 $\operatorname{rank}(A) \cdot \operatorname{rank}(B)$。

**8. [Vec运算] 对于 $2 \times 2$ 阵 $X = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$，写出 $\operatorname{vec}(X)$。**

??? success "参考答案"
    **步骤：**
    1. 提取第一列：$(a, c)^T$。
    2. 提取第二列：$(b, d)^T$。
    3. 纵向堆叠。
    **结论**：$\operatorname{vec}(X) = (a, c, b, d)^T$。注意：是按列堆叠，而非按行。

**9. [Lyapunov] 将 Lyapunov 方程 $AX + XA^T = Q$ 向量化。**

??? success "参考答案"
    **步骤：**
    1. 向量化每一项：$\operatorname{vec}(AXI) + \operatorname{vec}(IXA^T) = \operatorname{vec}(Q)$。
    2. 应用公式：$(I \otimes A) \operatorname{vec}(X) + (A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(Q)$。
    3. 提取公因子：$(I \otimes A + A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(Q)$。
    **结论**：$(A \oplus A) \operatorname{vec}(X) = \operatorname{vec}(Q)$。

**10. [应用] 为什么在量子力学中两个粒子的联合态要用张量积（Kronecker 积）表示？**

??? success "参考答案"
    **物理逻辑：**
    1. 每个粒子的状态由一个向量空间描述。
    2. 联合系统的自由度是子系统自由度的组合。
    3. 如果系统 1 有 $n$ 个基态，系统 2 有 $m$ 个基态，那么联合系统就有 $n \times m$ 个基态。
    **代数映射**：Kronecker 积正好通过一种“全排列”的方式构建了这样一个 $nm$ 维空间，能够完美刻画粒子间的独立性以及**纠缠态**（即不能分解为 $v_1 \otimes v_2$ 形式的向量）。

## 本章小结

Kronecker 积与 Vec 算子提供了矩阵代数的一套“升维与降维”方案：

1.  **维度的乘法**：Kronecker 积通过一种平铺嵌套的方式，将两个独立算子的作用整合为一个巨大的复合算子，是处理多体系统互作用的唯一语言。
2.  **算子的解构**：Vec 算子消除了矩阵的二维拓扑，将其转化为线性空间中最基本的向量形式，从而释放了经典线性方程组求解器的全部威力。
3.  **方程的统一**：矩阵方程向量化恒等式是连接矩阵方程理论与数值线性代数的纽带，它证明了所有的线性矩阵方程在本质上都是同一个线性系统在不同基下的表现。
