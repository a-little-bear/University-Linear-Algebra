# 第 28 章 线性代数在量子信息中的应用

<div class="context-flow" markdown>

**前置**：内积空间 (Ch08) · 特征值 (Ch06) · 酉矩阵 (Ch55A) · 张量积 (Ch21)

**本章脉络**：量子世界的代数公理 $\to$ 状态：单位范数向量（Bra-ket 记号） $\to$ 叠加原理与基 $\to$ 演化：酉变换 $U$ $\to$ 测量：厄米算子 $M$ 的特征值与投影 $\to$ 复合系统：张量积与量子纠缠 (Entanglement) $\to$ 纯态与混合态：密度矩阵 (Density Matrix) $\to$ 量子门电路（Hadamard, CNOT）的矩阵表示 $\to$ 应用：量子比特 (Qubits)、量子隐形传态、Shor 分解算法

**延伸**：量子信息是线性代数最前卫、最纯粹的应用场；它将物理实在完全等同于复 Hilbert 空间中的算子运算，证明了微观世界的本质是高维线性叠加，是理解下一代计算革命的唯一数学门票

</div>

在经典计算机中，信息以 0 或 1 存储。而在量子计算机中，信息以复向量的形式存在。**量子信息**（Quantum Information）完全建立在复线性代数的公理之上。通过利用**叠加**（线性组合）和**纠缠**（张量积的不可约性），量子系统展现出了超越经典逻辑的计算能力。本章将介绍如何用线性代数的语言描述这个充满概率与波动的微观世界。

---

## 28.1 量子态与 Bra-ket 记号

!!! definition "定义 28.1 (量子比特 Qubit)"
    量子比特是 $\mathbb{C}^2$ 空间中的单位向量：
    $$|\psi\rangle = \alpha |0\rangle + \beta |1\rangle, \quad |\alpha|^2 + |\beta|^2 = 1$$
    - $|0\rangle = (1, 0)^T, |1\rangle = (0, 1)^T$ 构成计算基。
    - $\langle \phi | \psi \rangle$ 表示复内积。

---

## 28.2 演化与测量

!!! theorem "定理 28.1 (量子演化)"
    量子系统的封闭演化由**酉矩阵** $U$ 描述：$|\psi'\rangle = U |\psi\rangle$。
    **物理意义**：酉性保证了概率的总和始终为 1。

!!! technique "公理：测量"
    对物理量 $M$（厄米阵）进行测量，所得结果必然是其特征值 $\lambda_i$。
    测量后，状态坍缩到对应的特征空间。

---

## 28.3 纠缠与张量积

!!! theorem "定理 28.2 (纠缠态判定)"
    对于复合系统 $H_A \otimes H_B$ 中的态 $|\Psi\rangle$，若它不能分解为单个系统态的张量积 $|\psi_A\rangle \otimes |\psi_B\rangle$，则称其为**纠缠态**。
    **代数本质**：纠缠态对应于秩大于 1 的张量。

---

## 练习题

**1. [基础] 判定量子比特 $|\psi\rangle = \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$ 是否归一化。**

??? success "参考答案"
    **验证：**
    计算范数的平方：$\| \psi \|^2 = |1/\sqrt{2}|^2 + |1/\sqrt{2}|^2 = 1/2 + 1/2 = 1$。
    **结论**：是的。该状态代表了 0 和 1 的**均匀叠加态**。

**2. [计算] 计算 Hadamard 门 $H = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ 作用于 $|0\rangle$ 的结果。**

??? success "参考答案"
    **矩阵乘法：**
    $H \begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \cdot 1 + 1 \cdot 0 \\ 1 \cdot 1 - 1 \cdot 0 \end{pmatrix} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \end{pmatrix}$。
    **结论**：Hadamard 门将基态 $|0\rangle$ 转化为了叠加态 $\frac{|0\rangle + |1\rangle}{\sqrt{2}}$。

**3. [测量概率] 对状态 $|\psi\rangle = \frac{\sqrt{3}}{2} |0\rangle + \frac{1}{2} |1\rangle$ 进行计算基测量，得到 1 的概率是多少？**

??? success "参考答案"
    **波恩定则 (Born Rule)：**
    概率 $P(1) = |\langle 1 | \psi \rangle|^2$。
    内积为 $1 \cdot (1/2) = 1/2$。
    **结论**：$P(1) = (1/2)^2 = 1/4$ 或 25%。

**4. [纠缠判定] 判定 Bell 态 $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$ 是否纠缠。**

??? success "参考答案"
    **分析：**
    1. 写成系数矩阵形式：$M = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。
    2. 计算矩阵的秩：$\operatorname{rank}(M) = 2$。
    3. 根据张量积性质，可分态（非纠缠）的系数矩阵秩必为 1。
    **结论**：由于秩大于 1，这是一个**最大纠缠态**。

**5. [密度矩阵] 计算纯态 $|\psi\rangle$ 的密度矩阵 $\rho$。**

??? success "参考答案"
    **公式：**
    $\rho = |\psi\rangle \langle \psi|$（外积）。
    **性质**：纯态密度矩阵满足 $\rho^2 = \rho$ 且 $\operatorname{tr}(\rho) = 1$。

**6. [酉性验证] 验证 Pauli-X 门 $X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 是酉矩阵。**

??? success "参考答案"
    **计算：**
    $X^* X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$。
    **结论**：是酉矩阵。在量子逻辑中，它对应于经典电路的非门 (NOT)。

**7. [测量影响] 为什么测量算子必须是厄米的？**

??? success "参考答案"
    **物理理由：**
    测量结果（观测值）必须是**实数**。
    在线性代数中，保证所有特征值均为实数的算子类正是**厄米算子**（Hermitian）。

**8. [计算] 两个 Qubit 系统的复合基 $|0\rangle \otimes |1\rangle$ 写成向量是什么？**

??? success "参考答案"
    **利用 Kronecker 积：**
    $\begin{pmatrix} 1 \\ 0 \end{pmatrix} \otimes \begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 1 \cdot 0 \\ 1 \cdot 1 \\ 0 \cdot 0 \\ 0 \cdot 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}$。
    这是 4 维空间中的第二个标准基向量，常记为 $|01\rangle$。

**9. [不确定性] 简述为什么 $[A, B] \neq 0$ 导致无法同时测量。**

??? success "参考答案"
    **代数解释：**
    如果不交换，则 $A$ 和 $B$ 没有共同的特征基（见 Ch63A）。
    测量 $A$ 会使状态坍缩到 $A$ 的特征空间，但这通常不是 $B$ 的特征空间。因此，测量 $A$ 之后立即测量 $B$ 会引入新的随机性，使得两个物理量无法同时具有确定的值。

**10. [应用] 什么是量子比特的“布洛赫球”（Bloch Sphere）表示？**

??? success "参考答案"
    由于量子比特是一个模为 1 的 2 维复向量，且整体相位无意义。
    其状态可以唯一映射到 3 维空间的单位球面上：
    $|\psi\rangle = \cos(\theta/2)|0\rangle + e^{i\phi}\sin(\theta/2)|1\rangle$。
    这证明了单一量子比特的代数结构与 3 维球面的几何旋转是同构的。

## 本章小结

线性代数是描述量子实在的唯一语法：

1.  **叠加的代数本质**：量子世界的并行性源于向量空间的线性组合，通过基变换，我们实现了不同物理属性之间的灵活转换。
2.  **纠缠的张量特性**：纠缠态揭示了张量空间中不可约元素的强大关联，证明了整体信息量可以超越部分信息的简单加和。
3.  **酉性的守恒律**：酉矩阵演化确立了量子计算的逻辑可逆性与概率守恒，是设计量子算法与纠错协议的最高行为准则。
