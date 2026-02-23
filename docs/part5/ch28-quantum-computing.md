# 第 28 章 线性代数在量子信息中的应用

<div class="context-flow" markdown>

**前置**：内积空间与酉变换(Ch7-8) · 张量积(Ch21) · 正定矩阵(Ch16) · SVD(Ch11) · 矩阵不等式(Ch18) · Kronecker 积(Ch19)

**脉络**：量子态(Hilbert 空间中的单位向量) → 量子门(酉矩阵) → 张量积(多体系统) → 纠缠(Bell 态/Schmidt 分解) → 测量(投影/POVM) → 密度矩阵与量子信道(Kraus 算子) → 量子纠错(稳定子码) → 矩阵不等式(von Neumann 熵/强次可加性)

</div>

量子信息科学完全建立在线性代数的基础之上。量子态是 Hilbert 空间中的向量，量子演化是酉变换，量子测量是投影运算，量子纠缠是张量积空间中的非可分结构，量子信道是完全正保迹映射。本章从量子态与 Hilbert 空间出发，系统展示线性代数工具在量子计算、量子通信和量子纠错中的核心作用，并以量子信息中的矩阵不等式作为结束。

---

## 28.1 量子态与 Hilbert 空间

<div class="context-flow" markdown>

**基本对应**：经典比特 $\{0, 1\}$ → 量子比特 $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle \in \mathbb{C}^2$ → Dirac 记号 → 纯态 vs 混合态

**链接**：Ch8 内积空间的直接应用

</div>

量子比特（qubit）是量子信息的基本单元，其数学描述是二维复 Hilbert 空间中的单位向量。Dirac 记号为量子力学提供了一种简洁而强大的线性代数语言。

!!! definition "定义 28.1 (量子态)"
    一个**量子态**（quantum state）是 Hilbert 空间 $\mathcal{H}$ 中的单位向量 $|\psi\rangle$，满足 $\langle\psi|\psi\rangle = 1$。在 Dirac 记号中：

    - **右矢**（ket）$|\psi\rangle$ 表示列向量。
    - **左矢**（bra）$\langle\psi|$ 表示 $|\psi\rangle$ 的共轭转置，即行向量。
    - **内积** $\langle\phi|\psi\rangle$ 为复数。
    - **外积** $|\psi\rangle\langle\phi|$ 为矩阵（线性算子）。

    全局相位 $e^{i\theta}|\psi\rangle$ 与 $|\psi\rangle$ 表示相同的物理状态，因此量子态的物理空间实际上是射影 Hilbert 空间。

!!! definition "定义 28.2 (量子比特)"
    一个**量子比特**（qubit）的状态为 $\mathbb{C}^2$ 中的单位向量

    $$
    |\psi\rangle = \alpha|0\rangle + \beta|1\rangle, \quad |\alpha|^2 + |\beta|^2 = 1,
    $$

    其中 $|0\rangle = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$，$|1\rangle = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ 为**计算基**（computational basis），$\alpha, \beta \in \mathbb{C}$。$|\alpha|^2$ 和 $|\beta|^2$ 分别是测量得到 $|0\rangle$ 和 $|1\rangle$ 的概率。

    $n$ 个量子比特的联合状态空间为张量积 $\mathcal{H} = (\mathbb{C}^2)^{\otimes n} \cong \mathbb{C}^{2^n}$，维数随比特数指数增长。

!!! theorem "定理 28.1 (密度矩阵的刻画)"
    算子 $\rho$ 是一个合法的密度矩阵（表示量子态）当且仅当满足以下两个条件：

    1. $\rho \succeq 0$（半正定）。
    2. $\operatorname{tr}(\rho) = 1$。

    进一步地：

    - $\rho$ 表示**纯态** $\iff$ $\rho^2 = \rho$（幂等）$\iff$ $\operatorname{rank}(\rho) = 1$ $\iff$ $\operatorname{tr}(\rho^2) = 1$。
    - $\rho$ 表示**混合态** $\iff$ $\operatorname{tr}(\rho^2) < 1$，此时 $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$（$p_i > 0$，$\sum p_i = 1$）。

??? proof "证明"
    **必要性**：纯态 $\rho = |\psi\rangle\langle\psi|$ 显然半正定（对任意 $|\phi\rangle$，$\langle\phi|\rho|\phi\rangle = |\langle\phi|\psi\rangle|^2 \ge 0$），且 $\operatorname{tr}(\rho) = \langle\psi|\psi\rangle = 1$。混合态 $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ 是半正定矩阵的凸组合，仍半正定，且 $\operatorname{tr}(\rho) = \sum_i p_i = 1$。

    **充分性**：若 $\rho \succeq 0$ 且 $\operatorname{tr}(\rho) = 1$，则 $\rho$ 有谱分解 $\rho = \sum_i \lambda_i |e_i\rangle\langle e_i|$，其中 $\lambda_i \ge 0$（半正定），$\sum_i \lambda_i = 1$（迹为 $1$）。因此 $\lambda_i$ 构成概率分布，$\rho$ 确实是纯态的概率混合。

    **纯态判据**：若 $\rho = |\psi\rangle\langle\psi|$，则 $\rho^2 = |\psi\rangle\langle\psi|\psi\rangle\langle\psi| = |\psi\rangle\langle\psi| = \rho$。反之，设 $\operatorname{tr}(\rho^2) = 1$。由 $\rho$ 的特征值 $\lambda_i \ge 0$、$\sum \lambda_i = 1$、$\sum \lambda_i^2 = 1$，以及 Cauchy-Schwarz 不等式 $\sum \lambda_i^2 \le (\max \lambda_i)\sum \lambda_i = \max \lambda_i$，得 $\max \lambda_i \ge 1$。又 $\lambda_i \le \sum \lambda_i = 1$，故恰有一个 $\lambda_i = 1$，其余为 $0$。$\blacksquare$

!!! example "例 28.1"
    **Bloch 球表示。** 单量子比特纯态（去除全局相位后）可参数化为

    $$
    |\psi\rangle = \cos\frac{\theta}{2}|0\rangle + e^{i\phi}\sin\frac{\theta}{2}|1\rangle, \quad \theta \in [0, \pi],\, \phi \in [0, 2\pi).
    $$

    对应密度矩阵 $\rho = |\psi\rangle\langle\psi| = \frac{1}{2}(I + \mathbf{r} \cdot \boldsymbol{\sigma})$，其中 Bloch 向量 $\mathbf{r} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$，$|\mathbf{r}| = 1$。

    混合态对应 $|\mathbf{r}| < 1$（球内部）。最大混合态 $\rho = I/2$ 对应 $\mathbf{r} = \mathbf{0}$（球心）。

    | 态 | Bloch 向量 | 位置 |
    |:---:|:---:|:---:|
    | $\|0\rangle$ | $(0, 0, 1)$ | 北极 |
    | $\|1\rangle$ | $(0, 0, -1)$ | 南极 |
    | $\|+\rangle = \frac{1}{\sqrt{2}}(\|0\rangle + \|1\rangle)$ | $(1, 0, 0)$ | $x$ 轴正向 |
    | $\|-\rangle = \frac{1}{\sqrt{2}}(\|0\rangle - \|1\rangle)$ | $(-1, 0, 0)$ | $x$ 轴负向 |

---

## 28.2 量子门与酉变换

<div class="context-flow" markdown>

**核心原理**：量子演化 = 酉变换 $U^\dagger U = I$ → 保持内积（概率守恒）→ 单比特门（Pauli, Hadamard, S, T）→ 多比特门（CNOT, Toffoli）→ 通用门集

**链接**：Ch7 酉矩阵理论的直接应用

</div>

封闭量子系统的演化由酉矩阵描述。量子门是量子计算的基本操作单元。

!!! definition "定义 28.3 (量子门)"
    **量子门**（quantum gate）是作用在量子比特上的酉矩阵 $U$（$U^\dagger U = UU^\dagger = I$）。

    **单比特 Pauli 门**：

    $$
    X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\  Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix},\  Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
    $$

    **Hadamard 门**：$H = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，满足 $H|0\rangle = |+\rangle$，$H|1\rangle = |-\rangle$。

    **相位门**：$S = \begin{pmatrix} 1 & 0 \\ 0 & i \end{pmatrix}$，$T = \begin{pmatrix} 1 & 0 \\ 0 & e^{i\pi/4} \end{pmatrix}$。注意 $S = T^2$，$Z = S^2$。

!!! definition "定义 28.4 (多比特门)"
    **CNOT（受控非）门**是最重要的双比特门：

    $$
    \text{CNOT} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix} = |0\rangle\langle 0| \otimes I + |1\rangle\langle 1| \otimes X.
    $$

    CNOT 以第一个比特（控制位）控制第二个比特（目标位）的翻转。

    **Toffoli 门**（CCNOT）是三比特门：当且仅当两个控制位均为 $|1\rangle$ 时翻转目标位。Toffoli 门是经典可逆计算的通用门。

!!! theorem "定理 28.2 (量子门集的通用性)"
    门集 $\{H, T, \text{CNOT}\}$ 是**通用的**：对任意 $n$ 比特酉矩阵 $U \in U(2^n)$ 和任意 $\varepsilon > 0$，存在由 $H$、$T$、$\text{CNOT}$ 组成的有限量子电路 $\tilde{U}$，使得

    $$
    \|U - \tilde{U}\| < \varepsilon.
    $$

    更精确地，所需门数为 $O(4^n \cdot \operatorname{poly}(\log(1/\varepsilon)))$。

??? proof "证明"
    证明分三步：

    **第 1 步**：任意 $n$ 比特酉矩阵可分解为双比特门的乘积。通过 QR 分解的量子类比（Givens 旋转），$U(2^n)$ 中的任意酉矩阵可以写成 $O(4^n)$ 个作用在两比特上的酉门的乘积。

    **第 2 步**：任意双比特门可用 CNOT 和单比特门实现。这通过 KAK 分解（Cartan 分解）实现：任意 $SU(4)$ 元素可以写成 $(A_1 \otimes A_2) \cdot \exp(i(a X \otimes X + b Y \otimes Y + c Z \otimes Z)) \cdot (A_3 \otimes A_4)$ 的形式，其中 $A_i \in SU(2)$。中间的纠缠部分可用至多 $3$ 个 CNOT 实现。

    **第 3 步**：$\{H, T\}$ 在 $SU(2)$ 中稠密（Solovay-Kitaev 定理）。$H$ 和 $T$ 生成的群在 $SU(2)$ 中稠密，且任何 $SU(2)$ 元素可以用 $O(\log^c(1/\varepsilon))$ 个 $H$ 和 $T$ 门近似到精度 $\varepsilon$，其中 $c \approx 3.97$。

    综合三步，总门数为 $O(4^n \cdot \operatorname{poly}(\log(1/\varepsilon)))$。$\blacksquare$

!!! example "例 28.2"
    **Bell 态的制备电路。** 利用 Hadamard 门和 CNOT 门从 $|00\rangle$ 制备 Bell 态：

    $$
    |00\rangle \xrightarrow{H \otimes I} \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)|0\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |10\rangle) \xrightarrow{\text{CNOT}} \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle) = |\Phi^+\rangle.
    $$

    矩阵表示为

    $$
    \text{CNOT} \cdot (H \otimes I) = \begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&0&1\\0&0&1&0\end{pmatrix} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1&0&1&0\\0&1&0&1\\1&0&-1&0\\0&1&0&-1\end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix}1&0&1&0\\0&1&0&1\\0&1&0&-1\\1&0&-1&0\end{pmatrix}.
    $$

    作用在 $|00\rangle = (1,0,0,0)^T$ 上得 $\frac{1}{\sqrt{2}}(1,0,0,1)^T = |\Phi^+\rangle$。

---

## 28.3 张量积与多体系统

<div class="context-flow" markdown>

**核心**：$\mathcal{H}_A \otimes \mathcal{H}_B$ = 复合系统的状态空间 → Kronecker 积给出矩阵表示 → 维数指数增长 = 量子优势的根源

**链接**：Ch19 Kronecker 积和 Ch21 张量代数的直接应用

</div>

量子多体系统的数学描述依赖于张量积。复合系统的状态空间是各子系统状态空间的张量积，这一结构导致了指数级的维度增长。

!!! definition "定义 28.5 (复合量子系统)"
    由子系统 $A$（状态空间 $\mathcal{H}_A$，$\dim = d_A$）和 $B$（状态空间 $\mathcal{H}_B$，$\dim = d_B$）组成的**复合量子系统**的状态空间为

    $$
    \mathcal{H}_{AB} = \mathcal{H}_A \otimes \mathcal{H}_B, \quad \dim \mathcal{H}_{AB} = d_A \cdot d_B.
    $$

    若 $\{|i\rangle_A\}$ 和 $\{|j\rangle_B\}$ 分别为 $\mathcal{H}_A$ 和 $\mathcal{H}_B$ 的正交基，则 $\{|i\rangle_A \otimes |j\rangle_B\}$ 为 $\mathcal{H}_{AB}$ 的正交基。

    对于算子，$A$ 上的算子 $O_A$ 和 $B$ 上的算子 $O_B$ 的联合作用为 Kronecker 积 $O_A \otimes O_B$。

!!! theorem "定理 28.3 (张量积的性质)"
    设 $A \in \mathbb{C}^{m \times m}$，$B \in \mathbb{C}^{n \times n}$，则 Kronecker 积 $A \otimes B \in \mathbb{C}^{mn \times mn}$ 满足：

    1. **混合乘积性质**：$(A \otimes B)(C \otimes D) = AC \otimes BD$。
    2. **谱性质**：若 $A$ 的特征值为 $\{\alpha_i\}$，$B$ 的特征值为 $\{\beta_j\}$，则 $A \otimes B$ 的特征值为 $\{\alpha_i \beta_j\}$。
    3. **迹**：$\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B)$。
    4. **行列式**：$\det(A \otimes B) = (\det A)^n (\det B)^m$。
    5. **酉性保持**：若 $A$ 和 $B$ 都是酉矩阵，则 $A \otimes B$ 也是酉矩阵。

??? proof "证明"
    (1)：$(A \otimes B)(C \otimes D)$ 的 $(i,j)$-块（$A$ 的 $(i,k)$-元素和 $C$ 的 $(k,j)$-元素）为 $\sum_k A_{ik}C_{kj} \cdot BD = (AC)_{ij} BD$，即 $AC \otimes BD$ 的 $(i,j)$-块。

    (2)：若 $A\mathbf{u} = \alpha\mathbf{u}$，$B\mathbf{v} = \beta\mathbf{v}$，则 $(A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = A\mathbf{u} \otimes B\mathbf{v} = \alpha\beta(\mathbf{u} \otimes \mathbf{v})$。

    (3)：$\operatorname{tr}(A \otimes B) = \sum_{i,j}(A \otimes B)_{(i,j),(i,j)} = \sum_i A_{ii} \sum_j B_{jj} = \operatorname{tr}(A)\operatorname{tr}(B)$。

    (5)：$(A \otimes B)^\dagger(A \otimes B) = (A^\dagger \otimes B^\dagger)(A \otimes B) = A^\dagger A \otimes B^\dagger B = I \otimes I = I$。$\blacksquare$

!!! example "例 28.3"
    **两比特系统的矩阵表示。** 考虑两比特态 $|\psi\rangle = \frac{1}{\sqrt{2}}|00\rangle + \frac{1}{2}|01\rangle + \frac{1}{2}|10\rangle$。

    向量表示：$|\psi\rangle = \frac{1}{\sqrt{2}}(1, 0)^T \otimes (1, 0)^T + \frac{1}{2}(1, 0)^T \otimes (0, 1)^T + \frac{1}{2}(0, 1)^T \otimes (1, 0)^T = (\frac{1}{\sqrt{2}}, \frac{1}{2}, \frac{1}{2}, 0)^T$。

    验证归一化：$\frac{1}{2} + \frac{1}{4} + \frac{1}{4} + 0 = 1$。

    对此态施加 $X \otimes Z$：

    $$
    (X \otimes Z)|\psi\rangle = \begin{pmatrix}0&0&1&0\\0&0&0&-1\\1&0&0&0\\0&-1&0&0\end{pmatrix}\begin{pmatrix}1/\sqrt{2}\\1/2\\1/2\\0\end{pmatrix} = \begin{pmatrix}1/2\\0\\1/\sqrt{2}\\-1/2\end{pmatrix}.
    $$

---

## 28.4 纠缠与 Bell 态

<div class="context-flow" markdown>

**核心概念**：乘积态 $|\alpha\rangle \otimes |\beta\rangle$ vs 纠缠态（不可分）→ Bell 态（最大纠缠）→ Schmidt 分解 = SVD 的量子版本 → Schmidt 秩判定纠缠

**链接**：Ch11 SVD 和 Ch21 张量积的核心应用

</div>

量子纠缠是量子信息与经典信息最本质的区别，其数学本质是张量积空间中的不可分结构。

!!! definition "定义 28.6 (纠缠态与可分态)"
    双系统纯态 $|\psi\rangle_{AB} \in \mathcal{H}_A \otimes \mathcal{H}_B$ 称为**可分**（separable），若存在 $|\alpha\rangle \in \mathcal{H}_A$ 和 $|\beta\rangle \in \mathcal{H}_B$ 使得

    $$
    |\psi\rangle_{AB} = |\alpha\rangle \otimes |\beta\rangle.
    $$

    否则称为**纠缠态**（entangled state）。

    对混合态，$\rho_{AB}$ 可分当且仅当存在分解 $\rho_{AB} = \sum_i p_i \rho_A^{(i)} \otimes \rho_B^{(i)}$（$p_i > 0$，$\sum p_i = 1$）。

!!! definition "定义 28.7 (Bell 态)"
    四个 **Bell 态**构成 $\mathbb{C}^2 \otimes \mathbb{C}^2$ 的正交基：

    $$
    |\Phi^{\pm}\rangle = \frac{1}{\sqrt{2}}(|00\rangle \pm |11\rangle), \quad |\Psi^{\pm}\rangle = \frac{1}{\sqrt{2}}(|01\rangle \pm |10\rangle).
    $$

    Bell 态是**最大纠缠态**：每个子系统的约化密度矩阵为 $\rho_A = \rho_B = \frac{I}{2}$（最大混合态），von Neumann 熵取最大值 $S = \log 2 = 1$ bit。

!!! theorem "定理 28.4 (Schmidt 分解定理)"
    任意双系统纯态 $|\psi\rangle \in \mathcal{H}_A \otimes \mathcal{H}_B$（$\dim \mathcal{H}_A = d_A$，$\dim \mathcal{H}_B = d_B$）可以写为

    $$
    |\psi\rangle = \sum_{i=1}^{r} \sqrt{\lambda_i}\, |u_i\rangle_A \otimes |v_i\rangle_B,
    $$

    其中 $r \le \min(d_A, d_B)$ 为 **Schmidt 秩**，$\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_r > 0$ 为 **Schmidt 系数**（$\sum_{i=1}^r \lambda_i = 1$），$\{|u_i\rangle\}$ 和 $\{|v_i\rangle\}$ 分别为 $\mathcal{H}_A$ 和 $\mathcal{H}_B$ 中的正交集。Schmidt 系数唯一确定，当 Schmidt 系数互不相同时向量也唯一（到相位因子）。

    $|\psi\rangle$ 为乘积态 $\iff$ Schmidt 秩 $r = 1$。

??? proof "证明"
    将 $|\psi\rangle = \sum_{j,k} c_{jk} |j\rangle_A |k\rangle_B$ 的系数排列为矩阵 $C \in \mathbb{C}^{d_A \times d_B}$，其中 $C_{jk} = c_{jk}$。对 $C$ 进行奇异值分解（SVD，Ch11）：$C = U\Sigma V^\dagger$，其中 $U \in \mathbb{C}^{d_A \times d_A}$ 和 $V \in \mathbb{C}^{d_B \times d_B}$ 为酉矩阵，$\Sigma$ 为对角矩阵，对角元素 $\sigma_i \ge 0$。

    则 $c_{jk} = \sum_{i=1}^{r} \sigma_i U_{ji} V_{ki}^*$。定义

    $$
    |u_i\rangle_A = \sum_j U_{ji} |j\rangle_A, \quad |v_i\rangle_B = \sum_k V_{ki}^* |k\rangle_B.
    $$

    由 $U$ 和 $V$ 的酉性，$\{|u_i\rangle\}$ 和 $\{|v_i\rangle\}$ 分别正交。归一化 $\sum_i \sigma_i^2 = \|C\|_F^2 = \langle\psi|\psi\rangle = 1$，令 $\lambda_i = \sigma_i^2$，得

    $$
    |\psi\rangle = \sum_{i=1}^r \sigma_i |u_i\rangle_A |v_i\rangle_B = \sum_{i=1}^r \sqrt{\lambda_i} |u_i\rangle_A |v_i\rangle_B.
    $$

    唯一性来自 SVD 奇异值的唯一性。$\blacksquare$

!!! example "例 28.4"
    **Bell 态的 Schmidt 分解与纠缠判定。**

    $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$。系数矩阵 $C = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = \frac{1}{\sqrt{2}}I$。SVD：$U = V = I$，$\sigma_1 = \sigma_2 = \frac{1}{\sqrt{2}}$。Schmidt 秩 $r = 2 > 1$，因此 $|\Phi^+\rangle$ 是纠缠态。Schmidt 系数 $\lambda_1 = \lambda_2 = \frac{1}{2}$（均匀），确认最大纠缠。

    比较：$|\psi\rangle = |01\rangle$。系数矩阵 $C = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。SVD 给出 $\sigma_1 = 1$，Schmidt 秩 $r = 1$，因此 $|01\rangle = |0\rangle \otimes |1\rangle$ 是乘积态（可分）。

---

## 28.5 量子测量与 POVM

<div class="context-flow" markdown>

**投影测量**：正交分解 $I = \sum_m P_m$ → Born 规则 $p(m) = \langle\psi|P_m|\psi\rangle$ → 测后态塌缩 → POVM 推广到非正交测量 → Naimark 膨胀定理

**链接**：Ch7 正交投影和谱分解的直接应用

</div>

量子测量将线性代数的投影算子理论与概率论结合。

!!! definition "定义 28.8 (投影测量与 Born 规则)"
    一个**投影测量**由一组正交投影算子 $\{P_m\}$ 定义，满足

    $$
    P_m P_{m'} = \delta_{mm'} P_m, \quad \sum_m P_m = I, \quad P_m = P_m^\dagger.
    $$

    **Born 规则**：对态 $|\psi\rangle$，测量结果为 $m$ 的概率为

    $$
    p(m) = \langle\psi|P_m|\psi\rangle = \operatorname{tr}(P_m |\psi\rangle\langle\psi|).
    $$

    测量后态塌缩为 $\frac{P_m|\psi\rangle}{\sqrt{p(m)}}$。对密度矩阵 $\rho$，$p(m) = \operatorname{tr}(P_m \rho)$。

!!! definition "定义 28.9 (POVM)"
    **正算子值测度**（Positive Operator-Valued Measure, POVM）是一般化的测量框架。POVM 由一组**正算子** $\{E_m\}$ 定义，满足

    $$
    E_m \succeq 0, \quad \sum_m E_m = I.
    $$

    测量结果为 $m$ 的概率为 $p(m) = \operatorname{tr}(E_m \rho)$。POVM 元素 $E_m$ 不必是投影，因此可以有比维数更多的测量结果。POVM 描述了信息提取，但不指定测后态。

!!! theorem "定理 28.5 (Naimark 膨胀定理)"
    任何作用在 $\mathcal{H}$ 上的 POVM $\{E_m\}_{m=1}^M$ 都可以实现为扩展空间 $\mathcal{H} \otimes \mathcal{H}_{\text{aux}}$ 上的投影测量。即存在辅助空间 $\mathcal{H}_{\text{aux}}$、初始辅助态 $|0\rangle_{\text{aux}}$ 和扩展空间上的正交投影 $\{\Pi_m\}$ 使得

    $$
    E_m = \langle 0|_{\text{aux}} \Pi_m |0\rangle_{\text{aux}}, \quad \forall m.
    $$

    等价地，$\operatorname{tr}(E_m \rho) = \operatorname{tr}(\Pi_m (\rho \otimes |0\rangle\langle 0|))$。

??? proof "证明"
    构造 $\mathcal{H}_{\text{aux}} = \mathbb{C}^M$。对每个 POVM 元素 $E_m$，由于 $E_m \succeq 0$，存在 $A_m$ 使得 $E_m = A_m^\dagger A_m$（取 $A_m = \sqrt{E_m}$）。定义等距映射 $V: \mathcal{H} \to \mathcal{H} \otimes \mathcal{H}_{\text{aux}}$，

    $$
    V|\psi\rangle = \sum_{m=1}^M (A_m|\psi\rangle) \otimes |m\rangle.
    $$

    验证 $V$ 是等距映射：$\langle\psi|V^\dagger V|\psi\rangle = \sum_m \langle\psi|A_m^\dagger A_m|\psi\rangle = \langle\psi|\sum_m E_m|\psi\rangle = \langle\psi|\psi\rangle$。

    定义 $\Pi_m = I_{\mathcal{H}} \otimes |m\rangle\langle m|$（扩展空间上的投影）。则

    $$
    \langle\psi|V^\dagger \Pi_m V|\psi\rangle = \langle\psi|A_m^\dagger A_m|\psi\rangle = \langle\psi|E_m|\psi\rangle = p(m).
    $$

    因此 POVM 测量等价于先嵌入扩展空间再做投影测量。$\blacksquare$

!!! example "例 28.5"
    **三元 POVM 区分非正交态。** 考虑以等概率 $1/2$ 给出 $|0\rangle$ 或 $|+\rangle = \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$。由于 $\langle 0|+\rangle = \frac{1}{\sqrt{2}} \ne 0$，两态不正交，不能被投影测量完美区分。

    构造无歧义态区分（USD）的三元 POVM $\{E_0, E_+, E_?\}$：

    - $E_0 = \frac{1}{1 + 1/\sqrt{2}} |1\rangle\langle 1|$：当触发时确定态为 $|+\rangle$（因为 $\langle 1|0\rangle = 0$，故 $\operatorname{tr}(E_0 |0\rangle\langle 0|) = 0$）。
    - $E_+ = \frac{1}{1 + 1/\sqrt{2}} |-\rangle\langle -|$：当触发时确定态为 $|0\rangle$（因为 $\langle -|+\rangle = 0$）。
    - $E_? = I - E_0 - E_+$：不确定结果。

    最优无歧义区分成功概率为 $p_{\text{succ}} = 1 - |\langle 0|+\rangle| = 1 - \frac{1}{\sqrt{2}} \approx 0.293$。

---

## 28.6 密度矩阵与量子信道

<div class="context-flow" markdown>

**开放系统**：混合态 $\rho$ → 偏迹(约化密度矩阵) → 量子信道 $\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger$（Kraus 表示）→ CPTP 映射 → Choi-Kraus 表示定理

**链接**：Ch16 正定矩阵/半正定锥和 Ch21 张量积的推广

</div>

量子信道描述开放量子系统的演化，其数学本质是完全正保迹（CPTP）映射。

!!! definition "定义 28.10 (偏迹与约化密度矩阵)"
    对复合系统态 $\rho_{AB} \in \mathcal{B}(\mathcal{H}_A \otimes \mathcal{H}_B)$，子系统 $A$ 的**约化密度矩阵**为

    $$
    \rho_A = \operatorname{tr}_B(\rho_{AB}) = \sum_j (I_A \otimes \langle j|_B) \rho_{AB} (I_A \otimes |j\rangle_B),
    $$

    其中 $\{|j\rangle_B\}$ 为 $\mathcal{H}_B$ 的任意正交基。偏迹是唯一满足以下性质的运算：对所有 $A$ 上的算子 $O_A$，

    $$
    \operatorname{tr}(O_A \rho_A) = \operatorname{tr}((O_A \otimes I_B)\rho_{AB}).
    $$

!!! definition "定义 28.11 (量子信道与 Kraus 表示)"
    **量子信道**是完全正保迹（CPTP）映射 $\mathcal{E}: \mathcal{B}(\mathcal{H}_{\text{in}}) \to \mathcal{B}(\mathcal{H}_{\text{out}})$。其 **Kraus 表示**为

    $$
    \mathcal{E}(\rho) = \sum_{k=1}^{r} K_k \rho K_k^\dagger, \quad \sum_{k=1}^{r} K_k^\dagger K_k = I.
    $$

    $\{K_k\}$ 称为 **Kraus 算子**。条件 $\sum_k K_k^\dagger K_k = I$ 保证迹守恒（$\operatorname{tr}(\mathcal{E}(\rho)) = \operatorname{tr}(\rho) = 1$）。

    常见量子信道包括：

    - **去极化信道**：$\mathcal{E}(\rho) = (1-p)\rho + \frac{p}{d}I$，以概率 $p$ 将态替换为最大混合态。
    - **振幅阻尼信道**：$K_0 = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-\gamma} \end{pmatrix}$，$K_1 = \begin{pmatrix} 0 & \sqrt{\gamma} \\ 0 & 0 \end{pmatrix}$，模拟自发辐射。

!!! theorem "定理 28.6 (Choi-Kraus 表示定理)"
    线性映射 $\mathcal{E}: \mathcal{B}(\mathcal{H}_{\text{in}}) \to \mathcal{B}(\mathcal{H}_{\text{out}})$ 是完全正映射当且仅当存在算子 $\{K_k\}_{k=1}^r$ 使得

    $$
    \mathcal{E}(\rho) = \sum_{k=1}^r K_k \rho K_k^\dagger.
    $$

    等价地，$\mathcal{E}$ 是完全正的当且仅当其 **Choi 矩阵**

    $$
    J(\mathcal{E}) = \sum_{i,j} |i\rangle\langle j| \otimes \mathcal{E}(|i\rangle\langle j|) \in \mathcal{B}(\mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}})
    $$

    是半正定的。$\mathcal{E}$ 还是保迹的当且仅当 $\operatorname{tr}_{\text{out}}(J(\mathcal{E})) = I_{\text{in}}$。

??? proof "证明"
    **Kraus $\Rightarrow$ CP**：对任意扩展空间 $\mathcal{H}_R$ 和 $\rho_{RA} \succeq 0$，

    $$
    (\operatorname{id}_R \otimes \mathcal{E})(\rho_{RA}) = \sum_k (I_R \otimes K_k) \rho_{RA} (I_R \otimes K_k)^\dagger \succeq 0.
    $$

    每一项 $X\rho X^\dagger \succeq 0$（保持半正定性），故和也半正定。

    **CP $\Rightarrow$ Kraus**：设 $\mathcal{E}$ 完全正，则 Choi 矩阵 $J(\mathcal{E}) \succeq 0$。对 $J(\mathcal{E})$ 做谱分解 $J(\mathcal{E}) = \sum_k |\phi_k\rangle\langle\phi_k|$，将 $|\phi_k\rangle \in \mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}}$ 重塑为算子 $K_k: \mathcal{H}_{\text{in}} \to \mathcal{H}_{\text{out}}$（通过自然同构 $\mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}} \cong \operatorname{Hom}(\mathcal{H}_{\text{in}}, \mathcal{H}_{\text{out}})$）。则可以验证 $\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger$。

    保迹条件 $\sum_k K_k^\dagger K_k = I$ 等价于 $\operatorname{tr}_{\text{out}}(J(\mathcal{E})) = I_{\text{in}}$，这可通过直接计算验证。$\blacksquare$

!!! example "例 28.6"
    **退相干信道。** 单比特退相干（dephasing）信道 $\mathcal{E}(\rho) = (1-p)\rho + pZ\rho Z$，Kraus 算子为 $K_0 = \sqrt{1-p}\, I$，$K_1 = \sqrt{p}\, Z$。

    作用在一般态 $\rho = \begin{pmatrix} a & b \\ b^* & 1-a \end{pmatrix}$ 上：

    $$
    \mathcal{E}(\rho) = (1-p)\begin{pmatrix} a & b \\ b^* & 1-a \end{pmatrix} + p\begin{pmatrix} a & -b \\ -b^* & 1-a \end{pmatrix} = \begin{pmatrix} a & (1-2p)b \\ (1-2p)b^* & 1-a \end{pmatrix}.
    $$

    对角元不变，非对角元衰减为 $(1-2p)$ 倍。当 $p = 1/2$ 时，非对角元消失，相干性完全丧失。Choi 矩阵为

    $$
    J(\mathcal{E}) = \begin{pmatrix} 1 & 0 & 0 & 1-2p \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 1-2p & 0 & 0 & 1 \end{pmatrix},
    $$

    特征值为 $\{1 + (1-2p), 1 - (1-2p), 0, 0\} = \{2-2p, 2p, 0, 0\} \ge 0$，确认完全正性。

---

## 28.7 量子纠错

<div class="context-flow" markdown>

**核心思想**：逻辑比特编码到物理比特的子空间 → 错误 = 线性算子 → Knill-Laflamme 条件 = 码空间上的正交性条件 → 稳定子码 → CSS 码

**链接**：Ch4 子空间和 Ch7 投影的核心应用

</div>

量子纠错是容错量子计算的基础。

!!! definition "定义 28.12 (量子纠错码)"
    一个 $[[n, k, d]]$ **量子纠错码**将 $k$ 个逻辑量子比特编码到 $n$ 个物理量子比特中。码空间 $\mathcal{C}$ 是 $(\mathbb{C}^2)^{\otimes n}$ 的 $2^k$ 维子空间。码距 $d$ 是可以检测的最大错误权重加 $1$；码可纠正 $\lfloor(d-1)/2\rfloor$ 个任意单比特错误。

!!! theorem "定理 28.7 (Knill-Laflamme 量子纠错条件)"
    量子码 $\mathcal{C}$（码空间投影 $P$）可以纠正错误集 $\{E_a\}$ 当且仅当存在 Hermite 矩阵 $(\alpha_{ab})$ 使得

    $$
    P E_a^\dagger E_b P = \alpha_{ab} P, \quad \forall a, b.
    $$

    等价地，对码空间的任意正交基 $\{|\psi_i\rangle\}$，

    $$
    \langle\psi_i|E_a^\dagger E_b|\psi_j\rangle = \alpha_{ab}\delta_{ij}, \quad \forall a, b, i, j.
    $$

    直观含义：不同错误要么将码空间映射到正交子空间（可区分），要么对码空间的作用成比例（等价错误）。

??? proof "证明"
    **充分性**：Knill-Laflamme 条件保证 $\{E_a P\}$ 将码空间映射到可区分的子空间。对 Hermite 矩阵 $\alpha = (\alpha_{ab})$ 做酉对角化 $\alpha = W D W^\dagger$，定义"标准错误" $\tilde{E}_c = \sum_a W_{ac} E_a$。则

    $$
    P \tilde{E}_c^\dagger \tilde{E}_{c'} P = d_c \delta_{cc'} P,
    $$

    即标准错误将码空间映射到两两正交的子空间。纠错操作 = 投影到错误子空间（识别哪个错误发生）+ 逆向旋转（恢复原始态）。

    **必要性**：若可以纠错，存在恢复操作 $\mathcal{R}$ 使得对码空间上任意 $\rho$，$\mathcal{R}(\sum_a E_a \rho E_a^\dagger) = \rho$。设 $\mathcal{R}$ 的 Kraus 算子为 $\{R_l\}$。对码空间的正交基向量 $|\psi_i\rangle$ 和 $|\psi_j\rangle$，恢复条件意味着 $\langle\psi_i|E_a^\dagger E_b|\psi_j\rangle$ 不依赖于 $i, j$（当 $i = j$）且当 $i \ne j$ 时为零。$\blacksquare$

!!! definition "定义 28.13 (稳定子码与 CSS 码)"
    **稳定子码**由 $n$ 比特 Pauli 群的一个 Abel 子群 $\mathcal{S}$（稳定子群）定义。码空间为

    $$
    \mathcal{C} = \{|\psi\rangle : S|\psi\rangle = |\psi\rangle,\, \forall S \in \mathcal{S}\}.
    $$

    若 $\mathcal{S}$ 由 $n - k$ 个独立生成元生成，则 $\dim \mathcal{C} = 2^k$。

    **CSS 码**（Calderbank-Shor-Steane）是一类特殊的稳定子码，由两个经典线性码 $C_1 \supset C_2$ 构造，使得 $X$ 错误和 $Z$ 错误可以分别纠正。

!!! example "例 28.7"
    **Shor 9 比特码和 Steane 7 比特码。**

    **Shor $[[9,1,3]]$ 码**将 $1$ 个逻辑比特编码为 $9$ 个物理比特：

    $$
    |0_L\rangle = \frac{1}{2\sqrt{2}}(|000\rangle + |111\rangle)^{\otimes 3}, \quad |1_L\rangle = \frac{1}{2\sqrt{2}}(|000\rangle - |111\rangle)^{\otimes 3}.
    $$

    此码可纠正任意单比特错误。结构：外层 3 比特重复码纠正位翻转（$X$）错误，内层 3 比特相位码纠正相位翻转（$Z$）错误。

    **Steane $[[7,1,3]]$ 码**是 CSS 码，基于经典 $[7,4,3]$ Hamming 码。由于 Hamming 码自对偶（$C_2 = C_1^\perp \subset C_1$），Steane 码仅用 $7$ 个物理比特即可纠正任意单比特错误，比 Shor 码更紧凑。

    验证 Knill-Laflamme 条件：对 Steane 码，任意两个权重 $\le 1$ 的 Pauli 错误 $E_a, E_b$，$E_a^\dagger E_b$ 的权重 $\le 2 < d = 3$，因此 $P E_a^\dagger E_b P = \alpha_{ab} P$。

!!! example "例 28.8"
    **稳定子码的结构。** 5 比特码 $[[5,1,3]]$ 是最小的能纠正任意单比特错误的量子码。其稳定子群由以下 $4 = n - k$ 个生成元生成：

    $$
    g_1 = XZZXI, \quad g_2 = IXZZX, \quad g_3 = XIXZZ, \quad g_4 = ZXIXZ.
    $$

    码空间由所有满足 $g_i|\psi\rangle = |\psi\rangle$（$i = 1, 2, 3, 4$）的态组成，维数 $2^1 = 2$。

    码距 $d = 3$：任何权重 $\le 2$ 的 Pauli 算子要么在稳定子群中，要么与某个生成元反对易，因此可以被检测。5 比特码达到了量子 Singleton 界 $n - k \ge 2(d-1)$，即 $4 \ge 4$。

---

## 28.8 量子信息中的矩阵不等式

<div class="context-flow" markdown>

**核心工具**：von Neumann 熵 $S(\rho) = -\operatorname{tr}(\rho \log \rho)$ → 量子相对熵 → 强次可加性 → 数据处理不等式

**链接**：Ch18 矩阵不等式的量子信息推广

</div>

量子信息论中的核心不等式都可以表述为矩阵函数的不等式。

!!! definition "定义 28.14 (von Neumann 熵)"
    密度矩阵 $\rho$ 的 **von Neumann 熵**定义为

    $$
    S(\rho) = -\operatorname{tr}(\rho \log \rho) = -\sum_i \lambda_i \log \lambda_i,
    $$

    其中 $\{\lambda_i\}$ 为 $\rho$ 的特征值，约定 $0 \log 0 = 0$。

    - $S(\rho) = 0$ 当且仅当 $\rho$ 为纯态。
    - $S(\rho) \le \log d$（$d = \dim \mathcal{H}$），等号当且仅当 $\rho = I/d$（最大混合态）。

!!! definition "定义 28.15 (量子相对熵)"
    两个密度矩阵 $\rho$ 和 $\sigma$ 之间的**量子相对熵**为

    $$
    S(\rho \| \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma),
    $$

    当 $\ker(\sigma) \cap \operatorname{supp}(\rho) = \{0\}$ 时有限，否则定义为 $+\infty$。

!!! theorem "定理 28.8 (强次可加性)"
    对三体系统 $ABC$ 的任意密度矩阵 $\rho_{ABC}$，von Neumann 熵满足**强次可加性**（strong subadditivity, SSA）：

    $$
    S(\rho_{ABC}) + S(\rho_B) \le S(\rho_{AB}) + S(\rho_{BC}),
    $$

    其中 $\rho_{AB} = \operatorname{tr}_C(\rho_{ABC})$ 等为约化密度矩阵。

    等价形式：

    - **条件熵递减**：$S(A|BC) \le S(A|B)$，其中 $S(A|B) = S(\rho_{AB}) - S(\rho_B)$。
    - **互信息非负**：$I(A:C|B) = S(A|B) - S(A|BC) \ge 0$。

??? proof "证明"
    强次可加性的标准证明使用 Lieb 凹性定理。定义映射 $f(X) = \operatorname{tr}(\exp(\log M + \log X))$（$M$ 固定正定矩阵），Lieb 证明了 $f$ 是凹函数。

    另一种证明基于量子相对熵的单调性（数据处理不等式）。量子相对熵在量子信道（CPTP 映射）下单调递减：

    $$
    S(\mathcal{E}(\rho) \| \mathcal{E}(\sigma)) \le S(\rho \| \sigma).
    $$

    取 $\mathcal{E} = \operatorname{tr}_C$（偏迹是 CPTP 映射），$\rho = \rho_{ABC}$，$\sigma = I_A \otimes \rho_{BC} / d_A$，代入并化简即可得到强次可加性。

    具体地，$S(\rho_{ABC} \| I_A/d_A \otimes \rho_{BC}) = \log d_A - S(A|BC)$ 和 $S(\rho_{AB} \| I_A/d_A \otimes \rho_B) = \log d_A - S(A|B)$。偏迹单调性给出 $\log d_A - S(A|BC) \ge \log d_A - S(A|B)$，即 $S(A|B) \ge S(A|BC)$，这就是强次可加性。$\blacksquare$

!!! theorem "定理 28.9 (Klein 不等式与量子相对熵非负性)"
    对任意密度矩阵 $\rho$ 和 $\sigma$，

    $$
    S(\rho \| \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma) \ge 0,
    $$

    等号当且仅当 $\rho = \sigma$。

??? proof "证明"
    这是 Klein 不等式的推论。Klein 不等式断言：对凸函数 $f$ 和自伴算子 $A, B$，

    $$
    \operatorname{tr}(f(A) - f(B) - f'(B)(A - B)) \ge 0.
    $$

    取 $f(t) = t \log t$（凸函数），$A = \rho$，$B = \sigma$，$f'(t) = \log t + 1$。则

    $$
    \operatorname{tr}(\rho \log \rho - \sigma \log \sigma - (\log \sigma + I)(\rho - \sigma)) \ge 0.
    $$

    展开并利用 $\operatorname{tr}(\rho) = \operatorname{tr}(\sigma) = 1$：

    $$
    \operatorname{tr}(\rho \log \rho - \sigma \log \sigma - \rho \log \sigma - \rho + \sigma \log \sigma + \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma) = S(\rho \| \sigma) \ge 0.
    $$

    等号条件由 $f$ 的严格凸性给出：$\rho = \sigma$。$\blacksquare$

!!! example "例 28.9"
    **验证强次可加性：GHZ 态。** 三比特 GHZ 态 $|\text{GHZ}\rangle = \frac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$。

    约化密度矩阵：

    - $\rho_{ABC} = |\text{GHZ}\rangle\langle\text{GHZ}|$（纯态），$S(\rho_{ABC}) = 0$。
    - $\rho_{AB} = \frac{1}{2}(|00\rangle\langle 00| + |11\rangle\langle 11|)$，$S(\rho_{AB}) = \log 2 = 1$。
    - $\rho_{BC} = \frac{1}{2}(|00\rangle\langle 00| + |11\rangle\langle 11|)$，$S(\rho_{BC}) = 1$。
    - $\rho_B = \frac{I}{2}$，$S(\rho_B) = 1$。

    SSA 验证：$S(\rho_{ABC}) + S(\rho_B) = 0 + 1 = 1 \le S(\rho_{AB}) + S(\rho_{BC}) = 1 + 1 = 2$。不等式成立。

    条件互信息 $I(A:C|B) = S(A|B) - S(A|BC) = (S(\rho_{AB}) - S(\rho_B)) - (S(\rho_{ABC}) - S(\rho_{BC})) = (1-1) - (0-1) = 1 \ge 0$。

!!! note "注"
    本章的线性代数工具贯穿了量子信息的各个领域：Hilbert 空间提供状态空间的框架（Ch8），酉矩阵描述量子演化（Ch7），SVD 给出 Schmidt 分解（Ch11），Kronecker 积描述多体系统（Ch19），正定矩阵理论刻画密度矩阵和 POVM（Ch16），矩阵不等式约束信息处理的极限（Ch18）。量子信息科学是线性代数最深刻、最活跃的应用领域之一。

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

本章全景式地展示了线性代数作为量子力学与量子计算唯一“母语”的绝对统治力：


****：以复向量和酉矩阵重构了量子力学的公理体系。薛定谔方程不过是酉群 $U(N)$ 在希尔伯特空间上的连续李群作用。

****：通过张量积（Kronecker 积）将孤立的量子比特拼合成指数级爆炸的组合空间。借助 SVD（Schmidt 分解），我们获得了数学上极其精确的“量子纠缠”判据。

****：超越了传统的正交投影测量，引入了基于正定矩阵的 POVM 测量和刻画混合态的密度矩阵 $\rho$。展示了偏迹运算如何自然地描述子系统信息的丢失。

****：将量子计算机视作一台极其庞大的硬件矩阵乘法器，介绍了 Pauli、Hadamard、CNOT 等基本酉门，及其在构建诸如量子隐形传态和 Grover 搜索等算法中的基石作用。

****：利用有限域 $\mathbb{F}_2$ 上的线性代数和 Pauli 群的代数结构，引出了稳定子码（Stabilizer Codes），证明了不可观测的“综合征测量”如何保护脆弱的量子叠加态。

****：将经典香农熵自然过渡为基于矩阵对数函数的冯·诺依曼熵，并在矩阵不等式的武装下，证明了决定一切量子通信极限的强次可加性（SSA）定理。
