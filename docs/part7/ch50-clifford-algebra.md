# 第 50 章 Clifford 代数与几何代数

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 二次型 (Ch9) · 外代数 (Ch49)

**本章脉络**：Clifford 代数 $\mathrm{Cl}(V,Q)$ 定义 → 几何积 = 内积 + 楔积 → 低维例子（复数、四元数）→ 旋量（rotor）与旋转 → Spin 群 → 矩阵表示 → Bott 周期性 → 几何代数在物理中的应用

**延伸**：Clifford 代数统一了复数、四元数和外代数；Spin 群是正交群的二重覆盖，是量子场论中旋量场的数学基础；几何代数在计算机图形学中提供了比矩阵方法更优雅的旋转和反射表示

</div>

外代数 $\Lambda(V)$ 将向量空间 $V$ 的反对称乘法系统化，但它完全忽略了 $V$ 上的**度量结构**（内积或更一般的二次型）。William Kingdon Clifford 在 1878 年提出的 Clifford 代数巧妙地将度量信息编码到乘法规则中：向量的"平方"不再为零（如外代数中 $v \wedge v = 0$），而是等于其二次型值 $Q(v)$。

这一看似微小的修改带来了深远的后果。Clifford 代数统一了实数、复数、四元数——它们不过是 $\mathbb{R}^0, \mathbb{R}^1, \mathbb{R}^2$ 上特定二次型的 Clifford 代数。更重要的是，Clifford 代数为旋转和反射提供了一种纯代数的、坐标无关的表示，催生了物理学中**旋量**的概念和计算机图形学中的**几何代数**方法。

---

## 50.1 Clifford 代数的定义

<div class="context-flow" markdown>

**核心问题**：如何构建一个代数，使得向量的乘积编码内积信息？

</div>

!!! definition "定义 50.1 (二次型与对称双线性形式)"
    回顾：设 $V$ 是 $\mathbb{F}$-向量空间。**二次型** $Q: V \to \mathbb{F}$ 满足 $Q(\lambda v) = \lambda^2 Q(v)$，且

    $$B(u, v) = \frac{1}{2}[Q(u + v) - Q(u) - Q(v)]$$

    是对称双线性形式（假设 $\operatorname{char}(\mathbb{F}) \neq 2$）。$Q$ 由 $B$ 确定：$Q(v) = B(v, v)$。

    典型例子：$V = \mathbb{R}^n$，$Q(v) = v_1^2 + \cdots + v_p^2 - v_{p+1}^2 - \cdots - v_{p+q}^2$（符号 $(p, q)$）。

!!! definition "定义 50.2 (Clifford 代数)"
    设 $(V, Q)$ 是配备二次型 $Q$ 的 $\mathbb{F}$-向量空间。$(V, Q)$ 的 **Clifford 代数** $\mathrm{Cl}(V, Q)$ 是满足以下万有性质的结合 $\mathbb{F}$-代数：

    存在线性映射 $\iota: V \to \mathrm{Cl}(V, Q)$ 使得

    $$\iota(v)^2 = Q(v) \cdot 1, \quad \forall v \in V,$$

    且对任意结合代数 $A$ 和线性映射 $f: V \to A$ 满足 $f(v)^2 = Q(v) \cdot 1_A$，存在唯一的代数同态 $\tilde{f}: \mathrm{Cl}(V, Q) \to A$ 使得 $\tilde{f} \circ \iota = f$。

!!! theorem "定理 50.1 (Clifford 代数的构造与维数)"
    Clifford 代数可以构造为张量代数的商：

    $$\mathrm{Cl}(V, Q) = T(V) / \langle v \otimes v - Q(v) \cdot 1 : v \in V \rangle.$$

    若 $\dim V = n$，则 $\dim \mathrm{Cl}(V, Q) = 2^n$。

??? proof "证明"
    **存在性：** 令 $\mathcal{I}_Q = \langle v \otimes v - Q(v) \cdot 1 : v \in V \rangle$。定义 $\mathrm{Cl}(V,Q) = T(V)/\mathcal{I}_Q$，$\iota: V \to \mathrm{Cl}(V,Q)$ 为自然映射。万有性质由张量代数的万有性质继承。

    **维数：** 选取 $V$ 的正交基 $\{e_1, \ldots, e_n\}$（即 $B(e_i, e_j) = 0$（$i \neq j$））。在 $\mathrm{Cl}(V,Q)$ 中，$\iota(e_i)$（简记 $e_i$）满足

    $$e_i^2 = Q(e_i) = q_i, \quad e_i e_j = -e_j e_i \quad (i \neq j).$$

    后者由 $\iota(e_i + e_j)^2 = Q(e_i + e_j) = q_i + q_j$ 和 $\iota(e_i + e_j)^2 = e_i^2 + e_i e_j + e_j e_i + e_j^2$ 比较得到。

    因此 $\mathrm{Cl}(V,Q)$ 由单项式 $e_{i_1} e_{i_2} \cdots e_{i_k}$（$i_1 < i_2 < \cdots < i_k$）生成（利用 $e_i e_j = -e_j e_i$ 排序，$e_i^2 = q_i$ 消除重复）。这些单项式共 $2^n$ 个。

    线性无关性可通过构造一个 $2^n$ 维的忠实表示来证明（例如对 $\Lambda(V)$ 的作用，见后文）。

!!! definition "定义 50.3 (Clifford 代数的基本关系)"
    设 $\{e_1, \ldots, e_n\}$ 是 $(V, Q)$ 的正交基，$q_i = Q(e_i)$。$\mathrm{Cl}(V,Q)$ 的乘法完全由以下关系确定：

    $$e_i^2 = q_i, \quad e_i e_j = -e_j e_i \quad (i \neq j).$$

    $\mathrm{Cl}(V,Q)$ 的基为

    $$\{e_I = e_{i_1} e_{i_2} \cdots e_{i_k} : I = \{i_1, \ldots, i_k\} \subseteq \{1, \ldots, n\}, \, i_1 < \cdots < i_k\} \cup \{1\},$$

    共 $2^n$ 个。

!!! example "例 50.1 (退化情形：外代数)"
    当 $Q = 0$（零二次型）时，$e_i^2 = 0$ 对所有 $i$，关系变为 $e_i e_j = -e_j e_i$ 和 $e_i^2 = 0$。这正是外代数 $\Lambda(V)$ 的定义关系。因此

    $$\mathrm{Cl}(V, 0) \cong \Lambda(V).$$

    Clifford 代数是外代数的"带度量版本"。

---

## 50.2 几何积

<div class="context-flow" markdown>

**核心问题**：Clifford 代数中两个向量的乘积如何分解为对称和反对称部分？

</div>

!!! definition "定义 50.4 (几何积)"
    在 $\mathrm{Cl}(V,Q)$ 中，两个向量 $a, b \in V$ 的乘积 $ab$（称为**几何积**，geometric product）可以分解为对称部分和反对称部分：

    $$ab = \frac{1}{2}(ab + ba) + \frac{1}{2}(ab - ba) = a \cdot b + a \wedge b,$$

    其中

    - **内积**（标量部分）：$a \cdot b = \frac{1}{2}(ab + ba) = B(a, b) \in \mathbb{F}$；
    - **楔积**（二向量部分）：$a \wedge b = \frac{1}{2}(ab - ba) \in \Lambda^2(V) \subset \mathrm{Cl}(V,Q)$。

    特别地，$a^2 = a \cdot a = Q(a)$。

!!! theorem "定理 50.2 (几何积的性质)"
    设 $a, b \in V$。则：

    1. $ab = a \cdot b + a \wedge b$；
    2. $a \cdot b = B(a, b)$（关联双线性形式）；
    3. $a \parallel b$（即 $b = \lambda a$）$\Leftrightarrow$ $ab = ba$（$\Leftrightarrow$ $a \wedge b = 0$）；
    4. $a \perp b$（即 $B(a,b) = 0$）$\Leftrightarrow$ $ab = -ba$（$\Leftrightarrow$ $a \cdot b = 0$）；
    5. 若 $a$ 可逆（即 $Q(a) \neq 0$），则 $a^{-1} = \frac{a}{Q(a)}$。

??? proof "证明"
    **(1)-(2)** 由定义直接计算。$(a+b)^2 = Q(a+b) = Q(a) + 2B(a,b) + Q(b)$。展开左边：$(a+b)^2 = a^2 + ab + ba + b^2 = Q(a) + ab + ba + Q(b)$。比较得 $ab + ba = 2B(a,b)$。

    **(3)** 若 $b = \lambda a$，则 $ab = \lambda a^2 = \lambda Q(a) = ba$。反之，$ab = ba$ 意味着 $a \wedge b = 0$，即 $a, b$ 在 $\Lambda^2(V)$ 中的楔积为零，故 $a, b$ 线性相关。

    **(4)** $B(a,b) = 0 \Leftrightarrow ab + ba = 0 \Leftrightarrow ab = -ba$。

    **(5)** 若 $Q(a) \neq 0$，$a \cdot \frac{a}{Q(a)} = \frac{a^2}{Q(a)} = \frac{Q(a)}{Q(a)} = 1$。

!!! example "例 50.2 (欧氏平面中的几何积)"
    取 $V = \mathbb{R}^2$，$Q(v) = v_1^2 + v_2^2$（欧氏度量）。正交基 $\{e_1, e_2\}$ 满足 $e_1^2 = e_2^2 = 1$，$e_1 e_2 = -e_2 e_1$。

    设 $a = a_1 e_1 + a_2 e_2$，$b = b_1 e_1 + b_2 e_2$。则

    $$ab = (a_1 b_1 + a_2 b_2) + (a_1 b_2 - a_2 b_1) e_1 e_2.$$

    标量部分 $a_1 b_1 + a_2 b_2 = a \cdot b$（点积），二向量部分 $(a_1 b_2 - a_2 b_1) e_1 e_2$（其系数是叉积的"$z$-分量"或有向面积）。

    **特别注意**元素 $I = e_1 e_2$（称为**赝标量**）：$I^2 = e_1 e_2 e_1 e_2 = -e_1 e_1 e_2 e_2 = -1$。因此 $I$ 的行为类似于虚数单位！

!!! example "例 50.3 (欧氏三维空间中的几何积)"
    取 $V = \mathbb{R}^3$，$Q(v) = v_1^2 + v_2^2 + v_3^2$。$\mathrm{Cl}(\mathbb{R}^3, Q)$ 的维数为 $2^3 = 8$，基为

    $$\{1, \, e_1, e_2, e_3, \, e_1 e_2, e_1 e_3, e_2 e_3, \, e_1 e_2 e_3\}.$$

    分次结构：标量（1维）、向量（3维）、二向量（3维）、赝标量（1维），共 $1+3+3+1 = 8$。

    赝标量 $I = e_1 e_2 e_3$ 满足 $I^2 = -1$，且 $I$ 与所有元素对易（即 $I$ 在代数的中心）。

    向量叉积可通过几何积表达：

    $$a \times b = -I(a \wedge b) = -\frac{I}{2}(ab - ba).$$

---

## 50.3 低维 Clifford 代数

<div class="context-flow" markdown>

**核心问题**：低维情形的 Clifford 代数同构于哪些熟知的代数？

</div>

!!! theorem "定理 50.3 (低维 Clifford 代数的分类)"
    以下是实向量空间上主要的低维 Clifford 代数（$(p,q)$ 表示符号为 $(p,q)$ 的二次型）：

    | $(p,q)$ | $\mathrm{Cl}(p,q)$ | 维数 | 说明 |
    |:---:|:---:|:---:|:---|
    | $(0,0)$ | $\mathbb{R}$ | $1$ | 实数 |
    | $(1,0)$ | $\mathbb{R} \oplus \mathbb{R}$ | $2$ | 分裂复数 |
    | $(0,1)$ | $\mathbb{C}$ | $2$ | 复数 |
    | $(2,0)$ | $M_2(\mathbb{R})$ | $4$ | $2 \times 2$ 实矩阵 |
    | $(1,1)$ | $M_2(\mathbb{R})$ | $4$ | $2 \times 2$ 实矩阵 |
    | $(0,2)$ | $\mathbb{H}$ | $4$ | 四元数 |
    | $(3,0)$ | $M_2(\mathbb{C})$ | $8$ | $2 \times 2$ 复矩阵 |
    | $(0,3)$ | $\mathbb{H} \oplus \mathbb{H}$ | $8$ | 两份四元数 |
    | $(1,3)$ | $M_2(\mathbb{H})$ | $16$ | $2 \times 2$ 四元数矩阵 |

??? proof "验证主要情形"
    **$\mathrm{Cl}(0,1) \cong \mathbb{C}$：** $V = \mathbb{R}$，$Q(v) = -v^2$。$\mathrm{Cl}(0,1)$ 由 $1$ 和 $e_1$ 生成，$e_1^2 = -1$。这正是 $\mathbb{C} = \{a + bi : i^2 = -1\}$。

    **$\mathrm{Cl}(0,2) \cong \mathbb{H}$：** $V = \mathbb{R}^2$，$Q(v) = -(v_1^2 + v_2^2)$。基 $\{e_1, e_2\}$ 满足 $e_1^2 = e_2^2 = -1$，$e_1 e_2 = -e_2 e_1$。

    令 $i = e_1$，$j = e_2$，$k = e_1 e_2$。则 $k^2 = e_1 e_2 e_1 e_2 = -e_1 e_1 e_2 e_2 = -(-1)(-1) = -1$。且 $ij = k$，$ji = -k$，$jk = i$，$ki = j$，等等。这正是四元数 $\mathbb{H}$ 的定义关系。

    **$\mathrm{Cl}(2,0) \cong M_2(\mathbb{R})$：** $e_1^2 = e_2^2 = 1$，$e_1 e_2 = -e_2 e_1$。对应

    $$e_1 \mapsto \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}, \quad e_2 \mapsto \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}.$$

    验证：$e_1^2 = I$，$e_2^2 = I$，$e_1 e_2 = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix} = -e_2 e_1$。四个基元素 $I, e_1, e_2, e_1 e_2$ 线性无关，恰好张成 $M_2(\mathbb{R})$。

    **$\mathrm{Cl}(1,0) \cong \mathbb{R} \oplus \mathbb{R}$：** $e_1^2 = 1$。令 $e_+ = \frac{1+e_1}{2}$，$e_- = \frac{1-e_1}{2}$，则 $e_+^2 = e_+$，$e_-^2 = e_-$，$e_+ e_- = 0$，$e_+ + e_- = 1$。故 $\mathrm{Cl}(1,0) = \mathbb{R} e_+ \oplus \mathbb{R} e_- \cong \mathbb{R} \oplus \mathbb{R}$。

!!! example "例 50.4 (复数作为 Clifford 代数)"
    $\mathbb{C} \cong \mathrm{Cl}(0,1)$。复数乘法 $(a+bi)(c+di) = (ac-bd) + (ad+bc)i$ 正是 Clifford 积。

    复数的共轭 $\overline{a+bi} = a-bi$ 对应 Clifford 代数的**主反自同构**（reversion）。$|z|^2 = z\bar{z} = a^2+b^2 = Q(-v) \cdot 1$（这里 $Q$ 是负定的）。

!!! example "例 50.5 (四元数作为 Clifford 代数)"
    $\mathbb{H} \cong \mathrm{Cl}(0,2)$。四元数乘法规则 $i^2 = j^2 = k^2 = ijk = -1$ 编码了 $\mathbb{R}^2$ 上的负定二次型。

    四元数的共轭 $\overline{a+bi+cj+dk} = a-bi-cj-dk$ 对应 Clifford 代数的反自同构。$|q|^2 = q\bar{q} = a^2+b^2+c^2+d^2$。

---

## 50.4 Clifford 代数的结构

<div class="context-flow" markdown>

**核心问题**：Clifford 代数有什么内在的代数结构？分次、滤过、周期性？

</div>

!!! definition "定义 50.5 (偶子代数)"
    $\mathrm{Cl}(V,Q)$ 有自然的 $\mathbb{Z}/2\mathbb{Z}$-分次：

    $$\mathrm{Cl}(V,Q) = \mathrm{Cl}^0(V,Q) \oplus \mathrm{Cl}^1(V,Q),$$

    其中 $\mathrm{Cl}^0$ 由偶数个向量的乘积张成（标量、二向量、四向量、...），$\mathrm{Cl}^1$ 由奇数个向量的乘积张成。$\mathrm{Cl}^0(V,Q)$ 是一个子代数，称为**偶子代数**（even subalgebra），维数为 $2^{n-1}$。

!!! theorem "定理 50.4 (偶子代数的同构)"
    $$\mathrm{Cl}^0(p,q) \cong \mathrm{Cl}(p, q-1) \quad \text{（若 } q \geq 1 \text{）},$$
    $$\mathrm{Cl}^0(p,q) \cong \mathrm{Cl}(q, p-1) \quad \text{（若 } p \geq 1 \text{）}.$$

    例如：$\mathrm{Cl}^0(3,0) \cong \mathrm{Cl}(0,2) \cong \mathbb{H}$。这意味着三维欧氏空间的 Clifford 代数的偶子代数同构于四元数。

!!! definition "定义 50.6 (Clifford 代数的反自同构)"
    $\mathrm{Cl}(V,Q)$ 上有三个重要的对合/反自同构：

    1. **分次对合** $\hat{\phantom{x}}$（grade involution）：$\hat{v} = -v$（$v \in V$），$\widehat{ab} = \hat{a}\hat{b}$。作用于 $k$-向量：$\hat{a}_k = (-1)^k a_k$。

    2. **反转** $\tilde{\phantom{x}}$（reversion）：$\tilde{v} = v$（$v \in V$），$\widetilde{ab} = \tilde{b}\tilde{a}$。作用于 $k$-向量：$\tilde{a}_k = (-1)^{k(k-1)/2} a_k$。

    3. **Clifford 共轭** $\bar{\phantom{x}}$：$\bar{a} = \widehat{\tilde{a}}$。作用于 $k$-向量：$\bar{a}_k = (-1)^{k(k+1)/2} a_k$。

    | $k$ | $\hat{a}_k$ | $\tilde{a}_k$ | $\bar{a}_k$ |
    |:---:|:---:|:---:|:---:|
    | 0 | $+a_0$ | $+a_0$ | $+a_0$ |
    | 1 | $-a_1$ | $+a_1$ | $-a_1$ |
    | 2 | $+a_2$ | $-a_2$ | $-a_2$ |
    | 3 | $-a_3$ | $-a_3$ | $+a_3$ |

!!! theorem "定理 50.5 (Bott 周期性)"
    实 Clifford 代数满足 **8-周期性**（Bott periodicity）：

    $$\mathrm{Cl}(n+8, 0) \cong \mathrm{Cl}(n, 0) \otimes M_{16}(\mathbb{R}),$$
    $$\mathrm{Cl}(0, n+8) \cong \mathrm{Cl}(0, n) \otimes M_{16}(\mathbb{R}).$$

    更一般地，$\mathrm{Cl}(p+1, q+1) \cong \mathrm{Cl}(p, q) \otimes M_2(\mathbb{R})$。

    完整的 8-周期分类表：

    | $n \bmod 8$ | $\mathrm{Cl}(n, 0)$ |
    |:---:|:---|
    | 0 | $M_{2^{n/2}}(\mathbb{R})$ |
    | 1 | $M_{2^{(n-1)/2}}(\mathbb{C})$ |
    | 2 | $M_{2^{(n-2)/2}}(\mathbb{H})$ |
    | 3 | $M_{2^{(n-3)/2}}(\mathbb{H}) \oplus M_{2^{(n-3)/2}}(\mathbb{H})$ |
    | 4 | $M_{2^{(n-2)/2}}(\mathbb{H})$ |
    | 5 | $M_{2^{(n-1)/2}}(\mathbb{C})$ |
    | 6 | $M_{2^{n/2}}(\mathbb{R})$ |
    | 7 | $M_{2^{(n-1)/2}}(\mathbb{R}) \oplus M_{2^{(n-1)/2}}(\mathbb{R})$ |

??? proof "证明"
    我们分步骤完整证明。

    **第一步：$\mathrm{Cl}(p+1,q+1) \cong \mathrm{Cl}(p,q) \otimes M_2(\mathbb{R})$。**

    设 $V = \mathbb{R}^{p+q+2}$ 带二次型 $Q$ 且符号为 $(p+1, q+1)$。取正交基 $\{e_1, \ldots, e_p, f_1, \ldots, f_q, e, f\}$，其中 $e_i^2 = 1$，$f_j^2 = -1$，$e^2 = 1$，$f^2 = -1$。

    令 $\varepsilon = ef$。则 $\varepsilon^2 = efef = -effe = -e^2 f^2(-1) \ldots$ 更仔细地：$\varepsilon^2 = (ef)(ef) = e(fe)f = -e(ef)f = -e^2 f^2 = -(1)(-1) = 1$。此外 $\varepsilon$ 与每个 $e_i$ 和 $f_j$ 反交换：$e_i \varepsilon = e_i ef = -ee_i f = -e(e_i f) = -ef e_i \cdot (-1)^0 \ldots$ 准确地说，$e_i e = -ee_i$（正交向量反交换），$e_i f = -f e_i$，故 $e_i \varepsilon = e_i ef = -ee_i f = e(-e_i)f = -e(-f e_i) = ef e_i = \varepsilon e_i$... 需更小心。由于 $e_i$ 与 $e$ 和 $f$ 都正交，$e_i e = -ee_i$，$e_i f = -fe_i$，故

    $$e_i \varepsilon = e_i(ef) = -(ee_i)f = -e(e_i f) = -e(-fe_i) = efe_i = \varepsilon e_i.$$

    等等，这说明 $e_i$ 与 $\varepsilon$ **对易**。类似地 $f_j$ 与 $\varepsilon$ 对易。

    现在定义映射 $\Phi: \mathrm{Cl}(p,q) \otimes M_2(\mathbb{R}) \to \mathrm{Cl}(p+1,q+1)$ 如下。令

    $$e_i \mapsto e_i \otimes I_2, \quad f_j \mapsto f_j \otimes I_2$$

    对应 $\mathrm{Cl}(p,q)$ 的生成元，以及

    $$E_{11} - E_{22} \mapsto \varepsilon, \quad E_{12} + E_{21} \mapsto e, \quad E_{12} - E_{21} \mapsto f$$

    对应 $M_2(\mathbb{R})$ 的生成元。通过验证这些映射保持 Clifford 关系和矩阵关系，可以建立代数同构。维数验证：两边都是 $2^{p+q+2}$ 维。

    **第二步：辅助同构。**

    $$\mathrm{Cl}(n+2, 0) \cong \mathrm{Cl}(0, n) \otimes \mathrm{Cl}(2, 0) \cong \mathrm{Cl}(0, n) \otimes M_2(\mathbb{R}).$$

    设 $\{e_1, \ldots, e_{n+2}\}$ 为 $\mathrm{Cl}(n+2,0)$ 的正交基（$e_i^2 = 1$）。令 $\eta = e_1 e_2$，则 $\eta^2 = -1$。定义 $f_i = e_{i+2}\eta$（$i = 1, \ldots, n$），则 $f_i^2 = e_{i+2}\eta e_{i+2}\eta = -e_{i+2}^2 \eta^2 = -(1)(-1) = 1$... 但我们需要 $f_i^2 = -1$。令 $f_i = \eta e_{i+2}$，则 $f_i^2 = \eta e_{i+2}\eta e_{i+2} = -\eta^2 e_{i+2}^2 = -(-1)(1) = 1$。

    采用不同的构造：$\mathrm{Cl}(0, n+2) \cong \mathrm{Cl}(n, 0) \otimes \mathrm{Cl}(0, 2) \cong \mathrm{Cl}(n, 0) \otimes \mathbb{H}$。

    **第三步：8-周期性。**

    由低维分类（定理 50.3），$\mathrm{Cl}(0,0) = \mathbb{R}$，$\mathrm{Cl}(1,0) = \mathbb{R} \oplus \mathbb{R}$，...，$\mathrm{Cl}(8,0)$ 可通过反复应用第一步计算：

    $$\mathrm{Cl}(2,0) \cong M_2(\mathbb{R}), \quad \mathrm{Cl}(4,0) \cong M_2(\mathbb{H}), \quad \mathrm{Cl}(8,0) \cong M_{16}(\mathbb{R}).$$

    最后一个可验证：$\mathrm{Cl}(8,0) \cong \mathrm{Cl}(6,0) \otimes M_2(\mathbb{R}) \cong \cdots \cong M_{16}(\mathbb{R})$。

    因此 $\mathrm{Cl}(n+8,0) \cong \mathrm{Cl}(n,0) \otimes \mathrm{Cl}(8,0) \cong \mathrm{Cl}(n,0) \otimes M_{16}(\mathbb{R})$。

    这里用到了 Clifford 代数的张量积分解：若 $V = V_1 \perp V_2$（正交直和分解），则在分次意义下 $\mathrm{Cl}(V, Q) \cong \mathrm{Cl}(V_1, Q|_{V_1}) \hat{\otimes} \mathrm{Cl}(V_2, Q|_{V_2})$（分次张量积），这对正定二次型简化为普通张量积。

    对 $\mathrm{Cl}(0, n)$ 的 8-周期性完全类似：$\mathrm{Cl}(0, 8) \cong M_{16}(\mathbb{R})$，故 $\mathrm{Cl}(0, n+8) \cong \mathrm{Cl}(0, n) \otimes M_{16}(\mathbb{R})$。

---

## 50.5 旋量与旋转

<div class="context-flow" markdown>

**核心问题**：如何用 Clifford 代数元素表示旋转和反射？Spin 群如何自然出现？

</div>

!!! theorem "定理 50.6 (Clifford 代数中的反射)"
    设 $v \in V$ 满足 $Q(v) \neq 0$。映射

    $$\rho_v: V \to V, \quad \rho_v(w) = -v w v^{-1} = -\frac{v w v}{Q(v)}$$

    是关于 $v^{\perp}$（$v$ 的正交补超平面）的**反射**。即

    $$\rho_v(w) = w - 2\frac{B(w,v)}{Q(v)} v.$$

??? proof "证明"
    $vwv^{-1} = vw \frac{v}{Q(v)}$。利用 $vw = 2B(v,w) - wv$：

    $$vwv = v(2B(v,w) - wv) \cdot \frac{v}{Q(v)} \cdot Q(v) = 2B(v,w)v - wv^2 = 2B(v,w)v - Q(v)w.$$

    （更直接地：利用 $vw + wv = 2B(v,w)$。）

    故 $-vwv^{-1} = -\frac{vwv}{Q(v)} = w - 2\frac{B(v,w)}{Q(v)}v$，正好是关于 $v^{\perp}$ 的反射公式。

!!! definition "定义 50.7 (旋量 / Rotor)"
    **旋量**（rotor）是 $\mathrm{Cl}^0(V,Q)$ 中的可逆元素 $R$，满足

    $$R \tilde{R} = 1,$$

    其中 $\tilde{R}$ 是 $R$ 的反转。旋量通过**夹心积**（sandwich product）作用于向量：

    $$v \mapsto R v \tilde{R}.$$

    每个旋量定义 $V$ 上的一个旋转（保持 $Q$ 的正交变换，且行列式为 $+1$）。

!!! theorem "定理 50.7 (Cartan-Dieudonné 定理的 Clifford 版本)"
    每个正交变换 $\varphi \in O(V, Q)$ 可以写成有限个反射的复合：

    $$\varphi = \rho_{v_1} \circ \rho_{v_2} \circ \cdots \circ \rho_{v_m}, \quad m \leq n.$$

    对应地，在 Clifford 代数中：

    $$\varphi(w) = (\pm 1)^m (v_1 v_2 \cdots v_m) w (v_1 v_2 \cdots v_m)^{-1}.$$

    若 $m$ 为偶数，$\varphi \in SO(V,Q)$（旋转），$R = v_1 v_2 \cdots v_m \in \mathrm{Cl}^0$ 是旋量。

!!! definition "定义 50.8 (Pin 群和 Spin 群)"
    **Pin 群**：

    $$\mathrm{Pin}(V,Q) = \{v_1 v_2 \cdots v_m \in \mathrm{Cl}(V,Q) : v_i \in V, \, Q(v_i) = \pm 1\}.$$

    **Spin 群**：

    $$\mathrm{Spin}(V,Q) = \mathrm{Pin}(V,Q) \cap \mathrm{Cl}^0(V,Q)$$

    $$= \{v_1 v_2 \cdots v_{2l} : v_i \in V, \, Q(v_i) = \pm 1\}.$$

!!! theorem "定理 50.8 (Spin 群是 SO 的二重覆盖)"
    映射 $\rho: \mathrm{Spin}(V,Q) \to SO(V,Q)$，$\rho(R)(w) = R w \tilde{R}$，是**满的群同态**，其核为 $\{+1, -1\}$。因此

    $$SO(V,Q) \cong \mathrm{Spin}(V,Q) / \{\pm 1\}.$$

    $\mathrm{Spin}(V,Q)$ 是 $SO(V,Q)$ 的**二重覆盖**。

!!! example "例 50.6 (三维旋转)"
    在 $\mathrm{Cl}(3,0)$ 中，绕单位向量 $\hat{n}$ 旋转角度 $\theta$ 的旋量为

    $$R = \exp\left(-\frac{\theta}{2} \hat{n} I_3\right) = \cos\frac{\theta}{2} - I_3 \hat{n} \sin\frac{\theta}{2},$$

    其中 $I_3 = e_1 e_2 e_3$ 是赝标量。注意 $I_3 \hat{n}$ 是一个二向量。

    更具体地，绕 $e_3$ 轴旋转角度 $\theta$：

    $$R = \cos\frac{\theta}{2} - e_1 e_2 \sin\frac{\theta}{2}.$$

    验证：$R e_1 \tilde{R} = e_1 \cos\theta + e_2 \sin\theta$，$R e_2 \tilde{R} = -e_1 \sin\theta + e_2 \cos\theta$，$R e_3 \tilde{R} = e_3$。

    这与四元数旋转 $q = \cos\frac{\theta}{2} + \sin\frac{\theta}{2}(n_1 i + n_2 j + n_3 k)$，$v \mapsto qvq^{-1}$ 完全一致。事实上 $\mathrm{Spin}(3) \cong SU(2) \cong S^3$（三维球面/单位四元数群）。

!!! theorem "定理 50.13 ($\mathrm{Spin}(3) \cong SU(2)$)"
    存在群同构 $\mathrm{Spin}(3) \cong SU(2)$。

??? proof "证明"
    **第一步：构造同态。** $\mathrm{Cl}(3,0) \cong M_2(\mathbb{C})$（由 Pauli 矩阵给出）。在此同构下，$\mathrm{Cl}^0(3,0)$ 的偶子代数由 $\{I_2, \sigma_1\sigma_2, \sigma_1\sigma_3, \sigma_2\sigma_3\} = \{I_2, i\sigma_3, -i\sigma_2, i\sigma_1\}$ 张成，同构于 $\mathbb{H}$，也同构于 $M_2(\mathbb{C})$ 中形如 $\begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix}$ 的矩阵全体。

    $\mathrm{Spin}(3)$ 由 $\mathrm{Cl}^0(3,0)$ 中满足 $R\tilde{R} = 1$ 的元素构成。在 Pauli 矩阵表示下，反转 $\tilde{R}$ 对应 $R^\dagger$（Hermite 共轭），故条件 $R\tilde{R} = 1$ 变为 $RR^\dagger = I_2$，即 $R$ 是酉矩阵。

    又 $R \in \mathrm{Cl}^0$，故 $R = \begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix}$，$|a|^2 + |b|^2 = 1$。这正是 $SU(2)$ 的定义。

    **第二步：验证群同态。** $\mathrm{Spin}(3)$ 中的乘法对应矩阵乘法，故嵌入 $\mathrm{Spin}(3) \hookrightarrow SU(2)$ 是群同态。

    **第三步：满射性。** 每个 $SU(2)$ 元素 $U = \begin{pmatrix} a & -\bar{b} \\ b & \bar{a} \end{pmatrix}$（$|a|^2+|b|^2=1$）可以写成 $U = \cos\frac{\theta}{2} I_2 - i\sin\frac{\theta}{2}(\hat{n} \cdot \vec{\sigma})$ 的形式（其中 $\hat{n}$ 是单位向量），这对应于 $\mathrm{Spin}(3)$ 中的旋量。故映射是满射。

    维数和连通性都匹配（$\mathrm{Spin}(3)$ 和 $SU(2)$ 都是 $3$ 维紧致连通 Lie 群），加上双射性质，这给出了群同构 $\mathrm{Spin}(3) \cong SU(2)$。

!!! example "例 50.7 (旋转的复合)"
    两个旋转的复合在 Clifford 代数中表现为旋量的乘积：

    $$R_{12} = R_1 R_2.$$

    这比旋转矩阵的乘法更紧凑（$4$ 个分量 vs. $9$ 个矩阵元素），且避免了万向节锁（gimbal lock）问题。

---

## 50.6 矩阵表示

<div class="context-flow" markdown>

**核心问题**：Clifford 代数如何用矩阵实现？物理学中的 Pauli 矩阵和 Dirac 矩阵从何而来？

</div>

!!! definition "定义 50.9 (Pauli 矩阵)"
    **Pauli 矩阵**是 $\mathrm{Cl}(3,0)$ 的 $2 \times 2$ 复矩阵表示中 $e_1, e_2, e_3$ 的像：

    $$\sigma_1 = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_2 = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \sigma_3 = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.$$

    验证：$\sigma_i^2 = I_2$，$\sigma_i \sigma_j = -\sigma_j \sigma_i$（$i \neq j$），$\sigma_1 \sigma_2 \sigma_3 = iI_2$。

    这些矩阵连同 $I_2$ 在 $\mathbb{C}$ 上线性无关，张成 $M_2(\mathbb{C})$。这给出同构 $\mathrm{Cl}(3,0) \cong M_2(\mathbb{C})$。

!!! theorem "定理 50.9 (Pauli 矩阵的性质)"
    Pauli 矩阵满足：

    1. $\sigma_i \sigma_j + \sigma_j \sigma_i = 2\delta_{ij} I_2$（Clifford 关系）；
    2. $\sigma_i \sigma_j = \delta_{ij} I_2 + i\varepsilon_{ijk} \sigma_k$（乘法公式）；
    3. $\operatorname{tr}(\sigma_i) = 0$；
    4. $\det(\sigma_i) = -1$；
    5. $\sigma_i^{\dagger} = \sigma_i$（Hermite 性）。

!!! definition "定义 50.10 (Dirac 矩阵)"
    **Dirac 矩阵**（或 $\gamma$ 矩阵）是 $\mathrm{Cl}(1,3)$（闵可夫斯基时空）的 $4 \times 4$ 复矩阵表示：

    $$\gamma^0 = \begin{pmatrix} I_2 & 0 \\ 0 & -I_2 \end{pmatrix}, \quad \gamma^i = \begin{pmatrix} 0 & \sigma_i \\ -\sigma_i & 0 \end{pmatrix}, \quad i = 1, 2, 3.$$

    满足

    $$\gamma^{\mu} \gamma^{\nu} + \gamma^{\nu} \gamma^{\mu} = 2 \eta^{\mu\nu} I_4,$$

    其中 $\eta = \operatorname{diag}(1, -1, -1, -1)$ 是闵可夫斯基度量。

    此外定义

    $$\gamma^5 = i\gamma^0 \gamma^1 \gamma^2 \gamma^3 = \begin{pmatrix} 0 & I_2 \\ I_2 & 0 \end{pmatrix},$$

    满足 $(\gamma^5)^2 = I_4$，$\gamma^5 \gamma^{\mu} = -\gamma^{\mu} \gamma^5$。$\gamma^5$ 定义了手征性（chirality）。

!!! example "例 50.8 (旋量表示)"
    $\mathrm{Spin}(3) \cong SU(2)$ 的自然表示空间 $\mathbb{C}^2$ 称为**旋量空间**。旋量 $\psi = \begin{pmatrix} \psi_1 \\ \psi_2 \end{pmatrix} \in \mathbb{C}^2$ 在旋转 $R = \exp(-\frac{\theta}{2} \hat{n} \cdot \vec{\sigma})$ 下变换为

    $$\psi \mapsto R \psi = \exp\left(-\frac{i\theta}{2} \hat{n} \cdot \vec{\sigma}\right) \psi.$$

    注意旋转 $2\pi$（$\theta = 2\pi$）给出 $R = -I_2$，即 $\psi \mapsto -\psi$。旋量在 $2\pi$ 旋转下变号！这是 Spin 群二重覆盖的物理表现。

    $\mathrm{Spin}(1,3) \cong SL_2(\mathbb{C})$ 的自然表示空间 $\mathbb{C}^4$ 是 Dirac 旋量空间，$\gamma^5$ 的特征子空间给出 Weyl 旋量。

---

## 50.7 几何代数的应用

<div class="context-flow" markdown>

**核心问题**：几何代数如何简化旋转、反射和其他几何运算的表达？它在物理学中有哪些应用？

</div>

!!! theorem "定理 50.10 (旋转和反射的统一表示)"
    在 $\mathrm{Cl}(n,0)$（欧氏 Clifford 代数）中：

    - **反射**关于超平面 $v^{\perp}$：$w \mapsto -vwv^{-1}$（$v$ 为超平面法向量）；
    - **旋转**绕平面 $B = e_i \wedge e_j$：$w \mapsto RwR^{-1}$，$R = \exp(-\frac{\theta}{2}B)$；
    - **两个反射的复合**是旋转：$\rho_{v_1} \circ \rho_{v_2}$ 对应 $R = v_1 v_2$。

    这一框架对任意维数都适用，不需要选择坐标系或欧拉角。

!!! example "例 50.9 (用几何代数做投影)"
    向量 $w$ 在向量 $v$（$Q(v) \neq 0$）上的正交投影：

    $$w_{\parallel} = (w \cdot v) v^{-1} = \frac{B(w,v)}{Q(v)} v,$$

    $$w_{\perp} = (w \wedge v) v^{-1} = w - w_{\parallel}.$$

    **注意**：这里利用了 $wv = w \cdot v + w \wedge v$，故 $w \cdot v = \frac{1}{2}(wv + vw)$，$w \wedge v = \frac{1}{2}(wv - vw)$。

!!! example "例 50.10 (Maxwell 方程的几何代数形式)"
    在时空代数 $\mathrm{Cl}(1,3)$ 中，电磁场张量 $F$ 可以写成二向量：

    $$F = \vec{E} + Ic\vec{B},$$

    其中 $\vec{E} = E_1 \gamma_1 \gamma_0 + E_2 \gamma_2 \gamma_0 + E_3 \gamma_3 \gamma_0$ 是电场二向量，$I = \gamma_0 \gamma_1 \gamma_2 \gamma_3$ 是赝标量，$\vec{B}$ 类似。

    全部四个 Maxwell 方程统一为一个方程：

    $$\nabla F = J,$$

    其中 $\nabla = \gamma^{\mu} \partial_{\mu}$ 是时空向量导数，$J$ 是电流密度四维向量。

    标量部分给出 $\nabla \cdot \vec{E} = \rho$（Gauss 定律），向量部分给出 $\nabla \times \vec{B} - \partial_t \vec{E} = \vec{J}$（Ampere 定律），等等。

!!! example "例 50.11 (共形几何代数简介)"
    **共形几何代数**（CGA）将 $\mathbb{R}^n$ 嵌入到 $\mathrm{Cl}(n+1,1)$ 中，使得平移、旋转、缩放和特殊共形变换都可以用旋量表示。

    基本思想：$\mathbb{R}^n$ 中的点 $x$ 表示为空向量

    $$X = x + \frac{1}{2}|x|^2 e_{\infty} + e_0,$$

    其中 $e_0, e_{\infty}$ 是额外的两个基向量，$e_0 \cdot e_{\infty} = -1$。

    在 CGA 中：

    - **球面**（含平面）表示为向量；
    - **圆**（含直线）表示为二向量（两球面的交）；
    - **旋转**：$X \mapsto RXR^{-1}$，$R = \exp(-\frac{\theta}{2}B)$；
    - **平移**：$X \mapsto TXT^{-1}$，$T = 1 - \frac{1}{2}t e_{\infty}$。

    这在计算机图形学和机器人学中极为有用。

!!! theorem "定理 50.11 (Clifford 代数的分类定理)"
    实 Clifford 代数 $\mathrm{Cl}(p,q)$ 完全由 $(p-q) \bmod 8$ 决定（Bott 周期性），同构于以下形式之一：

    $$M_N(\mathbb{R}), \quad M_N(\mathbb{C}), \quad M_N(\mathbb{H}), \quad M_N(\mathbb{R}) \oplus M_N(\mathbb{R}), \quad M_N(\mathbb{H}) \oplus M_N(\mathbb{H}).$$

    具体地，$N = 2^{\lfloor (p+q)/2 \rfloor}$ 或其适当因子，类型由 $(p-q) \bmod 8$ 决定：

    | $(p-q) \bmod 8$ | $\mathrm{Cl}(p,q)$ |
    |:---:|:---|
    | $0$ | $M_N(\mathbb{R})$ |
    | $1$ | $M_N(\mathbb{C})$ |
    | $2$ | $M_N(\mathbb{H})$ |
    | $3$ | $M_N(\mathbb{H}) \oplus M_N(\mathbb{H})$ |
    | $4$ | $M_N(\mathbb{H})$ |
    | $5$ | $M_N(\mathbb{C})$ |
    | $6$ | $M_N(\mathbb{R})$ |
    | $7$ | $M_N(\mathbb{R}) \oplus M_N(\mathbb{R})$ |

    复 Clifford 代数更简单（2-周期性）：

    $$\mathrm{Cl}_{\mathbb{C}}(2m) \cong M_{2^m}(\mathbb{C}), \quad \mathrm{Cl}_{\mathbb{C}}(2m+1) \cong M_{2^m}(\mathbb{C}) \oplus M_{2^m}(\mathbb{C}).$$

??? proof "证明"
    **归纳法证明**，利用以下递推关系：

    **(a)** $\mathrm{Cl}(p+1, q+1) \cong \mathrm{Cl}(p, q) \otimes M_2(\mathbb{R})$（定理 50.5 已证）。

    **(b)** $\mathrm{Cl}(p+1, q) \cong \mathrm{Cl}(q, p) \otimes \mathrm{Cl}(1, 0)$（符号翻转张量积）。

    证明 (b)：设 $\{e_1, \ldots, e_p, f_1, \ldots, f_q, g\}$ 为 $\mathrm{Cl}(p+1, q)$ 的正交基（$e_i^2 = 1$，$f_j^2 = -1$，$g^2 = 1$）。令 $e_i' = e_i g$，$f_j' = f_j g$。则 $(e_i')^2 = e_i g e_i g = -e_i^2 g^2 = -1$，$(f_j')^2 = f_j g f_j g = -f_j^2 g^2 = 1$。且 $e_i', f_j'$ 之间反交换，它们与 $g$ 对易（经计算 $e_i' g = e_i g^2 = e_i$，$g e_i' = g e_i g = -e_i g^2 = -e_i$... 故实际上 $e_i'$ 与 $g$ 反交换）。$\{e_i', f_j'\}$ 生成与 $\mathrm{Cl}(q, p)$ 同构的子代数，$g$ 生成 $\mathrm{Cl}(1, 0) \cong \mathbb{R} \oplus \mathbb{R}$。

    **(c)** 类似地，$\mathrm{Cl}(p, q+1) \cong \mathrm{Cl}(q, p) \otimes \mathrm{Cl}(0, 1)$。

    **基础**：低维情形 $\mathrm{Cl}(0,0) = \mathbb{R}$，$\mathrm{Cl}(1,0) = \mathbb{R} \oplus \mathbb{R}$，$\mathrm{Cl}(0,1) = \mathbb{C}$（定理 50.3 已验证）。

    **归纳步骤**：对 $p + q$ 归纳。利用关系 (a)，每当 $p \geq 1$ 且 $q \geq 1$，可将 $(p+q)$ 维的分类归结为 $(p+q-2)$ 维。利用 (b) 和 (c)，可处理仅增加 $p$ 或 $q$ 的情形。

    关键的矩阵代数张量积公式：$M_a(\mathbb{K}) \otimes M_b(\mathbb{K}') \cong M_{ab}(\mathbb{K}'')$，其中 $\mathbb{K}''$ 由 Brauer 群中 $[\mathbb{K}] \cdot [\mathbb{K}']$ 决定。在实情形：

    - $\mathbb{R} \otimes \mathbb{R} \cong \mathbb{R}$，$\mathbb{R} \otimes \mathbb{C} \cong \mathbb{C}$，$\mathbb{R} \otimes \mathbb{H} \cong \mathbb{H}$；
    - $\mathbb{C} \otimes \mathbb{C} \cong \mathbb{C} \oplus \mathbb{C}$（作为实代数）；
    - $\mathbb{H} \otimes \mathbb{H} \cong M_4(\mathbb{R})$；
    - $(\mathbb{K} \oplus \mathbb{K}) \otimes \mathbb{K}' \cong \mathbb{K} \otimes \mathbb{K}' \oplus \mathbb{K} \otimes \mathbb{K}'$。

    将这些规则与 $\mathrm{Cl}(1,0) \cong \mathbb{R} \oplus \mathbb{R}$，$\mathrm{Cl}(0,1) \cong \mathbb{C}$ 结合，可从低维逐步推出所有 $(p,q)$ 的分类。8-周期性来自 Brauer 群 $\mathrm{Br}(\mathbb{R}) = \{[\mathbb{R}], [\mathbb{H}]\} \cong \mathbb{Z}/2\mathbb{Z}$ 和上述张量积规则。

!!! theorem "定理 50.12 (Chevalley 同构)"
    设 $(V, Q)$ 是配备二次型的 $n$ 维 $\mathbb{F}$-向量空间。存在**典范的向量空间同构**

    $$\sigma: \Lambda(V) \xrightarrow{\sim} \mathrm{Cl}(V, Q),$$

    即 Clifford 代数作为**向量空间**（非代数）同构于外代数，维数都是 $2^n$。

    更精确地，$\mathrm{Cl}(V, Q)$ 具有自然的**滤过**（filtration）

    $$F_0 \subset F_1 \subset F_2 \subset \cdots \subset F_n = \mathrm{Cl}(V, Q),$$

    其中 $F_k$ 由至多 $k$ 个向量乘积的线性组合张成。关联分次（associated graded）为

    $$\mathrm{gr}(\mathrm{Cl}(V, Q)) = \bigoplus_{k=0}^n F_k / F_{k-1} \cong \Lambda(V).$$

    Chevalley 同构 $\sigma$ 可具体构造为：对 $v_1, \ldots, v_k \in V$，

    $$\sigma(v_1 \wedge \cdots \wedge v_k) = \frac{1}{k!}\sum_{\tau \in S_k} \operatorname{sgn}(\tau)\, v_{\tau(1)} v_{\tau(2)} \cdots v_{\tau(k)},$$

    即全反对称化。这个映射是向量空间同构，但**不是代数同构**（除非 $Q = 0$，此时 $\mathrm{Cl}(V, 0) = \Lambda(V)$）。

??? proof "证明"
    **构造映射**：定义 $\sigma_k: \Lambda^k(V) \to \mathrm{Cl}(V, Q)$ 为上述全反对称化。首先验证 $\sigma_k$ 良定义：它保持交替性，因为交换两个变元后 $\operatorname{sgn}(\tau)$ 变号。

    **单射性**：选取正交基 $\{e_1, \ldots, e_n\}$（$B(e_i, e_j) = 0$（$i \neq j$））。在正交基下，

    $$\sigma(e_{i_1} \wedge \cdots \wedge e_{i_k}) = e_{i_1} e_{i_2} \cdots e_{i_k} \quad (i_1 < \cdots < i_k),$$

    因为在 Clifford 代数中正交向量反交换，反对称化后只剩递增序的乘积。这些元素 $\{e_I\}$ 在 $\mathrm{Cl}(V, Q)$ 中线性无关（由维数计数），故 $\sigma$ 是单射。

    **满射性**：$\mathrm{Cl}(V, Q)$ 由 $\{e_I : I \subseteq \{1, \ldots, n\}\}$ 张成，每个都在 $\sigma$ 的像中。

    **滤过结构**：$F_k = \sigma(\bigoplus_{j=0}^k \Lambda^j(V))$。关系 $e_i e_j = -e_j e_i + 2B(e_i, e_j)$ 表明 Clifford 乘积将两个 $1$-向量映到 $F_2$ 中的元素，但模 $F_0$（标量部分 $2B(e_i, e_j)$）后给出 $e_i \wedge e_j$。这正是关联分次同构于外代数的原因。

---

### 本章总结

Clifford 代数 $\mathrm{Cl}(V,Q)$ 是将向量空间的度量结构（二次型 $Q$）编码到代数乘法中的万有构造。它的核心思想极为简洁：**向量的平方等于其"长度"**（$v^2 = Q(v)$），而正交向量反交换（$vw = -wv$）。从这个简单的出发点，自然地产生了：

| 构造 | 来源 |
|:---|:---|
| 复数 $\mathbb{C}$ | $\mathrm{Cl}(0,1)$ |
| 四元数 $\mathbb{H}$ | $\mathrm{Cl}(0,2)$ |
| 外代数 $\Lambda(V)$ | $\mathrm{Cl}(V, 0)$ |
| Spin 群 | $\mathrm{Cl}^0$ 中的特定子群 |
| 旋量表示 | Clifford 代数的矩阵表示 |

---

### 习题

!!! exercise "习题 50.1"
    在 $\mathrm{Cl}(2,0)$ 中，验证 $e_1 \mapsto \sigma_3$，$e_2 \mapsto \sigma_1$ 给出同构 $\mathrm{Cl}(2,0) \cong M_2(\mathbb{R})$。

!!! exercise "习题 50.2"
    在 $\mathrm{Cl}(3,0)$ 中，计算绕 $\hat{n} = \frac{1}{\sqrt{3}}(e_1 + e_2 + e_3)$ 旋转 $120°$ 的旋量 $R$，并验证 $R e_1 \tilde{R} = e_2$。

!!! exercise "习题 50.3"
    证明 $\mathrm{Cl}(0,2) \cong \mathbb{H}$ 中，偶子代数 $\mathrm{Cl}^0(0,2) \cong \mathbb{C}$。

!!! exercise "习题 50.4"
    验证 Dirac 矩阵满足 $\{\gamma^{\mu}, \gamma^{\nu}\} = 2\eta^{\mu\nu}I_4$。

!!! exercise "习题 50.5"
    设 $a, b$ 是 $\mathrm{Cl}(n,0)$ 中的单位向量（$a^2 = b^2 = 1$）。证明 $R = ab$ 是旋量（$R\tilde{R} = 1$），并说明 $v \mapsto RvR^{-1}$ 的几何意义。

!!! exercise "习题 50.6"
    证明：在 $\mathrm{Cl}(3,0)$ 中，赝标量 $I = e_1 e_2 e_3$ 满足 $Iv = vI$（对所有 $v \in \mathbb{R}^3$），即 $I$ 在代数中心。
