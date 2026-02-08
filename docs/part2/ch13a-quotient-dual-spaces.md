# 第 13A 章 商空间与对偶空间

<div class="context-flow" markdown>

**前置**：向量空间与子空间（第 4 章）· 线性映射（第 5 章）· 维数公式（第 5 章）· **本章脉络**：等价类与商空间 → 维数公式与同构定理 → 线性泛函与对偶空间 → 对偶基 → 零化子 → 转置映射 → 双对偶与自然同构 → 有限维应用
**一句话本质**：商空间将子空间"折叠为零"从而简化结构，对偶空间将线性映射的信息编码进线性泛函——二者是理解线性代数深层对称性的关键工具

</div>

商空间和对偶空间是线性代数中两个深刻而优雅的构造。商空间通过将子空间"坍缩为一点"来简化向量空间的结构，是代数中等价关系和同余思想的具体体现。对偶空间则从"测量"的角度审视向量空间，将线性泛函本身组织成一个新的向量空间。两者的交互——零化子、转置映射、双对偶——揭示了线性代数中深层的对称性。

---

## 13A.1 商空间的构造

<div class="context-flow" markdown>

**核心问题**：给定子空间 $W$，能否将 $V$ 中"模 $W$ 等价"的向量视为一个元素？ → 陪集 → 商空间 $V/W$

</div>

### 陪集与等价关系

!!! definition "定义 13A.1 (陪集)"
    设 $V$ 为 $\mathbb{F}$ 上的向量空间，$W$ 为 $V$ 的子空间。对 $v \in V$，集合

    $$
    v + W = \{v + w : w \in W\}
    $$

    称为 $v$ 关于 $W$ 的**陪集**（coset）或**仿射子集**（affine subset）。元素 $v$ 称为该陪集的一个**代表元**（representative）。

!!! theorem "定理 13A.1 (陪集的等价刻画)"
    设 $u, v \in V$，$W$ 为 $V$ 的子空间。以下条件等价：

    1. $u + W = v + W$；
    2. $u - v \in W$；
    3. $u \in v + W$。

??? proof "证明"
    $(1) \Rightarrow (3)$：$u = u + 0 \in u + W = v + W$。

    $(3) \Rightarrow (2)$：$u \in v + W$ 意味着 $u = v + w$（某个 $w \in W$），故 $u - v = w \in W$。

    $(2) \Rightarrow (1)$：设 $u - v = w_0 \in W$。对任意 $u + w \in u + W$，$u + w = v + w_0 + w \in v + W$（因 $w_0 + w \in W$）。同理 $v + W \subseteq u + W$。$\blacksquare$

!!! note "注"
    陪集的本质是等价关系：定义 $u \sim v$ 当且仅当 $u - v \in W$，则 $\sim$ 是 $V$ 上的等价关系，$v + W$ 恰好是 $v$ 所在的等价类。

!!! definition "定义 13A.2 (商空间)"
    设 $W$ 为 $V$ 的子空间。所有陪集的集合

    $$
    V/W = \{v + W : v \in V\}
    $$

    配上以下运算：

    - **加法**：$(u + W) + (v + W) = (u + v) + W$；
    - **标量乘法**：$\lambda(v + W) = \lambda v + W$，$\lambda \in \mathbb{F}$，

    构成 $\mathbb{F}$ 上的向量空间，称为 $V$ 关于 $W$ 的**商空间**（quotient space）。$V/W$ 的零元素为 $0 + W = W$。

!!! theorem "定理 13A.2 (商空间运算的良定义性)"
    商空间 $V/W$ 的加法和标量乘法与代表元的选取无关。

??? proof "证明"
    设 $u + W = u' + W$，$v + W = v' + W$。则 $u - u' \in W$，$v - v' \in W$。

    **加法**：$(u + v) - (u' + v') = (u - u') + (v - v') \in W$，故 $(u + v) + W = (u' + v') + W$。

    **标量乘法**：$\lambda u - \lambda u' = \lambda(u - u') \in W$，故 $\lambda u + W = \lambda u' + W$。$\blacksquare$

!!! example "例 13A.1"
    设 $V = \mathbb{R}^3$，$W = \{(x, y, 0) : x, y \in \mathbb{R}\}$（$xy$-平面）。

    陪集 $(a, b, c) + W = \{(x, y, c) : x, y \in \mathbb{R}\}$，即高度为 $c$ 的水平平面。两个陪集相等当且仅当它们的第三坐标相同。因此 $V/W \cong \mathbb{R}$，其中同构映射为 $(a, b, c) + W \mapsto c$。

!!! example "例 13A.2"
    设 $V = \mathbb{R}^2$，$W = \{(t, t) : t \in \mathbb{R}\}$（过原点斜率为 $1$ 的直线）。

    陪集 $(a, b) + W = \{(a + t, b + t) : t \in \mathbb{R}\}$，这是斜率为 $1$ 的直线族。两个陪集相同当且仅当 $b - a$ 相等。故 $V/W \cong \mathbb{R}$。

### 商映射

!!! definition "定义 13A.3 (商映射)"
    **商映射**（quotient map）$\pi: V \to V/W$ 定义为

    $$
    \pi(v) = v + W.
    $$

    $\pi$ 是满射线性映射，$\ker \pi = W$。

!!! theorem "定理 13A.3 (商映射的性质)"
    商映射 $\pi: V \to V/W$ 满足：

    1. $\pi$ 是线性映射；
    2. $\pi$ 是满射；
    3. $\ker \pi = W$。

??? proof "证明"
    1. $\pi(u + v) = (u + v) + W = (u + W) + (v + W) = \pi(u) + \pi(v)$；$\pi(\lambda v) = \lambda v + W = \lambda(v + W) = \lambda \pi(v)$。

    2. 对任意 $v + W \in V/W$，$\pi(v) = v + W$。

    3. $\pi(v) = 0_{V/W} = W$ 当且仅当 $v + W = W$ 当且仅当 $v \in W$。$\blacksquare$

---

## 13A.2 商空间的维数与同构定理

<div class="context-flow" markdown>

**核心问题**：$V/W$ 的维数是多少？ → 维数公式 → 第一同构定理——线性代数中最重要的结构定理之一

</div>

!!! theorem "定理 13A.4 (商空间的维数)"
    设 $V$ 为有限维向量空间，$W$ 为 $V$ 的子空间。则

    $$
    \dim(V/W) = \dim V - \dim W.
    $$

??? proof "证明"
    商映射 $\pi: V \to V/W$ 是满射且 $\ker \pi = W$。由维数公式（秩-零化度定理）：

    $$
    \dim V = \dim \ker \pi + \dim \operatorname{im} \pi = \dim W + \dim(V/W).
    $$

    因此 $\dim(V/W) = \dim V - \dim W$。$\blacksquare$

!!! note "注"
    $\dim(V/W)$ 也称为 $W$ 在 $V$ 中的**余维数**（codimension），记作 $\operatorname{codim} W$。

!!! theorem "定理 13A.5 (第一同构定理)"
    设 $T: V \to U$ 为线性映射。则存在唯一的同构

    $$
    \widetilde{T}: V/\ker T \xrightarrow{\;\sim\;} \operatorname{im} T,
    $$

    使得 $\widetilde{T}(v + \ker T) = T(v)$。即 $V/\ker T \cong \operatorname{im} T$。

??? proof "证明"
    **良定义性**：若 $v + \ker T = v' + \ker T$，则 $v - v' \in \ker T$，故 $T(v) = T(v')$。

    **线性性**：$\widetilde{T}((u + \ker T) + (v + \ker T)) = \widetilde{T}((u+v) + \ker T) = T(u+v) = T(u) + T(v) = \widetilde{T}(u + \ker T) + \widetilde{T}(v + \ker T)$。标量乘法类似。

    **单射**：$\widetilde{T}(v + \ker T) = 0$ 意味着 $T(v) = 0$，即 $v \in \ker T$，故 $v + \ker T = 0 + \ker T$。

    **满射**：对任意 $T(v) \in \operatorname{im} T$，$\widetilde{T}(v + \ker T) = T(v)$。

    **唯一性**：由 $\widetilde{T} \circ \pi = T$ 唯一确定。$\blacksquare$

!!! example "例 13A.3"
    设 $T: \mathbb{R}^3 \to \mathbb{R}^2$ 定义为 $T(x, y, z) = (x + y, y + z)$。

    $\ker T = \{(x, y, z) : x + y = 0, y + z = 0\} = \{(t, -t, t) : t \in \mathbb{R}\}$，维数为 $1$。

    $\operatorname{im} T = \mathbb{R}^2$（验证 $T$ 是满射）。

    第一同构定理给出 $\mathbb{R}^3/\ker T \cong \mathbb{R}^2$，且 $\dim(\mathbb{R}^3/\ker T) = 3 - 1 = 2 = \dim \mathbb{R}^2$。

!!! theorem "定理 13A.6 (第二同构定理)"
    设 $U, W$ 为 $V$ 的子空间。则

    $$
    (U + W)/W \cong U/(U \cap W).
    $$

    特别地，$\dim(U + W) - \dim W = \dim U - \dim(U \cap W)$。

??? proof "证明"
    定义 $\varphi: U \to (U + W)/W$，$\varphi(u) = u + W$。$\varphi$ 是线性映射。

    **满射**：对 $(u + w) + W \in (U+W)/W$，$(u + w) + W = u + W = \varphi(u)$。

    **核**：$\varphi(u) = W$ 当且仅当 $u \in W$，即 $u \in U \cap W$。故 $\ker \varphi = U \cap W$。

    由第一同构定理，$U/(U \cap W) \cong (U + W)/W$。$\blacksquare$

!!! theorem "定理 13A.7 (第三同构定理)"
    设 $W \subseteq U$ 为 $V$ 的子空间。则

    $$
    (V/W)/(U/W) \cong V/U.
    $$

??? proof "证明"
    定义 $\varphi: V/W \to V/U$，$\varphi(v + W) = v + U$。

    **良定义性**：若 $v + W = v' + W$，则 $v - v' \in W \subseteq U$，故 $v + U = v' + U$。

    **线性性**和**满射性**显然。

    **核**：$\varphi(v + W) = U$ 当且仅当 $v \in U$，即 $v + W \in U/W$。故 $\ker \varphi = U/W$。

    由第一同构定理，$(V/W)/(U/W) \cong V/U$。$\blacksquare$

---

## 13A.3 线性泛函与对偶空间

<div class="context-flow" markdown>

**核心问题**：线性映射 $V \to \mathbb{F}$ 有何特殊性？ → 线性泛函 → 所有线性泛函构成对偶空间 $V^*$

</div>

!!! definition "定义 13A.4 (线性泛函)"
    设 $V$ 为 $\mathbb{F}$ 上的向量空间。线性映射 $\varphi: V \to \mathbb{F}$ 称为 $V$ 上的一个**线性泛函**（linear functional）。

!!! example "例 13A.4"
    以下均为线性泛函：

    1. $\varphi: \mathbb{R}^n \to \mathbb{R}$，$\varphi(x_1, \ldots, x_n) = a_1 x_1 + \cdots + a_n x_n$（坐标的线性组合）。
    2. $\varphi: C[0,1] \to \mathbb{R}$，$\varphi(f) = \int_0^1 f(t)\, dt$（定积分）。
    3. $\varphi: \mathbb{F}[x]_{\leq n} \to \mathbb{F}$，$\varphi(p) = p(c)$（在点 $c$ 处求值）。

!!! definition "定义 13A.5 (对偶空间)"
    $V$ 上所有线性泛函的集合，配上逐点加法和标量乘法，构成向量空间

    $$
    V^* = \mathcal{L}(V, \mathbb{F}) = \operatorname{Hom}(V, \mathbb{F}),
    $$

    称为 $V$ 的**对偶空间**（dual space）或**代数对偶空间**。

!!! theorem "定理 13A.8 (对偶空间的维数)"
    若 $V$ 为有限维向量空间，则

    $$
    \dim V^* = \dim V.
    $$

??? proof "证明"
    $V^* = \mathcal{L}(V, \mathbb{F})$。由线性映射空间的维数公式，$\dim \mathcal{L}(V, \mathbb{F}) = \dim V \cdot \dim \mathbb{F} = \dim V \cdot 1 = \dim V$。$\blacksquare$

!!! theorem "定理 13A.9 (线性泛函的核)"
    设 $\varphi \in V^*$ 为非零线性泛函。则 $\ker \varphi$ 是 $V$ 的余维数为 $1$ 的子空间（超平面），即

    $$
    \dim \ker \varphi = \dim V - 1.
    $$

    反之，$V$ 的每个余维数为 $1$ 的子空间都是某个非零线性泛函的核。

??? proof "证明"
    $\varphi \neq 0$ 意味着 $\operatorname{im} \varphi = \mathbb{F}$（因为 $\operatorname{im} \varphi$ 是 $\mathbb{F}$ 的子空间，非零则为 $\mathbb{F}$ 本身）。由维数公式，$\dim \ker \varphi = \dim V - \dim \operatorname{im} \varphi = \dim V - 1$。

    反之，设 $\dim W = \dim V - 1$。取 $v_0 \in V \setminus W$，则 $V = W \oplus \operatorname{span}(v_0)$。定义 $\varphi(w + cv_0) = c$，则 $\varphi$ 是非零线性泛函，$\ker \varphi = W$。$\blacksquare$

!!! example "例 13A.5"
    在 $\mathbb{R}^3$ 中，线性泛函 $\varphi(x, y, z) = 2x - y + 3z$ 的核为

    $$
    \ker \varphi = \{(x, y, z) : 2x - y + 3z = 0\},
    $$

    这是 $\mathbb{R}^3$ 中过原点的一个平面（$2$ 维子空间，余维数为 $1$）。

---

## 13A.4 对偶基

<div class="context-flow" markdown>

**核心问题**：$V$ 的基如何确定 $V^*$ 的基？ → 对偶基的定义与构造 → 对偶基是坐标提取映射

</div>

!!! definition "定义 13A.6 (对偶基)"
    设 $\mathcal{B} = \{e_1, e_2, \ldots, e_n\}$ 为有限维向量空间 $V$ 的一组基。$V^*$ 的**对偶基**（dual basis）$\mathcal{B}^* = \{e_1^*, e_2^*, \ldots, e_n^*\}$ 定义为满足

    $$
    e_i^*(e_j) = \delta_{ij} = \begin{cases} 1 & \text{若 } i = j, \\ 0 & \text{若 } i \neq j \end{cases}
    $$

    的线性泛函。

!!! theorem "定理 13A.10 (对偶基的存在性与唯一性)"
    对偶基 $\mathcal{B}^*$ 存在、唯一，且构成 $V^*$ 的一组基。

??? proof "证明"
    **存在性与唯一性**：对每个 $i$，$e_i^*$ 由它在基 $\{e_1, \ldots, e_n\}$ 上的值唯一确定。定义 $e_i^*(e_j) = \delta_{ij}$ 并线性延拓，这给出了良定义的线性泛函。

    **是基**：设 $\sum c_i e_i^* = 0$（零泛函）。对 $e_j$ 求值：$\sum c_i e_i^*(e_j) = c_j = 0$。故 $e_1^*, \ldots, e_n^*$ 线性无关。由于 $\dim V^* = n$，它们构成基。$\blacksquare$

!!! theorem "定理 13A.11 (对偶基的坐标解释)"
    设 $v = \sum_{i=1}^n a_i e_i \in V$。则

    $$
    e_i^*(v) = a_i.
    $$

    即 $e_i^*$ 是提取 $v$ 在基 $\mathcal{B}$ 下第 $i$ 个坐标的映射。

    进而，$V^*$ 中任意线性泛函 $\varphi = \sum_{i=1}^n \varphi(e_i) e_i^*$。

??? proof "证明"
    $e_i^*(v) = e_i^*\left(\sum_j a_j e_j\right) = \sum_j a_j e_i^*(e_j) = \sum_j a_j \delta_{ij} = a_i$。

    对于第二个断言，设 $\varphi \in V^*$，令 $\psi = \sum_i \varphi(e_i) e_i^*$。对基向量 $e_j$：$\psi(e_j) = \sum_i \varphi(e_i) \delta_{ij} = \varphi(e_j)$。故 $\psi = \varphi$。$\blacksquare$

!!! example "例 13A.6"
    $V = \mathbb{R}^3$，标准基 $\{e_1, e_2, e_3\}$。对偶基为

    $$
    e_1^*(x, y, z) = x, \quad e_2^*(x, y, z) = y, \quad e_3^*(x, y, z) = z.
    $$

    即三个坐标函数。

!!! example "例 13A.7"
    $V = \mathbb{R}^2$，基 $\{v_1, v_2\} = \{(1, 1), (1, -1)\}$。求对偶基。

    设 $v_1^*(x, y) = ax + by$。由 $v_1^*(1, 1) = 1$ 和 $v_1^*(1, -1) = 0$，得 $a + b = 1$，$a - b = 0$，解得 $a = b = \frac{1}{2}$。

    设 $v_2^*(x, y) = cx + dy$。由 $v_2^*(1, 1) = 0$ 和 $v_2^*(1, -1) = 1$，得 $c + d = 0$，$c - d = 1$，解得 $c = \frac{1}{2}$，$d = -\frac{1}{2}$。

    因此 $v_1^*(x, y) = \frac{x + y}{2}$，$v_2^*(x, y) = \frac{x - y}{2}$。

---

## 13A.5 零化子

<div class="context-flow" markdown>

**核心问题**：$V^*$ 中哪些泛函在给定子空间 $U$ 上为零？ → 零化子 $U^0$ → 零化子的维数公式 → 零化子建立了 $V$ 的子空间与 $V^*$ 的子空间之间的对偶关系

</div>

!!! definition "定义 13A.7 (零化子)"
    设 $U$ 为 $V$ 的子空间。$U$ 的**零化子**（annihilator）定义为

    $$
    U^0 = \{\varphi \in V^* : \varphi(u) = 0 \text{ 对所有 } u \in U\}.
    $$

    $U^0$ 是 $V^*$ 的子空间。

!!! theorem "定理 13A.12 (零化子的维数)"
    设 $V$ 为有限维向量空间，$U$ 为 $V$ 的子空间。则

    $$
    \dim U^0 = \dim V - \dim U.
    $$

??? proof "证明"
    取 $U$ 的基 $\{e_1, \ldots, e_k\}$ 并扩充为 $V$ 的基 $\{e_1, \ldots, e_k, e_{k+1}, \ldots, e_n\}$。设 $\{e_1^*, \ldots, e_n^*\}$ 为对偶基。

    断言 $U^0 = \operatorname{span}\{e_{k+1}^*, \ldots, e_n^*\}$。

    一方面，对 $i > k$ 和 $j \leq k$：$e_i^*(e_j) = 0$，故 $e_i^* \in U^0$。

    另一方面，设 $\varphi = \sum_{i=1}^n c_i e_i^* \in U^0$。对 $j \leq k$：$\varphi(e_j) = c_j = 0$。故 $\varphi = \sum_{i=k+1}^n c_i e_i^*$。

    因此 $U^0 = \operatorname{span}\{e_{k+1}^*, \ldots, e_n^*\}$，$\dim U^0 = n - k = \dim V - \dim U$。$\blacksquare$

!!! theorem "定理 13A.13 (零化子的对偶性质)"
    设 $V$ 为有限维向量空间，$U, W$ 为 $V$ 的子空间。则：

    1. $U \subseteq W \Leftrightarrow W^0 \subseteq U^0$。
    2. $(U + W)^0 = U^0 \cap W^0$。
    3. $(U \cap W)^0 = U^0 + W^0$。
    4. $(U^0)^0 = U$（在自然同构 $V \cong V^{**}$ 下）。

??? proof "证明"
    1. 若 $U \subseteq W$ 且 $\varphi \in W^0$，则 $\varphi$ 在 $W$ 上为零，从而在 $U$ 上也为零，故 $\varphi \in U^0$。反之用维数计算。

    2. $\varphi \in (U + W)^0$ 当且仅当 $\varphi(u + w) = 0$ 对所有 $u \in U, w \in W$，当且仅当 $\varphi(u) = 0$ 对所有 $u \in U$ 且 $\varphi(w) = 0$ 对所有 $w \in W$，当且仅当 $\varphi \in U^0 \cap W^0$。

    3. 由 (2)，$(U^0 \cap W^0)^0 = ((U+W)^0)^0 = U + W$（利用 (4)）。两边取零化子再利用 (4)，或直接用维数计算：

    $$
    \dim(U^0 + W^0) = \dim U^0 + \dim W^0 - \dim(U^0 \cap W^0).
    $$

    由 (2)，$\dim(U^0 \cap W^0) = \dim(U + W)^0 = n - \dim(U+W)$。代入上式并化简得 $\dim(U^0 + W^0) = n - \dim(U \cap W) = \dim(U \cap W)^0$。再验证包含关系即得。

    4. 见 13A.7 节双对偶的讨论。$\blacksquare$

!!! example "例 13A.8"
    $V = \mathbb{R}^3$，$U = \operatorname{span}\{(1, 0, 1), (0, 1, 1)\}$。求 $U^0$。

    设 $\varphi(x, y, z) = ax + by + cz \in U^0$。由 $\varphi(1, 0, 1) = a + c = 0$ 和 $\varphi(0, 1, 1) = b + c = 0$，得 $a = -c$，$b = -c$。取 $c = 1$：$\varphi(x, y, z) = -x - y + z$。

    $$
    U^0 = \operatorname{span}\{(-1, -1, 1)\} = \{\lambda(-x - y + z) : \lambda \in \mathbb{R}\}.
    $$

    验证：$\dim U^0 = 1 = 3 - 2 = \dim V - \dim U$。

!!! example "例 13A.9"
    商空间与零化子的关系：$V/U$ 的对偶空间 $(V/U)^*$ 与 $U^0$ 自然同构。

    定义 $\Phi: U^0 \to (V/U)^*$，$\Phi(\varphi)(v + U) = \varphi(v)$。这是良定义的（$\varphi \in U^0$ 保证了与代表元无关），且是同构。这可以由维数相等（$\dim U^0 = \dim V - \dim U = \dim(V/U) = \dim(V/U)^*$）及单射性得到。

---

## 13A.6 转置映射（对偶映射）

<div class="context-flow" markdown>

**核心问题**：线性映射 $T: V \to W$ 如何诱导 $W^*$ 到 $V^*$ 的映射？ → 转置映射 $T^*$ → 转置映射反转方向 → 与矩阵转置的关系

</div>

!!! definition "定义 13A.8 (转置映射)"
    设 $T: V \to W$ 为线性映射。$T$ 的**转置映射**（transpose / dual map）$T^*: W^* \to V^*$ 定义为

    $$
    T^*(\psi) = \psi \circ T, \quad \psi \in W^*.
    $$

    即对 $v \in V$，$(T^*\psi)(v) = \psi(T(v))$。

!!! note "注"
    注意方向的反转：$T: V \to W$ 诱导 $T^*: W^* \to V^*$（逆变函子）。

!!! theorem "定理 13A.14 (转置映射的性质)"
    设 $S, T: V \to W$，$R: W \to U$ 为线性映射，$\lambda \in \mathbb{F}$。则：

    1. $T^*$ 是线性映射。
    2. $(S + T)^* = S^* + T^*$。
    3. $(\lambda T)^* = \lambda T^*$。
    4. $(R \circ T)^* = T^* \circ R^*$（注意顺序反转）。
    5. $(\operatorname{id}_V)^* = \operatorname{id}_{V^*}$。

??? proof "证明"
    1. $T^*(\psi_1 + \psi_2) = (\psi_1 + \psi_2) \circ T = \psi_1 \circ T + \psi_2 \circ T = T^*\psi_1 + T^*\psi_2$。标量乘法类似。

    4. 对 $\varphi \in U^*$：$(R \circ T)^*(\varphi) = \varphi \circ (R \circ T) = (\varphi \circ R) \circ T = T^*(R^*(\varphi)) = (T^* \circ R^*)(\varphi)$。

    其余类似。$\blacksquare$

!!! theorem "定理 13A.15 (转置映射的核与像)"
    设 $T: V \to W$ 为有限维向量空间之间的线性映射。则：

    1. $\ker T^* = (\operatorname{im} T)^0$。
    2. $\operatorname{im} T^* = (\ker T)^0$。
    3. $\operatorname{rank} T^* = \operatorname{rank} T$。

??? proof "证明"
    1. $\psi \in \ker T^*$ 当且仅当 $T^*\psi = 0$，即 $\psi \circ T = 0$，即 $\psi(T(v)) = 0$ 对所有 $v \in V$，即 $\psi$ 在 $\operatorname{im} T$ 上为零，即 $\psi \in (\operatorname{im} T)^0$。

    2. 先证 $\operatorname{im} T^* \subseteq (\ker T)^0$：若 $\varphi = T^*\psi$，则对 $v \in \ker T$，$\varphi(v) = \psi(T(v)) = \psi(0) = 0$，故 $\varphi \in (\ker T)^0$。

    维数计算：$\dim \operatorname{im} T^* = \dim W^* - \dim \ker T^* = \dim W - \dim(\operatorname{im} T)^0 = \dim W - (\dim W - \dim \operatorname{im} T) = \operatorname{rank} T$。而 $\dim(\ker T)^0 = \dim V - \dim \ker T = \operatorname{rank} T$。维数相等且有包含关系，故相等。

    3. 由 (2)，$\operatorname{rank} T^* = \dim \operatorname{im} T^* = \operatorname{rank} T$。$\blacksquare$

!!! theorem "定理 13A.16 (转置映射与矩阵转置)"
    设 $\mathcal{B}$ 为 $V$ 的基，$\mathcal{C}$ 为 $W$ 的基，$\mathcal{B}^*$ 和 $\mathcal{C}^*$ 为对应的对偶基。若 $T: V \to W$ 关于 $\mathcal{B}, \mathcal{C}$ 的矩阵为 $A$，则 $T^*: W^* \to V^*$ 关于 $\mathcal{C}^*, \mathcal{B}^*$ 的矩阵为 $A^T$（$A$ 的转置）。

??? proof "证明"
    设 $A = (a_{ij})$，即 $T(e_j) = \sum_i a_{ij} f_i$。则

    $$
    (T^* f_i^*)(e_j) = f_i^*(T(e_j)) = f_i^*\left(\sum_k a_{kj} f_k\right) = a_{ij}.
    $$

    另一方面，若 $T^* f_i^* = \sum_j b_{ji} e_j^*$，则 $(T^* f_i^*)(e_j) = b_{ji}$。

    因此 $b_{ji} = a_{ij}$，即 $T^*$ 的矩阵为 $A^T$。$\blacksquare$

!!! example "例 13A.10"
    设 $T: \mathbb{R}^2 \to \mathbb{R}^3$，$T(x, y) = (x + y, 2x, y)$。关于标准基，$T$ 的矩阵为

    $$
    A = \begin{pmatrix} 1 & 1 \\ 2 & 0 \\ 0 & 1 \end{pmatrix}.
    $$

    $T^*: (\mathbb{R}^3)^* \to (\mathbb{R}^2)^*$ 关于对偶标准基的矩阵为

    $$
    A^T = \begin{pmatrix} 1 & 2 & 0 \\ 1 & 0 & 1 \end{pmatrix}.
    $$

    验证：$T^*(\varphi)(x, y) = \varphi(x + y, 2x, y)$。若 $\varphi(a, b, c) = \alpha a + \beta b + \gamma c$，则 $T^*\varphi(x, y) = \alpha(x+y) + 2\beta x + \gamma y = (\alpha + 2\beta)x + (\alpha + \gamma)y$，与 $A^T$ 的矩阵乘法一致。

---

## 13A.7 双对偶空间与自然同构

<div class="context-flow" markdown>

**核心问题**：$V^{**} = (V^*)^*$ 与 $V$ 是什么关系？ → 自然同构（canonical isomorphism） → "自然"意味着不依赖基的选取

</div>

!!! definition "定义 13A.9 (求值映射)"
    对 $v \in V$，定义 $\hat{v}: V^* \to \mathbb{F}$ 为

    $$
    \hat{v}(\varphi) = \varphi(v), \quad \varphi \in V^*.
    $$

    $\hat{v}$ 是 $V^*$ 上的线性泛函，即 $\hat{v} \in V^{**}$。

!!! theorem "定理 13A.17 (自然同构)"
    映射 $\iota: V \to V^{**}$，$\iota(v) = \hat{v}$，是**单射线性映射**。当 $V$ 为有限维时，$\iota$ 是**同构**，称为**自然同构**（canonical isomorphism）。

??? proof "证明"
    **线性性**：$\widehat{u + v}(\varphi) = \varphi(u + v) = \varphi(u) + \varphi(v) = \hat{u}(\varphi) + \hat{v}(\varphi)$，故 $\widehat{u+v} = \hat{u} + \hat{v}$。标量乘法类似。

    **单射**：设 $\iota(v) = 0$，即 $\hat{v} = 0$。则 $\varphi(v) = 0$ 对所有 $\varphi \in V^*$。若 $v \neq 0$，将 $v$ 扩充为 $V$ 的基的一部分 $\{v, e_2, \ldots, e_n\}$，取对偶基中的 $v^*$（第一个对偶基向量）则 $v^*(v) = 1 \neq 0$，矛盾。故 $v = 0$。

    **有限维时是同构**：$\dim V^{**} = \dim V^* = \dim V$，单射线性映射在维数相等时自动是同构。$\blacksquare$

!!! note "注"
    "自然"一词有精确的数学含义（范畴论中的自然变换）：$\iota$ 的定义不依赖于基的选取，且与线性映射相容。具体地，若 $T: V \to W$，则 $T^{**} \circ \iota_V = \iota_W \circ T$（交换图）。

    相比之下，$V$ 与 $V^*$ 之间虽然同构（有限维时维数相等），但不存在自然同构——任何同构都依赖于基的选取或内积的引入。

!!! theorem "定理 13A.18 (对偶基的对偶基)"
    设 $\{e_1, \ldots, e_n\}$ 为 $V$ 的基，$\{e_1^*, \ldots, e_n^*\}$ 为对偶基。在自然同构 $\iota: V \to V^{**}$ 下，

    $$
    \iota(e_i) = e_i^{**}
    $$

    其中 $\{e_1^{**}, \ldots, e_n^{**}\}$ 是 $\{e_1^*, \ldots, e_n^*\}$ 的对偶基。即 $V$ 的基在自然同构下恰好映射为 $V^*$ 的对偶基的对偶基。

??? proof "证明"
    需要验证 $\hat{e}_i(e_j^*) = \delta_{ij}$。由定义，$\hat{e}_i(e_j^*) = e_j^*(e_i) = \delta_{ji} = \delta_{ij}$。因此 $\hat{e}_i$ 正是 $\{e_1^*, \ldots, e_n^*\}$ 的对偶基中的第 $i$ 个元素。$\blacksquare$

!!! example "例 13A.11"
    自然同构的应用：证明零化子的性质 (4)——$(U^0)^0 = U$。

    在自然同构 $\iota: V \to V^{**}$ 下，$(U^0)^0 = \{\Phi \in V^{**} : \Phi(\varphi) = 0 \text{ 对所有 } \varphi \in U^0\}$。

    若 $v \in U$，则对任意 $\varphi \in U^0$，$\hat{v}(\varphi) = \varphi(v) = 0$，故 $\iota(U) \subseteq (U^0)^0$。

    维数计算：$\dim(U^0)^0 = \dim V^* - \dim U^0 = \dim V - (\dim V - \dim U) = \dim U = \dim \iota(U)$。

    故 $\iota(U) = (U^0)^0$，即在自然同构下 $(U^0)^0$ 对应于 $U$。

---

## 13A.8 对偶空间在有限维中的应用

<div class="context-flow" markdown>

**核心问题**：对偶空间理论如何应用于具体的线性代数问题？ → 线性方程组的对偶描述 → 基变换的对偶视角 → 秩的对偶刻画

</div>

### 线性方程组的对偶描述

!!! theorem "定理 13A.19 (解空间与零化子)"
    齐次线性方程组 $A\mathbf{x} = \mathbf{0}$（$A$ 为 $m \times n$ 矩阵）的解空间 $S$ 满足

    $$
    S = \{v \in \mathbb{F}^n : \varphi_i(v) = 0, \; i = 1, \ldots, m\},
    $$

    其中 $\varphi_i \in (\mathbb{F}^n)^*$ 是由 $A$ 的第 $i$ 行定义的线性泛函。设 $W = \operatorname{span}\{\varphi_1, \ldots, \varphi_m\} \subseteq (\mathbb{F}^n)^*$，则 $S = W^0$（将 $W$ 视为 $(\mathbb{F}^n)^*$ 的子空间，$W^0$ 在自然同构下对应 $\mathbb{F}^n$ 的子空间）。

    因此 $\dim S = n - \dim W = n - \operatorname{rank} A$。

??? proof "证明"
    $v \in S$ 当且仅当 $Av = 0$，即 $A$ 的每一行与 $v$ 的内积为零，即 $\varphi_i(v) = 0$ 对所有 $i$。$v$ 在所有 $\varphi_i$ 上取值为零等价于 $v$ 在 $\operatorname{span}\{\varphi_i\}$ 上取值为零。

    通过自然同构 $\iota: \mathbb{F}^n \to (\mathbb{F}^n)^{**}$，$S$ 恰好是 $W$ 在 $(\mathbb{F}^n)^{**}$ 中的零化子在 $\iota$ 下的原像，即 $S = W^0$（广义意义下）。$\dim S = n - \dim W = n - \operatorname{rank} A$。$\blacksquare$

### 基变换的对偶视角

!!! theorem "定理 13A.20 (对偶基的变换公式)"
    设 $\mathcal{B} = \{e_1, \ldots, e_n\}$ 和 $\mathcal{B}' = \{e_1', \ldots, e_n'\}$ 为 $V$ 的两组基，过渡矩阵为 $P$（即 $e_j' = \sum_i p_{ij} e_i$）。则对偶基的过渡矩阵为 $(P^{-1})^T$：

    $$
    e_j'^* = \sum_i \left[(P^{-1})^T\right]_{ij} e_i^* = \sum_i (P^{-1})_{ji} e_i^*.
    $$

??? proof "证明"
    设 $e_j'^* = \sum_i q_{ij} e_i^*$。由 $e_j'^*(e_k') = \delta_{jk}$：

    $$
    \delta_{jk} = e_j'^*(e_k') = \sum_i q_{ij} e_i^*\left(\sum_l p_{lk} e_l\right) = \sum_i q_{ij} p_{ik} = (Q^T P)_{jk}.
    $$

    故 $Q^T P = I$，即 $Q^T = P^{-1}$，$Q = (P^{-1})^T$。$\blacksquare$

### 秩的对偶刻画

!!! theorem "定理 13A.21 (行秩等于列秩的对偶证明)"
    $m \times n$ 矩阵 $A$ 的行秩等于列秩。

??? proof "证明"
    将 $A$ 视为线性映射 $T: \mathbb{F}^n \to \mathbb{F}^m$。列秩 $= \dim \operatorname{im} T = \operatorname{rank} T$。

    行秩 $= \dim \operatorname{span}\{\text{$A$ 的行}\}$。$A$ 的行定义了 $(\mathbb{F}^n)^*$ 中的线性泛函 $\varphi_1, \ldots, \varphi_m$。注意 $A^T: (\mathbb{F}^m)^* \to (\mathbb{F}^n)^*$ 的像恰好是 $\operatorname{span}\{\varphi_1, \ldots, \varphi_m\}$（因为 $T^*(f_i^*) = f_i^* \circ T = \varphi_i$，其中 $\{f_1^*, \ldots, f_m^*\}$ 是 $(\mathbb{F}^m)^*$ 的标准基）。故行秩 $= \operatorname{rank} T^*$。

    由定理 13A.15 (3)，$\operatorname{rank} T^* = \operatorname{rank} T$。$\blacksquare$

!!! example "例 13A.12"
    对偶理论在 Lagrange 插值中的体现。

    设 $V = \mathbb{F}[x]_{\leq n}$（次数 $\leq n$ 的多项式空间），取 $n + 1$ 个不同点 $c_0, c_1, \ldots, c_n \in \mathbb{F}$。定义求值泛函 $\varepsilon_i: V \to \mathbb{F}$，$\varepsilon_i(p) = p(c_i)$。

    $\{\varepsilon_0, \varepsilon_1, \ldots, \varepsilon_n\}$ 是 $V^*$ 的一组基（因为若 $\sum a_i \varepsilon_i = 0$，则多项式 $\sum a_i \ell_i(x) = 0$，其中 $\ell_i$ 为 Lagrange 基多项式，故 $a_i = 0$）。

    这组基的对偶基（在 $V^{**} \cong V$ 下）恰好是 **Lagrange 基多项式** $\ell_0, \ell_1, \ldots, \ell_n$，因为 $\varepsilon_i(\ell_j) = \ell_j(c_i) = \delta_{ij}$。

    这正是第 0 章多项式插值的对偶空间解释。

!!! example "例 13A.13"
    设 $T: V \to V$ 为线性算子（$V$ 有限维），且 $T^2 = T$（$T$ 是幂等算子）。证明 $(T^*)^2 = T^*$。

    由定理 13A.14 (4)：$(T^*)^2 = T^* \circ T^* = (T \circ T)^* = (T^2)^* = T^*$。

    进一步，$V = \ker T \oplus \operatorname{im} T$，而对偶空间中 $V^* = \ker T^* \oplus \operatorname{im} T^*$。由定理 13A.15：

    - $\ker T^* = (\operatorname{im} T)^0$；
    - $\operatorname{im} T^* = (\ker T)^0$。

    这展示了幂等分解在对偶空间中的完美对称性。
