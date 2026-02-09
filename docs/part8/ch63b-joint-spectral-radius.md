# 第 63B 章 联合谱半径

<div class="context-flow" markdown>

**前置**：特征值与谱半径(Ch6) · 矩阵范数(Ch15) · 同时三角化(Ch63A)

**本章脉络**：联合谱半径定义(Rota-Strang) → 良定义性与范数无关性 → 基本性质 → 广义谱半径 → Berger-Wang 定理 → 有限性猜想 → 不可判定性 → Barabanov 范数与极值范数理论 → 多面体范数与 SOS/SDP 方法 → 下谱半径 → 切换线性系统稳定性 → CQLF → 多重 Lyapunov 函数 → 约束切换 → 小波正则性

**延伸**：联合谱半径在自动控制（切换线性系统稳定性）、小波理论（细分格式收敛性与光滑性）、编码理论（码的容量）中有核心应用

</div>

当我们从单一矩阵的谱理论走向**矩阵集合**的谱理论时，单一矩阵的谱半径 $\rho(A)$ 自然推广为一组矩阵的**联合谱半径**。对于单个矩阵，Gelfand 公式 $\rho(A) = \lim_{k\to\infty} \|A^k\|^{1/k}$ 刻画了 $A^k$ 的增长率；而对于一组矩阵 $\Sigma = \{A_1, \ldots, A_m\}$，联合谱半径 $\rho(\Sigma)$ 刻画了所有可能的乘积 $A_{i_1} A_{i_2} \cdots A_{i_k}$ 的**最大增长率**。

这一概念由 Rota 和 Strang 于 1960 年引入，此后在动力系统、自动控制、小波理论和编码理论中获得了广泛应用。Berger-Wang 定理（1992）建立了联合谱半径与广义谱半径的相等性，Barabanov 范数理论提供了深层的几何刻画，而有限性猜想的否定回答揭示了联合谱半径计算的本质困难。

---

## 63B.1 联合谱半径的定义

### 动机与定义

对于单个矩阵 $A$，谱半径 $\rho(A) = \max |\lambda_i(A)|$ 通过 Gelfand 公式等于 $\lim_{k\to\infty} \|A^k\|^{1/k}$。当我们有一组矩阵时，需要考虑所有可能的乘积序列的增长行为。

!!! definition "定义 63B.1 (联合谱半径, Rota-Strang 1960)"
    设 $\Sigma = \{A_1, A_2, \ldots, A_m\} \subset M_n(\mathbb{C})$ 是一个有限矩阵集合，$\|\cdot\|$ 是 $M_n(\mathbb{C})$ 上的某个矩阵范数。定义 $\Sigma$ 的**联合谱半径**（joint spectral radius, JSR）为
    $$
    \rho(\Sigma) = \lim_{k \to \infty} \max_{(i_1, \ldots, i_k) \in \{1,\ldots,m\}^k} \|A_{i_1} A_{i_2} \cdots A_{i_k}\|^{1/k}
    $$

    记 $\Sigma^k = \{A_{i_1} \cdots A_{i_k} : i_j \in \{1,\ldots,m\}\}$ 为所有长度 $k$ 的乘积的集合，则等价地
    $$
    \rho(\Sigma) = \lim_{k \to \infty} \left(\max_{A \in \Sigma^k} \|A\|\right)^{1/k}
    $$

### 良定义性与范数无关性

!!! theorem "定理 63B.1 (联合谱半径的良定义性)"
    定义 63B.1 中的极限存在，且与所选矩阵范数 $\|\cdot\|$ 无关。

??? proof "证明"
    **极限存在性（Fekete 引理）**：定义 $a_k = \log \max_{A \in \Sigma^k} \|A\|$。对于任意 $p, q \geq 1$，长度为 $p + q$ 的乘积可以拆分为长度 $p$ 和长度 $q$ 的乘积的组合：
    $$
    \max_{A \in \Sigma^{p+q}} \|A\| \leq \max_{B \in \Sigma^p} \|B\| \cdot \max_{C \in \Sigma^q} \|C\|
    $$
    这是因为对任意 $A = B \cdot C$（$B \in \Sigma^p$，$C \in \Sigma^q$），有 $\|A\| = \|BC\| \leq \|B\| \cdot \|C\|$。

    取对数得 $a_{p+q} \leq a_p + a_q$，即序列 $\{a_k\}$ 是**次可加的**（subadditive）。

    由 Fekete 引理，次可加序列满足
    $$
    \lim_{k\to\infty} \frac{a_k}{k} = \inf_{k \geq 1} \frac{a_k}{k}
    $$
    且该极限存在（可能为 $-\infty$，但由于 $\max_{A \in \Sigma^k} \|A\| \geq 0$，实际上 $a_k$ 有下界，取值在 $[-\infty, +\infty)$ 中，当 $\Sigma$ 非空时极限有限）。

    因此 $\lim_{k\to\infty} \left(\max_{A \in \Sigma^k} \|A\|\right)^{1/k} = \exp\left(\lim_{k\to\infty} \frac{a_k}{k}\right)$ 存在。

    **范数无关性**：设 $\|\cdot\|_a$ 和 $\|\cdot\|_b$ 是两个矩阵范数。由有限维空间上范数的等价性，存在常数 $c, C > 0$ 使得
    $$
    c\|A\|_a \leq \|A\|_b \leq C\|A\|_a, \quad \forall A \in M_n(\mathbb{C})
    $$

    设 $\hat{\rho}_k^{(a)} = \left(\max_{A \in \Sigma^k} \|A\|_a\right)^{1/k}$，类似定义 $\hat{\rho}_k^{(b)}$。则
    $$
    c^{1/k} \hat{\rho}_k^{(a)} \leq \hat{\rho}_k^{(b)} \leq C^{1/k} \hat{\rho}_k^{(a)}
    $$

    令 $k \to \infty$，由于 $c^{1/k} \to 1$ 和 $C^{1/k} \to 1$，两个范数给出相同的极限值。$\blacksquare$

!!! example "例 63B.1"
    设 $\Sigma = \{A_1, A_2\}$，其中
    $$
    A_1 = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}, \quad
    A_2 = \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}
    $$

    $\rho(A_1) = \rho(A_2) = 1$（均为 Jordan 块），但乘积
    $$
    A_1 A_2 = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix}
    $$
    的特征值为 $\frac{3 \pm \sqrt{5}}{2}$，谱半径 $\rho(A_1 A_2) = \frac{3+\sqrt{5}}{2} \approx 2.618$。

    因此 $\rho(\Sigma) \geq \rho(A_1 A_2)^{1/2} = \phi \approx 1.618$（黄金比例）。实际上可以证明 $\rho(\Sigma) = \phi$。

---

## 63B.2 基本性质

!!! theorem "定理 63B.2 (联合谱半径的基本性质)"
    设 $\Sigma, \Sigma' \subset M_n(\mathbb{C})$ 为有限矩阵集合。则：

    1. **非负性**：$\rho(\Sigma) \geq 0$，且 $\rho(\Sigma) = 0$ 当且仅当 $\Sigma$ 生成的矩阵半群中每个元素都是幂零的；
    2. **齐次性**：$\rho(\alpha \Sigma) = |\alpha| \rho(\Sigma)$，其中 $\alpha\Sigma = \{\alpha A : A \in \Sigma\}$；
    3. **次可乘性**：$\rho(\Sigma \cdot \Sigma') \leq \rho(\Sigma) \cdot \rho(\Sigma')$；
    4. **单调性**：若 $\Sigma \subseteq \Sigma'$，则 $\rho(\Sigma) \leq \rho(\Sigma')$；
    5. **连续性**：$\rho(\cdot)$ 作为 $(M_n(\mathbb{C}))^m$ 到 $\mathbb{R}_{\geq 0}$ 的函数是连续的；
    6. **退化为经典谱半径**：当 $\Sigma = \{A\}$ 时，$\rho(\Sigma) = \rho(A)$。

??? proof "证明"
    **性质 6（退化性）**：当 $\Sigma = \{A\}$ 时，$\Sigma^k = \{A^k\}$，因此
    $$
    \rho(\Sigma) = \lim_{k\to\infty} \|A^k\|^{1/k} = \rho(A)
    $$
    最后一步是 Gelfand 公式。

    **性质 2（齐次性）**：$(\alpha\Sigma)^k = \{\alpha^k B : B \in \Sigma^k\}$，因此
    $$
    \max_{B \in (\alpha\Sigma)^k} \|B\|^{1/k} = |\alpha| \cdot \max_{B \in \Sigma^k} \|B\|^{1/k}
    $$
    取极限即得。

    **性质 4（单调性）**：$\Sigma \subseteq \Sigma'$ 蕴含 $\Sigma^k \subseteq (\Sigma')^k$，因此最大值只可能增大。

    **性质 5（连续性）**：利用范数的连续性和极限运算的性质。更精确的证明需要利用次可加函数极限的一致收敛性。$\blacksquare$

!!! example "例 63B.2"
    设 $\Sigma = \{\alpha A : \alpha \in [0, 1]\}$，其中 $A$ 为固定矩阵。由单调性，$\rho(\Sigma) \geq \rho(\{A\}) = \rho(A)$（因 $A \in \Sigma$）。另一方面，长度 $k$ 的乘积为 $\alpha_1 \cdots \alpha_k A^k$，$\alpha_i \in [0,1]$，最大值在 $\alpha_i = 1$ 时达到。因此 $\rho(\Sigma) = \rho(A)$。

---

## 63B.3 广义谱半径与 Berger-Wang 定理

### 广义谱半径

!!! definition "定义 63B.2 (广义谱半径)"
    设 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{C})$。$\Sigma$ 的**广义谱半径**（generalized spectral radius）定义为
    $$
    \hat{\rho}(\Sigma) = \limsup_{k \to \infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$
    其中 $\rho(A)$ 是矩阵 $A$ 的经典谱半径。

!!! note "注"
    广义谱半径 $\hat{\rho}(\Sigma)$ 完全由矩阵乘积的**特征值**决定，不涉及范数选取。定义中使用 $\limsup$ 是因为极限不一定存在（但 Berger-Wang 定理将表明它等于 JSR，后者极限存在）。

### Berger-Wang 定理

!!! theorem "定理 63B.3 (Berger-Wang 定理, 1992)"
    设 $\Sigma \subset M_n(\mathbb{C})$ 是有界矩阵集合。则联合谱半径等于广义谱半径：
    $$
    \rho(\Sigma) = \hat{\rho}(\Sigma) = \limsup_{k \to \infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$

??? proof "证明"
    **方向一：$\hat{\rho}(\Sigma) \leq \rho(\Sigma)$。**

    对任意 $A \in \Sigma^k$，由谱半径与范数的关系 $\rho(A) \leq \|A\|$：
    $$
    \max_{A \in \Sigma^k} \rho(A)^{1/k} \leq \left(\max_{A \in \Sigma^k} \|A\|\right)^{1/k}
    $$
    取 $\limsup$ 得 $\hat{\rho}(\Sigma) \leq \rho(\Sigma)$。

    **方向二：$\rho(\Sigma) \leq \hat{\rho}(\Sigma)$（深刻方向）。**

    设 $\hat{\rho} = \hat{\rho}(\Sigma)$。目标：对任意 $\varepsilon > 0$，证明 $\rho(\Sigma) \leq \hat{\rho} + \varepsilon$。

    **步骤 1**：由 $\hat{\rho}$ 的定义，存在 $k_0$ 使得对所有 $k \geq k_0$：
    $$
    \max_{A \in \Sigma^k} \rho(A) \leq (\hat{\rho} + \varepsilon/2)^k
    $$

    **步骤 2**：利用 Gelfand 公式的矩阵半群推广，构造一个适当的范数。对于任意矩阵 $A$，$\rho(A) \leq r$ 蕴含对任意 $\delta > 0$，存在范数 $\|\cdot\|_\delta$ 使得 $\|A\|_\delta \leq r + \delta$。

    将此思想推广到矩阵集合：定义
    $$
    \|X\|_\varepsilon = \sum_{k=0}^{\infty} \frac{1}{(\hat{\rho} + \varepsilon)^k} \max_{A \in \Sigma^k} \|AX\|
    $$
    （约定 $\Sigma^0 = \{I\}$）。需要验证该级数收敛。

    由步骤 1，对 $k \geq k_0$，$\max_{A \in \Sigma^k} \|AX\| \leq C \cdot (\hat{\rho} + \varepsilon/2)^k \|X\|$（$C$ 是与有限个低阶项有关的常数，利用 Gelfand 公式），因此级数中的项以 $\left(\frac{\hat{\rho}+\varepsilon/2}{\hat{\rho}+\varepsilon}\right)^k$ 的速度衰减，级数收敛。

    **步骤 3**：验证 $\|\cdot\|_\varepsilon$ 是一个范数，且对所有 $A_i \in \Sigma$：
    $$
    \|A_i X\|_\varepsilon \leq (\hat{\rho} + \varepsilon) \|X\|_\varepsilon
    $$
    这直接由定义的移位性质得到。

    因此 $\rho(\Sigma) \leq \hat{\rho} + \varepsilon$。由 $\varepsilon$ 的任意性得 $\rho(\Sigma) \leq \hat{\rho}(\Sigma)$。$\blacksquare$

!!! example "例 63B.3"
    回到例 63B.1 中的 $\Sigma = \{A_1, A_2\}$。Berger-Wang 定理保证
    $$
    \rho(\Sigma) = \limsup_{k\to\infty} \max_{A \in \Sigma^k} \rho(A)^{1/k}
    $$
    由于 $\rho(A_1 A_2) = \phi^2$，取 $k=2$ 得 $\rho(\Sigma) \geq \phi$。进一步分析所有长度 $k$ 的乘积的特征值，可以严格证明 $\rho(\Sigma) = \phi$。

---

## 63B.4 有限性猜想及其否定

!!! definition "定义 63B.3 (有限性猜想)"
    **有限性猜想**（Finiteness Conjecture）断言：对任意有限矩阵集 $\Sigma$，存在有限长度的乘积 $A_{i_1} \cdots A_{i_k}$ 使得
    $$
    \rho(\Sigma) = \rho(A_{i_1} \cdots A_{i_k})^{1/k}
    $$
    即联合谱半径可由某个有限长度乘积的谱半径精确实现。

!!! theorem "定理 63B.4 (有限性猜想的否定)"
    有限性猜想不成立。存在一对 $2 \times 2$ 实矩阵 $\Sigma = \{A_0, A_1\}$，使得 $\rho(\Sigma)$ 不能由任何有限长度乘积的谱半径实现。

!!! note "注"
    反例由 Bousch-Mairesse（2002）和 Blondel-Theys-Vladimirov（2003）独立构造。然而对于许多重要的特殊情形，有限性性质仍然成立：

    - **非负矩阵对**：有限性猜想对 $\Sigma = \{A_1, A_2\}$（$A_1, A_2 \geq 0$）成立。
    - **可交换矩阵集**：$\rho(\Sigma) = \max_{A \in \Sigma} \rho(A)$，由 $k=1$ 的乘积实现。
    - **单参数矩阵族**的某些情形。

---

## 63B.5 不可判定性

!!! theorem "定理 63B.5 (JSR 的不可判定性)"
    以下问题是算法不可判定的：给定一组整数矩阵 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{Z})$，判断 $\rho(\Sigma) \leq 1$ 还是 $\rho(\Sigma) > 1$。

??? proof "证明"
    证明通过将**有界性问题**归约到经典的不可判定问题。具体地：

    **步骤 1**：矩阵集 $\Sigma$ 生成的半群 $\mathcal{S}(\Sigma) = \bigcup_{k \geq 1} \Sigma^k$ 有界（即 $\sup_{A \in \mathcal{S}(\Sigma)} \|A\| < \infty$）当且仅当 $\rho(\Sigma) \leq 1$（实际上等价于 $\rho(\Sigma) < 1$ 或 $\rho(\Sigma) = 1$ 且满足额外条件）。

    更精确地说，$\rho(\Sigma) \leq 1$ 等价于 $\sup_k \max_{A \in \Sigma^k} \|A\|^{1/k} \leq 1$，即乘积不指数增长。

    **步骤 2**：将 Post 对应问题（PCP）——一个已知不可判定的问题——编码为矩阵半群的有界性判定问题。PCP 的一个实例可以转化为一组整数矩阵 $\Sigma$，使得 PCP 有解当且仅当 $\rho(\Sigma) > 1$。

    因此判定 $\rho(\Sigma) \leq 1$ 是不可判定的。$\blacksquare$

!!! note "注"
    不可判定性意味着**不存在通用算法**在有限步内精确计算联合谱半径。然而，JSR 可以被任意精度地**近似**，这在实际应用中已经足够。

---

## 63B.6 Barabanov 范数与极值范数理论

### Barabanov 范数

!!! definition "定义 63B.4 (极值范数)"
    矩阵范数 $\|\cdot\|$ 称为矩阵集 $\Sigma$ 的**极值范数**（extremal norm），若
    $$
    \max_{A \in \Sigma} \|Ax\| = \rho(\Sigma) \|x\|, \quad \forall x \in \mathbb{R}^n
    $$
    即 $\Sigma$ 中的矩阵在此范数下的算子范数恰好等于 $\rho(\Sigma)$。

!!! theorem "定理 63B.6 (Barabanov 范数存在性)"
    设 $\Sigma \subset M_n(\mathbb{R})$ 是紧矩阵集，且 $\Sigma$ 不可约（即 $\Sigma$ 中的矩阵没有公共的非平凡不变子空间）。则 $\Sigma$ 存在极值范数。

    更具体地，存在一个范数 $\|\cdot\|_*$（称为 **Barabanov 范数**），使得
    $$
    \max_{A \in \Sigma} \|Ax\|_* = \rho(\Sigma) \|x\|_*, \quad \forall x \in \mathbb{R}^n
    $$

??? proof "证明"
    **构造方法**：定义序列
    $$
    \|x\|_k = \frac{1}{\rho(\Sigma)^k} \max_{A \in \Sigma^k} \|Ax\|
    $$
    其中 $\|\cdot\|$ 是任意初始范数。由 JSR 的定义，$\|x\|_k$ 有界。

    取 $\|x\|_* = \limsup_{k \to \infty} \|x\|_k$。可以验证 $\|\cdot\|_*$ 是一个范数（齐次性和三角不等式由逐点 $\limsup$ 保持），且满足
    $$
    \max_{A \in \Sigma} \|Ax\|_* = \rho(\Sigma) \|x\|_*
    $$

    不可约性条件保证了 $\|x\|_* > 0$ 对所有 $x \neq 0$ 成立（若存在 $x \neq 0$ 使得 $\|x\|_* = 0$，则 $V = \{x : \|x\|_* = 0\}$ 是 $\Sigma$ 的公共不变子空间，矛盾）。$\blacksquare$

### Protasov 定理

!!! theorem "定理 63B.7 (Protasov 定理)"
    设 $\Sigma \subset M_n(\mathbb{R})$ 是有限不可约矩阵集。则：

    1. Barabanov 范数的单位球 $B_* = \{x : \|x\|_* \leq 1\}$ 是严格凸的紧凸体；
    2. 对单位球面上的每个点 $x$（$\|x\|_* = 1$），存在 $A \in \Sigma$ 使得 $\|Ax\|_* = \rho(\Sigma)$；
    3. Barabanov 范数在 $C^0$ 等价类意义下不一定唯一，但其几何结构受 $\Sigma$ 的不可约结构严格约束。

!!! note "注"
    Barabanov 范数理论的一个重要推论是：若 $\Sigma$ 不可约且 $\rho(\Sigma) = 1$，则 Barabanov 范数的单位球是 $\Sigma$ 的**不变凸体**，即 $A B_* \subseteq B_*$ 对所有 $A \in \Sigma$ 成立，且边界上每个点都有至少一个 $A$ 使其映射回边界。

---

## 63B.7 计算方法

### 多面体范数方法

!!! definition "定义 63B.5 (多面体范数方法)"
    **多面体范数方法**通过寻找一个多面体范数 $\|\cdot\|_P$（其单位球是多面体 $P$）来逼近联合谱半径。若存在 $\gamma > 0$ 使得对所有 $A \in \Sigma$ 和所有 $x$，$\|Ax\|_P \leq \gamma \|x\|_P$，则 $\rho(\Sigma) \leq \gamma$。

    核心思想是寻找一个**不变多面体**（invariant polytope）$P$ 使得 $A_i P \subseteq \gamma P$ 对所有 $i$ 成立。

!!! theorem "定理 63B.8 (多面体逼近收敛性)"
    设 $\Sigma$ 是有限矩阵集。对任意 $\varepsilon > 0$，存在一个多面体范数 $\|\cdot\|_P$，使得
    $$
    \rho(\Sigma) \leq \max_{A \in \Sigma} \|A\|_P \leq \rho(\Sigma) + \varepsilon
    $$

### SOS/SDP 方法

!!! theorem "定理 63B.9 (SOS/SDP 上界方法)"
    **平方和**（Sum of Squares, SOS）方法利用半定规划（SDP）构造 Lyapunov 函数。对给定次数 $2d$，寻找正定矩阵 $Q \succ 0$ 使得
    $$
    A_i^T Q A_i \preceq \gamma^2 Q, \quad \forall i = 1, \ldots, m
    $$
    这是线性矩阵不等式（LMI）。若此 LMI 可行，则 $\rho(\Sigma) \leq \gamma$。

    对于更高次的 SOS Lyapunov 函数（次数 $2d > 2$），利用提升技术：令 $\tilde{\Sigma} = \{A \otimes \cdots \otimes A : A \in \Sigma^d\}$（$d$ 次 Kronecker 积），则二次 LMI 应用于 $\tilde{\Sigma}$ 等价于 $2d$ 次 SOS Lyapunov 函数。

    随着 $d \to \infty$，SOS 上界收敛到 $\rho(\Sigma)$。

!!! example "例 63B.4"
    对于 $\Sigma = \{A_1, A_2\}$，其中 $A_1 = \begin{pmatrix} 0.5 & 0.6 \\ 0 & 0.5 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0.5 & 0 \\ 0.6 & 0.5 \end{pmatrix}$。

    使用 $d=1$（二次 Lyapunov 函数），求解 SDP $A_i^T Q A_i \preceq \gamma^2 Q$，$Q \succ 0$。取 $Q = I$ 检验：
    $$
    A_1^T A_1 = \begin{pmatrix} 0.25 & 0.30 \\ 0.30 & 0.61 \end{pmatrix}
    $$
    最大特征值约为 $0.786$，因此 $\gamma_1 \approx 0.887$。类似地 $\gamma_2 \approx 0.887$。可得上界 $\rho(\Sigma) \leq 0.887$。

    下界：$\rho(A_1 A_2)^{1/2} \approx 0.781$。

---

## 63B.8 下谱半径与联合谱子半径

!!! definition "定义 63B.6 (下谱半径 / 联合谱子半径)"
    设 $\Sigma \subset M_n(\mathbb{C})$ 是有限矩阵集合。$\Sigma$ 的**下谱半径**（lower spectral radius）或**联合谱子半径**（joint spectral subradius）定义为
    $$
    \check{\rho}(\Sigma) = \lim_{k \to \infty} \min_{A \in \Sigma^k} \|A\|^{1/k}
    $$

    等价地，$\check{\rho}(\Sigma) = \lim_{k \to \infty} \left(\min_{A \in \Sigma^k} \|A\|\right)^{1/k}$。

!!! theorem "定理 63B.10 (下谱半径的性质)"
    1. $\check{\rho}(\Sigma)$ 的定义良好且与范数无关；
    2. $0 \leq \check{\rho}(\Sigma) \leq \rho(\Sigma)$；
    3. $\check{\rho}(\Sigma) = 0$ 当且仅当 $\Sigma$ 中存在有限乘积为零矩阵；
    4. 下谱半径的 Berger-Wang 类型定理也成立：$\check{\rho}(\Sigma) = \liminf_{k\to\infty} \min_{A \in \Sigma^k} \rho(A)^{1/k}$。

!!! note "注"
    下谱半径刻画了矩阵乘积的**最小增长率**。在切换系统语境中，$\check{\rho}(\Sigma) > 1$ 意味着无论如何切换，状态都会发散（系统**绝对不稳定**）。$\check{\rho}(\Sigma) < 1$ 意味着存在使系统稳定的切换策略。

!!! example "例 63B.5"
    设 $A_1 = \begin{pmatrix} 2 & 0 \\ 0 & 0.5 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0.5 & 0 \\ 0 & 2 \end{pmatrix}$。

    $\rho(\Sigma) = 2$（取乘积 $A_1^k$，谱半径为 $2^k$）。
    $\check{\rho}(\Sigma) = 1$（乘积 $A_1 A_2 = I$，谱半径为 $1$，无法更小）。

---

## 63B.9 切换线性系统的稳定性

### 基本模型

!!! definition "定义 63B.7 (离散时间切换线性系统)"
    给定矩阵集 $\Sigma = \{A_1, \ldots, A_m\} \subset M_n(\mathbb{R})$ 和切换信号 $\sigma: \mathbb{N} \to \{1, \ldots, m\}$，**离散时间切换线性系统**定义为
    $$
    x_{k+1} = A_{\sigma(k)} x_k, \quad k = 0, 1, 2, \ldots
    $$

!!! definition "定义 63B.8 (绝对渐近稳定性)"
    切换系统称为**绝对渐近稳定**（absolutely asymptotically stable），若对**所有**初始条件 $x_0$ 和**所有**切换信号 $\sigma$，$x_k \to 0$（$k \to \infty$）。

### 稳定性判据

!!! theorem "定理 63B.11 (切换系统稳定性与 JSR)"
    离散时间切换线性系统绝对渐近稳定，当且仅当 $\rho(\Sigma) < 1$。

??? proof "证明"
    **充分性**：若 $\rho(\Sigma) < 1$，取 $\gamma \in (\rho(\Sigma), 1)$。由极值范数理论（或直接由 JSR 定义），存在矩阵范数 $\|\cdot\|_*$ 使得 $\|A\|_* \leq \gamma$ 对所有 $A \in \Sigma$ 成立。

    则对任意切换信号 $\sigma$ 和初始条件 $x_0$：
    $$
    \|x_k\|_* = \|A_{\sigma(k-1)} \cdots A_{\sigma(0)} x_0\|_* \leq \gamma^k \|x_0\|_* \to 0
    $$

    **必要性**：若 $\rho(\Sigma) \geq 1$，由 JSR 定义，存在长度 $k$ 的乘积序列使得 $\|A_{i_1^{(k)}} \cdots A_{i_k^{(k)}}\|^{1/k} \to \rho(\Sigma) \geq 1$。选取使范数最大化的切换序列，取合适的初始条件 $x_0$，可以构造不收敛到零的轨道。$\blacksquare$

### 公共二次 Lyapunov 函数

!!! definition "定义 63B.9 (CQLF)"
    矩阵集 $\Sigma = \{A_1, \ldots, A_m\}$ 具有**公共二次 Lyapunov 函数**（Common Quadratic Lyapunov Function, CQLF），若存在 $P \succ 0$ 使得
    $$
    A_i^T P A_i \prec P, \quad \forall i = 1, \ldots, m
    $$
    此时 $V(x) = x^T P x$ 是公共 Lyapunov 函数。

!!! theorem "定理 63B.12 (CQLF 与稳定性的关系)"
    1. 若 $\Sigma$ 具有 CQLF，则 $\rho(\Sigma) < 1$；
    2. 反之不然：存在 $\rho(\Sigma) < 1$ 但不存在 CQLF 的矩阵集。
    3. 对于两个矩阵的情形 $\Sigma = \{A_1, A_2\}$，CQLF 存在当且仅当 $A_1^T P A_1 \prec P$ 和 $A_2^T P A_2 \prec P$ 对某个 $P \succ 0$ 同时成立。当 $A_1, A_2$ 均为正矩阵时，CQLF 存在性可以用乘积 $A_1 A_2^{-1}$ 的特征值条件刻画。

### 多重 Lyapunov 函数

!!! definition "定义 63B.10 (分段二次 Lyapunov 函数)"
    **多重 Lyapunov 函数**（multiple Lyapunov functions）方法为每个模式 $i$ 分配一个二次函数 $V_i(x) = x^T P_i x$（$P_i \succ 0$），要求
    $$
    A_i^T P_j A_i \prec P_i, \quad \forall (i, j) \text{ 允许的切换对}
    $$

    更一般地，**分段二次 Lyapunov 函数**允许 Lyapunov 函数在状态空间的不同区域取不同的二次形式。

!!! theorem "定理 63B.13 (多重 Lyapunov 函数的必要性)"
    若 $\rho(\Sigma) < 1$，则对充分精细的划分，存在分段二次 Lyapunov 函数证明稳定性。即多重/分段二次 Lyapunov 函数方法在理论上是**充要的**（尽管可能需要非常多的片段）。

!!! example "例 63B.6"
    设 $A_1 = \begin{pmatrix} 0.5 & 0.3 \\ 0 & 0.4 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0.4 & 0 \\ 0.3 & 0.5 \end{pmatrix}$。

    取 $P = I$ 检验 CQLF：
    $$
    A_1^T A_1 = \begin{pmatrix} 0.25 & 0.15 \\ 0.15 & 0.25 \end{pmatrix} \prec I, \quad
    A_2^T A_2 = \begin{pmatrix} 0.25 & 0.15 \\ 0.15 & 0.34 \end{pmatrix} \prec I
    $$
    $P = I$ 是 CQLF，因此 $\rho(\Sigma) < 1$，系统在所有切换下稳定。

!!! example "例 63B.7"
    设 $A_1 = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，$A_2 = \begin{pmatrix} 0 & 0 \\ 1 & 0 \end{pmatrix}$。

    $\rho(A_1) = \rho(A_2) = 0$（幂零矩阵），但 $A_1 A_2 = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$，$\rho(A_1 A_2) = 1$。

    交替切换 $\sigma = (1, 2, 1, 2, \ldots)$ 给出 $x_{2k} = (A_1 A_2)^k x_0$，不收敛到零。虽然每个矩阵谱半径为零，但 $\rho(\Sigma) = 1$，系统不稳定。

---

## 63B.10 约束切换

### 驻留时间约束

!!! definition "定义 63B.11 (驻留时间切换)"
    **驻留时间**（dwell-time）约束要求切换信号 $\sigma$ 在每个模式上至少停留 $\tau$ 步才能切换到下一个模式。即若 $\sigma(k) \neq \sigma(k-1)$，则 $\sigma(k) = \sigma(k+1) = \cdots = \sigma(k+\tau-1)$。

    驻留时间 $\tau$ 下的联合谱半径定义为
    $$
    \rho_\tau(\Sigma) = \rho(\Sigma^\tau) = \rho(\{A_1^\tau, \ldots, A_m^\tau\})^{1/\tau}
    $$

!!! theorem "定理 63B.14 (驻留时间稳定性)"
    对任意矩阵集 $\Sigma$，若每个 $A_i$ 都是 Schur 稳定的（$\rho(A_i) < 1$），则存在足够大的驻留时间 $\tau_0$ 使得 $\rho_\tau(\Sigma) < 1$ 对所有 $\tau \geq \tau_0$ 成立。

### 图约束切换

!!! definition "定义 63B.12 (图约束切换)"
    给定有向图 $G = (V, E)$，$V = \{1, \ldots, m\}$。**图约束切换**要求切换信号 $\sigma$ 满足 $(\sigma(k), \sigma(k+1)) \in E$ 对所有 $k$ 成立。

    图约束下的联合谱半径 $\rho_G(\Sigma)$ 定义为仅考虑 $G$ 中允许的路径对应的乘积的增长率。

!!! theorem "定理 63B.15 (图约束 JSR 的计算)"
    图约束联合谱半径可以通过构造**提升矩阵**来计算。定义 $mn \times mn$ 的分块矩阵
    $$
    \tilde{A} = (B_{ij})_{i,j=1}^m, \quad B_{ij} = \begin{cases} A_j & \text{若 } (i,j) \in E \\ 0 & \text{否则} \end{cases}
    $$
    则 $\rho_G(\Sigma) = \rho(\tilde{A})$，归结为单个矩阵的谱半径计算。

---

## 63B.11 小波正则性中的联合谱半径

### 细分格式与转移矩阵

!!! definition "定义 63B.13 (二进细分格式)"
    给定细分系数 $\{c_0, c_1, \ldots, c_N\} \subset \mathbb{R}$，定义两个**转移矩阵** $T_0$ 和 $T_1$：
    $$
    (T_\epsilon)_{ij} = c_{2i - j + \epsilon}, \quad \epsilon = 0, 1
    $$

### 收敛性与光滑性

!!! theorem "定理 63B.16 (细分格式的收敛性, Daubechies-Lagarias)"
    二进细分格式收敛（即存在连续极限函数），当且仅当
    $$
    \rho(\{T_0|_V, T_1|_V\}) < 1
    $$
    其中 $V$ 是差分空间（由向量 $e_i - e_{i+1}$ 张成的子空间）。

!!! theorem "定理 63B.17 (尺度函数的 Sobolev/Holder 正则性)"
    设细分格式的转移矩阵在 $s$ 阶差分空间上的限制为 $\{T_0^{(s)}, T_1^{(s)}\}$。尺度函数 $\phi$ 的 Holder 正则性指数为
    $$
    \alpha = s - \log_2 \rho(\{T_0^{(s)}, T_1^{(s)}\})
    $$
    其中 $s$ 是格式的逼近阶数。

    更一般地，$\phi$ 属于 Sobolev 空间 $W^{s,p}$，当且仅当
    $$
    \rho_p(\{T_0^{(s)}, T_1^{(s)}\}) < 2^{s - 1 + 1/p}
    $$

!!! example "例 63B.8"
    线性 B 样条的细分系数 $c = (1/2, 1, 1/2)$。转移矩阵在差分空间上的限制为
    $$
    T_0|_V = \begin{pmatrix} 1/2 \end{pmatrix}, \quad T_1|_V = \begin{pmatrix} 1/2 \end{pmatrix}
    $$
    $\rho(\{T_0|_V, T_1|_V\}) = 1/2 < 1$，格式收敛。Holder 正则性 $\alpha = 1 - \log_2(1/2) = 2$，但实际的 Holder 指数是 $1$（帽子函数是 Lipschitz 但不是 $C^2$），需要更仔细地区分 $s$ 的选取。

!!! example "例 63B.9"
    Daubechies 小波 $D_4$（4 系数）的细分系数为
    $$
    c_k = \frac{1}{4\sqrt{2}} \{1+\sqrt{3}, \; 3+\sqrt{3}, \; 3-\sqrt{3}, \; 1-\sqrt{3}\}
    $$
    计算限制转移矩阵在一阶差分空间上的联合谱半径，确定尺度函数 $\phi$ 的 Holder 正则性约为 $\alpha \approx 1.09$。

!!! note "注"
    联合谱半径在小波理论中的应用框架由 Daubechies 和 Lagarias（1991, 1992）建立。后续推广包括：

    - **多元小波**：膨胀矩阵确定多个转移矩阵；
    - **向量值细分**：控制点为向量，转移矩阵更大；
    - **非平稳细分**：每层使用不同规则，需计算无限矩阵集的 JSR。

---

## 习题

!!! question "习题 63B.1"
    计算 $\Sigma = \left\{\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}, \begin{pmatrix} 0 & 0 \\ 1 & 1 \end{pmatrix}\right\}$ 的联合谱半径。

!!! question "习题 63B.2"
    设 $\Sigma = \{\alpha A : \alpha \in [0, 1]\}$，其中 $A$ 为固定矩阵。证明 $\rho(\Sigma) = \rho(A)$。

!!! question "习题 63B.3"
    设 $A_1, A_2$ 是两个 $2 \times 2$ 幂零矩阵。证明 $\rho(\{A_1, A_2\}) = \sqrt{\rho(A_1 A_2)}$。

!!! question "习题 63B.4"
    证明：若 $\Sigma$ 中所有矩阵两两可交换，则 $\rho(\Sigma) = \max_{A \in \Sigma} \rho(A)$。

!!! question "习题 63B.5"
    给出一个切换系统的例子，其中每个 $A_i$ 都满足 $\rho(A_i) < 1$，但 $\rho(\Sigma) \geq 1$。

!!! question "习题 63B.6"
    证明 Berger-Wang 定理中 $\hat{\rho}(\Sigma) \leq \rho(\Sigma)$ 的方向。

!!! question "习题 63B.7"
    对于系数 $c = (1/4, 1/2, 1/4)$ 的细分格式，计算转移矩阵并判断收敛性。

!!! question "习题 63B.8"
    设 $\Sigma = \{A_1, A_2\}$，$A_1 = \begin{pmatrix} a & 0 \\ 0 & b \end{pmatrix}$，$A_2 = \begin{pmatrix} 0 & c \\ d & 0 \end{pmatrix}$。计算 $\rho(\Sigma)$。

!!! question "习题 63B.9"
    证明：若 $\rho(\Sigma) < 1$，则级数 $\sum_{k=0}^{\infty} \sum_{A \in \Sigma^k} A$ 绝对收敛。

!!! question "习题 63B.10"
    设 $\Sigma = \{A_1, A_2\} \subset M_2(\mathbb{R})$。构造一个例子说明 $\rho(\Sigma) < 1$ 但不存在公共二次 Lyapunov 函数。（提示：利用旋转矩阵。）

!!! question "习题 63B.11"
    证明图约束联合谱半径满足 $\rho_G(\Sigma) \leq \rho(\Sigma)$，且当 $G$ 是完全图时等号成立。

!!! question "习题 63B.12"
    设 $\Sigma$ 中每个矩阵 $A_i$ 的谱半径满足 $\rho(A_i) < 1$。证明存在驻留时间 $\tau$ 使得 $\rho(\{A_1^\tau, \ldots, A_m^\tau\})^{1/\tau} < 1$。
