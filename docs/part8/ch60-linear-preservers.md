# 第 60 章 线性保持问题

<div class="context-flow" markdown>

**前置**：线性变换(Ch5) · 行列式(Ch3) · 特征值(Ch6) · 矩阵运算(Ch2)

**本章脉络**：线性保持问题的一般框架 → 行列式保持（Frobenius 定理） → 秩保持 → 谱保持 → 正定保持 → 可逆性保持 → 交换性保持 → 现代方向

**延伸**：线性保持问题从 1897 年 Frobenius 开创至今仍是活跃的研究方向；量子信息中量子信道的完全正性保持、算子代数中的自同构刻画都是线性保持问题的推广

</div>

线性保持问题（Linear Preserver Problem, LPP）是矩阵理论中最古老且最持久的研究方向之一。其核心问题极为自然：给定矩阵空间上的某个性质 $\mathcal{P}$（如行列式的值、秩、谱、正定性等），刻画所有保持该性质的**线性映射** $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$。

这一方向的起源可以追溯到 1897 年，当时 Frobenius 证明了保持行列式的线性映射必然具有非常特殊的形式。此后一个多世纪中，线性保持问题发展为矩阵理论的一个主要分支，与算子代数、量子信息、多线性代数等领域产生了深刻的联系。

线性保持问题的哲学意义在于：它揭示了矩阵的各种性质之间的**内在刚性**——看似不同的性质往往由相同的保持映射所保持，这些映射的形式高度受限。

---

## 60.1 线性保持问题的一般框架

<div class="context-flow" markdown>

**核心问题**：什么是线性保持问题？如何系统地分类和研究各种保持问题？

</div>

!!! definition "定义 60.1 (线性保持映射)"
    设 $\mathbb{F}$ 为域（通常为 $\mathbb{R}$ 或 $\mathbb{C}$），$M_n(\mathbb{F})$ 为 $n \times n$ 矩阵空间。设 $\mathcal{P}$ 是 $M_n(\mathbb{F})$ 上的某个性质或函数。**线性 $\mathcal{P}$-保持映射**是线性映射 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$，使得

    $$\mathcal{P}(\phi(A)) = \mathcal{P}(A) \quad \forall\, A \in M_n(\mathbb{F}).$$

    或者更一般地，$A$ 满足 $\mathcal{P}$ 当且仅当 $\phi(A)$ 满足 $\mathcal{P}$。

!!! definition "定义 60.2 (标准形式)"
    最常出现的保持映射形式是以下两类**标准形式**：

    **(a) 类型 I**（相似变换型）：$\phi(A) = MAN$，其中 $M, N \in M_n(\mathbb{F})$ 是固定的可逆矩阵。

    **(b) 类型 II**（含转置型）：$\phi(A) = MA^\top N$ 或 $\phi(A) = MA^* N$（复数域上取共轭转置）。

    大量线性保持问题的结论都是：保持映射必然是类型 I 或类型 II 的某种特殊形式。

!!! definition "定义 60.3 (保持问题的分类)"
    线性保持问题可按保持的对象分为几大类：

    - **标量值函数保持**：行列式、迹、谱半径等。
    - **集合保持**：秩-$k$ 矩阵集合、可逆矩阵集合、幂等矩阵集合等。
    - **子空间保持**：对称矩阵、上三角矩阵等子空间。
    - **关系保持**：交换性、相似性等二元关系。
    - **序保持**：正定序（Loewner 序）、谱序等。

!!! example "例 60.1"
    一些具体的线性保持问题实例：

    - **行列式保持**：$\det(\phi(A)) = \det(A)$，$\forall A$。
    - **秩保持**：$\mathrm{rank}(\phi(A)) = \mathrm{rank}(A)$，$\forall A$。
    - **可逆性保持**：$\phi(A)$ 可逆 $\Leftrightarrow$ $A$ 可逆。
    - **幂等保持**：$A^2 = A \Rightarrow \phi(A)^2 = \phi(A)$。
    - **正定保持**：$A \succ 0 \Rightarrow \phi(A) \succ 0$。

---

## 60.2 Frobenius 行列式保持定理

<div class="context-flow" markdown>

**核心问题**：什么样的线性映射保持行列式？Frobenius 在 1897 年得到了什么结果？

</div>

!!! theorem "定理 60.1 (Frobenius 行列式保持定理, 1897)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 是线性映射，满足

    $$\det(\phi(A)) = \det(A), \quad \forall\, A \in M_n(\mathbb{C}).$$

    则存在可逆矩阵 $M, N \in M_n(\mathbb{C})$，$\det(MN) = 1$，使得

    $$\phi(A) = MAN \quad \text{或} \quad \phi(A) = MA^\top N.$$

??? proof "证明"
    这是一个深刻的定理，完整证明较长。我们给出主要步骤。

    **第一步：$\phi$ 保持秩-1 矩阵。** 设 $A$ 是秩-1 矩阵。对任意 $B$，考虑多项式

    $$p(t) = \det(B + tA) = \det(B) + t \cdot \mathrm{adj}(B) \cdot A \cdot \text{(某些项)} + \cdots$$

    利用行列式的多线性性和 $\det(\phi(B + tA)) = \det(B + tA)$，可以证明 $\phi(A)$ 也是秩-1 矩阵。

    具体地，$\mathrm{rank}(A) = 1$ 意味着 $\det(A + tB)$ 作为 $t$ 的多项式次数至多为 1（对适当选择的 $B$）。由 $\det(\phi(A + tB)) = \det(A + tB)$，$\phi(A)$ 也满足相同的次数约束，故 $\mathrm{rank}(\phi(A)) = 1$。

    **第二步：秩-1 保持映射的刻画。** 由 Marcus-Moyls 定理（定理 60.3），秩-1 保持的线性映射必然具有以下形式之一：

    $$\phi(xy^\top) = (Mx)(Ny)^\top = Mxy^\top N^\top,$$

    或

    $$\phi(xy^\top) = (My)(Nx)^\top = Myx^\top N^\top = M(xy^\top)^\top N^\top.$$

    因此 $\phi(A) = MAN^\top$ 或 $\phi(A) = MA^\top N^\top$。

    **第三步：确定 $\det(MN^\top) = 1$。** 由 $\det(\phi(I)) = \det(I) = 1$，得 $\det(MN^\top) = 1$（第一种情形）或 $\det(M)\det(N) = 1$（第二种情形）。

    将 $N^\top$ 重新记为 $N$（使得 $\det(MN) = 1$），即得 Frobenius 定理。$\blacksquare$

!!! theorem "定理 60.2 (行列式乘法保持)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 线性，满足 $\det(\phi(A)) = c \cdot \det(A)$，$\forall A$，其中 $c \neq 0$ 是常数。则

    $$\phi(A) = MAN \quad \text{或} \quad \phi(A) = MA^\top N,$$

    其中 $\det(MN) = c$。

!!! example "例 60.2"
    **转置映射保持行列式。** $\phi(A) = A^\top$ 满足 $\det(A^\top) = \det(A)$。这对应 Frobenius 定理中的第二种形式，取 $M = I$，$N = I$。

    **相似变换保持行列式。** $\phi(A) = PAP^{-1}$ 满足 $\det(PAP^{-1}) = \det(A)$。这对应第一种形式，$M = P$，$N = P^{-1}$，$\det(MN) = 1$。

    **非标准形式不可能。** 映射 $\phi(A) = \mathrm{tr}(A) \cdot I$ 是线性的，但 $\det(\mathrm{tr}(A) \cdot I) = (\mathrm{tr}(A))^n$，这一般不等于 $\det(A)$（$n \geq 2$ 时）。

---

## 60.3 秩保持

<div class="context-flow" markdown>

**核心问题**：保持矩阵秩的线性映射有怎样的结构？秩-1 保持是否足以确定映射的形式？

</div>

!!! theorem "定理 60.3 (Marcus-Moyls 定理, 1959)"
    设 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$（$\mathbb{F}$ 为无限域，$n \geq 2$）是线性映射，满足

    $$\mathrm{rank}(A) = 1 \implies \mathrm{rank}(\phi(A)) = 1.$$

    则存在可逆矩阵 $M, N$ 使得

    $$\phi(A) = MAN \quad \text{或} \quad \phi(A) = MA^\top N.$$

??? proof "证明"
    **核心思路。** 秩-1 矩阵恰好是形如 $uv^\top$ 的矩阵（$u \neq 0$，$v \neq 0$）。

    **第一步：$\phi$ 保持秩-1 矩阵空间的结构。** 固定非零 $u$，考虑 $\{uv^\top : v \in \mathbb{F}^n\}$，这是一个 $n$ 维子空间，其非零元素都是秩-1。$\phi$ 将其映到一个子空间，非零元素仍为秩-1。

    可以证明：一个子空间的非零元素全是秩-1 矩阵，当且仅当该子空间形如 $\{u'v^\top : v \in \mathbb{F}^n\}$（固定列空间）或 $\{uv'^\top : u \in \mathbb{F}^n\}$（固定行空间）中的某种。

    **第二步：确定映射对秩-1 矩阵的作用形式。** 利用第一步的结构，可以证明 $\phi(e_i e_j^\top) = M e_{\sigma(i)} e_{\tau(j)}^\top N'$（某种行列置换），进而组合论证得到全局形式。

    **第三步：线性扩展。** 每个矩阵都是秩-1 矩阵的线性组合，因此在秩-1 矩阵上确定的形式唯一地扩展到全体矩阵。$\blacksquare$

!!! theorem "定理 60.4 (秩保持的完整刻画)"
    设 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$ 线性，$n \geq 2$。以下条件等价：

    (a) $\phi$ 保持秩：$\mathrm{rank}(\phi(A)) = \mathrm{rank}(A)$，$\forall A$。

    (b) $\phi$ 保持秩-1：$\mathrm{rank}(A) = 1 \implies \mathrm{rank}(\phi(A)) = 1$，且 $\phi$ 是可逆的。

    (c) 存在可逆矩阵 $M, N$ 使得 $\phi(A) = MAN$ 或 $\phi(A) = MA^\top N$。

??? proof "证明"
    $(c) \Rightarrow (a)$：$\mathrm{rank}(MAN) = \mathrm{rank}(A)$（因 $M, N$ 可逆）。$\mathrm{rank}(MA^\top N) = \mathrm{rank}(A^\top) = \mathrm{rank}(A)$。

    $(a) \Rightarrow (b)$：秩保持显然蕴含秩-1 保持。$\phi$ 的可逆性由 $\phi$ 保持秩的事实推出（若 $\phi(A) = 0$，则 $\mathrm{rank}(\phi(A)) = 0 = \mathrm{rank}(A)$，故 $A = 0$）。

    $(b) \Rightarrow (c)$：由 Marcus-Moyls 定理（定理 60.3）。$\blacksquare$

!!! example "例 60.3"
    **秩-$k$ 保持（$k > 1$）不蕴含标准形式。** 映射 $\phi(A) = A + c\,\mathrm{tr}(A)\, I$（$c$ 为适当选择的常数）保持秩-$n$ 矩阵（可逆矩阵）但不保持秩-1。这说明秩-1 保持的特殊地位。

---

## 60.4 谱保持

<div class="context-flow" markdown>

**核心问题**：保持矩阵谱（特征值集合）的线性映射是什么形式？

</div>

!!! definition "定义 60.4 (谱保持映射)"
    线性映射 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 称为**谱保持的**，如果

    $$\sigma(\phi(A)) = \sigma(A), \quad \forall\, A \in M_n(\mathbb{C}),$$

    其中 $\sigma(A)$ 是 $A$ 的特征值多重集。

!!! theorem "定理 60.5 (谱保持映射的刻画)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 是**单射**线性映射，满足 $\sigma(\phi(A)) = \sigma(A)$，$\forall A$。则存在可逆矩阵 $M$ 使得

    $$\phi(A) = MAM^{-1} \quad \text{或} \quad \phi(A) = MA^\top M^{-1}.$$

??? proof "证明"
    **证明策略。** 分三步。

    **第一步：$\phi$ 保持迹和行列式。** $\mathrm{tr}(A) = \sum \lambda_i = \sum \sigma_i(\phi(A)) = \mathrm{tr}(\phi(A))$。类似地 $\det(A) = \prod \lambda_i = \det(\phi(A))$。更一般地，$\phi$ 保持特征多项式的所有系数。

    **第二步：$\phi$ 保持幂零矩阵。** $A$ 幂零当且仅当 $\sigma(A) = \{0\}$。因此 $\phi$ 保持幂零矩阵集合。

    **第三步：利用幂零保持推出标准形式。** 这需要精细的结构分析。幂零矩阵的保持映射已被完全刻画（需要用到 Jordan 结构理论），最终得到 $\phi$ 必须是 $A \mapsto MAM^{-1}$ 或 $A \mapsto MA^\top M^{-1}$。

    关键的技术引理是：如果 $\phi$ 保持幂零性且保持迹，那么对任意两个矩阵 $A, B$，$\mathrm{tr}(\phi(A)\phi(B)) = \mathrm{tr}(AB)$ 或 $\mathrm{tr}(\phi(A)\phi(B)) = \mathrm{tr}(AB^\top)$。由此通过 Killing 形式的非退化性推出映射的全局形式。$\blacksquare$

!!! theorem "定理 60.6 (谱半径保持)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 是满射线性映射，满足 $\rho(\phi(A)) = \rho(A)$（$\rho$ 为谱半径），$\forall A$。则

    $$\phi(A) = e^{i\theta} MAM^{-1} \quad \text{或} \quad \phi(A) = e^{i\theta} MA^\top M^{-1},$$

    其中 $\theta \in \mathbb{R}$，$M$ 可逆。

!!! example "例 60.4"
    **转置保持谱但不是相似变换。** 映射 $\phi(A) = A^\top$ 满足 $\sigma(A^\top) = \sigma(A)$（因为 $A$ 和 $A^\top$ 有相同的特征多项式）。但 $A^\top$ 一般不相似于 $A$（在实矩阵中可以有不同的 Jordan 结构），不过特征值作为多重集确实相同。

    这说明谱保持映射允许"转置"这种非相似变换型。

---

## 60.5 正定保持与正映射

<div class="context-flow" markdown>

**核心问题**：什么样的线性映射将正半定矩阵映射到正半定矩阵？正映射与完全正映射有什么区别？

</div>

!!! definition "定义 60.5 (正映射)"
    线性映射 $\phi: M_n(\mathbb{C}) \to M_m(\mathbb{C})$ 称为**正映射**（positive map），如果

    $$A \succeq 0 \implies \phi(A) \succeq 0.$$

!!! example "例 60.5"
    以下都是正映射：

    (a) $\phi(A) = MAM^*$（$M$ 为任意矩阵）。验证：$A \succeq 0 \Rightarrow x^* MAM^* x = (M^*x)^* A(M^*x) \geq 0$。

    (b) $\phi(A) = A^\top$（转置映射）。$A \succeq 0 \Rightarrow A^\top \succeq 0$。

    (c) $\phi(A) = \mathrm{tr}(A) \cdot I_m$。$A \succeq 0 \Rightarrow \mathrm{tr}(A) \geq 0 \Rightarrow \mathrm{tr}(A) \cdot I_m \succeq 0$。

!!! definition "定义 60.6 (完全正映射)"
    线性映射 $\phi: M_n \to M_m$ 称为**完全正的**（completely positive, CP），如果对所有 $k \geq 1$，扩展映射

    $$\phi \otimes \mathrm{id}_k: M_n \otimes M_k \to M_m \otimes M_k, \quad A \otimes B \mapsto \phi(A) \otimes B$$

    是正映射。等价地，$(\phi \otimes \mathrm{id}_k)(X) \succeq 0$，只要 $X \in M_{nk}$ 是正半定的（将 $M_{nk}$ 视为 $M_n \otimes M_k$）。

!!! theorem "定理 60.7 (转置不是完全正映射)"
    转置映射 $\phi(A) = A^\top$ 是正映射但**不是**完全正映射。

??? proof "证明"
    考虑 $n = 2$，$k = 2$。构造 $M_4$ 中的正半定矩阵

    $$X = \begin{pmatrix}1 & 0 & 0 & 1\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\1 & 0 & 0 & 1\end{pmatrix} = \begin{pmatrix}E_{11} & E_{11}\\E_{11} & E_{11}\end{pmatrix}_{\text{block}}$$

    实际上取 $X = |\psi\rangle\langle\psi|$，其中 $|\psi\rangle = (1, 0, 0, 1)^\top / \sqrt{2}$（Bell 态的变体），则 $X \succeq 0$。

    对 $X$ 的每个 $2 \times 2$ 块施加转置：

    $$(\phi \otimes \mathrm{id}_2)(X) = \begin{pmatrix}1 & 0 & 0 & 0\\0 & 0 & 1 & 0\\0 & 1 & 0 & 0\\0 & 0 & 0 & 1\end{pmatrix}.$$

    这个矩阵的特征值为 $\{1, 1, 1, -1\}$，存在负特征值，故不是正半定的。

    因此 $\phi \otimes \mathrm{id}_2$ 不保持正半定性，$\phi$ 不是完全正映射。$\blacksquare$

!!! theorem "定理 60.8 (Choi 定理, 1975)"
    线性映射 $\phi: M_n \to M_m$ 是完全正的，当且仅当存在矩阵 $V_1, \ldots, V_r \in M_{m \times n}$ 使得

    $$\phi(A) = \sum_{i=1}^r V_i A V_i^*, \quad \forall\, A \in M_n.$$

    这称为 **Kraus 表示**（或算子和表示）。

??? proof "证明"
    **"$\Leftarrow$"**：若 $\phi(A) = \sum V_i A V_i^*$ 且 $A \succeq 0$，则每个 $V_i A V_i^* \succeq 0$，故 $\phi(A) \succeq 0$。对 $\phi \otimes \mathrm{id}_k$ 同理验证，因此 $\phi$ 完全正。

    **"$\Rightarrow$"（Choi 矩阵方法）**：定义 **Choi 矩阵** $C_\phi \in M_{mn}$：

    $$C_\phi = \sum_{i,j=1}^n \phi(E_{ij}) \otimes E_{ij} = (\phi \otimes \mathrm{id}_n)\!\Bigl(\sum_{i,j} E_{ij} \otimes E_{ij}\Bigr).$$

    注意到 $\sum_{i,j} E_{ij} \otimes E_{ij}$ 是 $M_{n^2}$ 中的正半定矩阵（它是 $|\Omega\rangle\langle\Omega|$ 的非归一化版本，$|\Omega\rangle = \sum_i e_i \otimes e_i$）。

    如果 $\phi$ 完全正，则 $C_\phi \succeq 0$。对 $C_\phi$ 做谱分解 $C_\phi = \sum_i \lambda_i u_i u_i^*$，$\lambda_i \geq 0$。将每个 $u_i \in \mathbb{C}^{mn}$ 重塑为 $m \times n$ 矩阵 $V_i$（乘以 $\sqrt{\lambda_i}$），可以验证 $\phi(A) = \sum V_i A V_i^*$。$\blacksquare$

!!! theorem "定理 60.9 (正映射的结构——Storevig 分解)"
    每个正映射 $\phi: M_n \to M_n$ 可以分解为

    $$\phi = \phi_1 + \phi_2 \circ \tau,$$

    其中 $\phi_1, \phi_2$ 是完全正映射，$\tau$ 是转置映射。但这种分解**不总是可行的**（仅在 $n = 2$ 时成立）。对一般 $n$，正映射的完整结构描述是开放问题。

---

## 60.6 可逆性保持

<div class="context-flow" markdown>

**核心问题**：什么样的线性映射保持矩阵的可逆性？

</div>

!!! definition "定义 60.7 (可逆性保持映射)"
    线性映射 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$ 称为**可逆性保持的**，如果

    $$A \text{ 可逆} \implies \phi(A) \text{ 可逆}.$$

!!! theorem "定理 60.10 (Dieudonne 定理, 1949)"
    设 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$（$\mathbb{F}$ 为无限域）是线性映射，保持可逆性。如果 $\phi$ 是满射的，则存在可逆矩阵 $M, N$ 使得

    $$\phi(A) = MAN \quad \text{或} \quad \phi(A) = MA^\top N.$$

??? proof "证明"
    **核心思路。** $\phi$ 保持可逆性等价于 $\phi$ 保持奇异性的补集。更有用的是等价条件：对所有 $A, B$，如果 $\det(A + tB) \not\equiv 0$（作为 $t$ 的多项式），则 $\det(\phi(A) + t\phi(B)) \not\equiv 0$。

    **第一步。** 证明 $\phi$ 保持秩-$(n-1)$ 矩阵的零空间维度不增（即 $\phi$ 不将秩 $n-1$ 矩阵映到更低秩的矩阵）。

    **第二步。** 利用第一步，证明 $\phi$ 是可逆的（$\phi$ 的核是零），因此 $\phi$ 是双射。

    **第三步。** 双射且保持可逆性的线性映射也保持奇异性。由此可证明 $\phi$ 保持秩-1（通过分析 $\det(A + tB) = 0$ 的根的重数），然后由 Marcus-Moyls 定理得到标准形式。$\blacksquare$

!!! example "例 60.6"
    **单位映射保持可逆性的推广。** $\phi(A) = A + f(A) \cdot I$，其中 $f: M_n \to \mathbb{F}$ 是线性泛函。如果 $f$ 足够"小"，$\phi$ 可能保持可逆性。

    例如 $f(A) = c\,\mathrm{tr}(A)$。$\phi(A) = A + c\,\mathrm{tr}(A)\,I$。其行列式

    $$\det(A + c\,\mathrm{tr}(A)\,I) = \det(A) \cdot \det(I + c\,\mathrm{tr}(A)\,A^{-1}),$$

    一般不等于 $\det(A)$。当 $c$ 足够小时，$\phi$ 保持可逆性，但不具有标准形式。

    然而 Dieudonne 定理要求 $\phi$ 是**满射**的。$\phi(A) = A + c\,\mathrm{tr}(A)I$ 是满射的（对 $n \geq 2$，$1 + cn \neq 0$ 时），因此必须具有标准形式——这意味着当 $c \neq 0$ 时，它不保持可逆性。

---

## 60.7 交换性保持

<div class="context-flow" markdown>

**核心问题**：如果线性映射保持矩阵的交换性（$AB = BA$），它必须是什么形式？

</div>

!!! definition "定义 60.8 (交换性保持映射)"
    线性映射 $\phi: M_n(\mathbb{F}) \to M_n(\mathbb{F})$ 称为**交换性保持的**，如果

    $$AB = BA \implies \phi(A)\phi(B) = \phi(B)\phi(A).$$

!!! theorem "定理 60.11 (Watkins 定理, 1976)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$（$n \geq 3$）是**满射**线性映射，保持交换性。则 $\phi$ 具有以下形式之一：

    (a) $\phi(A) = \alpha MAM^{-1} + f(A) I$；

    (b) $\phi(A) = \alpha MA^\top M^{-1} + f(A) I$；

    其中 $\alpha \in \mathbb{C}^*$（非零标量），$M$ 可逆，$f: M_n \to \mathbb{C}$ 是线性泛函。

??? proof "证明"
    **证明要点。**

    **第一步：$\phi$ 将标量矩阵映到标量矩阵。** 标量矩阵 $cI$ 与所有矩阵交换，因此 $\phi(cI)$ 也与所有 $\phi(A)$ 交换。由 $\phi$ 的满射性，$\phi(cI)$ 与所有矩阵交换，故 $\phi(cI) = g(c)I$（$g$ 是某个标量函数）。

    **第二步：模去标量部分。** 定义 $\psi(A) = \phi(A) - f(A)I$，其中 $f(A) = \frac{1}{n}\mathrm{tr}(\phi(A))$。则 $\psi$ 保持交换性（模标量矩阵），且 $\mathrm{tr}(\psi(A)) = 0$。

    **第三步：分析 $\psi$ 在迹零矩阵上的作用。** 利用迹零矩阵空间 $\mathfrak{sl}_n$ 的 Lie 代数结构，$\psi$ 保持 Lie 括号 $[A,B] = AB - BA = 0$ 的条件。这与 $\mathfrak{sl}_n$ 的自同构密切相关。

    **第四步：$\mathfrak{sl}_n$ 的自同构。** 对 $n \geq 3$，简单 Lie 代数 $\mathfrak{sl}_n(\mathbb{C})$ 的自同构只有内自同构 $A \mapsto MAM^{-1}$ 和图自同构 $A \mapsto -A^\top$。因此 $\psi$ 的限制（到 $\mathfrak{sl}_n$）必为 $A \mapsto \alpha MAM^{-1}$ 或 $A \mapsto \alpha MA^\top M^{-1}$。$\blacksquare$

!!! example "例 60.7"
    **$n = 2$ 的特殊情况。** 当 $n = 2$ 时，$\mathfrak{sl}_2(\mathbb{C})$ 是 3 维 Lie 代数，其自同构群更大（$SO(3, \mathbb{C})$），因此交换性保持映射有更多形式。Watkins 定理需要 $n \geq 3$ 的条件。

!!! theorem "定理 60.12 (强交换性保持)"
    如果 $\phi: M_n \to M_n$ 是线性双射，满足更强的条件

    $$\phi(A)\phi(B) = \phi(B)\phi(A) \iff AB = BA,$$

    则 $\phi(A) = \alpha MAM^{-1} + f(A)I$ 或 $\phi(A) = \alpha MA^\top M^{-1} + f(A)I$（$n \geq 3$）。

---

## 60.8 现代方向

<div class="context-flow" markdown>

**核心问题**：线性保持问题有哪些现代推广和开放问题？

</div>

### 60.8.1 数值域保持

!!! definition "定义 60.9 (数值域)"
    矩阵 $A \in M_n(\mathbb{C})$ 的**数值域**（numerical range）为

    $$W(A) = \{x^* A x : x \in \mathbb{C}^n,\, \|x\| = 1\}.$$

    Toeplitz-Hausdorff 定理保证 $W(A)$ 是 $\mathbb{C}$ 的紧凸子集。

!!! theorem "定理 60.13 (数值域保持映射)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 是线性映射，满足 $W(\phi(A)) = W(A)$，$\forall A$。则存在酉矩阵 $U$ 使得

    $$\phi(A) = UAU^* \quad \text{或} \quad \phi(A) = UA^\top U^*.$$

### 60.8.2 算子代数上的保持问题

!!! definition "定义 60.10 (Jordan 同态)"
    线性映射 $\phi: \mathcal{A} \to \mathcal{B}$（$\mathcal{A}, \mathcal{B}$ 为代数）称为 **Jordan 同态**，如果

    $$\phi(A^2) = \phi(A)^2, \quad \forall\, A \in \mathcal{A}.$$

    等价条件：$\phi(AB + BA) = \phi(A)\phi(B) + \phi(B)\phi(A)$。

!!! theorem "定理 60.14 (Jordan 同态的刻画)"
    $M_n(\mathbb{C})$ 上的满射 Jordan 同态 $\phi$ 是以下形式之一：

    (a) **同态**：$\phi(AB) = \phi(A)\phi(B)$，即 $\phi(A) = MAM^{-1}$；

    (b) **反同态**：$\phi(AB) = \phi(B)\phi(A)$，即 $\phi(A) = MA^\top M^{-1}$。

### 60.8.3 非线性保持问题

!!! definition "定义 60.11 (非线性保持问题)"
    现代研究的一个方向是将"线性"条件放松为更弱的条件：

    - **加法保持**：$\phi(A + B) = \phi(A) + \phi(B)$（不要求齐次性）。
    - **乘法保持**：$\phi(AB) = \phi(A)\phi(B)$。
    - **谱加法保持**：$\sigma(A + B) = \sigma(\phi(A) + \phi(B))$。

    在这些更弱的条件下，保持映射的刻画往往更困难但也更有趣。

!!! theorem "定理 60.15 (Kaplansky 猜想的解决)"
    设 $\phi: M_n(\mathbb{C}) \to M_n(\mathbb{C})$ 是满射映射（不一定线性），满足

    $$\sigma(\phi(A) - \phi(B)) = \sigma(A - B), \quad \forall\, A, B.$$

    则 $\phi(A) = MAM^{-1} + C$ 或 $\phi(A) = MA^\top M^{-1} + C$（$M$ 可逆，$C$ 为常矩阵）。

### 60.8.4 量子信息中的保持问题

!!! example "例 60.8"
    在量子信息理论中，**量子信道**是密度矩阵空间上的完全正迹保持（CPTP）映射。Choi 定理（定理 60.8）为量子信道提供了 Kraus 表示。判断一个映射是否为合法的量子信道，本质上是一个线性保持问题：保持正半定性（量子态的合法性）和迹（概率归一化）。

    **纠缠见证**（entanglement witness）利用了正映射与完全正映射的差异：若 $\phi$ 是正映射但非完全正映射，则 $(\phi \otimes \mathrm{id})(\rho) \not\succeq 0$ 可以检测纠缠态 $\rho$。

---

**本章要点总结：**

1. 线性保持问题研究保持矩阵某种性质的线性映射的完整刻画。
2. Frobenius 定理（1897）：行列式保持映射的形式为 $A \mapsto MAN$ 或 $A \mapsto MA^\top N$。
3. Marcus-Moyls 定理：秩-1 保持是许多保持问题的核心中间步骤。
4. 谱保持映射必为相似变换或转置+相似变换。
5. 正映射与完全正映射之间的鸿沟（转置是正映射但非完全正映射）有深刻的量子信息含义。
6. Choi 定理给出完全正映射的 Kraus 算子和表示。
7. 交换性保持、可逆性保持等问题都指向相同的标准形式，体现了矩阵性质的内在刚性。
8. 现代方向包括非线性保持、算子代数保持和量子信息应用。
