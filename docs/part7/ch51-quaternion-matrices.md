# 第 51 章 四元数矩阵与除环上的线性代数

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch2) · 特征值 (Ch6) · SVD (Ch11)

**本章脉络**：四元数回顾 → 除环上的线性代数 → 四元数矩阵 → 左/右特征值 → 四元数行列式 → 四元数 SVD → 复表示 → 旋转表示

**延伸**：四元数在计算机图形学（避免万向节锁）、航天工程（姿态表示）、机器人学（旋转插值 SLERP）中是标准工具；Frobenius 定理指出 $\mathbb{R}, \mathbb{C}, \mathbb{H}$ 是仅有的有限维实结合除代数

</div>

线性代数的经典理论建立在域（field）上——乘法交换且每个非零元素可逆的代数结构。然而，某些重要的应用——特别是三维旋转的表示——自然地引出了一种乘法不交换的"域"，即**四元数** $\mathbb{H}$。四元数是一个**除环**（division ring / skew field）：每个非零元素可逆，但乘法不满足交换律。

在除环上发展线性代数，与域上的情形既有深刻的相似性，也有令人惊讶的差异。例如，四元数矩阵的特征值理论比实矩阵或复矩阵复杂得多——甚至"特征值"这个概念本身就需要仔细定义（左特征值 vs. 右特征值）。而行列式的推广也面临根本困难（交换律的缺失使得 Leibniz 公式无意义）。

尽管如此，四元数矩阵理论不仅在纯数学中有其内在价值，还在航天工程、计算机图形学、机器人学等领域有着不可替代的应用。

---

## 51.1 四元数代数回顾

<div class="context-flow" markdown>

**核心问题**：四元数的基本代数性质是什么？它与实数和复数有什么异同？

</div>

!!! definition "定义 51.1 (四元数)"
    **四元数体** $\mathbb{H}$ 是实数域 $\mathbb{R}$ 上的四维代数，由基 $\{1, i, j, k\}$ 生成，满足

    $$i^2 = j^2 = k^2 = ijk = -1.$$

    等价地，乘法规则为：

    $$ij = k, \quad jk = i, \quad ki = j, \quad ji = -k, \quad kj = -i, \quad ik = -j.$$

    一般四元数写成 $q = a + bi + cj + dk$，其中 $a, b, c, d \in \mathbb{R}$。

    - **实部**：$\operatorname{Re}(q) = a$；
    - **虚部（向量部分）**：$\operatorname{Im}(q) = bi + cj + dk$；
    - **共轭**：$\bar{q} = a - bi - cj - dk$；
    - **模**：$|q| = \sqrt{q\bar{q}} = \sqrt{a^2 + b^2 + c^2 + d^2}$。

!!! theorem "定理 51.1 (四元数的基本性质)"
    1. **除环**：$\mathbb{H}$ 的每个非零元素可逆，$q^{-1} = \frac{\bar{q}}{|q|^2}$。
    2. **非交换**：一般地 $pq \neq qp$。事实上 $pq = qp$ 对所有 $p$ 当且仅当 $q \in \mathbb{R}$（即 $q$ 是实数）。
    3. **共轭的性质**：$\overline{pq} = \bar{q}\bar{p}$（注意顺序反转），$\overline{p+q} = \bar{p} + \bar{q}$。
    4. **模的乘法性**：$|pq| = |p||q|$。
    5. **$\mathbb{H}$ 作为 $\mathbb{R}$-代数**：$\mathbb{H}$ 是 $4$ 维实代数，中心为 $\mathbb{R}$。
    6. **$\mathbb{H}$ 不能被排序**：$\mathbb{H}$ 上不存在与环运算相容的全序。

??? proof "证明（部分）"
    **(1)** $q\bar{q} = (a+bi+cj+dk)(a-bi-cj-dk) = a^2+b^2+c^2+d^2 = |q|^2 > 0$（当 $q \neq 0$）。故 $q^{-1} = \bar{q}/|q|^2$。

    **(3)** 对基元素验证后利用线性性推广。$\overline{ij} = \overline{k} = -k$，$\bar{j}\bar{i} = (-j)(-i) = ji = -k$。

    **(4)** $|pq|^2 = (pq)\overline{(pq)} = pq\bar{q}\bar{p} = p|q|^2\bar{p} = |q|^2 p\bar{p} = |q|^2|p|^2$。（利用 $|q|^2 \in \mathbb{R}$ 可与 $p$ 交换。）

!!! theorem "定理 51.2 (Frobenius 定理)"
    $\mathbb{R}$，$\mathbb{C}$，$\mathbb{H}$ 是仅有的（在同构意义下）有限维实结合除代数。

    特别地，不存在"八元数域"——八元数 $\mathbb{O}$ 是除代数但**不满足结合律**。

!!! example "例 51.1 (四元数的极坐标形式)"
    任何单位四元数 $q$（$|q| = 1$）可以写成

    $$q = \cos\theta + \hat{u}\sin\theta,$$

    其中 $\hat{u}$ 是纯虚单位四元数（$\hat{u}^2 = -1$，$|\hat{u}| = 1$），$\theta \in [0, \pi]$。

    更一般地，任何非零四元数可以写成 $q = |q|(\cos\theta + \hat{u}\sin\theta) = |q|e^{\hat{u}\theta}$。

    **指数映射**：$e^{\hat{u}\theta} = \cos\theta + \hat{u}\sin\theta$（类比 Euler 公式 $e^{i\theta} = \cos\theta + i\sin\theta$）。

!!! definition "定义 51.2 (四元数的共轭类)"
    两个四元数 $p, q$ 称为**相似的**（或**共轭的**），若存在非零 $s \in \mathbb{H}$ 使得 $q = sps^{-1}$。

    **关键事实**：$p, q \in \mathbb{H}$ 相似当且仅当 $\operatorname{Re}(p) = \operatorname{Re}(q)$ 且 $|p| = |q|$。

    特别地，$\mathbb{H}$ 的每个共轭类或者是一个实数 $\{a\}$，或者是以 $(a, 0, 0, 0)$ 为中心、半径 $\sqrt{b^2+c^2+d^2}$ 的二维球面。

---

## 51.2 除环上的线性代数

<div class="context-flow" markdown>

**核心问题**：乘法不交换时，如何定义向量空间和矩阵理论？

</div>

!!! definition "定义 51.3 (除环上的模 / 向量空间)"
    设 $D$ 是除环（如 $\mathbb{H}$）。**右 $D$-模**（或右 $D$-向量空间）$V$ 是一个 Abel 群，配备右标量乘法 $V \times D \to V$，满足

    $$(v_1 + v_2)d = v_1 d + v_2 d, \quad v(d_1 + d_2) = vd_1 + vd_2, \quad v(d_1 d_2) = (vd_1)d_2, \quad v \cdot 1 = v.$$

    类似地定义**左 $D$-模**。注意：由于 $D$ 不交换，左模和右模是不同的概念。

!!! theorem "定理 51.3 (除环上向量空间的维数)"
    设 $D$ 是除环。右 $D$-模 $V$ 的任意两个基的基数相同。这个基数称为 $V$ 在 $D$ 上的**维数**，记为 $\dim_D V$。

??? proof "证明思路"
    虽然 $D$ 不交换，但基的唯一替换性质（Steinitz 交换定理）在除环上仍然成立。关键步骤：若 $v = v_1 d_1 + \cdots + v_n d_n$，$d_1 \neq 0$，则可以解出 $v_1$（利用 $D$ 的可除性），从而进行基的替换。

!!! definition "定义 51.4 (四元数矩阵)"
    $M_{m \times n}(\mathbb{H})$ 是所有 $m \times n$ 四元数矩阵的集合。矩阵乘法以通常方式定义：

    $$(AB)_{ij} = \sum_{k=1}^{p} A_{ik} B_{kj}.$$

    **注意**：由于 $\mathbb{H}$ 不交换，矩阵乘法的性质需要仔细处理。例如：

    - $(AB)^T \neq B^T A^T$ 一般成立；
    - 但可以定义**共轭转置** $A^* = \bar{A}^T$（$(A^*)_{ij} = \overline{A_{ji}}$），则 $(AB)^* = B^* A^*$。

!!! definition "定义 51.5 (四元数线性映射)"
    设 $V, W$ 是右 $\mathbb{H}$-模。映射 $T: V \to W$ 称为**右 $\mathbb{H}$-线性**，若

    $$T(v_1 + v_2) = T(v_1) + T(v_2), \quad T(vq) = T(v)q, \quad \forall v, v_1, v_2 \in V,\, q \in \mathbb{H}.$$

    注意标量在右边。选取 $V$ 和 $W$ 的基后，$T$ 对应于一个四元数矩阵 $A$，**左乘** $v$：$T(v) = Av$（这里 $v$ 是列向量，标量乘法在右，矩阵乘法在左）。

!!! example "例 51.2 (左线性 vs. 右线性)"
    映射 $T: \mathbb{H} \to \mathbb{H}$，$T(q) = iq$ 是**左** $\mathbb{H}$-线性的：$T(qr) = iqr = T(q)r$。但 $T$ 不是右 $\mathbb{H}$-线性的。

    映射 $S: \mathbb{H} \to \mathbb{H}$，$S(q) = qi$ 是**右** $\mathbb{H}$-线性的：$S(q + p) = (q+p)i = qi + pi$，但 $S(qr) = qri \neq S(q)r = qir$（一般）。

    **结论**：对右 $\mathbb{H}$-模，线性映射矩阵从左边作用。

---

## 51.3 四元数矩阵的基本运算

<div class="context-flow" markdown>

**核心问题**：四元数矩阵的特征值如何定义？为什么需要区分左右特征值？

</div>

!!! definition "定义 51.6 (右特征值与左特征值)"
    设 $A \in M_n(\mathbb{H})$。

    - **右特征值**：$\lambda \in \mathbb{H}$ 称为 $A$ 的右特征值，若存在非零 $v \in \mathbb{H}^n$ 使得 $Av = v\lambda$。
    - **左特征值**：$\mu \in \mathbb{H}$ 称为 $A$ 的左特征值，若存在非零 $v \in \mathbb{H}^n$ 使得 $Av = \mu v$。

    **关键区别**：在域上，$Av = \lambda v$ 和 $Av = v\lambda$ 是一回事（因为标量可交换）。但在 $\mathbb{H}$ 上不然！

!!! theorem "定理 51.4 (右特征值的基本性质)"
    设 $A \in M_n(\mathbb{H})$。

    1. 若 $\lambda$ 是右特征值，$Av = v\lambda$，则对任何非零 $s \in \mathbb{H}$，$s^{-1}\lambda s$ 也是右特征值（特征向量为 $vs$）。
    2. 因此右特征值以**共轭类**（相似类）出现。
    3. 每个共轭类恰好包含一个复数（实部非负，虚部在 $\mathbb{R}i$ 中且非负）。这个复数称为**标准右特征值**。
    4. $A$ 恰有 $n$ 个标准右特征值（含重数）。

??? proof "证明"
    **(1)** $Av = v\lambda$。令 $w = vs$（$s \neq 0$），则 $Aw = A(vs) = (Av)s = (v\lambda)s = v(s(s^{-1}\lambda s)) = w(s^{-1}\lambda s)$。

    **(2)** 由 (1) 直接得出。

    **(3)** 设 $\lambda = a + bi + cj + dk$。取 $\hat{u} = \frac{bi+cj+dk}{|bi+cj+dk|}$（若虚部非零），则存在 $s$ 使得 $s^{-1}\hat{u}s = i$（所有纯虚单位四元数共轭），故 $s^{-1}\lambda s = a + |bi+cj+dk| \cdot i \in \mathbb{C}$。

    **(4)** 通过复表示（51.7节）证明。

!!! theorem "定理 51.5 (左特征值的复杂性)"
    左特征值的行为远比右特征值复杂：

    1. $2 \times 2$ 四元数矩阵可以有无穷多个左特征值（甚至不可数个）。
    2. 左特征值不一定存在（某些矩阵没有左特征值）。
    3. 左特征值不具有共轭类封闭性。

!!! example "例 51.3 (左特征值的反例)"
    考虑 $A = \begin{pmatrix} i & 0 \\ 0 & j \end{pmatrix} \in M_2(\mathbb{H})$。

    **右特征值**：$Av = v\lambda$ 要求 $iv_1 = v_1\lambda$，$jv_2 = v_2\lambda$。取 $v_1 \neq 0, v_2 = 0$：$iv_1 = v_1\lambda$，即 $\lambda = v_1^{-1}iv_1$，这是 $i$ 的共轭类中的任意元素。标准右特征值为 $i$。类似地取 $v_2 \neq 0, v_1 = 0$ 得标准右特征值也是 $i$（因为 $j$ 与 $i$ 共轭）。

    **左特征值**：$Av = \mu v$ 要求 $iv_1 = \mu v_1$，$jv_2 = \mu v_2$。若 $v_1, v_2$ 都非零，则 $\mu = iv_1 v_1^{-1} = jv_2 v_2^{-1}$。第一个等式给出 $\mu$ 是 $i$ 的左共轭（不是通常的共轭！），这更加复杂。

!!! definition "定义 51.7 (四元数矩阵的相似)"
    $A, B \in M_n(\mathbb{H})$ 称为**相似的**，若存在可逆 $P \in \operatorname{GL}_n(\mathbb{H})$ 使得 $B = P^{-1}AP$。

    相似矩阵有相同的标准右特征值（含重数）。

!!! example "例 51.4 (四元数矩阵的特征值计算)"
    设 $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \in M_2(\mathbb{H})$。

    求右特征值：$Av = v\lambda$，即

    $$\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} \lambda.$$

    得 $-v_2 = v_1\lambda$，$v_1 = v_2\lambda$。代入：$-v_2 = v_2\lambda^2$，即 $\lambda^2 = -1$。

    在 $\mathbb{H}$ 中，$\lambda^2 = -1$ 的解是所有纯虚单位四元数 $\lambda = ai + bj + ck$（$a^2+b^2+c^2=1$）。这是一个二维球面 $S^2$。

    标准右特征值为 $i$（重数 $2$）。

---

## 51.4 四元数行列式

<div class="context-flow" markdown>

**核心问题**：如何为四元数矩阵定义有意义的行列式？

</div>

四元数的非交换性使得 Leibniz 行列式公式（$\det A = \sum_{\sigma} \operatorname{sgn}(\sigma) \prod a_{i\sigma(i)}$）失去意义——积的顺序不同得到不同的值。历史上出现了几种不同的四元数行列式定义。

!!! definition "定义 51.8 (Study 行列式)"
    设 $A \in M_n(\mathbb{H})$。将 $A$ 视为 $2n \times 2n$ 复矩阵（通过 51.7 节的复表示 $\chi_A$）。**Study 行列式**定义为

    $$\operatorname{det}_S(A) = \sqrt{\det_{\mathbb{C}}(\chi_A)},$$

    其中 $\det_{\mathbb{C}}$ 是通常的复数行列式。$\det_{\mathbb{C}}(\chi_A)$ 总是非负实数，故 $\operatorname{det}_S(A) \in \mathbb{R}_{\geq 0}$。

!!! definition "定义 51.9 (Moore 行列式)"
    对**四元数自伴矩阵**（$A = A^*$，即 $A_{ij} = \overline{A_{ji}}$），Moore 定义了一种递归行列式。但此定义仅适用于自伴矩阵。

!!! definition "定义 51.10 (Dieudonné 行列式)"
    **Dieudonné 行列式**是群同态

    $$\operatorname{det}_D: \operatorname{GL}_n(\mathbb{H}) \to \mathbb{R}_{>0},$$

    由 $\operatorname{det}_D(A) = |Nrd(A)|$ 定义，其中 $Nrd$ 是约化范数。等价地，$\operatorname{det}_D(A) = \operatorname{det}_S(A)^{2/n}$（归一化使得 $\operatorname{det}_D(qI) = |q|^{2n}$）。

    Dieudonné 行列式满足：

    1. $\operatorname{det}_D(AB) = \operatorname{det}_D(A) \operatorname{det}_D(B)$；
    2. $\operatorname{det}_D(A) > 0$ 当且仅当 $A$ 可逆；
    3. $\operatorname{det}_D(I) = 1$。

!!! theorem "定理 51.6 (四元数行列式的困难)"
    不存在满足以下全部条件的映射 $\det: M_n(\mathbb{H}) \to \mathbb{H}$：

    1. $\det(AB) = \det(A) \det(B)$；
    2. $\det$ 是行（或列）的多线性函数；
    3. 交换两行改变 $\det$ 的符号；
    4. $\det(I) = 1$。

    **原因**：条件 2 和 3 要求 $\det$ 对行的置换是交替的，但在非交换环中乘积的顺序不确定，使得条件不相容。

??? proof "证明"
    考虑 $n = 1$。条件 2 要求 $\det(qa) = q\det(a)$（或 $= \det(a)q$）。条件 1 要求 $\det(ab) = \det(a)\det(b)$。但 $\det(ab) = qab$（若 $\det(a) = qa$），而 $\det(a)\det(b) = qa \cdot qb$，两者一般不等。

    更精确地，对 $2 \times 2$ 矩阵，无论乘积顺序如何选取，都无法同时满足多线性和乘法性。

!!! example "例 51.5 (Study 行列式的计算)"
    设 $A = \begin{pmatrix} 1+i & j \\ -j & 1-i \end{pmatrix} \in M_2(\mathbb{H})$。

    复表示（每个四元数 $a+bi+cj+dk \mapsto \begin{pmatrix} a+bi & c+di \\ -c+di & a-bi \end{pmatrix}$）：

    $$\chi_A = \begin{pmatrix} 1+i & 0 & 0 & 1 \\ 0 & 1-i & -1 & 0 \\ 0 & -1 & 1-i & 0 \\ 1 & 0 & 0 & 1+i \end{pmatrix}.$$

    $\det_{\mathbb{C}}(\chi_A) = \ldots = 4$（经计算）。故 $\operatorname{det}_S(A) = \sqrt{4} = 2$。

---

## 51.5 四元数矩阵的谱定理

<div class="context-flow" markdown>

**核心问题**：四元数自伴矩阵是否可以对角化？特征值是什么？

</div>

!!! definition "定义 51.11 (四元数 Hermite 矩阵)"
    $A \in M_n(\mathbb{H})$ 称为**四元数 Hermite 矩阵**（或自伴矩阵），若 $A = A^*$，即 $A_{ij} = \overline{A_{ji}}$。

    四元数 Hermite 矩阵的对角元素必为实数（因为 $A_{ii} = \overline{A_{ii}}$ 意味着 $A_{ii} \in \mathbb{R}$）。

!!! theorem "定理 51.7 (四元数谱定理)"
    设 $A \in M_n(\mathbb{H})$ 是四元数 Hermite 矩阵。则：

    1. $A$ 的所有右特征值都是**实数**。
    2. 存在**酉矩阵** $U \in M_n(\mathbb{H})$（$U^* U = I$）使得 $U^* A U = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$，$\lambda_i \in \mathbb{R}$。
    3. 不同特征值对应的特征向量正交（关于四元数内积 $\langle u, v \rangle = u^* v = \sum_i \bar{u}_i v_i$）。

??? proof "证明"
    **(1)** 设 $Av = v\lambda$，$v \neq 0$。计算

    $$v^* A v = v^* (v\lambda) = (v^* v)\lambda = \|v\|^2 \lambda.$$

    另一方面，$(v^* A v)^* = v^* A^* v = v^* A v$（因为 $A = A^*$，且 $v^* Av$ 是 $1 \times 1$ 矩阵）。因此 $v^* Av \in \mathbb{R}$，故 $\lambda = \frac{v^* Av}{\|v\|^2} \in \mathbb{R}$。

    **(2)** 存在性的证明与实/复情形类似，通过归纳法。设 $\lambda_1$ 是 $A$ 的一个特征值（实数），$v_1$ 是对应的单位特征向量。将 $v_1$ 扩展为 $\mathbb{H}^n$ 的正交归一基，得到酉矩阵 $U_1 = [v_1, v_2, \ldots, v_n]$。

    $$U_1^* A U_1 = \begin{pmatrix} \lambda_1 & w^* \\ w & A' \end{pmatrix}.$$

    由 Hermite 性，$w = 0$（与复情形的论证相同）。对 $A'$（$(n-1) \times (n-1)$ 四元数 Hermite 矩阵）递归。

    **(3)** 设 $Av = v\lambda$，$Aw = w\mu$，$\lambda \neq \mu$。

    $$\lambda \langle v, w \rangle = \lambda v^* w = (\bar{\lambda}v)^* w = (Av)^* w = v^* A^* w = v^* Aw = v^* (w\mu) = (v^*w)\mu = \langle v, w \rangle \mu.$$

    由于 $\lambda, \mu \in \mathbb{R}$，$\lambda \langle v, w \rangle = \langle v, w \rangle \mu$，即 $(\lambda - \mu)\langle v, w \rangle = 0$。$\lambda \neq \mu$ 故 $\langle v, w \rangle = 0$。

!!! example "例 51.6 (四元数 Hermite 矩阵的对角化)"
    $$A = \begin{pmatrix} 2 & j \\ -j & 3 \end{pmatrix}.$$

    验证 $A = A^*$：$A_{12} = j$，$\overline{A_{21}} = \overline{-j} = j$。

    特征方程（由复表示或直接计算）：$(2-\lambda)(3-\lambda) - j(-j) = (2-\lambda)(3-\lambda) - 1 = 0$，即 $\lambda^2 - 5\lambda + 5 = 0$。解为 $\lambda = \frac{5 \pm \sqrt{5}}{2}$（实数）。

---

## 51.6 四元数 SVD

<div class="context-flow" markdown>

**核心问题**：四元数矩阵是否有奇异值分解？奇异值仍然是实数吗？

</div>

!!! theorem "定理 51.8 (四元数 SVD)"
    设 $A \in M_{m \times n}(\mathbb{H})$。则存在四元数酉矩阵 $U \in M_m(\mathbb{H})$（$U^*U = I_m$）和 $V \in M_n(\mathbb{H})$（$V^*V = I_n$），以及非负实对角矩阵 $\Sigma$ 使得

    $$A = U \Sigma V^*.$$

    对角元素 $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_{\min(m,n)} \geq 0$ 是**实的非负奇异值**。

??? proof "证明"
    考虑四元数 Hermite 矩阵 $A^* A \in M_n(\mathbb{H})$。由谱定理（定理 51.7），存在酉 $V$ 使得 $V^*(A^*A)V = \operatorname{diag}(\sigma_1^2, \ldots, \sigma_n^2)$，$\sigma_i \geq 0$。

    对非零 $\sigma_i$，定义 $u_i = Av_i / \sigma_i$，则 $\{u_i\}$ 是正交归一的（$u_i^* u_j = v_i^* A^* A v_j / (\sigma_i \sigma_j) = \sigma_i^2 \delta_{ij} / (\sigma_i \sigma_j) = \delta_{ij}$）。将 $\{u_i\}$ 扩展为 $\mathbb{H}^m$ 的正交归一基得 $U$。

    则 $A = U\Sigma V^*$。

!!! example "例 51.7 (四元数 SVD 的计算)"
    设 $A = \begin{pmatrix} i \\ j \end{pmatrix} \in M_{2 \times 1}(\mathbb{H})$。

    $A^* A = \begin{pmatrix} -i & -j \end{pmatrix} \begin{pmatrix} i \\ j \end{pmatrix} = (-i)(i) + (-j)(j) = 1 + 1 = 2$。

    故奇异值 $\sigma_1 = \sqrt{2}$，$V = (1)$（$1 \times 1$ 酉矩阵）。

    $u_1 = A / \sqrt{2} = \frac{1}{\sqrt{2}}\begin{pmatrix} i \\ j \end{pmatrix}$。扩展为正交归一基：$u_2 = \frac{1}{\sqrt{2}}\begin{pmatrix} j \\ -i \end{pmatrix}$（验证 $u_1^* u_2 = \frac{1}{2}((-i)(j) + (-j)(-i)) = \frac{1}{2}(-k + k) = 0$）。

    $$A = \frac{1}{\sqrt{2}}\begin{pmatrix} i & j \\ j & -i \end{pmatrix} \begin{pmatrix} \sqrt{2} \\ 0 \end{pmatrix} (1).$$

---

## 51.7 复表示

<div class="context-flow" markdown>

**核心问题**：如何将四元数矩阵"翻译"为复矩阵？这保持了哪些信息？

</div>

!!! definition "定义 51.12 (四元数到复矩阵的嵌入)"
    每个四元数 $q = a + bi + cj + dk$ 可以表示为 $q = z_1 + z_2 j$，其中 $z_1 = a + bi$，$z_2 = c + di \in \mathbb{C}$（注意 $j z = \bar{z} j$ 对 $z \in \mathbb{C}$）。

    定义**复表示映射** $\chi: \mathbb{H} \to M_2(\mathbb{C})$：

    $$\chi(z_1 + z_2 j) = \begin{pmatrix} z_1 & z_2 \\ -\bar{z}_2 & \bar{z}_1 \end{pmatrix}.$$

    推广到矩阵：$\chi: M_n(\mathbb{H}) \to M_{2n}(\mathbb{C})$，将 $n \times n$ 四元数矩阵 $A = (a_{ij})$ 映为 $2n \times 2n$ 复矩阵，通过对每个元素 $a_{ij}$ 替换为其 $2 \times 2$ 复表示。

!!! theorem "定理 51.9 (复表示的性质)"
    $\chi: M_n(\mathbb{H}) \to M_{2n}(\mathbb{C})$ 满足：

    1. **加法性**：$\chi(A + B) = \chi(A) + \chi(B)$；
    2. **乘法性**：$\chi(AB) = \chi(A)\chi(B)$；
    3. **共轭**：$\chi(\bar{A}^T) = \chi(A)^*$（Hermite 共轭）；
    4. **保秩**：$\operatorname{rank}_{\mathbb{H}}(A) = \frac{1}{2}\operatorname{rank}_{\mathbb{C}}(\chi(A))$；
    5. **特征值**：$A$ 的标准右特征值 $\lambda_1, \ldots, \lambda_n \in \mathbb{C}$ 恰好是 $\chi(A)$ 的特征值（每个出现两次，与 $\bar{\lambda}_i$ 配对）；
    6. **$\chi$ 的像的刻画**：$C \in M_{2n}(\mathbb{C})$ 在 $\chi$ 的像中当且仅当 $J_n \bar{C} J_n^{-1} = C$，其中 $J_n = \operatorname{diag}(J, \ldots, J)$，$J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$。

??? proof "证明（部分）"
    **(2) 乘法性：** 验证 $\chi(q_1 q_2) = \chi(q_1)\chi(q_2)$ 对 $1 \times 1$ 情形成立后，矩阵情形自动成立。

    设 $q_1 = z_1 + z_2 j$，$q_2 = w_1 + w_2 j$。则

    $$q_1 q_2 = (z_1 w_1 - z_2 \bar{w}_2) + (z_1 w_2 + z_2 \bar{w}_1) j$$

    （利用 $jw = \bar{w}j$ 对 $w \in \mathbb{C}$）。另一方面，

    $$\chi(q_1)\chi(q_2) = \begin{pmatrix} z_1 & z_2 \\ -\bar{z}_2 & \bar{z}_1 \end{pmatrix} \begin{pmatrix} w_1 & w_2 \\ -\bar{w}_2 & \bar{w}_1 \end{pmatrix} = \begin{pmatrix} z_1 w_1 - z_2\bar{w}_2 & z_1 w_2 + z_2\bar{w}_1 \\ -\bar{z}_2 w_1 - \bar{z}_1\bar{w}_2 & -\bar{z}_2 w_2 + \bar{z}_1\bar{w}_1 \end{pmatrix}.$$

    与 $\chi(q_1 q_2)$ 比较一致。

    **(5)** $Av = v\lambda$ 在复表示下变为 $\chi(A)\chi(v) = \chi(v)\chi(\lambda)$。$\chi(\lambda) = \begin{pmatrix} \lambda & 0 \\ 0 & \bar{\lambda} \end{pmatrix}$（当 $\lambda \in \mathbb{C}$）。因此 $\chi(A)$ 的特征值包含 $\lambda$ 和 $\bar{\lambda}$。

!!! example "例 51.8 (复表示的计算)"
    $A = \begin{pmatrix} i & j \\ k & 1 \end{pmatrix} \in M_2(\mathbb{H})$。

    $i = i + 0 \cdot j \mapsto \begin{pmatrix} i & 0 \\ 0 & -i \end{pmatrix}$，$j = 0 + 1 \cdot j \mapsto \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$，

    $k = 0 + i \cdot j \mapsto \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}$，$1 \mapsto \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$。

    $$\chi(A) = \begin{pmatrix} i & 0 & 0 & 1 \\ 0 & -i & -1 & 0 \\ 0 & i & 1 & 0 \\ i & 0 & 0 & 1 \end{pmatrix}.$$

---

## 51.8 旋转表示

<div class="context-flow" markdown>

**核心问题**：四元数如何表示三维旋转？与其他表示方法相比有何优势？

</div>

!!! theorem "定理 51.10 (单位四元数与三维旋转)"
    设 $q \in \mathbb{H}$，$|q| = 1$。映射

    $$R_q: \mathbb{R}^3 \to \mathbb{R}^3, \quad R_q(v) = qv\bar{q}$$

    （将纯虚四元数 $v = xi + yj + zk$ 视为 $\mathbb{R}^3$ 的向量）是一个三维旋转。

    映射 $\Phi: S^3 \to SO(3)$，$\Phi(q) = R_q$，是**满的群同态**，其核为 $\{1, -1\}$。因此

    $$SO(3) \cong S^3 / \{\pm 1\} \cong SU(2) / \{\pm I\}.$$

??? proof "证明"
    **$R_q$ 保持纯虚性：** 若 $v$ 是纯虚的（$\bar{v} = -v$），则 $\overline{qv\bar{q}} = q\bar{v}\bar{q} = -qv\bar{q}$，故 $R_q(v)$ 也是纯虚的。

    **$R_q$ 是正交变换：** $|R_q(v)| = |qv\bar{q}| = |q||v||\bar{q}| = |v|$。

    **$R_q$ 是旋转（$\det R_q = +1$）：** $\Phi$ 连续，$\Phi(1) = I$（$\det = 1$），$S^3$ 连通，故 $\det R_q = +1$ 恒成立。

    **满射性：** 设 $q = \cos\frac{\theta}{2} + \hat{u}\sin\frac{\theta}{2}$（$\hat{u}$ 为纯虚单位四元数）。可以验证 $R_q$ 是绕 $\hat{u}$ 轴旋转角度 $\theta$。由于任意旋转都可以这样表示（轴角表示），$\Phi$ 满射。

    **核：** $R_q = I \Leftrightarrow qv = vq$ 对所有纯虚 $v \Leftrightarrow q \in \mathbb{R} \Leftrightarrow q = \pm 1$（因为 $|q| = 1$）。

!!! definition "定义 51.13 (SLERP 插值)"
    **球面线性插值**（SLERP, Spherical Linear intERPolation）：给定两个单位四元数 $q_0, q_1$（代表两个旋转），参数 $t \in [0, 1]$，SLERP 定义为

    $$\operatorname{SLERP}(q_0, q_1, t) = q_0 (q_0^{-1} q_1)^t = q_0 \frac{\sin((1-t)\Omega)}{\sin\Omega} + q_1 \frac{\sin(t\Omega)}{\sin\Omega},$$

    其中 $\Omega = \arccos(\operatorname{Re}(q_0^{-1} q_1))$ 是 $q_0, q_1$ 之间的角度。

    SLERP 在单位四元数球面 $S^3$ 上给出**大圆弧**插值，对应于恒定角速度的旋转插值。这是计算机图形学和动画中的标准旋转插值方法。

!!! theorem "定理 51.11 (四元数旋转 vs. 旋转矩阵 vs. 欧拉角)"
    | 表示方法 | 参数数 | 万向节锁 | 插值 | 复合 | 归一化 |
    |:---|:---:|:---:|:---:|:---:|:---:|
    | 旋转矩阵 $R \in SO(3)$ | 9（6 约束） | 否 | 困难 | 矩阵乘法 | 需 Gram-Schmidt |
    | 欧拉角 $(\alpha, \beta, \gamma)$ | 3 | **有** | 线性（非均匀） | 困难 | 自动 |
    | 轴角 $(\hat{n}, \theta)$ | 4（1 约束） | 否 | 困难 | 困难 | 归一化轴 |
    | 单位四元数 $q$ | 4（1 约束） | **否** | SLERP（最优） | 四元数乘法 | 归一化 |

    四元数表示的优势在于：

    1. **无万向节锁**：不存在欧拉角在 $\beta = \pm 90°$ 时丧失一个自由度的问题；
    2. **高效插值**：SLERP 提供恒角速度插值；
    3. **紧凑表示**：仅 4 个参数（vs. 矩阵的 9 个）；
    4. **数值稳定**：归一化只需除以模长（vs. 矩阵的 Gram-Schmidt）。

!!! example "例 51.9 (四元数到旋转矩阵的转换)"
    设 $q = a + bi + cj + dk$（$|q| = 1$）。对应的旋转矩阵为

    $$R_q = \begin{pmatrix} 1-2(c^2+d^2) & 2(bc-ad) & 2(bd+ac) \\ 2(bc+ad) & 1-2(b^2+d^2) & 2(cd-ab) \\ 2(bd-ac) & 2(cd+ab) & 1-2(b^2+c^2) \end{pmatrix}.$$

    **例**：$q = \frac{1}{\sqrt{2}}(1 + k) = \frac{1}{\sqrt{2}} + \frac{1}{\sqrt{2}}k$（绕 $z$ 轴旋转 $90°$）。$a = d = 1/\sqrt{2}$，$b = c = 0$。

    $$R_q = \begin{pmatrix} 1-2 \cdot 1/2 & 0-2 \cdot 1/2 & 0 \\ 0+2 \cdot 1/2 & 1-2 \cdot 1/2 & 0 \\ 0 & 0 & 1 \end{pmatrix} = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 1 \end{pmatrix}.$$

    这确实是绕 $z$ 轴逆时针旋转 $90°$ 的矩阵。

!!! example "例 51.10 (SLERP 的应用)"
    设初始姿态 $q_0 = 1$（无旋转），目标姿态 $q_1 = \frac{1}{\sqrt{2}}(1+i)$（绕 $x$ 轴旋转 $90°$）。

    $q_0^{-1} q_1 = q_1$，$\Omega = \arccos(\operatorname{Re}(q_1)) = \arccos(1/\sqrt{2}) = \pi/4$。

    在 $t = 0.5$ 处：

    $$\operatorname{SLERP}(q_0, q_1, 0.5) = q_0 \frac{\sin(\pi/8)}{\sin(\pi/4)} + q_1 \frac{\sin(\pi/8)}{\sin(\pi/4)} = \frac{\sin(\pi/8)}{\sin(\pi/4)}(1 + q_1).$$

    这对应绕 $x$ 轴旋转 $45°$，正好是 $0°$ 和 $90°$ 的中间旋转。

---

### 本章总结

四元数矩阵理论展示了当标量环从域变为除环时，线性代数的哪些定理需要修改：

| 概念 | 域上 | 四元数上 |
|:---|:---|:---|
| 特征值 | 唯一定义 | 左/右之分；右特征值以共轭类出现 |
| 行列式 | Leibniz 公式 | 无自然定义；Study/Dieudonné 行列式 |
| 谱定理 | Hermite $\Rightarrow$ 实特征值 | 四元数 Hermite $\Rightarrow$ 实特征值（类似） |
| SVD | 存在 | 存在，奇异值为实数 |
| 维数 | 良定义 | 良定义（Steinitz 定理推广） |

---

### 习题

!!! exercise "习题 51.1"
    验证 $q = \frac{1}{2}(1+i+j+k)$ 是单位四元数，计算 $q^{-1}$，并求 $R_q$ 对应的旋转轴和角度。

!!! exercise "习题 51.2"
    求 $A = \begin{pmatrix} 0 & i \\ i & 0 \end{pmatrix}$ 的所有标准右特征值。

!!! exercise "习题 51.3"
    设 $A = \begin{pmatrix} 1 & i \\ -i & 2 \end{pmatrix}$。验证 $A = A^*$，求其特征值并将其对角化。

!!! exercise "习题 51.4"
    计算 $A = \begin{pmatrix} j \end{pmatrix}$（$1 \times 1$ 矩阵）的复表示 $\chi(A)$，并验证 $\det_{\mathbb{C}}(\chi(A)) = |j|^2 = 1$。

!!! exercise "习题 51.5"
    证明 SLERP 满足 $\operatorname{SLERP}(q_0, q_1, 0) = q_0$，$\operatorname{SLERP}(q_0, q_1, 1) = q_1$。

!!! exercise "习题 51.6"
    给出一个 $2 \times 2$ 四元数矩阵的例子，使其具有无穷多个左特征值。
