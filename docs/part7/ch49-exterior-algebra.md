# 第 49 章 外代数与 Grassmannian

<div class="context-flow" markdown>

**前置**：向量空间 (Ch4) · 线性变换 (Ch5) · 多线性代数 (Ch21) · 行列式 (Ch3)

**本章脉络**：外积（楔积）→ 外幂空间 $\Lambda^k(V)$ → 交替多线性映射 → 行列式作为顶外幂 → 复合矩阵 → Grassmannian → Plücker 坐标 → 外代数的万有性质

**延伸**：外代数是微分形式（现代微分几何和物理学的语言）的代数基础；Grassmannian 是代数几何和优化（子空间跟踪）中的核心对象；复合矩阵给出了特征值乘积的不等式

</div>

向量空间中的标量积 $\langle u, v \rangle$ 衡量两个向量的"对齐程度"。但很多几何对象——面积、体积、通量——本质上需要一种**反对称**的乘法运算。两个向量 $u, v$ 张成的平行四边形的（有向）面积，不能用点积或矩阵乘法自然表达，却恰好由**楔积**（wedge product）$u \wedge v$ 完美捕捉。

外代数（exterior algebra）系统地构建了这种反对称乘法。它不仅为行列式提供了最自然的定义，还引出了 Grassmannian——参数化所有 $k$ 维子空间的几何对象。从微分形式到代数几何中的 Schubert 演算，外代数的影响无处不在。

---

## 49.1 外积的定义

<div class="context-flow" markdown>

**核心问题**：如何定义一种反对称的向量乘法？它与行列式有什么关系？

</div>

!!! definition "定义 49.1 (外积 / 楔积)"
    设 $V$ 是域 $\mathbb{F}$ 上的向量空间。两个向量 $u, v \in V$ 的**外积**（或**楔积**，wedge product）$u \wedge v$ 是一个满足以下性质的抽象对象：

    1. **双线性**：$(au_1 + bu_2) \wedge v = a(u_1 \wedge v) + b(u_2 \wedge v)$，$u \wedge (av_1 + bv_2) = a(u \wedge v_1) + b(u \wedge v_2)$；
    2. **反对称性**：$u \wedge u = 0$，对所有 $u \in V$。

    由反对称性立即推出**交替性**（当 $\operatorname{char}(\mathbb{F}) \neq 2$ 时两者等价）：

    $$u \wedge v = -(v \wedge u).$$

    这是因为 $0 = (u+v) \wedge (u+v) = u \wedge u + u \wedge v + v \wedge u + v \wedge v = u \wedge v + v \wedge u$。

!!! example "例 49.1 (平面中的楔积)"
    设 $V = \mathbb{R}^2$，$e_1 = (1,0)$，$e_2 = (0,1)$。对 $u = a_1 e_1 + a_2 e_2$，$v = b_1 e_1 + b_2 e_2$：

    $$u \wedge v = (a_1 e_1 + a_2 e_2) \wedge (b_1 e_1 + b_2 e_2) = a_1 b_2 (e_1 \wedge e_2) + a_2 b_1 (e_2 \wedge e_1)$$

    $$= (a_1 b_2 - a_2 b_1)(e_1 \wedge e_2).$$

    系数 $a_1 b_2 - a_2 b_1 = \det \begin{pmatrix} a_1 & b_1 \\ a_2 & b_2 \end{pmatrix}$ 正好是 $u, v$ 张成的平行四边形的有向面积。

!!! example "例 49.2 (楔积的几何意义)"
    在 $\mathbb{R}^3$ 中，设 $u = (1, 0, 0)$，$v = (0, 1, 0)$，$w = (0, 0, 1)$。

    $$u \wedge v = e_1 \wedge e_2, \quad u \wedge w = e_1 \wedge e_3, \quad v \wedge w = e_2 \wedge e_3.$$

    $$u \wedge v \wedge w = e_1 \wedge e_2 \wedge e_3,$$

    这个三重楔积代表 $u, v, w$ 围成的平行六面体的有向体积。

    对一般的 $u = (u_1, u_2, u_3)$，$v = (v_1, v_2, v_3)$：

    $$u \wedge v = (u_1 v_2 - u_2 v_1)(e_1 \wedge e_2) + (u_1 v_3 - u_3 v_1)(e_1 \wedge e_3) + (u_2 v_3 - u_3 v_2)(e_2 \wedge e_3).$$

    三个系数正好是叉积 $u \times v$ 的分量（带适当符号）。事实上，在 $\mathbb{R}^3$ 中，二重楔积与叉积通过 Hodge 对偶联系：$\star(u \wedge v) = u \times v$。

---

## 49.2 外幂空间

<div class="context-flow" markdown>

**核心问题**：$k$ 个向量的楔积生成什么空间？这个空间的维数是多少？

</div>

!!! definition "定义 49.2 ($k$ 次外幂空间)"
    设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间。**$k$ 次外幂空间** $\Lambda^k(V)$ 是由所有 $k$ 重楔积

    $$v_1 \wedge v_2 \wedge \cdots \wedge v_k, \quad v_1, \ldots, v_k \in V$$

    的 $\mathbb{F}$-线性组合构成的向量空间。$\Lambda^k(V)$ 的元素称为 **$k$-向量**（$k$-vector）或 **$k$-形式**。

    特别约定：$\Lambda^0(V) = \mathbb{F}$，$\Lambda^1(V) = V$。

!!! theorem "定理 49.1 (外幂空间的基与维数)"
    设 $\{e_1, \ldots, e_n\}$ 是 $V$ 的基。则 $\Lambda^k(V)$ 有基

    $$\{e_{i_1} \wedge e_{i_2} \wedge \cdots \wedge e_{i_k} : 1 \leq i_1 < i_2 < \cdots < i_k \leq n\},$$

    因此

    $$\dim \Lambda^k(V) = \binom{n}{k}.$$

    特别地，$\dim \Lambda^n(V) = 1$（**顶外幂**），$\Lambda^k(V) = 0$（$k > n$ 时）。

??? proof "证明"
    **生成性：** 任意 $v_1 \wedge \cdots \wedge v_k$ 中，将每个 $v_j$ 展开为 $e_i$ 的线性组合，利用多线性和反对称性，可将其表示为 $e_{i_1} \wedge \cdots \wedge e_{i_k}$（$i_1 < \cdots < i_k$）的线性组合。

    具体地，若 $v_j = \sum_i a_{ij} e_i$，则

    $$v_1 \wedge \cdots \wedge v_k = \sum_{i_1, \ldots, i_k} a_{i_1 1} a_{i_2 2} \cdots a_{i_k k} \, e_{i_1} \wedge \cdots \wedge e_{i_k}.$$

    当某两个指标相同时，$e_{i_1} \wedge \cdots \wedge e_{i_k} = 0$（反对称性）。对不同的指标，可通过交换排序为 $i_1 < \cdots < i_k$，每次交换引入一个负号。

    **线性无关性：** 通过构造对偶基（交替多线性函数）来证明。对每个递增指标组 $I = \{i_1, \ldots, i_k\}$，定义 $\varphi_I: V^k \to \mathbb{F}$ 为

    $$\varphi_I(v_1, \ldots, v_k) = \det(v_j \text{ 的第 } i_l \text{ 分量})_{l,j}.$$

    $\varphi_I$ 是交替 $k$-线性函数。可以验证 $\varphi_I(e_{j_1}, \ldots, e_{j_k}) = \delta_{I,J}$（$J = \{j_1, \ldots, j_k\}$），这证明了线性无关性。

!!! definition "定义 49.3 (可分解与不可分解 $k$-向量)"
    $\Lambda^k(V)$ 中的元素 $\omega$ 称为**可分解的**（decomposable），若 $\omega = v_1 \wedge \cdots \wedge v_k$ 对某些 $v_i \in V$。否则称为**不可分解的**。

    当 $k = 1$ 或 $k = n-1$ 时，每个非零 $k$-向量都是可分解的。但对 $2 \leq k \leq n-2$，一般的 $k$-向量不可分解。

!!! example "例 49.3 (不可分解的 2-向量)"
    在 $\mathbb{R}^4$ 中取标准基 $\{e_1, e_2, e_3, e_4\}$。考虑

    $$\omega = e_1 \wedge e_2 + e_3 \wedge e_4 \in \Lambda^2(\mathbb{R}^4).$$

    **断言：** $\omega$ 不可分解。

    **证明：** 若 $\omega = u \wedge v$，其中 $u = \sum a_i e_i$，$v = \sum b_i e_i$，则

    $$\omega \wedge \omega = (u \wedge v) \wedge (u \wedge v) = 0$$

    （因为 $u \wedge v \wedge u \wedge v$ 中 $u$ 出现两次）。但

    $$\omega \wedge \omega = (e_1 \wedge e_2 + e_3 \wedge e_4) \wedge (e_1 \wedge e_2 + e_3 \wedge e_4) = 2 \, e_1 \wedge e_2 \wedge e_3 \wedge e_4 \neq 0.$$

    矛盾，故 $\omega$ 不可分解。

    这个判据可以推广：$\omega \in \Lambda^2(V)$ 可分解当且仅当 $\omega \wedge \omega = 0$。

---

## 49.3 外代数

<div class="context-flow" markdown>

**核心问题**：如何将所有外幂空间组合成一个带有乘法结构的代数？

</div>

!!! definition "定义 49.4 (外代数)"
    $V$ 的**外代数**（exterior algebra）定义为所有外幂空间的直和：

    $$\Lambda(V) = \bigoplus_{k=0}^{n} \Lambda^k(V) = \Lambda^0(V) \oplus \Lambda^1(V) \oplus \cdots \oplus \Lambda^n(V),$$

    其维数为 $\sum_{k=0}^n \binom{n}{k} = 2^n$。

    $\Lambda(V)$ 配备楔积乘法 $\wedge: \Lambda^p(V) \times \Lambda^q(V) \to \Lambda^{p+q}(V)$，使之成为一个**分次代数**（graded algebra）。

!!! theorem "定理 49.2 (外代数的乘法规则)"
    设 $\alpha \in \Lambda^p(V)$，$\beta \in \Lambda^q(V)$。则：

    1. **分次反交换性**（graded commutativity）：$\alpha \wedge \beta = (-1)^{pq} \beta \wedge \alpha$；
    2. **结合律**：$(\alpha \wedge \beta) \wedge \gamma = \alpha \wedge (\beta \wedge \gamma)$；
    3. **单位元**：$1 \in \Lambda^0(V) = \mathbb{F}$ 满足 $1 \wedge \alpha = \alpha \wedge 1 = \alpha$。

    特别地，当 $p$ 和 $q$ 都是奇数时，$\alpha \wedge \beta = -\beta \wedge \alpha$。

??? proof "证明"
    对可分解元素 $\alpha = u_1 \wedge \cdots \wedge u_p$，$\beta = v_1 \wedge \cdots \wedge v_q$，有

    $$\beta \wedge \alpha = v_1 \wedge \cdots \wedge v_q \wedge u_1 \wedge \cdots \wedge u_p.$$

    要将此排列回 $\alpha \wedge \beta = u_1 \wedge \cdots \wedge u_p \wedge v_1 \wedge \cdots \wedge v_q$ 的顺序，需将 $q$ 个 $v_i$ 依次穿过 $p$ 个 $u_j$，共需 $pq$ 次相邻交换，每次交换引入因子 $(-1)$，故 $\beta \wedge \alpha = (-1)^{pq} \alpha \wedge \beta$。

    结合律和单位元的证明直接由定义给出。

!!! definition "定义 49.5 (外代数的张量构造)"
    更严格地，外代数可通过张量代数的商来定义。设 $T(V) = \bigoplus_{k=0}^{\infty} V^{\otimes k}$ 是 $V$ 的张量代数，$\mathcal{I}$ 是由所有 $v \otimes v$（$v \in V$）生成的双边理想。则

    $$\Lambda(V) = T(V) / \mathcal{I}.$$

    楔积 $v_1 \wedge \cdots \wedge v_k$ 就是 $v_1 \otimes \cdots \otimes v_k$ 在商代数中的像。

!!! example "例 49.4 (低维外代数)"
    **(a) $\dim V = 2$：** $\Lambda(\mathbb{R}^2)$ 的维数为 $2^2 = 4$。基为

    $$\{1, \, e_1, \, e_2, \, e_1 \wedge e_2\}.$$

    $\Lambda^0 = \mathbb{R}$（标量），$\Lambda^1 = \mathbb{R}^2$（向量），$\Lambda^2 = \mathbb{R}$（面积元素）。

    **(b) $\dim V = 3$：** $\Lambda(\mathbb{R}^3)$ 的维数为 $2^3 = 8$。基为

    $$\{1, \, e_1, e_2, e_3, \, e_1 \wedge e_2, e_1 \wedge e_3, e_2 \wedge e_3, \, e_1 \wedge e_2 \wedge e_3\}.$$

    维数分布：$1, 3, 3, 1$（对称的！这是因为 $\binom{3}{k} = \binom{3}{3-k}$）。

    **(c) $\dim V = 4$：** $\Lambda(\mathbb{R}^4)$ 的维数为 $2^4 = 16$。维数分布：$1, 4, 6, 4, 1$（Pascal 三角第 4 行）。

---

## 49.4 行列式与顶外幂

<div class="context-flow" markdown>

**核心问题**：行列式如何自然地从外代数中产生？这给出了行列式性质的哪些优雅证明？

</div>

!!! theorem "定理 49.3 (行列式的外代数定义)"
    设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$T: V \to V$ 是线性变换。$T$ 自然诱导映射 $\Lambda^n(T): \Lambda^n(V) \to \Lambda^n(V)$：

    $$\Lambda^n(T)(v_1 \wedge \cdots \wedge v_n) = T(v_1) \wedge T(v_2) \wedge \cdots \wedge T(v_n).$$

    由于 $\dim \Lambda^n(V) = 1$，$\Lambda^n(T)$ 是标量乘法。这个标量就是 $\det(T)$：

    $$T(v_1) \wedge T(v_2) \wedge \cdots \wedge T(v_n) = \det(T) \cdot (v_1 \wedge v_2 \wedge \cdots \wedge v_n).$$

??? proof "证明"
    选取 $V$ 的基 $\{e_1, \ldots, e_n\}$，$\Lambda^n(V)$ 的基为 $\{e_1 \wedge \cdots \wedge e_n\}$。设 $T(e_j) = \sum_i a_{ij} e_i$。则

    $$T(e_1) \wedge \cdots \wedge T(e_n) = \left(\sum_{i_1} a_{i_1 1} e_{i_1}\right) \wedge \cdots \wedge \left(\sum_{i_n} a_{i_n n} e_{i_n}\right)$$

    $$= \sum_{i_1, \ldots, i_n} a_{i_1 1} a_{i_2 2} \cdots a_{i_n n} \, e_{i_1} \wedge \cdots \wedge e_{i_n}.$$

    只有当 $(i_1, \ldots, i_n)$ 是 $(1, \ldots, n)$ 的置换时非零：

    $$= \sum_{\sigma \in S_n} a_{\sigma(1), 1} a_{\sigma(2), 2} \cdots a_{\sigma(n), n} \, e_{\sigma(1)} \wedge \cdots \wedge e_{\sigma(n)}$$

    $$= \left(\sum_{\sigma \in S_n} \operatorname{sgn}(\sigma) \, a_{\sigma(1), 1} \cdots a_{\sigma(n), n}\right) \, e_1 \wedge \cdots \wedge e_n$$

    $$= \det(A) \cdot e_1 \wedge \cdots \wedge e_n.$$

!!! theorem "定理 49.4 (行列式性质的外代数证明)"
    外代数定义使行列式的基本性质变得显然：

    1. **$\det(ST) = \det(S)\det(T)$：** $\Lambda^n(ST) = \Lambda^n(S) \circ \Lambda^n(T)$，故 $\det(ST) \cdot \omega = \det(S) \cdot \det(T) \cdot \omega$。

    2. **$\det(I) = 1$：** $\Lambda^n(I) = \operatorname{id}_{\Lambda^n(V)}$。

    3. **$T$ 可逆 $\Leftrightarrow$ $\det(T) \neq 0$：** $T$ 可逆 $\Leftrightarrow$ $T$ 将基映为基 $\Leftrightarrow$ $T(e_1) \wedge \cdots \wedge T(e_n) \neq 0$。

    4. **行/列相同则行列式为零：** 若 $T$ 的两列相同，则 $T(e_i) = T(e_j)$（$i \neq j$），楔积中出现相同的向量，故为 $0$。

!!! example "例 49.5 (用外积计算行列式)"
    设 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$。列向量 $v_1 = e_1 + 3e_2$，$v_2 = 2e_1 + 4e_2$。

    $$v_1 \wedge v_2 = (e_1 + 3e_2) \wedge (2e_1 + 4e_2) = 4(e_1 \wedge e_2) + 6(e_2 \wedge e_1)$$

    $$= 4(e_1 \wedge e_2) - 6(e_1 \wedge e_2) = -2(e_1 \wedge e_2).$$

    故 $\det(A) = -2$。

!!! definition "定义 49.6 (定向与体积形式)"
    $\Lambda^n(V)$ 的一个非零元素 $\omega$ 确定了 $V$ 的一个**定向**（orientation）。两个非零元素 $\omega, \omega'$ 确定相同的定向当且仅当 $\omega' = c\omega$（$c > 0$，当 $\mathbb{F} = \mathbb{R}$ 时）。

    配备了内积的 $n$ 维实向量空间 $V$ 有一个自然的体积形式：选取正交归一基 $\{e_1, \ldots, e_n\}$，则 $\omega = e_1 \wedge \cdots \wedge e_n$ 是**标准体积形式**。对任意 $v_1, \ldots, v_n \in V$，

    $$v_1 \wedge \cdots \wedge v_n = \det([v_1 \cdots v_n]) \cdot \omega,$$

    其中 $[v_1 \cdots v_n]$ 是以 $v_i$ 为列的矩阵在正交归一基下的坐标。

---

## 49.5 复合矩阵

<div class="context-flow" markdown>

**核心问题**：线性变换 $T$ 对 $k$-向量的诱导作用是什么？如何用矩阵表示？

</div>

!!! definition "定义 49.7 (诱导映射与复合矩阵)"
    设 $T: V \to V$ 是线性变换。$T$ 自然诱导**$k$ 次外幂映射** $\Lambda^k(T): \Lambda^k(V) \to \Lambda^k(V)$：

    $$\Lambda^k(T)(v_1 \wedge \cdots \wedge v_k) = T(v_1) \wedge \cdots \wedge T(v_k).$$

    选取 $V$ 的基 $\{e_1, \ldots, e_n\}$，$\Lambda^k(V)$ 的基 $\{e_I : I \in \binom{[n]}{k}\}$（$I = \{i_1, \ldots, i_k\}$，$e_I = e_{i_1} \wedge \cdots \wedge e_{i_k}$）。$\Lambda^k(T)$ 在此基下的矩阵称为 $T$（或其矩阵 $A$）的 **$k$ 次复合矩阵**（$k$-th compound matrix），记为 $C_k(A)$。

    $C_k(A)$ 是 $\binom{n}{k} \times \binom{n}{k}$ 矩阵，其 $(I, J)$ 元素为 $A$ 的由行 $I$、列 $J$ 确定的 $k \times k$ 子式 $\det(A[I|J])$。

!!! theorem "定理 49.5 (复合矩阵的性质)"
    设 $A, B$ 是 $n \times n$ 矩阵。则：

    1. **乘法性**：$C_k(AB) = C_k(A) C_k(B)$（因为 $\Lambda^k(ST) = \Lambda^k(S) \circ \Lambda^k(T)$）；
    2. **$C_1(A) = A$**，$C_n(A) = \det(A)$；
    3. **$C_k(I) = I_{\binom{n}{k}}$**；
    4. 若 $A$ 可逆，则 $C_k(A)$ 可逆，$C_k(A^{-1}) = C_k(A)^{-1}$。

??? proof "证明"
    性质 1 是函子性的直接推论。$\Lambda^k$ 是一个从向量空间范畴到自身的函子，保持复合。

    **Cauchy-Binet 公式**：$C_k(AB)$ 的 $(I,J)$ 元素为

    $$\det((AB)[I|J]) = \sum_{K \in \binom{[n]}{k}} \det(A[I|K]) \det(B[K|J]),$$

    这正是 $C_k(A) C_k(B)$ 的 $(I,J)$ 元素。因此 Cauchy-Binet 公式不过是 $\Lambda^k(AB) = \Lambda^k(A) \circ \Lambda^k(B)$ 的坐标表达。

!!! theorem "定理 49.6 (复合矩阵的特征值)"
    设 $A$ 是 $n \times n$ 矩阵，特征值为 $\lambda_1, \ldots, \lambda_n$（含重数，在代数闭包中）。则 $C_k(A)$ 的特征值恰好是

    $$\{\lambda_{i_1} \lambda_{i_2} \cdots \lambda_{i_k} : 1 \leq i_1 < i_2 < \cdots < i_k \leq n\},$$

    共 $\binom{n}{k}$ 个。

??? proof "证明"
    不妨设 $A$ 可上三角化（在代数闭域上总是可以的）：$A = PTP^{-1}$，$T$ 上三角，对角元素为 $\lambda_1, \ldots, \lambda_n$。

    $C_k(A) = C_k(PTP^{-1}) = C_k(P) C_k(T) C_k(P)^{-1}$，故 $C_k(A)$ 与 $C_k(T)$ 相似。

    $C_k(T)$ 也是上三角矩阵（这需要适当排列基序），其对角元素为 $\det(T[\{i_1,\ldots,i_k\}|\{i_1,\ldots,i_k\}]) = \lambda_{i_1} \cdots \lambda_{i_k}$（$T$ 上三角故子式也上三角）。

!!! example "例 49.6 (复合矩阵计算)"
    设 $A = \begin{pmatrix} 1 & 2 & 0 \\ 0 & 3 & 1 \\ 0 & 0 & 2 \end{pmatrix}$，特征值 $\lambda_1=1, \lambda_2=3, \lambda_3=2$。

    $C_2(A)$ 是 $3 \times 3$ 矩阵，基为 $e_1 \wedge e_2, e_1 \wedge e_3, e_2 \wedge e_3$。

    $C_2(A)$ 的特征值为 $\lambda_1\lambda_2 = 3$，$\lambda_1\lambda_3 = 2$，$\lambda_2\lambda_3 = 6$。

    验证：$\operatorname{tr}(C_2(A)) = 3 + 2 + 6 = 11$，$\det(C_2(A)) = 3 \cdot 2 \cdot 6 = 36 = (\det A)^2 = 6^2$。

    （一般地，$\det(C_k(A)) = (\det A)^{\binom{n-1}{k-1}}$。）

---

## 49.6 Grassmannian

<div class="context-flow" markdown>

**核心问题**：如何参数化向量空间中所有 $k$ 维子空间的集合？它有什么几何结构？

</div>

!!! definition "定义 49.8 (Grassmannian)"
    设 $V$ 是 $n$ 维向量空间。**Grassmannian** $\operatorname{Gr}(k, V)$（或 $\operatorname{Gr}(k, n)$）是 $V$ 中所有 $k$ 维子空间的集合：

    $$\operatorname{Gr}(k, V) = \{W \subseteq V : W \text{ 是 } k \text{ 维子空间}\}.$$

    特别地，$\operatorname{Gr}(1, V) = \mathbb{P}(V)$ 是**射影空间**。

!!! theorem "定理 49.7 (Grassmannian 的维数)"
    $\operatorname{Gr}(k, n)$ 是一个光滑（实或复）流形，维数为

    $$\dim \operatorname{Gr}(k, n) = k(n - k).$$

??? proof "证明"
    每个 $k$ 维子空间 $W$ 可以用一个 $n \times k$ 矩阵 $X$（列构成 $W$ 的基）来表示。两个矩阵 $X, X'$ 表示同一子空间当且仅当 $X' = XG$，$G \in \operatorname{GL}_k(\mathbb{F})$。

    因此 $\operatorname{Gr}(k,n)$ 可以看作 **Stiefel 流形**（所有 $n \times k$ 列满秩矩阵的集合，维数 $nk - k^2/2 \ldots$ ）的商：

    $$\operatorname{Gr}(k,n) = V_{k,n} / \operatorname{GL}_k(\mathbb{F}),$$

    其中 $V_{k,n}$ 是秩为 $k$ 的 $n \times k$ 矩阵的集合。$V_{k,n}$ 的维数为 $nk$（作为 $\mathbb{F}^{n \times k}$ 的开子集），$\operatorname{GL}_k$ 的维数为 $k^2$，故

    $$\dim \operatorname{Gr}(k,n) = nk - k^2 = k(n-k).$$

    **局部坐标（仿射坐标卡）：** 对 $k$ 元子集 $I \subseteq \{1, \ldots, n\}$，定义开集 $U_I = \{W \in \operatorname{Gr}(k,n) : \pi_I|_W \text{ 是同构}\}$，其中 $\pi_I$ 是投射到第 $I$ 坐标的映射。在 $U_I$ 中，$W$ 可以唯一地用 $k \times (n-k)$ 矩阵参数化（选取 $I$ 坐标为单位矩阵部分，非 $I$ 坐标自由）。$U_I$ 同胚于 $\mathbb{F}^{k(n-k)}$。

!!! example "例 49.7 (低维 Grassmannian)"
    **(a)** $\operatorname{Gr}(1, n) = \mathbb{P}^{n-1}$（射影空间），维数 $1 \cdot (n-1) = n-1$。

    **(b)** $\operatorname{Gr}(2, 4)$：$\mathbb{R}^4$ 中所有 $2$ 维子空间（即平面）的集合，维数 $2 \cdot 2 = 4$。

    **(c)** $\operatorname{Gr}(k, n) \cong \operatorname{Gr}(n-k, n)$（通过取正交补，或对偶地通过 $\Lambda^k(V) \cong \Lambda^{n-k}(V^*)$）。

!!! definition "定义 49.9 (Schubert 胞腔)"
    固定 $V = \mathbb{F}^n$ 的标准完全旗 $\{0\} = F_0 \subset F_1 \subset \cdots \subset F_n = V$（$F_i = \langle e_1, \ldots, e_i \rangle$）。对递增序列 $1 \leq a_1 < a_2 < \cdots < a_k \leq n$，**Schubert 胞腔** $\Omega^\circ_{a_1, \ldots, a_k}$ 定义为

    $$\Omega^\circ_{a_1, \ldots, a_k} = \{W \in \operatorname{Gr}(k,n) : \dim(W \cap F_{a_i}) = i, \, \dim(W \cap F_{a_i - 1}) = i - 1\}.$$

    Schubert 胞腔给出 $\operatorname{Gr}(k,n)$ 的一个 CW-分解，这是研究 Grassmannian 拓扑和上同调的基本工具。

---

## 49.7 Plücker 坐标与嵌入

<div class="context-flow" markdown>

**核心问题**：如何将 Grassmannian 嵌入射影空间？嵌入的像由什么方程描述？

</div>

!!! definition "定义 49.10 (Plücker 嵌入)"
    **Plücker 嵌入**是映射

    $$\iota: \operatorname{Gr}(k, V) \to \mathbb{P}(\Lambda^k(V)),$$

    将 $k$ 维子空间 $W = \langle v_1, \ldots, v_k \rangle$ 映到

    $$\iota(W) = [v_1 \wedge v_2 \wedge \cdots \wedge v_k] \in \mathbb{P}(\Lambda^k(V)).$$

!!! theorem "定理 49.8 (Plücker 嵌入的良定义性和单射性)"
    1. $\iota$ 良定义：$v_1 \wedge \cdots \wedge v_k$ 在基的选取下至多差一个非零标量因子。
    2. $\iota$ 是单射。
    3. $\iota$ 的像恰好是 $\mathbb{P}(\Lambda^k(V))$ 中的**可分解** $k$-向量（构成的射影子簇）。

??? proof "证明"
    **(1) 良定义性：** 若 $\{v_1', \ldots, v_k'\}$ 是 $W$ 的另一组基，则 $v_j' = \sum_i g_{ij} v_i$，$G = (g_{ij}) \in \operatorname{GL}_k(\mathbb{F})$。

    $$v_1' \wedge \cdots \wedge v_k' = \det(G) \cdot v_1 \wedge \cdots \wedge v_k.$$

    由于 $\det(G) \neq 0$，两者在射影空间中定义同一个点。

    **(2) 单射：** 若 $[v_1 \wedge \cdots \wedge v_k] = [w_1 \wedge \cdots \wedge w_k]$，则 $v_1 \wedge \cdots \wedge v_k = c \cdot w_1 \wedge \cdots \wedge w_k$（$c \neq 0$）。可证对任意 $v \in V$：

    $$v \wedge v_1 \wedge \cdots \wedge v_k = 0 \quad \Leftrightarrow \quad v \in \langle v_1, \ldots, v_k \rangle.$$

    因此 $\langle v_1, \ldots, v_k \rangle = \langle w_1, \ldots, w_k \rangle$。

!!! definition "定义 49.11 (Plücker 坐标)"
    选取 $V$ 的基 $\{e_1, \ldots, e_n\}$。$W = \langle v_1, \ldots, v_k \rangle$ 的 **Plücker 坐标**是

    $$p_I = \det(X_I), \quad I \in \binom{[n]}{k},$$

    其中 $X$ 是以 $v_j$ 为列的 $n \times k$ 矩阵，$X_I$ 是 $X$ 的第 $I$ 行构成的 $k \times k$ 子矩阵。

    Plücker 坐标 $(p_I)$ 是齐次坐标，确定 $\mathbb{P}(\Lambda^k(V))$ 中的一个点。

!!! theorem "定理 49.9 (Plücker 关系)"
    Plücker 坐标 $(p_I)_{I \in \binom{[n]}{k}}$ 对应 $\operatorname{Gr}(k,n)$ 中的点当且仅当它满足以下**二次 Plücker 关系**：对任意 $I \in \binom{[n]}{k-1}$，$J \in \binom{[n]}{k+1}$，

    $$\sum_{l=1}^{k+1} (-1)^l \, p_{I \cup \{j_l\}} \, p_{J \setminus \{j_l\}} = 0,$$

    其中 $J = \{j_1, \ldots, j_{k+1}\}$（$j_1 < \cdots < j_{k+1}$）。

!!! example "例 49.8 ($\operatorname{Gr}(2, 4)$ 的 Plücker 嵌入)"
    $\operatorname{Gr}(2, 4)$ 嵌入到 $\mathbb{P}(\Lambda^2(\mathbb{F}^4)) = \mathbb{P}^5$。Plücker 坐标有 $\binom{4}{2} = 6$ 个：$p_{12}, p_{13}, p_{14}, p_{23}, p_{24}, p_{34}$。

    唯一的 Plücker 关系（取 $I = \{i\}$，$J = \{j_1, j_2, j_3\}$，或等价地）为：

    $$p_{12} p_{34} - p_{13} p_{24} + p_{14} p_{23} = 0.$$

    这是 $\mathbb{P}^5$ 中的一个二次超曲面方程。$\operatorname{Gr}(2,4)$ 同构于这个二次超曲面（称为 **Klein 二次曲面**）。

    **验证维数：** $\operatorname{Gr}(2,4)$ 的维数为 $2 \cdot 2 = 4$，而 $\mathbb{P}^5$ 中二次超曲面的维数为 $5 - 1 = 4$，一致。

    **具体例子：** 子空间 $W = \langle (1,0,1,0), (0,1,0,1) \rangle$ 的 Plücker 坐标：

    $$X = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 0 \\ 0 & 1 \end{pmatrix}, \quad p_{12} = 1, \, p_{13} = 0, \, p_{14} = 1, \, p_{23} = -1, \, p_{24} = 0, \, p_{34} = 1.$$

    验证 Plücker 关系：$1 \cdot 1 - 0 \cdot 0 + 1 \cdot (-1) = 0$。

---

## 49.8 万有性质

<div class="context-flow" markdown>

**核心问题**：外代数在什么意义上是"最好的"反对称代数？如何用范畴论语言刻画？

</div>

!!! theorem "定理 49.10 (外代数的万有性质)"
    设 $V$ 是 $\mathbb{F}$-向量空间。外代数 $(\Lambda(V), \iota)$（其中 $\iota: V \to \Lambda(V)$ 是自然嵌入 $V = \Lambda^1(V) \hookrightarrow \Lambda(V)$）满足以下万有性质：

    对任意结合 $\mathbb{F}$-代数 $A$ 和线性映射 $f: V \to A$ 满足 $f(v)^2 = 0$（$\forall v \in V$），存在唯一的代数同态 $\tilde{f}: \Lambda(V) \to A$ 使得 $\tilde{f} \circ \iota = f$：

    $$\begin{CD}
    V @>{\iota}>> \Lambda(V) \\
    @V{f}VV @VV{\exists!\, \tilde{f}}V \\
    & A
    \end{CD}$$

??? proof "证明"
    由 $\Lambda(V) = T(V)/\mathcal{I}$（$\mathcal{I}$ 由 $v \otimes v$ 生成），$f: V \to A$ 自然延拓为代数同态 $\hat{f}: T(V) \to A$（$\hat{f}(v_1 \otimes \cdots \otimes v_k) = f(v_1) \cdots f(v_k)$）。

    条件 $f(v)^2 = 0$ 保证 $\hat{f}(v \otimes v) = f(v)^2 = 0$，故 $\mathcal{I} \subseteq \ker \hat{f}$，$\hat{f}$ 过渡为 $\tilde{f}: \Lambda(V) = T(V)/\mathcal{I} \to A$。唯一性由 $\Lambda(V)$ 由 $\iota(V)$ 生成保证。

!!! theorem "定理 49.11 (外代数的函子性)"
    设 $f: V \to W$ 是线性映射。则存在唯一的分次代数同态

    $$\Lambda(f): \Lambda(V) \to \Lambda(W),$$

    使得 $\Lambda(f)|_V = f$，且对每个 $k$：

    $$\Lambda^k(f)(v_1 \wedge \cdots \wedge v_k) = f(v_1) \wedge \cdots \wedge f(v_k).$$

    进一步：

    1. $\Lambda(\operatorname{id}_V) = \operatorname{id}_{\Lambda(V)}$；
    2. $\Lambda(g \circ f) = \Lambda(g) \circ \Lambda(f)$。

    即 $\Lambda$ 是从 $\mathbb{F}$-向量空间范畴到分次 $\mathbb{F}$-代数范畴的**函子**（functor）。

!!! definition "定义 49.12 (对偶外幂与交替多线性映射)"
    $\Lambda^k(V)$ 的对偶空间 $\Lambda^k(V)^* \cong \Lambda^k(V^*)$ 与 $V$ 上的**交替 $k$-线性函数**空间同构：

    $$\operatorname{Alt}^k(V) = \{f: V^k \to \mathbb{F} : f \text{ 多线性且交替}\} \cong \Lambda^k(V^*).$$

    同构由 $\varphi_1 \wedge \cdots \wedge \varphi_k \mapsto ((v_1, \ldots, v_k) \mapsto \det(\varphi_i(v_j)))$ 给出。

    这正是微分形式的代数基础：$k$-形式是余切空间 $T^*_p M$ 上的 $k$ 次交替多线性函数，即 $\Lambda^k(T^*_p M)$ 的元素。

!!! example "例 49.9 (外积与行列式的联系)"
    设 $\varphi_1, \ldots, \varphi_k \in V^*$。对应的交替 $k$-线性函数为

    $$(\varphi_1 \wedge \cdots \wedge \varphi_k)(v_1, \ldots, v_k) = \det \begin{pmatrix} \varphi_1(v_1) & \cdots & \varphi_1(v_k) \\ \vdots & & \vdots \\ \varphi_k(v_1) & \cdots & \varphi_k(v_k) \end{pmatrix}.$$

    特别地，取 $V = \mathbb{F}^n$，$\varphi_i = e_i^*$（坐标函数），则 $e_1^* \wedge \cdots \wedge e_n^*$ 对应的就是行列式函数。

!!! theorem "定理 49.12 (Hodge 对偶)"
    设 $V$ 是 $n$ 维实内积空间，带正交归一基和定向。**Hodge 星算子** $\star: \Lambda^k(V) \to \Lambda^{n-k}(V)$ 定义为：

    $$\alpha \wedge \star\beta = \langle \alpha, \beta \rangle \, \omega, \quad \forall \alpha, \beta \in \Lambda^k(V),$$

    其中 $\omega = e_1 \wedge \cdots \wedge e_n$ 是体积形式，$\langle \cdot, \cdot \rangle$ 是 $\Lambda^k(V)$ 上由 $V$ 的内积诱导的内积。

    基本性质：$\star \star = (-1)^{k(n-k)} \operatorname{id}$（在实情形）。

!!! example "例 49.10 ($\mathbb{R}^3$ 中的 Hodge 对偶)"
    在 $\mathbb{R}^3$ 中：

    | $\alpha$ | $\star\alpha$ |
    |:---:|:---:|
    | $1$ | $e_1 \wedge e_2 \wedge e_3$ |
    | $e_1$ | $e_2 \wedge e_3$ |
    | $e_2$ | $-e_1 \wedge e_3 = e_3 \wedge e_1$ |
    | $e_3$ | $e_1 \wedge e_2$ |
    | $e_1 \wedge e_2$ | $e_3$ |
    | $e_1 \wedge e_3$ | $-e_2$ |
    | $e_2 \wedge e_3$ | $e_1$ |
    | $e_1 \wedge e_2 \wedge e_3$ | $1$ |

    叉积可以表示为 $u \times v = \star(u \wedge v)$。这解释了为什么叉积只在 $\mathbb{R}^3$ 中自然定义：它依赖于 Hodge 对偶将 2-向量映射为 1-向量，而这只在 $n = 3$ 时两者维数相同（$\binom{3}{2} = 3 = \binom{3}{1}$）。

---

## 49.9 内积（缩并）

<div class="context-flow" markdown>

**核心问题**：如何定义向量对 $k$-形式的"插入"运算？这在微分形式和物理学中为何如此重要？

</div>

!!! definition "定义 49.13 (内积 / 缩并)"
    设 $V$ 是 $n$ 维 $\mathbb{F}$-向量空间，$V^*$ 为其对偶。对 $v \in V$，**内积**（interior product / contraction）是线性映射

    $$\iota_v: \Lambda^k(V^*) \to \Lambda^{k-1}(V^*),$$

    定义为：对 $\alpha \in \Lambda^k(V^*)$ 和 $v_1, \ldots, v_{k-1} \in V$，

    $$(\iota_v \alpha)(v_1, \ldots, v_{k-1}) = \alpha(v, v_1, \ldots, v_{k-1}).$$

    即将 $v$ "插入"到 $\alpha$ 的第一个变元中。约定 $\iota_v: \Lambda^0(V^*) \to 0$（标量上的内积为零）。

    更一般地，若 $V$ 配备了非退化双线性形式，可将 $V$ 与 $V^*$ 等同，定义 $\iota_v: \Lambda^k(V) \to \Lambda^{k-1}(V)$。

!!! theorem "定理 49.13 (内积的性质)"
    设 $v, w \in V$，$\alpha \in \Lambda^p(V^*)$，$\beta \in \Lambda^q(V^*)$。则：

    1. **线性**：$\iota_{av+bw} = a\,\iota_v + b\,\iota_w$；
    2. **幂零性**：$\iota_v \circ \iota_v = 0$；
    3. **分次 Leibniz 法则**（反导子）：

    $$\iota_v(\alpha \wedge \beta) = (\iota_v \alpha) \wedge \beta + (-1)^p \alpha \wedge (\iota_v \beta);$$

    4. **反交换性**：$\iota_v \circ \iota_w = -\iota_w \circ \iota_v$；
    5. **降次**：$\iota_v$ 将 $\Lambda^k$ 映到 $\Lambda^{k-1}$，即 $\iota_v$ 是**次数 $-1$ 的反导子**。

??? proof "证明"
    **(2)** 对可分解元素 $\alpha = \varphi_1 \wedge \cdots \wedge \varphi_k$（$\varphi_i \in V^*$），

    $$\iota_v \alpha = \sum_{i=1}^{k} (-1)^{i-1} \varphi_i(v) \, \varphi_1 \wedge \cdots \wedge \widehat{\varphi_i} \wedge \cdots \wedge \varphi_k,$$

    其中 $\widehat{\varphi_i}$ 表示省略。再施加 $\iota_v$：

    $$\iota_v(\iota_v \alpha) = \sum_{i<j} (-1)^{i-1}(-1)^{j-2} \varphi_i(v)\varphi_j(v)(\cdots) + \sum_{j<i} (-1)^{i-1}(-1)^{j-1} \varphi_i(v)\varphi_j(v)(\cdots).$$

    交换 $i, j$ 的求和指标后，两部分逐项抵消，故 $\iota_v^2 = 0$。

    **(3)** 对可分解元素用归纳法。设 $\alpha = \varphi \in \Lambda^1(V^*) = V^*$。则

    $$\iota_v(\varphi \wedge \beta) = \varphi(v)\beta - \varphi \wedge \iota_v\beta = (\iota_v\varphi)\wedge\beta + (-1)^1 \varphi \wedge (\iota_v\beta).$$

    一般情形通过对 $p$ 归纳和线性扩展得到。

!!! example "例 49.11 (内积的具体计算)"
    在 $\mathbb{R}^3$ 中取标准基 $\{e_1, e_2, e_3\}$ 和对偶基 $\{e^1, e^2, e^3\}$。设 $v = e_1$。

    **(a)** $\iota_{e_1}(e^1 \wedge e^2) = e^1(e_1)\,e^2 - e^2(e_1)\,e^1 = 1 \cdot e^2 - 0 \cdot e^1 = e^2$。

    **(b)** $\iota_{e_1}(e^2 \wedge e^3) = e^2(e_1)\,e^3 - e^3(e_1)\,e^2 = 0$。

    **(c)** $\iota_{e_1}(e^1 \wedge e^2 \wedge e^3) = e^2 \wedge e^3$（将 $e_1$ 插入第一个位置）。

    **(d)** 设 $v = 2e_1 + 3e_2$，$\alpha = e^1 \wedge e^2 \wedge e^3$。则

    $$\iota_v \alpha = 2\,e^2 \wedge e^3 - 3\,e^1 \wedge e^3 = 2\,e^2 \wedge e^3 + 3\,e^3 \wedge e^1.$$

!!! theorem "定理 49.14 (Cartan 魔术公式)"
    设 $M$ 是光滑流形，$X$ 是向量场，$d$ 是外微分算子，$\iota_X$ 是关于 $X$ 的内积。则 **Lie 导数**满足

    $$\mathcal{L}_X = \iota_X \circ d + d \circ \iota_X.$$

    这一公式将三个基本运算——Lie 导数、外微分、内积——统一起来，是微分几何和理论物理中最重要的公式之一。

    进一步，$\mathcal{L}_X$ 满足：

    - $\mathcal{L}_X(\alpha \wedge \beta) = (\mathcal{L}_X \alpha) \wedge \beta + \alpha \wedge (\mathcal{L}_X \beta)$；
    - $[\mathcal{L}_X, \iota_Y] = \iota_{[X,Y]}$；
    - $[\mathcal{L}_X, d] = 0$。

---

## 49.10 Cauchy-Binet 公式

<div class="context-flow" markdown>

**核心问题**：对于非方阵 $A$ 和 $B$，$\det(AB)$ 如何用子式表达？

</div>

!!! theorem "定理 49.15 (Cauchy-Binet 公式)"
    设 $A \in \mathbb{F}^{m \times n}$，$B \in \mathbb{F}^{n \times m}$，$m \leq n$。则

    $$\det(AB) = \sum_{S \in \binom{[n]}{m}} \det(A_S)\det(B_S),$$

    其中 $\binom{[n]}{m}$ 是 $\{1, \ldots, n\}$ 的所有 $m$ 元子集，$A_S$ 是 $A$ 的列限制在 $S$ 上的 $m \times m$ 子矩阵，$B_S$ 是 $B$ 的行限制在 $S$ 上的 $m \times m$ 子矩阵。

    当 $m = n$ 时，唯一的 $S = \{1, \ldots, n\}$，退化为 $\det(AB) = \det(A)\det(B)$。

??? proof "证明"
    **外代数方法**：设 $a_1, \ldots, a_m$ 是 $A$ 的行向量（$a_i \in \mathbb{F}^n$），$b_1, \ldots, b_m$ 是 $B$ 的列向量（$b_j \in \mathbb{F}^n$）。则 $AB$ 的 $(i,j)$ 元素为 $\langle a_i, b_j \rangle$。

    在 $\Lambda^m(\mathbb{F}^n)$ 中，设 $\{e_1, \ldots, e_n\}$ 为标准基。将 $A$ 视为线性映射 $\mathbb{F}^n \to \mathbb{F}^m$，$B$ 视为 $\mathbb{F}^m \to \mathbb{F}^n$。

    设 $f_1, \ldots, f_m$ 是 $\mathbb{F}^m$ 的标准基。则

    $$B(f_j) = \sum_{i=1}^n B_{ij} e_i, \quad j = 1, \ldots, m.$$

    $$\Lambda^m(B)(f_1 \wedge \cdots \wedge f_m) = B(f_1) \wedge \cdots \wedge B(f_m) = \sum_{S \in \binom{[n]}{m}} \det(B_S)\, e_{s_1} \wedge \cdots \wedge e_{s_m}.$$

    类似地，$\Lambda^m(A)$ 作用于 $e_{s_1} \wedge \cdots \wedge e_{s_m}$ 得到 $\det(A_S)\, f_1 \wedge \cdots \wedge f_m$。

    由 $\Lambda^m(AB) = \Lambda^m(A) \circ \Lambda^m(B)$，

    $$\det(AB)\, f_1 \wedge \cdots \wedge f_m = \sum_{S} \det(A_S)\det(B_S)\, f_1 \wedge \cdots \wedge f_m.$$

    比较系数即得。

!!! example "例 49.12 (Cauchy-Binet 公式的应用)"
    设 $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \end{pmatrix}$，$B = A^T$。则 $AB$ 是 $2 \times 2$ 矩阵。

    $S$ 取遍 $\{1,2\}, \{1,3\}, \{2,3\}$：

    - $S = \{1,2\}$：$\det(A_S) = 1 \cdot 5 - 2 \cdot 4 = -3$，$\det(B_S) = -3$；
    - $S = \{1,3\}$：$\det(A_S) = 1 \cdot 6 - 3 \cdot 4 = -6$，$\det(B_S) = -6$；
    - $S = \{2,3\}$：$\det(A_S) = 2 \cdot 6 - 3 \cdot 5 = -3$，$\det(B_S) = -3$。

    $$\det(AA^T) = (-3)^2 + (-6)^2 + (-3)^2 = 9 + 36 + 9 = 54.$$

    验证：$AA^T = \begin{pmatrix} 14 & 32 \\ 32 & 77 \end{pmatrix}$，$\det = 14 \cdot 77 - 32^2 = 1078 - 1024 = 54$。

---

## 49.11 Koszul 复形

<div class="context-flow" markdown>

**核心问题**：外代数如何产生重要的链复形？这与交换代数中的正则序列有什么关系？

</div>

!!! definition "定义 49.14 (Koszul 复形)"
    设 $R$ 是交换环，$x_1, \ldots, x_n \in R$。设 $V$ 为秩 $n$ 的自由 $R$-模，基为 $e_1, \ldots, e_n$。定义映射 $\partial: \Lambda^k(V) \to \Lambda^{k-1}(V)$ 为

    $$\partial(e_{i_1} \wedge \cdots \wedge e_{i_k}) = \sum_{j=1}^{k} (-1)^{j-1} x_{i_j}\, e_{i_1} \wedge \cdots \wedge \widehat{e_{i_j}} \wedge \cdots \wedge e_{i_k}.$$

    即 $\partial = \sum_{i=1}^n x_i \, \iota_{e_i^*}$（用内积语言），其中 $\iota_{e_i^*}$ 是关于对偶基元素的缩并。

    **Koszul 复形** $K_\bullet(x_1, \ldots, x_n)$ 是

    $$0 \to \Lambda^n(V) \xrightarrow{\partial} \Lambda^{n-1}(V) \xrightarrow{\partial} \cdots \xrightarrow{\partial} \Lambda^1(V) \xrightarrow{\partial} \Lambda^0(V) = R \to 0.$$

!!! theorem "定理 49.16 ($\partial^2 = 0$)"
    上述 $\partial$ 满足 $\partial \circ \partial = 0$，即 $K_\bullet$ 确实构成链复形。

??? proof "证明"
    对基元素 $e_{i_1} \wedge \cdots \wedge e_{i_k}$ 计算 $\partial^2$：

    $$\partial^2(e_{i_1} \wedge \cdots \wedge e_{i_k}) = \partial\left(\sum_{j=1}^k (-1)^{j-1} x_{i_j} \, e_{i_1} \wedge \cdots \wedge \widehat{e_{i_j}} \wedge \cdots \wedge e_{i_k}\right)$$

    $$= \sum_{j=1}^k \sum_{\substack{l=1 \\ l \neq j}}^k (-1)^{j-1}(-1)^{l'-1} x_{i_j} x_{i_l}\, e_{i_1} \wedge \cdots \wedge \widehat{e_{i_j}} \wedge \cdots \wedge \widehat{e_{i_l}} \wedge \cdots \wedge e_{i_k},$$

    其中 $l'$ 是 $l$ 在删去第 $j$ 项后的重新编号。对每对 $(j, l)$（$j < l$），出现的两项为

    $$(-1)^{j-1}(-1)^{l-2} x_{i_j} x_{i_l}(\cdots) + (-1)^{l-1}(-1)^{j-1} x_{i_l} x_{i_j}(\cdots).$$

    由于 $(-1)^{j-1+l-2} + (-1)^{l-1+j-1} = (-1)^{j+l-3} + (-1)^{j+l-2} = 0$，故 $\partial^2 = 0$。

!!! theorem "定理 49.17 (Koszul 复形的正合性)"
    设 $R$ 是 Noetherian 交换环，$x_1, \ldots, x_n \in R$。Koszul 复形 $K_\bullet(x_1, \ldots, x_n)$ 是**正合**的（即所有同调群为零）当且仅当 $x_1, \ldots, x_n$ 构成 $R$ 的一个**正则序列**，即：

    1. $(x_1, \ldots, x_n) \neq R$；
    2. 对每个 $i = 1, \ldots, n$，$x_i$ 不是 $R/(x_1, \ldots, x_{i-1})$ 的零因子。

    特别地，当 $R$ 是局部环且 $x_1, \ldots, x_n$ 在极大理想中时，Koszul 复形正合等价于 $x_1, \ldots, x_n$ 是正则序列。

!!! example "例 49.13 (简单的 Koszul 复形)"
    设 $R = \mathbb{F}[x, y]$，$n = 2$，元素为 $x, y$。Koszul 复形为

    $$0 \to R \xrightarrow{\partial_2} R^2 \xrightarrow{\partial_1} R \to 0,$$

    其中 $\partial_2(1) = (-y, x)$（即 $\partial_2(e_1 \wedge e_2) = x \cdot e_2 - y \cdot e_1$），$\partial_1(a, b) = xa + yb$。

    验证 $\partial_1 \circ \partial_2 = 0$：$\partial_1(-y, x) = x(-y) + y(x) = 0$。

    由于 $x, y$ 是 $\mathbb{F}[x,y]$ 的正则序列，此复形正合。其同调计算了 $R/(x,y) \cong \mathbb{F}$ 的自由消解。

---

### 本章总结

外代数 $\Lambda(V)$ 是向量空间 $V$ 上的基本代数结构，它：

- 通过楔积 $\wedge$ 提供了反对称乘法，自然编码面积、体积等几何量；
- 通过顶外幂 $\Lambda^n(V)$ 给出行列式的最内在定义，使行列式的性质成为外代数函子性的直接推论；
- 通过复合矩阵 $C_k(A)$ 将线性变换的特征值乘积信息编码在外幂映射中；
- 通过 Plücker 嵌入将 Grassmannian 实现为射影空间中的代数簇；
- 通过 Hodge 对偶连接不同次的外幂空间，统一了梯度、旋度、散度等微分算子；
- 通过内积（缩并）$\iota_v$ 为微分形式理论提供基础运算，并通过 Cartan 魔术公式统一 Lie 导数、外微分和内积；
- 通过 Cauchy-Binet 公式给出非方阵乘积的行列式的子式展开；
- 通过 Koszul 复形连接外代数与同调代数，正合性等价于正则序列条件。

---

### 习题

!!! exercise "习题 49.1"
    在 $\mathbb{R}^4$ 中，判断以下 2-向量是否可分解：(a) $e_1 \wedge e_2 + e_1 \wedge e_3$；(b) $e_1 \wedge e_2 + e_3 \wedge e_4$。

!!! exercise "习题 49.2"
    计算 $3 \times 3$ 矩阵 $A = \begin{pmatrix} 2 & 1 & 0 \\ 0 & 2 & 1 \\ 0 & 0 & 2 \end{pmatrix}$ 的 $2$ 次复合矩阵 $C_2(A)$，并验证其特征值为 $4, 4, 4$。

!!! exercise "习题 49.3"
    设 $V$ 是 $n$ 维向量空间，$W \subseteq V$ 是 $k$ 维子空间。证明 $v \in W$ 当且仅当 $v \wedge w_1 \wedge \cdots \wedge w_k = 0$（其中 $\{w_1, \ldots, w_k\}$ 是 $W$ 的基）。

!!! exercise "习题 49.4"
    验证 $\operatorname{Gr}(2, 4)$ 的 Plücker 关系对子空间 $W = \langle (1,1,0,0), (0,0,1,1) \rangle$ 成立。

!!! exercise "习题 49.5"
    证明：$\det(C_k(A)) = (\det A)^{\binom{n-1}{k-1}}$。（提示：用特征值。）

!!! exercise "习题 49.6"
    设 $\omega = \sum_{i<j} a_{ij} \, e_i \wedge e_j \in \Lambda^2(\mathbb{F}^n)$。定义反对称矩阵 $A = (a_{ij})$（$a_{ji} = -a_{ij}$，$a_{ii} = 0$）。证明 $\omega$ 可分解当且仅当 $\operatorname{rank}(A) \leq 2$。
