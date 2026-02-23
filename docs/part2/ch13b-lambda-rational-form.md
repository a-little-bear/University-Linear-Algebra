# 第 13B 章 λ-矩阵与有理标准形

<div class="context-flow" markdown>

**前置**：多项式代数 (Ch00) · 矩阵基础 (Ch02) · Jordan 标准形 (Ch12) · 商空间 (Ch13A)

**本章脉络**：从数域到多项式环 $\to$ λ-矩阵定义与初等变换 $\to$ Smith 标准形及其唯一性定理 $\to$ 行列式因子与不变因子 (Invariant Factors) $\to$ 初等因子 $\to$ 矩阵相似的充要条件（特征矩阵等价） $\to$ 伴随矩阵 (Companion Matrix) 构造 $\to$ 有理标准形 (Rational Canonical Form) 定理 $\to$ 与 Jordan 标准形的关系 $\to$ 应用：一般域上的矩阵结构分析、控制系统的正则形式

**延伸**：有理标准形解决了在一般数域（如 $\mathbb{Q}$）上矩阵的标准形问题，而无需像 Jordan 形那样依赖代数闭域（如 $\mathbb{C}$）；它是主理想整环上模结构定理的具体应用

</div>

虽然 Jordan 标准形在理论上非常完美，但它的构造依赖于特征值在域内的存在性。如果我们在有理数域 $\mathbb{Q}$ 上工作，特征多项式可能无法分解。**有理标准形**（Rational Canonical Form）则克服了这一限制，通过多项式环的分解提供了任何域上都通用的矩阵标准形式。本章将利用 λ-矩阵的初等变换揭示这一深刻的代数构造。

---

## 13B.1 λ-矩阵与 Smith 标准形

!!! definition "定义 13B.1 (λ-矩阵)"
    矩阵元素为关于 $\lambda$ 的多项式的矩阵称为 **λ-矩阵**。
    例如特征矩阵 $\lambda I - A$ 就是最典型的 λ-矩阵。

!!! theorem "定理 13B.1 (Smith 标准形)"
    每一个 $n$ 阶 λ-矩阵 $A(\lambda)$ 都可以通过初等变换化为唯一的对角形：
    $$S(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0)$$
    其中 $d_i(\lambda)$ 是首一多项式，且满足 $d_i(\lambda) \mid d_{i+1}(\lambda)$。
    这些多项式称为 $A(\lambda)$ 的**不变因子**（Invariant Factors）。

---

## 13B.2 相似判定与初等因子

!!! theorem "定理 13B.2 (相似判别定理)"
    两个 $n$ 阶矩阵 $A$ 与 $B$ 相似的充要条件是：
    它们的特征矩阵 $\lambda I - A$ 与 $\lambda I - B$ 有相同的 Smith 标准形（即相同的不变因子）。

---

## 13B.3 伴随矩阵与有理标准形

!!! definition "定义 13B.2 (伴随矩阵 $C(p)$)"
    对于首一多项式 $p(\lambda) = \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_0$，其**伴随矩阵**定义为：
    $$C(p) = \begin{pmatrix} 0 & 0 & \cdots & 0 & -a_0 \\ 1 & 0 & \cdots & 0 & -a_1 \\ 0 & 1 & \cdots & 0 & -a_2 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & 1 & -a_{k-1} \end{pmatrix}$$
    **性质**：$C(p)$ 的特征多项式和最小多项式均为 $p(\lambda)$。

!!! theorem "定理 13B.3 (有理标准形)"
    每一个方阵 $A$ 都相似于一个分块对角矩阵，其对角块是不变因子对应的伴随矩阵：
    $$R = \operatorname{diag}(C(d_1), C(d_2), \ldots, C(d_k))$$
    这被称为 $A$ 的**有理标准形**。注意，$d_i$ 越往后次数越高。

---

## 练习题

**1. [Smith形] 计算 $\begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$ 的不变因子。**

??? success "参考答案"
    **步骤：**
    1. **行列式因子**：
       - $D_1(\lambda) = \gcd(\lambda, 1, 0, \lambda) = 1$。
       - $D_2(\lambda) = \det = \lambda^2$。
    2. **不变因子**：
       - $d_1(\lambda) = D_1 = 1$。
       - $d_2(\lambda) = D_2 / D_1 = \lambda^2$。
    **结论**：不变因子为 $1, \lambda^2$。

**2. [相似判定] 若 $A, B$ 的特征多项式相同且均为 $\lambda^2$，它们一定相似吗？**

??? success "参考答案"
    **结论**：不一定。
    **解析**：
    - 特征多项式相同仅意味着不变因子的乘积相同。
    - 考虑 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，不变因子为 $1, \lambda^2$。
    - 考虑 $B = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$，不变因子为 $\lambda, \lambda$。
    - 由于不变因子序列不同，它们不相似。这反映了 Jordan 块的存在差异。

**3. [伴随阵] 写出 $p(\lambda) = \lambda^2 - 3\lambda + 2$ 的伴随矩阵。**

??? success "参考答案"
    **构造：**
    根据定义，$a_1 = -3, a_0 = 2$。
    伴随矩阵为 $\begin{pmatrix} 0 & -a_0 \\ 1 & -a_1 \end{pmatrix} = \begin{pmatrix} 0 & -2 \\ 1 & 3 \end{pmatrix}$。
    验证：该矩阵的特征值为 1 和 2，符合多项式的根。

**4. [最小多项式] 有理标准形中，哪一个块对应于矩阵的最小多项式？**

??? success "参考答案"
    **结论：**
    最后一个非平凡不变因子 $d_k(\lambda)$ 对应的伴随矩阵块。
    **理由**：在 Smith 标准形中，$d_i \mid d_{i+1}$，因此 $d_k$ 包含了所有初等因子的最高幂次，这正是最小多项式的定义。

**5. [不变因子] 证明 $d_i(\lambda) \mid d_{i+1}(\lambda)$。**

??? success "参考答案"
    **证明思路：**
    1. 不变因子的定义是 $d_i = D_i / D_{i-1}$，其中 $D_i$ 是所有 $i$ 阶子式的最大公因式。
    2. 根据行列式按行展开的性质，任何 $i+1$ 阶子式都是 $i$ 阶子式的线性组合。
    3. 因此 $D_i$ 必然整除 $D_{i+1}$。
    4. 深入的代数引理证明了商序列也满足除法链关系。

**6. [对比] Jordan 形与有理标准形的主要区别是什么？**

??? success "参考答案"
    **核心差异：**
    1. **依赖域不同**：Jordan 形要求特征多项式必须在域内完全分解（通常需要复数域 $\mathbb{C}$）；有理标准形在任何数域（如 $\mathbb{Q}, \mathbb{R}$）上都存在且唯一。
    2. **分解深度不同**：Jordan 形将空间分解到不可约的线性项幂（特征空间）；有理标准形分解到不可约多项式的伴随块（循环空间）。

**7. [初等因子] 若不变因子为 $1, (\lambda-1)(\lambda-2)$，对应的初等因子是什么？**

??? success "参考答案"
    **计算：**
    初等因子是不变因子在复数域（或代数闭域）上分解得到的最高次幂项。
    $(\lambda-1)(\lambda-2)$ 分解为 $(\lambda-1)$ 和 $(\lambda-2)$。
    **结论**：初等因子为 $(\lambda-1), (\lambda-2)$。

**8. [秩] $\lambda I - A$ 的 Smith 标准形中，非零对角元的个数 $r$ 等于什么？**

??? success "参考答案"
    **结论：**
    $r = n$。
    **理由**：特征矩阵 $\lambda I - A$ 的行列式是一个 $n$ 次多项式，不恒为 0。因此该 λ-矩阵是满秩的，在 Smith 标准形中必然有 $n$ 个非零对角元。

**9. [计算] 求 $J_2(\lambda_0)$ 的不变因子。**

??? success "参考答案"
    **分析：**
    $A = \begin{pmatrix} \lambda_0 & 1 \\ 0 & \lambda_0 \end{pmatrix} \implies \lambda I - A = \begin{pmatrix} \lambda-\lambda_0 & -1 \\ 0 & \lambda-\lambda_0 \end{pmatrix}$。
    1. $D_1 = \gcd(\lambda-\lambda_0, -1, 0, \lambda-\lambda_0) = 1$。
    2. $D_2 = (\lambda-\lambda_0)^2$。
    **结论**：不变因子为 $1, (\lambda-\lambda_0)^2$。

**10. [应用] 为什么有理标准形在计算代数中很重要？**

??? success "参考答案"
    **工程意义：**
    1. **避免数值近似**：求特征值通常涉及求根运算，在计算机中会有精度损失。而计算有理标准形只需要进行多项式的加减乘除（精确运算）。
    2. **域的普适性**：对于无法求出精确特征值的矩阵（如在 $\mathbb{Q}$ 上），有理标准形提供了描述其结构的唯一精确方式。

## 本章小结

有理标准形提供了矩阵相似类在任意数域下的普适刻画：

1.  **域的独立性**：通过使用不可约多项式的伴随块，有理标准形消除了对复数域的依赖，成为抽象代数处理线性算子的核心。
2.  **多项式逻辑**：不变因子理论揭示了特征矩阵 $\lambda I - A$ 背后深层的模结构，确立了矩阵相似的终极代数判据。
3.  **计算精准度**：相比于不稳定的特征值计算，基于初等变换的 Smith 形构造为精确代数软件提供了稳健的算法基础。
