# 第 13B 章 λ-矩阵与有理标准形

<div class="context-flow" markdown>

**前置**：多项式代数 (Ch00) · 矩阵基础 (Ch02) · Jordan 标准形 (Ch12)

**本章脉络**：λ-矩阵定义 $\to$ 初等变换与等价 $\to$ Smith 标准形 $\to$ 不变因子 $\to$ 初等因子 $\to$ 矩阵相似的充要条件 $\to$ 伴随矩阵 (Companion Matrix) $\to$ 有理标准形 (Rational Canonical Form) $\to$ 结构对比（Jordan vs 有理）

**延伸**：有理标准形解决了在一般数域（如 $\mathbb{Q}$）上矩阵的标准形问题，而无需像 Jordan 形那样依赖代数闭域（如 $\mathbb{C}$）；它是主理想整环上模结构定理的具体应用

</div>

虽然 Jordan 标准形在理论上非常完美，但它的构造依赖于特征值在域内的存在性。如果我们在有理数域 $\mathbb{Q}$ 上工作，特征多项式可能无法分解。有理标准形（Rational Canonical Form）则克服了这一限制，通过多项式环的分解提供了任何域上都通用的矩阵标准形式。

---

## 13B.1 λ-矩阵与 Smith 标准形

!!! definition "定义 13B.1 (λ-矩阵)"
    矩阵元素为关于 $\lambda$ 的多项式的矩阵称为 **λ-矩阵**（或矩阵多项式）。

!!! theorem "定理 13B.1 (Smith 标准形)"
    每一个 $n$ 阶 λ-矩阵 $A(\lambda)$ 都可以通过初等变换化为唯一的对角形：
    $$S(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0)$$
    其中 $d_i(\lambda)$ 是首一多项式，且满足 $d_i(\lambda) \mid d_{i+1}(\lambda)$。这些多项式称为 $A(\lambda)$ 的**不变因子**。

---

## 13B.2 相似的充要条件

!!! theorem "定理 13B.2 (相似判别定理)"
    两个 $n$ 阶矩阵 $A$ 与 $B$ 相似的充要条件是：它们的特征矩阵 $\lambda I - A$ 与 $\lambda I - B$ 有相同的 Smith 标准形（即相同的不变因子）。

---

## 13B.3 伴随矩阵与有理标准形

!!! definition "定义 13B.2 (伴随矩阵 $C(p)$)"
    对于首一多项式 $p(\lambda) = \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_0$，其**伴随矩阵**定义为：
    $$C(p) = \begin{pmatrix} 0 & 0 & \cdots & -a_0 \\ 1 & 0 & \cdots & -a_1 \\ 0 & 1 & \cdots & -a_2 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & -a_{k-1} \end{pmatrix}$$
    **性质**：$C(p)$ 的特征多项式和最小多项式均为 $p(\lambda)$。

!!! theorem "定理 13B.3 (有理标准形)"
    每一个方阵 $A$ 都相似于一个分块对角矩阵，其对角块是不变因子对应的伴随矩阵：
    $$R = \operatorname{diag}(C(d_1), C(d_2), \ldots, C(d_k))$$
    这被称为 $A$ 的**有理标准形**。

---

## 练习题

1. **[Smith形] 计算 $\begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$ 的不变因子。**
   ??? success "参考答案"
       行列式为 $\lambda^2$。最大公因式 $D_1(\lambda) = \gcd(\lambda, 1, 0, \lambda) = 1$。
       故 $d_1 = D_1 = 1$，$d_2 = D_2/D_1 = \lambda^2$。不变因子为 $1, \lambda^2$。

2. **[相似判定] 若 $A, B$ 的特征多项式相同且均为 $\lambda^2$，它们一定相似吗？**
   ??? success "参考答案"
       不一定。例如 $\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ 和 $\begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$。前者不变因子为 $1, \lambda^2$，后者为 $\lambda, \lambda$。

3. **[伴随阵] 写出 $\lambda^2 - 3\lambda + 2$ 的伴随矩阵。**
   ??? success "参考答案"
       $C(p) = \begin{pmatrix} 0 & -2 \\ 1 & 3 \end{pmatrix}$。

4. **[最小多项式] 有理标准形中，哪一个块对应于矩阵的最小多项式？**
   ??? success "参考答案"
       最后一个非平凡不变因子 $d_k(\lambda)$ 对应的伴随矩阵块。

5. **[不变因子] 证明 $d_i(\lambda) \mid d_{i+1}(\lambda)$。**
   ??? success "参考答案"
       这是 Smith 标准形的结构性质，由 $k$ 阶子式最大公因式的定义推导得出。

6. **[对比] Jordan 形与有理标准形的主要区别是什么？**
   ??? success "参考答案"
       Jordan 形将多项式分解到一次因式的幂（初等因子），依赖于域的封闭性；有理标准形直接使用不变因子，不要求分解，适用于任何域。

7. **[初等因子] 若不变因子为 $1, (\lambda-1)(\lambda-2)$，对应的初等因子是什么？**
   ??? success "参考答案"
       $\lambda-1$ 和 $\lambda-2$。

8. **[秩] $\lambda I - A$ 的 Smith 标准形中，非零对角元的个数 $r$ 等于什么？**
   ??? success "参考答案"
       由于 $\lambda I - A$ 总是非奇异的（行列式为 $n$ 次多项式），故 $r = n$。

9. **[计算] 求 $J_2(\lambda_0)$ 的不变因子。**
   ??? success "参考答案"
       $1, (\lambda-\lambda_0)^2$。

10. **[应用] 为什么有理标准形在计算代数中很重要？**
    ??? success "参考答案"
        因为它避开了求特征根的过程（这通常涉及数值近似），只需进行矩阵的初等变换（精确的代数运算）。

## 本章小结

有理标准形提供了矩阵相似类在任意数域下的普适刻画：

1.  **域的独立性**：通过使用不可约多项式的伴随块，有理标准形消除了对复数域的依赖，成为抽象代数处理线性算子的核心。
2.  **多项式逻辑**：不变因子理论揭示了特征矩阵 $\lambda I - A$ 背后深层的模结构，确立了矩阵相似的终极代数判据。
3.  **计算精准度**：相比于不稳定的特征值计算，基于初等变换的 Smith 形构造为精确代数软件提供了稳健的算法基础。
