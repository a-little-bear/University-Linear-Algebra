# 第 31 章 Majorization 与双随机矩阵

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 凸性(Ch64A) · 双随机矩阵(Ch64A)

**本章脉络**：Majorization 定义 → $\mathbb{R}^n$ 上的偏序 → 双随机矩阵与 Majorization (Hardy-Littlewood-Pólya) → Schur-Horn 定理 → Ky Fan 最大值原理 → Lidskii 定理 → 对数 Majorization → 概率与物理应用

**延伸**：Majorization 是不等式论与不确定性（熵）比较的数学语言；它建立了矩阵对角元与谱之间的深刻联系

</div>

Majorization 提供了一种精确的数学方法来描述一个向量比另一个向量“更分散”。它在 $\mathbb{R}^n$ 向量上定义了一种偏序关系，这种关系在所有双随机变换下保持不变。在矩阵理论中，Majorization 是连接矩阵元素与其谱（特征值）的桥梁，其中最著名的代表是 Schur-Horn 定理。

---

## 31.1 定义与刻画

!!! definition "定义 31.1 (Majorization)"
    设 $x, y \in \mathbb{R}^n$。称 $x$ 被 $y$ **Majorized**（记作 $x \prec y$），如果：
    $$\sum_{i=1}^k x_i^\downarrow \le \sum_{i=1}^k y_i^\downarrow, \quad k=1, \dots, n-1, \quad \text{且 } \sum_{i=1}^n x_i = \sum_{i=1}^n y_i$$
    其中 $x_i^\downarrow$ 是将 $x$ 的分量按非增顺序排列后的结果。

!!! theorem "定理 31.1 (Hardy-Littlewood-Pólya)"
    $x \prec y$ 当且仅当存在双随机矩阵 $D$ 使得 $x = Dy$。

---

## 练习题

1. **[基础] 验证 $(1/2, 1/2) \prec (1, 0)$。**
   ??? success "参考答案"
       和的条件：$1/2 \le 1$ 且 $1/2+1/2 = 1+0 = 1$。条件满足。几何上，$(1/2, 1/2)$ 是 $(1, 0)$ 和 $(0, 1)$ 的平均值，因此更不“极端”。

2. **[凸性] 证明 $x \prec y$ 蕴含对于任何凸函数 $f$，都有 $\sum f(x_i) \le \sum f(y_i)$。**
   ??? success "参考答案"
       由于 $x = Dy$，每个 $x_i$ 都是 $y_j$ 的凸组合。根据 Jensen 不等式，$f(x_i) \le \sum_j d_{ij} f(y_j)$。对 $i$ 求和并利用双随机矩阵列和为 1 的性质 $\sum_i d_{ij} = 1$，即可得证。

3. **[Schur-Horn] 简述 Schur-Horn 定理的内容。**
   ??? success "参考答案"
       对于实对称矩阵 $A$，其对角元向量 $d$ 被其特征值向量 $\lambda$ 所 Majorized：$d \prec \lambda$。反之，给定任何满足 $d \prec \lambda$ 的向量，总能构造一个具有相应对角元和特征值的对称矩阵。

4. **[置换多面体] 描述集合 $\{x : x \prec y\}$ 的几何形状。**
   ??? success "参考答案"
       该集合是 $y$ 的**置换多面体**（Permutahedron）：即由向量 $y$ 的所有 $n!$ 个置换作为顶点构成的凸包。这是双随机矩阵 Birkhoff 定理的直接推论。

5. **[Ky Fan] Ky Fan 最大值原理如何与 Majorization 联系？**
   ??? success "参考答案"
       Ky Fan 原理指出前 $k$ 个特征值之和等于 $k$ 维正交投影下迹的最大值。这为证明矩阵对角元之和受到特征值之和的变分约束提供了基础。

6. **[熵] 证明若概率分布 $p \prec q$，则 Shannon 熵 $H(p) \ge H(q)$。**
   ??? success "参考答案"
       函数 $f(t) = t \log t$ 是凸的。$p \prec q \implies \sum p_i \log p_i \le \sum q_i \log q_i$，即 $-H(p) \le -H(q)$，从而 $H(p) \ge H(q)$。Majorization 过程对应于分布向均匀化（高熵）演向。

7. **[对数] 定义正向量的对数 Majorization $x \prec_{\log} y$。**
   ??? success "参考答案"
       $x \prec_{\log} y$ 指前 $k$ 个最大分量的乘积满足 $\prod_{i=1}^k x_i^\downarrow \le \prod_{i=1}^k y_i^\downarrow$，且在 $k=n$ 时取等号。这等价于 $\log x \prec \log y$。

8. **[奇异值] 两个矩阵乘积 $AB$ 的奇异值与原矩阵奇异值有何 Majorization 关系？**
   ??? success "参考答案"
       满足 $\sigma(AB) \prec_{\log} \sigma(A) \circ \sigma(B)$，其中 $\circ$ 代表分量对应乘积。这说明乘积操作会压缩奇异值的分布。

9. **[Lidskii] 简述 Lidskii 定理关于对称阵之和的特征值结论。**
   ??? success "参考答案"
       $\lambda(A+B) - \lambda(A) \prec \lambda(B)$。它刻画了谱扰动在 Majorization 意义下的范围。

10. **[量子] 解释 Majorization 在量子态演化（LOCC）中的意义。**
    ??? success "参考答案"
        量子态 $\rho$ 能通过局部操作和经典通信（LOCC）转化为 $\sigma$，当且仅当其谱满足 $\operatorname{spec}(\rho) \prec \operatorname{spec}(\sigma)$。这反映了量子纠缠转化中的守恒与耗散规律。

## 本章小结

本章形式化了向量“分散程度”的比较及其矩阵含义：

1. **多样性排序**：定义了 Majorization 这一刻画分量集中程度的偏序关系。
2. **矩阵对角元**：利用 Schur-Horn 定理确立了特征值是所有可能对角元中最“极端”的边界。
3. **变分原理**：将 Majorization 与 Ky Fan 迹极大化联系，提供了特征值和的计算微积分。
4. **信息度量**：论证了 Majorization 是概率与量子理论中比较熵与不确定性的范畴化工具。
