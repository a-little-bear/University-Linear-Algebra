# 第 40B 章 不变量与广义矩阵函数

<div class="context-flow" markdown>

**前置**：行列式(Ch3) · 积和式(Ch40A) · 群论基础 · 对称群的表示论

**本章脉络**：对称群 $S_n$ 的表示论入门 $\to$ Immanant 定义 $\to$ Merris 广义矩阵函数 $\to$ Schur 不等式 $\to$ Lieb 定理与积和式支配猜想 $\to$ Stembridge 猜想 $\to$ 计数的 Kasteleyn 定理 $\to$ VP 与 VNP 猜想

**延伸**：Immanant 通过对称群的表示论统一了行列式与积和式，是代数组合学的核心对象

</div>

行列式和积和式分别对应于对称群 $S_n$ 的两个“极端”特征标。**Immanant** 提供了一个统一的框架，允许使用 $S_n$ 的*任意*不可约特征标 $\chi^\lambda$。本章将形式化这些函数，并探讨它们与对称函数及计算复杂性的深层联系。

---

## 40B.1 对称群 $S_n$ 的表示论

!!! definition "定义 40B.1 (分拆与 Young 图)"
    分拆 $\lambda \vdash n$ 是一个非递增的正整数序列 $(\lambda_1, \dots, \lambda_k)$，其和为 $n$。每个分拆唯一对应 $S_n$ 的一个不可约特征标 $\chi^\lambda$。

!!! theorem "定理 40B.1 (Hook 长度公式)"
    不可约表示 $S^\lambda$ 的维数（记作 $f^\lambda$）由 $f^\lambda = n! / \prod h(i,j)$ 给出，其中 $h(i,j)$ 是 Young 图的钩长。

---

## 40B.2 Immanants

!!! definition "定义 40B.6 (Immanant)"
    对于 $n \times n$ 矩阵 $A$ 和 $S_n$ 的特征标 $\chi$，$\chi$-immanant 定义为：
    $$d_\chi(A) = \sum_{\sigma \in S_n} \chi(\sigma) \prod_{i=1}^n a_{i,\sigma(i)}$$
    - 若 $\chi = \operatorname{sgn}$，则 $d_\chi(A) = \det(A)$。
    - 若 $\chi = \mathbf{1}$，则 $d_\chi(A) = \operatorname{perm}(A)$。

---

## 练习题

1. **[基础] 利用特征标表计算 $3 \times 3$ 矩阵的 $d_{(2,1)}$ immanant。**
   ??? success "参考答案"
       对于 $S_3$，$\chi^{(2,1)}$ 的取值为：$\chi(\text{id})=2, \chi(对换)=0, \chi(3\text{-轮换})=-1$。
       故 $d_{(2,1)}(A) = 2 a_{11}a_{22}a_{33} - (a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32})$。

2. **[对角阵] 求对角矩阵 $D = \operatorname{diag}(x_1, \dots, x_n)$ 的 $d_\lambda(D)$。**
   ??? success "参考答案"
       $d_\lambda(D) = \chi^\lambda(\text{id}) \prod x_i = f^\lambda \prod x_i$。由于只有恒等置换项非零，故结果为维数乘以对角元乘积。

3. **[Schur不等式] 证明对于 $A \succeq 0$，有 $d_\lambda(A) \ge f^\lambda \det(A)$。**
   ??? success "参考答案"
       这是 Schur 定理 (1918)。它源于 immanant 可以表示为对称化子空间上正算子的迹，而行列式对应于该算子特征值的最小值。

4. **[积和式支配] 叙述关于 immanant 的 Lieb 定理。**
   ??? success "参考答案"
       对于 $A \succeq 0$，满足 $|d_\chi(A)| \le \chi(\text{id}) \operatorname{perm}(A)$。这证明了在所有 immanant 中，积和式在半正定矩阵上是取值最大的。

5. **[Stembridge] 什么是 Stembridge 猜想？**
   ??? success "参考答案"
       该猜想指出某些 Jacobi-Trudi 矩阵的 immanant 在展开为对称函数时，其 Schur 函数的系数均为非负。

6. **[拉丁方] 积和式如何用于拉丁方的计数？**
   ??? success "参考答案"
       在 $r \times n$ 拉丁矩形中增加一行的方式数，等于一个描述允许填入数字的 (0,1)-矩阵的积和式。

7. **[Kasteleyn] 为什么平面图的完美匹配数可以在 $O(n^3)$ 时间内计算？**
   ??? success "参考答案"
       Kasteleyn 定理证明可以为平面图定向，使得匹配数等于邻接矩阵的 Pfaffian，而 Pfaffian 是行列式的平方根，可在多项式时间内求解。

8. **[复杂度] 定义类 VP 与 VNP。**
   ??? success "参考答案"
       VP 是可由多项式大小的电路计算的多项式类（如行列式）。VNP 是系数易于计算的多项式类（如积和式）。

9. **[Valiant] 计算积和式是 VNP-完全的吗？**
   ??? success "参考答案"
       是的。Valiant (1979) 证明了积和式是 VNP 中最难的问题。证明它不属于 VP 相当于证明代数版的 $P \neq NP$。

10. **[对称性] 证明对于实特征标，有 $d_\chi(A^T) = d_\chi(A)$。**
    ??? success "参考答案"
        $d_\chi(A^T) = \sum \chi(\sigma) \prod a_{\sigma(i),i} = \sum \chi(\sigma) \prod a_{i,\sigma^{-1}(i)}$。由于 $S_n$ 的特征标满足 $\chi(\sigma) = \chi(\sigma^{-1})$，命题得证。

## 本章小结

本章通过群论的视角探讨了矩阵不变量：

1. **表示框架**：将行列式与积和式统一为 immanant 的特殊实例。
2. **谱界限**：确立了半正定算子在不同特征标下的取值层次（Schur 和 Lieb 定理）。
3. **组合计数**：建立了 immanant 与拉丁方、平面匹配理论的深刻联系。
4. **复杂度障碍**：引入了 VP 与 VNP 问题，突显了不同矩阵函数之间巨大的计算鸿沟。
