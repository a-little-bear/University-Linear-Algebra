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

****

??? success "参考答案"
    
       故 $d_{(2,1)}(A) = 2 a_{11}a_{22}a_{33} - (a_{12}a_{23}a_{31} + a_{13}a_{21}a_{32})$。

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

本章通过群论的视角探讨了矩阵不变量：

1. **表示框架**：将行列式与积和式统一为 immanant 的特殊实例。
2. **谱界限**：确立了半正定算子在不同特征标下的取值层次（Schur 和 Lieb 定理）。
3. **组合计数**：建立了 immanant 与拉丁方、平面匹配理论的深刻联系。
4. **复杂度障碍**：引入了 VP 与 VNP 问题，突显了不同矩阵函数之间巨大的计算鸿沟。
