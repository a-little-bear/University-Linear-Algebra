# 第 17 章 非负矩阵与 Perron-Frobenius 理论

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 矩阵分析(Ch14) · 图论初步(Ch27)

**本章脉络**：非负矩阵与正矩阵定义 → 不可约性 (Irreducibility) → Perron-Frobenius 定理 → Perron 根 $\rho(A)$ 的性质 → 随机矩阵 (Stochastic Matrices) → 矩阵的指数 (Exponent) → 幂法 (Power Method) → 应用（PageRank、人口模型、Markov 链）

**延伸**：Perron-Frobenius 理论揭示了系统的“主导增长模式”，是链接代数结构与系统长期行为的纽带

</div>

非负矩阵广泛存在于概率论、经济学和生物学中。Perron-Frobenius 理论是非负矩阵分析的皇冠，它保证了在一定连通性条件下，矩阵存在一个唯一的、正的主特征向量。

---

## 17.1 Perron-Frobenius 定理

!!! definition "定义 17.1 (不可约矩阵)"
    若非负矩阵 $A$ 的关联有向图是强连通的（任意两点间均有路径），则称 $A$ 是**不可约的**。

!!! theorem "定理 17.3 (Perron-Frobenius 定理)"
    若 $A$ 是不可约非负矩阵，则：
    1. 谱半径 $\lambda = \rho(A)$ 是 $A$ 的一个特征值（Perron 根）。
    2. 存在唯一的（经归一化）正特征向量 $v > 0$ 满足 $Av = \lambda v$。
    3. $\lambda$ 是单特征值。

---

## 练习题

1. **[基础] 矩阵 $A = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 是不可约的吗？它是本原的（Primitive）吗？**
   ??? success "参考答案"
       - 是不可约的，因为有通路 $1 \to 2$ 和 $2 \to 1$。
       - 不是本原的，因为它具有周期性。$A^2 = I, A^3 = A, \dots$ 永远不会变为全正矩阵。其特征值为 $1, -1$，模长均为 1。

2. **[Perron根] 若非负矩阵 $A$ 的每一行之和均为 $s$，证明 $\rho(A) = s$。**
   ??? success "参考答案"
       令 $\mathbf{1} = (1, \dots, 1)^T$。由行和条件得 $A \mathbf{1} = s \mathbf{1}$。故 $s$ 是一个特征值。由于 $A \ge 0$ 且所有特征值模长不超过最大行和（Gershgorin 定理），故 $\rho(A) = s$。

3. **[比较性质] 若 $0 \le A \le B$，证明 $\rho(A) \le \rho(B)$。**
   ??? success "参考答案"
       由于 $A, B$ 是非负的，随着指数增加，$A^k \le B^k$。由 Gelfand 公式 $\rho(A) = \lim \|A^k\|^{1/k} \le \lim \|B^k\|^{1/k} = \rho(B)$。这体现了非负矩阵谱半径的单调性。

4. **[本原性] 判定 $A = \begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix}$ 是否为本原矩阵。**
   ??? success "参考答案"
       计算 $A^2 = \begin{pmatrix} 2 & 1 \\ 1 & 1 \end{pmatrix} > 0$。由于存在某次幂为全正矩阵，故 $A$ 是本原的。

5. **[Collatz-Wielandt] 利用 Collatz-Wielandt 公式估计 $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$ 的谱半径范围。**
   ??? success "参考答案"
       行和分别为 3 和 7。故 $3 \le \rho(A) \le 7$。精确值为 $\frac{5+\sqrt{13}}{2} \approx 5.30$。

6. **[应用] 在 PageRank 算法中，为什么要在原始链接矩阵中加入阻尼因子（加上一个全正矩阵的倍数）？**
   ??? success "参考答案"
       为了使矩阵变为**本原的**。这保证了幂法能够收敛到唯一的正平稳分布（即网页排名向量），消除了原图中孤立节点或死循环导致的收敛问题。

7. **[特征值分布] 若 $A > 0$，除 Perron 根外，其他特征值的模长满足什么条件？**
   ??? success "参考答案"
       对于本原矩阵（特别是全正矩阵），所有其他特征值的模长严格小于 Perron 根：$|\lambda_i| < \rho(A)$（$i \ge 2$）。

8. **[随机矩阵] 证明随机矩阵（行和为 1）的谱半径等于 1。**
   ??? success "参考答案"
       由 $A \mathbf{1} = 1 \mathbf{1}$ 知 1 是特征值。由于行和均为 1，由 Gershgorin 定理知所有特征值模长 $\le 1$。故 $\rho(A) = 1$。

9. **[不可约性判定] 判定 $\begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$ 的不可约性。**
   ??? success "参考答案"
       可约。节点 1 无法到达节点 2（只有 $2 \to 1$ 的路径）。矩阵已经是下三角形式，反映了其可约性。

10. **[极限行为] 若 $A$ 是本原的且 $\rho(A)=1$，证明 $\lim_{k \to \infty} A^k = v w^T$，其中 $v, w$ 为左右 Perron 特征向量。**
    ??? success "参考答案"
        这是由于除 $\lambda=1$ 外所有特征值模长小于 1。在谱分解中，只有对应 $\lambda=1$ 的项（即秩-1 投影算子）在极限下存活。

## 本章小结

非负矩阵理论将分析与组合完美交织：

1. **增长主导**：谱半径不再仅仅是抽象数值，而是系统扩张的实际速率。
2. **正性保证**：不可约性是确保系统中每一个组件都能获得正向演化的结构前提。
3. **收敛必然**：本原性确立了离散动力系统向稳态演化的唯一终点。
