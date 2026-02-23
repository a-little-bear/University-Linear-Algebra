# 第 33 章 广义逆矩阵

<div class="context-flow" markdown>

**前置**：矩阵秩(Ch2) · 正交投影(Ch5) · 奇异值分解(Ch11) · 最小二乘法(Ch72B)

**本章脉络**：广义逆的定义 → Penrose 四条件 → Moore-Penrose 逆 $A^\dagger$ → 存在性与唯一性 → 通过 SVD 和 QR 计算 → 性质（$A^\dagger A$ 作为投影） → 求解超静定方程组 → 最小范数解

**延伸**：广义逆是数值线性代数的引擎，为秩亏损和超静定系统提供了稳定的解；在统计学（线性模型）和信号处理中不可或缺

</div>

广义逆矩阵将矩阵逆的概念推广到了长方形矩阵和奇异矩阵。虽然标准的逆 $A^{-1}$ 仅对方阵且非奇异矩阵存在，但 **Moore-Penrose 逆** $A^\dagger$ 是唯一的，且对*任何*矩阵都存在。它在最小二乘和最小范数意义下提供了线性系统的“最佳可能”解，建立了 $A^T$ 的列空间到 $A$ 的列空间之间的双射。

---

## 33.1 Moore-Penrose 条件

!!! definition "定义 33.1 (Moore-Penrose 逆)"
    对于任意 $m \times n$ 矩阵 $A$，Moore-Penrose 逆 $A^\dagger$ 是满足以下四个条件的唯一 $n \times m$ 矩阵：
    1. $A A^\dagger A = A$ （内逆）
    2. $A^\dagger A A^\dagger = A^\dagger$ （外逆）
    3. $(A A^\dagger)^* = A A^\dagger$ （左自伴）
    4. $(A^\dagger A)^* = A^\dagger A$ （右自伴）

!!! theorem "Theorem 33.1 (通过 SVD 构造)"
    若 $A = U \Sigma V^*$ 是 $A$ 的奇异值分解，则 $A^\dagger = V \Sigma^\dagger U^*$，其中 $\Sigma^\dagger$ 是将非零奇异值取倒数并转置得到的。

---

## 练习题

1. **[基础] 所有的矩阵（包括非方阵和奇异阵）都有 Moore-Penrose 逆吗？它是唯一的吗？**
   ??? success "参考答案"
       是的，复数域上的所有矩阵 $A$ 都存在唯一的 Moore-Penrose 逆 $A^\dagger$。虽然满足部分条件（如仅满足 $AGA=A$）的广义逆可能不唯一，但同时满足 Penrose 四条件的 $A^\dagger$ 是由矩阵 $A$ 唯一确定的。

2. **[投影] 证明 $P = A A^\dagger$ 是到列空间 $\operatorname{Im}(A)$ 的正交投影算子。**
   ??? success "参考答案"
       由条件 1 和 3，$P^2 = (A A^\dagger A) A^\dagger = A A^\dagger = P$ 且 $P^* = P$。因此 $P$ 是正交投影。由于 $P A = A$，其像空间包含 $\operatorname{Im}(A)$；又由于 $P = A(A^\dagger)$，其像空间包含于 $\operatorname{Im}(A)$。

3. **[满秩] 推导列满秩矩阵 $A$ 的广义逆 $A^\dagger$ 公式。**
   ??? success "参考答案"
       若 $A$ 列满秩，则 $A^* A$ 可逆。此时 $A^\dagger = (A^* A)^{-1} A^*$。这就是最小二乘法中标准的“左逆”。

4. **[线性系统] 证明 $x = A^\dagger b$ 是最小二乘问题 $\min \|Ax - b\|^2$ 的最小范数解。**
   ??? success "参考答案"
       最小二乘要求 $Ax = P_{\operatorname{Im}(A)} b = A A^\dagger b$。$x = A^\dagger b$ 显然满足。任何其他解 $x'$ 与 $x$ 之差属于 $\ker(A)$。由于 $A^\dagger b \in \operatorname{Im}(A^*) = (\ker A)^\perp$，根据勾股定理，$\|A^\dagger b\|$ 必然是所有解中范数最小的。

5. **[性质] $(AB)^\dagger = B^\dagger A^\dagger$ 总是成立吗？**
   ??? success "参考答案"
       不一定。该等式成立的充要条件包括：$A$ 列满秩且 $B$ 行满秩，或者 $A^* A B B^*$ 是 Hermitian 矩阵。通常情况下，由于列空间和行空间的相互作用，该恒等式不成立。

6. **[秩-1] 求秩-1 矩阵 $A = uv^*$ 的伪逆。**
   ??? success "参考答案"
       $A^\dagger = \frac{1}{\|u\|^2 \|v\|^2} v u^*$。

7. **[求逆] 证明 $(A^\dagger)^\dagger = A$。**
   ??? success "参考答案"
       $A^\dagger$ 满足作为 $A$ 的伪逆的四个条件，对称地，$A$ 也满足作为 $A^\dagger$ 的伪逆的四个条件。

8. **[计算] 求 $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ 的 $A^\dagger$。**
   ??? success "参考答案"
       $A$ 列满秩，$A^* A = (2)$。因此 $A^\dagger = (1/2) \begin{pmatrix} 1 & 1 \end{pmatrix} = \begin{pmatrix} 0.5 & 0.5 \end{pmatrix}$。

9. **[零空间] 广义逆 $A^\dagger$ 的零空间与 $A$ 的伴随矩阵有什么关系？**
   ??? success "参考答案"
       $\ker(A^\dagger) = \ker(A^*)$。伪逆将 $A$ 的像空间映射回行空间，并将与像空间正交的噪声部分（属于 $\ker(A^*)$）映射为零。

10. **[连续性] 映射 $A \mapsto A^\dagger$ 是连续的吗？**
    ??? success "参考答案"
        不是。在矩阵秩发生改变的地方，该映射是不连续的。微小的扰动可能产生极小的非零奇异值，其倒数 $1/\sigma$ 会趋于无穷大。这也是为什么在实践中经常使用截断 SVD（Truncated SVD）。

## 本章小结

本章确立了 Moore-Penrose 逆作为通用逆矩阵的地位：

1. **解析唯一性**：通过 Penrose 四条件定义了 $A^\dagger$，确保了任何线性算子都有唯一的“最佳”逆。
2. **几何映射**：展示了 $A^\dagger$ 如何将投影到值域与逆映射到行空间有效结合。
3. **优化求解**：确立了伪逆作为求解最小范数最小二乘问题的数学工具。
4. **数值特性**：关联了伪逆与 SVD，强调了秩亏损情况下的敏感性与稳定性挑战。
