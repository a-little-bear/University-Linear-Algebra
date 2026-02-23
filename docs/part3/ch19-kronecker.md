# 第 19 章 Kronecker 积与 Vec 算子

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 线性变换(Ch5) · 张量积初步(Ch21)

**本章脉络**：Kronecker 积定义 $A \otimes B$ → 代数性质（分配律、结合律、转置、逆） → 谱性质（特征值与奇异值） → Vec 算子定义 → 核心恒等式 $\operatorname{vec}(AXB) = (B^T \otimes A)\operatorname{vec}(X)$ → 应用（矩阵方程、系统辨识、张量计算初步）

**延伸**：Kronecker 积是处理矩阵方程的“线性化”工具，它将原本复杂的算子运算转化为标准的向量映射

</div>

Kronecker 积（张量积在矩阵上的体现）提供了一种将两个小矩阵组合成一个大矩阵的系统方法。结合 Vec 算子，它可以将矩阵形式的线性方程转化为标准的向量形式，从而直接应用现有的线性系统求解理论。

---

## 19.1 定义与核心恒等式

!!! definition "定义 19.1 (Kronecker 积)"
    设 $A$ 为 $m \times n$ 矩阵，$B$ 为 $p \times q$ 矩阵。$A \otimes B$ 是一个 $mp \times nq$ 的分块矩阵：
    $$A \otimes B = \begin{pmatrix} a_{11}B & \dots & a_{1n}B \\ \vdots & \ddots & \dots \\ a_{m1}B & \dots & a_{mn}B \end{pmatrix}$$

!!! theorem "定理 19.3 (Vec 算子恒等式)"
    对于矩阵方程 $AXB = C$，其向量化形式为：
    $$(B^T \otimes A) \operatorname{vec}(X) = \operatorname{vec}(C)$$

---

## 练习题

1. **[基础计算] 计算 $\begin{pmatrix} 1 & 2 \end{pmatrix} \otimes \begin{pmatrix} 3 \\ 4 \end{pmatrix}$。**
   ??? success "参考答案"
       根据定义：$\begin{pmatrix} 1\begin{pmatrix} 3 \\ 4 \end{pmatrix} & 2\begin{pmatrix} 3 \\ 4 \end{pmatrix} \end{pmatrix} = \begin{pmatrix} 3 & 6 \\ 4 & 8 \end{pmatrix}$。

2. **[行列式] 若 $A$ 是 $n \times n$ 矩阵，$B$ 是 $m \times m$ 矩阵，$\det(A \otimes B)$ 的公式是什么？**
   ??? success "参考答案"
       $\det(A \otimes B) = (\det A)^m (\det B)^n$。

3. **[特征值] 若 $A$ 的特征值为 $\lambda_i$，$B$ 的特征值为 $\mu_j$，证明 $A \otimes B$ 的特征值为 $\lambda_i \mu_j$。**
   ??? success "参考答案"
       设 $Ax = \lambda x, By = \mu y$。
       则 $(A \otimes B)(x \otimes y) = Ax \otimes By = (\lambda x) \otimes (\mu y) = (\lambda \mu)(x \otimes y)$。
       由于共有 $nm$ 个这样的组合，构成了 $A \otimes B$ 的完整谱。

4. **[求逆] 证明：$(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$（假设逆存在）。**
   ??? success "参考答案"
       $(A \otimes B)(A^{-1} \otimes B^{-1}) = (AA^{-1}) \otimes (BB^{-1}) = I \otimes I = I$。

5. **[Vec应用] 将矩阵方程 $AX + XA^T = C$ 转化为向量形式。**
   ??? success "参考答案"
       $\operatorname{vec}(AX) + \operatorname{vec}(XA^T) = \operatorname{vec}(C)$。
       利用恒等式：$(I \otimes A) \operatorname{vec}(X) + (A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$。
       故 $((I \otimes A) + (A \otimes I)) \operatorname{vec}(X) = \operatorname{vec}(C)$。这就是著名的 Lyapunov 方程的线性化形式。

6. **[秩] 证明 $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A) \operatorname{rank}(B)$。**
   ??? success "参考答案"
       利用奇异值分解：$A \otimes B$ 的奇异值是 $A$ 和 $B$ 奇异值的乘积。非零奇异值的个数即为秩，显然为两矩阵秩的乘积。

7. **[迹] 证明 $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \operatorname{tr}(B)$。**
   ??? success "参考答案"
       $\operatorname{tr}(A \otimes B) = \sum \lambda_{ij} = \sum_i \sum_j \lambda_i \mu_j = (\sum \lambda_i) (\sum \mu_j) = \operatorname{tr}(A) \operatorname{tr}(B)$。

8. **[交换矩阵] 存在置换矩阵 $K$ 使得 $B \otimes A = K (A \otimes B) K^T$ 吗？**
   ??? success "参考答案"
       是的。这个矩阵 $K$ 被称为**交换矩阵 (Commutation Matrix)**。它通过重排条目实现 Kronecker 积的非对易交换。

9. **[混合乘积] 证明 $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$。**
   ??? success "参考答案"
       这是 Kronecker 积最核心的代数性质，可通过分块矩阵乘法直接验证。

10. **[应用] 为什么在量子信息中经常用到 Kronecker 积？**
    ??? success "参考答案"
        因为复合量子系统的态空间（Hilbert 空间）是由各个子系统空间的张量积构成的。若系统 A 有 $n$ 个态，系统 B 有 $m$ 个态，则复合系统 A+B 的态由 $n \times m$ 维向量（或密度矩阵的 Kronecker 积）描述。

## 本章小结

Kronecker 积是高维线性代数的转换桥梁：

1. **维数爆炸与降维**：它提供了一种系统性增加系统维数的方法，同时通过 Vec 算子保持了线性性。
2. **谱的传承**：特征值和奇异值的乘积结构揭示了复合系统的物理本质。
3. **方程求解**：它是将复杂的算子方程（如 Lyapunov, Sylvester 方程）转化为标准线性问题的标准技术路线。
