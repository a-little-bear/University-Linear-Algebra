# 第 06 章 特征值与特征向量

<div class="context-flow" markdown>

**前置**：行列式 (Ch03) · 线性变换 (Ch05)

**本章脉络**：特征值与特征向量定义 $\to$ 特征方程与特征多项式 $\to$ 代数重数与几何重数 $\to$ 相似变换 $\to$ 矩阵可对角化的判定定理 $\to$ 特殊矩阵的谱（对称阵、三角阵） $\to$ Cayley-Hamilton 定理 $\to$ 矩阵幂的计算 $\to$ 谱半径与稳定性初步

**延伸**：特征值分析是理解动力系统稳定性、Google PageRank 算法及量子力学能级的关键；它是 Jordan 标准形 (Ch12) 的基础

</div>

特征值和特征向量揭示了线性变换最本质的“不变方向”。一个复杂的矩阵作用在特定向量上，可能仅仅表现为一种缩放。寻找这些缩放因子及其对应的方向，是简化矩阵运算、分析系统长期行为的核心手段。

---

## 06.1 定义与特征方程

!!! definition "定义 06.1 (特征值与特征向量)"
    设 $A$ 是 $n$ 阶方阵。若存在非零向量 $\mathbf{v}$ 和标量 $\lambda$，使得：
    $$A\mathbf{v} = \lambda\mathbf{v}$$
    则称 $\lambda$ 为 $A$ 的**特征值**（Eigenvalue），$\mathbf{v}$ 为对应于 $\lambda$ 的**特征向量**（Eigenvector）。

!!! definition "定义 06.2 (特征多项式)"
    方程 $\det(A - \lambda I) = 0$ 称为 $A$ 的**特征方程**。其左侧 $p(\lambda) = \det(A - \lambda I)$ 是关于 $\lambda$ 的 $n$ 次多项式，称为**特征多项式**。

---

## 06.2 重数与特征空间

!!! definition "定义 06.3 (代数重数与几何重数)"
    1.  **代数重数**：特征值 $\lambda_i$ 作为特征方程根的重数。
    2.  **几何重数**：特征空间 $E_{\lambda_i} = \ker(A - \lambda_i I)$ 的维数，即线性无关特征向量的最大个数。
    **性质**：几何重数 $\leq$ 代数重数。

---

## 06.3 对角化

!!! theorem "定理 06.1 (可对角化判定)"
    $n$ 阶方阵 $A$ 可对角化 $\iff$ $A$ 有 $n$ 个线性无关的特征向量 $\iff$ 每个特征值的几何重数等于代数重数。
    对角化形式为：$P^{-1}AP = \Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$。

---

## 06.4 Cayley-Hamilton 定理

!!! theorem "定理 06.2 (Cayley-Hamilton 定理)"
    每一个方阵都满足其自身的特征方程。即若 $p(\lambda) = \det(A - \lambda I)$，则 $p(A) = O$。
    **应用**：该定理可用于高效计算矩阵的高次幂以及逆矩阵 $A^{-1}$。

---

## 练习题

1. **[计算] 求 $A = \begin{pmatrix} 4 & -5 \\ 2 & -3 \end{pmatrix}$ 的特征值。**
   ??? success "参考答案"
       $\det(A-\lambda I) = (4-\lambda)(-3-\lambda) + 10 = \lambda^2 - \lambda - 2 = 0$。解得 $\lambda_1 = 2, \lambda_2 = -1$。

2. **[特征向量] 求上题中 $\lambda = 2$ 对应的特征向量。**
   ??? success "参考答案"
       解 $(A-2I)\mathbf{v} = 0$，即 $\begin{pmatrix} 2 & -5 \\ 2 & -5 \end{pmatrix}\begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = 0$。解得 $\mathbf{v} = k \begin{pmatrix} 5 \\ 2 \end{pmatrix}$。

3. **[对角化] 若 $A$ 的特征值各不相同，它是否一定可对角化？**
   ??? success "参考答案"
       是的。不同特征值对应的特征向量必然线性无关。

4. **[性质] 证明：若 $A$ 可逆，则其特征值均不为 0。**
   ??? success "参考答案"
       $A\mathbf{v} = \lambda\mathbf{v}$。若 $\lambda = 0$，则 $A\mathbf{v} = \mathbf{0}$。由于 $\mathbf{v} \neq \mathbf{0}$，这说明 $A$ 奇异，与可逆矛盾。

5. **[迹与行列式] 矩阵的迹与特征值之和有什么关系？**
   ??? success "参考答案"
       $\operatorname{tr}(A) = \sum \lambda_i$。

6. **[三角阵] 上三角矩阵的特征值是什么？**
   ??? success "参考答案"
       其主对角线上的元素。

7. **[相似] 证明相似矩阵具有相同的特征多项式。**
   ??? success "参考答案"
       $\det(P^{-1}AP - \lambda I) = \det(P^{-1}(A-\lambda I)P) = \det(P^{-1})\det(A-\lambda I)\det(P) = \det(A-\lambda I)$。

8. **[幂运算] 若 $A = PDP^{-1}$，计算 $A^{10}$。**
   ??? success "参考答案"
       $A^{10} = P D^{10} P^{-1}$。

9. **[C-H定理] 已知 $A^2 - 3A + 2I = O$，若 $A$ 可逆，求 $A^{-1}$。**
   ??? success "参考答案"
       两边同乘 $A^{-1}$：$A - 3I + 2A^{-1} = O \implies A^{-1} = \frac{1}{2}(3I - A)$。

10. **[重数] 举例说明几何重数小于代数重数的情况。**
    ??? success "参考答案"
        $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。特征值 0 的代数重数为 2，但特征向量只有 $\operatorname{span}\{e_1\}$，几何重数为 1。

## 本章小结

特征分析是对方阵最深刻的解构：

1.  **坐标系无关性**：特征值是矩阵的固有属性，不随基的变化而改变，这使其成为物理定律的理想载体。
2.  **解构与重组**：对角化过程实质上是找到了一组完美的基，使线性算子的作用变得纯粹且独立。
3.  **多项式约束**：Cayley-Hamilton 定理建立了矩阵算术与多项式代数之间的终极纽带，揭示了方阵作用的某种自我回归特性。
