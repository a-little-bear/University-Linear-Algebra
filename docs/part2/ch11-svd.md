# 第 11 章 奇异值分解 (SVD)

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 正交性(Ch7) · 矩阵分解(Ch10)

**本章脉络**：SVD 定义 $A = U \Sigma V^*$ → 几何意义（超椭圆映射） → 奇异值与特征值的联系 → 截断 SVD 与低秩逼近 (Eckart-Young 定理) → SVD 求解广义逆 → 应用（图像压缩、主成分分析 PCA、协同过滤）

**延伸**：SVD 是线性代数中最具通用性的分解定理，它对长方形矩阵、奇异矩阵均有效，被誉为“线性代数的基石”

</div>

如果特征分解是对方阵内禀结构的探索，那么 **奇异值分解 (SVD)** 就是对任意形状算子的全面解析。它揭示了任何线性映射都可以看作：旋转 → 各向异性缩放 → 旋转。在数据科学时代，SVD 是提取高维数据核心特征的终极利器。

---

## 11.1 定义与核心定理

!!! definition "定义 11.1 (奇异值分解)"
    对于任意 $m \times n$ 复矩阵 $A$，存在 $m$ 阶酉矩阵 $U$ 和 $n$ 阶酉矩阵 $V$，使得：
    $$A = U \Sigma V^*$$
    其中 $\Sigma$ 是对角线上为非负实数 $\sigma_1 \ge \sigma_2 \dots \ge \sigma_r > 0$ 的对角阵。

!!! theorem "定理 11.3 (Eckart-Young 定理)"
    在 Frobenius 范数下，$A$ 的最优秩-$k$ 逼近是由截断 SVD 给出的：
    $$A_k = \sum_{i=1}^k \sigma_i u_i v_i^*$$

---

## 练习题

1. **[基础计算] 计算 $A = \begin{pmatrix} 3 & 0 \\ 0 & -2 \end{pmatrix}$ 的奇异值。**
   ??? success "参考答案"
       奇异值是 $A^* A$ 特征值的平方根。
       $A^* A = \begin{pmatrix} 9 & 0 \\ 0 & 4 \end{pmatrix}$。特征值为 9 和 4。
       故奇异值为 $\sigma_1 = 3, \sigma_2 = 2$。注意奇异值总是非负的。

2. **[SVD与秩] 若 $A$ 有 3 个非零奇异值，其秩是多少？**
   ??? success "参考答案"
       秩为 3。非零奇异值的个数精确等于矩阵的秩。

3. **[几何意义] 线性变换 $A = U \Sigma V^*$ 将单位圆映射为什么形状？**
   ??? success "参考答案"
       映射为一个**超椭圆**。半轴长度由奇异值 $\sigma_i$ 决定，半轴方向由左奇异向量 $u_i$ 决定。

4. **[谱范数] 证明矩阵的谱范数 $\|A\|_2$ 等于其最大奇异值 $\sigma_{\max}$。**
   ??? success "参考答案"
       $\|A\|_2 = \max \frac{\|Ax\|}{\|x\|} = \max \sqrt{\frac{x^* A^* A x}{x^* x}}$。
       由瑞利商可知，最大值为 $\sqrt{\lambda_{\max}(A^* A)} = \sigma_{\max}$。

5. **[条件数] 定义利用奇异值表示的矩阵条件数 $\kappa(A)$。**
   ??? success "参考答案"
       $\kappa(A) = \sigma_{\max} / \sigma_{\min}$。它衡量了矩阵在求逆运算中对误差的放大程度。

6. **[低秩逼近] 若 $A$ 的奇异值为 $100, 50, 1, 0.1$，其秩-2 逼近的误差（Frobenius 范数）是多少？**
   ??? success "参考答案"
       误差为剩余奇异值的平方和的平方根：$\sqrt{1^2 + 0.1^2} \approx 1.005$。

7. **[特征值对比] 若 $A$ 是正定对称阵，其奇异值与特征值有何关系？**
   ??? success "参考答案"
       对于正定对称阵，奇异值精确等于特征值。

8. **[计算] 求 $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$ 的 SVD。**
   ??? success "参考答案"
       $A^* A = (2)$，特征值 2，故奇异值 $\sigma_1 = \sqrt{2}$。
       $V = (1)$。$u_1 = Av_1/\sigma_1 = \begin{pmatrix} 1/\sqrt{2} \\ 1/\sqrt{2} \end{pmatrix}$。
       $A = \begin{pmatrix} 1/\sqrt{2} \\ 1/\sqrt{2} \end{pmatrix} \begin{pmatrix} \sqrt{2} \end{pmatrix} (1)$。

9. **[应用] 在图像压缩中，保留前 $k$ 个奇异值如何实现压缩？**
   ??? success "参考答案"
       原始图像需要 $m \times n$ 个存储空间。保留前 $k$ 个奇异值及其对应的 $u, v$ 向量仅需 $k(m+n+1)$ 个空间。当 $k \ll \min(m,n)$ 时，压缩率极大。

10. **[伪逆] 利用 SVD 写出 $A$ 的 Moore-Penrose 伪逆 $A^\dagger$。**
    ??? success "参考答案"
        $A^\dagger = V \Sigma^\dagger U^*$，其中 $\Sigma^\dagger$ 是对非零奇异值取倒数并转置得到的。

## 本章小结

SVD 是现代应用数学的皇冠：

1. **普适性**：打破了方阵的局限，为一切线性系统提供了统一视图。
2. **能量集中**：通过奇异值的大小分布，揭示了数据中的主导成分。
3. **鲁棒性**：作为低秩逼近的核心，SVD 是噪声过滤和信息压缩的代数准则。
