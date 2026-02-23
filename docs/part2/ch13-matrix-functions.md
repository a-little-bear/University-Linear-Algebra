# 第 13 章 矩阵函数

<div class="context-flow" markdown>

**前置**：特征分解(Ch10) · Jordan 标准形(Ch12) · 幂级数

**本章脉络**：矩阵幂级数定义 → 收敛性 → 矩阵指数 $e^A$ → 矩阵对数与三角函数 → 通过对角化计算 → 通过 Jordan 形计算 → Cauchy 积分公式在矩阵上的推广 → 应用（线性微分方程组求解）

**延伸**：矩阵指数是现代控制理论与量子力学的核心算子，它将代数加法映射为群乘法

</div>

矩阵函数将标量函数 $f(z)$ 的定义域扩展到了方阵。这种扩展不是按元素计算，而是保持代数结构的一致性。最核心的工具是**矩阵指数** $e^A$，它是解决一切线性连续动力系统的关键。

---

## 13.1 定义与计算方法

!!! definition "定义 13.1 (矩阵幂级数)"
    若标量函数 $f(z) = \sum_{k=0}^\infty a_k z^k$ 在某个圆盘内收敛，则矩阵函数定义为：
    $$f(A) = \sum_{k=0}^\infty a_k A^k$$

!!! theorem "定理 13.3 (对角化计算法)"
    若 $A = PDP^{-1}$，则 $f(A) = P f(D) P^{-1} = P \operatorname{diag}(f(\lambda_i)) P^{-1}$。

---

## 练习题

1. **[基础计算] 若 $A = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$，计算 $e^A$。**
   ??? success "参考答案"
       $e^A = \begin{pmatrix} e^1 & 0 \\ 0 & e^2 \end{pmatrix}$。对角阵的函数只需对对角元分别求函数值。

2. **[指数性质] 证明：若 $AB = BA$，则 $e^{A+B} = e^A e^B$。**
   ??? success "参考答案"
       利用幂级数展开并应用二项式定理。由于 $A, B$ 对易，$(A+B)^k$ 可以像标量一样展开。若不对易，该等式通常不成立（需使用 BCH 公式）。

3. **[幂零矩阵] 计算 $e^{At}$ 其中 $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。**
   ??? success "参考答案"
       由于 $A^2 = 0$，级数展开为 $e^{At} = I + At + 0 = \begin{pmatrix} 1 & t \\ 0 & 1 \end{pmatrix}$。

4. **[迹恒等式] 证明 $\det(e^A) = e^{\operatorname{tr}(A)}$。**
   ??? success "参考答案"
       设 $A$ 的特征值为 $\lambda_i$，则 $e^A$ 的特征值为 $e^{\lambda_i}$。
       $\det(e^A) = \prod e^{\lambda_i} = e^{\sum \lambda_i} = e^{\operatorname{tr}(A)}$。

5. **[反余弦] 若 $A$ 是投影矩阵（$A^2=A$），计算 $\sin(\pi A)$。**
   ??? success "参考答案"
       $A^k = A$ 对所有 $k \ge 1$ 成立。
       $\sin(\pi A) = \sum \frac{(-1)^k}{(2k+1)!} (\pi A)^{2k+1} = A \sum \frac{(-1)^k \pi^{2k+1}}{(2k+1)!} = A \sin(\pi) = 0$。

6. **[Jordan块函数] 写出 $f(J_2(\lambda))$ 的公式。**
   ??? success "参考答案"
       $f \begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix} = \begin{pmatrix} f(\lambda) & f'(\lambda) \\ 0 & f(\lambda) \end{pmatrix}$。这反映了非对角化分量与导数的关系。

7. **[求逆] $e^A$ 总是可逆的吗？**
   ??? success "参考答案"
       是的。因为其行列式 $\det(e^A) = e^{\operatorname{tr}(A)}$ 永远不为 0。其逆矩阵为 $e^{-A}$。

8. **[微分] 证明 $\frac{d}{dt} e^{At} = A e^{At}$。**
   ??? success "参考答案"
       对级数 $\sum \frac{t^k A^k}{k!}$ 逐项求导得 $\sum \frac{k t^{k-1} A^k}{k!} = A \sum \frac{t^{k-1} A^{k-1}}{(k-1)!} = A e^{At}$。

9. **[周期性] 计算 $e^{Jt}$ 其中 $J = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$。**
   ??? success "参考答案"
       由于 $J^2 = -I$，级数按奇偶项分开即得旋转矩阵：$e^{Jt} = \begin{pmatrix} \cos t & -\sin t \\ \sin t & \cos t \end{pmatrix}$。

10. **[应用] 如何用矩阵指数求解 $\dot{x} = Ax, x(0)=x_0$？**
    ??? success "参考答案"
        其解为 $x(t) = e^{At} x_0$。矩阵指数将初始状态映射到任意时刻的状态。

## 本章小结

矩阵函数是线性代数的高级微积分：

1. **结构映射**：它将标量域的解析性质完美平移到算子空间。
2. **计算核心**：对角化和 Jordan 链是计算任意矩阵函数的标准化路径。
3. **动力学价值**：矩阵指数是连续系统的“时间算子”，是描述系统演化的终极工具。
