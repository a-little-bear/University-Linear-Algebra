# 第 47A 章 矩阵微积分基础

<div class="context-flow" markdown>

**前置**：矩阵运算(Ch2) · 多元微积分 · Kronecker 积(Ch19) · 矩阵范数(Ch15)

**本章脉络**：矩阵函数的微分 → 梯度、Jacobian 与 Hessian → 微分的迹技巧 → 逆矩阵与行列式的导数 → 矩阵链式法则 → Vec 算子与 Kronecker 积 → 矩阵流形上的优化 → 在机器学习（反向传播）中的应用

**延伸**：矩阵微积分是现代深度学习（神经网络权重优化）和高斯过程回归的数学引擎

</div>

矩阵微积分将微分规则扩展到了以矩阵为变量的函数。与其处理 $n^2$ 个独立的偏导数，我们不如将矩阵视为一个整体变量。其核心技术是**微分法**：利用恒等式 $df = \operatorname{tr}(G^T dX)$ 来识别梯度 $G$。这种方法与布局（Layout）无关，能够极大简化涉及逆矩阵、行列式和迹的复杂表达式。

---

## 47A.1 微分与梯度

!!! definition "定义 47A.1 (标量函数的梯度)"
    对于标量函数 $f(X)$（其中 $X \in M_{m \times n}$），**梯度** $\nabla_X f$ 是同型的矩阵：
    $$(\nabla_X f)_{ij} = \frac{\partial f}{\partial X_{ij}}$$
    微分与梯度的关系为：$df = \langle \nabla_X f, dX \rangle = \operatorname{tr}((\nabla_X f)^T dX)$。

!!! theorem "定理 47A.1 (三大基础微分公式)"
    1. $d(AXB) = A(dX)B$
    2. $d(X^{-1}) = -X^{-1}(dX)X^{-1}$
    3. $d(\det X) = \det X \cdot \operatorname{tr}(X^{-1} dX)$

---

## 练习题

1. **[基础] 求 $f(x) = x^T A x$ 的梯度。**
   ??? success "参考答案"
       $df = d(x^T A x) = (dx)^T A x + x^T A (dx) = x^T A^T dx + x^T A dx = x^T (A + A^T) dx$。因此 $\nabla_x f = (A + A^T) x$。若 $A$ 是对称的，则 $\nabla_x f = 2Ax$。

2. **[迹技巧] 求 $f(X) = \operatorname{tr}(AX)$ 的梯度。**
   ??? success "参考答案"
       $df = \operatorname{tr}(A dX)$。根据定义 $df = \operatorname{tr}(G^T dX)$，我们有 $G^T = A$，故 $\nabla_X f = A^T$。

3. **[逆矩阵] 计算 $f(X) = \operatorname{tr}(X^{-1}A)$ 的导数。**
   ??? success "参考答案"
       $df = \operatorname{tr}(d(X^{-1})A) = \operatorname{tr}(-X^{-1}(dX)X^{-1}A) = \operatorname{tr}(-X^{-1}AX^{-1} dX)$。因此 $\nabla_X f = -(X^{-1}AX^{-1})^T = -X^{-T}A^TX^{-T}$。

4. **[对数行列式] 推导 $f(X) = \log\det X$ 对 $X \succ 0$ 的梯度。**
   ??? success "参考答案"
       $df = d(\log\det X) = \frac{1}{\det X} d(\det X) = \frac{1}{\det X} (\det X \operatorname{tr}(X^{-1} dX)) = \operatorname{tr}(X^{-1} dX)$。因此 $\nabla_X f = X^{-T} = X^{-1}$（由于 $X$ 对称）。

5. **[链式法则] 求 $f(X) = \|AX-B\|_F^2$ 的梯度。**
   ??? success "参考答案"
       令 $R = AX-B$，则 $f = \operatorname{tr}(R^T R)$。$df = 2 \operatorname{tr}(R^T dR) = 2 \operatorname{tr}(R^T A dX)$。因此 $\nabla_X f = (2 R^T A)^T = 2 A^T (AX-B)$。令此梯度为零即得最小二乘的正规方程。

6. **[Vec算子] 利用恒等式 $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ 求 $f(X) = AXB$ 的 Jacobian。**
   ??? success "参考答案"
       $\frac{\partial \operatorname{vec}(AXB)}{\partial \operatorname{vec}(X)} = B^T \otimes A$。这允许我们将矩阵微积分转化为标准的向量微积分运算。

7. **[乘法规则] 证明 $d(XY) = (dX)Y + X(dY)$。**
   ??? success "参考答案"
       根据增量定义：$(X+dX)(Y+dY) - XY = XY + (dX)Y + X(dY) + (dX)(dY) - XY$。忽略二阶无穷小量 $(dX)(dY)$ 即得结果。

8. **[Hessian] 定义标量函数 $f(X)$ 的 Hessian 矩阵。**
   ??? success "参考答案"
       Hessian 是二阶偏导数矩阵。在矩阵微积分中，它通常表示为一个双线性算子，使得 $d^2 f = \operatorname{vec}(dX)^T \mathbf{H} \operatorname{vec}(dX)$。

9. **[特征值] 对于对称矩阵，第 $i$ 个特征值的微分 $d\lambda_i$ 是什么？**
   ??? success "参考答案"
       $d\lambda_i = v_i^T (dA) v_i$，其中 $v_i$ 是对应的单位特征向量。该结果是振动分析和统计学中灵敏度分析的基础。

10. **[复杂度] 为什么微分法在大规模矩阵运算中优于逐元素求导法？**
    ??? success "参考答案"
        它将矩阵视为单一的代数对象，保持了结构的完整性（如迹的循环性质），从而产生与坐标系和布局无关的结果，更易于验证并在高级编程语言中实现。

## 本章小结

本章建立了以矩阵为变量的微积分体系：

1. **微分法核心**：开发了微分法作为识别梯度和 Jacobian 的标准工具。
2. **结构恒等式**：推导了包括逆矩阵和行列式在内的基础矩阵算子的导数。
3. **Kronecker 桥梁**：利用 vec-Kronecker 恒等式将矩阵微积分与向量优化联系起来。
4. **优化引擎**：展示了矩阵导数在最小二乘、似然极大化及神经网络训练中的应用。
