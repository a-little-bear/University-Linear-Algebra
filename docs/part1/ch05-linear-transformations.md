# 第 05 章 线性变换

<div class="context-flow" markdown>

**前置**：向量空间(Ch4) · 矩阵运算(Ch2)

**本章脉络**：线性变换定义 → 核 (Kernel) 与像 (Image) → 矩阵表示 $T \leftrightarrow A$ → 相似变换与基变换 → 同态与同构 → 投影、反射、旋转的算子形式

**延伸**：线性变换是“动”的视角，它描述了空间如何被拉伸、扭转或压缩

</div>

如果向量空间是静态的舞台，那么线性变换就是在这个舞台上进行的演变。每一个线性变换在给定一组基后，都可以由一个矩阵唯一代表。

---

## 05.1 定义与核心性质

!!! definition "定义 05.1 (线性变换)"
    映射 $T: V \to W$ 是线性的，若对所有 $u, v \in V$ 和标量 $c$ 满足：
    1. $T(u+v) = T(u) + T(v)$
    2. $T(cv) = cT(v)$

!!! theorem "定理 05.3 (基变换公式)"
    设 $T$ 在基 $B$ 下的矩阵为 $A$，在基 $B'$ 下的矩阵为 $A'$。若 $P$ 是从 $B'$ 到 $B$ 的过渡矩阵，则：
    $$A' = P^{-1} A P$$
    这类矩阵被称为**相似矩阵**。

---

## 练习题

1. **[判定] 映射 $T(x, y) = (x+1, y)$ 是线性变换吗？说明理由。**
   ??? success "参考答案"
       不是。线性变换必须满足 $T(0) = 0$。而此处 $T(0, 0) = (1, 0) \neq (0, 0)$。另外它也不满足数乘封闭性。

2. **[矩阵表示] 设 $T: \mathbb{R}^2 \to \mathbb{R}^2$ 是将向量投影到 $x$ 轴的变换。写出其在标准基下的矩阵。**
   ??? success "参考答案"
       $T(1, 0) = (1, 0)$；$T(0, 1) = (0, 0)$。
       故矩阵 $A = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$。

3. **[核与像] 设 $T$ 的矩阵为 $A = \begin{pmatrix} 1 & 2 \\ 2 & 2 \end{pmatrix}$。计算 $\dim(\ker T)$。**
   ??? success "参考答案"
       计算 $\det A = 2-4 = -2 \neq 0$。矩阵满秩，故核空间只包含零向量。
       $\dim(\ker T) = 0$。

4. **[基变换] 已知 $A = \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix}$。若基变换矩阵 $P = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$，求相似矩阵 $A'$。**
   ??? success "参考答案"
       $P^{-1} = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$。
       $A' = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix} \begin{pmatrix} 2 & 0 \\ 0 & 3 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -3 \\ 0 & 3 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix} = \begin{pmatrix} 2 & -1 \\ 0 & 3 \end{pmatrix}$。

5. **[微分算子] 证明求导运算 $D: P_2 \to P_1$（$D(p) = p'$）是一个线性变换。**
   ??? success "参考答案"
       根据微积分法则：$(f+g)' = f' + g'$ 且 $(cf)' = cf'$。满足加法和数乘的线性性，故 $D$ 是线性变换。

6. **[同构] 证明：维数相同的向量空间是同构的。**
   ??? success "参考答案"
       设 $V, W$ 维数均为 $n$。取各自的一组基，定义将一个基映射到另一个基的线性映射。该映射是双射且保持运算，故 $V \cong W$。

7. **[旋转矩阵] 写出在 $\mathbb{R}^2$ 中逆时针旋转 $90^\circ$ 的变换矩阵。**
   ??? success "参考答案"
       $T(1, 0) = (0, 1)$；$T(0, 1) = (-1, 0)$。
       矩阵 $R = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$。

8. **[幂等性] 证明：若 $P$ 是投影算子（$P^2=P$），则 $P$ 的特征值只能是 0 或 1。**
   ??? success "参考答案"
       设 $Pv = \lambda v$，则 $P^2 v = \lambda^2 v$。
       由 $P^2=P$ 得 $\lambda^2 v = \lambda v \implies (\lambda^2-\lambda)v = 0$。
       对于非零向量 $v$，必有 $\lambda^2-\lambda=0$，即 $\lambda=0$ 或 $1$。

9. **[相似性] 相似矩阵具有相同的特征多项式吗？**
   ??? success "参考答案"
       是的。$\det(\lambda I - P^{-1}AP) = \det(P^{-1}(\lambda I - A)P) = \det(P^{-1}) \det(\lambda I - A) \det(P) = \det(\lambda I - A)$。

10. **[秩-零化度应用] 变换 $T: \mathbb{R}^n \to \mathbb{R}^m$ 是单射的条件是什么？**
    ??? success "参考答案"
        条件是核空间只包含零向量，即 $\dim(\ker T) = 0$。由秩-零化度定理，这等价于 $\operatorname{rank}(T) = n$（满列秩）。

## 本章小结

线性变换统一了算子与矩阵：

1. **表现形式**：矩阵是线性变换在特定观测坐标系（基）下的数值化表现。
2. **内在逻辑**：相似变换反映了同一算子在不同视角下的描述，其特征属性保持不变。
3. **结构分类**：核与像的维数分布定义了变换对空间信息的保留或压缩程度。
