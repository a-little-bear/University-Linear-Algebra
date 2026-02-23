# 第 64A 章 矩阵空间中的凸集

<div class="context-flow" markdown>

**前置**：正定矩阵(Ch16) · 矩阵不等式(Ch18) · 优化基础(Ch25) · 共正矩阵(Ch45b)

**本章脉络**：维数与迹内积 → PSD 锥 → Birkhoff 多面体 → 相关矩阵集 (椭圆体) → CP 与 COP 锥 → 谱面体 → 分离定理 → S-引理 → 极小极大定理 → SDP 对偶性

**延伸**：矩阵空间中的凸集理论是半定规划 (SDP) 和量子信息论（纠缠见证、状态空间结构）的数学基础

</div>

凸性是优化理论的核心支柱。当从 $\mathbb{R}^n$ 推广到矩阵空间时，凸性揭示了丰富的几何结构。半正定锥 $S_n^+$ 是最基础的矩阵凸锥，其自对偶性确立了半定规划的强大框架。

---

## 64A.1 核心矩阵凸集

!!! definition "定义 64A.3 (PSD 锥)"
    $S_n^+ = \{A \in S_n(\mathbb{R}) : A \succeq 0\}$ 是实对称矩阵空间中的一个闭凸锥。其自对偶性意味着对所有 $A, B \in S_n^+$，满足 $\langle A, B \rangle = \operatorname{tr}(AB) \ge 0$。

!!! theorem "定理 64A.6 (Birkhoff 定理)"
    $n \times n$ 双随机矩阵的集合 $\mathcal{B}_n$ 是一个凸多面体，其顶点恰好是 $n!$ 个置换矩阵。

!!! definition "定义 64A.10 (谱面体)"
    谱面体是线性矩阵不等式 (LMI) 的可行集：$\{x \in \mathbb{R}^d : A_0 + \sum x_i A_i \succeq 0\}$。它是 SDP 研究的基础对象。

---

## 练习题

1. **[基础] 证明 PSD 锥 $S_n^+$ 是凸的。**
   ??? success "参考答案"
       设 $A, B \in S_n^+$ 且 $t \in [0, 1]$。对于任意向量 $x$，有 $x^T(tA + (1-t)B)x = t(x^T Ax) + (1-t)(x^T Bx)$。由于 $x^T Ax \ge 0$ 且 $x^T Bx \ge 0$，且权重 $t, 1-t \ge 0$，故该和项保持非负。因此凸组合仍属于 $S_n^+$。

2. **[迹内积] 计算 $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ 与 $B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ 的内积 $\langle A, B \rangle$。**
   ??? success "参考答案"
       $\langle A, B \rangle = \operatorname{tr}(A^T B) = \operatorname{tr}(AB) = 1(0) + 2(1) + 2(1) + 1(0) = 4$。

3. **[Birkhoff] 描述所有 $2 \times 2$ 双随机矩阵的一般形式，并验证其维数。**
   ??? success "参考答案"
       约束条件为：$a_{11}+a_{12}=1, a_{21}+a_{22}=1, a_{11}+a_{21}=1, a_{12}+a_{22}=1$。设 $a_{11}=t \in [0,1]$，则 $A = \begin{pmatrix} t & 1-t \\ 1-t & t \end{pmatrix}$。维数为 $(2-1)^2 = 1$。

4. **[分离性] 若对称矩阵 $A$ 具有负特征值 $-1$，找出一个 $H \succeq 0$ 使得 $\operatorname{tr}(AH) < 0$。**
   ??? success "参考答案"
       令 $v$ 为与 $-1$ 关联的单位特征向量。取 $H = vv^T \succeq 0$。则 $\operatorname{tr}(AH) = v^T Av = -1 < 0$。在几何上，$H$ 代表了一个将 $A$ 与 PSD 锥分离的超平面。

5. **[S-引理] 叙述 S-引理在控制理论中的意义。**
   ??? success "参考答案"
       S-引理允许将“一个二次约束蕴含另一个二次约束”的条件转化为单个 LMI。这使得利用高效的凸优化求解器来验证非线性系统的稳定性成为可能。

6. **[谱面体] 识别集合 $\{ (x, y) : \begin{pmatrix} 1+x & y \\ y & 1-x \end{pmatrix} \succeq 0 \}$ 的几何形状。**
   ??? success "参考答案"
       PSD 条件要求 $1+x \ge 0, 1-x \ge 0$ 且行列式满足 $(1+x)(1-x) - y^2 \ge 0 \implies 1 - x^2 - y^2 \ge 0$。这是 $\mathbb{R}^2$ 中的单位圆盘。

7. **[对偶性] SDP 的强对偶性通常需要什么条件？**
   ??? success "参考答案"
       通常需要严格可行性（Slater 条件）：即必须存在一个位于锥内部的可行点（即 $X \succ 0$ 或对偶松弛变量 $Z \succ 0$）。

8. **[自对偶] 证明 $S_n^+$ 在迹内积下是自对偶的。**
   ??? success "参考答案"
       对偶锥为 $K^* = \{Y : \operatorname{tr}(XY) \ge 0, \forall X \succeq 0\}$。取 $X=vv^T$ 可知 $v^T Y v \ge 0$，故 $Y \succeq 0$。反之，两个 PSD 矩阵乘积的迹总是非负的，故 $S_n^+ \subseteq K^*$。

9. **[极小极大] Von Neumann 极小极大定理为矩阵博弈保证了什么？**
   ??? success "参考答案"
       它保证了对于任何二人零和博弈，都存在一个值 $V$ 和一组混合策略，使得任何一方都无法仅通过改变自己的策略来改善期望结果。

10. **[复杂度] 对比在 PSD 锥与共正 (Copositive) 锥中检查成员身份的复杂度。**
    ??? success "参考答案"
        检查 $A \succeq 0$ 是多项式时间问题（通过特征值）。而检查 $A$ 是否共正（$x^T Ax \ge 0, \forall x \ge 0$）是 NP-困难的，这展示了在向量上增加非负性约束所引入的计算困难性。

## 本章小结

本章通过凸性的视角审视矩阵空间的几何：

1. **锥支配**：确立了 PSD 锥作为矩阵分析中主要正则锥的地位，详述了其自对偶性与极端射线。
2. **组合凸性**：通过 Birkhoff 定理建立了双随机矩阵与置换的联系，桥接了离散与连续对称性。
3. **谱面体几何**：定义了 LMI 的可行域，展示了谱面体作为现代优化中广义“多面体”的角色。
4. **SDP 基础**：阐述了半定规划的对偶理论，为解决矩阵值约束提供了分析框架。
