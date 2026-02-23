# 第 28 章 线性代数在量子信息中的应用

<div class="context-flow" markdown>

**前置**：内积空间(Ch8) · Kronecker 积(Ch19) · 正规算子(Ch8) · SVD(Ch11)

**本章脉络**：量子态与 Hilbert 空间 → 算符与可观测量（Hermitian 算子） → 算符的谱分解 → 测量与投影算子 → 态叠加与张量积 → 量子纠缠与 Schmidt 分解 → 密度矩阵与迹类算子 → 应用（量子逻辑门、量子搜索算法、贝尔不等式）

**延伸**：量子力学是建立在复内积空间之上的线性代数；量子信息的每一条公理本质上都是对线性算子性质的物理描述

</div>

量子信息科学将经典逻辑位（0 或 1）扩展为量子位（态空间的单位向量）。在这个微观世界中，所有的物理过程——从态的演化到信息的读取——都严格遵循复线性空间上的线性算子理论。

---

## 28.1 核心概念与公理

!!! definition "定义 28.1 (量子态与算子)"
    - 一个孤立量子系统的态由 Hilbert 空间中的单位向量 $|\psiangle$ 表示。
    - 物理可观测量由 **Hermitian 算子** $H$（满足 $H = H^*$）表示。

!!! theorem "定理 28.3 (Schmidt 分解)"
    任何复合系统 $AB$ 的态 $|\Psiangle$ 均可分解为 $|\Psiangle = \sum \sqrt{\lambda_i} |u_iangle |v_iangle$，其中 $\lambda_i$ 是通过部分迹得到的密度矩阵的特征值。

---

## 练习题

1. **[基础] 证明两个量子态的叠加态 $|\psiangle = a|0angle + b|1angle$ 是单位向量的条件。**
   ??? success "参考答案"
       $\langle \psi | \psi angle = 1$。
       展开得 $(a^* \langle 0| + b^* \langle 1|)(a|0angle + b|1angle) = |a|^2 \langle 0|0angle + |b|^2 \langle 1|1angle = |a|^2 + |b|^2 = 1$。这就是概率归一化条件。

2. **[Hermitian算子] 若 $H$ 代表能量（Hamiltonian），证明其特征值必为实数。**
   ??? success "参考答案"
       这是 Hermitian 算子的基本性质。设 $H|vangle = E|vangle$，则 $\langle v|H|vangle = E \langle v|vangle$。
       取伴随得 $\langle v|H^*|vangle = E^* \langle v|vangle$。
       由于 $H^* = H$，故 $E = E^*$，即 $E$ 为实数。这保证了物理测量值（能量）是可观测的实数。

3. **[量子门] 验证 Hadamard 门 $H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix}$ 是酉矩阵。**
   ??? success "参考答案"
       计算 $HH^*$：$\frac{1}{2} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} = \frac{1}{2} \begin{pmatrix} 2 & 0 \ 0 & 2 \end{pmatrix} = I$。
       酉性保证了量子演化过程是保持概率幅长度的，即信息不丢失。

4. **[张量积] 计算复合态 $|+angle \otimes |-angle$，其中 $|+angle = \frac{1}{\sqrt{2}}(|0angle+|1angle), |-angle = \frac{1}{\sqrt{2}}(|0angle-|1angle)$。**
   ??? success "参考答案"
       $|+angle \otimes |-angle = \frac{1}{2} \begin{pmatrix} 1 \ 1 \end{pmatrix} \otimes \begin{pmatrix} 1 \ -1 \end{pmatrix} = \frac{1}{2} [1, -1, 1, -1]^T$。
       展开为 $|00angle - |01angle + |10angle - |11angle$。

5. **[量子纠缠] 判定贝尔态 $|\Phi^+angle = \frac{1}{\sqrt{2}}(|00angle + |11angle)$ 是否可分离。**
   ??? success "参考答案"
       不可分离。若可分离，则 $|\Phi^+angle = (a|0angle+b|1angle) \otimes (c|0angle+d|1angle) = ac|00angle + ad|01angle + bc|10angle + bd|11angle$。
       由系数得 $ad=0$ 且 $bc=0$。但这会推导出 $ac=0$ 或 $bd=0$，矛盾。故该态是纠缠态。

6. **[测量] 给定态 $|\psiangle = \frac{\sqrt{3}}{2}|0angle + \frac{1}{2}|1angle$，测量基为 $\{|0angle, |1angle\}$。求得到 0 的概率。**
   ??? success "参考答案"
       $P(0) = |\langle 0|\psiangle|^2 = | \frac{\sqrt{3}}{2} \langle 0|0angle + \frac{1}{2} \langle 0|1angle |^2 = (\frac{\sqrt{3}}{2})^2 = 3/4$。

7. **[密度矩阵] 纯态 $|\psiangle$ 的密度矩阵 $ho = |\psiangle\langle \psi|$ 满足什么幂等性？**
   ??? success "参考答案"
       $ho^2 = (|\psiangle\langle \psi|)(|\psiangle\langle \psi|) = |\psiangle (\langle \psi|\psiangle) \langle \psi| = |\psiangle \langle \psi| = ho$。
       纯态密度矩阵是幂等的投影算子，且 $\operatorname{tr}(ho) = 1$。

8. **[部分迹] 为什么部分迹运算（Partial Trace）对应于“忽略一个子系统”？**
   ??? success "参考答案"
       部分迹 $ho_A = \operatorname{tr}_B(ho_{AB})$ 是一个线性映射，它将全系统的密度矩阵映射到子系统的密度矩阵。代数上，它通过对子系统 B 的基向量求内积和，消除了 B 的所有信息，仅保留 A 的边际概率分布。

9. **[Pauli矩阵] 计算 Pauli-X 和 Pauli-Z 矩阵的对易子 $[X, Z]$。**
   ??? success "参考答案"
       $X = \begin{pmatrix} 0 & 1 \ 1 & 0 \end{pmatrix}, Z = \begin{pmatrix} 1 & 0 \ 0 & -1 \end{pmatrix}$。
       $XZ = \begin{pmatrix} 0 & -1 \ 1 & 0 \end{pmatrix}, ZX = \begin{pmatrix} 0 & 1 \ -1 & 0 \end{pmatrix}$。
       $[X, Z] = XZ - ZX = \begin{pmatrix} 0 & -2 \ 2 & 0 \end{pmatrix} = -2iY$。非零对易子体现了量子力学中的不确定性。

10. **[Schmidt秩] Schmidt 分解中非零系数的数量代表什么？**
    ??? success "参考答案"
        非零系数的数量称为 **Schmidt 秩**。若 Schmidt 秩 $> 1$，则该复合系统处于纠缠态。Schmidt 秩是量化纠缠程度最基础的代数指标。

## 本章小结

量子信息是线性代数最前沿的应用场景：

1. **叠加与张量**：量子位通过线性叠加占据空间，复合系统通过张量积爆炸式扩张。
2. **算子即逻辑**：量子计算的所有门操作都是复空间上的酉旋转。
3. **特征值即结果**：所有的物理测量本质上都是在投影算子定义的特征空间上提取特征值及其概率。
