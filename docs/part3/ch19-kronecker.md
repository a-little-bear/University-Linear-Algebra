# 第 19 章 Kronecker 积与 Vec 算子

<div class="context-flow" markdown>

**前置**：矩阵运算 (Ch02) · 矩阵方程初步 (Ch20) · 张量积初步 (Ch21)

**本章脉络**：Kronecker 积的定义 $\to$ 基本代数性质（乘积、转置、逆） $\to$ 谱性质（特征值与迹） $\to$ Vec 算子及其性质 $\to$ 核心恒等式 $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ $\to$ Kronecker 和 ($A \oplus B$) $\to$ 应用：Sylvester 与 Lyapunov 矩阵方程的向量化求解 $\to$ 信号处理中的张量结构

**延伸**：Kronecker 积是张量积 (Tensor Product) 在矩阵维度的具体实现；它是处理多变量系统、多维信号以及量子纠缠态 (Ch28) 的基础代数语言

</div>

当我们处理涉及多个矩阵相互作用的复杂方程（如 $AX + XB = C$）时，传统的矩阵算术往往难以直接给出闭式解。Kronecker 积与 Vec 算子提供了一种将“矩阵方程”转化为“向量方程”的标准工具，将高维的算子相互作用降维打击为经典的线性系统。

---

## 19.1 Kronecker 积

!!! definition "定义 19.1 (Kronecker 积)"
    设 $A$ 为 $m \times n$ 矩阵，$B$ 为 $p \times q$ 矩阵。它们的 **Kronecker 积** $A \otimes B$ 是一个 $mp \times nq$ 的分块矩阵：
    $$A \otimes B = \begin{pmatrix} a_{11}B & \cdots & a_{1n}B \\ \vdots & \ddots & \vdots \\ a_{m1}B & \cdots & a_{mn}B \end{pmatrix}$$

!!! theorem "定理 19.1 (核心性质)"
    1.  **混合乘积性质**：$(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$。
    2.  **谱性质**：若 $A$ 的特征值为 $\{\lambda_i\}$，$B$ 的特征值为 $\{\mu_j\}$，则 $A \otimes B$ 的特征值为 $\{\lambda_i \mu_j\}$。
    3.  **迹与行列式**：
        - $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$。
        - $\det(A \otimes B) = (\det A)^p (\det B)^m$。

---

## 19.2 Vec 算子与向量化

!!! definition "定义 19.2 (Vec 算子)"
    对于 $m \times n$ 矩阵 $X$，$\operatorname{vec}(X)$ 是将 $X$ 的列按顺序堆叠成的一个 $mn \times 1$ 向量。

!!! theorem "定理 19.2 (矩阵方程向量化恒等式)"
    对于适当维数的矩阵 $A, X, B$：
    $$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$$
    **意义**：这是本章最重要的公式，它将矩阵变量 $X$ 提取出来，使方程化为 $My = f$ 的标准形式。

---

## 19.3 Kronecker 和

!!! definition "定义 19.3 (Kronecker 和)"
    对于 $n$ 阶方阵 $A$ 和 $m$ 阶方阵 $B$：
    $$A \oplus B = A \otimes I_m + I_n \otimes B$$
    **性质**：$A \oplus B$ 的特征值为 $\{\lambda_i + \mu_j\}$。
    **应用**：求解 Sylvester 方程 $AX + XB = C$ 等价于求解系统 $(I \otimes A + B^T \otimes I) \operatorname{vec}(X) = \operatorname{vec}(C)$。

---

## 练习题

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

****

??? success "参考答案"
    

## 本章小结

Kronecker 积与 Vec 算子提供了矩阵代数的一套“升维与降维”方案：

1.  **维度的乘法**：Kronecker 积通过一种平铺嵌套的方式，将两个独立算子的作用整合为一个巨大的复合算子，是处理多体系统互作用的唯一语言。
2.  **算子的解构**：Vec 算子消除了矩阵的二维拓扑，将其转化为线性空间中最基本的向量形式，从而释放了经典线性方程组求解器的全部威力。
3.  **方程的统一**：矩阵方程向量化恒等式是连接矩阵方程理论与数值线性代数的纽带，它证明了所有的线性矩阵方程在本质上都是同一个线性系统在不同基下的表现。
