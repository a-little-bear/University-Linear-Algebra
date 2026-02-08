# 线性代数全览

<div style="text-align: center; margin: 2em 0;">
<strong style="font-size: 1.3em;">从大一到博士——全面系统的线性代数知识库</strong>
</div>

---

## 关于本站

本站是一个**全面、系统、自包含**的线性代数知识库，涵盖从本科大一基础课程到博士研究生阶段的所有核心线性代数内容。无论你是初学者还是研究者，都可以在这里找到所需的知识。

本站内容按照难度递进、逻辑连贯的方式组织为四个部分，共 25 章。每章包含完整的定义、定理、证明和例题，力求做到内容严谨、叙述清晰。

---

## 内容导航

### 第一部分：基础线性代数 <span class="difficulty-tag beginner">本科基础</span>

适合大一、大二学生，涵盖线性代数入门课程的全部核心内容。

| 章节 | 内容概要 |
|------|---------|
| [第 1 章 线性方程组](part1/ch01-linear-equations.md) | 线性方程组的求解、高斯消元法、解的结构 |
| [第 2 章 矩阵与矩阵运算](part1/ch02-matrices.md) | 矩阵运算、逆矩阵、分块矩阵、初等矩阵、矩阵的秩 |
| [第 3 章 行列式](part1/ch03-determinants.md) | 行列式的定义与性质、展开定理、Cramer 法则 |
| [第 4 章 向量空间](part1/ch04-vector-spaces.md) | 向量空间公理、子空间、基与维数、秩-零化度定理 |
| [第 5 章 线性变换](part1/ch05-linear-transformations.md) | 线性映射、核与像、矩阵表示、基变换 |
| [第 6 章 特征值与特征向量](part1/ch06-eigenvalues.md) | 特征多项式、对角化、Cayley-Hamilton 定理 |
| [第 7 章 正交性与最小二乘](part1/ch07-orthogonality.md) | 正交集、Gram-Schmidt、正交投影、最小二乘法 |

### 第二部分：中级线性代数 <span class="difficulty-tag intermediate">本科进阶</span>

适合大二至大四学生以及研究生初级阶段，深入探讨线性代数的核心理论。

| 章节 | 内容概要 |
|------|---------|
| [第 8 章 内积空间](part2/ch08-inner-product-spaces.md) | 一般内积空间、正交补、伴随算子、谱定理 |
| [第 9 章 二次型](part2/ch09-quadratic-forms.md) | 二次型的标准形、惯性定理、正定性判别 |
| [第 10 章 矩阵分解](part2/ch10-matrix-decompositions.md) | LU、Cholesky、QR、Schur 分解 |
| [第 11 章 奇异值分解](part2/ch11-svd.md) | SVD 的理论与应用、低秩近似、伪逆 |
| [第 12 章 Jordan 标准形](part2/ch12-jordan-form.md) | 广义特征向量、Jordan 块、最小多项式 |
| [第 13 章 矩阵函数](part2/ch13-matrix-functions.md) | 矩阵指数、矩阵对数、矩阵幂级数 |

### 第三部分：高级线性代数 <span class="difficulty-tag advanced">研究生</span>

适合硕士及博士研究生，涵盖矩阵分析与高级理论。

| 章节 | 内容概要 |
|------|---------|
| [第 14 章 矩阵分析](part3/ch14-matrix-analysis.md) | 矩阵序列与级数、谱半径、Gershgorin 定理 |
| [第 15 章 范数与扰动理论](part3/ch15-norms-perturbation.md) | 矩阵范数、条件数、特征值扰动 |
| [第 16 章 正定矩阵](part3/ch16-positive-definite.md) | 正定矩阵的等价条件、Schur 补、Löwner 偏序 |
| [第 17 章 非负矩阵与 Perron-Frobenius 理论](part3/ch17-nonnegative-matrices.md) | Perron-Frobenius 定理、不可约矩阵、随机矩阵 |
| [第 18 章 矩阵不等式](part3/ch18-matrix-inequalities.md) | 特征值不等式、迹不等式、行列式不等式、Majorization |
| [第 19 章 Kronecker 积与 Vec 算子](part3/ch19-kronecker.md) | Kronecker 积、Vec 算子及其在矩阵方程中的应用 |
| [第 20 章 矩阵方程](part3/ch20-matrix-equations.md) | Sylvester 方程、Lyapunov 方程、Riccati 方程 |

### 第四部分：专题研究 <span class="difficulty-tag research">博士</span>

面向博士研究生和研究者，介绍线性代数的前沿专题。

| 章节 | 内容概要 |
|------|---------|
| [第 21 章 多线性代数与张量](part4/ch21-multilinear-algebra.md) | 对偶空间、张量积、外代数、张量分解 |
| [第 22 章 数值线性代数](part4/ch22-numerical-linear-algebra.md) | 迭代法、Krylov 子空间、数值稳定性 |
| [第 23 章 随机矩阵初步](part4/ch23-random-matrices.md) | Wigner 半圆律、Marchenko-Pastur 律、特征值分布 |
| [第 24 章 矩阵流形](part4/ch24-matrix-manifolds.md) | Stiefel 流形、Grassmann 流形、矩阵 Lie 群 |
| [第 25 章 线性代数在优化中的应用](part4/ch25-optimization.md) | 半定规划、矩阵补全、压缩感知、PCA |

---

## 符号约定

本站统一使用以下数学符号约定：

| 符号 | 含义 |
|------|------|
| $\mathbb{R}, \mathbb{C}, \mathbb{F}$ | 实数域、复数域、一般数域 |
| $\mathbb{R}^n, \mathbb{R}^{m \times n}$ | $n$ 维实向量空间、$m \times n$ 实矩阵空间 |
| $A, B, C$ | 矩阵（大写字母） |
| $\mathbf{v}, \mathbf{u}, \mathbf{w}$ | 向量（粗体小写字母） |
| $a, b, \lambda, \alpha$ | 标量（小写字母或希腊字母） |
| $V, W, U$ | 向量空间 |
| $T, S$ | 线性变换 |
| $A^T, A^H$ | 转置、共轭转置（Hermitian 转置） |
| $A^{-1}$ | 逆矩阵 |
| $\det(A)$ | 行列式 |
| $\operatorname{tr}(A)$ | 迹 |
| $\operatorname{rank}(A)$ | 秩 |
| $\dim(V)$ | 维数 |
| $\ker(T), \operatorname{im}(T)$ | 核（零空间）、像（值域） |
| $\langle \mathbf{u}, \mathbf{v} \rangle$ | 内积 |
| $\|\mathbf{v}\|$ | 范数 |
| $\sigma_i(A)$ | 第 $i$ 个奇异值 |
| $\lambda_i(A)$ | 第 $i$ 个特征值 |
| $I_n$ 或 $I$ | $n$ 阶单位矩阵 |
| $O$ | 零矩阵 |
| $\mathbf{0}$ | 零向量 |
| $\oplus$ | 直和 |
| $\otimes$ | 张量积 / Kronecker 积 |

---

## 如何使用本站

- **搜索**：使用页面顶部的搜索栏，可以快速查找任何概念、定理或关键词。支持中英文搜索。
- **导航**：通过左侧导航栏按章节浏览，或使用右侧目录快速定位到当前页面的某一节。
- **深浅模式**：点击顶部的亮度图标可在浅色和深色模式之间切换。
- **数学公式**：本站所有公式使用 MathJax 渲染，支持复制公式的 LaTeX 源码（右键点击公式）。

---

## 参考文献

本站内容参考了以下经典教材：

1. **Strang, G.** *Introduction to Linear Algebra*. Wellesley-Cambridge Press.
2. **Lay, D.C.** *Linear Algebra and Its Applications*. Pearson.
3. **Axler, S.** *Linear Algebra Done Right*. Springer.
4. **Hoffman, K. & Kunze, R.** *Linear Algebra*. Prentice Hall.
5. **Horn, R.A. & Johnson, C.R.** *Matrix Analysis*. Cambridge University Press.
6. **Horn, R.A. & Johnson, C.R.** *Topics in Matrix Analysis*. Cambridge University Press.
7. **Meyer, C.D.** *Matrix Analysis and Applied Linear Algebra*. SIAM.
8. **Golub, G.H. & Van Loan, C.F.** *Matrix Computations*. Johns Hopkins University Press.
9. **Halmos, P.R.** *Finite-Dimensional Vector Spaces*. Springer.
10. **Lax, P.D.** *Linear Algebra and Its Applications*. Wiley.
