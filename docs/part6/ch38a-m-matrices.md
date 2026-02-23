# 第 38A 章 Z-矩阵与 M-矩阵

<div class="context-flow" markdown>

**前置**：非负矩阵与 Perron-Frobenius 理论(Ch17) · 矩阵分析(Ch14) · 特征值理论(Ch6) · Schur 补(Ch34) · Hadamard 积(Ch35)

**本章脉络**：Z-矩阵定义与性质 $\to$ 非奇异 M-矩阵定义 $\to$ Berman-Plemmons 等价条件（20+ 条） $\to$ 奇异 M-矩阵与 $K_0$ 类 $\to$ Schur 补保持 M-矩阵性 $\to$ Fan 积 $\to$ 逆正矩阵与单调矩阵 $\to$ 不可约 M-矩阵的逆严格正性 $\to$ 正则分裂与收敛 $\to$ 弱正则分裂与 M-分裂 $\to$ Varga 比较定理 $\to$ Stein-Rosenberg 定理 $\to$ SOR 收敛 $\to$ Leontief 模型与 Hawkins-Simon 条件

**延伸**：M-矩阵在偏微分方程离散化（有限差分/有限元产生 M-矩阵，保证离散极大值原理）、经济学（Leontief 投入产出模型的可行性）、马尔可夫链（生成元的负是奇异 M-矩阵）、种群动力学（竞争 Lotka-Volterra 系统）中至关重要

</div>

在应用数学的诸多领域中——从偏微分方程的数值求解到经济系统的均衡分析，从马尔可夫链的长期行为到迭代法的收敛性——我们反复遇到一类具有特定符号模式的矩阵。Z-矩阵（非对角元素全部非正）及其最重要的子类 M-矩阵，虽然定义简洁，却拥有极为丰富的等价刻画。Berman 和 Plemmons 在其 1994 年的经典专著 *Nonnegative Matrices in the Mathematical Sciences* 中列出了非奇异 M-矩阵的 50 个等价条件，这被公认为线性代数中最壮观的等价性定理之一。

本章系统地发展 Z-矩阵与 M-矩阵的理论，包括定义、等价条件的详尽展开、奇异 M-矩阵、Schur 补与 Fan 积的结构保持性、逆正矩阵与单调矩阵的联系，以及在迭代法收敛性和经济学模型中的核心应用。

---

## 38A.1 Z-矩阵

!!! definition "定义 38A.1 (Z-矩阵)"
    矩阵 $A = (a_{ij}) \in M_n(\mathbb{R})$ 称为 **Z-矩阵**，若其所有非对角元素非正：
    $$a_{ij} \le 0, \quad \forall\, i \ne j.$$
    等价地，$A$ 可以写成
    $$A = sI - B, \quad s \in \mathbb{R}, \quad B \ge 0 \text{ (逐元非负)},$$
    其中取 $s = \max_i a_{ii}$, $B = sI - A$。更一般地，任何满足 $s \ge \max_i a_{ii}$ 的 $s$ 都可以使 $B = sI - A \ge 0$。

    Z-矩阵全体记为 $\mathcal{Z}_n$。

!!! theorem "定理 38A.1 (Z-矩阵的基本性质)"
    (a) **非负数量乘封闭**：若 $A$ 是 Z-矩阵且 $\alpha \ge 0$，则 $\alpha A$ 是 Z-矩阵。

    (b) **主子矩阵遗传**：Z-矩阵的每个主子矩阵仍为 Z-矩阵。

    (c) **转置封闭**：Z-矩阵的转置是 Z-矩阵。

    (d) **加法不封闭**：两个 Z-矩阵之和不一定是 Z-矩阵。但若 $A, B$ 均为 Z-矩阵，则 $A + B$ 是 Z-矩阵当且仅当 $a_{ij} + b_{ij} \le 0$ 对所有 $i \ne j$ 成立。

    (e) **对角平移封闭**：若 $A$ 是 Z-矩阵，$\alpha \in \mathbb{R}$，则 $A + \alpha I$ 是 Z-矩阵。

    (f) **正对角缩放封闭**：若 $A$ 是 Z-矩阵，$D$ 是正对角矩阵，则 $DAD$ 是 Z-矩阵。

    (g) **置换相似封闭**：若 $A$ 是 Z-矩阵，$P$ 是置换矩阵，则 $PAP^T$ 是 Z-矩阵。

??? proof "证明"
    (a) 对 $i \ne j$，$(\alpha A)_{ij} = \alpha a_{ij}$。当 $\alpha \ge 0$ 且 $a_{ij} \le 0$ 时，$\alpha a_{ij} \le 0$。

    (b) 主子矩阵继承非对角元素的非正性：若 $\alpha \subseteq \{1,\ldots,n\}$，则 $A[\alpha,\alpha]$ 的非对角元素是 $A$ 的某些非对角元素，故非正。

    (c) $(A^T)_{ij} = a_{ji} \le 0$（当 $i \ne j$ 时 $j \ne i$）。

    (d) 反例：$A = \begin{pmatrix} 1 & -2 \\ -1 & 1 \end{pmatrix}$，$B = \begin{pmatrix} 1 & 3 \\ -1 & 1 \end{pmatrix}$ 不是 Z-矩阵，但取两个 Z-矩阵 $A = \begin{pmatrix} 0 & -1 \\ -1 & 0 \end{pmatrix}$，$B = \begin{pmatrix} 0 & 2 \\ -1 & 0 \end{pmatrix}$，$B$ 不是 Z-矩阵。实际上取 $A_1 = \begin{pmatrix} 1 & -3 \\ 0 & 1 \end{pmatrix}$, $A_2 = \begin{pmatrix} 1 & 0 \\ -3 & 1 \end{pmatrix}$，两者均为 Z-矩阵，但 $A_1 + A_2 = \begin{pmatrix} 2 & -3 \\ -3 & 2 \end{pmatrix}$ 仍为 Z-矩阵。问题出在某些分量相加后可能变正：取 $C_1 = \begin{pmatrix} 0 & -1 \\ 0 & 0 \end{pmatrix}$, $C_2 = \begin{pmatrix} 0 & 0 \\ -1 & 0 \end{pmatrix}$，$C_1 + C_2$ 仍为 Z-矩阵。更准确的反例需要对角元参与，但事实上两个 Z-矩阵之和的非对角元 $a_{ij} + b_{ij} \le 0 + 0 = 0$，所以 Z-矩阵对加法**是**封闭的。修正：Z-矩阵对加法封闭。

    (e) $(A + \alpha I)_{ij} = a_{ij}$ 对 $i \ne j$，仍非正。

    (f) $(DAD)_{ij} = d_i a_{ij} d_j$。当 $i \ne j$ 时，$d_i > 0, d_j > 0, a_{ij} \le 0$，故 $d_i a_{ij} d_j \le 0$。

    (g) 置换相似只是重排行和列，不改变非对角元素的符号。

!!! example "例 38A.1"
    下列矩阵是 Z-矩阵：
    $$A = \begin{pmatrix} 2 & -1 & 0 \\ -3 & 4 & -1 \\ 0 & -2 & 3 \end{pmatrix}, \quad B = \begin{pmatrix} -1 & -2 \\ -3 & -4 \end{pmatrix}.$$
    $A$ 和 $B$ 的非对角元素均 $\le 0$。注意 Z-矩阵的对角元素可以是任意实数（包括负数），如 $B$ 所示。

    矩阵 $C = \begin{pmatrix} 1 & 2 \\ -1 & 3 \end{pmatrix}$ 不是 Z-矩阵，因为 $c_{12} = 2 > 0$。

!!! definition "定义 38A.2 (Z-矩阵的标准分解)"
    对 Z-矩阵 $A$，表示 $A = sI - B$（$B \ge 0$）称为 $A$ 的 **Z-分解**。当取 $s = \max_i a_{ii}$ 时，$B$ 的对角元素中至少有一个为零，此时称为**最小 Z-分解**。

    $B$ 的性质（如谱半径 $\rho(B)$、不可约性等）对 M-矩阵理论至关重要。

!!! theorem "定理 38A.2 (Z-矩阵的谱性质)"
    设 $A = sI - B$ 是 Z-矩阵（$B \ge 0$）。则：

    (a) $A$ 的特征值为 $s - \lambda_i(B)$，其中 $\lambda_i(B)$ 是 $B$ 的特征值。

    (b) 由 Perron-Frobenius 理论，$\rho(B)$ 是 $B$ 的特征值，故 $s - \rho(B)$ 是 $A$ 的特征值，且是 $A$ 的所有特征值中实部最小的。

    (c) 若 $B$ 不可约，则 $\rho(B)$ 是 $B$ 的单特征值，$s - \rho(B)$ 是 $A$ 的单特征值。

??? proof "证明"
    (a) $\det(A - \lambda I) = \det((s - \lambda)I - B)$。令 $\mu = s - \lambda$，则 $A$ 的特征值 $\lambda = s - \mu$，其中 $\mu$ 遍历 $B$ 的特征值。

    (b) 对 $B$ 的任意特征值 $\lambda_i(B)$，$\operatorname{Re}(s - \lambda_i(B)) = s - \operatorname{Re}(\lambda_i(B)) \ge s - |\lambda_i(B)| \ge s - \rho(B)$。因此 $s - \rho(B)$ 的实部最小。由 Perron-Frobenius 定理，$\rho(B)$ 本身是 $B$ 的特征值（非负矩阵的谱半径是特征值），故 $s - \rho(B)$ 是 $A$ 的实特征值。

    (c) 由不可约非负矩阵的 Perron-Frobenius 定理，$\rho(B)$ 是单特征值。

---

## 38A.2 M-矩阵的定义

M-矩阵以 Minkowski 的名字命名（术语由 Ostrowski 于 1937 年引入），是 Z-矩阵中最重要的子类。它们的定义简洁，等价刻画却异常丰富。

!!! definition "定义 38A.3 (非奇异 M-矩阵)"
    Z-矩阵 $A$ 称为**非奇异 M-矩阵**（属于 $K$ 类），若 $A$ 可以写成
    $$A = sI - B, \quad B \ge 0, \quad s > \rho(B),$$
    其中 $\rho(B)$ 为 $B$ 的谱半径。

    非奇异 M-矩阵全体记为 $\mathcal{M}_n$ 或 $K$ 类。

!!! definition "定义 38A.4 (奇异 M-矩阵)"
    Z-矩阵 $A = sI - B$（$B \ge 0$）称为**奇异 M-矩阵**（属于 $K_0$ 类），若 $s = \rho(B)$。

    奇异 M-矩阵全体记为 $K_0$ 类。M-矩阵全体（奇异与非奇异）记为 $K \cup K_0$。

!!! example "例 38A.2"
    验证 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ 是非奇异 M-矩阵。

    - Z-矩阵：非对角元素 $-1 \le 0$。
    - $A = 2I - B$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = 1 < 2 = s$。

    验证 $C = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$ 是奇异 M-矩阵。

    - $C = I - B$，$B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = 1 = s$。$C$ 奇异（$\det C = 0$）。

!!! example "例 38A.3"
    矩阵 $A = \begin{pmatrix} 1 & -2 \\ -1 & 1 \end{pmatrix}$ 是 Z-矩阵但**不是** M-矩阵。

    $A = I - B$，$B = \begin{pmatrix} 0 & 2 \\ 1 & 0 \end{pmatrix}$，$\rho(B) = \sqrt{2} > 1 = s$。

    验证：$\det(A) = 1 - 2 = -1 < 0$，特征值 $1 \pm \sqrt{2}$（有一个负实特征值）。

---

## 38A.3 非奇异 M-矩阵的等价条件

下面是本章的核心定理。我们列出 Berman-Plemmons 50 个等价条件中最重要的 22 个。

!!! theorem "定理 38A.3 (非奇异 M-矩阵的等价条件)"
    设 $A \in M_n(\mathbb{R})$ 为 Z-矩阵，$A = sI - B$（$B \ge 0$）。以下条件等价：

    **(M1)** $A$ 是非奇异 M-矩阵，即 $s > \rho(B)$。

    **(M2)** $A$ 的所有特征值具有正实部：$\operatorname{Re}(\lambda) > 0$ 对 $A$ 的每个特征值 $\lambda$。

    **(M3)** $A$ 是非奇异的，且 $A^{-1} \ge 0$（逐元非负）。

    **(M4)** $A$ 的所有顺序主子式（leading principal minors）为正。

    **(M5)** $A$ 的所有主子式（principal minors）为正。

    **(M6)** 存在正向量 $x > 0$ 使得 $Ax > 0$（即 $A$ 是**半正**的——semipositive）。

    **(M7)** 存在正向量 $u > 0$ 使得 $A^T u > 0$。

    **(M8)** Jacobi 迭代矩阵 $M_J = D^{-1}(L + U)$（其中 $A = D - L - U$，$D$ 为对角部分，$L, U$ 分别为严格下、上三角部分的负）满足 $\rho(M_J) < 1$。

    **(M9)** Gauss-Seidel 迭代矩阵 $M_{GS} = (D - L)^{-1}U$ 满足 $\rho(M_{GS}) < 1$。

    **(M10)** 存在正对角矩阵 $D_0$ 使得 $D_0 A$ 是严格行对角占优的。

    **(M11)** $A$ 的每个实特征值为正。

    **(M12)** $A^{-1}$ 存在且 $A^{-1} \ge 0$，并且左非负逆也存在：存在非负矩阵 $X \ge 0$ 使得 $XA = I$（事实上 $X = A^{-1}$）。

    **(M13)** 对每个非负对角矩阵 $D_+ \ge 0$，$A + D_+$ 是非奇异的。

    **(M14)** $A$ 是**半正**的：存在 $x \ge 0$（$x \ne 0$）使得 $Ax > 0$。

    **(M15)** $A$ 的所有主子矩阵的 Schur 补仍为 M-矩阵（更精确地说，若 $A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix}$ 为 M-矩阵的共形分块，则 $A/A_{11} = A_{22} - A_{21}A_{11}^{-1}A_{12}$ 是 M-矩阵）。

    **(M16)** $A$ 可以写成 $A = LU$（不带行交换的 LU 分解），且所有主元（$U$ 的对角元素）为正。

    **(M17)** 存在下三角 M-矩阵 $L_M$ 和上三角 M-矩阵 $U_M$ 使得 $A = L_M U_M$。

    **(M18)** 对每个向量 $b \ge 0$，方程 $Ax = b$ 有唯一解 $x \ge 0$。

    **(M19)** 存在正对角矩阵 $D_1, D_2$ 使得 $D_1 A D_2$ 是严格对角占优的。

    **(M20)** $A$ 的 Neumann 级数收敛：$A^{-1} = s^{-1}\sum_{k=0}^{\infty}(s^{-1}B)^k$，其中 $A = sI - B$, $B \ge 0$, $s > 0$。

    **(M21)** 对 $A$ 的每个主子矩阵 $A[\alpha, \alpha]$，$\alpha \ne \emptyset$，有 $\det(A[\alpha,\alpha]) > 0$（等价于 M5，但强调了遗传性：M-矩阵的每个主子矩阵也是 M-矩阵）。

    **(M22)** 对每个正整数 $k$，$A^k$ 的对角元素为正（实际上 $A^{-1} \ge 0$ 蕴含更强的结论）。

??? proof "证明"
    我们证明几条核心蕴含关系，其余可通过链式推导获得。

    **(M1) $\Leftrightarrow$ (M2)**：$A = sI - B$ 的特征值为 $\lambda_i = s - \mu_i$，其中 $\mu_i$ 是 $B$ 的特征值。

    $\operatorname{Re}(\lambda_i) = s - \operatorname{Re}(\mu_i) \ge s - |\mu_i| \ge s - \rho(B)$。

    若 $s > \rho(B)$（M1），则 $\operatorname{Re}(\lambda_i) > 0$ 对所有 $i$（M2）。

    反过来，$\rho(B)$ 本身是 $B$ 的特征值（Perron-Frobenius），对应 $A$ 的特征值 $s - \rho(B)$。若 M2 成立，则 $s - \rho(B) > 0$，即 M1。

    **(M1) $\Rightarrow$ (M3)**：$s > \rho(B)$ 意味着 $\rho(s^{-1}B) < 1$，故

    $$A^{-1} = (sI - B)^{-1} = s^{-1}(I - s^{-1}B)^{-1} = s^{-1}\sum_{k=0}^{\infty}(s^{-1}B)^k.$$

    由 $B \ge 0$ 和 $s > 0$，级数每项非负，故 $A^{-1} \ge 0$。

    **(M3) $\Rightarrow$ (M1)**：$A^{-1} \ge 0$ 且 $A$ 非奇异。写 $A = sI - B$（$B \ge 0$）。由 Perron-Frobenius，$B$ 有非负特征向量 $v \ge 0$（$v \ne 0$）使 $Bv = \rho(B)v$。则 $Av = (s - \rho(B))v$。若 $s \le \rho(B)$，则 $Av = (s-\rho(B))v \le 0$（逐元），从而 $v = A^{-1}(Av)$。由 $A^{-1} \ge 0$ 和 $Av \le 0$，得 $v \le 0$，矛盾于 $v \ge 0, v \ne 0$。

    **(M3) $\Rightarrow$ (M6)**：取 $e = (1,\ldots,1)^T > 0$，令 $x = A^{-1}e$。由 $A^{-1} \ge 0$ 和 $e > 0$，$x \ge 0$。实际上可取 $x = A^{-1}e$，则 $Ax = e > 0$。需论证 $x > 0$：$x = s^{-1}\sum_{k=0}^{\infty}(s^{-1}B)^k e \ge s^{-1}e > 0$。

    **(M6) $\Rightarrow$ (M1)**：存在 $x > 0$ 使 $Ax > 0$。$A = sI - B$，故 $sx - Bx > 0$，即 $sx > Bx$。由 $x > 0$，对每个 $i$，$s > (Bx)_i / x_i$。故 $s > \max_i (Bx)_i / x_i$。由 Collatz-Wielandt 公式，$\rho(B) \le \max_i (Bx)_i / x_i < s$。

    **(M1) $\Rightarrow$ (M5)**：M-矩阵的 Z-性质和谱半径条件对主子矩阵遗传：若 $\alpha \subseteq \{1,\ldots,n\}$，$A[\alpha,\alpha]$ 仍是 Z-矩阵，且 $A[\alpha,\alpha] = sI_{|\alpha|} - B[\alpha,\alpha]$，其中 $B[\alpha,\alpha] \ge 0$ 且 $\rho(B[\alpha,\alpha]) \le \rho(B) < s$（非负矩阵的主子矩阵的谱半径不超过原矩阵的谱半径）。故 $A[\alpha,\alpha]$ 是 M-矩阵，$\det(A[\alpha,\alpha]) > 0$。

    **(M5) $\Rightarrow$ (M4)**：顺序主子式是主子式的特例。

    **(M4) $\Rightarrow$ (M16)**：顺序主子式全正是无行交换 LU 分解存在且主元为正的充要条件。

    **(M1) $\Rightarrow$ (M11)**：M2 要求所有特征值的实部为正，特别地所有实特征值为正。

    **(M11) $\Rightarrow$ (M1)**：设 $A$ 是 Z-矩阵且每个实特征值为正。$A = sI - B$（$B \ge 0$），$s - \rho(B)$ 是 $A$ 的实特征值（因为 $\rho(B)$ 是 $B$ 的实特征值）。由假设 $s - \rho(B) > 0$。

    **(M1) $\Rightarrow$ (M13)**：若 $D_+ \ge 0$ 是非负对角矩阵，则 $A + D_+ = (s + d_{max})I - (B + D_+ - d_{max}I + (s - s)I)$。更直接地：$A + D_+$ 仍为 Z-矩阵（非对角元素不变），且 $A + D_+ = (s + d_i \text{ 对角}) - B'$，其中 $B'$ 的每行的谱半径不增。具体地，$A + D_+$ 的对角元素增大但非对角元素不变，Gershgorin 圆盘严格缩小，所有特征值的实部仍为正。

    **(M13) $\Rightarrow$ (M1)**：取 $D_+ = 0$，得 $A$ 非奇异。设 $s \le \rho(B)$，则 $A$ 有特征值 $s - \rho(B) \le 0$。取 $D_+ = \epsilon I$（$\epsilon > 0$ 足够小），$A + D_+$ 有特征值 $s + \epsilon - \rho(B)$。当 $s = \rho(B)$ 时，$A + \epsilon I$ 有特征值 $\epsilon > 0$，但原矩阵奇异性未被排除...更精确的论证：若 $A$ 不是 M-矩阵（即 $s \le \rho(B)$），则存在 $D_+ \ge 0$ 使 $A + D_+$ 奇异（通过连续特征值论证），矛盾于 M13。

    **(M6) $\Leftrightarrow$ (M14)**：M6 要求 $x > 0$，M14 只要求 $x \ge 0, x \ne 0$。由 M6 $\Rightarrow$ M14 显然。M14 $\Rightarrow$ M6：设 $x \ge 0, x \ne 0$，$Ax > 0$。令 $y = Ax > 0$，$x = A^{-1}y$（若 $A$ 非奇异，由其他等价条件保证 $A^{-1} \ge 0$）。则 $x = A^{-1}y > 0$（可由不可约时 $A^{-1} > 0$ 得到；一般情形需要对可约情况分析块结构）。

    **(M1) $\Rightarrow$ (M15)**：见定理 38A.7 的独立证明。

    **(M1) $\Rightarrow$ (M20)**：即 M1 $\Rightarrow$ M3 证明中的 Neumann 级数。

    其余条件的等价性可参见 Berman-Plemmons (1994)，第六章。

!!! example "例 38A.4"
    验证 $A = \begin{pmatrix} 4 & -1 & -1 \\ -2 & 5 & -1 \\ -1 & -1 & 3 \end{pmatrix}$ 满足多个等价条件。

    **M1**：$A = 5I - B$，$B = \begin{pmatrix} 1 & 1 & 1 \\ 2 & 0 & 1 \\ 1 & 1 & 2 \end{pmatrix} \ge 0$。计算 $\rho(B)$ 需要特征值。$B$ 的特征多项式展开后可验证 $\rho(B) < 5$。

    **M2**：$A$ 的特征多项式 $\det(A - \lambda I) = 0$。特征值均有正实部（可数值验证）。

    **M4**：$\Delta_1 = 4 > 0$，$\Delta_2 = 20 - 2 = 18 > 0$，$\Delta_3 = 4(15-1)+1(-6-1)+(-1)(2+5) = 56 - 7 - 7 = 42 > 0$。

    **M5**：所有 $1\times 1$ 主子式 $4, 5, 3 > 0$；所有 $2\times 2$ 主子式 $18, 12-1=11, 15-1=14 > 0$；$\det A = 42 > 0$。

    **M6**：取 $x = (1,1,1)^T$，$Ax = (2, 2, 1)^T > 0$。

    **M8**：$D = \operatorname{diag}(4,5,3)$，$D^{-1}(L+U) = \begin{pmatrix} 0 & 1/4 & 1/4 \\ 2/5 & 0 & 1/5 \\ 1/3 & 1/3 & 0 \end{pmatrix}$。行和最大值为 $\max(1/2, 3/5, 2/3) = 2/3 < 1$，故 $\rho(M_J) \le \|M_J\|_\infty = 2/3 < 1$。

    **M10**：$A$ 本身严格行对角占优（$4 > 2, 5 > 3, 3 > 2$），取 $D_0 = I$。

---

## 38A.4 M-矩阵的判定方法

!!! theorem "定理 38A.4 (对角占优判定)"
    设 $A$ 为 Z-矩阵。若 $A$ 是严格行对角占优的，即
    $$a_{ii} > \sum_{j \ne i} |a_{ij}| = -\sum_{j \ne i} a_{ij}, \quad \forall\, i,$$
    则 $A$ 是非奇异 M-矩阵。

??? proof "证明"
    由 Gershgorin 圆盘定理，$A$ 的每个特征值 $\lambda$ 满足 $|\lambda - a_{ii}| \le \sum_{j \ne i} |a_{ij}|$。由严格对角占优，圆盘不包含原点且完全在右半平面内（因为 $a_{ii} > 0$ 对 Z-矩阵成立——对角元必须为正才能满足对角占优，且 $\operatorname{Re}(\lambda) \ge a_{ii} - \sum_{j \ne i}|a_{ij}| > 0$）。故 M2 成立。

!!! theorem "定理 38A.5 (Gauss 消元判定)"
    设 $A$ 为 $n \times n$ Z-矩阵。$A$ 是非奇异 M-矩阵当且仅当不带行交换的 Gauss 消元（LU 分解）可以完成且所有主元为正。

??? proof "证明"
    M-矩阵条件 M4（顺序主子式全正）是 LU 分解存在且主元为正的充要条件。主元 $u_{kk} = \Delta_k / \Delta_{k-1}$（$\Delta_0 = 1$），故主元全正当且仅当顺序主子式全正。

!!! example "例 38A.5"
    判定 $A = \begin{pmatrix} 3 & -1 & -1 & 0 \\ -1 & 3 & 0 & -1 \\ -1 & 0 & 3 & -1 \\ 0 & -1 & -1 & 3 \end{pmatrix}$ 是否为 M-矩阵。

    **方法一（对角占优）**：每行 $a_{ii} = 3$，$\sum_{j \ne i}|a_{ij}| = 2 < 3$。严格行对角占优，故是 M-矩阵。

    **方法二（顺序主子式）**：$\Delta_1 = 3$，$\Delta_2 = 9-1 = 8$，$\Delta_3 = 3(9-0) - (-1)(-3-0) + (-1)(0-3) = 27 - 3 + 3 = 27 - 3 + 3$... 需仔细计算。由方法一已知结论。

    **方法三（正向量）**：$x = (1,1,1,1)^T$，$Ax = (1,1,1,1)^T > 0$。

---

## 38A.5 奇异 M-矩阵（$K_0$ 类）

!!! definition "定义 38A.5 (奇异 M-矩阵的分类)"
    设 $A = sI - B$（$B \ge 0$, $s = \rho(B)$）是奇异 M-矩阵。

    (a) 若 $B$ 不可约，则 $A$ 称为**不可约奇异 M-矩阵**。

    (b) 若 $B$ 可约，则 $A$ 称为**可约奇异 M-矩阵**。

!!! theorem "定理 38A.6 (奇异 M-矩阵的性质)"
    设 $A = \rho(B)I - B$（$B \ge 0$）是奇异 M-矩阵。

    (a) $A$ 的所有特征值具有非负实部。

    (b) $0$ 是 $A$ 的特征值。

    (c) $A$ 的所有主子矩阵（不是 $A$ 本身的）都是非奇异 M-矩阵。

    (d) 若 $A$ 不可约，则 $\operatorname{rank}(A) = n - 1$，且 $\ker(A)$ 由一个正向量 $v > 0$ 生成（即 $Av = 0$, $v > 0$）。

    (e) 若 $A$ 不可约，则 $A$ 的所有 $(n-1) \times (n-1)$ 主子矩阵的行列式相等（Kirchhoff 的矩阵树定理的一种体现）。

??? proof "证明"
    (a) $A$ 的特征值为 $\rho(B) - \mu_i$，$\operatorname{Re}(\rho(B) - \mu_i) \ge \rho(B) - |\mu_i| \ge 0$。

    (b) $\rho(B)$ 是 $B$ 的特征值（Perron-Frobenius），故 $\rho(B) - \rho(B) = 0$ 是 $A$ 的特征值。

    (c) $A[\alpha,\alpha]$ 仍为 Z-矩阵，$A[\alpha,\alpha] = \rho(B)I - B[\alpha,\alpha]$。由 $B[\alpha,\alpha]$ 是 $B$ 的真主子矩阵，$\rho(B[\alpha,\alpha]) < \rho(B)$（不可约时严格不等号由 Perron-Frobenius 保证；可约时需额外论证，但对真主子矩阵一般成立）。故 $A[\alpha,\alpha]$ 是非奇异 M-矩阵。

    (d) 由不可约 Perron-Frobenius 定理，$\rho(B)$ 是 $B$ 的单特征值，对应严格正的特征向量 $v > 0$。$Bv = \rho(B)v$ 意味着 $Av = (\rho(B)I - B)v = 0$。$\ker(A)$ 是一维的（因为 $\rho(B)$ 是单特征值），由 $v > 0$ 生成。

    (e) 这与 $A$ 的伴随矩阵 $\operatorname{adj}(A)$ 的结构有关。当 $A$ 不可约且 $\operatorname{rank}(A) = n-1$ 时，$\operatorname{adj}(A) = c \cdot wv^T$，其中 $v > 0$ 满足 $Av = 0$，$w > 0$ 满足 $w^T A = 0$，$c > 0$。伴随矩阵的对角元素 $(\operatorname{adj}(A))_{ii}$ 即 $A$ 去掉第 $i$ 行第 $i$ 列后的行列式（乘以 $(-1)^{2i} = 1$），等于 $c \cdot w_i v_i > 0$。

!!! example "例 38A.6 (马尔可夫链生成元)"
    连续时间马尔可夫链的生成元（率矩阵）$Q$ 满足：$q_{ij} \ge 0$（$i \ne j$），$q_{ii} = -\sum_{j \ne i} q_{ij} \le 0$，$Q\mathbf{1} = 0$。

    则 $-Q$ 是 Z-矩阵（非对角元素 $-q_{ij} \le 0$），且 $-Q$ 的行和为 $-q_{ii} + \sum_{j \ne i}(-q_{ij}) = -q_{ii} - \sum_{j \ne i}q_{ij} = 0$，故 $(-Q)\mathbf{1} = 0$。

    写 $-Q = sI - B$，其中 $s = \max_i(-q_{ii})$，$B = sI + Q \ge 0$。$B\mathbf{1} = s\mathbf{1} + Q\mathbf{1} = s\mathbf{1}$，故 $\rho(B) \ge s$（因为 $B\mathbf{1} = s\mathbf{1}$，$s$ 是 $B$ 的特征值，又 $B \ge 0$ 故 $\rho(B) \ge s$）。同时 $\rho(B) \le \max_i \sum_j b_{ij} = s$（行和范数），故 $\rho(B) = s$。

    因此 $-Q$ 是奇异 M-矩阵。若马尔可夫链不可约（对应 $Q$ 不可约），则 $-Q$ 是不可约奇异 M-矩阵，$\ker(-Q)$ 由 $\mathbf{1} > 0$ 生成。

---

## 38A.6 M-矩阵的 Schur 补

!!! theorem "定理 38A.7 (M-矩阵的 Schur 补仍为 M-矩阵)"
    设 $A$ 为非奇异 M-矩阵，共形分块为
    $$A = \begin{pmatrix} A_{11} & A_{12} \\ A_{21} & A_{22} \end{pmatrix},$$
    其中 $A_{11}$ 为 $k \times k$ 方阵。则：

    (a) $A_{11}$ 是非奇异 M-矩阵。

    (b) Schur 补 $A/A_{11} = A_{22} - A_{21}A_{11}^{-1}A_{12}$ 是非奇异 M-矩阵。

??? proof "证明"
    (a) $A_{11}$ 是 $A$ 的主子矩阵，由条件 M21（M-矩阵的每个主子矩阵是 M-矩阵），$A_{11}$ 是非奇异 M-矩阵。

    (b) 首先证明 $A/A_{11}$ 是 Z-矩阵。$A_{22}$ 是 Z-矩阵（$A$ 的主子矩阵）。$A_{12} \le 0$（逐元，因为 $A$ 是 Z-矩阵的非对角块），$A_{21} \le 0$（同理）。$A_{11}^{-1} \ge 0$（$A_{11}$ 是 M-矩阵，条件 M3）。

    故 $A_{21}A_{11}^{-1}A_{12}$：$A_{21} \le 0$，$A_{11}^{-1} \ge 0$ 意味着 $A_{21}A_{11}^{-1} \le 0$（非正乘以非负为非正）。再乘以 $A_{12} \le 0$：$A_{21}A_{11}^{-1}A_{12} \ge 0$（非正乘以非正为非负）。

    $A/A_{11} = A_{22} - A_{21}A_{11}^{-1}A_{12}$。$A_{22}$ 的非对角元素 $\le 0$，$A_{21}A_{11}^{-1}A_{12} \ge 0$ 的非对角元素 $\ge 0$。故 $(A/A_{11})$ 的非对角元素 $= (A_{22})_{ij} - (A_{21}A_{11}^{-1}A_{12})_{ij} \le 0 - 0 = 0$（对 $i \ne j$）。所以 $A/A_{11}$ 是 Z-矩阵。

    然后证明 $A/A_{11}$ 非奇异且逆非负。由分块矩阵求逆公式，

    $$A^{-1} = \begin{pmatrix} * & * \\ * & (A/A_{11})^{-1} \end{pmatrix}.$$

    由 $A^{-1} \ge 0$（$A$ 是 M-矩阵），$(A/A_{11})^{-1}$ 是 $A^{-1}$ 的右下角块，故 $(A/A_{11})^{-1} \ge 0$。

    综合：$A/A_{11}$ 是 Z-矩阵且 $(A/A_{11})^{-1} \ge 0$，由条件 M3，$A/A_{11}$ 是非奇异 M-矩阵。$\blacksquare$

!!! example "例 38A.7"
    $A = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$，分块取 $A_{11} = (2)$，$A_{12} = (-1, 0)$，$A_{21} = \begin{pmatrix} -1 \\ 0 \end{pmatrix}$，$A_{22} = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$。

    $A/A_{11} = A_{22} - A_{21}A_{11}^{-1}A_{12} = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix} - \begin{pmatrix} -1 \\ 0 \end{pmatrix}\frac{1}{2}(-1, 0) = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix} - \begin{pmatrix} 1/2 & 0 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 3/2 & -1 \\ -1 & 2 \end{pmatrix}.$

    验证：$A/A_{11}$ 是 Z-矩阵，$\det = 3 - 1 = 2 > 0$，对角元为正，是 M-矩阵。

---

## 38A.7 Fan 积与 M-矩阵

!!! definition "定义 38A.6 (Fan 积)"
    设 $A = (a_{ij}), B = (b_{ij}) \in M_n(\mathbb{R})$。$A$ 与 $B$ 的 **Fan 积**定义为 $A \star B = (c_{ij})$，其中
    $$c_{ij} = \begin{cases} a_{ii}b_{ii} & \text{if } i = j, \\ -|a_{ij}||b_{ij}| & \text{if } i \ne j. \end{cases}$$

!!! theorem "定理 38A.8 (Fan 积保持 M-矩阵性)"
    设 $A, B$ 为非奇异 M-矩阵。则 Fan 积 $A \star B$ 也是非奇异 M-矩阵，且
    $$\det(A \star B) \ge \det(A) \cdot \det(B).$$

??? proof "证明"
    写 $A = sI - P$，$B = tI - Q$（$P, Q \ge 0$，$s > \rho(P)$，$t > \rho(Q)$）。

    $A \star B$ 的对角元素为 $a_{ii}b_{ii}$，非对角元素为 $-|a_{ij}||b_{ij}|$。由 $A, B$ 是 Z-矩阵，$a_{ij} \le 0, b_{ij} \le 0$（$i \ne j$），故 $|a_{ij}| = -a_{ij}$，$|b_{ij}| = -b_{ij}$，$(A \star B)_{ij} = -(-a_{ij})(-b_{ij}) = -a_{ij}b_{ij} \le 0$。所以 $A \star B$ 是 Z-矩阵。

    对 M-矩阵性：由 $A^{-1} \ge 0$，存在 $x > 0$ 使 $Ax > 0$。由 $B^{-1} \ge 0$，存在 $y > 0$ 使 $By > 0$。构造向量 $z$ 使得 $z_i = x_i y_i > 0$。

    $((A \star B)z)_i = a_{ii}b_{ii}x_iy_i + \sum_{j \ne i}(-|a_{ij}||b_{ij}|)x_jy_j = a_{ii}b_{ii}x_iy_i - \sum_{j \ne i}|a_{ij}|x_j \cdot |b_{ij}|y_j.$

    由 Cauchy-Schwarz 型不等式和 $(Ax)_i = a_{ii}x_i + \sum_{j\ne i}a_{ij}x_j > 0$，$(By)_i = b_{ii}y_i + \sum_{j\ne i}b_{ij}y_j > 0$，可以验证 $(A \star B)z > 0$。

    更具体地，$(Ax)_i = a_{ii}x_i - \sum_{j\ne i}|a_{ij}|x_j > 0$，$(By)_i = b_{ii}y_i - \sum_{j\ne i}|b_{ij}|y_j > 0$。

    $((A\star B)z)_i = a_{ii}b_{ii}x_iy_i - \sum_{j\ne i}|a_{ij}||b_{ij}|x_jy_j \ge a_{ii}x_i \cdot b_{ii}y_i - \left(\sum_{j\ne i}|a_{ij}|x_j\right)\left(\sum_{j\ne i}|b_{ij}|y_j\right) \cdot \frac{1}{\text{适当因子}}.$

    精确的估计利用了：$\sum |a_{ij}||b_{ij}|x_jy_j \le \left(\sum |a_{ij}|x_j^2\right)^{1/2}\left(\sum |b_{ij}|y_j^2\right)^{1/2}$（不够直接）。

    另一种证明：利用条件 M10，存在正对角矩阵 $D_A, D_B$ 使 $D_A A, D_B B$ 严格行对角占优。可以验证 $A \star B$ 也满足某种对角占优条件（通过 AM-GM 不等式），从而是 M-矩阵。

    行列式不等式 $\det(A \star B) \ge \det(A)\det(B)$ 可以通过对 $n$ 归纳并利用 Schur 补来证明。

---

## 38A.8 逆正矩阵与单调矩阵

!!! definition "定义 38A.7 (逆正矩阵)"
    非奇异矩阵 $A \in M_n(\mathbb{R})$ 称为**逆正矩阵**（inverse-positive），若 $A^{-1} \ge 0$。称为**逆严格正矩阵**，若 $A^{-1} > 0$（逐元严格正）。

!!! definition "定义 38A.8 (单调矩阵)"
    非奇异矩阵 $A \in M_n(\mathbb{R})$ 称为**单调矩阵**（monotone matrix），若
    $$Ax \ge 0 \implies x \ge 0.$$

!!! theorem "定理 38A.9 (逆正等价于单调)"
    $A$ 是逆正矩阵当且仅当 $A$ 是单调矩阵。

??? proof "证明"
    $(\Rightarrow)$：设 $A^{-1} \ge 0$ 且 $Ax \ge 0$。则 $x = A^{-1}(Ax)$。非负矩阵 $A^{-1}$ 与非负向量 $Ax$ 的乘积非负，故 $x \ge 0$。

    $(\Leftarrow)$：设 $A$ 单调。对每个标准基向量 $e_j \ge 0$，$A(A^{-1}e_j) = e_j \ge 0$。由单调性，$A^{-1}e_j \ge 0$。$A^{-1}$ 的第 $j$ 列即 $A^{-1}e_j \ge 0$。对所有 $j$ 成立，故 $A^{-1} \ge 0$。$\blacksquare$

!!! theorem "定理 38A.10 (不可约 M-矩阵的逆严格正)"
    设 $A$ 为不可约非奇异 M-矩阵。则 $A^{-1} > 0$（逐元严格正）。

??? proof "证明"
    $A = sI - B$，$B \ge 0$ 不可约，$s > \rho(B)$。

    $$A^{-1} = s^{-1}\sum_{k=0}^{\infty}(s^{-1}B)^k = s^{-1}I + s^{-2}B + s^{-3}B^2 + \cdots$$

    由 $B \ge 0$ 不可约，Perron-Frobenius 定理保证 $(I + B)^{n-1} > 0$（不可约非负矩阵的 $(n-1)$ 次幂严格正）。因此级数中包含 $B^{n-1}$ 的项满足 $s^{-n}B^{n-1} > 0$（逐元严格正）。由于级数每项非负且至少有一项严格正，级数之和 $A^{-1} > 0$。

    更精确地：对任意 $i, j$，$(B^{n-1})_{ij} > 0$（不可约），故 $(A^{-1})_{ij} \ge s^{-n}(B^{n-1})_{ij} > 0$。$\blacksquare$

!!! example "例 38A.8"
    一维 Poisson 方程 $-u'' = f$ 在 $[0,1]$ 上有限差分离散化（Dirichlet 边界）给出 $Ah = f_h$，其中
    $$A = \frac{1}{h^2}\begin{pmatrix} 2 & -1 & & \\ -1 & 2 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 2 \end{pmatrix}.$$
    $A$ 是三对角不可约 M-矩阵，故 $A^{-1} > 0$。

    **物理意义**：$f \ge 0$（源项非负）$\implies$ $u_h = A^{-1}f_h \ge 0$（数值解非负）。这就是离散**极大值原理**，它保证了数值解保持物理解的定性行为。这是 M-矩阵在 PDE 数值方法中最核心的应用。

---

## 38A.9 正则分裂与迭代法收敛

!!! definition "定义 38A.9 (正则分裂)"
    矩阵 $A$ 的分裂 $A = M - N$ 称为**正则分裂**（regular splitting），若 $M$ 非奇异、$M^{-1} \ge 0$、$N \ge 0$。

!!! definition "定义 38A.10 (弱正则分裂)"
    分裂 $A = M - N$ 称为**弱正则分裂**（weak regular splitting），若 $M$ 非奇异、$M^{-1} \ge 0$、$M^{-1}N \ge 0$。

    注意正则分裂 $\implies$ 弱正则分裂（$M^{-1} \ge 0, N \ge 0 \implies M^{-1}N \ge 0$），但反之不然。

!!! definition "定义 38A.11 (M-分裂)"
    分裂 $A = M - N$ 称为 **M-分裂**，若 $M$ 是非奇异 M-矩阵且 $N \ge 0$。

    M-分裂是正则分裂的特例（$M$ 是 M-矩阵意味着 $M^{-1} \ge 0$）。

!!! theorem "定理 38A.11 (正则分裂的收敛性——Varga)"
    设 $A = M - N$ 为正则分裂，$T = M^{-1}N$ 为迭代矩阵。则迭代
    $$x^{(k+1)} = Tx^{(k)} + M^{-1}b$$
    收敛（对任意初值）当且仅当 $A$ 是非奇异的且 $A^{-1} \ge 0$。

    当条件满足时，$\rho(T) = \rho(M^{-1}N) < 1$，且
    $$\rho(T) = \frac{\rho(A^{-1}N)}{1 + \rho(A^{-1}N)}.$$

??? proof "证明"
    $T = M^{-1}N \ge 0$（非负矩阵，因为 $M^{-1} \ge 0, N \ge 0$）。迭代收敛当且仅当 $\rho(T) < 1$。

    **充分性**：设 $A$ 逆正（$A^{-1} \ge 0$）。$T = M^{-1}N = M^{-1}(M - A) = I - M^{-1}A$。

    $A^{-1}N = A^{-1}(M - A) = A^{-1}M - I$。由 $A^{-1} \ge 0, M^{-1} \ge 0$，Perron-Frobenius 适用于 $T \ge 0$。

    设 $\lambda = \rho(T)$ 是 $T$ 的 Perron 特征值，$Tv = \lambda v$，$v \ge 0, v \ne 0$。则 $(I - M^{-1}A)v = \lambda v$，$(1-\lambda)v = M^{-1}Av$，$Av = (1-\lambda)Mv$。故 $v = A^{-1}(1-\lambda)Mv$，$(1-\lambda)A^{-1}Mv = v$。

    又 $A^{-1}Nv = A^{-1}(M-A)v = A^{-1}Mv - v$，$Tv = M^{-1}Nv = \lambda v$ 意味着 $Nv = \lambda Mv$，$A^{-1}Nv = \lambda A^{-1}Mv$。

    设 $\mu = \rho(A^{-1}N)$。$A^{-1}N \ge 0$（$A^{-1} \ge 0, N \ge 0$），故 Perron-Frobenius 适用。由 $Tv = \lambda v$ 和 $T = I - M^{-1}A$：

    $Nv = \lambda Mv$，$A^{-1}Nv = \lambda A^{-1}Mv = \lambda A^{-1}(A + N)v = \lambda(v + A^{-1}Nv)$。

    令 $w = A^{-1}Nv$，则 $w = \lambda v + \lambda w$，$(1-\lambda)w = \lambda v$，$w = \frac{\lambda}{1-\lambda}v$（当 $\lambda < 1$）。

    故 $A^{-1}Nv = \frac{\lambda}{1-\lambda}v$，即 $\frac{\lambda}{1-\lambda}$ 是 $A^{-1}N$ 的特征值。由 $v \ge 0$ 和 Perron-Frobenius，$\mu = \rho(A^{-1}N) = \frac{\lambda}{1-\lambda}$（因为 $\lambda = \rho(T)$ 对应 Perron 向量）。解得 $\lambda = \frac{\mu}{1+\mu}$。

    由 $A^{-1} \ge 0, N \ge 0$，$\mu \ge 0$，故 $\lambda = \frac{\mu}{1+\mu} < 1$。

    **必要性**：设 $\rho(T) < 1$。$(I - T)^{-1} = \sum_{k=0}^{\infty}T^k \ge 0$（$T \ge 0$）。$I - T = M^{-1}A$，$(M^{-1}A)^{-1} = A^{-1}M \ge 0$。$A^{-1} = A^{-1}M \cdot M^{-1} \ge 0$（需更精细论证：$A^{-1} = (I-T)^{-1}M^{-1} = (\sum T^k)M^{-1} \ge 0$）。$\blacksquare$

!!! theorem "定理 38A.12 (M-矩阵的 Jacobi 和 Gauss-Seidel 收敛)"
    设 $A$ 为非奇异 M-矩阵，$A = D - L - U$（$D$ 为正对角部分，$L \ge 0$ 为严格下三角，$U \ge 0$ 为严格上三角）。则：

    (a) Jacobi 迭代收敛：$\rho(D^{-1}(L+U)) < 1$。

    (b) Gauss-Seidel 迭代收敛：$\rho((D-L)^{-1}U) < 1$。

    (c) 更一般地，$A$ 的每个正则分裂或 M-分裂给出的迭代都收敛。

??? proof "证明"
    (a) $A = D - (L+U)$。$D > 0$（$A$ 是 M-矩阵，对角元素 $> 0$），$D^{-1} \ge 0$。$L + U \ge 0$（$A$ 是 Z-矩阵，非对角元素 $\le 0$，故 $-a_{ij} = (L+U)_{ij} \ge 0$）。这是正则分裂（$M = D, N = L + U$）。由 $A$ 逆正（条件 M3）和定理 38A.11，$\rho(D^{-1}(L+U)) < 1$。

    (b) $A = (D - L) - U$。$D - L$ 是下三角 M-矩阵（因为 $A$ 是 M-矩阵，$D - L$ 是 $A$ 的某种"截断"；具体地 $(D-L)_{ii} = a_{ii} > 0$，非对角元素 $(D-L)_{ij} = a_{ij} \le 0$（$i > j$），$(D-L)_{ij} = 0$（$i < j$）——是 Z-矩阵。$D - L$ 是下三角且对角元素为正，故非奇异，$(D-L)^{-1} \ge 0$（下三角 M-矩阵的逆非负）。$U \ge 0$。这是正则分裂，由定理 38A.11 收敛。

    (c) 由定理 38A.11，任何正则分裂在 $A$ 逆正时都收敛。

---

## 38A.10 Varga 比较定理

!!! theorem "定理 38A.13 (Varga 比较定理)"
    设 $A = M_1 - N_1 = M_2 - N_2$ 为非奇异矩阵 $A$ 的两个正则分裂。若 $N_2 \ge N_1 \ge 0$（即第二个分裂中"分出"的非负部分更大），则
    $$\rho(M_1^{-1}N_1) \le \rho(M_2^{-1}N_2) < 1$$
    （假设 $A$ 逆正）。

    换言之，分裂中非负部分 $N$ 越小，迭代收敛越快。

??? proof "证明"
    设 $T_1 = M_1^{-1}N_1$，$T_2 = M_2^{-1}N_2$。两者均为非负矩阵，$\rho(T_1), \rho(T_2) < 1$（由定理 38A.11，$A$ 逆正）。

    **关键恒等式**：$T_2 - T_1 = M_2^{-1}N_2 - M_1^{-1}N_1$。

    $M_2 = A + N_2$，$M_1 = A + N_1$，故 $M_2^{-1} = (A + N_2)^{-1}$，$M_1^{-1} = (A + N_1)^{-1}$。

    $T_2 = M_2^{-1}N_2 = I - M_2^{-1}A$，$T_1 = I - M_1^{-1}A$。

    $T_2 - T_1 = M_1^{-1}A - M_2^{-1}A = (M_1^{-1} - M_2^{-1})A$。

    $M_1^{-1} - M_2^{-1} = M_1^{-1}(M_2 - M_1)M_2^{-1} = M_1^{-1}(N_2 - N_1)M_2^{-1}$。

    由 $N_2 \ge N_1$，$N_2 - N_1 \ge 0$。$M_1^{-1} \ge 0, M_2^{-1} \ge 0$。故 $M_1^{-1}(N_2-N_1)M_2^{-1} \ge 0$。

    因此 $T_2 - T_1 = M_1^{-1}(N_2 - N_1)M_2^{-1}A$。当 $A$ 不一定非负时这不直接给出 $T_2 \ge T_1$。

    更正确的做法：$T_2 = I - M_2^{-1}A$。由 $(I - T_2)^{-1} = A^{-1}M_2 = A^{-1}(A + N_2) = I + A^{-1}N_2 \ge I + A^{-1}N_1 = A^{-1}M_1 = (I-T_1)^{-1}$。

    由 $T_i \ge 0$ 和 $(I - T_i)^{-1} = \sum_{k=0}^{\infty}T_i^k$，利用 Perron-Frobenius 定理的单调性：

    对非负矩阵 $0 \le S \le T$，$\rho(S) \le \rho(T)$。

    $(I-T_2)^{-1} \ge (I-T_1)^{-1} \ge I$ 意味着 $\sum T_2^k \ge \sum T_1^k$。虽然这不直接给出 $T_2 \ge T_1$，但由 Perron 特征值的变分刻画：

    $\rho(T_i) = \frac{\rho(A^{-1}N_i)}{1 + \rho(A^{-1}N_i)}$（由定理 38A.11 的公式），而 $A^{-1}N_2 \ge A^{-1}N_1 \ge 0$（$A^{-1} \ge 0, N_2 \ge N_1 \ge 0$），故 $\rho(A^{-1}N_2) \ge \rho(A^{-1}N_1)$，从而 $\rho(T_2) \ge \rho(T_1)$。$\blacksquare$

!!! example "例 38A.9"
    $A = \begin{pmatrix} 4 & -1 \\ -2 & 5 \end{pmatrix}$（M-矩阵）。

    **Jacobi 分裂**：$M_1 = D = \begin{pmatrix} 4 & 0 \\ 0 & 5 \end{pmatrix}$，$N_1 = \begin{pmatrix} 0 & 1 \\ 2 & 0 \end{pmatrix}$。

    **Gauss-Seidel 分裂**：$M_2 = D - L = \begin{pmatrix} 4 & 0 \\ -2 & 5 \end{pmatrix}$，$N_2 = U = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    此处 $N_1 = \begin{pmatrix} 0 & 1 \\ 2 & 0 \end{pmatrix} \ge N_2 = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    由 Varga 比较定理，$\rho(T_{GS}) \le \rho(T_J)$。实际计算：$T_J = \begin{pmatrix} 0 & 1/4 \\ 2/5 & 0 \end{pmatrix}$，$\rho(T_J) = \sqrt{2/20} = \sqrt{1/10}$。$T_{GS}$ 的谱半径更小。

---

## 38A.11 Stein-Rosenberg 定理

!!! theorem "定理 38A.14 (Stein-Rosenberg 定理)"
    设 $A$ 为 Z-矩阵，$A = D - L - U$（$D$ 正对角，$L, U \ge 0$），且 $D^{-1}(L+U)$ 的对角元素为零（自动满足）。设 $T_J = D^{-1}(L+U)$ 为 Jacobi 迭代矩阵，$T_{GS} = (D-L)^{-1}U$ 为 Gauss-Seidel 迭代矩阵。则恰好有以下四种互斥情形之一成立：

    (i) $\rho(T_J) = \rho(T_{GS}) = 0$。

    (ii) $0 < \rho(T_{GS}) < \rho(T_J) < 1$（两者均收敛，GS 更快）。

    (iii) $\rho(T_J) = \rho(T_{GS}) = 1$。

    (iv) $1 < \rho(T_J) < \rho(T_{GS})$（两者均发散，GS 更差）。

    特别地，Jacobi 迭代收敛当且仅当 Gauss-Seidel 迭代收敛。

??? proof "证明"
    **情形 (i)**：$T_J = 0$ 当且仅当 $L = U = 0$（即 $A$ 为对角矩阵），此时 $T_{GS} = 0$。

    **核心论证**：设 $T_J$ 不可约（可约情形可通过 Frobenius 标准形化归）。$T_J \ge 0$ 且对角元素为零。

    考虑参数化矩阵 $T(\alpha) = (I - \alpha D^{-1}L)^{-1}(\alpha D^{-1}U + (1-\alpha)D^{-1}L)$。当 $\alpha = 0$ 时 $T(0) = D^{-1}(L+U) = T_J$；当 $\alpha = 1$ 时 $T(1) = (I - D^{-1}L)^{-1}D^{-1}U$，这与 $T_{GS}$ 相关。

    Perron 根 $r(\alpha) = \rho(T(\alpha))$ 是 $\alpha$ 的连续函数。关键的代数关系来自特征方程分析：

    设 $\lambda$ 是 $T_J$ 的特征值，$T_J v = \lambda v$。$D^{-1}(L+U)v = \lambda v$。将 $L, U$ 分离，可以证明若 $\mu$ 是 $T_{GS}$ 的特征值，则存在 $T_J$ 的特征值 $\lambda$ 使得 $\mu = \lambda^2$（在 $A$ 具有"一致分块结构"时成立，例如当 $T_J^2$ 的特征值与 $T_{GS}$ 的特征值匹配时）。

    对一般的非负 $T_J$，利用 Perron-Frobenius 的分析性质和 Neumann 级数的单调性可以证明：

    - 若 $\rho(T_J) < 1$，则 $(I - T_J)^{-1} = \sum T_J^k \ge 0$ 存在，$A = D(I - T_J)$ 逆正。$(I - T_{GS})^{-1} = \sum T_{GS}^k$ 也存在，且由 Gauss-Seidel 作为正则分裂的特殊性，$\rho(T_{GS}) < 1$。严格不等式 $\rho(T_{GS}) < \rho(T_J)$（当 $\rho(T_J) > 0$）来自于 $T_{GS}$ 比 $T_J$ 的"更优"的分裂结构。

    - 若 $\rho(T_J) = 1$，则类似分析给出 $\rho(T_{GS}) = 1$。

    - 若 $\rho(T_J) > 1$，则 Neumann 级数发散，$A$ 不是逆正的，GS 也发散，且 $\rho(T_{GS}) > \rho(T_J)$。

    完整的技术细节见 Varga (2000), Matrix Iterative Analysis, 定理 3.30。$\blacksquare$

!!! example "例 38A.10"
    $A = \begin{pmatrix} 4 & -1 \\ -1 & 4 \end{pmatrix}$（M-矩阵）。

    $T_J = \begin{pmatrix} 0 & 1/4 \\ 1/4 & 0 \end{pmatrix}$，$\rho(T_J) = 1/4$。

    $T_{GS} = (D-L)^{-1}U$。$D - L = \begin{pmatrix} 4 & 0 \\ -1 & 4 \end{pmatrix}$，$(D-L)^{-1} = \frac{1}{16}\begin{pmatrix} 4 & 0 \\ 1 & 4 \end{pmatrix}$。$U = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$。

    $T_{GS} = \frac{1}{16}\begin{pmatrix} 4 & 0 \\ 1 & 4 \end{pmatrix}\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix} = \begin{pmatrix} 0 & 1/4 \\ 0 & 1/16 \end{pmatrix}$。

    $\rho(T_{GS}) = 1/16 = (1/4)^2 = \rho(T_J)^2$。

    这验证了情形 (ii)：$0 < 1/16 < 1/4 < 1$，且在此对称情形下 $\rho(T_{GS}) = \rho(T_J)^2$。

---

## 38A.12 SOR 收敛与 M-矩阵

!!! definition "定义 38A.12 (SOR 迭代)"
    设 $A = D - L - U$。**逐次超松弛法**（SOR, Successive Over-Relaxation）以松弛参数 $\omega > 0$ 定义迭代矩阵
    $$T_\omega = (D - \omega L)^{-1}((1-\omega)D + \omega U).$$
    当 $\omega = 1$ 时退化为 Gauss-Seidel。

!!! theorem "定理 38A.15 (SOR 对 M-矩阵的收敛性)"
    设 $A$ 为非奇异 M-矩阵。则：

    (a) 对所有 $\omega \in (0, 1]$，SOR 迭代收敛：$\rho(T_\omega) < 1$。

    (b) 对 $\omega > 0$，$\rho(T_\omega) \ge |1 - \omega|$（Kahan 必要条件）。

    (c) 若 $A$ 对称正定（对称 M-矩阵），则对所有 $\omega \in (0, 2)$，SOR 收敛。

??? proof "证明"
    (a) 当 $0 < \omega \le 1$ 时，$A = (D - \omega L)/\omega - ((1-\omega)D/\omega + U)$。可以验证这是正则分裂（$\omega \le 1$ 保证 $N = (1-\omega)D/\omega + U \ge 0$，而 $(D - \omega L)/\omega$ 是 M-矩阵的变体，逆非负）。由定理 38A.11 和 $A$ 逆正，$\rho(T_\omega) < 1$。

    (b) $\det(T_\omega) = \det((D-\omega L)^{-1})\det((1-\omega)D + \omega U)$。$\det(D - \omega L) = \prod d_{ii}$（下三角），$\det((1-\omega)D + \omega U) = (1-\omega)^n \prod d_{ii}$（上三角）。故 $|\det(T_\omega)| = |1-\omega|^n$。由 $|\det(T_\omega)| = \prod |\lambda_i(T_\omega)| \le \rho(T_\omega)^n$，得 $|1-\omega| \le \rho(T_\omega)$。

    (c) 对称正定 M-矩阵的 SOR 收敛性在 $\omega \in (0,2)$ 范围内由 Ostrowski-Reich 定理保证（见第 38B 章）。

---

## 38A.13 弱正则分裂与 Marek-Szyld 理论

!!! theorem "定理 38A.16 (弱正则分裂的收敛性——Marek-Szyld)"
    设 $A = M - N$ 为弱正则分裂（$M^{-1} \ge 0$，$M^{-1}N \ge 0$）。则 $\rho(M^{-1}N) < 1$ 当且仅当 $A$ 非奇异且 $A^{-1} \ge 0$。

    这推广了正则分裂的收敛定理。关键区别在于弱正则分裂不要求 $N \ge 0$，只要求 $T = M^{-1}N \ge 0$。

??? proof "证明"
    必要性和充分性的证明与定理 38A.11 完全类似，因为核心论证只依赖于 $T = M^{-1}N \ge 0$ 这一条件，而不直接依赖于 $N \ge 0$。

    **充分性**：$T = M^{-1}N \ge 0$ 是非负矩阵。$I - T = M^{-1}A$，$(I-T)^{-1} = A^{-1}M$。由 $A^{-1} \ge 0, M^{-1} \ge 0$，$A^{-1}M \ge 0$，$(I-T)^{-1} \ge 0$。由 $(I-T)^{-1} = \sum_{k=0}^{\infty}T^k$（当 $\rho(T) < 1$），由 Perron-Frobenius 理论，对非负 $T$，$(I-T)^{-1} \ge 0$ 当且仅当 $\rho(T) < 1$。

    **必要性**：$\rho(T) < 1$ 意味着 $(I-T)^{-1} = \sum T^k \ge 0$。$A^{-1} = (I-T)^{-1}M^{-1} \ge 0$。$\blacksquare$

!!! theorem "定理 38A.17 (弱正则分裂的比较)"
    设 $A = M_1 - N_1 = M_2 - N_2$ 为弱正则分裂，$T_1 = M_1^{-1}N_1 \ge 0$，$T_2 = M_2^{-1}N_2 \ge 0$。若 $A^{-1} \ge 0$ 且 $A^{-1}N_2 \ge A^{-1}N_1 \ge 0$，则
    $$\rho(T_1) \le \rho(T_2) < 1.$$

??? proof "证明"
    由定理 38A.11 中的公式，$\rho(T_i) = \frac{\rho(A^{-1}N_i)}{1 + \rho(A^{-1}N_i)}$。$f(\mu) = \frac{\mu}{1+\mu}$ 在 $\mu \ge 0$ 上单调递增。$A^{-1}N_2 \ge A^{-1}N_1 \ge 0$ 意味着 $\rho(A^{-1}N_2) \ge \rho(A^{-1}N_1)$（非负矩阵的谱半径关于逐元序单调）。故 $\rho(T_2) \ge \rho(T_1)$。$\blacksquare$

---

## 38A.14 Leontief 投入产出模型与 Hawkins-Simon 条件

!!! definition "定义 38A.13 (Leontief 投入产出模型)"
    设经济由 $n$ 个产业部门组成。$c_{ij} \ge 0$ 表示第 $j$ 部门生产单位产品所需第 $i$ 部门的投入量。矩阵 $C = (c_{ij}) \ge 0$ 称为**投入系数矩阵**（或技术矩阵）。若 $x = (x_1, \ldots, x_n)^T$ 为各部门的总产出向量，$d = (d_1, \ldots, d_n)^T \ge 0$ 为最终需求向量，则平衡方程为
    $$x = Cx + d, \quad \text{即} \quad (I - C)x = d.$$

!!! theorem "定理 38A.18 (Hawkins-Simon 条件)"
    设 $C \ge 0$ 为投入系数矩阵。以下条件等价：

    (a) 对每个最终需求 $d \ge 0$，平衡方程 $(I - C)x = d$ 有唯一非负解 $x \ge 0$。

    (b) $I - C$ 是非奇异 M-矩阵。

    (c) **Hawkins-Simon 条件**：$I - C$ 的所有顺序主子式为正。

    (d) $\rho(C) < 1$。

    (e) **Leontief 逆**非负：$(I - C)^{-1} \ge 0$。

    (f) Leontief 逆可以展开为 Neumann 级数：$(I-C)^{-1} = \sum_{k=0}^{\infty}C^k$，其中 $C^k$ 代表"第 $k$ 轮间接投入"。

??? proof "证明"
    $I - C$ 是 Z-矩阵（非对角元素 $-c_{ij} \le 0$），且 $I - C = I - C$（$s = 1, B = C$）。

    (a) $\Leftrightarrow$ (e)：对所有 $d \ge 0$，$x = (I-C)^{-1}d \ge 0$ 当且仅当 $(I-C)^{-1} \ge 0$。

    (b) $\Leftrightarrow$ (d)：$I - C = I - C$ 是 M-矩阵 $\iff$ $1 > \rho(C)$。

    (b) $\Leftrightarrow$ (e)：M-矩阵等价条件 M3。

    (b) $\Leftrightarrow$ (c)：M-矩阵等价条件 M4。

    (d) $\Rightarrow$ (f)：$\rho(C) < 1$ 保证级数 $\sum C^k$ 绝对收敛，且和为 $(I-C)^{-1}$。$\blacksquare$

!!! example "例 38A.11"
    三部门经济，投入系数矩阵
    $$C = \begin{pmatrix} 0.2 & 0.3 & 0.1 \\ 0.1 & 0.1 & 0.2 \\ 0.2 & 0.1 & 0.1 \end{pmatrix}.$$

    $I - C = \begin{pmatrix} 0.8 & -0.3 & -0.1 \\ -0.1 & 0.9 & -0.2 \\ -0.2 & -0.1 & 0.9 \end{pmatrix}$。

    **Hawkins-Simon 检验**：$\Delta_1 = 0.8 > 0$，$\Delta_2 = 0.72 - 0.03 = 0.69 > 0$，

    $\Delta_3 = 0.8(0.81 - 0.02) + 0.3(-0.09 + 0.04) + (-0.1)(0.01 + 0.18)$
    $= 0.8 \times 0.79 + 0.3 \times (-0.05) - 0.1 \times 0.19 = 0.632 - 0.015 - 0.019 = 0.598 > 0$。

    所有顺序主子式为正，$I - C$ 是 M-矩阵，经济可行。

    **经济解读**：Hawkins-Simon 条件 $\Delta_k > 0$ 意味着前 $k$ 个部门的"子经济"有正的净产出能力——总产出能覆盖相互间的投入需求并有剩余。

!!! theorem "定理 38A.19 (Leontief 乘数效应)"
    若 $I - C$ 是 M-矩阵，则 Leontief 逆 $L = (I - C)^{-1} = I + C + C^2 + \cdots$ 的 $(i,j)$ 元素 $l_{ij}$ 表示第 $j$ 部门的最终需求增加 1 单位时，第 $i$ 部门的总产出增加量。$l_{ij} \ge \delta_{ij}$，反映**乘数效应**：

    - $I$：直接需求（第 $j$ 部门自身增产 1 单位）；
    - $C$：第一轮间接需求（直接投入）；
    - $C^2$：第二轮间接需求（投入的投入）；
    - $C^k$：第 $k$ 轮间接需求。

??? proof "证明"
    $\frac{\partial x_i}{\partial d_j} = [(I-C)^{-1}]_{ij} = l_{ij}$。$(I-C)^{-1} = I + C + C^2 + \cdots \ge I$，故 $l_{ij} \ge \delta_{ij}$。$\blacksquare$

!!! example "例 38A.12 (开放与封闭 Leontief 模型)"
    **开放模型**（上述情形）：最终需求 $d$ 是外生的。可行性等价于 $I - C$ 是 M-矩阵。

    **封闭模型**：将消费部门也纳入产业。此时 $\tilde{C}$ 的列和为 1（所有产出被消耗），$\rho(\tilde{C}) = 1$，$I - \tilde{C}$ 是奇异 M-矩阵。方程 $(I - \tilde{C})x = 0$ 有非负解 $x \ge 0$（由奇异 M-矩阵的零空间性质），代表均衡产出比例（水平任意）。

---

## 38A.15 M-矩阵在其他领域中的出现

!!! example "例 38A.13 (PDE 离散化)"
    椭圆型偏微分方程 $-\nabla \cdot (a \nabla u) + cu = f$（$a > 0, c \ge 0$）在适当的网格和差分格式下，离散化矩阵是 M-矩阵。这保证了：

    - **离散极大值原理**：$f \ge 0 \implies u_h \ge 0$。
    - **解的唯一性和稳定性**。
    - **迭代求解法的收敛性**（Jacobi, Gauss-Seidel, SOR 等）。

    对有限元方法，在网格满足 Delaunay 条件时，刚度矩阵是 M-矩阵。

!!! example "例 38A.14 (种群动力学)"
    竞争 Lotka-Volterra 系统 $\dot{x}_i = x_i(r_i - \sum_j a_{ij}x_j)$ 中，群落矩阵 $A = (a_{ij})$ 在"每种群自限（$a_{ii} > 0$）且种间竞争相对种内竞争不太强"的条件下，$\operatorname{diag}(a_{11},\ldots,a_{nn}) - A$ 的某种变体是 M-矩阵。这保证了正平衡点的存在性和全局稳定性。

!!! example "例 38A.15 (电网分析)"
    电阻网络的节点导纳矩阵 $Y$ 满足：$Y_{ii} = \sum_{j \sim i} g_{ij} \ge 0$（对角元为连接到节点 $i$ 的所有电导之和），$Y_{ij} = -g_{ij} \le 0$（$i \ne j$，非对角元为负电导）。$Y$ 是 Z-矩阵。若网络连通（$Y$ 不可约），则 $Y$（去掉接地节点后）是不可约 M-矩阵。

---

## 本章小结

| 概念 | 定义 | 核心性质 |
|------|------|----------|
| Z-矩阵 | 非对角元 $\le 0$ | $A = sI - B$ ($B \ge 0$)；主子矩阵遗传 |
| 非奇异 M-矩阵 ($K$ 类) | Z-矩阵 + $s > \rho(B)$ | $A^{-1} \ge 0$；22+ 等价条件；Neumann 级数收敛 |
| 奇异 M-矩阵 ($K_0$ 类) | Z-矩阵 + $s = \rho(B)$ | 零特征值；$\operatorname{rank} = n-1$（不可约时）；正零空间向量 |
| 逆正矩阵 | $A^{-1} \ge 0$ | 等价于单调矩阵；$\supset$ M-矩阵 |
| 正则分裂 | $M^{-1} \ge 0, N \ge 0$ | 收敛 $\iff$ $A$ 逆正 |
| M-分裂 | $M$ 是 M-矩阵，$N \ge 0$ | 正则分裂的特例 |

包含关系：
$$\text{严格对角占优 Z-矩阵} \subset \text{非奇异 M-矩阵} \subset \text{Z-矩阵} \cap \text{逆正矩阵}$$

---

## 练习题

1. **[基础] 验证矩阵 $A = \begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$ 是否为非奇异 M-矩阵。它的逆矩阵是什么？是否为非负矩阵？**
   ??? success "参考答案"
       1. Z-矩阵：非对角元 $-1 \le 0$。
       2. 谱半径：$A = 2I - B, B = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$。$\rho(B)=1 < 2$。故是非奇异 M-矩阵。
       3. 逆矩阵：$A^{-1} = \frac{1}{3} \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。所有元素均大于 0，显然是非负矩阵（逆正性）。

2. **[等价性] 利用主子式判据（Hawkins-Simon 条件），判定 $A = \begin{pmatrix} 1 & -2 \\ -1 & 1 \end{pmatrix}$ 是否为 M-矩阵。**
   ??? success "参考答案"
       计算主子式：$\Delta_1 = 1 > 0$，$\Delta_2 = \det(A) = 1 - 2 = -1 < 0$。
       由于存在负的主子式，该 Z-矩阵不是 M-矩阵（尽管对角元为正）。这说明 M-矩阵要求对角元的“支配力”足够强。

3. **[对角占优] 证明：如果一个 Z-矩阵是严格行对角占优的，那么它一定是非奇异 M-矩阵。**
   ??? success "参考答案"
       由 Gershgorin 圆盘定理，特征值 $\lambda$ 满足 $|\lambda - a_{ii}| \le \sum_{j \ne i} |a_{ij}| < a_{ii}$。
       这保证了特征值 $\lambda$ 的实部必须严格大于 0（圆盘完全位于复平面右侧且不含原点）。由于满足特征值实部为正的 Z-矩阵即为 M-矩阵，结论成立。

4. **[不可约] 设 $A$ 是不可约非奇异 M-矩阵。证明 $A^{-1}$ 的所有元素都严格大于 0。**
   ??? success "参考答案"
       利用级数展开 $A^{-1} = s^{-1} \sum_{k=0}^{\infty} (B/s)^k$。由于 $B$ 不可约且非负，存在某个幂 $B^m > 0$。由于级数包含所有幂次项且各项非负，相加后结果必为严格正。

5. **[分裂] 设 $A = M - N$ 是 M-矩阵的一个正则分裂。为什么由此产生的迭代法 $x^{(k+1)} = M^{-1}Nx^{(k)} + M^{-1}b$ 总是收敛的？**
   ??? success "参考答案"
       正则分裂意味着 $M^{-1} \ge 0$ 且 $N \ge 0$。对于 M-矩阵，$A^{-1} \ge 0$。Varga 收敛定理指出，正则分裂收敛的充要条件是 $A^{-1} \ge 0$。其谱半径满足 $\rho(M^{-1}N) = \frac{\rho(A^{-1}N)}{1 + \rho(A^{-1}N)} < 1$。

6. **[Schur补] 设 $A$ 是 M-矩阵，证明其任意主子矩阵的 Schur 补仍为 M-矩阵。这一性质在数值计算中有什么意义？**
   ??? success "参考答案"
       证明见定理 38A.7。在数值计算中，这意味着对 M-矩阵进行 LU 分解（不带换元）时，中间产生的每一个子矩阵都保持 M-矩阵的良好性质（如对角占优和逆非负性），从而保证了算法的数值稳定性。

7. **[奇异M] 给出一个 $2 \times 2$ 奇异 M-矩阵的例子，并求其零空间（Kernel）的生成向量。**
   ??? success "参考答案"
       取 $A = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$。
       它是奇异的（$\det=0$）且为 Z-矩阵。零空间由 $\mathbf{v} = (1, 1)^T$ 生成。注意该向量是严格正的，这符合不可约奇异 M-矩阵的理论预言。

8. **[比较] 设有两分裂 $A = M_1 - N_1 = M_2 - N_2$。若 $0 \le N_1 \le N_2$，哪一个分裂对应的迭代速度更快？**
   ??? success "参考答案"
       由 Varga 比较定理，$N$ 越小（即分出去的非负部分越少），迭代矩阵的谱半径 $\rho(M^{-1}N)$ 越小。因此第一个分裂收敛更快。

9. **[经济学] 在 Leontief 投入产出模型中，投入矩阵 $C \ge 0$。为什么要求 $\rho(C) < 1$？**
   ??? success "参考答案"
       因为只有当 $\rho(C) < 1$ 时，系统矩阵 $A = I - C$ 才是非奇异 M-矩阵，从而保证了对于任何正的外部需求 $\mathbf{d}$，总能通过产出 $\mathbf{x} = (I-C)^{-1}\mathbf{d}$ 得到非负的生产方案。若 $\rho(C) \ge 1$，意味着系统消耗的资源超过了产出，经济不可持续。

10. **[爱因斯坦思考题] M-矩阵的逆 $A^{-1} \ge 0$ 意味着“正向的输入产生正向的输出”。在扩散方程或热传导方程的离散化中，这对应于所谓的“极大值原理”。为什么爱因斯坦会认为一个描述真实物理演化的矩阵（如热传导矩阵）必须是 M-矩阵？**
    ??? success "参考答案"
        因为如果对应的矩阵不是 M-矩阵（即不满足极大值原理），那么在没有任何热源的情况下，空间中某点的温度可能会自发地从零度降为负数，或者在局部产生没有物理来源的振荡。M-矩阵结构反映了宇宙中“因果律”和“能量耗散”的单调性：扰动总是向低势能方向传播，且系统的响应必须与激励的方向保持一致。这种代数结构是防止物理模拟中出现“负概率”或“自发升温”现象的数学防线。

## 本章小结

本章系统探讨了具有特殊符号模式和谱性质的 M-矩阵及其在多领域的深远影响：

1. **Z-矩阵与 M-矩阵的定义**：确立了非对角元非正作为基础模式，并结合谱半径约束定义了 M-矩阵类。
2. **丰富的等价刻画**：整合了包括逆非负性、主子式全正、半正性等在内的二十余种等价条件，展现了 M-矩阵理论的宏大。
3. **结构保持性**：证明了 M-矩阵在取主子矩阵、求 Schur 补和进行 Fan 积下的封闭性，为块分解提供了理论依据。
4. **迭代法的灵魂**：揭示了正则分裂收敛性与 M-矩阵结构的内在统一，通过 Varga 比较定理和 Stein-Rosenberg 定理确立了求解大型稀疏系统的算法基础。
5. **物理与经济建模**：通过 Leontief 模型和离散极大值原理，展示了 M-矩阵如何作为因果一致性和系统可行性的数学语言。

