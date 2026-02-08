# 第 70 章 线性代数在生物学中的应用

<div class="context-flow" markdown>

**前置**：特征值(Ch6) · 非负矩阵(Ch17) · 矩阵指数(Ch13) · 微分方程(Ch26)

**本章脉络**：Leslie 矩阵(种群模型) $\to$ Perron-Frobenius 的生态学解释 $\to$ Lotka-Volterra 系统线性化 $\to$ 房室模型(药代动力学) $\to$ SEIR 流行病模型 $\to$ 序列比对评分矩阵 $\to$ 系统发育 $\to$ 基因调控网络

**延伸**：Leslie 矩阵模型是种群生态学和渔业管理的标准工具；房室模型是药代动力学和放射性示踪剂动力学的基础；基因调控网络分析推动了系统生物学的发展

</div>

生物学在 20 世纪下半叶经历了深刻的数学化转变。从种群动力学到分子生物学，线性代数提供了不可或缺的分析工具。Leslie 矩阵将年龄结构种群的增长编码为矩阵乘法，Perron-Frobenius 定理直接给出长期增长率和稳定年龄分布。在流行病学中，基本再生数 $R_0$ 可以表示为次代矩阵的谱半径。在基因组学中，序列比对的评分矩阵本质上是突变概率矩阵的对数。

本章将系统地展示矩阵方法如何贯穿生物学的各个层次——从种群到细胞到分子。

---

## 70.1 Leslie 矩阵与种群动力学

<div class="context-flow" markdown>

**核心问题**：一个具有年龄结构的种群，各年龄组有不同的繁殖力和存活率，其长期增长行为如何？

</div>

P. H. Leslie 于 1945 年提出了以矩阵形式描述年龄结构种群动力学的方法。这个模型至今仍是种群生态学和渔业管理的标准工具。

!!! definition "定义 70.1 (Leslie 矩阵)"
    将种群按年龄分为 $n$ 个年龄组。设 $\mathbf{n}_t = (n_1(t), n_2(t), \ldots, n_n(t))^T$ 为第 $t$ 个时间步各年龄组的个体数。定义 **Leslie 矩阵**：

    $$L = \begin{pmatrix} f_1 & f_2 & f_3 & \cdots & f_{n-1} & f_n \\ s_1 & 0 & 0 & \cdots & 0 & 0 \\ 0 & s_2 & 0 & \cdots & 0 & 0 \\ \vdots & & \ddots & & & \vdots \\ 0 & 0 & \cdots & s_{n-1} & 0 & 0 \end{pmatrix}$$

    其中：

    - $f_i \geq 0$ 为第 $i$ 个年龄组的**繁殖力**（每个个体在一个时间步内产生的后代数）
    - $s_i \in (0, 1]$ 为从第 $i$ 个年龄组存活到第 $i+1$ 个年龄组的**存活率**

    种群动力学方程为：

    $$\mathbf{n}_{t+1} = L \mathbf{n}_t$$

    因此 $\mathbf{n}_t = L^t \mathbf{n}_0$。

!!! theorem "定理 70.1 (Leslie 矩阵的特征方程)"
    Leslie 矩阵 $L$ 的特征方程为：

    $$\lambda^n = f_1 \lambda^{n-1} + f_2 s_1 \lambda^{n-2} + f_3 s_1 s_2 \lambda^{n-3} + \cdots + f_n s_1 s_2 \cdots s_{n-1}$$

    等价地，令 $\ell_i = s_1 s_2 \cdots s_{i-1}$（$\ell_1 = 1$）为存活到第 $i$ 个年龄组的概率，则特征方程为：

    $$1 = \sum_{i=1}^{n} f_i \ell_i \lambda^{-i}$$

??? proof "证明"
    计算 $\det(L - \lambda I)$。沿第一列展开：

    $$\det(L - \lambda I) = (f_1 - \lambda) \det(M_{11}) - s_1 \det(M_{21})$$

    其中 $M_{11}$ 和 $M_{21}$ 是相应的余子式矩阵。$M_{11}$ 是 $(n-1) \times (n-1)$ 的下三角矩阵（对角线上的 Leslie 矩阵去掉第一行第一列），其行列式涉及次对角线元素。

    通过归纳法或直接展开，可以得到特征多项式为：

    $$p(\lambda) = \lambda^n - f_1\lambda^{n-1} - f_2 s_1 \lambda^{n-2} - \cdots - f_n s_1 s_2 \cdots s_{n-1}$$

    将 $p(\lambda) = 0$ 两端除以 $\lambda^n$（$\lambda \neq 0$），得 $1 = \sum_{i=1}^n f_i \ell_i \lambda^{-i}$。

    $\blacksquare$

!!! theorem "定理 70.2 (Leslie 矩阵的主特征值)"
    若 Leslie 矩阵 $L$ 至少有两个非零繁殖力 $f_i, f_j > 0$（$i < j$），且 $\gcd$ 条件满足使得 $L$ 为本原矩阵，则 $L$ 有唯一的正主特征值 $\lambda_1 > 0$，且 $\lambda_1 > |\lambda_k|$ 对所有其他特征值 $\lambda_k$。

??? proof "证明"
    Leslie 矩阵 $L$ 是非负矩阵。若 $L$ 不可约且非周期（本原），则由 Perron-Frobenius 定理，$L$ 有唯一的正主特征值 $\lambda_1 = \rho(L)$，且对应的特征向量 $\mathbf{v}_1 > 0$。

    $L$ 不可约的充分条件是所有 $s_i > 0$ 且至少 $f_n > 0$（或更一般地，存在从最后一个年龄组到第一个年龄组的路径）。

    $L$ 非周期的充分条件是 $\gcd\{i : f_i > 0\} = 1$。例如，若 $f_1 > 0$ 和 $f_2 > 0$，则 $\gcd(1,2) = 1$。

    $\blacksquare$

!!! example "例 70.1"
    一个三龄组种群的 Leslie 矩阵：

    $$L = \begin{pmatrix} 0 & 4 & 3 \\ 0.5 & 0 & 0 \\ 0 & 0.25 & 0 \end{pmatrix}$$

    即第 1 龄组不繁殖（$f_1 = 0$），第 2 龄组繁殖力为 4，第 3 龄组为 3；从第 1 到第 2 龄组的存活率为 0.5，从第 2 到第 3 龄组的存活率为 0.25。

    特征方程：$\lambda^3 = 4 \times 0.5 \lambda + 3 \times 0.5 \times 0.25 = 2\lambda + 0.375$

    即 $\lambda^3 - 2\lambda - 0.375 = 0$。数值求解得主特征值 $\lambda_1 \approx 1.489$。

    这意味着种群长期按每个时间步增长 $48.9\%$ 的速率增长。对应的特征向量给出稳定年龄分布。

---

## 70.2 Perron-Frobenius 在生态学中的应用

<div class="context-flow" markdown>

**核心问题**：Leslie 矩阵的特征值和特征向量有什么生态学含义？

</div>

Perron-Frobenius 定理在种群生态学中有三个核心应用：确定长期增长率、稳定年龄分布和繁殖价值。

!!! theorem "定理 70.3 (长期增长率与稳定年龄分布)"
    设 Leslie 矩阵 $L$ 本原，主特征值为 $\lambda_1$，对应的右特征向量为 $\mathbf{v}_1$（$L\mathbf{v}_1 = \lambda_1 \mathbf{v}_1$），左特征向量为 $\mathbf{w}_1$（$\mathbf{w}_1^T L = \lambda_1 \mathbf{w}_1^T$）。则：

    1. **渐近增长率**：$\lambda_1$ 是种群的渐近增长率，即 $\mathbf{n}_t \sim c \lambda_1^t \mathbf{v}_1$
    2. **稳定年龄分布**：$\mathbf{v}_1 / \|\mathbf{v}_1\|_1$ 是种群的稳定年龄分布——无论初始年龄结构如何，最终各年龄组的比例趋向 $\mathbf{v}_1$ 的归一化形式
    3. **繁殖价值**：$\mathbf{w}_1$ 是各年龄组的**繁殖价值**——$w_i$ 度量第 $i$ 个年龄组的一个个体对未来种群增长的贡献

??? proof "证明"
    设 $L$ 的 Jordan 标准形为 $L = S \Lambda S^{-1}$（由于 $L$ 本原，可以对角化或使用谱分解）。则：

    $$L^t = S \Lambda^t S^{-1} = \lambda_1^t \mathbf{v}_1 \mathbf{w}_1^T / (\mathbf{w}_1^T \mathbf{v}_1) + \sum_{k \geq 2} \lambda_k^t \mathbf{v}_k \mathbf{w}_k^T / (\mathbf{w}_k^T \mathbf{v}_k)$$

    由于 $|\lambda_k / \lambda_1| < 1$ 对 $k \geq 2$，

    $$\frac{L^t}{\lambda_1^t} \to \frac{\mathbf{v}_1 \mathbf{w}_1^T}{\mathbf{w}_1^T \mathbf{v}_1}$$

    因此 $\mathbf{n}_t = L^t \mathbf{n}_0 \sim \lambda_1^t \frac{\mathbf{w}_1^T \mathbf{n}_0}{\mathbf{w}_1^T \mathbf{v}_1} \mathbf{v}_1$。

    这表明：(1) 增长率为 $\lambda_1$；(2) 年龄结构趋向 $\mathbf{v}_1$；(3) 初始条件通过 $\mathbf{w}_1^T \mathbf{n}_0$ 影响总量——$\mathbf{w}_1$ 的各分量给出了各年龄组个体的权重，即繁殖价值。

    $\blacksquare$

!!! definition "定义 70.2 (弹性分析与灵敏度分析)"
    主特征值 $\lambda_1$ 对 Leslie 矩阵元素 $a_{ij}$ 的**灵敏度**为：

    $$\frac{\partial \lambda_1}{\partial a_{ij}} = \frac{w_i v_j}{\mathbf{w}^T \mathbf{v}}$$

    **弹性**（相对灵敏度）为：

    $$e_{ij} = \frac{a_{ij}}{\lambda_1} \cdot \frac{\partial \lambda_1}{\partial a_{ij}} = \frac{a_{ij} w_i v_j}{\lambda_1 \mathbf{w}^T \mathbf{v}}$$

    弹性矩阵的元素之和等于 1：$\sum_{ij} e_{ij} = 1$。

!!! example "例 70.2"
    对例 70.1 中的 Leslie 矩阵，主特征值 $\lambda_1 \approx 1.489$。

    右特征向量（稳定年龄分布）：$\mathbf{v}_1 \propto (1, 0.336, 0.0564)^T$

    归一化后约为 $(0.718, 0.241, 0.041)^T$，即稳定状态下约 71.8% 的个体在第 1 龄组，24.1% 在第 2 龄组，4.1% 在第 3 龄组。

    左特征向量（繁殖价值）：$\mathbf{w}_1 \propto (1, 2.978, 4.020)^T$

    第 3 龄组个体的繁殖价值最高（约为第 1 龄组的 4 倍），因为它们马上就要繁殖。

---

## 70.3 Lotka-Volterra 与线性化

<div class="context-flow" markdown>

**核心问题**：非线性种群动力学系统（如捕食者-猎物模型）的稳定性如何通过线性化分析来判定？

</div>

!!! definition "定义 70.3 (Lotka-Volterra 捕食者-猎物模型)"
    经典 Lotka-Volterra 系统描述猎物 ($x$) 与捕食者 ($y$) 的动态：

    $$\frac{dx}{dt} = \alpha x - \beta xy = x(\alpha - \beta y)$$

    $$\frac{dy}{dt} = \delta xy - \gamma y = y(\delta x - \gamma)$$

    其中 $\alpha, \beta, \gamma, \delta > 0$ 为正参数。

!!! theorem "定理 70.4 (Lotka-Volterra 的均衡与稳定性)"
    Lotka-Volterra 系统有两个均衡点：

    1. 平凡均衡 $(0, 0)$：Jacobi 矩阵 $J_0 = \begin{pmatrix} \alpha & 0 \\ 0 & -\gamma \end{pmatrix}$，特征值 $\alpha > 0, -\gamma < 0$，是鞍点（不稳定）

    2. 非平凡均衡 $(x^*, y^*) = \left(\frac{\gamma}{\delta}, \frac{\alpha}{\beta}\right)$：

    $$J^* = \begin{pmatrix} 0 & -\frac{\beta\gamma}{\delta} \\ \frac{\delta\alpha}{\beta} & 0 \end{pmatrix}$$

    特征值 $\lambda = \pm i\sqrt{\alpha\gamma}$，纯虚数，是中心（线性分析不能确定非线性稳定性）

??? proof "证明"
    设 $f(x,y) = x(\alpha - \beta y)$，$g(x,y) = y(\delta x - \gamma)$。

    Jacobi 矩阵为：

    $$J = \begin{pmatrix} \partial f/\partial x & \partial f/\partial y \\ \partial g/\partial x & \partial g/\partial y \end{pmatrix} = \begin{pmatrix} \alpha - \beta y & -\beta x \\ \delta y & \delta x - \gamma \end{pmatrix}$$

    在 $(0,0)$：$J_0 = \begin{pmatrix} \alpha & 0 \\ 0 & -\gamma \end{pmatrix}$，特征值 $\alpha, -\gamma$，一正一负，鞍点。

    在 $(x^*, y^*) = (\gamma/\delta, \alpha/\beta)$：

    $$J^* = \begin{pmatrix} 0 & -\beta \gamma/\delta \\ \delta\alpha/\beta & 0 \end{pmatrix}$$

    特征方程 $\lambda^2 + \alpha\gamma = 0$，特征值 $\lambda = \pm i\sqrt{\alpha\gamma}$，纯虚数。

    迹 $\text{tr}(J^*) = 0$，行列式 $\det(J^*) = \alpha\gamma > 0$。线性化给出中心，但经典 Lotka-Volterra 实际上有一个守恒量 $H(x,y) = \delta x - \gamma \ln x + \beta y - \alpha \ln y$，因此轨道是封闭曲线，均衡确实是（非线性）中心。

    $\blacksquare$

!!! definition "定义 70.4 (竞争 Lotka-Volterra 模型)"
    两个物种竞争同一资源的模型：

    $$\frac{dx_1}{dt} = r_1 x_1 \left(1 - \frac{x_1 + \alpha_{12} x_2}{K_1}\right)$$

    $$\frac{dx_2}{dt} = r_2 x_2 \left(1 - \frac{\alpha_{21} x_1 + x_2}{K_2}\right)$$

    其中 $\alpha_{12}, \alpha_{21}$ 为竞争系数。共存均衡 $(x_1^*, x_2^*)$ 的稳定性由 Jacobi 矩阵的特征值决定。

!!! example "例 70.3"
    设 $\alpha = 1, \beta = 0.1, \gamma = 1.5, \delta = 0.075$。

    非平凡均衡：$x^* = 1.5/0.075 = 20$，$y^* = 1/0.1 = 10$。

    $$J^* = \begin{pmatrix} 0 & -0.1 \times 20 \\ 0.075 \times 10 & 0 \end{pmatrix} = \begin{pmatrix} 0 & -2 \\ 0.75 & 0 \end{pmatrix}$$

    特征值 $\lambda = \pm i\sqrt{1.5} \approx \pm 1.225i$。

    振荡周期 $T = 2\pi / \sqrt{1.5} \approx 5.13$ 时间单位。

---

## 70.4 房室模型

<div class="context-flow" markdown>

**核心问题**：药物在体内的吸收、分布和消除过程如何用线性微分方程组来描述？

</div>

!!! definition "定义 70.5 (房室模型)"
    **房室模型**将生物体分为 $n$ 个房室（如血液、组织、器官等），物质在房室间按一级动力学传输。系统方程为：

    $$\dot{\mathbf{x}}(t) = A\mathbf{x}(t) + \mathbf{u}(t)$$

    其中 $\mathbf{x}(t)$ 为各房室中的物质量，$\mathbf{u}(t)$ 为外部输入，$A$ 为**房室矩阵**，满足：

    - $a_{ij} \geq 0$ 对 $i \neq j$（物质从房室 $j$ 流向房室 $i$ 的速率常数）
    - $a_{ii} \leq 0$（从房室 $i$ 流出的总速率）
    - $\sum_i a_{ij} \leq 0$ 对每个 $j$（物质守恒或有外排）

!!! theorem "定理 70.5 (房室矩阵的特征值)"
    房室矩阵 $A$ 的所有特征值具有非正实部：$\text{Re}(\lambda_k) \leq 0$。若系统是封闭的（$\sum_i a_{ij} = 0$），则 $\lambda = 0$ 是特征值。若系统是开放的（至少存在外排），则所有特征值实部严格为负。

??? proof "证明"
    $A$ 的列和 $\leq 0$（即 $\mathbf{e}^T A \leq \mathbf{0}^T$，其中 $\mathbf{e} = (1,\ldots,1)^T$）。$-A$ 是 Z-矩阵（非对角元素非正）且列对角占优（$|a_{ii}| \geq \sum_{j \neq i} a_{ji}$）。

    由 Gershgorin 圆盘定理，$A$ 的特征值满足：

    $$|\lambda - a_{ii}| \leq \sum_{j \neq i} |a_{ji}| = \sum_{j \neq i} a_{ji}$$

    由于 $a_{ii} \leq 0$ 且 $|a_{ii}| \geq \sum_{j \neq i} a_{ji}$，Gershgorin 圆盘完全在左半平面（含虚轴）。

    对于开放系统，至少存在一个 $j$ 使 $\sum_i a_{ij} < 0$，此时 $-A$ 是不可约列严格对角占优矩阵（若系统连通），所有特征值实部严格为负。

    $\blacksquare$

!!! example "例 70.4"
    **二房室药代动力学模型**。药物从给药部位（房室 1）吸收到血液（房室 2），并从血液中消除。

    $$A = \begin{pmatrix} -k_a & 0 \\ k_a & -k_e \end{pmatrix}$$

    其中 $k_a$ 为吸收速率常数，$k_e$ 为消除速率常数。解为：

    $$\mathbf{x}(t) = e^{At}\mathbf{x}(0)$$

    特征值为 $\lambda_1 = -k_a$，$\lambda_2 = -k_e$。设 $k_a \neq k_e$，解为：

    $$x_1(t) = x_1(0) e^{-k_a t}$$

    $$x_2(t) = \frac{k_a x_1(0)}{k_a - k_e}(e^{-k_e t} - e^{-k_a t}) + x_2(0) e^{-k_e t}$$

    血药浓度 $x_2(t)$ 先上升（吸收阶段）后下降（消除阶段），达到峰值的时间为 $t_{\max} = \frac{\ln(k_a/k_e)}{k_a - k_e}$。

!!! example "例 70.5"
    **三房室模型**。血浆（房室 1）与快速分布组织（房室 2）和慢速分布组织（房室 3）交换，并从血浆消除。

    $$A = \begin{pmatrix} -(k_{12}+k_{13}+k_e) & k_{21} & k_{31} \\ k_{12} & -k_{21} & 0 \\ k_{13} & 0 & -k_{31} \end{pmatrix}$$

    静脉注射后，血药浓度呈三指数衰减：

    $$C(t) = A_1 e^{-\alpha t} + A_2 e^{-\beta t} + A_3 e^{-\gamma t}$$

    其中 $\alpha > \beta > \gamma > 0$ 是 $-A$ 的三个正特征值。

---

## 70.5 SEIR 流行病模型

<div class="context-flow" markdown>

**核心问题**：传染病能否在种群中传播？传播速度有多快？线性代数如何帮助计算基本再生数 $R_0$？

</div>

!!! definition "定义 70.6 (SEIR 模型)"
    SEIR 模型将种群分为四类：

    - $S$：易感者（Susceptible）
    - $E$：潜伏者（Exposed）
    - $I$：感染者（Infectious）
    - $R$：康复者（Recovered）

    动力学方程为：

    $$\frac{dS}{dt} = -\beta SI, \quad \frac{dE}{dt} = \beta SI - \sigma E, \quad \frac{dI}{dt} = \sigma E - \gamma I, \quad \frac{dR}{dt} = \gamma I$$

    其中 $\beta$ 为传播率，$\sigma$ 为潜伏期倒数，$\gamma$ 为恢复率。

!!! definition "定义 70.7 (次代矩阵与基本再生数)"
    在无病均衡点 $(S_0, 0, 0, 0)$ 附近，感染子系统 $(E, I)$ 的线性化为：

    $$\frac{d}{dt}\begin{pmatrix} E \\ I \end{pmatrix} = (F - V)\begin{pmatrix} E \\ I \end{pmatrix}$$

    其中：

    $$F = \begin{pmatrix} 0 & \beta S_0 \\ 0 & 0 \end{pmatrix} \quad (\text{新感染矩阵}), \quad V = \begin{pmatrix} \sigma & 0 \\ -\sigma & \gamma \end{pmatrix} \quad (\text{转移矩阵})$$

    **次代矩阵**（Next-Generation Matrix）定义为：

    $$K = FV^{-1}$$

    **基本再生数**为 $R_0 = \rho(K)$，即 $K$ 的谱半径。

!!! theorem "定理 70.6 (基本再生数判据)"
    无病均衡是局部渐近稳定的当且仅当 $R_0 < 1$。当 $R_0 > 1$ 时，疾病可以侵入种群。

??? proof "证明"
    感染子系统的 Jacobi 矩阵为 $J = F - V$。无病均衡稳定当且仅当 $J$ 的所有特征值实部为负。

    由次代矩阵方法（Diekmann-Heesterbeek-Metz, van den Driessche-Watmough），$F - V$ 的特征值实部全为负，当且仅当 $\rho(FV^{-1}) < 1$。

    直觉上，$K = FV^{-1}$ 的 $(i,j)$ 元素表示一个处于感染状态 $j$ 的个体在其整个感染期内产生的处于感染状态 $i$ 的新感染者数。$R_0 = \rho(K)$ 是"一代感染的放大倍数"，$R_0 < 1$ 意味着感染链最终消亡。

    更严格地，$V$ 是 M-矩阵（因为 $V$ 的非对角元素非正且 $V^{-1} \geq 0$），而 $F \geq 0$。由 M-矩阵理论，$s(F-V) < 0$（其中 $s$ 表示谱横坐标）当且仅当 $\rho(FV^{-1}) < 1$。

    $\blacksquare$

!!! example "例 70.6"
    对 SEIR 模型：

    $$V^{-1} = \begin{pmatrix} 1/\sigma & 0 \\ 1/\gamma & 1/\gamma \end{pmatrix}$$

    $$K = FV^{-1} = \begin{pmatrix} 0 & \beta S_0 \\ 0 & 0 \end{pmatrix}\begin{pmatrix} 1/\sigma & 0 \\ 1/\gamma & 1/\gamma \end{pmatrix} = \begin{pmatrix} \beta S_0/\gamma & \beta S_0/\gamma \\ 0 & 0 \end{pmatrix}$$

    $R_0 = \rho(K) = \beta S_0 / \gamma$。

    这是经典结果：$R_0 = \beta S_0 / \gamma$，即传播率乘以感染期（$1/\gamma$）乘以易感人口比例（$S_0$）。注意潜伏期 $1/\sigma$ 不影响 $R_0$（它影响疫情的时间尺度，但不影响最终是否传播）。

    例如 COVID-19 早期，若 $\beta = 0.3$/天，$\gamma = 0.1$/天（平均感染期 10 天），$S_0 \approx 1$，则 $R_0 = 3$。

---

## 70.6 序列比对评分矩阵

<div class="context-flow" markdown>

**核心问题**：比较两条蛋白质序列时，氨基酸之间的替换评分应该如何确定？评分矩阵的数学本质是什么？

</div>

!!! definition "定义 70.8 (PAM 矩阵)"
    **PAM**（Point Accepted Mutation）矩阵由 Dayhoff 等人于 1978 年提出。设 $M$ 为 $20 \times 20$ 的氨基酸突变概率矩阵（$M_{ij}$ 为氨基酸 $i$ 在一个进化时间单位内突变为氨基酸 $j$ 的条件概率），$\mathbf{q}$ 为背景氨基酸频率。则 PAM-$k$ 评分矩阵定义为：

    $$S_{ij}^{(k)} = \log_2 \frac{(M^k)_{ij}}{q_j}$$

    其中 $M^k$ 是 $M$ 的 $k$ 次方（$k$ 个进化时间单位后的突变概率），$q_j$ 是随机配对下观察到氨基酸 $j$ 的概率。

!!! theorem "定理 70.7 (评分矩阵的矩阵对数形式)"
    PAM-1 评分矩阵可以写为：

    $$S^{(1)} = \log_2 M - \mathbf{e}\log_2(\mathbf{q})^T$$

    PAM-$k$ 评分矩阵满足：

    $$S^{(k)} = \log_2 M^k - \mathbf{e}\log_2(\mathbf{q})^T$$

    其中对数逐元素取，$\mathbf{e}$ 为全 1 向量。

??? proof "证明"
    由定义，$S_{ij}^{(k)} = \log_2 (M^k)_{ij} - \log_2 q_j$，写成矩阵形式即得。

    注意 $M^k$ 的计算可以通过 $M$ 的谱分解来实现：若 $M = P \Lambda P^{-1}$，则 $M^k = P \Lambda^k P^{-1}$。

    $\blacksquare$

!!! definition "定义 70.9 (BLOSUM 矩阵)"
    **BLOSUM**（BLOcks SUbstitution Matrix）由 Henikoff 和 Henikoff 于 1992 年直接从序列比对数据库中估计，不依赖于进化模型的外推。BLOSUM-$n$ 矩阵基于序列相似度阈值为 $n\%$ 的聚类：

    $$S_{ij} = \log_2 \frac{p_{ij}}{q_i q_j}$$

    其中 $p_{ij}$ 是在比对中观察到氨基酸对 $(i,j)$ 的频率，$q_i q_j$ 是随机情况下的期望频率。

!!! example "例 70.7"
    BLOSUM62 矩阵的部分条目（以 1/2 bit 为单位）：

    |   | A  | R  | D  | C  |
    |---|----|----|----|----|
    | A |  4 | -1 | -2 |  0 |
    | R | -1 |  5 | -2 | -3 |
    | D | -2 | -2 |  6 | -3 |
    | C |  0 | -3 | -3 |  9 |

    对角线元素（自身替换）总是最高的。半胱氨酸（C）的自替换得分最高（9），因为半胱氨酸高度保守（通过二硫键固定）。

    该矩阵的正特征值和负特征值的分布决定了比对算法区分同源序列和随机序列的能力。有效的评分矩阵应该有少数几个大的正特征值（对应于保守的物理化学性质维度）。

---

## 70.7 系统发育树

<div class="context-flow" markdown>

**核心问题**：如何从物种之间的遗传距离矩阵推断进化关系（系统发育树）？

</div>

!!! definition "定义 70.10 (距离矩阵)"
    **距离矩阵** $D \in \mathbb{R}^{n \times n}$ 是对称非负矩阵，$D_{ij}$ 表示物种 $i$ 和 $j$ 之间的进化距离，$D_{ii} = 0$。若 $D$ 满足三角不等式 $D_{ij} \leq D_{ik} + D_{kj}$，则称 $D$ 为度量矩阵。

!!! definition "定义 70.11 (超度量矩阵)"
    若对所有 $i, j, k$：

    $$D_{ij} \leq \max(D_{ik}, D_{kj})$$

    则 $D$ 称为**超度量矩阵**。超度量条件等价于：对任意三个物种，三个距离中最大的两个相等。

    超度量矩阵对应于分子钟假设下的系统发育树（所有叶节点到根的距离相等）。

!!! theorem "定理 70.8 (UPGMA 算法的正确性)"
    若真实距离矩阵 $D$ 是超度量的，则 UPGMA（Unweighted Pair Group Method with Arithmetic Mean）算法能正确重建系统发育树。

??? proof "证明"
    UPGMA 算法的步骤：

    1. 找到 $D$ 中最小的非对角元素 $D_{ij}$
    2. 将物种 $i$ 和 $j$ 合并为一个簇，分支长度为 $D_{ij}/2$
    3. 更新距离矩阵：新簇与其他物种的距离取平均
    4. 重复直到所有物种合并

    正确性证明：在超度量条件下，若 $D_{ij}$ 是最小距离，则 $i$ 和 $j$ 在真实树中是最近邻（姊妹类群）。这是因为：假设 $i$ 的真正最近邻是 $k \neq j$。在超度量树中，$D_{ik} < D_{ij}$，矛盾。

    合并后的距离更新保持超度量性质（可用归纳法证明），因此递归正确。

    $\blacksquare$

!!! definition "定义 70.12 (邻接法)"
    **邻接法**（Neighbor-Joining，NJ）不要求超度量条件。它选择使总支长最小的合并对。对于每个物种 $i$，定义：

    $$u_i = \frac{1}{n-2} \sum_{k=1}^{n} D_{ik}$$

    选择使 $D_{ij} - u_i - u_j$ 最小的配对 $(i,j)$。分支长度为：

    $$v_i = \frac{1}{2}(D_{ij} + u_i - u_j), \quad v_j = D_{ij} - v_i$$

!!! example "例 70.8"
    四个物种的距离矩阵：

    $$D = \begin{pmatrix} 0 & 2 & 4 & 4 \\ 2 & 0 & 4 & 4 \\ 4 & 4 & 0 & 2 \\ 4 & 4 & 2 & 0 \end{pmatrix}$$

    验证超度量性：对任意三元组，最大的两个距离相等。例如 $(1,2,3)$：$D_{12}=2, D_{13}=4, D_{23}=4$，最大两个相等。

    UPGMA 首先合并最近的一对。$D_{12} = D_{34} = 2$ 是最小距离，合并 $(1,2)$ 和 $(3,4)$，分支长度为 1。更新距离后，两个簇之间的距离为 4，最终树为 $((1,2),(3,4))$，内部分支长度为 $4/2 - 1 = 1$。

---

## 70.8 基因调控网络

<div class="context-flow" markdown>

**核心问题**：基因之间的调控关系（激活/抑制）如何用矩阵来表示和分析？网络的稳态有什么生物学意义？

</div>

!!! definition "定义 70.13 (线性基因调控网络)"
    在简化模型中，$n$ 个基因的表达水平 $\mathbf{x}(t) = (x_1(t), \ldots, x_n(t))^T$ 满足线性微分方程：

    $$\dot{\mathbf{x}}(t) = W\mathbf{x}(t) + \mathbf{b}$$

    其中 $W \in \mathbb{R}^{n \times n}$ 为**调控矩阵**（$W_{ij} > 0$ 表示基因 $j$ 激活基因 $i$，$W_{ij} < 0$ 表示抑制），$\mathbf{b}$ 为基础表达率。

!!! theorem "定理 70.9 (基因网络的稳定性)"
    若调控矩阵 $W$ 的所有特征值具有负实部（$W$ 是 Hurwitz 稳定的），则系统有唯一稳态：

    $$\mathbf{x}^* = -W^{-1}\mathbf{b}$$

    且 $\mathbf{x}^*$ 是全局渐近稳定的。

??? proof "证明"
    稳态条件 $\dot{\mathbf{x}} = 0$ 给出 $W\mathbf{x}^* + \mathbf{b} = 0$，即 $\mathbf{x}^* = -W^{-1}\mathbf{b}$（$W$ 可逆，因为 0 不是特征值）。

    令 $\mathbf{y} = \mathbf{x} - \mathbf{x}^*$，则 $\dot{\mathbf{y}} = W\mathbf{y}$，解为 $\mathbf{y}(t) = e^{Wt}\mathbf{y}(0)$。

    由于 $W$ 的所有特征值实部为负，$\|e^{Wt}\| \to 0$（$t \to \infty$），故 $\mathbf{y}(t) \to 0$，即 $\mathbf{x}(t) \to \mathbf{x}^*$。

    $\blacksquare$

!!! definition "定义 70.14 (Boolean 网络模型)"
    在 **Boolean 网络**中，每个基因的状态为 $x_i \in \{0, 1\}$（关闭/开启），更新规则为布尔函数：

    $$x_i(t+1) = f_i(x_1(t), x_2(t), \ldots, x_n(t))$$

    网络的状态空间有 $2^n$ 个状态。网络最终会进入吸引子——固定点（稳态）或极限环（振荡）。

!!! theorem "定理 70.10 (Boolean 网络的吸引子)"
    一个 $n$ 基因的 Boolean 网络有有限的状态空间 $\{0,1\}^n$，因此其动力学最终必进入长度有限的吸引子。

    - **固定点**对应于网络的稳定表达模式（如细胞类型）
    - **极限环**对应于振荡表达模式（如细胞周期）

!!! example "例 70.9"
    一个简单的三基因抑制网络（Repressilator）：

    $$x_1(t+1) = \overline{x_3(t)}, \quad x_2(t+1) = \overline{x_1(t)}, \quad x_3(t+1) = \overline{x_2(t)}$$

    用调控矩阵表示线性化版本：

    $$W = \begin{pmatrix} -1 & 0 & -\alpha \\ -\alpha & -1 & 0 \\ 0 & -\alpha & -1 \end{pmatrix}$$

    其中 $\alpha > 0$ 为抑制强度，对角线 $-1$ 为自降解。

    $W$ 的特征值为 $\lambda_k = -1 - \alpha \omega^k$（$k = 0, 1, 2$，$\omega = e^{2\pi i/3}$）。当 $\alpha$ 足够大时，特征值可以有正实部，导致不稳定（对应于 Boolean 模型中的振荡吸引子）。

!!! example "例 70.10"
    **从表达数据推断网络**。给定 $m$ 个时间点的基因表达数据 $X \in \mathbb{R}^{n \times m}$，推断调控矩阵 $W$。

    假设 $\dot{\mathbf{x}} \approx (\mathbf{x}_{t+1} - \mathbf{x}_t)/\Delta t$，则：

    $$\frac{X_{\text{next}} - X_{\text{curr}}}{\Delta t} \approx W X_{\text{curr}} + \mathbf{b}\mathbf{e}^T$$

    令 $\dot{X} = (X_{\text{next}} - X_{\text{curr}})/\Delta t$，这是一个矩阵回归问题：

    $$\dot{X} = W X_{\text{curr}} + \mathbf{b}\mathbf{e}^T$$

    最小二乘解为 $W = \dot{X} X_{\text{curr}}^T (X_{\text{curr}} X_{\text{curr}}^T)^{-1}$（当 $\mathbf{b}$ 的影响已去除时）。

    推断出的 $W$ 的非零元素揭示了基因之间的调控关系。

---

## 本章小结

本章展示了线性代数在生物学各个层次中的应用：

1. **种群层次**：Leslie 矩阵将年龄结构种群的增长编码为矩阵迭代，Perron-Frobenius 定理给出长期增长率（主特征值）、稳定年龄分布（右特征向量）和繁殖价值（左特征向量）。

2. **生态系统层次**：Lotka-Volterra 模型的线性化分析通过 Jacobi 矩阵的特征值来判定均衡的稳定性。

3. **生理层次**：房室模型（药代动力学）是线性微分方程组 $\dot{\mathbf{x}} = A\mathbf{x}$ 的直接应用，矩阵指数 $e^{At}$ 给出解。

4. **流行病学层次**：SEIR 模型的基本再生数 $R_0$ 是次代矩阵 $K = FV^{-1}$ 的谱半径，连接了非负矩阵理论和公共卫生决策。

5. **分子层次**：序列比对的评分矩阵是突变概率矩阵的对数，系统发育树的重建依赖于距离矩阵的超度量性质，基因调控网络的分析用到了矩阵指数和谱分析。

贯穿这些应用的核心数学工具是特征值分析、Perron-Frobenius 定理、矩阵指数和 M-矩阵理论。
