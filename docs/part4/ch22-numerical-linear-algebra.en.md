# Chapter 22  Numerical Linear Algebra

<div class="context-flow" markdown>

**Prerequisites**: LU/QR/SVD/eigenvalue decompositions (Ch5-8) · **Arc**: Floating-point errors are the root cause → Condition numbers determine accuracy → Direct methods (LU, $O(n^3)$) vs iterative methods (Krylov, $O(n \cdot \text{nnz})$) → Sparsity structure dictates algorithm choice
**Essence**: Theoretical linear algebra assumes exact arithmetic — this chapter answers "how accurately and how fast can we compute?"

</div>

Numerical linear algebra is the core of scientific computing. In practical applications, matrix computations are inevitably affected by floating-point rounding errors, making algorithmic stability and efficiency critical concerns. Starting from floating-point arithmetic fundamentals, this chapter systematically introduces direct and iterative methods for linear systems, Krylov subspace methods, numerical eigenvalue computation, and sparse matrix techniques. These methods form the foundation of modern scientific computing and engineering simulation.

## 22.1 Floating-Point Arithmetic and Rounding Errors

<div class="context-flow" markdown>

**Starting point**: Real numbers → Finite-precision floating-point numbers · Machine epsilon $\epsilon_{\text{mach}} \approx 10^{-16}$ (double precision) is the atomic unit of all errors → **Catastrophic cancellation** is the most common pitfall

</div>

### IEEE 754 Floating-Point Numbers

!!! definition "Definition 22.1 (Floating-point number system)"
    A **floating-point number system** $\mathbb{F}(\beta, t, L, U)$ is determined by the following parameters:

    - $\beta$: base, typically $\beta = 2$;
    - $t$: number of mantissa digits;
    - $[L, U]$: exponent range.

    A floating-point number in this system is represented as

    $$x = \pm d_0.d_1 d_2 \cdots d_{t-1} \times \beta^e = \pm \beta^e \sum_{i=0}^{t-1} d_i \beta^{-i},$$

    where $0 \leq d_i \leq \beta - 1$, $d_0 \neq 0$ (normalized), $L \leq e \leq U$.

    The IEEE 754 standard defines two commonly used formats:

    | Format | Bits | Mantissa digits ($t$) | Exponent range |
    |--------|------|----------------------|----------------|
    | Single precision (float32) | 32 | 24 | $[-126, 127]$ |
    | Double precision (float64) | 64 | 53 | $[-1022, 1023]$ |

!!! definition "Definition 22.2 (Machine epsilon)"
    The **machine epsilon** $\epsilon_{\text{mach}}$ is defined as the smallest positive floating-point number $\epsilon$ such that $\text{fl}(1 + \epsilon) > 1$, i.e., the smallest $\epsilon$ that the floating-point system can distinguish from $1 + \epsilon$ versus $1$. Equivalently,

    $$\epsilon_{\text{mach}} = \frac{1}{2} \beta^{1-t}.$$

    For IEEE double precision, $\epsilon_{\text{mach}} \approx 1.11 \times 10^{-16}$.

    The floating-point representation $\text{fl}(x)$ of any real number $x$ satisfies

    $$\text{fl}(x) = x(1 + \delta), \quad |\delta| \leq \epsilon_{\text{mach}}.$$

!!! theorem "Theorem 22.1 (Fundamental error bound for floating-point operations)"
    Let $\oplus$ denote floating-point addition, $\ominus$ floating-point subtraction, $\otimes$ floating-point multiplication, and $\oslash$ floating-point division. For any floating-point numbers $a, b$,

    $$a \circledast b = (a * b)(1 + \delta), \quad |\delta| \leq \epsilon_{\text{mach}},$$

    where $\circledast \in \{\oplus, \ominus, \otimes, \oslash\}$ and $* \in \{+, -, \times, /\}$ is the corresponding exact operation.

??? proof "Proof"
    This is the fundamental guarantee of the IEEE 754 standard. The standard requires that the result of each floating-point operation must be the correctly rounded result (round to nearest) of the exact result, i.e., $\text{fl}(a * b)$ is the floating-point number closest to $a * b$. The relative error bound $\epsilon_{\text{mach}}$ follows directly.

!!! example "Example 22.1"
    **Catastrophic cancellation.** Consider computing $f(x) = \sqrt{x+1} - \sqrt{x}$ when $x$ is large.

    Take $x = 10^{16}$. The exact value is approximately $f(x) \approx \frac{1}{2\sqrt{x}} \approx 5 \times 10^{-9}$.

    But in double-precision floating-point, $\text{fl}(\sqrt{10^{16}+1}) = \text{fl}(\sqrt{x})$ (since $10^{16}+1$ and $10^{16}$ may be the same floating-point number in double precision), leading to $\text{fl}(f(x)) = 0$, a relative error of $100\%$.

    Using the numerically stable equivalent form $f(x) = \frac{1}{\sqrt{x+1} + \sqrt{x}}$ avoids this problem.

## 22.2 Numerical Stability

<div class="context-flow" markdown>

**Core framework**: Forward error (output deviation) $\leq$ **condition number** $\times$ backward error (equivalent input perturbation) → Backward stable algorithm + well-conditioned problem = accurate result
**Key**: $\kappa(A) = \sigma_{\max}/\sigma_{\min}$ (SVD, Ch8) — the condition number is inherent to the problem; the algorithm only controls backward error

</div>

!!! definition "Definition 22.3 (Forward error and backward error)"
    Let algorithm $\hat{f}$ be a numerical approximation of function $f$. For input $x$:

    - **Forward error**: $\|\hat{f}(x) - f(x)\|$, measuring the error in the output;
    - **Backward error**: The smallest $\|\Delta x\|$ satisfying $\hat{f}(x) = f(x + \Delta x)$, interpreting the output error as an input perturbation.

!!! definition "Definition 22.4 (Condition number)"
    The **condition number** of problem $f$ is defined as

    $$\kappa = \lim_{\delta \to 0} \sup_{\|\Delta x\| \leq \delta} \frac{\|f(x + \Delta x) - f(x)\| / \|f(x)\|}{\|\Delta x\| / \|x\|}.$$

    The condition number measures the inherent sensitivity of the problem to input perturbations. Problems with large condition numbers are called **ill-conditioned**.

    For the linear system $A\mathbf{x} = \mathbf{b}$, the condition number is

    $$\kappa(A) = \|A\| \cdot \|A^{-1}\|.$$

    Using the $2$-norm, $\kappa_2(A) = \sigma_{\max}(A) / \sigma_{\min}(A)$, where $\sigma_{\max}$ and $\sigma_{\min}$ are the largest and smallest singular values.

!!! definition "Definition 22.5 (Numerical stability)"
    - **Forward stable**: The forward error of the algorithm is small, i.e., $\|\hat{f}(x) - f(x)\| / \|f(x)\| = O(\epsilon_{\text{mach}})$.
    - **Backward stable**: The output can be interpreted as the exact output of a slightly perturbed input, i.e., $\hat{f}(x) = f(x + \Delta x)$ where $\|\Delta x\| / \|x\| = O(\epsilon_{\text{mach}})$.
    - **Mixed stable**: $\hat{f}(x) + \Delta y = f(x + \Delta x)$ where both $\|\Delta x\|/\|x\|$ and $\|\Delta y\|/\|f(x)\|$ are $O(\epsilon_{\text{mach}})$.

!!! theorem "Theorem 22.2 (Relationship between backward stability and forward error)"
    If algorithm $\hat{f}$ is backward stable, then the forward relative error satisfies

    $$\frac{\|\hat{f}(x) - f(x)\|}{\|f(x)\|} \leq \kappa(f, x) \cdot O(\epsilon_{\text{mach}}).$$

    That is, a backward stable algorithm gives accurate results for well-conditioned problems (small $\kappa$).

??? proof "Proof"
    By backward stability, $\hat{f}(x) = f(x + \Delta x)$, $\|\Delta x\|/\|x\| = O(\epsilon_{\text{mach}})$. By the definition of condition number,

    $$\frac{\|f(x + \Delta x) - f(x)\|}{\|f(x)\|} \leq \kappa(f, x) \cdot \frac{\|\Delta x\|}{\|x\|} = \kappa(f, x) \cdot O(\epsilon_{\text{mach}}).$$

!!! example "Example 22.2"
    **Examples of matrix condition numbers.** The Hilbert matrix $H_n$ with entries $H_{ij} = 1/(i+j-1)$ is a famously ill-conditioned matrix:

    | $n$ | $\kappa_2(H_n)$ |
    |-----|----------------|
    | 5   | $4.8 \times 10^5$ |
    | 10  | $1.6 \times 10^{13}$ |
    | 15  | $3.7 \times 10^{17}$ |

    When $n \geq 13$, $\kappa_2(H_n) > 1/\epsilon_{\text{mach}}$, so solving $H_n \mathbf{x} = \mathbf{b}$ in double precision may yield completely unreliable results.

## 22.3 Direct Methods for Linear Systems

<div class="context-flow" markdown>

**Direct methods**: $PA = LU$ ($\frac{2}{3}n^3$) → Partial pivoting ensures $|l_{ij}| \le 1$ → Backward stable · Growth factor $g(n)$ can theoretically reach $2^{n-1}$, but is $O(n)$ in practice

</div>

### Gaussian Elimination with Partial Pivoting

!!! definition "Definition 22.6 (Partial pivoting)"
    At step $k$ of Gaussian elimination, the **partial pivoting** strategy selects the element with the largest absolute value in column $k$ from row $k$ and below as the pivot, then swaps rows to make it the pivot row.

    Specifically, select $p = \arg\max_{k \leq i \leq n} |a^{(k)}_{ik}|$, then swap row $k$ and row $p$.

!!! theorem "Theorem 22.3 (LU decomposition)"
    Let $A$ be an $n \times n$ nonsingular matrix. Gaussian elimination with partial pivoting produces the decomposition

    $$PA = LU,$$

    where $P$ is a permutation matrix, $L$ is a unit lower triangular matrix (diagonal entries are $1$, and $|l_{ij}| \leq 1$), and $U$ is an upper triangular matrix.

    The computational cost is $\frac{2}{3}n^3 + O(n^2)$ floating-point operations.

??? proof "Proof"
    Step $k$ of Gaussian elimination ($k = 1, \ldots, n-1$):

    1. Pivot selection and row swap (encoded as permutation matrix $P_k$);
    2. For row $i$ ($i > k$), compute multiplier $l_{ik} = a^{(k)}_{ik} / a^{(k)}_{kk}$, then perform row operation $R_i \leftarrow R_i - l_{ik} R_k$.

    Since pivoting ensures $|a^{(k)}_{kk}| \geq |a^{(k)}_{ik}|$ ($\forall i > k$), we have $|l_{ik}| \leq 1$.

    The elimination process is equivalent to $M_{n-1} P_{n-1} \cdots M_1 P_1 A = U$, where $M_k$ are elimination matrices. Rearranging gives $PA = LU$.

    Operation count: Step $k$ requires $2(n-k)^2$ floating-point operations (multiplications and additions), totaling $\sum_{k=1}^{n-1} 2(n-k)^2 = \frac{2}{3}n^3 + O(n^2)$.

### Error Analysis

!!! theorem "Theorem 22.4 (Backward error of LU decomposition)"
    LU decomposition with partial pivoting is backward stable: the computed solution $\hat{\mathbf{x}}$ satisfies

    $$(A + \Delta A)\hat{\mathbf{x}} = \mathbf{b}, \quad \|\Delta A\| \leq c \, n \, g(n) \, \epsilon_{\text{mach}} \|A\|,$$

    where $c$ is a modest constant and $g(n)$ is the **growth factor**:

    $$g(n) = \frac{\max_{i,j,k} |a^{(k)}_{ij}|}{\max_{i,j} |a_{ij}|}.$$

    Theoretically $g(n)$ can reach $2^{n-1}$, but in practice it is almost always $g(n) = O(n)$.

??? proof "Proof"
    This is the result of the classical Wilkinson backward error analysis. The rounding error introduced at each elimination step can equivalently be viewed as a perturbation to the original matrix. Partial pivoting ensures $|l_{ij}| \leq 1$, preventing excessive error amplification. The growth factor $g(n)$ controls the growth of matrix entries during elimination. For a detailed step-by-step analysis, see Higham (2002).

!!! example "Example 22.3"
    Solving with Gaussian elimination with partial pivoting:

    $$\begin{pmatrix} 0.001 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} x_1 \\ x_2 \end{pmatrix} = \begin{pmatrix} 1 \\ 2 \end{pmatrix}.$$

    The exact solution is $x_1 \approx 1.001$, $x_2 \approx 0.999$.

    **Without pivoting**: Using $0.001$ as pivot, multiplier $l_{21} = 1000$. In finite precision, the large multiplier causes severe rounding errors.

    **With partial pivoting**: After swapping rows, pivot is $1$, multiplier $l_{21} = 0.001$, and the computation is stable.

## 22.4 Iterative Methods for Linear Systems

<div class="context-flow" markdown>

**Turning point**: Direct methods at $O(n^3)$ are infeasible for large sparse systems → Iterative methods only require matrix-vector products → Convergence is determined by the **spectral radius** $\rho(G) < 1$ (eigenvalue theory, Ch6)
**Hierarchy**: Jacobi/GS (classical) → SOR (accelerated) → Krylov (optimal polynomial approximation, S22.5)

</div>

For large sparse systems, the $O(n^3)$ complexity of direct methods is prohibitive, making iterative methods the better choice.

### Basic Iterative Methods

!!! definition "Definition 22.7 (Matrix splitting)"
    Let $A = M - N$ where $M$ is invertible. Then the linear system $A\mathbf{x} = \mathbf{b}$ is equivalent to

    $$\mathbf{x} = M^{-1}N\mathbf{x} + M^{-1}\mathbf{b}.$$

    This gives the iteration $\mathbf{x}^{(k+1)} = M^{-1}N\mathbf{x}^{(k)} + M^{-1}\mathbf{b}$. The matrix $G = M^{-1}N$ is called the **iteration matrix**.

    Common splittings: Decompose $A = D - L - U$ ($D$ is the diagonal part, $-L$ is the strict lower triangular part, $-U$ is the strict upper triangular part):

    | Method | $M$ | Iteration matrix $G$ |
    |--------|-----|---------------------|
    | **Jacobi iteration** | $D$ | $D^{-1}(L + U)$ |
    | **Gauss-Seidel iteration** | $D - L$ | $(D - L)^{-1}U$ |
    | **SOR ($\omega$)** | $\frac{1}{\omega}D - L$ | $(\frac{1}{\omega}D - L)^{-1}[(\frac{1}{\omega}-1)D + U]$ |

!!! theorem "Theorem 22.5 (Necessary and sufficient condition for convergence)"
    The iteration $\mathbf{x}^{(k+1)} = G\mathbf{x}^{(k)} + \mathbf{c}$ converges to the fixed point for any initial vector $\mathbf{x}^{(0)}$ if and only if

    $$\rho(G) < 1,$$

    where $\rho(G) = \max_i |\lambda_i(G)|$ is the **spectral radius** of $G$.

??? proof "Proof"
    Let the fixed point be $\mathbf{x}^* = G\mathbf{x}^* + \mathbf{c}$. Let $\mathbf{e}^{(k)} = \mathbf{x}^{(k)} - \mathbf{x}^*$, then $\mathbf{e}^{(k+1)} = G\mathbf{e}^{(k)}$, so $\mathbf{e}^{(k)} = G^k \mathbf{e}^{(0)}$.

    $G^k \to 0$ ($k \to \infty$) for all $\mathbf{e}^{(0)}$ if and only if $\rho(G) < 1$.

    **Necessity**: If $\rho(G) \geq 1$, take $\mathbf{e}^{(0)}$ to be the eigenvector corresponding to the eigenvalue $\lambda$ of largest modulus, then $\|G^k \mathbf{e}^{(0)}\| = |\lambda|^k \|\mathbf{e}^{(0)}\| \not\to 0$.

    **Sufficiency**: If $\rho(G) < 1$, by Gelfand's formula $\rho(G) = \lim_{k \to \infty} \|G^k\|^{1/k}$, for any $\rho(G) < r < 1$, there exists $C > 0$ such that $\|G^k\| \leq C r^k$, so $\|G^k\| \to 0$.

!!! theorem "Theorem 22.6 (Convergence for diagonally dominant matrices)"
    If $A$ is **strictly diagonally dominant**, i.e.,

    $$|a_{ii}| > \sum_{j \neq i} |a_{ij}|, \quad \forall i,$$

    then both Jacobi and Gauss-Seidel iterations converge.

??? proof "Proof"
    For Jacobi iteration, let $G_J = D^{-1}(L + U)$. By Gershgorin's disc theorem, each eigenvalue $\lambda$ of $G_J$ satisfies

    $$|\lambda - 0| \leq \sum_{j \neq i} \left|\frac{a_{ij}}{a_{ii}}\right| < 1$$

    (using the strict diagonal dominance condition). Hence $\rho(G_J) < 1$.

    The convergence proof for Gauss-Seidel is more technical and can use the Stein-Rosenberg theorem or a direct Lyapunov function construction.

!!! example "Example 22.4"
    Consider the system $A\mathbf{x} = \mathbf{b}$ with

    $$A = \begin{pmatrix} 4 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 4 \end{pmatrix}, \quad \mathbf{b} = \begin{pmatrix} 1 \\ 5 \\ 0 \end{pmatrix}.$$

    $A$ is strictly diagonally dominant. The Jacobi iteration matrix is

    $$G_J = \begin{pmatrix} 0 & 1/4 & 0 \\ 1/4 & 0 & 1/4 \\ 0 & 1/4 & 0 \end{pmatrix}, \quad \rho(G_J) = \frac{1}{2\sqrt{2}} \approx 0.354.$$

    The spectral radius of the Gauss-Seidel iteration matrix is $\rho(G_{GS}) = \rho(G_J)^2 \approx 0.125$ (the Stein-Rosenberg relation for symmetric positive definite case). Hence Gauss-Seidel converges faster.

## 22.5 Krylov Subspace Methods

<div class="context-flow" markdown>

**Core idea**: $\mathcal{K}_k(A, \mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}$ = the maximal information space built from $k$ matrix-vector products → Arnoldi orthogonalization → Lanczos (tridiagonalization in the symmetric case)
**Unified framework**: CG (S22.6) = symmetric positive definite Krylov optimization · GMRES (S22.7) = general case residual minimization

</div>

### Krylov Subspace

!!! definition "Definition 22.8 (Krylov subspace)"
    Given an $n \times n$ matrix $A$ and vector $\mathbf{b}$, the $k$-th **Krylov subspace** is defined as

    $$\mathcal{K}_k(A, \mathbf{b}) = \text{span}\{\mathbf{b}, A\mathbf{b}, A^2\mathbf{b}, \ldots, A^{k-1}\mathbf{b}\}.$$

    Krylov subspaces satisfy the nesting relation $\mathcal{K}_1 \subseteq \mathcal{K}_2 \subseteq \cdots$, and there exists a smallest $m \leq n$ such that $\mathcal{K}_m = \mathcal{K}_{m+1} = \cdots$.

!!! proposition "Proposition 22.1 (Properties of Krylov subspaces)"
    1. $\dim \mathcal{K}_k(A, \mathbf{b}) \leq k$, with equality not necessarily holding.
    2. Let $p(\lambda)$ be the minimal polynomial of $A$, $\deg p = d$. Then $\dim \mathcal{K}_k(A, \mathbf{b}) = \min(k, d_{\mathbf{b}})$, where $d_{\mathbf{b}}$ is the smallest integer such that $\mathcal{K}_{d_{\mathbf{b}}}(A, \mathbf{b}) = \mathcal{K}_{d_{\mathbf{b}}+1}(A, \mathbf{b})$.
    3. $\mathcal{K}_k(A, \mathbf{b}) = \{p(A)\mathbf{b} : p \text{ is a polynomial of degree at most } k-1\}$.

??? proof "Proof"
    Property 3: $\mathcal{K}_k(A, \mathbf{b})$ is spanned by $\mathbf{b}, A\mathbf{b}, \ldots, A^{k-1}\mathbf{b}$, and $p(A)\mathbf{b} = c_0\mathbf{b} + c_1 A\mathbf{b} + \cdots + c_{k-1}A^{k-1}\mathbf{b}$ is exactly a linear combination of these. Conversely, any such linear combination can be written as $p(A)\mathbf{b}$.

### Arnoldi Process

!!! definition "Definition 22.9 (Arnoldi process)"
    The **Arnoldi process** is an algorithm for constructing an orthonormal basis of $\mathcal{K}_k(A, \mathbf{b})$:

    **Input**: Matrix $A$, vector $\mathbf{b}$, number of steps $k$.

    1. $\mathbf{q}_1 = \mathbf{b} / \|\mathbf{b}\|$
    2. **For** $j = 1, 2, \ldots, k$:
        - $\mathbf{v} = A\mathbf{q}_j$
        - **For** $i = 1, \ldots, j$: $h_{ij} = \mathbf{q}_i^T \mathbf{v}$, $\mathbf{v} = \mathbf{v} - h_{ij}\mathbf{q}_i$
        - $h_{j+1,j} = \|\mathbf{v}\|$
        - If $h_{j+1,j} = 0$ then stop (Krylov subspace is invariant)
        - $\mathbf{q}_{j+1} = \mathbf{v} / h_{j+1,j}$

    **Output**: Orthogonal matrix $Q_k = [\mathbf{q}_1, \ldots, \mathbf{q}_k]$ and upper Hessenberg matrix $H_k = (h_{ij})_{k \times k}$.

!!! theorem "Theorem 22.7 (Arnoldi relation)"
    The matrices produced by the Arnoldi process satisfy

    $$AQ_k = Q_k H_k + h_{k+1,k} \, \mathbf{q}_{k+1} \mathbf{e}_k^T = Q_{k+1} \tilde{H}_k,$$

    where $\tilde{H}_k$ is a $(k+1) \times k$ matrix. Equivalently, $Q_k^T A Q_k = H_k$.

??? proof "Proof"
    From step $j$ of the Arnoldi process:

    $$A\mathbf{q}_j = \sum_{i=1}^{j} h_{ij} \mathbf{q}_i + h_{j+1,j} \mathbf{q}_{j+1}.$$

    Combining the equations for $j = 1, \ldots, k$ gives $AQ_k = Q_k H_k + h_{k+1,k} \mathbf{q}_{k+1} \mathbf{e}_k^T$.

    By column orthogonality of $Q_k$, left-multiplying by $Q_k^T$ gives $Q_k^T A Q_k = H_k$ (using $Q_k^T \mathbf{q}_{k+1} = \mathbf{0}$).

### Lanczos Process

!!! theorem "Theorem 22.8 (Lanczos process)"
    When $A$ is a **symmetric matrix**, $H_k$ in the Arnoldi process reduces to a **tridiagonal matrix** $T_k$, and the Arnoldi process simplifies to the **Lanczos process**:

    $$\beta_{j+1} \mathbf{q}_{j+1} = A\mathbf{q}_j - \alpha_j \mathbf{q}_j - \beta_j \mathbf{q}_{j-1},$$

    where $\alpha_j = \mathbf{q}_j^T A \mathbf{q}_j$, $\beta_{j+1} = \|A\mathbf{q}_j - \alpha_j \mathbf{q}_j - \beta_j \mathbf{q}_{j-1}\|$.

    Each step only requires storing the current and previous vectors and one matrix-vector product.

??? proof "Proof"
    When $A = A^T$, $H_k = Q_k^T A Q_k$ is also symmetric. Since $H_k$ is an upper Hessenberg matrix and symmetric, it must be tridiagonal. Setting $\alpha_j = h_{jj}$, $\beta_j = h_{j,j-1} = h_{j-1,j}$, the Arnoldi recurrence simplifies to a three-term recurrence.

!!! example "Example 22.5"
    For the symmetric matrix $A = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix}$ and $\mathbf{b} = \begin{pmatrix} 1 \\ 0 \\ 0 \end{pmatrix}$, the Lanczos process gives:

    - $\mathbf{q}_1 = (1, 0, 0)^T$, $\alpha_1 = \mathbf{q}_1^T A \mathbf{q}_1 = 2$
    - $\mathbf{r}_1 = A\mathbf{q}_1 - 2\mathbf{q}_1 = (0, 1, 0)^T$, $\beta_2 = 1$, $\mathbf{q}_2 = (0, 1, 0)^T$
    - $\alpha_2 = \mathbf{q}_2^T A \mathbf{q}_2 = 3$
    - $\mathbf{r}_2 = A\mathbf{q}_2 - 3\mathbf{q}_2 - 1 \cdot \mathbf{q}_1 = (0, 0, 1)^T$, $\beta_3 = 1$, $\mathbf{q}_3 = (0, 0, 1)^T$

    The resulting tridiagonal matrix $T_3 = \begin{pmatrix} 2 & 1 & 0 \\ 1 & 3 & 1 \\ 0 & 1 & 2 \end{pmatrix} = A$ (exact recovery).

## 22.6 Conjugate Gradient Method

<div class="context-flow" markdown>

**Essence of CG**: Find the $A$-norm best approximation in the Krylov subspace → Convergence rate $O(\sqrt{\kappa})$ (Chebyshev polynomial optimality) → **Preconditioning** to reduce effective $\kappa$ is the key in practice

</div>

### Algorithm Derivation

!!! definition "Definition 22.10 ($A$-conjugacy)"
    Let $A$ be a symmetric positive definite matrix. Vectors $\mathbf{p}$ and $\mathbf{q}$ are called **$A$-conjugate** if $\mathbf{p}^T A \mathbf{q} = 0$.

The **Conjugate Gradient method** (CG) solves the symmetric positive definite system $A\mathbf{x} = \mathbf{b}$, equivalently minimizing the quadratic function

$$\phi(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T A \mathbf{x} - \mathbf{b}^T \mathbf{x}.$$

!!! theorem "Theorem 22.9 (Optimality of CG)"
    The CG iterate $\mathbf{x}_k$ at step $k$ satisfies

    $$\mathbf{x}_k = \arg\min_{\mathbf{x} \in \mathbf{x}_0 + \mathcal{K}_k(A, \mathbf{r}_0)} \|\mathbf{x} - \mathbf{x}^*\|_A,$$

    where $\mathbf{r}_0 = \mathbf{b} - A\mathbf{x}_0$ is the initial residual, $\mathbf{x}^*$ is the exact solution, and $\|\mathbf{v}\|_A = \sqrt{\mathbf{v}^T A \mathbf{v}}$ is the $A$-norm.

    Equivalently, CG finds the best approximation in the $A$-norm within the Krylov subspace.

??? proof "Proof"
    CG constructs $A$-conjugate search directions $\{\mathbf{p}_0, \mathbf{p}_1, \ldots\}$. It can be shown that $\text{span}\{\mathbf{p}_0, \ldots, \mathbf{p}_{k-1}\} = \mathcal{K}_k(A, \mathbf{r}_0)$.

    Since $\mathbf{x}_k = \mathbf{x}_0 + \sum_{i=0}^{k-1} \alpha_i \mathbf{p}_i$, and the $\alpha_i$ are chosen so that the error $\mathbf{e}_k = \mathbf{x}_k - \mathbf{x}^*$ satisfies $\mathbf{e}_k^T A \mathbf{p}_i = 0$ ($i = 0, \ldots, k-1$), this is exactly the orthogonal projection in $\mathbf{x}_0 + \mathcal{K}_k$ with respect to the $A$-inner product, equivalent to $A$-norm minimization.

### Convergence Analysis

<div class="context-flow" markdown>

**Insight**: CG convergence rate $\propto \sqrt{\kappa}$ rather than $\kappa$ — because the optimal approximation property of Chebyshev polynomials on $[\lambda_{\min}, \lambda_{\max}]$ reduces the number of iterations from $O(\kappa)$ to $O(\sqrt{\kappa})$

</div>

!!! theorem "Theorem 22.10 (Convergence rate of CG)"
    Let $\kappa = \kappa_2(A) = \lambda_{\max}/\lambda_{\min}$ be the condition number of $A$. The CG error satisfies

    $$\|\mathbf{x}_k - \mathbf{x}^*\|_A \leq 2 \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|\mathbf{x}_0 - \mathbf{x}^*\|_A.$$

    Therefore, the number of iterations needed to reduce the error by a factor of $\epsilon$ is $O(\sqrt{\kappa} \log(1/\epsilon))$.

??? proof "Proof"
    By the optimality of CG,

    $$\|\mathbf{e}_k\|_A = \min_{p \in \mathcal{P}_k, \, p(0)=1} \max_{\lambda \in \sigma(A)} |p(\lambda)| \cdot \|\mathbf{e}_0\|_A,$$

    where $\mathcal{P}_k$ is the set of polynomials of degree at most $k$. Choosing the Chebyshev polynomial

    $$p_k(\lambda) = \frac{T_k\left(\frac{\lambda_{\max}+\lambda_{\min}-2\lambda}{\lambda_{\max}-\lambda_{\min}}\right)}{T_k\left(\frac{\lambda_{\max}+\lambda_{\min}}{\lambda_{\max}-\lambda_{\min}}\right)},$$

    the estimate follows from properties of $T_k$.

!!! example "Example 22.6"
    For a symmetric positive definite system with $\kappa = 100$, CG reduces the error to $10^{-6}$ in approximately

    $$k \approx \frac{\sqrt{100}}{2} \ln(2/10^{-6}) \approx 5 \times 14.5 \approx 73$$

    iterations. Using a good preconditioner $M$ that reduces the effective condition number to $\kappa_{\text{eff}} = 10$ requires only about $23$ steps.

### Preconditioned Conjugate Gradient Method

!!! definition "Definition 22.11 (Preconditioned conjugate gradient method)"
    The **Preconditioned CG** (PCG) method applies CG to the system $M^{-1}A\mathbf{x} = M^{-1}\mathbf{b}$, where $M$ is the **preconditioner**, satisfying $M \approx A$ and $M^{-1}$ is easy to compute.

    The core steps of the PCG algorithm are:

    1. $\mathbf{r}_k = \mathbf{b} - A\mathbf{x}_k$
    2. Solve $M\mathbf{z}_k = \mathbf{r}_k$
    3. $\beta_k = \mathbf{r}_k^T \mathbf{z}_k / \mathbf{r}_{k-1}^T \mathbf{z}_{k-1}$
    4. $\mathbf{p}_k = \mathbf{z}_k + \beta_k \mathbf{p}_{k-1}$
    5. $\alpha_k = \mathbf{r}_k^T \mathbf{z}_k / (\mathbf{p}_k^T A \mathbf{p}_k)$
    6. $\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k$

    Common preconditioners include: incomplete Cholesky factorization, symmetric Gauss-Seidel, multigrid methods, etc.

## 22.7 GMRES Method

<div class="context-flow" markdown>

**Generalization**: CG is limited to symmetric positive definite → GMRES handles general square matrices · The Arnoldi relation converts the large problem into a small Hessenberg least-squares problem → Convergence depends on the clustering of eigenvalues

</div>

!!! definition "Definition 22.12 (GMRES method)"
    **GMRES** (Generalized Minimal Residual) is used to solve general (possibly nonsymmetric) linear systems $A\mathbf{x} = \mathbf{b}$. At step $k$, GMRES finds the approximation in the Krylov subspace that minimizes the residual $2$-norm:

    $$\mathbf{x}_k = \arg\min_{\mathbf{x} \in \mathbf{x}_0 + \mathcal{K}_k(A, \mathbf{r}_0)} \|\mathbf{b} - A\mathbf{x}\|_2.$$

!!! theorem "Theorem 22.11 (Relationship between GMRES and Arnoldi)"
    GMRES step $k$ is equivalent to solving the least-squares problem

    $$\min_{\mathbf{y} \in \mathbb{R}^k} \|\beta \mathbf{e}_1 - \tilde{H}_k \mathbf{y}\|_2,$$

    where $\beta = \|\mathbf{r}_0\|$ and $\tilde{H}_k$ is the $(k+1) \times k$ upper Hessenberg matrix from the Arnoldi process. The solution is $\mathbf{x}_k = \mathbf{x}_0 + Q_k \mathbf{y}_k$.

??? proof "Proof"
    By the Arnoldi relation $AQ_k = Q_{k+1}\tilde{H}_k$. For $\mathbf{x} = \mathbf{x}_0 + Q_k\mathbf{y}$,

    $$\mathbf{b} - A\mathbf{x} = \mathbf{r}_0 - AQ_k\mathbf{y} = \beta Q_{k+1}\mathbf{e}_1 - Q_{k+1}\tilde{H}_k\mathbf{y} = Q_{k+1}(\beta\mathbf{e}_1 - \tilde{H}_k\mathbf{y}).$$

    By column orthogonality of $Q_{k+1}$,

    $$\|\mathbf{b} - A\mathbf{x}\|_2 = \|\beta\mathbf{e}_1 - \tilde{H}_k\mathbf{y}\|_2.$$

    Therefore the GMRES optimization reduces to a small-scale ($(k+1) \times k$) least-squares problem.

!!! proposition "Proposition 22.2 (Convergence of GMRES)"
    1. GMRES converges to the exact solution in at most $n$ steps (in exact arithmetic).
    2. If the eigenvalues of $A$ are far from zero and clustered, GMRES converges fast. Specifically, if $A$ is diagonalizable, $A = X \Lambda X^{-1}$, then

    $$\frac{\|\mathbf{r}_k\|}{\|\mathbf{r}_0\|} \leq \kappa(X) \min_{p \in \mathcal{P}_k, \, p(0)=1} \max_{\lambda \in \sigma(A)} |p(\lambda)|.$$

??? proof "Proof"
    Property 1: $\mathcal{K}_n(A, \mathbf{r}_0) = \mathbb{R}^n$ (when $A$ is nonsingular), so $\mathbf{x}^* \in \mathbf{x}_0 + \mathcal{K}_n$.

    Property 2: By the optimality of GMRES, $\|\mathbf{r}_k\| \leq \|p_k(A)\mathbf{r}_0\|$ for any $k$-th degree polynomial $p_k$ with $p_k(0) = 1$. If $A = X\Lambda X^{-1}$, then $\|p_k(A)\| \leq \|X\| \cdot \|p_k(\Lambda)\| \cdot \|X^{-1}\| = \kappa(X) \max_\lambda |p_k(\lambda)|$.

!!! example "Example 22.7"
    For the upper triangular matrix $A = I + N$ ($N$ strictly upper triangular, $N^n = 0$), GMRES converges exactly in at most $n$ steps (since the minimal polynomial of $A$ has degree at most $n$). In fact, if $N^m = 0$ ($m < n$), then GMRES converges in $m$ steps.

## 22.8 Numerical Eigenvalue Computation

<div class="context-flow" markdown>

**Hierarchy**: Power iteration (simplest, $|\lambda_2/\lambda_1|^k$) → Inverse iteration + Rayleigh quotient (cubic convergence) → QR iteration (full spectrum, first reduce to Hessenberg then iterate)
**Link**: QR decomposition (Ch8) is not a decomposition tool here but an iteration engine · Ch25 Rayleigh quotient optimization perspective

</div>

### Power Iteration

!!! definition "Definition 22.13 (Power iteration)"
    **Power iteration** computes the eigenvalue of largest modulus and the corresponding eigenvector of matrix $A$:

    1. Choose initial vector $\mathbf{q}_0$ ($\|\mathbf{q}_0\| = 1$)
    2. **For** $k = 0, 1, 2, \ldots$:
        - $\mathbf{z}_{k+1} = A\mathbf{q}_k$
        - $\mathbf{q}_{k+1} = \mathbf{z}_{k+1} / \|\mathbf{z}_{k+1}\|$
        - $\lambda_{k+1} = \mathbf{q}_{k+1}^T A \mathbf{q}_{k+1}$ (Rayleigh quotient)

!!! theorem "Theorem 22.12 (Convergence of power iteration)"
    Let $A$ have eigenvalues $|\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|$, and suppose the initial vector $\mathbf{q}_0$ has a nonzero component in the direction of the eigenvector of $\lambda_1$. Then

    $$|\lambda_{k} - \lambda_1| = O\left(\left|\frac{\lambda_2}{\lambda_1}\right|^{2k}\right), \quad \text{dist}(\mathbf{q}_k, \text{span}\{\mathbf{v}_1\}) = O\left(\left|\frac{\lambda_2}{\lambda_1}\right|^k\right).$$

    When using the Rayleigh quotient, the eigenvalue convergence rate is twice that of the eigenvector.

??? proof "Proof"
    Let $\mathbf{q}_0 = \sum_{i=1}^n c_i \mathbf{v}_i$, $c_1 \neq 0$. Then

    $$A^k \mathbf{q}_0 = \sum_{i=1}^n c_i \lambda_i^k \mathbf{v}_i = c_1 \lambda_1^k \left[\mathbf{v}_1 + \sum_{i=2}^n \frac{c_i}{c_1}\left(\frac{\lambda_i}{\lambda_1}\right)^k \mathbf{v}_i\right].$$

    Since $|\lambda_i/\lambda_1| < 1$ ($i \geq 2$), as $k \to \infty$, $A^k \mathbf{q}_0 / \|A^k \mathbf{q}_0\| \to \mathbf{v}_1/\|\mathbf{v}_1\|$ at rate $|\lambda_2/\lambda_1|^k$.

    The Rayleigh quotient error $\lambda_k = \mathbf{q}_k^T A \mathbf{q}_k$ is $O(\|\mathbf{q}_k - \mathbf{v}_1\|^2) = O(|\lambda_2/\lambda_1|^{2k})$.

### QR Algorithm

!!! definition "Definition 22.14 (QR iteration)"
    **Basic QR iteration**:

    1. Set $A_0 = A$
    2. **For** $k = 0, 1, 2, \ldots$:
        - Compute QR decomposition $A_k = Q_k R_k$
        - Set $A_{k+1} = R_k Q_k$

    **QR iteration with shifts**:

    1. Set $A_0 = A$
    2. **For** $k = 0, 1, 2, \ldots$:
        - Choose shift $\mu_k$ (typically the Wilkinson shift of $A_k$)
        - Compute QR decomposition $A_k - \mu_k I = Q_k R_k$
        - Set $A_{k+1} = R_k Q_k + \mu_k I$

!!! theorem "Theorem 22.13 (Convergence of QR iteration)"
    Let the eigenvalues of $A$ satisfy $|\lambda_1| > |\lambda_2| > \cdots > |\lambda_n|$ (distinct moduli). Then in basic QR iteration, $A_k$ converges to an upper triangular matrix with diagonal entries converging to $\lambda_1, \ldots, \lambda_n$.

    Convergence rate: The $(i,j)$ entry of $A_k$ ($i > j$) tends to zero at rate $|\lambda_i/\lambda_j|^k$.

    QR iteration with Wilkinson shifts achieves **cubic convergence** (for symmetric matrices).

??? proof "Proof"
    Note that $A_{k+1} = R_k Q_k = Q_k^T (Q_k R_k) Q_k = Q_k^T A_k Q_k$, i.e., $A_{k+1}$ is orthogonally similar to $A_k$. Let $\hat{Q}_k = Q_0 Q_1 \cdots Q_{k-1}$, then $A_k = \hat{Q}_k^T A \hat{Q}_k$.

    On the other hand, $A^k = (Q_0 R_0)(Q_1 R_1) \cdots = \hat{Q}_k \hat{R}_k$ (where $\hat{R}_k = R_0 R_1 \cdots R_{k-1}$). This is closely related to the QR decomposition of $A^k$. Using the Schur decomposition $A = USU^*$ and asymptotic analysis of the QR decomposition of $A^k = US^kU^*$, one can prove that the lower triangular part of $A_k$ tends to zero.

### Hessenberg Reduction

!!! theorem "Theorem 22.14 (Hessenberg preprocessing)"
    Any matrix $A$ can be transformed to upper Hessenberg form $H = Q^T A Q$ (where $h_{ij} = 0$ for $i > j + 1$) by an orthogonal similarity transformation using $O(\frac{2}{3}n^3)$ operations.

    Performing one step of QR iteration on a Hessenberg matrix requires only $O(n^2)$ operations (instead of $O(n^3)$), making the QR algorithm practically efficient.

??? proof "Proof"
    Use Householder transformations to eliminate column by column. At step $k$, choose a Householder matrix $P_k$ to zero out entries below row $k+2$ in column $k$, performing the similarity transformation $A \leftarrow P_k A P_k$. Maintaining similarity means multiplying by $P_k$ on the right as well, which introduces nonzero elements to the right of column $k$ but does not destroy the Hessenberg structure of previous columns.

    QR decomposition of a Hessenberg matrix requires only $n-1$ Givens rotations (each $O(n)$), totaling $O(n^2)$.

!!! example "Example 22.8"
    Transform the matrix $A = \begin{pmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 0 \end{pmatrix}$ to Hessenberg form first, then use QR iteration for eigenvalues.

    Hessenberg reduction: Use a Householder transformation to eliminate $a_{31}$. Choose $P_1$ to map $(a_{21}, a_{31})^T = (4, 7)^T$ to $(\pm\sqrt{65}, 0)^T$, then

    $$H = P^T A P = \begin{pmatrix} 1 & * & * \\ -\sqrt{65} & * & * \\ 0 & * & * \end{pmatrix}$$

    is already in Hessenberg form (a $3 \times 3$ matrix requires only one step). Subsequent QR iterations operate on this Hessenberg matrix, each step costing only $O(n^2)$ operations.

## 22.9 Sparse Matrices

<div class="context-flow" markdown>

**Practice**: In scientific computing $n \sim 10^6$ but $\text{nnz} \sim O(n)$ (e.g., finite elements/differences) → COO/CSR/CSC formats → Matrix-vector product $O(\text{nnz})$ makes Krylov methods feasible

</div>

### Storage Formats

!!! definition "Definition 22.15 (Sparse matrix storage formats)"
    A **sparse matrix** is a matrix in which most elements are zero. Let an $n \times n$ matrix have $\text{nnz}$ nonzero elements, $\text{nnz} \ll n^2$. Common storage formats:

    - **COO (Coordinate format)**: Stores three arrays: row indices `row[]`, column indices `col[]`, values `val[]`, each of length $\text{nnz}$.
    - **CSR (Compressed Sparse Row)**: Stores arrays `val[]` (length $\text{nnz}$), `col_idx[]` (length $\text{nnz}$), `row_ptr[]` (length $n+1$). Here `row_ptr[i]` to `row_ptr[i+1]-1` indicate the positions of nonzero elements of row $i$ in `val[]` and `col_idx[]`.
    - **CSC (Compressed Sparse Column)**: Similar to CSR but compressed by column.

    | Format | Storage | Matrix-vector product | Access $a_{ij}$ |
    |--------|---------|----------------------|-----------------|
    | COO  | $3 \cdot \text{nnz}$ | $O(\text{nnz})$ | $O(\text{nnz})$ |
    | CSR  | $2 \cdot \text{nnz} + n + 1$ | $O(\text{nnz})$ (row-major) | $O(\log d_i)$ |
    | CSC  | $2 \cdot \text{nnz} + n + 1$ | $O(\text{nnz})$ (column-major) | $O(\log d_j)$ |

    where $d_i$ is the number of nonzero elements in row $i$.

!!! example "Example 22.9"
    Consider the matrix

    $$A = \begin{pmatrix} 5 & 0 & 0 & 3 \\ 0 & 8 & 0 & 0 \\ 0 & 0 & 3 & 0 \\ 0 & 6 & 0 & 1 \end{pmatrix}.$$

    **COO format**:

    - `val = [5, 3, 8, 3, 6, 1]`
    - `row = [0, 0, 1, 2, 3, 3]`
    - `col = [0, 3, 1, 2, 1, 3]`

    **CSR format**:

    - `val = [5, 3, 8, 3, 6, 1]`
    - `col_idx = [0, 3, 1, 2, 1, 3]`
    - `row_ptr = [0, 2, 3, 4, 6]`

    Storage: $\text{nnz} = 6$, COO requires $18$ numbers, CSR requires $17$ numbers, while dense storage requires $16$ numbers. When $n$ is large and $\text{nnz} \ll n^2$, the advantage of sparse storage is significant.

!!! proposition "Proposition 22.3 (Efficiency of sparse matrix-vector multiplication)"
    The algorithm for matrix-vector multiplication $\mathbf{y} = A\mathbf{x}$ in CSR format is:

    **For** $i = 0, \ldots, n-1$:

    $$y_i = \sum_{k=\text{row\_ptr}[i]}^{\text{row\_ptr}[i+1]-1} \text{val}[k] \cdot x_{\text{col\_idx}[k]}.$$

    The computational cost is $2 \cdot \text{nnz}$ floating-point operations, while dense matrix-vector multiplication requires $2n^2$. For structured sparse matrices such as five-point stencil matrices, $\text{nnz} = O(n)$, so sparse methods are far more efficient than dense methods.

??? proof "Proof"
    By the definition of CSR format, the nonzero elements of row $i$ are stored in `val[row_ptr[i]:row_ptr[i+1]]`, with corresponding column indices `col_idx[row_ptr[i]:row_ptr[i+1]]`. The matrix-vector product $y_i = \sum_j a_{ij} x_j$ only traverses the nonzero $a_{ij}$, so the total operation count is $\sum_i (\text{number of nonzero elements in row } i) \times 2 = 2 \cdot \text{nnz}$.

!!! example "Example 22.10"
    **Five-point stencil.** Discretizing the Laplace equation $-\Delta u = f$ on an $N \times N$ grid yields an $n \times n$ ($n = N^2$) sparse matrix with at most $5$ nonzero elements per row:

    $$A = \begin{pmatrix} T & -I & & \\ -I & T & -I & \\ & \ddots & \ddots & \ddots \\ & & -I & T \end{pmatrix}, \quad T = \begin{pmatrix} 4 & -1 & & \\ -1 & 4 & -1 & \\ & \ddots & \ddots & \ddots \\ & & -1 & 4 \end{pmatrix}.$$

    Here $\text{nnz} \leq 5n$. Using CSR storage and sparse matrix-vector multiplication enables efficient implementation of iterative methods like CG. For $N = 1000$ ($n = 10^6$), dense storage requires $10^{12}$ floating-point numbers (about $8$ TB), while sparse storage requires only about $5 \times 10^6$ (about $40$ MB).
