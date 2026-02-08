# Chapter 25  Linear Algebra in Optimization

<div class="context-flow" markdown>

**Prerequisites**: SVD/eigenvalues (Ch6-8) · positive definiteness (Ch7) · manifold optimization (Ch24) · **Chapter arc**: LP (basis = column selection) → least squares (QR/SVD) → SDP (positive semidefinite cone) → matrix completion/compressed sensing (nuclear norm/$\ell_1$) → PCA/Rayleigh quotient
**Essence**: The three major decompositions of linear algebra (LU/QR/SVD) and eigenvalue theory form the computational backbone of modern optimization

</div>

Optimization theory and linear algebra share deep and extensive connections. Linear algebra provides not only the modeling language and analytical tools for optimization problems, but its core decomposition methods (such as QR decomposition, SVD, and eigenvalue decomposition) directly constitute the computational backbone of many optimization algorithms. This chapter systematically demonstrates the central role of linear algebra in optimization, covering linear programming, least squares, semidefinite programming, matrix completion, compressed sensing, principal component analysis, low-rank approximation, and eigenvalue optimization.

---

## 25.1 Linear Programming Fundamentals

<div class="context-flow" markdown>

**Linear algebra perspective**: LP optimal solution = **basic feasible solution** = selecting $m$ columns to form an invertible $A_B$ → each step of the simplex method = solving the linear system $A_B^{-1}\mathbf{b}$ + column exchange (LU update)
**Link**: Ch22 LU decomposition is directly applied here

</div>

Linear Programming (LP) is the most fundamental type of optimization problem, with its theory and algorithms rooted in linear algebra.

!!! definition "Definition 25.1 (Standard form of linear programming)"
    The standard form of a **linear program** is
    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \mathbf{c}^T \mathbf{x}, \quad \text{s.t.} \quad A\mathbf{x} = \mathbf{b}, \quad \mathbf{x} \ge \mathbf{0},
    $$
    where $A \in \mathbb{R}^{m \times n}$ ($m < n$) is the constraint matrix, $\mathbf{b} \in \mathbb{R}^m$ is the right-hand side vector, and $\mathbf{c} \in \mathbb{R}^n$ is the objective vector. The inequality $\mathbf{x} \ge \mathbf{0}$ means componentwise nonnegativity.

!!! definition "Definition 25.2 (Basic feasible solution)"
    Let $A \in \mathbb{R}^{m \times n}$, $\operatorname{rank}(A) = m$. A **basis** $B$ of $A$ is an index set of $m$ columns such that $A_B$ (the submatrix formed by these $m$ columns) is invertible. The corresponding **basic feasible solution** (BFS) is
    $$
    \mathbf{x}_B = A_B^{-1} \mathbf{b}, \quad \mathbf{x}_N = \mathbf{0},
    $$
    where $N$ is the index set of nonbasic variables. If $\mathbf{x}_B \ge \mathbf{0}$, the basis is called a feasible basis.

!!! theorem "Theorem 25.1 (Fundamental theorem of linear programming)"
    If a linear program has an optimal solution, then there exists a basic feasible solution that is optimal.

??? proof "Proof"
    Let $\mathbf{x}^*$ be an optimal solution. If $\mathbf{x}^*$ is not a basic feasible solution, then the number $k$ of its positive components is less than $m$, or the columns of $A$ corresponding to the positive components are linearly dependent.

    **Case 1**: If the columns corresponding to positive components are linearly independent and $k < m$, we can augment with other columns to form a basis, so $\mathbf{x}^*$ is a degenerate basic feasible solution.

    **Case 2**: If the columns corresponding to positive components are linearly dependent, there exists a nonzero $\mathbf{d}$ such that $A\mathbf{d} = \mathbf{0}$, $d_j = 0$ ($j$ not in the positive component set). Consider $\mathbf{x}^* + t\mathbf{d}$ and $\mathbf{x}^* - t\mathbf{d}$, both satisfying $A\mathbf{x} = \mathbf{b}$. Adjust $t$ so that some positive component becomes zero while maintaining nonnegativity, and the objective value does not increase. Repeat this process until the columns corresponding to positive components are linearly independent, yielding a basic feasible solution. $\blacksquare$

!!! theorem "Theorem 25.2 (Linear algebra foundation of the simplex method)"
    Each step of the simplex method performs the following linear algebra operations:

    1. **Solve for the basic feasible solution**: $\mathbf{x}_B = A_B^{-1}\mathbf{b}$.
    2. **Compute reduced costs**: $\bar{c}_j = c_j - \mathbf{c}_B^T A_B^{-1} \mathbf{a}_j$, $j \in N$.
    3. **Determine the entering variable**: select a nonbasic variable $j$ with $\bar{c}_j < 0$.
    4. **Compute the direction**: $\mathbf{d}_B = A_B^{-1} \mathbf{a}_j$.
    5. **Determine the leaving variable**: $\theta^* = \min_{i: d_{B_i} > 0} \frac{(x_B)_i}{d_{B_i}}$ (minimum ratio rule).
    6. **Pivot**: update $B$ and $A_B^{-1}$.

    The core operations are solving the linear systems $A_B^{-1}\mathbf{b}$ and $A_B^{-1}\mathbf{a}_j$. In practice, LU decomposition is used with column exchange updates.

??? proof "Proof"
    The reduced cost $\bar{c}_j$ represents the rate of change of the objective value when the nonbasic variable $x_j$ increases by one unit. Let the current basis be $B$ with basic feasible solution $\mathbf{x}_B = A_B^{-1}\mathbf{b}$. If $x_j$ increases from $0$ to $\theta$, to maintain $A\mathbf{x} = \mathbf{b}$, we need $\mathbf{x}_B \to \mathbf{x}_B - \theta A_B^{-1}\mathbf{a}_j$. The change in objective value is
    $$
    \Delta z = c_j \theta - \mathbf{c}_B^T (A_B^{-1}\mathbf{a}_j) \theta = \bar{c}_j \theta.
    $$
    When $\bar{c}_j < 0$, increasing $x_j$ reduces the objective value. $\theta^*$ is determined by the tightest constraint from the nonnegativity requirement $\mathbf{x}_B - \theta A_B^{-1}\mathbf{a}_j \ge \mathbf{0}$. $\blacksquare$

!!! example "Example 25.1"
    **Solving a simple linear program.**

    $$
    \min -x_1 - 2x_2, \quad \text{s.t.} \quad x_1 + x_2 + s_1 = 4, \quad x_1 + 3x_2 + s_2 = 6, \quad \mathbf{x}, \mathbf{s} \ge \mathbf{0}.
    $$
    Initial basis $B = \{s_1, s_2\}$, $A_B = I_2$. Basic feasible solution $(s_1, s_2) = (4, 6)$, $z = 0$.

    Reduced costs $\bar{c}_{x_1} = -1$, $\bar{c}_{x_2} = -2$. Select $x_2$ to enter the basis, $\mathbf{d}_B = A_B^{-1}\mathbf{a}_{x_2} = (1, 3)^T$. Ratios: $4/1 = 4$, $6/3 = 2$, so $\theta^* = 2$, $s_2$ leaves the basis. New basic feasible solution $(x_2, s_1) = (2, 2)$, $z = -4$. Continue iterating until all reduced costs are nonnegative.

---

## 25.2 Least Squares Problems

<div class="context-flow" markdown>

**Three solution methods**: normal equations ($A^TA$, squared condition number) → QR decomposition (stable, Ch8) → SVD (most general, pseudoinverse + regularization) · **Tikhonov regularization** = spectral filtering: large $\sigma_i$ retained, small $\sigma_i$ suppressed

</div>

The least squares problem is a classical area at the intersection of linear algebra and optimization.

!!! definition "Definition 25.3 (Least squares problem)"
    Given $A \in \mathbb{R}^{m \times n}$ ($m \ge n$) and $\mathbf{b} \in \mathbb{R}^m$, the **least squares problem** is
    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|A\mathbf{x} - \mathbf{b}\|_2^2.
    $$

!!! definition "Definition 25.4 (Tikhonov regularization)"
    The **Tikhonov regularization** or ridge regression problem is
    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|A\mathbf{x} - \mathbf{b}\|_2^2 + \lambda \|\mathbf{x}\|_2^2, \quad \lambda > 0.
    $$
    The regularization parameter $\lambda$ controls the trade-off between solution smoothness and data fitting.

!!! theorem "Theorem 25.3 (Normal equations)"
    The solution of the least squares problem $\min \|A\mathbf{x} - \mathbf{b}\|_2^2$ satisfies the **normal equations**:
    $$
    A^T A \mathbf{x} = A^T \mathbf{b}.
    $$
    When $A$ has full column rank, the solution is unique: $\mathbf{x}^* = (A^T A)^{-1} A^T \mathbf{b}$. The solution of the Tikhonov regularization problem is
    $$
    \mathbf{x}_\lambda = (A^T A + \lambda I)^{-1} A^T \mathbf{b}.
    $$

??? proof "Proof"
    Let $f(\mathbf{x}) = \|A\mathbf{x} - \mathbf{b}\|_2^2 = \mathbf{x}^T A^T A \mathbf{x} - 2\mathbf{b}^T A\mathbf{x} + \|\mathbf{b}\|^2$. The necessary condition $\nabla f(\mathbf{x}) = 2A^TA\mathbf{x} - 2A^T\mathbf{b} = \mathbf{0}$ gives the normal equations. $A^TA$ is positive semidefinite; if $A$ has full column rank, then $A^TA$ is positive definite and the solution is unique.

    For Tikhonov regularization, $g(\mathbf{x}) = f(\mathbf{x}) + \lambda\|\mathbf{x}\|^2$, $\nabla g = 2(A^TA + \lambda I)\mathbf{x} - 2A^T\mathbf{b} = \mathbf{0}$. Since $\lambda > 0$, $A^TA + \lambda I$ is positive definite and the solution is unique. $\blacksquare$

!!! theorem "Theorem 25.4 (Solving least squares via QR decomposition)"
    Let $A = QR$, where $Q \in \mathbb{R}^{m \times n}$, $Q^TQ = I_n$, and $R \in \mathbb{R}^{n \times n}$ is upper triangular. Then the least squares solution is
    $$
    \mathbf{x}^* = R^{-1} Q^T \mathbf{b}.
    $$
    Numerically, QR decomposition is more stable than the normal equations because $\kappa(A^TA) = \kappa(A)^2$.

??? proof "Proof"
    $\|A\mathbf{x} - \mathbf{b}\|^2 = \|QR\mathbf{x} - \mathbf{b}\|^2$. Let $Q_\perp$ be the orthogonal complement of $Q$, so that $[Q \ Q_\perp]$ is an orthogonal matrix. Then
    $$
    \|QR\mathbf{x} - \mathbf{b}\|^2 = \|R\mathbf{x} - Q^T\mathbf{b}\|^2 + \|Q_\perp^T \mathbf{b}\|^2.
    $$
    The second term is independent of $\mathbf{x}$. Minimizing the first term: set $R\mathbf{x} = Q^T\mathbf{b}$, i.e., $\mathbf{x} = R^{-1}Q^T\mathbf{b}$. $\blacksquare$

!!! theorem "Theorem 25.5 (SVD solution of least squares and pseudoinverse)"
    Let $A = U\Sigma V^T$ be the SVD, $\Sigma = \operatorname{diag}(\sigma_1, \ldots, \sigma_r, 0, \ldots, 0)$. The minimum-norm solution of the least squares problem is
    $$
    \mathbf{x}^+ = A^+ \mathbf{b} = V \Sigma^+ U^T \mathbf{b} = \sum_{i=1}^{r} \frac{\mathbf{u}_i^T \mathbf{b}}{\sigma_i} \mathbf{v}_i,
    $$
    where $A^+ = V\Sigma^+ U^T$ is the Moore-Penrose pseudoinverse, $\Sigma^+ = \operatorname{diag}(\sigma_1^{-1}, \ldots, \sigma_r^{-1}, 0, \ldots, 0)$.

    The SVD expression of the Tikhonov regularization solution is
    $$
    \mathbf{x}_\lambda = \sum_{i=1}^{r} \frac{\sigma_i}{\sigma_i^2 + \lambda} (\mathbf{u}_i^T \mathbf{b}) \mathbf{v}_i.
    $$

??? proof "Proof"
    Substituting $A = U\Sigma V^T$:
    $$
    \|A\mathbf{x} - \mathbf{b}\|^2 = \|\Sigma V^T\mathbf{x} - U^T\mathbf{b}\|^2.
    $$
    Let $\mathbf{y} = V^T\mathbf{x}$, $\mathbf{c} = U^T\mathbf{b}$, then $\|\Sigma\mathbf{y} - \mathbf{c}\|^2 = \sum_{i=1}^r (\sigma_i y_i - c_i)^2 + \sum_{i=r+1}^m c_i^2$. Minimizing gives $y_i = c_i/\sigma_i$ ($i \le r$), $y_i$ arbitrary ($i > r$). The minimum-norm solution takes $y_i = 0$ ($i > r$), i.e., $\mathbf{x}^+ = V\Sigma^+ U^T\mathbf{b}$.

    For Tikhonov regularization, $(A^TA + \lambda I)\mathbf{x} = A^T\mathbf{b}$ becomes $(\Sigma^2 + \lambda I)\mathbf{y} = \Sigma \mathbf{c}$, so $y_i = \sigma_i c_i / (\sigma_i^2 + \lambda)$. $\blacksquare$

!!! example "Example 25.2"
    **Effect of condition number on least squares.**

    Let $A = \begin{pmatrix} 1 & 1 \\ \epsilon & 0 \\ 0 & \epsilon \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 2 \\ 0 \\ 0 \end{pmatrix}$. When $\epsilon$ is very small, $\kappa(A) \approx \sqrt{2}/\epsilon$. The normal equations give $A^TA = \begin{pmatrix} 1+\epsilon^2 & 1 \\ 1 & 1+\epsilon^2 \end{pmatrix}$, $\kappa(A^TA) \approx 2/\epsilon^2$. Using QR decomposition or SVD avoids the accuracy loss from squaring the condition number.

!!! example "Example 25.3"
    **Spectral filtering effect of Tikhonov regularization.**

    The SVD expression $\mathbf{x}_\lambda = \sum_i \frac{\sigma_i}{\sigma_i^2 + \lambda} (\mathbf{u}_i^T\mathbf{b})\mathbf{v}_i$ shows: for large singular values $\sigma_i \gg \sqrt{\lambda}$, the filter factor $\sigma_i/(\sigma_i^2 + \lambda) \approx 1/\sigma_i$ (same as the pseudoinverse); for small singular values $\sigma_i \ll \sqrt{\lambda}$, the filter factor $\approx \sigma_i/\lambda \approx 0$ (suppressed). This effectively truncates the amplification of noise in the directions of small singular values.

---

## 25.3 Semidefinite Programming (SDP)

<div class="context-flow" markdown>

**Generalization chain**: LP ($\mathbf{x} \ge 0$) → SDP ($X \succeq 0$) — generalizing nonnegativity constraints to **positive semidefinite cone** constraints · strong duality + complementary slackness $X^*S^* = 0$
**Power**: Goemans-Williamson relaxation for MAX-CUT (0.878 approximation ratio) · nuclear norm minimization = SDP (Section 25.4)

</div>

Semidefinite programming is a natural generalization of linear programming to matrix spaces.

!!! definition "Definition 25.5 (Semidefinite programming)"
    The standard form of **Semidefinite Programming** (SDP) is
    $$
    \min_{X \in \operatorname{Sym}(n)} \langle C, X \rangle, \quad \text{s.t.} \quad \langle A_i, X \rangle = b_i \; (i = 1, \ldots, m), \quad X \succeq 0,
    $$
    where $\langle A, B \rangle = \operatorname{tr}(A^T B)$ is the matrix inner product, $C, A_i \in \operatorname{Sym}(n)$, and $X \succeq 0$ means $X$ is positive semidefinite.

!!! definition "Definition 25.6 (SDP dual problem)"
    The **dual problem** of the standard form SDP is
    $$
    \max_{\mathbf{y} \in \mathbb{R}^m} \mathbf{b}^T \mathbf{y}, \quad \text{s.t.} \quad \sum_{i=1}^{m} y_i A_i + S = C, \quad S \succeq 0.
    $$
    The primal objective value (minimization) is always no less than the dual objective value (maximization): $\langle C, X \rangle \ge \mathbf{b}^T\mathbf{y}$ (weak duality).

!!! theorem "Theorem 25.6 (SDP strong duality)"
    If the primal problem (or the dual problem) satisfies the **Slater condition** (there exists a strictly feasible solution $X \succ 0$ such that $\langle A_i, X \rangle = b_i$), then strong duality holds:
    $$
    \min_X \langle C, X \rangle = \max_{\mathbf{y}} \mathbf{b}^T \mathbf{y},
    $$
    and the dual optimum is attained. At the optimal solution $(X^*, \mathbf{y}^*, S^*)$, the complementary slackness condition $X^* S^* = 0$ holds.

??? proof "Proof"
    **Proof sketch (conic programming duality theory).** SDP is a special case of conic programming, with the constraint cone being the positive semidefinite cone $\mathcal{S}_+^n$. The Slater condition ensures constraint regularity.

    Proof of weak duality: let $X$ be primal feasible and $(\mathbf{y}, S)$ be dual feasible, then
    $$
    \langle C, X \rangle - \mathbf{b}^T\mathbf{y} = \langle C, X \rangle - \sum_i y_i \langle A_i, X \rangle = \langle C - \sum_i y_i A_i, X \rangle = \langle S, X \rangle \ge 0,
    $$
    where the last inequality is guaranteed by $S \succeq 0$, $X \succeq 0$ ($\langle S, X \rangle = \operatorname{tr}(SX) \ge 0$, since the trace of the product of two positive semidefinite matrices is nonnegative).

    Strong duality requires the Slater condition and the separating hyperplane theorem from convex analysis to show that the duality gap is zero. Complementary slackness $\langle S^*, X^* \rangle = 0$ means $\operatorname{tr}(S^*X^*) = 0$; since $S^*, X^* \succeq 0$, this is equivalent to $S^*X^* = 0$. $\blacksquare$

!!! theorem "Theorem 25.7 (LP as a special case of SDP)"
    The linear program $\min \mathbf{c}^T\mathbf{x}$, s.t. $A\mathbf{x} = \mathbf{b}$, $\mathbf{x} \ge \mathbf{0}$ can be formulated as an SDP: take $X = \operatorname{diag}(\mathbf{x})$, $C = \operatorname{diag}(\mathbf{c})$, with the constraint $X \succeq 0$ (i.e., $X$ is a nonnegative diagonal matrix).

??? proof "Proof"
    Let $X = \operatorname{diag}(x_1, \ldots, x_n)$. Then $\langle C, X \rangle = \operatorname{tr}(\operatorname{diag}(\mathbf{c}) \operatorname{diag}(\mathbf{x})) = \mathbf{c}^T\mathbf{x}$. The constraints $\langle A_i, X \rangle = b_i$ can encode each row of $A\mathbf{x} = \mathbf{b}$. $X \succeq 0$ for a diagonal matrix is equivalent to $\mathbf{x} \ge \mathbf{0}$. $\blacksquare$

!!! example "Example 25.4"
    **SDP relaxation of the maximum cut problem.**

    The maximum cut (MAX-CUT) problem on a graph $G = (V, E)$ can be formulated as
    $$
    \max \frac{1}{4} \sum_{(i,j) \in E} w_{ij}(1 - x_i x_j), \quad x_i \in \{-1, +1\}.
    $$
    Let $X = \mathbf{x}\mathbf{x}^T$, then the constraints are equivalent to $X \succeq 0$, $\operatorname{rank}(X) = 1$, $X_{ii} = 1$. Dropping the rank-one constraint gives the SDP relaxation:
    $$
    \max \frac{1}{4} \langle L, X \rangle, \quad \text{s.t.} \quad X_{ii} = 1, \quad X \succeq 0,
    $$
    where $L$ is the Laplacian matrix of the graph. The Goemans-Williamson theorem guarantees that this relaxation achieves an approximation ratio of at least $0.878$.

!!! example "Example 25.5"
    **Matrix norm minimization.**

    Minimizing $\|A\|_{\text{op}}$ (operator norm) is equivalent to the SDP:
    $$
    \min t, \quad \text{s.t.} \quad \begin{pmatrix} tI & A \\ A^T & tI \end{pmatrix} \succeq 0.
    $$
    This is because the Schur complement condition $tI - A^T(tI)^{-1}A \succeq 0$ is equivalent to $t^2 I \succeq A^TA$, i.e., $t \ge \sigma_{\max}(A)$.

---

## 25.4 Matrix Completion

<div class="context-flow" markdown>

**Idea**: rank constraint (NP-hard) → **nuclear norm** relaxation (tightest convex envelope of rank) = SDP → incoherence condition + $O(\mu^2 r \log^2 n)$ observations → exact recovery
**Link**: Matrix version of Ch21 tensor decomposition · mathematical foundation of Netflix recommendation/collaborative filtering

</div>

Matrix completion is the problem of recovering a low-rank matrix from partial observations.

!!! definition "Definition 25.7 (Matrix completion problem)"
    Given partial entries $\{M_{ij} : (i,j) \in \Omega\}$ of $M \in \mathbb{R}^{m \times n}$, the **matrix completion problem** is
    $$
    \text{find } X \in \mathbb{R}^{m \times n}, \quad \operatorname{rank}(X) \le r, \quad X_{ij} = M_{ij} \; \forall (i,j) \in \Omega.
    $$
    Since the rank constraint is nonconvex, in practice nuclear norm relaxation is typically used.

!!! definition "Definition 25.8 (Nuclear norm)"
    The **nuclear norm** of a matrix $X \in \mathbb{R}^{m \times n}$ is defined as
    $$
    \|X\|_* = \sum_{i=1}^{\min(m,n)} \sigma_i(X) = \operatorname{tr}\!\left(\sqrt{X^T X}\right),
    $$
    where $\sigma_i(X)$ are the singular values. The nuclear norm is the tightest convex relaxation of the rank function (on the operator norm ball).

!!! theorem "Theorem 25.8 (SDP representation of nuclear norm minimization)"
    The nuclear norm minimization
    $$
    \min_{X} \|X\|_*, \quad \text{s.t.} \quad X_{ij} = M_{ij} \; \forall (i,j) \in \Omega
    $$
    is equivalent to the SDP:
    $$
    \min \frac{1}{2}(\operatorname{tr}(W_1) + \operatorname{tr}(W_2)), \quad \text{s.t.} \quad \begin{pmatrix} W_1 & X \\ X^T & W_2 \end{pmatrix} \succeq 0, \quad X_{ij} = M_{ij}.
    $$

??? proof "Proof"
    For any $X$ and positive semidefinite block matrix $\begin{pmatrix} W_1 & X \\ X^T & W_2 \end{pmatrix} \succeq 0$, by the Schur complement condition, $W_1 \succeq XW_2^{-1}X^T$ (when $W_2 \succ 0$). Therefore
    $$
    \operatorname{tr}(W_1) + \operatorname{tr}(W_2) \ge \operatorname{tr}(XW_2^{-1}X^T) + \operatorname{tr}(W_2).
    $$
    Optimizing over $W_2$, taking $W_2 = (X^TX)^{1/2}$, gives
    $$
    \operatorname{tr}(X(X^TX)^{-1/2}X^T) + \operatorname{tr}((X^TX)^{1/2}) = 2\operatorname{tr}((X^TX)^{1/2}) = 2\|X\|_*.
    $$
    Therefore $\frac{1}{2}(\operatorname{tr}(W_1) + \operatorname{tr}(W_2)) \ge \|X\|_*$, and equality is attainable. $\blacksquare$

!!! theorem "Theorem 25.9 (Exact recovery conditions for matrix completion)"
    Let $M \in \mathbb{R}^{m \times n}$ be a rank-$r$ matrix satisfying the **incoherence condition**: let $M = U\Sigma V^T$,
    $$
    \max_i \|U^T \mathbf{e}_i\|^2 \le \frac{\mu r}{m}, \quad \max_j \|V^T \mathbf{e}_j\|^2 \le \frac{\mu r}{n},
    $$
    where $\mu \ge 1$ is the incoherence parameter. If each entry in the observation set $\Omega$ is observed independently with probability $p$, and
    $$
    p \ge C \frac{\mu^2 r \log^2(m+n)}{m \wedge n}
    $$
    ($C$ is an absolute constant), then nuclear norm minimization recovers $M$ exactly with high probability.

??? proof "Proof"
    **Proof sketch (Candes-Recht 2009, Gross 2011).** One needs to show that $M$ is the unique solution of nuclear norm minimization. A dual certificate $Y$ is constructed satisfying: $\mathcal{P}_\Omega(Y) = \mathcal{P}_\Omega(\text{sgn}(M))$, $\mathcal{P}_T(Y) = UV^T$ (subgradient condition), $\|\mathcal{P}_{T^\perp}(Y)\|_{\text{op}} < 1$, where $T$ is the tangent space of $M$ and $\mathcal{P}_\Omega$ is the projection onto the observation set.

    The dual certificate is constructed via the **golfing scheme**: randomly partition $\Omega$ into multiple subsets $\Omega_1, \ldots, \Omega_L$, and iteratively correct $Y$ to approach the desired properties. The error at each step decays at a geometric rate, and the incoherence condition ensures good behavior (RIP property) of the projection operator $\mathcal{P}_{\Omega_j}\mathcal{P}_T$. $\blacksquare$

!!! example "Example 25.6"
    **The Netflix recommendation problem.**

    A user-movie rating matrix $M \in \mathbb{R}^{m \times n}$ ($m$ users, $n$ movies) is typically approximately low-rank (user preferences are determined by a small number of latent factors). The observed ratings form $\Omega$. Nuclear norm minimization can recover the complete rating matrix from sparse observations, enabling recommendations.

!!! example "Example 25.7"
    **Alternating minimization for low-rank recovery.**

    In practice, alternating minimization is often used instead of nuclear norm minimization: let $X = LR^T$ ($L \in \mathbb{R}^{m \times r}$, $R \in \mathbb{R}^{n \times r}$), and alternately optimize
    $$
    L \leftarrow \arg\min_L \sum_{(i,j) \in \Omega} (L_i^T R_j - M_{ij})^2, \quad
    R \leftarrow \arg\min_R \sum_{(i,j) \in \Omega} (L_i^T R_j - M_{ij})^2.
    $$
    Each step is a least squares problem and is computationally efficient.

---

## 25.5 Compressed Sensing

<div class="context-flow" markdown>

**Core condition**: **RIP** — the measurement matrix $A$ approximately preserves distances on sparse vectors → $\delta_{2s} < \sqrt{2}-1$ guarantees $\ell_1$ minimization exactly recovers $s$-sparse signals
**Random matrix connection**: Gaussian random $A$ with $m = O(s\log(n/s))$ rows satisfies RIP (Ch23 concentration inequalities) → far fewer than $n$ measurements suffice for reconstruction

</div>

Compressed sensing exploits signal sparsity to recover signals from far fewer measurements than required by the Nyquist sampling theorem.

!!! definition "Definition 25.9 (Restricted Isometry Property, RIP)"
    A matrix $A \in \mathbb{R}^{m \times n}$ satisfies the $s$-th order **Restricted Isometry Property** (RIP) with constant $\delta_s \in (0, 1)$ if for all $s$-sparse vectors $\mathbf{x}$ (with at most $s$ nonzero entries),
    $$
    (1 - \delta_s) \|\mathbf{x}\|_2^2 \le \|A\mathbf{x}\|_2^2 \le (1 + \delta_s) \|\mathbf{x}\|_2^2.
    $$
    Intuitively, RIP requires $A$ to approximately preserve distances on sparse vectors.

!!! definition "Definition 25.10 (Basis pursuit)"
    **Basis Pursuit** (BP) recovers sparse signals via $\ell_1$ norm minimization:
    $$
    \min_{\mathbf{x} \in \mathbb{R}^n} \|\mathbf{x}\|_1, \quad \text{s.t.} \quad A\mathbf{x} = \mathbf{b}.
    $$
    The noisy version (BPDN): $\min \|\mathbf{x}\|_1$, s.t. $\|A\mathbf{x} - \mathbf{b}\|_2 \le \epsilon$.

!!! theorem "Theorem 25.10 (RIP guarantees exact recovery)"
    Let $A \in \mathbb{R}^{m \times n}$ satisfy $\delta_{2s} < \sqrt{2} - 1$. If $\mathbf{x}^0$ is an $s$-sparse vector and $\mathbf{b} = A\mathbf{x}^0$, then the basis pursuit solution $\hat{\mathbf{x}}$ is unique and $\hat{\mathbf{x}} = \mathbf{x}^0$.

??? proof "Proof"
    Let $\hat{\mathbf{x}}$ be the basis pursuit solution, $\mathbf{h} = \hat{\mathbf{x}} - \mathbf{x}^0$. Then $A\mathbf{h} = 0$ and $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}^0\|_1$.

    Let $S$ be the support of $\mathbf{x}^0$ ($|S| \le s$) and $S^c$ its complement. By the triangle inequality:
    $$
    \|\mathbf{x}^0 + \mathbf{h}\|_1 = \|\mathbf{x}^0_S + \mathbf{h}_S\|_1 + \|\mathbf{h}_{S^c}\|_1 \ge \|\mathbf{x}^0_S\|_1 - \|\mathbf{h}_S\|_1 + \|\mathbf{h}_{S^c}\|_1.
    $$
    From $\|\hat{\mathbf{x}}\|_1 \le \|\mathbf{x}^0\|_1 = \|\mathbf{x}^0_S\|_1$, we get $\|\mathbf{h}_{S^c}\|_1 \le \|\mathbf{h}_S\|_1$ (**cone constraint**).

    Now sort the entries on $S^c$ by decreasing absolute value and partition into blocks of size $s$: $S_1, S_2, \ldots$. Using the cone constraint and RIP, one can show that
    $$
    \|\mathbf{h}\|_2^2 \le \frac{2\delta_{2s}}{1 - \delta_{2s}} \|\mathbf{h}_S\|_2 \|\mathbf{h}\|_2 + \text{tail terms},
    $$
    and when $\delta_{2s} < \sqrt{2} - 1$, one can derive $\|\mathbf{h}\|_2 = 0$, i.e., exact recovery. $\blacksquare$

!!! theorem "Theorem 25.11 (Gaussian random matrices satisfy RIP)"
    Let $A \in \mathbb{R}^{m \times n}$ with entries $A_{ij} \sim N(0, 1/m)$ i.i.d. If
    $$
    m \ge C \delta^{-2} s \log(n/s),
    $$
    then $A$ satisfies the $s$-th order RIP with constant $\delta_s \le \delta$ with high probability.

??? proof "Proof"
    **Proof sketch.** Fix an $s$-sparse vector $\mathbf{x}$. $\|A\mathbf{x}\|^2 = \sum_{i=1}^m (\mathbf{a}_i^T\mathbf{x})^2$, where $\mathbf{a}_i$ are the rows of $A$. $\mathbf{a}_i^T\mathbf{x} \sim N(0, \|\mathbf{x}\|^2/m)$, so $\|A\mathbf{x}\|^2/\|\mathbf{x}\|^2$ follows $\chi^2(m)/m$. By concentration inequalities:
    $$
    P\!\left(\left|\frac{\|A\mathbf{x}\|^2}{\|\mathbf{x}\|^2} - 1\right| > \delta\right) \le 2e^{-cm\delta^2}.
    $$
    Taking a union bound over all $s$-sparse vectors: the unit sphere of $s$-sparse vectors can be covered by an $\epsilon$-net of size $(Cn/s)^s$. Taking $m \ge C'\delta^{-2}s\log(n/s)$ ensures the union bound probability tends to zero. $\blacksquare$

!!! example "Example 25.8"
    **Sparse signal recovery.**

    Let $\mathbf{x}^0 \in \mathbb{R}^{1000}$ have $s = 20$ nonzero entries. Take $m = 200$ Gaussian random measurements $\mathbf{b} = A\mathbf{x}^0$. By Theorem 25.11, $m = 200 \ge C \cdot 20 \cdot \log(50) \approx 78C$ satisfies the RIP condition (for a moderate constant $C$). Basis pursuit can exactly recover the $1000$-dimensional signal from these $200$ measurements.

---

## 25.6 Principal Component Analysis (PCA)

<div class="context-flow" markdown>

**SVD is PCA**: principal component directions = eigenvectors of $S = \frac{1}{n}\bar{X}^T\bar{X}$ = right singular vectors of $\bar{X}$ → Eckart-Young optimal low-rank approximation
**Robust PCA**: $M = L + S$ (low-rank + sparse) → $\|L\|_* + \lambda\|S\|_1$ convex relaxation → video foreground/background separation · link to Ch23 BBP phase transition (when can the signal be detected?)

</div>

Principal component analysis is a classical method for data dimensionality reduction and feature extraction, with its mathematical foundation entirely built on linear algebra.

!!! definition "Definition 25.11 (Principal component analysis)"
    Given a data matrix $X \in \mathbb{R}^{n \times p}$ ($n$ samples, $p$ features), after centering $\bar{X} = X - \frac{1}{n}\mathbf{1}\mathbf{1}^TX$, **Principal Component Analysis** (PCA) seeks orthogonal directions $\mathbf{v}_1, \ldots, \mathbf{v}_k$ that maximize projected variance:
    $$
    \mathbf{v}_1 = \arg\max_{\|\mathbf{v}\|=1} \frac{1}{n}\|\bar{X}\mathbf{v}\|^2 = \arg\max_{\|\mathbf{v}\|=1} \mathbf{v}^T S \mathbf{v},
    $$
    where $S = \frac{1}{n}\bar{X}^T\bar{X}$ is the sample covariance matrix.

!!! definition "Definition 25.12 (Robust PCA)"
    **Robust PCA** decomposes a data matrix into a low-rank part and a sparse part:
    $$
    \min_{L, S} \|L\|_* + \lambda \|S\|_1, \quad \text{s.t.} \quad L + S = M,
    $$
    where $\|L\|_*$ is the nuclear norm (promoting low rank), $\|S\|_1 = \sum_{ij}|S_{ij}|$ (promoting sparsity), and $\lambda > 0$ is the trade-off parameter.

!!! theorem "Theorem 25.12 (SVD derivation of PCA)"
    Let $\bar{X} = U\Sigma V^T$ be the SVD. Then:

    1. The $k$-th principal component direction $\mathbf{v}_k$ is the $k$-th column of $V$ (the $k$-th eigenvector of $S$);
    2. The $k$-th principal component scores are $\bar{X}\mathbf{v}_k = \sigma_k \mathbf{u}_k$;
    3. The proportion of variance explained by the first $k$ principal components is $\sum_{i=1}^k \sigma_i^2 / \sum_{i=1}^p \sigma_i^2$;
    4. The best rank-$k$ approximation to $\bar{X}$ (in Frobenius norm) is $\bar{X}_k = U_k \Sigma_k V_k^T = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$.

??? proof "Proof"
    $S = \frac{1}{n}\bar{X}^T\bar{X} = \frac{1}{n}V\Sigma^2 V^T$, so the eigenvalues of $S$ are $\sigma_i^2/n$ and the eigenvectors are the columns of $V$.

    First principal component: $\max_{\|\mathbf{v}\|=1} \mathbf{v}^T S \mathbf{v} = \sigma_1^2/n$, attained at $\mathbf{v}_1$.

    Projected variance: $\frac{1}{n}\|\bar{X}\mathbf{v}_k\|^2 = \frac{1}{n}\sigma_k^2$.

    Best rank-$k$ approximation: by the Eckart-Young-Mirsky theorem, $\min_{\operatorname{rank}(B)\le k}\|\bar{X} - B\|_F = \sqrt{\sum_{i=k+1}^p \sigma_i^2}$, with optimal $B = \bar{X}_k$. $\blacksquare$

!!! theorem "Theorem 25.13 (Exact recovery for robust PCA)"
    Let $M = L_0 + S_0$, where $L_0$ has rank $r$ (satisfying the incoherence condition) and $S_0$ is sparse (with nonzero entry fraction $\rho$). Taking $\lambda = 1/\sqrt{\max(m,n)}$, under conditions that $r$ is sufficiently small and $\rho$ is sufficiently small, the robust PCA solution $(\hat{L}, \hat{S})$ satisfies $\hat{L} = L_0$, $\hat{S} = S_0$ with high probability.

??? proof "Proof"
    **Proof sketch (Candes-Li-Ma-Wright 2011).** Construct a dual certificate satisfying subgradient conditions. Let $L_0 = U\Sigma V^T$; one needs to construct a matrix $W$ satisfying:
    $$
    W = UV^T + D, \quad \|W\|_{\text{op}} \le 1, \quad \mathcal{P}_\Omega(W) = \lambda \operatorname{sgn}(S_0), \quad \|\mathcal{P}_{\Omega^c}(W - UV^T)\|_\infty < \lambda,
    $$
    where $\Omega$ is the support of $S_0$. The construction uses an alternating projection method (similar to the golfing scheme in matrix completion) together with analysis of the incoherence condition. $\blacksquare$

!!! example "Example 25.9"
    **Foreground-background separation in video surveillance.**

    Each frame of a video sequence is flattened into a column vector, forming the matrix $M$. The static background corresponds to a low-rank matrix $L$ (background frames are highly correlated), while the moving foreground corresponds to a sparse matrix $S$ (only a small number of pixels change). Robust PCA $M = L + S$ automatically performs foreground-background separation.

---

## 25.7 Low-Rank Approximation and Dimensionality Reduction

<div class="context-flow" markdown>

**Two major theorems**: **Eckart-Young-Mirsky** (truncated SVD = optimal low-rank approximation) + **JL lemma** (random projection $\mathbb{R}^n \to \mathbb{R}^{O(\log N/\epsilon^2)}$ approximately preserves distances)
**Link**: The concentration properties of Gaussian random matrices from Ch23 are at the core of the JL lemma proof

</div>

Low-rank approximation and dimensionality reduction are core techniques for handling high-dimensional data.

!!! definition "Definition 25.13 (Johnson-Lindenstrauss transform)"
    The **Johnson-Lindenstrauss (JL) transform** is a random linear map $\Phi: \mathbb{R}^n \to \mathbb{R}^m$ ($m \ll n$) that approximately preserves pairwise distances of a set of high-dimensional points in low-dimensional space.

!!! theorem "Theorem 25.14 (Johnson-Lindenstrauss lemma)"
    For any $\epsilon \in (0, 1)$ and $N$ points $\mathbf{x}_1, \ldots, \mathbf{x}_N$ in $\mathbb{R}^n$, there exists a map $f: \mathbb{R}^n \to \mathbb{R}^m$, $m = O(\epsilon^{-2} \log N)$, such that
    $$
    (1 - \epsilon)\|\mathbf{x}_i - \mathbf{x}_j\|_2^2 \le \|f(\mathbf{x}_i) - f(\mathbf{x}_j)\|_2^2 \le (1 + \epsilon)\|\mathbf{x}_i - \mathbf{x}_j\|_2^2
    $$
    for all $i, j$. In particular, $f$ can be taken as a random matrix $\Phi \in \mathbb{R}^{m \times n}$ ($\Phi_{ij} \sim N(0, 1/m)$), with $f(\mathbf{x}) = \Phi\mathbf{x}$.

??? proof "Proof"
    **Step 1: Concentration for a single vector.** Let $\mathbf{x} \in \mathbb{R}^n$, $\|\mathbf{x}\| = 1$, $\Phi_{ij} \sim N(0, 1/m)$ independent. Then $\|\Phi\mathbf{x}\|^2 = \sum_{i=1}^m (\boldsymbol{\phi}_i^T\mathbf{x})^2$, where $\boldsymbol{\phi}_i^T\mathbf{x} \sim N(0, 1/m)$. Therefore $m\|\Phi\mathbf{x}\|^2 \sim \chi^2(m)$. By concentration inequalities for the $\chi^2$ distribution:
    $$
    P\!\left(\left|\|\Phi\mathbf{x}\|^2 - 1\right| > \epsilon\right) \le 2\exp\!\left(-\frac{m\epsilon^2}{8}\right).
    $$

    **Step 2: Union bound.** For the $\binom{N}{2}$ difference vectors $\mathbf{x}_i - \mathbf{x}_j$ of $N$ points, taking a union bound:
    $$
    P(\exists \, i,j : \text{distance distortion} > \epsilon) \le 2\binom{N}{2}\exp\!\left(-\frac{m\epsilon^2}{8}\right) \le N^2 \exp\!\left(-\frac{m\epsilon^2}{8}\right).
    $$
    Taking $m \ge \frac{16}{\epsilon^2}\ln N$, the failure probability $\le N^2 \cdot N^{-2} = 1$ (in practice, a larger constant is used to drive the probability to zero). $\blacksquare$

!!! theorem "Theorem 25.15 (Eckart-Young-Mirsky theorem)"
    Let $A \in \mathbb{R}^{m \times n}$, $\operatorname{rank}(A) = r$, with SVD $A = \sum_{i=1}^r \sigma_i \mathbf{u}_i \mathbf{v}_i^T$. Then for $k < r$,
    $$
    \min_{\operatorname{rank}(B) \le k} \|A - B\|_F = \sqrt{\sum_{i=k+1}^r \sigma_i^2}, \quad
    \min_{\operatorname{rank}(B) \le k} \|A - B\|_{\text{op}} = \sigma_{k+1},
    $$
    and the optimal solution is the truncated SVD $A_k = \sum_{i=1}^k \sigma_i \mathbf{u}_i \mathbf{v}_i^T$.

??? proof "Proof"
    **Frobenius norm case.** Let $B$ be a matrix of rank $\le k$. $\|A - B\|_F^2 = \operatorname{tr}((A-B)^T(A-B))$. Let the column space of $B$ be $\mathcal{V}$ ($\dim \mathcal{V} \le k$) and the row space be $\mathcal{W}$ ($\dim \mathcal{W} \le k$).

    By the Courant-Fischer minimax principle, $\mathcal{V}^\perp$ is an $(n - k)$-dimensional subspace. The restriction of $A$ to $\mathcal{V}^\perp$ satisfies
    $$
    \|A - B\|_F^2 \ge \sum_{i=k+1}^r \sigma_i^2(A|_{\mathcal{V}^\perp}) \ge \sum_{i=k+1}^r \sigma_i^2(A),
    $$
    where the latter inequality follows from the minimax characterization of singular values. Equality is attained when $B = A_k$.

    **Operator norm case.** $\ker(B)$ has dimension at least $n - k$. Choose $\mathbf{v} \in \ker(B) \cap \operatorname{span}(\mathbf{v}_1, \ldots, \mathbf{v}_{k+1})$ (a dimension argument guarantees this intersection is nontrivial), $\|\mathbf{v}\| = 1$, then
    $$
    \|A - B\|_{\text{op}} \ge \|(A - B)\mathbf{v}\| = \|A\mathbf{v}\| \ge \sigma_{k+1}.
    $$
    Equality is attained when $B = A_k$. $\blacksquare$

!!! example "Example 25.10"
    **Application of the JL lemma: approximate nearest neighbor search.**

    Given $N = 10^6$ data points in $\mathbb{R}^{10000}$, the JL lemma guarantees that they can be projected to $m = O(\epsilon^{-2} \cdot 6\log 10) \approx O(\epsilon^{-2} \cdot 14)$ dimensions while preserving distances within a factor of $1 \pm \epsilon$. Taking $\epsilon = 0.1$, $m \approx 1400$, which significantly reduces the computational complexity of nearest neighbor search.

!!! example "Example 25.11"
    **Low-rank approximation in image compression.**

    A grayscale image can be represented as a matrix $A \in \mathbb{R}^{m \times n}$. If $A$ has rapidly decaying singular values, the rank-$k$ truncated SVD $A_k$ approximates $mn$ pixel values using $k(m + n + 1)$ parameters. For example, when $m = n = 1000$ and $k = 50$, the compression ratio is approximately $1000000 / 100050 \approx 10:1$.

---

## 25.8 Eigenvalue Optimization

<div class="context-flow" markdown>

**Eigenvalues = extrema**: Rayleigh quotient $R_A(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}/\|\mathbf{x}\|^2$ → $\lambda_{\min} \le R_A \le \lambda_{\max}$ · **Courant-Fischer** minimax → **Weyl perturbation inequality** $|\gamma_i - \alpha_i| \le \|B\|$
**Convergence**: Rayleigh quotient iteration (cubic convergence, Ch22) · graph Laplacian $\lambda_2$ = connectivity (Fiedler) · link to Ch24 eigenvalue problems on Stiefel manifolds

</div>

Many optimization problems have solutions that can be characterized through eigenvalues.

!!! definition "Definition 25.14 (Rayleigh quotient)"
    Let $A \in \operatorname{Sym}(n)$. The **Rayleigh quotient** is defined as
    $$
    R_A(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T \mathbf{x}}, \quad \mathbf{x} \ne \mathbf{0}.
    $$
    The Rayleigh quotient is a homogeneous function of degree zero, so it can be restricted to the unit sphere $\|\mathbf{x}\| = 1$.

!!! definition "Definition 25.15 (Generalized Rayleigh quotient)"
    Let $A, B \in \operatorname{Sym}(n)$, $B \succ 0$. The **generalized Rayleigh quotient** is
    $$
    R_{A,B}(\mathbf{x}) = \frac{\mathbf{x}^T A \mathbf{x}}{\mathbf{x}^T B \mathbf{x}}, \quad \mathbf{x} \ne \mathbf{0}.
    $$
    Its extrema are related to the generalized eigenvalue problem $A\mathbf{x} = \lambda B\mathbf{x}$.

!!! theorem "Theorem 25.16 (Extrema of the Rayleigh quotient)"
    Let $A \in \operatorname{Sym}(n)$ with eigenvalues $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$ and corresponding eigenvectors $\mathbf{v}_1, \ldots, \mathbf{v}_n$. Then
    $$
    \lambda_1 = \min_{\mathbf{x} \ne \mathbf{0}} R_A(\mathbf{x}), \quad \lambda_n = \max_{\mathbf{x} \ne \mathbf{0}} R_A(\mathbf{x}),
    $$
    with the extrema attained at $\mathbf{x} = \mathbf{v}_1$ and $\mathbf{x} = \mathbf{v}_n$, respectively.

??? proof "Proof"
    Let $\mathbf{x} = \sum_{i=1}^n c_i \mathbf{v}_i$, $\|\mathbf{x}\| = 1$, i.e., $\sum c_i^2 = 1$. Then
    $$
    R_A(\mathbf{x}) = \mathbf{x}^T A \mathbf{x} = \sum_{i=1}^n \lambda_i c_i^2.
    $$
    Since $\sum c_i^2 = 1$ and $\lambda_1 \le \cdots \le \lambda_n$,
    $$
    \lambda_1 = \lambda_1 \sum c_i^2 \le \sum \lambda_i c_i^2 \le \lambda_n \sum c_i^2 = \lambda_n.
    $$
    The lower bound is attained when $c_1 = 1$ (i.e., $\mathbf{x} = \mathbf{v}_1$), and the upper bound is attained when $c_n = 1$. $\blacksquare$

<div class="context-flow" markdown>

**Insight**: Courant-Fischer transforms eigenvalues from "algebraic objects" ($\det(A-\lambda I)=0$) into "optimization objects" (extrema over subspaces) — from this, one derives Weyl perturbation, Lidskii inequalities, graph partitioning, and all other eigenvalue inequalities

</div>

!!! theorem "Theorem 25.17 (Courant-Fischer minimax theorem)"
    Let $A \in \operatorname{Sym}(n)$ with eigenvalues $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$. Then
    $$
    \lambda_k = \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} \mathbf{x}^T A \mathbf{x} = \max_{\dim W = n-k+1} \min_{\mathbf{x} \in W, \|\mathbf{x}\|=1} \mathbf{x}^T A \mathbf{x}.
    $$

??? proof "Proof"
    **Proof of the first equality.** Let $V_k = \operatorname{span}(\mathbf{v}_1, \ldots, \mathbf{v}_k)$.

    **Upper bound**: Take $V = V_k$, then $\max_{\mathbf{x} \in V_k, \|\mathbf{x}\|=1} \mathbf{x}^TA\mathbf{x} = \lambda_k$ (since the largest eigenvalue of $A$ restricted to $V_k$ is $\lambda_k$).

    **Lower bound**: Let $V$ be any $k$-dimensional subspace. Consider $W = \operatorname{span}(\mathbf{v}_k, \ldots, \mathbf{v}_n)$, $\dim W = n - k + 1$. By the dimension formula, $\dim(V \cap W) \ge k + (n-k+1) - n = 1$. Choose $\mathbf{x} \in V \cap W$, $\|\mathbf{x}\| = 1$, then
    $$
    \max_{\mathbf{y} \in V, \|\mathbf{y}\|=1} \mathbf{y}^T A \mathbf{y} \ge \mathbf{x}^T A \mathbf{x} \ge \lambda_k,
    $$
    where the latter inequality holds because $\mathbf{x} \in W$ implies $\mathbf{x}^T A \mathbf{x} \ge \lambda_k$ (the smallest eigenvalue of $A$ on $W$ is $\lambda_k$).

    Combining gives $\lambda_k = \min_{\dim V = k} \max_{\mathbf{x} \in V} \mathbf{x}^TA\mathbf{x}$. The proof of the second equality is similar. $\blacksquare$

!!! theorem "Theorem 25.18 (Weyl eigenvalue perturbation inequality)"
    Let $A, B \in \operatorname{Sym}(n)$ with eigenvalues $\alpha_1 \le \cdots \le \alpha_n$ and $\beta_1 \le \cdots \le \beta_n$ respectively, and let $A + B$ have eigenvalues $\gamma_1 \le \cdots \le \gamma_n$. Then
    $$
    \alpha_i + \beta_1 \le \gamma_i \le \alpha_i + \beta_n, \quad i = 1, \ldots, n.
    $$
    In particular, $|\gamma_i - \alpha_i| \le \|B\|_{\text{op}}$.

??? proof "Proof"
    By the Courant-Fischer theorem,
    $$
    \gamma_k = \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} \mathbf{x}^T(A+B)\mathbf{x} \le \min_{\dim V = k} \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1} (\mathbf{x}^TA\mathbf{x} + \beta_n) = \alpha_k + \beta_n.
    $$
    Similarly, $\gamma_k \ge \alpha_k + \beta_1$. Taking $B' = -B$ gives symmetrically $|\gamma_k - \alpha_k| \le \max(|\beta_1|, |\beta_n|) = \|B\|_{\text{op}}$. $\blacksquare$

!!! theorem "Theorem 25.19 (Lidskii inequality)"
    Let $A, B \in \operatorname{Sym}(n)$ with eigenvalues $\alpha_1 \le \cdots \le \alpha_n$ and $\beta_1 \le \cdots \le \beta_n$ respectively. Let $C = A + B$ with eigenvalues $\gamma_1 \le \cdots \le \gamma_n$. Then for any $1 \le i_1 < i_2 < \cdots < i_k \le n$,
    $$
    \sum_{j=1}^{k} \gamma_{i_j} \le \sum_{j=1}^{k} \alpha_{i_j} + \sum_{j=1}^{k} \beta_{n-k+j}.
    $$

??? proof "Proof"
    Using a generalization of the Courant-Fischer theorem. For the index set $I = \{i_1, \ldots, i_k\}$, construct an appropriate chain of subspaces $V_1 \subset V_2 \subset \cdots \subset V_k$ with $\dim V_j = i_j$. Use the Courant-Fischer characterization of each $\gamma_{i_j}$, then through dimension arguments on subspace intersections, decompose the Rayleigh quotient contributions of $\gamma_{i_j}$ into those from $A$ and $B$ to obtain the inequality. For the complete proof, see Bhatia, *Matrix Analysis*, Theorem III.4.1. $\blacksquare$

!!! example "Example 25.12"
    **Rayleigh quotient iteration.**

    To compute eigenvalues and eigenvectors of $A \in \operatorname{Sym}(n)$: given initial $\mathbf{x}_0$, $\|\mathbf{x}_0\| = 1$, iterate
    $$
    \rho_k = \mathbf{x}_k^T A \mathbf{x}_k, \quad (A - \rho_k I)\mathbf{y}_{k+1} = \mathbf{x}_k, \quad \mathbf{x}_{k+1} = \mathbf{y}_{k+1}/\|\mathbf{y}_{k+1}\|.
    $$
    Rayleigh quotient iteration has cubic convergence: $|\rho_{k+1} - \lambda| = O(|\rho_k - \lambda|^3)$, making it one of the fastest iterative methods for computing eigenvalues.

!!! example "Example 25.13"
    **Application of the Courant-Fischer theorem in graph theory.**

    Let $L$ be the Laplacian matrix of a graph $G$. By the Courant-Fischer theorem,
    $$
    \lambda_2(L) = \min_{\mathbf{x} \perp \mathbf{1}, \|\mathbf{x}\|=1} \mathbf{x}^T L \mathbf{x} = \min_{\mathbf{x} \perp \mathbf{1}} \frac{\sum_{(i,j)\in E}(x_i - x_j)^2}{\sum_i x_i^2}.
    $$
    $\lambda_2(L)$ (the Fiedler value) measures graph connectivity: $\lambda_2 > 0$ if and only if $G$ is connected, and larger $\lambda_2$ indicates stronger connectivity. The Fiedler vector (the eigenvector corresponding to $\lambda_2$) is used for graph partitioning.

!!! example "Example 25.14"
    **Application of Weyl's inequality: stability of eigenvalues.**

    Let $A$ be a symmetric matrix and $E$ a symmetric perturbation with $\|E\|_{\text{op}} = \epsilon$. By Weyl's inequality, the difference between the $k$-th eigenvalue of $A + E$ and the $k$-th eigenvalue of $A$ does not exceed $\epsilon$. Therefore, to achieve eigenvalue accuracy of $10^{-6}$, it suffices to control the matrix entry precision to $O(10^{-6})$. This provides a theoretical guarantee for numerical eigenvalue computation.

---

## Chapter Summary

This chapter demonstrated the central role of linear algebra in optimization theory and algorithms:

1. The simplex method for **linear programming** is essentially an iterative process of solving linear systems, with basic feasible solutions corresponding to column selection.
2. **Least squares** problems are solved via normal equations, QR decomposition, or SVD; Tikhonov regularization suppresses noise amplification through spectral filtering.
3. **Semidefinite programming** generalizes linear programming to matrix spaces, with strong duality as its theoretical cornerstone.
4. **Matrix completion** uses nuclear norm relaxation to convert the NP-hard rank constraint problem into convex optimization.
5. In **compressed sensing**, the RIP condition (restricted isometry property from linear algebra) guarantees exact recovery via $\ell_1$ minimization.
6. **PCA** achieves optimal dimensionality reduction through SVD; robust PCA separates low-rank and sparse components.
7. The **JL lemma** and the **Eckart-Young-Mirsky theorem** provide theoretical foundations for dimensionality reduction.
8. The **Courant-Fischer theorem** and **Weyl's inequality** are fundamental tools for eigenvalue optimization and perturbation analysis.
