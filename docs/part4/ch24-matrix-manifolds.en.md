# Chapter 24  Matrix Manifolds

<div class="context-flow" markdown>

**Prerequisites**: Orthogonal matrices (Ch7) · SVD (Ch8) · Positive definite matrices (Ch7)

**Arc**: Matrix constraint sets → Smooth manifolds → Riemannian geometry (geodesics/gradients) → Manifold optimization

**Core objects**: $O(n)$ / $\text{St}(k,n)$ / $\text{Gr}(k,n)$ / $\mathcal{P}(n)$ — turning constrained optimization into unconstrained optimization on manifolds → Links to Ch25

**Further connections**：Matrix manifolds are widely applied in computer vision (essential matrix estimation, shape analysis), robotics (pose estimation on $SO(3)$), medical imaging (diffusion tensor MRI), and machine learning (low-rank optimization, subspace tracking)

</div>

Matrix manifolds are smooth manifolds formed by matrices satisfying certain constraints. These manifolds are ubiquitous in control theory, signal processing, computer vision, and machine learning. This chapter systematically introduces several important matrix manifolds — the general linear group, orthogonal group, Stiefel manifold, Grassmann manifold, and the manifold of positive definite matrices — along with their differential geometric structures and optimization methods.

---

## 24.1 Basic Concepts of Matrix Manifolds

<div class="context-flow" markdown>

**Basic toolbox**: Embedded submanifold ($h^{-1}(0)$, implicit function theorem) → Tangent space ($\ker Dh$) → Riemannian metric (inner product on tangent space) → Geodesic (locally shortest path) → Exponential map

</div>

We first review the basic language of manifold theory and explain how matrix sets naturally acquire manifold structure.

!!! definition "Definition 24.1 (Smooth manifold)"
    A **smooth manifold** $\mathcal{M}$ is a Hausdorff, second-countable topological space equipped with a maximal smooth atlas $\{(U_\alpha, \varphi_\alpha)\}$, where each $U_\alpha$ is an open subset of $\mathcal{M}$, $\varphi_\alpha: U_\alpha \to \mathbb{R}^d$ is a homeomorphism, and the transition maps $\varphi_\beta \circ \varphi_\alpha^{-1}$ are smooth. $d$ is called the dimension of $\mathcal{M}$.

!!! definition "Definition 24.2 (Embedded submanifold)"
    Let $\mathcal{M}$ be a subset of $\mathbb{R}^{n \times p}$. If $\mathcal{M}$ can be expressed as

    $$
    \mathcal{M} = \{ X \in \mathbb{R}^{n \times p} : h(X) = 0 \},
    $$

    where $h: \mathbb{R}^{n \times p} \to \mathbb{R}^m$ is a smooth map and $Dh(X)$ has full rank ($\operatorname{rank}(Dh(X)) = m$) for all $X \in \mathcal{M}$, then $\mathcal{M}$ is an **embedded submanifold** of $\mathbb{R}^{n \times p}$ of dimension $np - m$.

!!! definition "Definition 24.3 (Tangent space)"
    Let $\mathcal{M}$ be a smooth manifold and $X \in \mathcal{M}$. The **tangent space** $T_X \mathcal{M}$ of $\mathcal{M}$ at $X$ is the set of velocity vectors at $X$ of all smooth curves passing through $X$:

    $$
    T_X \mathcal{M} = \left\{ \gamma'(0) : \gamma: (-\epsilon, \epsilon) \to \mathcal{M} \text{ smooth}, \, \gamma(0) = X \right\}.
    $$

    If $\mathcal{M} = h^{-1}(0)$ is an embedded submanifold, then $T_X \mathcal{M} = \ker(Dh(X))$.

!!! definition "Definition 24.4 (Riemannian metric)"
    A **Riemannian metric** on manifold $\mathcal{M}$ is a smoothly varying family of inner products $\langle \cdot, \cdot \rangle_X: T_X\mathcal{M} \times T_X\mathcal{M} \to \mathbb{R}$, $X \in \mathcal{M}$. A manifold equipped with a Riemannian metric is called a Riemannian manifold.

!!! definition "Definition 24.5 (Geodesic)"
    A **geodesic** on a Riemannian manifold $(\mathcal{M}, \langle \cdot, \cdot \rangle)$ is a locally shortest curve. Formally, $\gamma: [0, 1] \to \mathcal{M}$ is a geodesic if and only if $\nabla_{\dot\gamma} \dot\gamma = 0$, where $\nabla$ is the Levi-Civita connection. The **exponential map** $\operatorname{Exp}_X: T_X\mathcal{M} \to \mathcal{M}$ maps a tangent vector $\xi$ to the endpoint at $t = 1$ of the geodesic starting from $X$ with initial velocity $\xi$.

!!! theorem "Theorem 24.1 (Implicit function theorem and submanifolds)"
    Let $h: \mathbb{R}^{n \times p} \to \mathbb{R}^m$ be a smooth map and $\mathbf{0}$ a regular value of $h$ (i.e., $Dh(X)$ has full rank for all $X \in h^{-1}(\mathbf{0})$). Then $\mathcal{M} = h^{-1}(\mathbf{0})$ is a smooth embedded submanifold of $\mathbb{R}^{n \times p}$ of dimension $np - m$.

??? proof "Proof"
    This is a direct consequence of the implicit function theorem. For $X_0 \in \mathcal{M}$, full rank of $Dh(X_0)$ implies there exists a neighborhood $U$ of $X_0$ and a smooth map $\psi$ such that $\mathcal{M} \cap U$ can be parameterized as the graph of a smooth function of $np - m$ free variables. These local parameterizations form a smooth atlas for $\mathcal{M}$. $\blacksquare$

!!! example "Example 24.1"
    **The unit sphere as a submanifold.**

    $S^{n-1} = \{\mathbf{x} \in \mathbb{R}^n : \|\mathbf{x}\|^2 = 1\} = h^{-1}(0)$, where $h(\mathbf{x}) = \mathbf{x}^T\mathbf{x} - 1$. $Dh(\mathbf{x}) = 2\mathbf{x}^T$, which clearly has full rank for $\mathbf{x} \in S^{n-1}$. The tangent space $T_{\mathbf{x}} S^{n-1} = \{\mathbf{v} : \mathbf{x}^T \mathbf{v} = 0\}$, i.e., vectors orthogonal to $\mathbf{x}$.

---

## 24.2 General Linear Group $GL(n)$

<div class="context-flow" markdown>

**Parent group**: $GL(n) = \{\det \ne 0\}$ is an open subset of $\mathbb{R}^{n^2}$ → Lie algebra $\mathfrak{gl}(n) = \mathbb{R}^{n \times n}$ → Matrix exponential $e^A$ connects Lie algebra and Lie group

**All matrix Lie groups are closed subgroups of $GL(n)$** (Cartan's theorem)

</div>

!!! definition "Definition 24.6 (General linear group)"
    The **general linear group** $GL(n, \mathbb{R})$ (abbreviated $GL(n)$) is the set of all $n \times n$ invertible real matrices:

    $$
    GL(n) = \{ A \in \mathbb{R}^{n \times n} : \det(A) \ne 0 \}.
    $$

    $GL(n)$ is a group under matrix multiplication and is an open subset of $\mathbb{R}^{n \times n}$, hence an $n^2$-dimensional smooth manifold. Similarly, $GL(n, \mathbb{C})$ is the group of all $n \times n$ invertible complex matrices.

!!! definition "Definition 24.7 (Lie algebra $\mathfrak{gl}(n)$)"
    The **Lie algebra** $\mathfrak{gl}(n)$ of $GL(n)$ is the tangent space of $GL(n)$ at the identity $I$, i.e., all $n \times n$ real matrices:

    $$
    \mathfrak{gl}(n) = T_I GL(n) = \mathbb{R}^{n \times n}.
    $$

    The Lie bracket on $\mathfrak{gl}(n)$ is defined as the matrix commutator $[A, B] = AB - BA$.

!!! theorem "Theorem 24.2 (Matrix exponential map)"
    The map $\exp: \mathfrak{gl}(n) \to GL(n)$, $A \mapsto e^A = \sum_{k=0}^{\infty} \frac{A^k}{k!}$ is a smooth map whose differential at $A = 0$ is the identity. Therefore $\exp$ is a diffeomorphism in a neighborhood of $0$.

??? proof "Proof"
    The matrix exponential series $e^A = \sum_{k=0}^\infty \frac{A^k}{k!}$ converges absolutely for all $A \in \mathbb{R}^{n \times n}$ (since $\|A^k/k!\| \le \|A\|^k/k!$ and the sum is bounded by $e^{\|A\|}$). Smoothness follows from term-by-term differentiation. At $A = 0$, $D\exp(0)(H) = \frac{d}{dt}\Big|_{t=0} e^{tH} = H$, so the differential is the identity. By the inverse function theorem, $\exp$ is a diffeomorphism near $0$. $\blacksquare$

!!! theorem "Theorem 24.3 (Left-invariant vector fields on $GL(n)$)"
    Every left-invariant vector field on $GL(n)$ is uniquely determined by its value at $I$. Specifically, for $A \in \mathfrak{gl}(n)$, the left-invariant vector field $\tilde{A}$ is defined by $\tilde{A}_g = gA$ ($g \in GL(n)$), and $[\tilde{A}, \tilde{B}] = \widetilde{[A, B]}$.

??? proof "Proof"
    The left translation $L_g: GL(n) \to GL(n)$, $h \mapsto gh$ has differential $(dL_g)_h(V) = gV$. The vector field $\tilde{A}$ satisfies $(dL_g)\tilde{A}_h = g(hA) = (gh)A = \tilde{A}_{gh}$, hence it is left-invariant. The Lie bracket $[\tilde{A}, \tilde{B}]_I = \frac{d}{dt}\Big|_{t=0} \frac{d}{ds}\Big|_{s=0} (e^{tA} e^{sB} e^{-tA}) = AB - BA = [A, B]$. $\blacksquare$

!!! example "Example 24.2"
    **Exponential map on $GL(2)$.**

    Let $A = \begin{pmatrix} 0 & -\theta \\ \theta & 0 \end{pmatrix}$. Then

    $$
    e^A = \begin{pmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{pmatrix} \in SO(2).
    $$

    This shows that skew-symmetric matrices generate rotation matrices via the exponential map.

!!! example "Example 24.3"
    **Exponential of a diagonal matrix.**

    Let $D = \operatorname{diag}(d_1, \ldots, d_n)$. Then $e^D = \operatorname{diag}(e^{d_1}, \ldots, e^{d_n})$. The eigenvalues of $e^D$ are all positive, so $e^D$ is a positive definite diagonal matrix.

---

## 24.3 Orthogonal Group and Unitary Group

<div class="context-flow" markdown>

**Constraint → Manifold**: $Q^TQ = I$ → $O(n)$, $\dim = \frac{n(n-1)}{2}$ · Lie algebra = **skew-symmetric matrices** $\text{Skew}(n)$ → Geodesic $\gamma(t) = Qe^{t\Omega}$ · The invariance of GOE from Ch23 is precisely $O(n)$-conjugation invariance

</div>

The orthogonal and unitary groups are the most important matrix Lie groups, with broad applications in geometry, physics, and engineering.

!!! definition "Definition 24.8 (Orthogonal group)"
    The **orthogonal group** $O(n)$ is the set of all $n \times n$ orthogonal matrices:

    $$
    O(n) = \{ Q \in \mathbb{R}^{n \times n} : Q^T Q = I_n \}.
    $$

    $O(n)$ is a compact Lie group of dimension $\frac{n(n-1)}{2}$. The **special orthogonal group** $SO(n) = \{Q \in O(n) : \det(Q) = 1\}$ is the identity component of $O(n)$.

!!! definition "Definition 24.9 (Unitary group)"
    The **unitary group** $U(n)$ is the set of all $n \times n$ unitary matrices:

    $$
    U(n) = \{ U \in \mathbb{C}^{n \times n} : U^* U = I_n \}.
    $$

    $U(n)$ is a compact Lie group of (real) dimension $n^2$. The **special unitary group** $SU(n) = \{U \in U(n) : \det(U) = 1\}$ has dimension $n^2 - 1$.

!!! theorem "Theorem 24.4 (Tangent space and Lie algebra of $O(n)$ and $SO(n)$)"
    The tangent space of $O(n)$ at $I$ (i.e., the Lie algebra) is the space of skew-symmetric matrices:

    $$
    \mathfrak{o}(n) = T_I O(n) = \operatorname{Skew}(n) = \{ \Omega \in \mathbb{R}^{n \times n} : \Omega^T = -\Omega \}.
    $$

    More generally, at $Q \in O(n)$, $T_Q O(n) = \{ Q\Omega : \Omega \in \operatorname{Skew}(n) \}$. The Lie algebra of $SO(n)$ is also $\operatorname{Skew}(n)$ (the determinant condition does not affect the tangent space).

??? proof "Proof"
    $O(n)$ is defined by the constraint $h(Q) = Q^T Q - I = 0$. $Dh(Q)(V) = V^T Q + Q^T V$. At $Q = I$, $Dh(I)(V) = V^T + V$. Therefore

    $$
    T_I O(n) = \ker(Dh(I)) = \{ V : V^T + V = 0 \} = \operatorname{Skew}(n).
    $$

    The dimension of $\operatorname{Skew}(n)$ is $\frac{n(n-1)}{2}$, matching the dimension of $O(n)$.

    At a general point $Q$, let $\gamma(t) = Qe^{t\Omega}$ ($\Omega \in \operatorname{Skew}(n)$). Then $\gamma(0) = Q$, $\gamma(t)^T\gamma(t) = e^{t\Omega^T}Q^TQe^{t\Omega} = e^{-t\Omega}e^{t\Omega} = I$, so $\gamma(t) \in O(n)$ and $\gamma'(0) = Q\Omega$. This shows $T_Q O(n) \supseteq \{Q\Omega : \Omega \in \operatorname{Skew}(n)\}$. Dimension counting shows equality. $\blacksquare$

!!! theorem "Theorem 24.5 (Lie algebra of $U(n)$)"
    The Lie algebra of $U(n)$ is the space of skew-Hermitian matrices:

    $$
    \mathfrak{u}(n) = T_I U(n) = \{ \Omega \in \mathbb{C}^{n \times n} : \Omega^* = -\Omega \}.
    $$

    The Lie algebra of $SU(n)$ is $\mathfrak{su}(n) = \{ \Omega \in \mathfrak{u}(n) : \operatorname{tr}(\Omega) = 0 \}$.

??? proof "Proof"
    Completely analogous to the proof for $O(n)$. The constraint is $U^* U = I$; differentiating gives $V^* U + U^* V = 0$, which at $U = I$ yields $V^* + V = 0$, i.e., $V$ is skew-Hermitian. For $SU(n)$, the additional condition $\det(U) = 1$ differentiates to $\operatorname{tr}(U^{-1}V) = \operatorname{tr}(V) = 0$ (at $U = I$). $\blacksquare$

!!! theorem "Theorem 24.6 (Geodesics on $O(n)$)"
    Equip $O(n)$ with the bi-invariant metric $\langle \xi, \eta \rangle_Q = \operatorname{tr}(\xi^T \eta)$ ($\xi, \eta \in T_Q O(n)$). Then the geodesic from $Q$ with initial velocity $Q\Omega$ ($\Omega \in \operatorname{Skew}(n)$) is

    $$
    \gamma(t) = Q e^{t\Omega}.
    $$

    The geodesic distance between two points $Q_1, Q_2 \in O(n)$ is

    $$
    d(Q_1, Q_2) = \|\log(Q_1^T Q_2)\|_F = \left(\sum_{k=1}^{\lfloor n/2 \rfloor} \theta_k^2\right)^{1/2},
    $$

    where $\theta_k$ are the rotation angles of $Q_1^T Q_2$.

??? proof "Proof"
    Under the bi-invariant metric, the Levi-Civita connection on $O(n)$ is $\nabla_{\tilde{A}}\tilde{B} = \frac{1}{2}[\tilde{A}, \tilde{B}]$ (where $\tilde{A}, \tilde{B}$ are left-invariant vector fields). The velocity of $\gamma(t) = Qe^{t\Omega}$ is $\dot\gamma = Qe^{t\Omega}\Omega$, and the corresponding left-invariant velocity $\gamma^{-1}\dot\gamma = \Omega$ is constant, hence $\nabla_{\dot\gamma}\dot\gamma = 0$, i.e., $\gamma$ is a geodesic. The distance formula follows from the geodesic length $L = \int_0^1 \|\dot\gamma\|dt = \|\Omega\|_F$, where $e^\Omega = Q_1^T Q_2$. $\blacksquare$

!!! example "Example 24.4"
    **Rotations in $SO(3)$.**

    The Lie algebra $\mathfrak{so}(3)$ consists of $3 \times 3$ skew-symmetric matrices, with dimension $3$. The standard basis is

    $$
    L_1 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & 1 & 0 \end{pmatrix}, \quad
    L_2 = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ -1 & 0 & 0 \end{pmatrix}, \quad
    L_3 = \begin{pmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}.
    $$

    Every $\Omega \in \mathfrak{so}(3)$ can be written as $\Omega = \theta_1 L_1 + \theta_2 L_2 + \theta_3 L_3$, corresponding to rotation about axis $(\theta_1, \theta_2, \theta_3)$ by angle $\|\boldsymbol{\theta}\|$. The Rodrigues formula gives

    $$
    e^\Omega = I + \frac{\sin\theta}{\theta}\Omega + \frac{1 - \cos\theta}{\theta^2}\Omega^2, \quad \theta = \|\boldsymbol{\theta}\|.
    $$

!!! example "Example 24.5"
    **Dimension computation for $O(n)$.**

    $O(n)$ is defined by the constraint $Q^T Q = I$, with $n^2$ matrix entries and $\frac{n(n+1)}{2}$ independent constraints (since $Q^T Q$ is symmetric). Therefore

    $$
    \dim O(n) = n^2 - \frac{n(n+1)}{2} = \frac{n(n-1)}{2}.
    $$

    For example, $\dim O(3) = 3$, $\dim O(4) = 6$.

---

## 24.4 Stiefel Manifold

<div class="context-flow" markdown>

**Generalization chain**: $S^{n-1} = \text{St}(1,n) \subset \text{St}(k,n) \subset O(n) = \text{St}(n,n)$ · Tangent space: $X^TZ$ is skew-symmetric → Projection $\Pi_X(Z) = Z - X\text{sym}(X^TZ)$ → Links to S24.8 Riemannian gradient descent

</div>

The Stiefel manifold generalizes the orthogonal group and is an important object in constrained optimization.

!!! definition "Definition 24.10 (Stiefel manifold)"
    The **Stiefel manifold** $V_k(\mathbb{R}^n)$ ($1 \le k \le n$) is the set of all orthonormal $k$-frames in $\mathbb{R}^n$:

    $$
    V_k(\mathbb{R}^n) = \operatorname{St}(k, n) = \{ X \in \mathbb{R}^{n \times k} : X^T X = I_k \}.
    $$

    When $k = n$, $V_n(\mathbb{R}^n) = O(n)$; when $k = 1$, $V_1(\mathbb{R}^n) = S^{n-1}$. $\operatorname{St}(k, n)$ is a compact manifold of dimension $nk - \frac{k(k+1)}{2}$.

!!! theorem "Theorem 24.7 (Tangent space of the Stiefel manifold)"
    The tangent space of $\operatorname{St}(k, n)$ at $X$ is

    $$
    T_X \operatorname{St}(k, n) = \{ Z \in \mathbb{R}^{n \times k} : X^T Z + Z^T X = 0 \}.
    $$

    Equivalently, $Z \in T_X \operatorname{St}(k, n)$ if and only if $X^T Z$ is a skew-symmetric matrix.

??? proof "Proof"
    The constraint is $h(X) = X^T X - I_k = 0$. Differentiating: $Dh(X)(Z) = X^T Z + Z^T X$. Therefore

    $$
    T_X \operatorname{St}(k, n) = \ker(Dh(X)) = \{ Z : X^T Z + Z^T X = 0 \}.
    $$

    Verifying the full-rank condition: $Dh(X)$ is a linear map from $\mathbb{R}^{n \times k}$ to $\operatorname{Sym}(k)$ ($k \times k$ symmetric matrices). For any symmetric matrix $S$, take $Z = \frac{1}{2}XS$; then $X^T Z + Z^T X = S$, so $Dh(X)$ is surjective. $\blacksquare$

!!! theorem "Theorem 24.8 (Canonical metric and geodesics on the Stiefel manifold)"
    Equip $\operatorname{St}(k, n)$ with the metric inherited from Euclidean space $\mathbb{R}^{n \times k}$ (canonical metric) $\langle Z_1, Z_2 \rangle = \operatorname{tr}(Z_1^T Z_2)$. The orthogonal projection onto $T_X \operatorname{St}(k, n)$ is given by

    $$
    \Pi_X(Z) = Z - X \operatorname{sym}(X^T Z),
    $$

    where $\operatorname{sym}(A) = \frac{A + A^T}{2}$.

    Under this metric, the geodesic from $X$ with initial velocity $Z \in T_X \operatorname{St}(k, n)$ can be computed as follows: let $A = X^T Z$ (skew-symmetric), $QR = (I - XX^T)Z$ be the QR decomposition ($Q \in \mathbb{R}^{n \times k}$, $R \in \mathbb{R}^{k \times k}$), then

    $$
    \gamma(t) = \begin{pmatrix} X & Q \end{pmatrix} \exp\!\left(t \begin{pmatrix} A & -R^T \\ R & 0 \end{pmatrix}\right) \begin{pmatrix} I_k \\ 0 \end{pmatrix}.
    $$

??? proof "Proof"
    Projection formula: For $Z \in \mathbb{R}^{n \times k}$, decompose $Z = \Pi_X(Z) + X \cdot \operatorname{sym}(X^T Z)$. Verify that $\Pi_X(Z) \in T_X\operatorname{St}(k,n)$: $X^T\Pi_X(Z) = X^TZ - \operatorname{sym}(X^TZ)$ is skew-symmetric. The residual $X \cdot \operatorname{sym}(X^T Z)$ is orthogonal to the tangent space (in the canonical inner product).

    Deriving the geodesic formula requires solving the connection equation $\nabla_{\dot\gamma}\dot\gamma = 0$. Decomposing the velocity into the $X$-direction and the orthogonal complement of $X$, using matrix exponentials to represent rotations, the geodesic equation reduces to a finite-dimensional matrix ODE whose solution is the formula above. For details, see Edelman-Arias-Smith (1998). $\blacksquare$

!!! example "Example 24.6"
    **Dimension and tangent space of $\operatorname{St}(2, 3)$.**

    $\operatorname{St}(2, 3) = \{X \in \mathbb{R}^{3 \times 2} : X^TX = I_2\}$, dimension $= 3 \times 2 - \frac{2 \times 3}{2} = 3$.

    Take $X = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$. The tangent space consists of $Z \in \mathbb{R}^{3 \times 2}$ with $X^TZ$ skew-symmetric:

    $$
    T_X \operatorname{St}(2, 3) = \left\{ \begin{pmatrix} 0 & -a \\ a & 0 \\ b & c \end{pmatrix} : a, b, c \in \mathbb{R} \right\},
    $$

    which is indeed $3$-dimensional.

---

## 24.5 Grassmann Manifold

<div class="context-flow" markdown>

**Quotient space**: $\text{Gr}(k,n) = \text{St}(k,n)/O(k)$ — parameterizing $k$-dimensional subspaces rather than bases · Distance = **principal angles** $\theta_i = \arccos\sigma_i(X_1^TX_2)$ (SVD, Ch8)

**Applications**: PCA = optimization on $\text{Gr}(k,n)$ (Ch25) · Subspace tracking / computer vision

</div>

The Grassmann manifold parameterizes subspaces of fixed dimension and arises naturally in PCA, subspace tracking, and related problems.

!!! definition "Definition 24.11 (Grassmann manifold)"
    The **Grassmann manifold** $\operatorname{Gr}(k, n)$ is the set of all $k$-dimensional subspaces of $\mathbb{R}^n$:

    $$
    \operatorname{Gr}(k, n) = \{ \mathcal{V} \subseteq \mathbb{R}^n : \dim \mathcal{V} = k \}.
    $$

    $\operatorname{Gr}(k, n)$ can be viewed as a quotient of the Stiefel manifold:

    $$
    \operatorname{Gr}(k, n) \cong \operatorname{St}(k, n) / O(k),
    $$

    where the equivalence relation is $X \sim XQ$ ($Q \in O(k)$), since $X$ and $XQ$ span the same subspace. The dimension of $\operatorname{Gr}(k, n)$ is $k(n - k)$.

!!! definition "Definition 24.12 (Projector representation of the Grassmann manifold)"
    Each $k$-dimensional subspace $\mathcal{V} \in \operatorname{Gr}(k, n)$ can be uniquely represented by its orthogonal projection matrix $P = XX^T$ ($X \in \operatorname{St}(k, n)$ is an orthonormal basis for $\mathcal{V}$). Therefore

    $$
    \operatorname{Gr}(k, n) \cong \{ P \in \mathbb{R}^{n \times n} : P^2 = P, \, P^T = P, \, \operatorname{tr}(P) = k \}.
    $$

!!! theorem "Theorem 24.9 (Tangent space of the Grassmann manifold)"
    In the projector representation $P = XX^T$, the tangent space of $\operatorname{Gr}(k, n)$ at $P$ is

    $$
    T_P \operatorname{Gr}(k, n) = \{ \Delta \in \mathbb{R}^{n \times n} : \Delta = \Delta^T, \, P\Delta + \Delta P = \Delta \}.
    $$

    Equivalently, in the orthonormal basis representation, the horizontal tangent vectors at $T_{[X]} \operatorname{Gr}(k, n)$ are

    $$
    \mathcal{H}_X = \{ Z \in \mathbb{R}^{n \times k} : X^T Z = 0 \},
    $$

    i.e., matrices whose column space is orthogonal to the column space of $X$.

??? proof "Proof"
    In the quotient space framework, tangent vectors $Z \in T_X \operatorname{St}(k, n)$ on the Stiefel manifold decompose into vertical components (along the $O(k)$-orbit direction, i.e., $XA$ with $A$ skew-symmetric) and horizontal components (orthogonal to the orbit). In the canonical metric, the horizontal space is

    $$
    \mathcal{H}_X = \{ Z \in T_X \operatorname{St}(k, n) : X^T Z = 0 \} = \{ Z \in \mathbb{R}^{n \times k} : X^T Z = 0 \}.
    $$

    Its dimension is $k(n - k) = \dim \operatorname{Gr}(k, n)$. $\blacksquare$

!!! theorem "Theorem 24.10 (Geodesics and distance on the Grassmann manifold)"
    Equip $\operatorname{Gr}(k, n)$ with the Riemannian metric inherited from the Stiefel quotient structure. The geodesic distance between two subspaces $\mathcal{V}_1, \mathcal{V}_2 \in \operatorname{Gr}(k, n)$ is

    $$
    d(\mathcal{V}_1, \mathcal{V}_2) = \|\boldsymbol{\theta}\|_2 = \left(\sum_{i=1}^{k} \theta_i^2\right)^{1/2},
    $$

    where $\theta_1, \ldots, \theta_k \in [0, \pi/2]$ are the **principal angles** between $\mathcal{V}_1$ and $\mathcal{V}_2$, defined by $\cos\theta_i = \sigma_i(X_1^T X_2)$, with $\sigma_i$ being the singular values.

    The geodesic from $[X_1]$ in the horizontal direction $Z$ is: let $Z = U\Sigma V^T$ be the compact SVD, then

    $$
    \gamma(t) = [X_1 V \cos(\Sigma t) + U \sin(\Sigma t)].
    $$

??? proof "Proof"
    Let $X_1, X_2 \in \operatorname{St}(k, n)$ be orthonormal bases for $\mathcal{V}_1, \mathcal{V}_2$. The SVD of $X_1^T X_2$ is $X_1^T X_2 = P \operatorname{diag}(\cos\theta_1, \ldots, \cos\theta_k) Q^T$. Take $\tilde{X}_1 = X_1 P$, $\tilde{X}_2 = X_2 Q$, so $\tilde{X}_1^T \tilde{X}_2 = \operatorname{diag}(\cos\theta_i)$.

    Construct the geodesic $\gamma(t) = \tilde{X}_1 \operatorname{diag}(\cos(t\theta_i)) + U_\perp \operatorname{diag}(\sin(t\theta_i))$, where $U_\perp = (\tilde{X}_2 - \tilde{X}_1 \operatorname{diag}(\cos\theta_i)) \operatorname{diag}(\sin\theta_i)^{-1}$. Verify that $\gamma(0) = [\tilde{X}_1] = [X_1]$, $\gamma(1) = [X_2]$, and $\gamma$ satisfies the geodesic equation. The length is $\int_0^1 \|\dot\gamma\| dt = \|\boldsymbol{\theta}\|_2$. $\blacksquare$

!!! example "Example 24.7"
    **$\operatorname{Gr}(1, n)$ — real projective space.**

    $\operatorname{Gr}(1, n)$ is the set of all one-dimensional subspaces (lines through the origin) in $\mathbb{R}^n$, i.e., the real projective space $\mathbb{RP}^{n-1}$. Dimension: $1 \times (n-1) = n - 1$. The geodesic distance between two lines is their angle $\theta \in [0, \pi/2]$.

!!! example "Example 24.8"
    **Computing principal angles between two subspaces.**

    Let $X_1 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix}$, $X_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}$ (already orthogonalized). Compute

    $$
    X_1^T X_2 = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix},
    $$

    with singular values $\sigma_1 = \sigma_2 = 1/\sqrt{2}$, principal angles $\theta_1 = \theta_2 = \pi/4$. Geodesic distance $d = \sqrt{(\pi/4)^2 + (\pi/4)^2} = \frac{\pi}{2\sqrt{2}}$.

---

## 24.6 Manifold of Positive Definite Matrices

<div class="context-flow" markdown>

**Special geometry**: Affine-invariant metric on $\mathcal{P}(n)$: $\langle \xi,\eta\rangle_P = \text{tr}(P^{-1}\xi P^{-1}\eta)$ → Geodesic $P^{1/2}e^{tP^{-1/2}\xi P^{-1/2}}P^{1/2}$ → **Hadamard manifold** (nonpositive curvature, unique Frechet mean)

**Applications**: Covariance estimation, diffusion tensor imaging, geometric means · $\mathcal{P}(1) = (0,\infty)$ with distance $= |\ln(p/q)|$

</div>

Symmetric positive definite matrices play important roles in statistics, diffusion tensor imaging, and covariance estimation.

!!! definition "Definition 24.13 (Manifold of positive definite matrices)"
    The **manifold of symmetric positive definite matrices** is defined as

    $$
    \mathcal{P}(n) = \operatorname{Sym}^+(n) = \{ P \in \mathbb{R}^{n \times n} : P = P^T, \, P \succ 0 \}.
    $$

    $\mathcal{P}(n)$ is an open subset of $\operatorname{Sym}(n)$, hence a $\frac{n(n+1)}{2}$-dimensional smooth manifold with tangent space $T_P \mathcal{P}(n) = \operatorname{Sym}(n)$.

!!! definition "Definition 24.14 (Affine-invariant Riemannian metric)"
    The **affine-invariant Riemannian metric** on $\mathcal{P}(n)$ is defined as

    $$
    \langle \xi, \eta \rangle_P = \operatorname{tr}(P^{-1} \xi P^{-1} \eta), \quad \xi, \eta \in T_P \mathcal{P}(n) = \operatorname{Sym}(n).
    $$

    This metric is invariant under the group action $P \mapsto APA^T$ ($A \in GL(n)$).

<div class="context-flow" markdown>

**Insight**: Under the affine-invariant metric, $d(P,Q) = \|\log(P^{-1/2}QP^{-1/2})\|_F$ — the distance is determined by the **log-ratios** of eigenvalues; the geometric mean replaces the arithmetic mean as the natural "center"

</div>

!!! theorem "Theorem 24.11 (Geodesics on $\mathcal{P}(n)$)"
    Under the affine-invariant metric, the geodesic on $\mathcal{P}(n)$ starting from $P$ with initial velocity $\xi \in \operatorname{Sym}(n)$ is

    $$
    \gamma(t) = P^{1/2} \exp(t P^{-1/2} \xi P^{-1/2}) P^{1/2}.
    $$

    The geodesic distance between $P, Q \in \mathcal{P}(n)$ is

    $$
    d(P, Q) = \left\| \log(P^{-1/2} Q P^{-1/2}) \right\|_F = \left( \sum_{i=1}^{n} \ln^2 \lambda_i \right)^{1/2},
    $$

    where $\lambda_1, \ldots, \lambda_n$ are the eigenvalues of $P^{-1}Q$ (or equivalently $P^{-1/2}QP^{-1/2}$).

??? proof "Proof"
    At $P = I$, the metric simplifies to $\langle \xi, \eta \rangle_I = \operatorname{tr}(\xi\eta)$, and the geodesic is $\gamma(t) = e^{t\xi}$. Verification: $\gamma(0) = I$, $\dot\gamma(0) = \xi$, and $\gamma(t)$ satisfies the geodesic equation $\ddot\gamma - \dot\gamma \gamma^{-1} \dot\gamma = 0$ on $\mathcal{P}(n)$ (the equation of motion is derived from the Levi-Civita connection $\nabla_\xi \eta = D_\xi \eta - \frac{1}{2}(\xi P^{-1}\eta + \eta P^{-1}\xi)$).

    For general $P$, use the isometry $\Phi_A(P) = APA^T$ to map $P$ to $I$ (take $A = P^{-1/2}$), yielding the general geodesic formula. The distance is

    $$
    d(P, Q) = d(I, P^{-1/2}QP^{-1/2}) = \|\log(P^{-1/2}QP^{-1/2})\|_F.
    $$

    $\blacksquare$

!!! theorem "Theorem 24.12 (Frechet mean on $\mathcal{P}(n)$)"
    Let $P_1, \ldots, P_m \in \mathcal{P}(n)$. Their **Frechet mean** with respect to the affine-invariant metric is

    $$
    \bar{P} = \arg\min_{P \in \mathcal{P}(n)} \sum_{i=1}^{m} d(P, P_i)^2.
    $$

    The Frechet mean exists and is unique (because $\mathcal{P}(n)$ under the affine-invariant metric is a nonpositively curved complete Riemannian manifold, i.e., a Hadamard manifold).

??? proof "Proof"
    $\mathcal{P}(n)$ under the affine-invariant metric is complete (geodesics exist for all time) and has sectional curvature $K \le 0$ (verifiable by computing the curvature tensor). On Hadamard manifolds, the Frechet mean exists and is unique, a consequence of the Cartan-Hadamard theorem and Karcher's theorem. Specifically, strictly nonpositive curvature ensures that the objective function $f(P) = \sum_i d(P, P_i)^2$ is strictly convex (in the geodesic sense), hence has a unique minimizer. $\blacksquare$

!!! example "Example 24.9"
    **$\mathcal{P}(1)$, i.e., positive real numbers.**

    $\mathcal{P}(1) = (0, \infty)$, with metric $ds^2 = dp^2/p^2$ (hyperbolic metric). The geodesic is $\gamma(t) = p_0 e^{vt}$, distance is $d(p, q) = |\ln(p/q)|$. The Frechet mean of two points $p, q > 0$ is the geometric mean $\bar{p} = (p \cdot q)^{1/2}$.

!!! example "Example 24.10"
    **Computing distance between $2 \times 2$ positive definite matrices.**

    Let $P = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$, $Q = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix}$. Then

    $$
    P^{-1}Q = \begin{pmatrix} 1/2 & 0 \\ 0 & 2 \end{pmatrix}, \quad \log(P^{-1}Q) = \begin{pmatrix} -\ln 2 & 0 \\ 0 & \ln 2 \end{pmatrix}.
    $$

    Distance $d(P, Q) = \sqrt{(\ln 2)^2 + (\ln 2)^2} = \ln 2 \cdot \sqrt{2} \approx 0.980$.

---

## 24.7 Matrix Lie Groups and Lie Algebras

<div class="context-flow" markdown>

**Unified theory**: BCH formula $e^Xe^Y = e^{X+Y+\frac{1}{2}[X,Y]+\cdots}$ encodes group multiplication as Lie brackets → Lie group homomorphism $\leftrightarrow$ Lie algebra homomorphism (full correspondence when simply connected)

**Practical**: Rodrigues formula in $SO(3)$ · Composition of small rotations in robotics / computer vision

</div>

This section systematically discusses the general theory of matrix Lie groups.

!!! definition "Definition 24.15 (Matrix Lie group)"
    A **matrix Lie group** is a closed subgroup $G$ of $GL(n, \mathbb{C})$. By Cartan's closed subgroup theorem, $G$ is automatically a smooth manifold and hence a Lie group. Common matrix Lie groups include $GL(n), SL(n), O(n), SO(n), U(n), SU(n), Sp(2n)$, etc.

!!! theorem "Theorem 24.13 (Baker-Campbell-Hausdorff formula)"
    Let $X, Y \in \mathfrak{g}$ (some Lie algebra) with $\|X\|, \|Y\|$ sufficiently small. Then there exists $Z \in \mathfrak{g}$ such that $e^X e^Y = e^Z$, and

    $$
    Z = X + Y + \frac{1}{2}[X, Y] + \frac{1}{12}\big([X, [X, Y]] - [Y, [X, Y]]\big) + \cdots
    $$

    This series is called the **Baker-Campbell-Hausdorff (BCH) formula**, where each term is a nested Lie bracket of $X$ and $Y$. In particular, when $[X, Y] = 0$, $Z = X + Y$.

??? proof "Proof"
    **Proof sketch.** Taking the logarithm of both sides of $e^X e^Y = e^Z$ (in the formal power series sense), using the series expansions

    $$
    e^X = \sum_k \frac{X^k}{k!}, \quad \log(I + W) = \sum_k \frac{(-1)^{k+1}}{k} W^k,
    $$

    one expands $Z = \log(e^X e^Y)$ as a power series in $X, Y$. By Dynkin's explicit formula, each term can be expressed as nested commutators.

    Verification of the first few terms:

    $$
    e^X e^Y = (I + X + \tfrac{X^2}{2} + \cdots)(I + Y + \tfrac{Y^2}{2} + \cdots) = I + (X+Y) + (XY + \tfrac{X^2 + Y^2}{2}) + \cdots
    $$

    $$
    \log(e^X e^Y) = (X+Y) + \frac{XY - YX}{2} + \cdots = X + Y + \frac{1}{2}[X, Y] + \cdots
    $$

    Convergence is guaranteed when $\|X\| + \|Y\|$ is sufficiently small by completeness of the Lie algebra. $\blacksquare$

!!! theorem "Theorem 24.14 (Correspondence between Lie group and Lie algebra homomorphisms)"
    Let $G, H$ be matrix Lie groups with Lie algebras $\mathfrak{g}, \mathfrak{h}$. If $\Phi: G \to H$ is a Lie group homomorphism (smooth group homomorphism), then $d\Phi_I: \mathfrak{g} \to \mathfrak{h}$ is a Lie algebra homomorphism, i.e., it preserves Lie brackets:

    $$
    d\Phi_I([X, Y]) = [d\Phi_I(X), d\Phi_I(Y)], \quad \forall X, Y \in \mathfrak{g}.
    $$

    Conversely, if $G$ is simply connected, then every Lie algebra homomorphism $\phi: \mathfrak{g} \to \mathfrak{h}$ lifts to a unique Lie group homomorphism $\Phi: G \to H$ with $d\Phi_I = \phi$.

??? proof "Proof"
    Let $\Phi: G \to H$ be a Lie group homomorphism. For $X, Y \in \mathfrak{g}$,

    $$
    d\Phi_I([X, Y]) = d\Phi_I\!\left(\frac{d}{dt}\Big|_{t=0} e^{tX} Y e^{-tX}\right) = \frac{d}{dt}\Big|_{t=0} \Phi(e^{tX}) d\Phi_I(Y) \Phi(e^{-tX}).
    $$

    Since $\Phi(e^{tX}) = e^{t \, d\Phi_I(X)}$, this equals $[d\Phi_I(X), d\Phi_I(Y)]$.

    The converse direction uses the simple connectivity of $G$ and the local diffeomorphism property of the exponential map to extend $\phi$ step by step to the entire group. $\blacksquare$

!!! example "Example 24.11"
    **Lie algebra map of $\det: GL(n) \to \mathbb{R}^*$.**

    The determinant $\det: GL(n) \to GL(1) = \mathbb{R}^*$ is a Lie group homomorphism. Its differential at $I$ is

    $$
    d(\det)_I(X) = \operatorname{tr}(X).
    $$

    This is because $\det(I + tX) = 1 + t\operatorname{tr}(X) + O(t^2)$. The corresponding Lie algebra homomorphism is $\operatorname{tr}: \mathfrak{gl}(n) \to \mathfrak{gl}(1) = \mathbb{R}$, preserving Lie brackets: $\operatorname{tr}([A, B]) = 0 = [\operatorname{tr}(A), \operatorname{tr}(B)]$.

!!! example "Example 24.12"
    **Application of the BCH formula: approximate products.**

    For small-angle rotations $R_1 = e^{\Omega_1}, R_2 = e^{\Omega_2} \in SO(3)$, the BCH formula gives

    $$
    R_1 R_2 \approx \exp\!\left(\Omega_1 + \Omega_2 + \frac{1}{2}[\Omega_1, \Omega_2]\right).
    $$

    This is used in robotics and computer vision for composing small rotations. When $\Omega_1, \Omega_2$ are sufficiently small, one can further approximate $R_1 R_2 \approx e^{\Omega_1 + \Omega_2}$.

---

## 24.8 Optimization on Matrix Manifolds

<div class="context-flow" markdown>

**Core transformation**: Constrained optimization ($X^TX=I$, etc.) → Unconstrained optimization on manifolds

**Riemannian gradient** = Euclidean gradient projected onto the tangent space → **Retraction** replaces the expensive geodesic step

**Links**: Eigenvalue problem = $\max \text{tr}(X^TAX)$ on $\text{St}(k,n)$ (Ch25) · PCA = optimization on $\text{Gr}(k,n)$

</div>

Many practical problems can be modeled as optimization on matrix manifolds. Manifold optimization transforms constrained optimization into unconstrained Riemannian optimization.

!!! definition "Definition 24.16 (Riemannian gradient)"
    Let $f: \mathcal{M} \to \mathbb{R}$ be a smooth function on a Riemannian manifold $(\mathcal{M}, \langle \cdot, \cdot \rangle)$. The **Riemannian gradient** of $f$ at $X \in \mathcal{M}$, $\operatorname{grad} f(X) \in T_X \mathcal{M}$, is the unique tangent vector satisfying

    $$
    \langle \operatorname{grad} f(X), \xi \rangle_X = Df(X)[\xi], \quad \forall \xi \in T_X \mathcal{M}.
    $$

!!! definition "Definition 24.17 (Retraction)"
    A **retraction** $R_X: T_X \mathcal{M} \to \mathcal{M}$ is a first-order approximation of the exponential map, satisfying:

    1. $R_X(0) = X$;
    2. $\frac{d}{dt}\Big|_{t=0} R_X(t\xi) = \xi$, $\forall \xi \in T_X \mathcal{M}$.

    Retractions are computationally cheaper than the exponential map and are commonly used in optimization algorithms as substitutes for geodesic steps.

!!! theorem "Theorem 24.15 (Convergence of Riemannian gradient descent)"
    Let $f: \mathcal{M} \to \mathbb{R}$ be a smooth function on a compact Riemannian manifold and $R$ a retraction. The Riemannian gradient descent iteration

    $$
    X_{k+1} = R_{X_k}(-\alpha_k \operatorname{grad} f(X_k))
    $$

    satisfies: under appropriate step size selection (e.g., Armijo line search), $\|\operatorname{grad} f(X_k)\| \to 0$, i.e., the iterates approach a critical point.

??? proof "Proof"
    **Proof sketch.** Using properties of the retraction and Taylor expansion:

    $$
    f(R_X(-\alpha \operatorname{grad} f)) = f(X) - \alpha \|\operatorname{grad} f(X)\|^2 + O(\alpha^2).
    $$

    Under the Armijo condition $f(X_{k+1}) \le f(X_k) - c \alpha_k \|\operatorname{grad} f(X_k)\|^2$ ($0 < c < 1$),

    $$
    \sum_{k=0}^{\infty} \alpha_k \|\operatorname{grad} f(X_k)\|^2 \le f(X_0) - \inf f < \infty.
    $$

    By the lower bound on step sizes (backtracking line search guarantees $\alpha_k \ge \alpha_{\min} > 0$), we obtain $\|\operatorname{grad} f(X_k)\| \to 0$. Compactness ensures the existence of lower bounds and Lipschitz constants. $\blacksquare$

!!! theorem "Theorem 24.16 (Riemannian gradient on the Stiefel manifold)"
    Let $f: \mathbb{R}^{n \times k} \to \mathbb{R}$ be a smooth function restricted to $\operatorname{St}(k, n)$. Under the canonical metric, the Riemannian gradient is

    $$
    \operatorname{grad} f(X) = \nabla f(X) - X \operatorname{sym}(X^T \nabla f(X)),
    $$

    where $\nabla f(X)$ is the Euclidean gradient of $f$ and $\operatorname{sym}(A) = (A + A^T)/2$.

??? proof "Proof"
    The Riemannian gradient is the orthogonal projection of the Euclidean gradient onto the tangent space: $\operatorname{grad} f(X) = \Pi_X(\nabla f(X))$. By the projection formula from Theorem 24.8, $\Pi_X(Z) = Z - X\operatorname{sym}(X^TZ)$. Substituting $Z = \nabla f(X)$ gives the result. $\blacksquare$

!!! example "Example 24.13"
    **Eigenvalue problem on the Stiefel manifold.**

    Finding the $k$ largest eigenvalues and corresponding eigenvectors of $A \in \operatorname{Sym}(n)$ is equivalent to

    $$
    \max_{X \in \operatorname{St}(k, n)} \operatorname{tr}(X^T A X).
    $$

    Euclidean gradient $\nabla f(X) = 2AX$, Riemannian gradient:

    $$
    \operatorname{grad} f(X) = 2AX - X \operatorname{sym}(X^T \cdot 2AX) = 2AX - X(X^TAX + (X^TAX)^T)/1 = 2(I - XX^T)AX,
    $$

    using the symmetry of $X^TAX$. Riemannian gradient ascent with a retraction (e.g., QR retraction $R_X(Z) = \operatorname{qf}(X + Z)$, where $\operatorname{qf}$ denotes the Q factor of the QR decomposition) can solve this problem.

!!! example "Example 24.14"
    **Subspace fitting on the Grassmann manifold.**

    Given a data matrix $Y \in \mathbb{R}^{n \times m}$, find the best $k$-dimensional subspace to maximize projected variance:

    $$
    \max_{[X] \in \operatorname{Gr}(k, n)} \|X^T Y\|_F^2 = \max_{[X] \in \operatorname{Gr}(k, n)} \operatorname{tr}(X^T Y Y^T X).
    $$

    This is the manifold optimization formulation of PCA. The Riemannian gradient is $\operatorname{grad} f([X]) = 2(I - XX^T)YY^TX$, and the retraction can be the polar decomposition retraction $R_X(Z) = (X + Z)(I + Z^TZ)^{-1/2}$.

!!! example "Example 24.15"
    **Geometric mean on the positive definite matrix manifold.**

    Given $P_1, \ldots, P_m \in \mathcal{P}(n)$, the geometric mean (Frechet mean) can be computed via Riemannian gradient descent:

    $$
    \operatorname{grad} f(P) = -\sum_{i=1}^{m} \log(P^{-1} P_i),
    $$

    where $f(P) = \frac{1}{2}\sum_i d(P, P_i)^2$ and $\log$ is the matrix logarithm. The iteration

    $$
    P_{k+1} = P_k^{1/2} \exp\!\left(\frac{\alpha}{m} P_k^{-1/2} \sum_{i=1}^{m} \log(P_k^{-1/2} P_i P_k^{-1/2}) P_k^{-1/2}\right) P_k^{1/2}
    $$

    converges to the geometric mean under appropriate step size $\alpha$.

---

## Chapter Summary

This chapter systematically introduced the theoretical framework and computational methods for matrix manifolds:

1. **Basic concepts**: Embedded submanifolds, tangent spaces, and Riemannian metrics provide geometric language for matrix constraint problems.
2. The **general linear group** $GL(n)$ is the parent group of all matrix Lie groups; its Lie algebra $\mathfrak{gl}(n)$ is the full matrix space.
3. The **orthogonal and unitary groups** $O(n), U(n)$ and their special subgroups are the most fundamental compact Lie groups, with Lie algebras being skew-symmetric and skew-Hermitian matrix spaces respectively.
4. The **Stiefel manifold** $\operatorname{St}(k, n)$ parameterizes orthonormal $k$-frames and is the natural space for orthogonal constraint optimization.
5. The **Grassmann manifold** $\operatorname{Gr}(k, n)$ parameterizes $k$-dimensional subspaces; principal angles provide a natural distance between subspaces.
6. The **positive definite matrix manifold** $\mathcal{P}(n)$ under the affine-invariant metric is a nonpositively curved space with a unique Frechet mean.
7. **Matrix Lie groups and Lie algebras** are connected through the exponential map and the BCH formula.
8. **Manifold optimization** transforms constrained optimization into unconstrained Riemannian optimization; the Riemannian gradient and retraction are core tools.
