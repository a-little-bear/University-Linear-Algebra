# Chapter 23  Introduction to Random Matrices

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues/SVD (Ch6-8) · Positive definite matrices (Ch7) · **Arc**: Statistical behavior of eigenvalues as $n \to \infty$ — Semicircle law (Wigner) → MP law (Wishart/sample covariance) → Tracy-Widom (edge fluctuations) → Universality
**Essence**: Deterministic matrix theory concerns individual eigenvalues; random matrix theory concerns the **collective behavior** of $n$ eigenvalues — a leap from algebra to probability

</div>

Random matrix theory (RMT) studies the statistical properties of eigenvalues and eigenvectors when matrix entries are random variables. The theory originated in the 1950s from Wigner's study of nuclear energy level statistics, and has since had profound impact in mathematical physics, number theory, wireless communications, high-dimensional statistics, and many other fields. This chapter systematically introduces the basic concepts, core limit theorems, and several frontier applications of random matrices.

---

## 23.1 Basic Concepts of Random Matrices

<div class="context-flow" markdown>

**Three main ensembles**: GOE ($\beta=1$, real symmetric) / GUE ($\beta=2$, Hermitian) / GSE ($\beta=4$, quaternion) → Symmetry determines $\beta$ → Joint density $\propto \prod_{i<j}|\lambda_i - \lambda_j|^\beta \cdot e^{-\sum \lambda_i^2}$
**Wishart matrix** $W = \frac{1}{n}X^TX$: prototype of the sample covariance matrix → Links to Ch25 PCA

</div>

A random matrix is a matrix whose entries are random variables. The core question is: when the matrix dimension $n \to \infty$, what deterministic limit does the empirical spectral distribution exhibit?

!!! definition "Definition 23.1 (Random matrix)"
    Let $(\Omega, \mathcal{F}, P)$ be a probability space. A **random matrix** is a measurable map $M: \Omega \to \mathbb{K}^{n \times n}$, where $\mathbb{K} = \mathbb{R}$ or $\mathbb{C}$, i.e., each entry $M_{ij}(\omega)$ is a random variable defined on this probability space.

!!! definition "Definition 23.2 (Gaussian Orthogonal Ensemble GOE)"
    The **Gaussian Orthogonal Ensemble** (GOE) is the probability distribution of $n \times n$ real symmetric random matrices $M$ with density

    $$
    f(M) = C_n \exp\!\left( -\frac{n}{4} \operatorname{tr}(M^2) \right),
    $$

    where $C_n$ is a normalization constant. Equivalently, the upper triangular entries of $M$ are independent, with diagonal entries $M_{ii} \sim N(0, 2/n)$, off-diagonal entries $M_{ij} \sim N(0, 1/n)$ ($i < j$), and $M_{ji} = M_{ij}$. The GOE distribution is invariant under orthogonal conjugation $M \mapsto O^T M O$ ($O \in O(n)$).

!!! definition "Definition 23.3 (Gaussian Unitary Ensemble GUE)"
    The **Gaussian Unitary Ensemble** (GUE) is the probability distribution of $n \times n$ Hermitian random matrices $M$ with density

    $$
    f(M) = \widetilde{C}_n \exp\!\left( -\frac{n}{2} \operatorname{tr}(M^2) \right).
    $$

    Equivalently, diagonal entries $M_{ii} \sim N(0, 1/n)$ are real; for off-diagonal entries ($i < j$), the real and imaginary parts are independent, each distributed as $N(0, 1/(2n))$, and $M_{ji} = \overline{M_{ij}}$. The GUE distribution is invariant under unitary conjugation $M \mapsto U^* M U$ ($U \in U(n)$).

!!! definition "Definition 23.4 (Gaussian Symplectic Ensemble GSE)"
    The **Gaussian Symplectic Ensemble** (GSE) is the probability distribution of $2n \times 2n$ self-dual quaternion Hermitian matrices. Its distribution is invariant under symplectic conjugation $M \mapsto S^* M S$ ($S \in Sp(2n)$), with density

    $$
    f(M) = \hat{C}_n \exp\!\left( -n \operatorname{tr}(M^2) \right).
    $$

    GOE, GUE, and GSE correspond to Dyson indices $\beta = 1, 2, 4$ respectively.

!!! definition "Definition 23.5 (Wishart matrix)"
    Let $X$ be an $n \times p$ matrix whose rows are independently distributed as $N(\mathbf{0}, \Sigma)$. The **Wishart matrix** is defined as

    $$
    W = \frac{1}{n} X^T X.
    $$

    When $\Sigma = I_p$, $W$ is called a white Wishart matrix. The Wishart distribution is denoted $W \sim \mathcal{W}_p(\Sigma, n)$.

!!! definition "Definition 23.6 (Empirical spectral distribution)"
    Let $M$ be an $n \times n$ Hermitian matrix with eigenvalues $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$. The **empirical spectral distribution** (ESD) is defined as

    $$
    F_n(x) = \frac{1}{n} \#\{ i : \lambda_i \le x \} = \frac{1}{n} \sum_{i=1}^{n} \mathbf{1}_{\{\lambda_i \le x\}}.
    $$

    The corresponding empirical spectral measure is $\mu_n = \frac{1}{n} \sum_{i=1}^{n} \delta_{\lambda_i}$.

!!! note "Note"
    The unified framework of the three Gaussian ensembles can be given by the **$\beta$-ensemble**: for $\beta > 0$, the joint eigenvalue density is

    $$
    p(\lambda_1, \ldots, \lambda_n) = Z_{n,\beta}^{-1} \prod_{i < j} |\lambda_i - \lambda_j|^\beta \prod_{i=1}^{n} e^{-\frac{n\beta}{4} \lambda_i^2}.
    $$

!!! theorem "Theorem 23.1 (Joint eigenvalue density of GOE/GUE)"
    Let $M$ be a GOE ($\beta=1$) or GUE ($\beta=2$) matrix. The joint probability density of its eigenvalues $\lambda_1, \ldots, \lambda_n$ is

    $$
    p(\lambda_1, \ldots, \lambda_n) = Z_{n,\beta}^{-1} \prod_{1 \le i < j \le n} |\lambda_i - \lambda_j|^\beta \cdot \prod_{i=1}^{n} e^{-\frac{n\beta}{4} \lambda_i^2},
    $$

    where $Z_{n,\beta}$ is the normalization constant.

??? proof "Proof"
    We illustrate with GUE ($\beta = 2$). Let $M = U \Lambda U^*$, where $U \in U(n)$, $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$. Performing a change of variables on the density $f(M) \propto e^{-\frac{n}{2}\operatorname{tr}(M^2)}$, by the Jacobi formula the Jacobian determinant is

    $$
    |J| = \prod_{i < j} (\lambda_i - \lambda_j)^2.
    $$

    Note that $\operatorname{tr}(M^2) = \operatorname{tr}(\Lambda^2) = \sum_i \lambda_i^2$. Integrating over the Haar measure on $U$ yields

    $$
    p(\lambda_1, \ldots, \lambda_n) \propto \prod_{i < j} (\lambda_i - \lambda_j)^2 \cdot \prod_i e^{-\frac{n}{2}\lambda_i^2}.
    $$

    This is the $\beta = 2$ case. The derivation for GOE ($\beta=1$) is similar, with the exponent in the Jacobian being $1$. $\blacksquare$

!!! example "Example 23.1"
    **Verifying the $2 \times 2$ GUE eigenvalue density.**

    Let $M = \begin{pmatrix} a & z \\ \bar{z} & b \end{pmatrix}$, where $a, b \sim N(0, 1/n)$, $z = x + iy$, $x, y \sim N(0, 1/(2n))$. The eigenvalues are

    $$
    \lambda_\pm = \frac{a+b}{2} \pm \sqrt{\left(\frac{a-b}{2}\right)^2 + |z|^2}.
    $$

    Performing the change of variables $(a, b, x, y) \to (\lambda_+, \lambda_-, \theta)$, where $\theta$ is the eigenvector parameter, one can verify that the factor $|\lambda_+ - \lambda_-|^2$ appears in the joint density, consistent with Theorem 23.1.

---

## 23.2 Wigner Matrices and the Semicircle Law

<div class="context-flow" markdown>

**Core theorem**: Wigner matrices (symmetric + independent + zero mean + variance 1) have empirical spectral distribution → $\rho_{sc}(x) = \frac{1}{2\pi}\sqrt{4-x^2}$ (semicircle)
**Proof route**: Method of moments — $\frac{1}{n}\mathbb{E}[\text{tr}(W^{2m})]$ → Closed path counting → **Catalan numbers** $C_m$ → Uniquely determine the semicircle distribution

</div>

The Wigner semicircle law is the most fundamental limit theorem in random matrix theory, describing the limiting behavior of the empirical spectral distribution of Wigner matrices.

!!! definition "Definition 23.7 (Wigner matrix)"
    A **Wigner matrix** is an $n \times n$ Hermitian (or real symmetric) random matrix $W_n = \frac{1}{\sqrt{n}}(X_{ij})_{1 \le i,j \le n}$, where:

    1. $\{X_{ij} : i \le j\}$ are independent;
    2. Diagonal entries $X_{ii}$ are i.i.d., $\mathbb{E}[X_{ii}] = 0$, $\mathbb{E}[X_{ii}^2] < \infty$;
    3. Off-diagonal entries $X_{ij}$ ($i < j$) are i.i.d., $\mathbb{E}[X_{ij}] = 0$, $\mathbb{E}[|X_{ij}|^2] = 1$;
    4. $X_{ji} = \overline{X_{ij}}$.

!!! theorem "Theorem 23.2 (Wigner semicircle law)"
    Let $W_n$ be a Wigner matrix, and let $\mu_n = \frac{1}{n}\sum_{i=1}^{n} \delta_{\lambda_i}$ be its empirical spectral measure. Then as $n \to \infty$, $\mu_n$ converges weakly almost surely to the **semicircle distribution** $\mu_{sc}$, whose density is

    $$
    \rho_{sc}(x) = \frac{1}{2\pi} \sqrt{4 - x^2} \cdot \mathbf{1}_{[-2, 2]}(x).
    $$

??? proof "Proof"
    The core idea of the **method of moments** is as follows.

    **Step 1: Compute moments.** We need to prove that for every positive integer $k$,

    $$
    \frac{1}{n} \mathbb{E}\!\left[\operatorname{tr}(W_n^k)\right] \to m_k = \int x^k \, \rho_{sc}(x) \, dx.
    $$

    Expanding $\operatorname{tr}(W_n^k) = \frac{1}{n^{k/2}} \sum_{i_1, \ldots, i_k} X_{i_1 i_2} X_{i_2 i_3} \cdots X_{i_k i_1}$, each term corresponds to a closed path $(i_1, i_2, \ldots, i_k, i_1)$.

    **Step 2: Combinatorial graph-theoretic analysis.** Since $\mathbb{E}[X_{ij}] = 0$, only paths where every edge is traversed at least twice contribute. For odd $k$, no such paths exist, so $m_k = 0$. For $k = 2m$, the main contribution comes from paths where each edge is traversed exactly twice, and such paths are in bijection with **Catalan numbers** $C_m$.

    **Step 3: Identify Catalan numbers.** One can show that

    $$
    m_{2m} = C_m = \frac{1}{m+1}\binom{2m}{m},
    $$

    and the $2m$-th moment of the semicircle distribution is exactly the Catalan number.

    **Step 4: Moments uniquely determine the distribution.** Since the semicircle distribution has bounded support ($[-2, 2]$), its moment sequence uniquely determines the distribution.

    **Step 5: Almost sure convergence.** Using the variance estimate $\operatorname{Var}\!\left(\frac{1}{n}\operatorname{tr}(W_n^k)\right) = O(n^{-2})$ and the Borel-Cantelli lemma, convergence in expectation can be upgraded to almost sure convergence. $\blacksquare$

!!! theorem "Theorem 23.3 (Stieltjes transform characterization of the semicircle law)"
    The Stieltjes transform of the semicircle distribution $\mu_{sc}$ is

    $$
    s(z) = \int \frac{\rho_{sc}(x)}{x - z} \, dx = \frac{-z + \sqrt{z^2 - 4}}{2}, \quad z \in \mathbb{C}^+,
    $$

    where the branch is chosen so that $\operatorname{Im}(s(z)) > 0$ when $\operatorname{Im}(z) > 0$. Equivalently, $s(z)$ satisfies the equation

    $$
    s(z)^2 + z \, s(z) + 1 = 0.
    $$

??? proof "Proof"
    Direct computation. Let $z \in \mathbb{C}^+$, then

    $$
    s(z) = \frac{1}{2\pi} \int_{-2}^{2} \frac{\sqrt{4 - x^2}}{x - z} \, dx.
    $$

    Substituting $x = 2\cos\theta$, $dx = -2\sin\theta \, d\theta$, $\sqrt{4 - x^2} = 2\sin\theta$, we get

    $$
    s(z) = \frac{1}{2\pi} \int_0^{\pi} \frac{4\sin^2\theta}{2\cos\theta - z} \, (-2) \, d\theta \cdot \frac{1}{-1} = \frac{2}{\pi} \int_0^{\pi} \frac{\sin^2\theta}{z - 2\cos\theta} \, d\theta.
    $$

    Using the residue theorem (setting $w = e^{i\theta}$ and converting to a contour integral), one computes

    $$
    s(z) = \frac{-z + \sqrt{z^2 - 4}}{2}.
    $$

    Verification: $s^2 + zs + 1 = \frac{z^2 - 2z\sqrt{z^2-4} + z^2 - 4}{4} + \frac{-z^2 + z\sqrt{z^2-4}}{2} + 1 = 0$. $\blacksquare$

!!! example "Example 23.2"
    **Numerical verification of the semicircle law.**

    Take $n = 1000$, generate a GOE matrix $M = \frac{1}{\sqrt{n}} A$, where $A$ is a symmetric matrix with independent standard normal upper triangular entries. Compute the eigenvalues and plot a histogram; observe the fit with $\rho_{sc}(x) = \frac{1}{2\pi}\sqrt{4 - x^2}$. Experiments show that even at $n = 1000$, the empirical spectral distribution is already very close to the semicircle distribution.

!!! example "Example 23.3"
    **Computing moments of the semicircle distribution.**

    Odd moments of the semicircle distribution are zero (by symmetry). Even moments:

    $$
    m_2 = \int_{-2}^{2} x^2 \cdot \frac{\sqrt{4 - x^2}}{2\pi} \, dx = 1, \quad m_4 = \int_{-2}^{2} x^4 \cdot \frac{\sqrt{4 - x^2}}{2\pi} \, dx = 2.
    $$

    In general, $m_{2k} = C_k = \frac{1}{k+1}\binom{2k}{k}$, where $C_k$ is the $k$-th Catalan number. First few values: $C_0=1, C_1=1, C_2=2, C_3=5, C_4=14$.

---

## 23.3 Sample Covariance Matrices and the Marchenko-Pastur Law

<div class="context-flow" markdown>

**From Wigner to Wishart**: Semicircle law = spectral limit of symmetric matrices → **MP law** = spectral limit of $\frac{1}{n}X^TX$, supported on $[(1-\sqrt{y})^2, (1+\sqrt{y})^2]$, $y = p/n$
**Cornerstone of high-dimensional statistics**: When $p/n \to y > 0$, the sample covariance matrix deviates enormously from the true covariance → Classical statistical theory breaks down → Links to Ch25 PCA

</div>

When we draw samples from a high-dimensional population, the spectral behavior of the sample covariance matrix is characterized by the Marchenko-Pastur law.

!!! definition "Definition 23.8 (Sample covariance matrix)"
    Let $\mathbf{x}_1, \ldots, \mathbf{x}_n \in \mathbb{R}^p$ be i.i.d. random vectors with $\mathbb{E}[\mathbf{x}_i] = \mathbf{0}$, $\operatorname{Cov}(\mathbf{x}_i) = \Sigma$. The **sample covariance matrix** is defined as

    $$
    S_n = \frac{1}{n} \sum_{i=1}^{n} \mathbf{x}_i \mathbf{x}_i^T = \frac{1}{n} X^T X,
    $$

    where $X = (\mathbf{x}_1, \ldots, \mathbf{x}_n)^T$ is the $n \times p$ data matrix.

!!! theorem "Theorem 23.4 (Marchenko-Pastur law)"
    Let $X$ be an $n \times p$ matrix with i.i.d. entries $X_{ij}$, $\mathbb{E}[X_{ij}] = 0$, $\mathbb{E}[X_{ij}^2] = 1$. Let $S_n = \frac{1}{n} X^T X$, $\gamma = p/n \to y \in (0, \infty)$. Then the empirical spectral distribution of $S_n$ converges weakly almost surely to the **Marchenko-Pastur distribution** $\mu_{MP}$, whose density is

    $$
    \rho_{MP}(x) = \frac{1}{2\pi xy} \sqrt{(\lambda_+ - x)(x - \lambda_-)} \cdot \mathbf{1}_{[\lambda_-, \lambda_+]}(x),
    $$

    where $\lambda_\pm = (1 \pm \sqrt{y})^2$. When $y > 1$, $\mu_{MP}$ has a point mass of $1 - 1/y$ at $x = 0$.

??? proof "Proof"
    **Stieltjes transform method.** Let $s_n(z)$ be the Stieltjes transform of the ESD of $S_n$. Using the identity

    $$
    s_n(z) = \frac{1}{p} \operatorname{tr}\!\left( S_n - zI \right)^{-1},
    $$

    and the matrix identity

    $$
    \frac{1}{n} X^T X - zI_p = -z \left( I_p - \frac{1}{nz} X^T X \right),
    $$

    combined with the Sherman-Morrison-Woodbury formula and column-by-column deletion technique, one can show that $s_n(z)$ converges to $s(z)$ satisfying the following equation:

    $$
    s(z) = \frac{1}{-z + y/(1 + s(z) \cdot z)} \cdot \frac{1}{1},
    $$

    which simplifies to the implicit equation

    $$
    y z s^2(z) + (z - 1 + y) s(z) + 1 = 0.
    $$

    Solving this equation and verifying the Stieltjes inversion formula $\rho(x) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \operatorname{Im} s(x + i\eta)$ recovers the density $\rho_{MP}$. $\blacksquare$

!!! theorem "Theorem 23.5 (Marchenko-Pastur law for general populations)"
    Let $X$ be an $n \times p$ matrix with i.i.d. entries, $\mathbb{E}[X_{ij}] = 0$, $\mathbb{E}[X_{ij}^2] = 1$. Let the ESD of the population covariance $\Sigma$ converge to $H$. Let $S_n = \frac{1}{n} X \Sigma X^T$, $p/n \to y$. Then the Stieltjes transform $s(z)$ of the limiting spectral distribution $F$ of $S_n$ satisfies

    $$
    s(z) = \int \frac{1}{\tau(1 - y - y z s(z)) - z} \, dH(\tau).
    $$

??? proof "Proof"
    The proof strategy is similar to Theorem 23.4, but the structure of $\Sigma$ must be taken into account during column deletion. Using $S_n = \frac{1}{n} X \Sigma X^T$ and the rank-one perturbation formula for resolvents, together with concentration inequalities and truncation arguments, one can show that the approximate equation satisfied by $s_n(z)$ converges to the above deterministic equation as $n \to \infty$. For a detailed proof, see Silverstein-Bai (1995). $\blacksquare$

!!! example "Example 23.4"
    **The Marchenko-Pastur distribution when $y = 1$.**

    When $p = n$ (i.e., $y = 1$), $\lambda_- = 0$, $\lambda_+ = 4$, and the density is

    $$
    \rho_{MP}(x) = \frac{1}{2\pi x} \sqrt{(4 - x) \cdot x} = \frac{1}{2\pi} \sqrt{\frac{4 - x}{x}}, \quad x \in (0, 4].
    $$

    Note that $\rho_{MP}(x) \to \infty$ as $x \to 0^+$, i.e., the eigenvalue density diverges near zero.

!!! example "Example 23.5"
    **Comparing MP distributions for different values of $y$.**

    - $y = 0.2$: $\lambda_- = (1 - \sqrt{0.2})^2 \approx 0.106$, $\lambda_+ = (1 + \sqrt{0.2})^2 \approx 2.294$, distribution concentrated near $1$.
    - $y = 1$: $\lambda_- = 0$, $\lambda_+ = 4$, distribution extends to $[0, 4]$.
    - $y = 5$: $\lambda_- = 0$, $\lambda_+ = (1 + \sqrt{5})^2 \approx 10.47$, point mass $1 - 1/5 = 0.8$ at $0$.

    As $y$ increases (sample size decreases relative to dimension), the spectral distribution spreads wider, reflecting the amplification of high-dimensional noise.

---

## 23.4 Empirical Spectral Distribution and the Stieltjes Transform Method

<div class="context-flow" markdown>

**Methodology**: Method of moments is suitable for proving existence → The **Stieltjes transform** $s(z) = \int \frac{d\mu(x)}{x-z}$ is the computational workhorse → Inversion formula $\rho(x) = \frac{1}{\pi}\text{Im}\,s(x+i0^+)$ recovers the density from the transform

</div>

The Stieltjes transform is the core analytical tool for studying limiting spectral distributions of random matrices.

!!! definition "Definition 23.9 (Stieltjes transform)"
    Let $\mu$ be a probability measure on $\mathbb{R}$. The **Stieltjes transform** of $\mu$ is defined as

    $$
    s_\mu(z) = \int_{\mathbb{R}} \frac{1}{x - z} \, d\mu(x), \quad z \in \mathbb{C}^+ = \{z : \operatorname{Im}(z) > 0\}.
    $$

    $s_\mu$ is a holomorphic function on $\mathbb{C}^+$ with $\operatorname{Im}(s_\mu(z)) > 0$.

!!! theorem "Theorem 23.6 (Stieltjes inversion formula)"
    Let $\mu$ be a probability measure and $s(z)$ its Stieltjes transform. If $\mu$ has a continuous density $\rho$ on $(a, b)$, then

    $$
    \rho(x) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \operatorname{Im}\, s(x + i\eta), \quad x \in (a, b).
    $$

    More generally, for continuity points $a < b$ of $\mu$,

    $$
    \mu\!\left((a, b)\right) = \frac{1}{\pi} \lim_{\eta \downarrow 0} \int_a^b \operatorname{Im}\, s(x + i\eta) \, dx.
    $$

??? proof "Proof"
    From $s(x + i\eta) = \int \frac{1}{t - x - i\eta} \, d\mu(t)$, taking the imaginary part gives

    $$
    \operatorname{Im}\, s(x + i\eta) = \int \frac{\eta}{(t - x)^2 + \eta^2} \, d\mu(t).
    $$

    Note that $\frac{\eta}{\pi((t-x)^2 + \eta^2)}$ is the Cauchy (Poisson) kernel centered at $x$ with half-width $\eta$, which tends to $\delta(t - x)$ as $\eta \to 0$. Therefore

    $$
    \frac{1}{\pi} \operatorname{Im}\, s(x + i\eta) = \int \frac{1}{\pi} \frac{\eta}{(t - x)^2 + \eta^2} \, d\mu(t) \to \rho(x),
    $$

    where convergence holds at continuity points of $\rho$. $\blacksquare$

!!! theorem "Theorem 23.7 (Continuity theorem for Stieltjes transforms)"
    Let $\{\mu_n\}$ be a sequence of probability measures with corresponding Stieltjes transforms $s_n(z)$. If for all $z \in \mathbb{C}^+$, $s_n(z) \to s(z)$, and $s(z)$ is the Stieltjes transform of some probability measure $\mu$, then $\mu_n \xrightarrow{w} \mu$ (weak convergence).

??? proof "Proof"
    There is a one-to-one correspondence between Stieltjes transforms and distribution functions (under appropriate conditions). Pointwise convergence $s_n(z) \to s(z)$ implies convergence of moments (via Laurent expansion), and uniqueness of the moment problem then yields weak convergence. A rigorous proof requires compactness arguments (Helly's selection theorem) and the property that the Stieltjes transform uniquely determines the measure. $\blacksquare$

!!! example "Example 23.6"
    **Verifying the semicircle law via Stieltjes transform.**

    For a Wigner matrix $W_n$, the normalized trace of its resolvent $s_n(z) = \frac{1}{n}\operatorname{tr}(W_n - zI)^{-1}$ satisfies the approximate equation

    $$
    s_n(z) \approx \frac{1}{-z - s_n(z)},
    $$

    and as $n \to \infty$, $s(z)$ satisfies $s = \frac{1}{-z - s}$, i.e., $s^2 + zs + 1 = 0$, with solution $s(z) = \frac{-z + \sqrt{z^2 - 4}}{2}$, exactly the Stieltjes transform of the semicircle distribution.

---

## 23.5 Eigenvalue Spacing and Repulsion

<div class="context-flow" markdown>

**Microscopic behavior**: Semicircle law / MP law = macroscopic (density) → Spacing statistics = microscopic · **Repulsion**: $p(s) \sim s^\beta$ ($s \to 0$) vs Poisson $p(s) = e^{-s}$ → Eigenvalues "push each other apart"; this is the essential difference between random matrices and independent random variables

</div>

A remarkable feature of random matrices is the repulsion effect between eigenvalues: eigenvalues tend to stay away from each other, and their spacing statistics differ fundamentally from those of independent random variables.

!!! definition "Definition 23.10 (Normalized spacing)"
    Let $\lambda_1 \le \lambda_2 \le \cdots \le \lambda_n$ be the ordered eigenvalues of a random matrix. At a bulk spectral point $E$, the local eigenvalue density is $\rho(E)$. The **normalized spacing** is defined as

    $$
    s_i = n \rho(E) (\lambda_{i+1} - \lambda_i),
    $$

    so that the mean of the normalized spacing is $1$.

!!! theorem "Theorem 23.8 (Wigner eigenvalue spacing statistics)"
    For GUE ($\beta = 2$) matrices, in the bulk of the spectrum, the distribution of normalized spacings approaches the **Gaudin distribution**. For small spacings $s \to 0$, the spacing probability density satisfies

    $$
    p(s) \sim c_\beta \, s^\beta, \quad s \to 0,
    $$

    where $\beta$ is the Dyson index (GOE: $\beta=1$, GUE: $\beta=2$, GSE: $\beta=4$). This means small spacings occur with vanishingly small probability — eigenvalues repel each other.

??? proof "Proof"
    We illustrate with GUE. Using the determinantal point process structure: GUE eigenvalues form a determinantal point process with the sine kernel

    $$
    K(x, y) = \frac{\sin(\pi(x - y))}{\pi(x - y)}
    $$

    as the correlation kernel. The gap probability is

    $$
    P(s > t) = \det(I - K_t),
    $$

    where $K_t$ is the restriction of the sine kernel operator to the interval $[0, t]$. By Fredholm determinant expansion, as $t \to 0$,

    $$
    P(s \le t) = 1 - \det(I - K_t) \sim \frac{(\pi t)^2}{2} + O(t^4),
    $$

    hence $p(t) \sim \pi^2 t$, i.e., $p(s) \propto s^2$ ($\beta = 2$). $\blacksquare$

!!! theorem "Theorem 23.9 (Wigner-Dyson-Mehta spacing distribution approximation)"
    In practice, the following approximate formulas (Wigner surmise) are commonly used for spacing distributions:

    - GOE ($\beta = 1$): $p(s) = \frac{\pi}{2} s \, e^{-\pi s^2 / 4}$;
    - GUE ($\beta = 2$): $p(s) = \frac{32}{\pi^2} s^2 \, e^{-4s^2 / \pi}$;
    - GSE ($\beta = 4$): $p(s) = \frac{2^{18}}{3^6 \pi^3} s^4 \, e^{-64 s^2 / (9\pi)}$.

    These formulas exactly describe the spacing distribution of $2 \times 2$ matrices and are excellent approximations for large matrices.

??? proof "Proof"
    We illustrate with the $2 \times 2$ GOE case. Let $M = \begin{pmatrix} a & b \\ b & c \end{pmatrix}$, $a, c \sim N(0,1)$, $b \sim N(0, 1/2)$. The eigenvalue gap is $s = \lambda_+ - \lambda_- = \sqrt{(a-c)^2 + 4b^2}$. Let $u = a - c$, $v = 2b$, then $u \sim N(0, 2)$, $v \sim N(0, 2)$, $s = \sqrt{u^2 + v^2}$. Converting to polar coordinates: $p(s) = \frac{s}{2} e^{-s^2/4}$. After appropriate normalization to make $\langle s \rangle = 1$, we obtain the Wigner surmise $p(s) = \frac{\pi}{2} s \, e^{-\pi s^2/4}$. $\blacksquare$

!!! example "Example 23.7"
    **Comparison of Poisson and GUE spacing statistics.**

    - Independent random eigenvalues (e.g., diagonal random matrices) have exponentially distributed spacings $p(s) = e^{-s}$ (Poisson statistics), $p(0) = 1$.
    - GUE spacings are $p(s) \approx \frac{32}{\pi^2} s^2 e^{-4s^2/\pi}$, $p(0) = 0$.

    GUE has zero density at $s = 0$, reflecting eigenvalue repulsion; Poisson statistics have maximum density at zero spacing, indicating that independent eigenvalues can come arbitrarily close.

---

## 23.6 Tracy-Widom Distribution

<div class="context-flow" markdown>

**Fine scale**: Semicircle law (macroscopic) → Spacing (microscopic) → **Largest eigenvalue fluctuations** (edge) · $\lambda_{\max} \approx 2 + n^{-2/3} \cdot F_\beta$ → Airy kernel + Painleve II equation
**Universality**: $F_\beta$ does not depend on the specific distribution of matrix entries, only on the symmetry class $\beta$ → Links to S23.7 statistical tests

</div>

The semicircle law describes the bulk distribution of eigenvalues, while the Tracy-Widom distribution characterizes the fluctuations of the largest eigenvalue.

!!! definition "Definition 23.11 (Tracy-Widom distribution)"
    The **Tracy-Widom distribution** $F_\beta$ ($\beta = 1, 2, 4$) describes the limiting distribution of the largest eigenvalue of a random matrix after appropriate centering and scaling. For $\beta = 2$ (GUE), the distribution function is

    $$
    F_2(s) = \exp\!\left( -\int_s^{\infty} (x - s) \, q(x)^2 \, dx \right),
    $$

    where $q(x)$ is the unique solution of the Painleve II equation $q''(x) = xq(x) + 2q(x)^3$ satisfying the Airy decay condition $q(x) \sim \operatorname{Ai}(x)$ ($x \to +\infty$).

!!! theorem "Theorem 23.10 (Tracy-Widom limit for GUE largest eigenvalue)"
    Let $M_n$ be an $n \times n$ GUE matrix and $\lambda_{\max}$ its largest eigenvalue. Then

    $$
    \frac{\lambda_{\max} - 2}{n^{-2/3}} \xrightarrow{d} F_2,
    $$

    i.e., $P\!\left( n^{2/3}(\lambda_{\max} - 2) \le s \right) \to F_2(s)$, where $F_2$ is the Tracy-Widom distribution.

??? proof "Proof"
    **Proof sketch.** GUE eigenvalues form a determinantal point process with correlation kernel

    $$
    K_n(x, y) = \sum_{k=0}^{n-1} \varphi_k(x) \varphi_k(y),
    $$

    where $\varphi_k$ are normalized Hermite functions. The distribution of the largest eigenvalue is

    $$
    P(\lambda_{\max} \le t) = \det(I - K_n)|_{L^2(t, \infty)}.
    $$

    At the spectral edge $t = 2 + s n^{-2/3}$, using the Plancherel-Rotach asymptotic formula, $K_n$ under appropriate scaling converges to the **Airy kernel**

    $$
    K_{\text{Airy}}(x, y) = \frac{\operatorname{Ai}(x)\operatorname{Ai}'(y) - \operatorname{Ai}'(x)\operatorname{Ai}(y)}{x - y}.
    $$

    Therefore $P(\lambda_{\max} \le 2 + sn^{-2/3}) \to \det(I - K_{\text{Airy}})|_{L^2(s, \infty)} = F_2(s)$. $\blacksquare$

!!! theorem "Theorem 23.11 (Tracy-Widom limit for GOE largest eigenvalue)"
    Let $M_n$ be an $n \times n$ GOE matrix. Then

    $$
    \frac{\lambda_{\max} - 2}{n^{-2/3}} \xrightarrow{d} F_1,
    $$

    where $F_1(s) = \exp\!\left( -\frac{1}{2}\int_s^\infty q(x) \, dx \right) \cdot F_2(s)^{1/2}$, and $q$ is the Painleve II solution above.

??? proof "Proof"
    GOE eigenvalues form a Pfaffian point process (rather than a determinantal one). After edge scaling, the correlation kernel converges to a symmetrized version of the Airy kernel. The distribution of the largest eigenvalue can be written as a Fredholm Pfaffian, whose limit gives $F_1$. $\blacksquare$

!!! example "Example 23.8"
    **Numerical characteristics of the Tracy-Widom distribution.**

    Numerical characteristics of $F_2$:
    - Mean $\approx -1.771$
    - Standard deviation $\approx 0.813$
    - Skewness $\approx 0.224$
    - Kurtosis $\approx 0.093$

    The distribution is left-skewed, meaning the largest eigenvalue tends to be slightly below its mean $2$. The $F_1$ distribution has larger fluctuations (standard deviation $\approx 1.268$) because real symmetric matrices have fewer degrees of freedom.

---

## 23.7 Applications of Random Matrices in Statistics

<div class="context-flow" markdown>

**Applications**: **BBP phase transition** — when signal strength $\theta > \sqrt{p/n}$, the largest eigenvalue "pops out" from the MP edge → Theoretical threshold for signal detection · Tracy-Widom distribution replaces classical $F$/$\chi^2$ distributions for high-dimensional hypothesis testing

</div>

Random matrix theory provides both theoretical foundations and practical tools for high-dimensional statistics.

!!! definition "Definition 23.12 (High-dimensional asymptotic framework)"
    In the **high-dimensional asymptotic framework**, the data dimension $p$ and sample size $n$ tend to infinity simultaneously with ratio $p/n \to y \in (0, \infty)$. This is fundamentally different from the classical framework where $p$ is fixed and $n \to \infty$.

<div class="context-flow" markdown>

**Insight**: The BBP phase transition reveals a profound "information-theoretic limit" — when signal strength $\theta$ is below $\sqrt{y}$, **no method** can detect the signal from the eigenvalues of the sample covariance

</div>

!!! theorem "Theorem 23.12 (Baik-Ben Arous-Peche phase transition)"
    **BBP phase transition**: Let the population covariance matrix be $\Sigma = I + \theta \mathbf{v}\mathbf{v}^T$ (rank-one perturbation), $p/n \to y$. Let $\ell_1$ be the largest eigenvalue of the sample covariance matrix $S_n$. Then:

    - If $\theta < \sqrt{y}$, then $\ell_1 \to (1 + \sqrt{y})^2$ (same as when $\theta = 0$);
    - If $\theta > \sqrt{y}$, then $\ell_1 \to (1 + \theta)(1 + y/\theta) > (1 + \sqrt{y})^2$.

    The critical value $\theta_c = \sqrt{y}$ marks the phase transition for signal detectability.

??? proof "Proof"
    **Proof sketch.** Using the Stieltjes transform method. Under $\Sigma = I + \theta \mathbf{v}\mathbf{v}^T$, the limiting spectral distribution of $S_n$ is still described by the Marchenko-Pastur law (rank-one perturbation does not affect the limiting spectral distribution), but the behavior of the largest eigenvalue depends on the magnitude of $\theta$.

    The key step is to analyze the behavior of the resolvent $(S_n - zI)^{-1}$ outside the spectral edge. When $z$ is outside $\lambda_+ = (1+\sqrt{y})^2$, $\frac{1}{n}\operatorname{tr}(S_n - zI)^{-1} \to s(z)$. Using the matrix perturbation formula, the largest sample eigenvalue satisfies

    $$
    1 + \theta \cdot \frac{1}{p}\operatorname{tr}\!\left(\frac{S_n^{(0)}}{S_n^{(0)} - \ell_1 I}\right) \approx 0,
    $$

    where $S_n^{(0)}$ is the matrix with the signal removed. When $\theta > \sqrt{y}$, this equation has a solution outside $(1+\sqrt{y})^2$, meaning the largest eigenvalue "pops out" from the Marchenko-Pastur support. $\blacksquare$

!!! example "Example 23.9"
    **Signal detection problem.**

    In wireless communications, the received signal model is $\mathbf{x} = \sqrt{\theta} \, s \, \mathbf{a} + \mathbf{n}$, where $s$ is the signal, $\mathbf{a}$ is the steering vector, and $\mathbf{n} \sim N(\mathbf{0}, I)$. After $n$ observations, the sample covariance matrix is $S_n = \frac{1}{n}XX^T$. By the BBP phase transition, when the signal-to-noise ratio $\theta > \sqrt{p/n}$, the largest eigenvalue of $S_n$ will significantly deviate from the upper edge of the Marchenko-Pastur distribution, enabling signal detection.

!!! example "Example 23.10"
    **Correction of Roy's largest root test.**

    The classical Roy's largest root test statistic is $\ell_1(S_1 S_2^{-1})$. Under the high-dimensional framework $p, n \to \infty$, $p/n \to y$, the limiting distribution of this statistic under the null hypothesis is no longer the classical Roy distribution but the Tracy-Widom distribution $F_1$. Therefore, the rejection region should be based on Tracy-Widom quantiles rather than traditional tables.

---

## 23.8 Introduction to Free Probability

<div class="context-flow" markdown>

**Algebraization**: Independent random matrices $A_n, B_n$ → Asymptotically **free** (Voiculescu) → The limiting spectrum of $A+B$ is computed via **$R$-transform** additivity → Free CLT: the limit is the semicircle distribution (compare classical CLT → normal distribution)
**Unified perspective**: Semicircle law = central limit theorem of free probability

</div>

Free probability theory studies noncommutative random variables and provides an algebraic framework for the asymptotic spectral behavior of random matrices.

!!! definition "Definition 23.13 (Noncommutative probability space)"
    A **noncommutative probability space** is a pair $(\mathcal{A}, \varphi)$, where $\mathcal{A}$ is a unital algebra (not necessarily commutative) and $\varphi: \mathcal{A} \to \mathbb{C}$ is a linear functional satisfying $\varphi(1) = 1$ (called a tracial state). Elements of $\mathcal{A}$ are called noncommutative random variables.

!!! definition "Definition 23.14 (Free independence)"
    In a noncommutative probability space $(\mathcal{A}, \varphi)$, subalgebras $\mathcal{A}_1, \ldots, \mathcal{A}_k$ are called **freely independent** if for any $a_j \in \mathcal{A}_{i_j}$ ($j = 1, \ldots, m$) satisfying $\varphi(a_j) = 0$ and with adjacent elements from different subalgebras ($i_1 \ne i_2 \ne \cdots \ne i_m$),

    $$
    \varphi(a_1 a_2 \cdots a_m) = 0.
    $$

    Free independence is the analogue of classical independence in the noncommutative setting, but the two are fundamentally different.

!!! theorem "Theorem 23.13 (Voiculescu's asymptotic freeness theorem)"
    Let $A_n, B_n$ be $n \times n$ independent random matrices, $A_n$ a Wigner matrix, $B_n = U_n D_n U_n^*$ where $D_n$ is a deterministic diagonal matrix (whose ESD converges to $\nu$) and $U_n$ is a Haar unitary matrix independent of $A_n$. Then as $n \to \infty$, $A_n$ and $B_n$ are **asymptotically free**, meaning their mixed moments with respect to the normalized trace $\varphi(\cdot) = \frac{1}{n}\operatorname{tr}(\cdot)$ satisfy the algebraic relations of free independence.

??? proof "Proof"
    **Proof sketch.** One needs to verify that for centered alternating products, the normalized trace tends to zero. That is, for $p(A_n) = A_n^k - \varphi(A_n^k)I$ and $q(B_n) = B_n^l - \varphi(B_n^l)I$, one must show

    $$
    \frac{1}{n}\operatorname{tr}\!\left(p_1(A_n) q_1(B_n) p_2(A_n) q_2(B_n) \cdots \right) \to 0.
    $$

    The key tool is the Weingarten integration formula, which gives the integral of products of Haar unitary matrix entries over $U_n$. Using this formula, the trace above can be expanded as a sum over permutations, and cancellations among the leading terms (guaranteed by centering) ensure that the whole expression tends to zero. $\blacksquare$

!!! definition "Definition 23.15 (Free convolution)"
    Let $\mu, \nu$ be probability measures on $\mathbb{R}$. If $a, b$ are self-adjoint elements in a noncommutative probability space with distributions $\mu, \nu$ respectively, and $a, b$ are freely independent, then the distribution of $a + b$ is called the **free (additive) convolution** of $\mu$ and $\nu$, denoted $\mu \boxplus \nu$.

!!! theorem "Theorem 23.14 ($R$-transform of free convolution)"
    Let $\mu, \nu$ be probability measures, $G_\mu(z) = \int \frac{d\mu(x)}{z - x}$ the Cauchy transform (note the sign difference from the Stieltjes transform), $K_\mu$ its functional inverse ($G_\mu(K_\mu(w)) = w$), and $R_\mu(w) = K_\mu(w) - 1/w$. Then

    $$
    R_{\mu \boxplus \nu}(w) = R_\mu(w) + R_\nu(w).
    $$

    That is, the $R$-transform is additive under free convolution, analogous to the log-additivity of characteristic functions under classical independence.

??? proof "Proof"
    Using combinatorial free probability theory (free cumulant theory). Define free cumulants $\kappa_n(\mu)$ via the moment-cumulant formula

    $$
    m_n = \sum_{\pi \in NC(n)} \prod_{V \in \pi} \kappa_{|V|},
    $$

    where the sum is over all non-crossing partitions $NC(n)$ of $\{1, \ldots, n\}$. One can show that free independence is equivalent to vanishing of mixed free cumulants. Therefore $\kappa_n(\mu \boxplus \nu) = \kappa_n(\mu) + \kappa_n(\nu)$. The Laurent expansion coefficients of the $R$-transform are exactly the free cumulants: $R_\mu(w) = \sum_{n=1}^\infty \kappa_n(\mu) w^{n-1}$. $\blacksquare$

!!! theorem "Theorem 23.15 (Free central limit theorem for the semicircle law)"
    Let $a_1, a_2, \ldots$ be freely independent identically distributed self-adjoint noncommutative random variables with $\varphi(a_i) = 0$, $\varphi(a_i^2) = 1$. Then

    $$
    \frac{a_1 + a_2 + \cdots + a_n}{\sqrt{n}} \xrightarrow{d} \mu_{sc},
    $$

    i.e., the normalized partial sum converges in distribution to the semicircle distribution. This is the free probability analogue of the classical central limit theorem (where the limit is the normal distribution).

??? proof "Proof"
    By the additivity of free cumulants, $S_n = \frac{1}{\sqrt{n}}(a_1 + \cdots + a_n)$ has free cumulants $\kappa_k(S_n) = n^{1 - k/2} \kappa_k(a_1)$. As $n \to \infty$, $\kappa_1(S_n) = 0$, $\kappa_2(S_n) = 1$, $\kappa_k(S_n) \to 0$ ($k \ge 3$). The semicircle distribution has free cumulants $\kappa_2 = 1$, $\kappa_k = 0$ ($k \ne 2$), because its $R$-transform is $R(w) = w$. $\blacksquare$

!!! example "Example 23.11"
    **Free convolution of two semicircle distributions.**

    Let $\mu = \mu_{sc}(0, 1)$ (standard semicircle, $R_\mu(w) = w$), $\nu = \mu_{sc}(0, 1)$. Then

    $$
    R_{\mu \boxplus \nu}(w) = 2w,
    $$

    corresponding to $\mu_{sc}(0, \sqrt{2})$, the semicircle distribution with radius $2\sqrt{2}$ and variance $2$. More generally, the free convolution of semicircle distributions with variances $\sigma_1^2$ and $\sigma_2^2$ is again a semicircle distribution with variance $\sigma_1^2 + \sigma_2^2$.

!!! example "Example 23.12"
    **Computing the limiting spectrum of $A + B$ using free probability.**

    Let $A_n$ be a Wigner matrix (limiting spectrum: semicircle $\mu_{sc}$), and $B_n$ be an independent deterministic matrix conjugated by a Haar unitary, whose ESD converges to the Bernoulli distribution $\nu = \frac{1}{2}\delta_{-1} + \frac{1}{2}\delta_1$. By asymptotic freeness, the limiting spectral distribution of $A_n + B_n$ is $\mu_{sc} \boxplus \nu$, which can be computed via the $R$-transform method. The $R$-transform of $\nu$ is $R_\nu(w) = \frac{w}{1 - w^2}$ (from the moment-cumulant relation), so $R_{\mu_{sc} \boxplus \nu}(w) = w + \frac{w}{1 - w^2}$, and the limiting density can be obtained numerically via the inverse function relation.

---

## Chapter Summary

This chapter introduced the basic framework and core results of random matrix theory:

1. **Gaussian ensembles** (GOE, GUE, GSE) serve as classical models of random matrices, with joint eigenvalue densities possessing elegant determinantal/Pfaffian structure.
2. The **Wigner semicircle law** describes the macroscopic limit of the empirical spectral distribution of Wigner matrices.
3. The **Marchenko-Pastur law** characterizes the spectral behavior of high-dimensional sample covariance matrices and is the theoretical cornerstone of high-dimensional statistics.
4. The **Stieltjes transform** is the core analytical tool for studying limiting spectral distributions.
5. Eigenvalue **repulsion** and Wigner-Dyson spacing statistics reveal the essential difference between random matrices and independent random variables.
6. The **Tracy-Widom distribution** describes fine fluctuations of the largest eigenvalue.
7. The **BBP phase transition** provides a theoretical threshold for high-dimensional signal detection.
8. **Free probability theory** provides algebraic tools for asymptotic spectral computations of random matrices.
