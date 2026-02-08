# Chapter 15  Norms and Perturbation Theory

<div class="context-flow" markdown>

**Prerequisites**: Ch14 spectral radius and norms · **Chapter arc**: Vector norms → Matrix norms (operator/Frobenius) → **Condition number** quantifies ill-conditioning → Bauer--Fike/Weyl quantify eigenvalue sensitivity

</div>

Norms provide a measure of "size" for vector spaces and matrix spaces, and are the cornerstone of matrix analysis and numerical linear algebra. Perturbation theory studies how eigenvalues, singular values, and solutions of linear systems change when a matrix undergoes small perturbations. These theories are essential for understanding the stability and accuracy of numerical computations. This chapter systematically develops the theory of vector norms and matrix norms, introduces the concept of condition number, and discusses perturbation bounds for eigenvalues and singular values in depth.

---

## 15.1 Vector norms

<div class="context-flow" markdown>

**Chapter arc**: $\ell_p$ norm family ($p=1,2,\infty$) → Holder inequality ($p,q$ conjugate) → finite-dimensional norm equivalence (unique topology)

</div>

!!! definition "Definition 15.1 (Vector norm)"
    A **vector norm** on $\mathbb{C}^n$ is a function $\|\cdot\| : \mathbb{C}^n \to \mathbb{R}$ satisfying for all $\mathbf{x}, \mathbf{y} \in \mathbb{C}^n$ and $\alpha \in \mathbb{C}$:

    (1) **Non-negativity**: $\|\mathbf{x}\| \geq 0$, with equality if and only if $\mathbf{x} = \mathbf{0}$;

    (2) **Homogeneity**: $\|\alpha\mathbf{x}\| = |\alpha|\|\mathbf{x}\|$;

    (3) **Triangle inequality**: $\|\mathbf{x} + \mathbf{y}\| \leq \|\mathbf{x}\| + \|\mathbf{y}\|$.

!!! definition "Definition 15.2 ($\ell_p$ norm)"
    For $1 \leq p \leq \infty$, the **$\ell_p$ norm** on $\mathbb{C}^n$ is defined as

    $$\|\mathbf{x}\|_p = \left(\sum_{i=1}^{n} |x_i|^p\right)^{1/p}, \quad 1 \leq p < \infty,$$

    $$\|\mathbf{x}\|_\infty = \max_{1 \leq i \leq n} |x_i|.$$

    Of particular importance are:

    - $\|\mathbf{x}\|_1 = \sum_{i=1}^{n}|x_i|$ (Manhattan norm);
    - $\|\mathbf{x}\|_2 = \sqrt{\sum_{i=1}^{n}|x_i|^2} = \sqrt{\mathbf{x}^*\mathbf{x}}$ (Euclidean norm);
    - $\|\mathbf{x}\|_\infty = \max_i |x_i|$ (Chebyshev norm).

!!! theorem "Theorem 15.1 (Holder inequality)"
    Let $p, q \geq 1$ satisfy $\frac{1}{p} + \frac{1}{q} = 1$ (conjugate exponents). Then for any $\mathbf{x}, \mathbf{y} \in \mathbb{C}^n$,

    $$|\mathbf{x}^*\mathbf{y}| \leq \|\mathbf{x}\|_p\|\mathbf{y}\|_q.$$

    In particular, $p = q = 2$ gives the Cauchy--Schwarz inequality.

??? proof "Proof"
    The case $\|\mathbf{x}\|_p = 0$ or $\|\mathbf{y}\|_q = 0$ is trivial. Otherwise assume $\|\mathbf{x}\|_p = \|\mathbf{y}\|_q = 1$. By Young's inequality $ab \leq \frac{a^p}{p} + \frac{b^q}{q}$ ($a, b \geq 0$), setting $a = |x_i|$, $b = |y_i|$ and summing gives

    $$\sum_{i=1}^{n}|x_i||y_i| \leq \frac{1}{p}\sum|x_i|^p + \frac{1}{q}\sum|y_i|^q = \frac{1}{p} + \frac{1}{q} = 1.$$

    Therefore $|\mathbf{x}^*\mathbf{y}| \leq \sum|x_i||y_i| \leq 1 = \|\mathbf{x}\|_p\|\mathbf{y}\|_q$. The general case is reduced to the above by dividing by $\|\mathbf{x}\|_p\|\mathbf{y}\|_q$. $\blacksquare$

!!! theorem "Theorem 15.2 (Norm equivalence theorem)"
    Any two norms $\|\cdot\|_\alpha$ and $\|\cdot\|_\beta$ on $\mathbb{C}^n$ are **equivalent**, i.e., there exist positive constants $c_1, c_2 > 0$ such that for all $\mathbf{x} \in \mathbb{C}^n$,

    $$c_1\|\mathbf{x}\|_\alpha \leq \|\mathbf{x}\|_\beta \leq c_2\|\mathbf{x}\|_\alpha.$$

??? proof "Proof"
    It suffices to show that any norm $\|\cdot\|$ is equivalent to $\|\cdot\|_2$. Let $\mathbf{e}_1, \ldots, \mathbf{e}_n$ be the standard basis. For $\mathbf{x} = \sum x_i\mathbf{e}_i$,

    $$\|\mathbf{x}\| \leq \sum |x_i|\|\mathbf{e}_i\| \leq \left(\max_i\|\mathbf{e}_i\|\right)\sum|x_i| \leq \left(\max_i\|\mathbf{e}_i\|\right)\sqrt{n}\|\mathbf{x}\|_2.$$

    Set $c_2 = \sqrt{n}\max_i\|\mathbf{e}_i\|$; then $\|\mathbf{x}\| \leq c_2\|\mathbf{x}\|_2$. This shows $\|\cdot\|$ is Lipschitz continuous with respect to $\|\cdot\|_2$.

    Consider the compact set $S = \{\mathbf{x} : \|\mathbf{x}\|_2 = 1\}$. The continuous function $\|\cdot\|$ attains its minimum $c_1 > 0$ on $S$ (since $\|\mathbf{x}\| > 0$ for $\mathbf{x} \in S$). Therefore for $\mathbf{x} \neq \mathbf{0}$, $\|\mathbf{x}/\|\mathbf{x}\|_2\| \geq c_1$, i.e., $\|\mathbf{x}\| \geq c_1\|\mathbf{x}\|_2$. $\blacksquare$

!!! proposition "Proposition 15.1 (Relations between $\ell_p$ norms)"
    For any $\mathbf{x} \in \mathbb{C}^n$ and $1 \leq p \leq q \leq \infty$,

    $$\|\mathbf{x}\|_q \leq \|\mathbf{x}\|_p \leq n^{1/p - 1/q}\|\mathbf{x}\|_q.$$

    In particular:

    - $\|\mathbf{x}\|_2 \leq \|\mathbf{x}\|_1 \leq \sqrt{n}\|\mathbf{x}\|_2$;
    - $\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_2 \leq \sqrt{n}\|\mathbf{x}\|_\infty$;
    - $\|\mathbf{x}\|_\infty \leq \|\mathbf{x}\|_1 \leq n\|\mathbf{x}\|_\infty$.

!!! example "Example 15.1"
    Let $\mathbf{x} = (3, -4, 0, 1)^T$. Compute the $\ell_p$ norms.

    - $\|\mathbf{x}\|_1 = |3| + |-4| + |0| + |1| = 8$;
    - $\|\mathbf{x}\|_2 = \sqrt{9 + 16 + 0 + 1} = \sqrt{26} \approx 5.10$;
    - $\|\mathbf{x}\|_\infty = \max\{3, 4, 0, 1\} = 4$.

    Verification: $\|\mathbf{x}\|_\infty = 4 \leq 5.10 = \|\mathbf{x}\|_2 \leq 8 = \|\mathbf{x}\|_1$.

    $\|\mathbf{x}\|_1 = 8 \leq \sqrt{4}\cdot 5.10 = 10.2$? Indeed $\sqrt{n}\|\mathbf{x}\|_2 = 2 \times 5.10 = 10.2$, and $\|\mathbf{x}\|_1 = 8 \leq 10.2$.

!!! example "Example 15.2"
    **Unit balls of $\ell_p$ norms.** The unit ball $B_p = \{\mathbf{x} \in \mathbb{R}^2 : \|\mathbf{x}\|_p \leq 1\}$ changes shape with $p$:

    - $p = 1$: diamond (square rotated $45°$), vertices at $(\pm 1, 0)$ and $(0, \pm 1)$;
    - $p = 2$: disc;
    - $p = \infty$: square $[-1, 1]^2$;
    - $p \to 0^+$: the unit ball degenerates to line segments on the coordinate axes.

    As $p$ increases, the unit ball expands from a diamond to a circle to a square, reflecting the relationship $\|\mathbf{x}\|_q \leq \|\mathbf{x}\|_p$ ($p \leq q$).

---

## 15.2 Matrix norms

<div class="context-flow" markdown>

**Chapter arc**: Matrix norm + **submultiplicativity** $\|AB\|\leq\|A\|\|B\|$ → Frobenius norm $=\sqrt{\sum\sigma_i^2}$ is a unitarily invariant "entrywise" norm

</div>

!!! definition "Definition 15.3 (Matrix norm)"
    A **matrix norm** on $\mathbb{C}^{m \times n}$ is a function $\|\cdot\| : \mathbb{C}^{m \times n} \to \mathbb{R}$ satisfying for all $A, B$ and $\alpha \in \mathbb{C}$:

    (1) **Non-negativity**: $\|A\| \geq 0$, with equality if and only if $A = O$;

    (2) **Homogeneity**: $\|\alpha A\| = |\alpha|\|A\|$;

    (3) **Triangle inequality**: $\|A + B\| \leq \|A\| + \|B\|$.

    If it additionally satisfies (4) **submultiplicativity**: $\|AB\| \leq \|A\|\|B\|$ (when the product is defined), it is called a **submultiplicative norm**.

!!! definition "Definition 15.4 (Frobenius norm)"
    Let $A = (a_{ij}) \in \mathbb{C}^{m \times n}$. The **Frobenius norm** is defined as

    $$\|A\|_F = \sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n}|a_{ij}|^2} = \sqrt{\operatorname{tr}(A^*A)} = \sqrt{\sum_{i=1}^{\min(m,n)}\sigma_i^2}.$$

!!! theorem "Theorem 15.3 (Properties of the Frobenius norm)"
    The Frobenius norm satisfies:

    (1) It is a matrix norm and is submultiplicative;

    (2) $\|A\|_F = \|A^*\|_F$;

    (3) For unitary (orthogonal) matrices $U, V$, $\|UAV\|_F = \|A\|_F$ (unitary invariance);

    (4) $\|A\|_F^2 = \sum_{i=1}^{r}\sigma_i^2$, where $\sigma_i$ are the singular values of $A$;

    (5) $\|A\|_2 \leq \|A\|_F \leq \sqrt{r}\|A\|_2$, where $r = \operatorname{rank}(A)$.

??? proof "Proof"
    **(1)** Submultiplicativity: $\|AB\|_F^2 = \sum_{i,k}\left|\sum_j a_{ij}b_{jk}\right|^2 \leq \sum_{i,k}\left(\sum_j|a_{ij}|^2\right)\left(\sum_j|b_{jk}|^2\right)$ (by Cauchy--Schwarz) $= \sum_i\left(\sum_j|a_{ij}|^2\right)\sum_k\left(\sum_j|b_{jk}|^2\right) = \|A\|_F^2\|B\|_F^2$.

    **(3)** $\|UAV\|_F^2 = \operatorname{tr}(V^*A^*U^*UAV) = \operatorname{tr}(V^*A^*AV) = \operatorname{tr}(A^*AVV^*) = \operatorname{tr}(A^*A) = \|A\|_F^2$.

    **(4)** By SVD, $A = U\Sigma V^*$, $\|A\|_F^2 = \|\Sigma\|_F^2 = \sum\sigma_i^2$.

    **(5)** $\|A\|_F^2 = \sum\sigma_i^2 \geq \sigma_1^2 = \|A\|_2^2$, and $\sum\sigma_i^2 \leq r\sigma_1^2$. $\blacksquare$

!!! example "Example 15.3"
    Let $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$.

    $\|A\|_F = \sqrt{1 + 4 + 9 + 16} = \sqrt{30} \approx 5.48$.

    $A^*A = \begin{pmatrix} 10 & 14 \\ 14 & 20 \end{pmatrix}$, $\operatorname{tr} = 30$, $\det = 200 - 196 = 4$. Eigenvalues $\mu = 15 \pm \sqrt{225 - 4} = 15 \pm \sqrt{221}$. $\sigma_1 = \sqrt{15 + \sqrt{221}} \approx \sqrt{29.87} \approx 5.47$, $\sigma_2 = \sqrt{15 - \sqrt{221}} \approx \sqrt{0.13} \approx 0.37$.

    Verification: $\|A\|_F = \sqrt{\sigma_1^2 + \sigma_2^2} = \sqrt{30} \approx 5.48$.

---

## 15.3 Operator norms (induced norms)

<div class="context-flow" markdown>

**Intuition**: $\|A\|_1$=max column sum, $\|A\|_\infty$=max row sum, $\|A\|_2=\sigma_1(A)$ · Operator norms are automatically submultiplicative and satisfy $\rho(A)\leq\|A\|$

</div>

!!! definition "Definition 15.5 (Operator norm)"
    Let $\|\cdot\|_\alpha$ and $\|\cdot\|_\beta$ be vector norms on $\mathbb{C}^n$ and $\mathbb{C}^m$ respectively. The **operator norm** **induced** by them is defined as

    $$\|A\|_{\alpha \to \beta} = \max_{\mathbf{x} \neq \mathbf{0}} \frac{\|A\mathbf{x}\|_\beta}{\|\mathbf{x}\|_\alpha} = \max_{\|\mathbf{x}\|_\alpha = 1} \|A\mathbf{x}\|_\beta.$$

    When $\alpha = \beta = p$, we write $\|A\|_p$.

!!! theorem "Theorem 15.4 (Computation formulas for operator norms)"
    Let $A = (a_{ij}) \in \mathbb{C}^{m \times n}$. Then:

    (1) $\|A\|_1 = \max_{1 \leq j \leq n} \sum_{i=1}^{m} |a_{ij}|$ (maximum column sum);

    (2) $\|A\|_\infty = \max_{1 \leq i \leq m} \sum_{j=1}^{n} |a_{ij}|$ (maximum row sum);

    (3) $\|A\|_2 = \sigma_1(A)$ (largest singular value).

??? proof "Proof"
    **(1)** Let $\mathbf{x}$ satisfy $\|\mathbf{x}\|_1 = 1$. Then

    $$\|A\mathbf{x}\|_1 = \sum_i\left|\sum_j a_{ij}x_j\right| \leq \sum_i\sum_j|a_{ij}||x_j| = \sum_j|x_j|\sum_i|a_{ij}| \leq \max_j\sum_i|a_{ij}|.$$

    Equality is attained at $\mathbf{x} = \mathbf{e}_{j^*}$ (the standard basis vector corresponding to the maximum column sum).

    **(2)** Similarly, $\|A\mathbf{x}\|_\infty = \max_i\left|\sum_j a_{ij}x_j\right| \leq \max_i\sum_j|a_{ij}|\|\mathbf{x}\|_\infty$. Equality is attained for a suitable choice of $\mathbf{x}$.

    **(3)** $\|A\mathbf{x}\|_2^2 = \mathbf{x}^*A^*A\mathbf{x}$. $\max_{\|\mathbf{x}\|_2=1}\mathbf{x}^*A^*A\mathbf{x} = \lambda_{\max}(A^*A) = \sigma_1^2(A)$ (by the Rayleigh quotient). Hence $\|A\|_2 = \sigma_1(A)$. $\blacksquare$

!!! theorem "Theorem 15.5 (Properties of operator norms)"
    Operator norms satisfy:

    (1) $\|I\| = 1$;

    (2) Submultiplicativity: $\|AB\| \leq \|A\|\|B\|$;

    (3) Consistency: $\|A\mathbf{x}\| \leq \|A\|\|\mathbf{x}\|$ (for the corresponding vector norm);

    (4) $\rho(A) \leq \|A\|$.

??? proof "Proof"
    **(2)** $\|AB\mathbf{x}\| \leq \|A\|\|B\mathbf{x}\| \leq \|A\|\|B\|\|\mathbf{x}\|$; taking the supremum over all $\|\mathbf{x}\| = 1$ gives $\|AB\| \leq \|A\|\|B\|$.

    **(3)** Follows directly from the definition.

    **(4)** Already proved in Theorem 14.6. $\blacksquare$

!!! example "Example 15.4"
    Let $A = \begin{pmatrix} 1 & -2 & 3 \\ 4 & 0 & -1 \end{pmatrix}$.

    - $\|A\|_1 = \max\{|1|+|4|,\; |-2|+|0|,\; |3|+|-1|\} = \max\{5, 2, 4\} = 5$;
    - $\|A\|_\infty = \max\{|1|+|-2|+|3|,\; |4|+|0|+|-1|\} = \max\{6, 5\} = 6$;
    - $\|A\|_2 = \sigma_1(A)$, which requires computing the largest eigenvalue of $A^TA$.

!!! example "Example 15.5"
    For normal matrices, $\|A\|_2 = \rho(A)$.

    Let $A = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}$ (rotation by $90°$). $A$ is a normal matrix (in fact orthogonal), with eigenvalues $\pm i$, $\rho(A) = 1$.

    $A^TA = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$, $\sigma_1 = 1$. Therefore $\|A\|_2 = 1 = \rho(A)$.

---

## 15.4 Relations between norms

<div class="context-flow" markdown>

**Key inequalities**: $\|A\|_2\leq\|A\|_F\leq\sqrt{r}\|A\|_2$ · $\|A\|_2\leq\sqrt{\|A\|_1\|A\|_\infty}$ connects three operator norms

</div>

!!! theorem "Theorem 15.6 (Inequalities between matrix norms)"
    Let $A \in \mathbb{C}^{m \times n}$. Then:

    (1) $\|A\|_2 \leq \|A\|_F \leq \sqrt{n}\|A\|_2$;

    (2) $\frac{1}{\sqrt{n}}\|A\|_\infty \leq \|A\|_2 \leq \sqrt{m}\|A\|_\infty$;

    (3) $\frac{1}{\sqrt{m}}\|A\|_1 \leq \|A\|_2 \leq \sqrt{n}\|A\|_1$;

    (4) $\|A\|_2 \leq \sqrt{\|A\|_1\|A\|_\infty}$;

    (5) $\frac{1}{\sqrt{mn}}\|A\|_F \leq \|A\|_\infty$, $\frac{1}{\sqrt{mn}}\|A\|_F \leq \|A\|_1$.

??? proof "Proof"
    **(1)** Already proved in Theorem 15.3 (5).

    **(4)** More directly: $\|A\|_2^2 = \rho(A^*A) \leq \|A^*A\|_\infty \leq \|A^*\|_\infty\|A\|_\infty = \|A\|_1\|A\|_\infty$. $\blacksquare$

!!! example "Example 15.6"
    For $A = \begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$:

    - $\|A\|_1 = \max\{4, 6\} = 6$;
    - $\|A\|_\infty = \max\{3, 7\} = 7$;
    - $\|A\|_F = \sqrt{30} \approx 5.48$;
    - $\|A\|_2 \approx 5.46$ (largest singular value).

    Verification of inequality (4): $\|A\|_2^2 \approx 29.87 \leq 42 = 6 \times 7 = \|A\|_1\|A\|_\infty$.

---

## 15.5 Condition number

<div class="context-flow" markdown>

**Core**: $\kappa(A)=\|A\|\|A^{-1}\|=\sigma_1/\sigma_n$ · The condition number is a "perturbation amplifier": $\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|}\leq\kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}$ · For normal matrices $\kappa_2=\rho(A)/|\lambda_{\min}|$

</div>

The condition number measures the inherent sensitivity of a linear system to input perturbations, and is one of the most important concepts in numerical linear algebra.

!!! definition "Definition 15.6 (Condition number)"
    Let $A \in \mathbb{C}^{n \times n}$ be invertible and $\|\cdot\|$ a matrix norm. The **condition number** of $A$ with respect to $\|\cdot\|$ is defined as

    $$\kappa(A) = \|A\|\cdot\|A^{-1}\|.$$

    For the $\ell_2$ norm, $\kappa_2(A) = \|A\|_2\|A^{-1}\|_2 = \frac{\sigma_1}{\sigma_n}$ (ratio of the largest to the smallest singular value).

!!! theorem "Theorem 15.7 (Basic properties of the condition number)"
    Let $A$ be invertible. Then:

    (1) $\kappa(A) \geq 1$;

    (2) $\kappa(\alpha A) = \kappa(A)$ for any $\alpha \neq 0$;

    (3) (Property regarding other invertible matrices of the same size);

    (4) If $U$ is unitary, then $\kappa_2(U) = 1$;

    (5) $\kappa_2(A) = \kappa_2(A^*)$.

??? proof "Proof"
    **(1)** $\kappa(A) = \|A\|\|A^{-1}\| \geq \|AA^{-1}\| = \|I\| = 1$.

    **(2)** $\kappa(\alpha A) = \|\alpha A\|\|(\alpha A)^{-1}\| = |\alpha|\|A\|\frac{1}{|\alpha|}\|A^{-1}\| = \kappa(A)$.

    **(4)** $\|U\|_2 = 1$ (all singular values of a unitary matrix are $1$), $\|U^{-1}\|_2 = \|U^*\|_2 = 1$.

    **(5)** $\kappa_2(A^*) = \sigma_1(A^*)/\sigma_n(A^*) = \sigma_1(A)/\sigma_n(A) = \kappa_2(A)$. $\blacksquare$

!!! theorem "Theorem 15.8 (Perturbation theorem for linear systems)"
    Let $A$ be invertible, $A\mathbf{x} = \mathbf{b}$ ($\mathbf{b} \neq \mathbf{0}$). If $\mathbf{b}$ is perturbed to $\mathbf{b} + \delta\mathbf{b}$, the corresponding solution becomes $\mathbf{x} + \delta\mathbf{x}$; then

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.$$

    If $A$ is perturbed to $A + \delta A$ ($\|\delta A\| < \|A^{-1}\|^{-1}$), the solution becomes $\mathbf{x} + \delta\mathbf{x}$; then

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x} + \delta\mathbf{x}\|} \leq \kappa(A)\frac{\|\delta A\|}{\|A\|}.$$

??? proof "Proof"
    **Right-hand side perturbation**: $A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$, so $A\delta\mathbf{x} = \delta\mathbf{b}$, $\delta\mathbf{x} = A^{-1}\delta\mathbf{b}$.

    $$\|\delta\mathbf{x}\| = \|A^{-1}\delta\mathbf{b}\| \leq \|A^{-1}\|\|\delta\mathbf{b}\|.$$

    From $\|\mathbf{b}\| = \|A\mathbf{x}\| \leq \|A\|\|\mathbf{x}\|$, i.e., $\frac{1}{\|\mathbf{x}\|} \leq \frac{\|A\|}{\|\mathbf{b}\|}$. Therefore

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\|\delta\mathbf{b}\|\frac{\|A\|}{\|\mathbf{b}\|} = \kappa(A)\frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}.$$

    **Coefficient matrix perturbation**: $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} = A\mathbf{x}$. Expanding gives $A\delta\mathbf{x} + \delta A(\mathbf{x} + \delta\mathbf{x}) = \mathbf{0}$, i.e., $\delta\mathbf{x} = -A^{-1}\delta A(\mathbf{x} + \delta\mathbf{x})$. Therefore

    $$\|\delta\mathbf{x}\| \leq \|A^{-1}\|\|\delta A\|\|\mathbf{x}+\delta\mathbf{x}\|,$$

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}+\delta\mathbf{x}\|} \leq \|A^{-1}\|\|\delta A\| = \kappa(A)\frac{\|\delta A\|}{\|A\|}. \quad \blacksquare$$

!!! example "Example 15.7"
    The famous Hilbert matrix $H_n = \left(\frac{1}{i+j-1}\right)_{i,j=1}^{n}$ is a classic example of an ill-conditioned matrix.

    - $H_3 = \begin{pmatrix} 1 & 1/2 & 1/3 \\ 1/2 & 1/3 & 1/4 \\ 1/3 & 1/4 & 1/5 \end{pmatrix}$, $\kappa_2(H_3) \approx 524$;
    - $\kappa_2(H_5) \approx 4.77 \times 10^5$;
    - $\kappa_2(H_{10}) \approx 1.60 \times 10^{13}$.

    The condition number grows exponentially with $n$, meaning that small perturbations in the input data lead to huge changes in the solution when solving linear systems with Hilbert matrices.

!!! example "Example 15.8"
    Let $A = \begin{pmatrix} 1 & 1 \\ 1 & 1.001 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 2 \\ 2.001 \end{pmatrix}$, with solution $\mathbf{x} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

    Perturbing $\mathbf{b}' = \begin{pmatrix} 2 \\ 2.002 \end{pmatrix}$ (relative perturbation about $0.05\%$), the new solution is $\mathbf{x}' = \begin{pmatrix} 0 \\ 2 \end{pmatrix}$.

    Relative change in the solution: $\frac{\|\mathbf{x}'-\mathbf{x}\|_2}{\|\mathbf{x}\|_2} = \frac{\sqrt{2}}{\sqrt{2}} = 1 = 100\%$.

    $\kappa_2(A) \approx 4002$, far greater than $1$, indicating that this problem is extremely ill-conditioned.

---

## 15.6 Perturbation of eigenvalues

<div class="context-flow" markdown>

**Chapter arc**: **Bauer--Fike** ($\kappa(X)\|E\|$ bound) → **Weyl** (Hermitian matrices: $|\Delta\lambda_i|\leq\|E\|_2$) → **Wielandt--Hoffman** (Frobenius norm bound) · Normal matrices are most stable ($\kappa=1$)

</div>

!!! theorem "Theorem 15.9 (Bauer--Fike theorem)"
    Let $A \in \mathbb{C}^{n \times n}$ be diagonalizable, $A = X\Lambda X^{-1}$ ($\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$). Let $\mu$ be any eigenvalue of $A + E$. Then there exists an eigenvalue $\lambda_j$ of $A$ such that

    $$|\mu - \lambda_j| \leq \kappa_p(X)\|E\|_p,$$

    where $\kappa_p(X) = \|X\|_p\|X^{-1}\|_p$ is the condition number of the eigenvector matrix.

??? proof "Proof"
    If $\mu = \lambda_j$ for some $j$, the inequality is trivial. Otherwise $\mu$ is not an eigenvalue of $A$, so $\mu I - A$ is invertible. From $(A+E)\mathbf{v} = \mu\mathbf{v}$ ($\mathbf{v} \neq \mathbf{0}$), we get $(\mu I - A)\mathbf{v} = E\mathbf{v}$, i.e.,

    $$\mathbf{v} = (\mu I - A)^{-1}E\mathbf{v}.$$

    Therefore $1 \leq \|(\mu I - A)^{-1}E\|_p \leq \|(\mu I - A)^{-1}\|_p\|E\|_p$.

    From $A = X\Lambda X^{-1}$, $(\mu I - A)^{-1} = X(\mu I - \Lambda)^{-1}X^{-1}$. Hence

    $$\|(\mu I - A)^{-1}\|_p \leq \|X\|_p\|(\mu I - \Lambda)^{-1}\|_p\|X^{-1}\|_p.$$

    $(\mu I - \Lambda)^{-1} = \operatorname{diag}\left(\frac{1}{\mu-\lambda_1}, \ldots, \frac{1}{\mu-\lambda_n}\right)$, whose norm is $\max_j\frac{1}{|\mu-\lambda_j|} = \frac{1}{\min_j|\mu-\lambda_j|}$.

    Therefore $1 \leq \frac{\kappa_p(X)}{\min_j|\mu-\lambda_j|}\|E\|_p$, i.e., $\min_j|\mu-\lambda_j| \leq \kappa_p(X)\|E\|_p$. $\blacksquare$

!!! note "Note"
    The Bauer--Fike theorem shows that the sensitivity of eigenvalues to perturbations is controlled by the condition number $\kappa(X)$ of the eigenvector matrix. For normal matrices, $X$ can be taken to be unitary, so $\kappa_2(X) = 1$, and eigenvalues are least sensitive to perturbations.

!!! theorem "Theorem 15.10 (Weyl inequality -- Hermitian matrices)"
    Let $A, B$ be $n \times n$ Hermitian matrices with eigenvalues in decreasing order $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$ and $\lambda_1(B) \geq \cdots \geq \lambda_n(B)$. Let the eigenvalues of $A + B$ be $\lambda_1(A+B) \geq \cdots \geq \lambda_n(A+B)$. Then for all $i + j - 1 \leq n$,

    $$\lambda_{i+j-1}(A + B) \leq \lambda_i(A) + \lambda_j(B).$$

    In particular, taking $j = 1$: $\lambda_i(A+B) \leq \lambda_i(A) + \lambda_1(B)$.

    Taking $i = 1$: $\lambda_j(A+B) \leq \lambda_1(A) + \lambda_j(B)$.

??? proof "Proof"
    Using the minimax principle (Courant--Fischer theorem):

    $$\lambda_k(M) = \min_{\dim V = n-k+1}\max_{\substack{\mathbf{x} \in V \\ \|\mathbf{x}\| = 1}}\mathbf{x}^*M\mathbf{x}.$$

    Let $U$ be the $(n-i+1)$-dimensional subspace achieving $\lambda_i(A)$, and $W$ the $(n-j+1)$-dimensional subspace achieving $\lambda_j(B)$. Set $V = U \cap W$; its dimension is $\geq (n-i+1) + (n-j+1) - n = n-i-j+2$, so $\dim V \geq n-(i+j-1)+1$.

    For $\mathbf{x} \in V$, $\|\mathbf{x}\| = 1$:

    $$\mathbf{x}^*(A+B)\mathbf{x} = \mathbf{x}^*A\mathbf{x} + \mathbf{x}^*B\mathbf{x} \leq \lambda_i(A) + \lambda_j(B).$$

    By the minimax principle, $\lambda_{i+j-1}(A+B) \leq \max_{\mathbf{x} \in V, \|\mathbf{x}\|=1}\mathbf{x}^*(A+B)\mathbf{x} \leq \lambda_i(A) + \lambda_j(B)$. $\blacksquare$

!!! theorem "Theorem 15.11 (Wielandt--Hoffman theorem)"
    Let $A, B$ be $n \times n$ Hermitian matrices with eigenvalues $\lambda_1(A) \geq \cdots \geq \lambda_n(A)$ and $\lambda_1(B) \geq \cdots \geq \lambda_n(B)$. Then

    $$\sum_{i=1}^{n}(\lambda_i(A) - \lambda_i(B))^2 \leq \|A - B\|_F^2.$$

??? proof "Proof"
    Let $A = U\Lambda_A U^*$, $B = V\Lambda_B V^*$ (spectral decompositions). Set $W = U^*V$; then $W$ is unitary.

    $$\|A - B\|_F^2 = \|\Lambda_A - W\Lambda_B W^*\|_F^2 = \operatorname{tr}(\Lambda_A^2) + \operatorname{tr}(\Lambda_B^2) - 2\operatorname{Re}\operatorname{tr}(\Lambda_A W\Lambda_B W^*).$$

    By the von Neumann trace inequality, $\operatorname{Re}\operatorname{tr}(\Lambda_A W\Lambda_B W^*) \leq \sum_i\lambda_i(A)\lambda_i(B)$ (equality when $W$ is a permutation matrix). Therefore

    $$\|A-B\|_F^2 \geq \sum\lambda_i(A)^2 + \sum\lambda_i(B)^2 - 2\sum\lambda_i(A)\lambda_i(B) = \sum(\lambda_i(A)-\lambda_i(B))^2. \quad \blacksquare$$

!!! example "Example 15.9"
    Let $A = \begin{pmatrix} 5 & 1 \\ 1 & 3 \end{pmatrix}$, $E = \begin{pmatrix} 0.1 & 0 \\ 0 & -0.1 \end{pmatrix}$.

    Eigenvalues of $A$: $\lambda = 4 \pm \sqrt{2}$, i.e., $\lambda_1 \approx 5.41$, $\lambda_2 \approx 2.59$.

    $A + E = \begin{pmatrix} 5.1 & 1 \\ 1 & 2.9 \end{pmatrix}$, eigenvalues $\lambda = 4 \pm \sqrt{2.21}$, i.e., $\lambda_1' \approx 5.49$, $\lambda_2' \approx 2.51$.

    $\max_i|\lambda_i' - \lambda_i| \approx 0.08 \leq \|E\|_2 = 0.1$ (Weyl inequality for Hermitian matrices gives a $\kappa_2 = 1$ estimate).

---

## 15.7 Perturbation of singular values

<div class="context-flow" markdown>

**Intuition**: Construct the Hermitian dilation $\hat{A}=\begin{pmatrix}0&A\\A^*&0\end{pmatrix}$ to reduce singular value perturbation to eigenvalue perturbation → Ch18 unified inequality framework

</div>

!!! theorem "Theorem 15.12 (Weyl inequality for singular values)"
    Let $A, B \in \mathbb{C}^{m \times n}$ with singular values in decreasing order $\sigma_1(A) \geq \cdots \geq \sigma_p(A)$ and $\sigma_1(B) \geq \cdots \geq \sigma_p(B)$ ($p = \min(m,n)$). Then for each $i$,

    $$|\sigma_i(A) - \sigma_i(B)| \leq \|A - B\|_2.$$

    More strongly,

    $$\sqrt{\sum_{i=1}^{p}(\sigma_i(A) - \sigma_i(B))^2} \leq \|A - B\|_F.$$

??? proof "Proof"
    Construct the $(m+n) \times (m+n)$ Hermitian matrices

    $$\hat{A} = \begin{pmatrix} O & A \\ A^* & O \end{pmatrix}, \quad \hat{B} = \begin{pmatrix} O & B \\ B^* & O \end{pmatrix}.$$

    The eigenvalues of $\hat{A}$ are $\pm\sigma_1(A), \ldots, \pm\sigma_p(A)$ and $|m-n|$ zeros (if $m \neq n$). Similarly for $\hat{B}$.

    Applying the Weyl inequality for Hermitian matrices to $\hat{A}$ and $\hat{B}$:

    $$|\sigma_i(A) - \sigma_i(B)| = |\lambda_i(\hat{A}) - \lambda_i(\hat{B})| \leq \|\hat{A} - \hat{B}\|_2 = \|A - B\|_2.$$

    The Frobenius norm version follows similarly from the Wielandt--Hoffman theorem. $\blacksquare$

!!! proposition "Proposition 15.2 (Subadditivity of singular values)"
    Let $A, B \in \mathbb{C}^{m \times n}$. Then for each $i + j - 1 \leq \min(m,n)$,

    $$\sigma_{i+j-1}(A + B) \leq \sigma_i(A) + \sigma_j(B).$$

!!! example "Example 15.10"
    Let $A = \begin{pmatrix} 3 & 0 \\ 0 & 1 \end{pmatrix}$, $E = \begin{pmatrix} 0 & 0.2 \\ 0.2 & 0 \end{pmatrix}$.

    $\sigma_1(A) = 3$, $\sigma_2(A) = 1$. $B = A + E = \begin{pmatrix} 3 & 0.2 \\ 0.2 & 1 \end{pmatrix}$.

    $B^TB = \begin{pmatrix} 9.04 & 0.8 \\ 0.8 & 1.04 \end{pmatrix}$, eigenvalues $\mu = 5.04 \pm \sqrt{16.64}$. $\mu_1 \approx 9.12$, $\mu_2 \approx 0.96$. $\sigma_1(B) \approx 3.02$, $\sigma_2(B) \approx 0.98$.

    $|\sigma_1(B) - \sigma_1(A)| \approx 0.02 \leq 0.2 = \|E\|_2$. $|\sigma_2(B) - \sigma_2(A)| \approx 0.02 \leq 0.2$.

---

## 15.8 Perturbation analysis of linear systems

<div class="context-flow" markdown>

**Chapter arc**: Small residual $\not\Rightarrow$ small error (condition number amplification) · $\frac{1}{\kappa}\cdot\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}\leq\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|}\leq\kappa\cdot\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}$ → two-sided bound

</div>

!!! definition "Definition 15.7 (Forward error and backward error)"
    Let $\hat{\mathbf{x}}$ be an approximate solution of the linear system $A\mathbf{x} = \mathbf{b}$.

    - **Forward error**: $\|\hat{\mathbf{x}} - \mathbf{x}\|$ (error in the solution);
    - **Backward error**: $\|\mathbf{r}\| = \|\mathbf{b} - A\hat{\mathbf{x}}\|$ (norm of the residual);
    - **Relative forward error**: $\frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|}$;
    - **Relative backward error**: $\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}$.

!!! theorem "Theorem 15.13 (Relation between forward and backward errors)"
    Let $A$ be invertible and $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}}$ be the residual. Then

    $$\frac{1}{\kappa(A)}\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \frac{\|\hat{\mathbf{x}} - \mathbf{x}\|}{\|\mathbf{x}\|} \leq \kappa(A)\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}.$$

??? proof "Proof"
    From $A(\hat{\mathbf{x}} - \mathbf{x}) = A\hat{\mathbf{x}} - \mathbf{b} = -\mathbf{r}$, i.e., $\hat{\mathbf{x}} - \mathbf{x} = -A^{-1}\mathbf{r}$.

    **Upper bound**: $\|\hat{\mathbf{x}}-\mathbf{x}\| \leq \|A^{-1}\|\|\mathbf{r}\|$, $\|\mathbf{b}\| \leq \|A\|\|\mathbf{x}\|$, therefore

    $$\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\frac{\|\mathbf{r}\|}{\|\mathbf{x}\|} \leq \|A^{-1}\|\|A\|\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} = \kappa(A)\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|}.$$

    **Lower bound**: $\|\mathbf{r}\| = \|A(\hat{\mathbf{x}}-\mathbf{x})\| \leq \|A\|\|\hat{\mathbf{x}}-\mathbf{x}\|$, $\|\mathbf{x}\| \leq \|A^{-1}\|\|\mathbf{b}\|$, therefore

    $$\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \|A\|\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{b}\|} \leq \|A\|\|A^{-1}\|\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|} = \kappa(A)\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|}.$$

    That is, $\frac{1}{\kappa(A)}\frac{\|\mathbf{r}\|}{\|\mathbf{b}\|} \leq \frac{\|\hat{\mathbf{x}}-\mathbf{x}\|}{\|\mathbf{x}\|}$. $\blacksquare$

!!! theorem "Theorem 15.14 (Simultaneous perturbation of $A$ and $\mathbf{b}$)"
    Let $(A + \delta A)(\mathbf{x} + \delta\mathbf{x}) = \mathbf{b} + \delta\mathbf{b}$, and $\kappa(A)\frac{\|\delta A\|}{\|A\|} < 1$. Then

    $$\frac{\|\delta\mathbf{x}\|}{\|\mathbf{x}\|} \leq \frac{\kappa(A)}{1 - \kappa(A)\frac{\|\delta A\|}{\|A\|}}\left(\frac{\|\delta A\|}{\|A\|} + \frac{\|\delta\mathbf{b}\|}{\|\mathbf{b}\|}\right).$$

!!! example "Example 15.11"
    Solving $A\mathbf{x} = \mathbf{b}$ in floating-point arithmetic, where $A = \begin{pmatrix} 1 & 1 \\ 1 & 1.0001 \end{pmatrix}$, $\mathbf{b} = \begin{pmatrix} 2 \\ 2.0001 \end{pmatrix}$.

    Exact solution $\mathbf{x} = \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.

    Suppose the computed solution is $\hat{\mathbf{x}} = \begin{pmatrix} 1.5 \\ 0.5 \end{pmatrix}$. Residual $\mathbf{r} = \mathbf{b} - A\hat{\mathbf{x}} = \begin{pmatrix} 0 \\ 0.00005 \end{pmatrix}$.

    Relative residual $\frac{\|\mathbf{r}\|_\infty}{\|\mathbf{b}\|_\infty} = \frac{0.00005}{2.0001} \approx 2.5 \times 10^{-5}$ (very small).

    But the relative forward error $\frac{\|\hat{\mathbf{x}}-\mathbf{x}\|_\infty}{\|\mathbf{x}\|_\infty} = \frac{0.5}{1} = 0.5 = 50\%$ (very large).

    This illustrates the amplification effect of the condition number. $\kappa_\infty(A) \approx 4 \times 10^4$, causing a small residual to correspond to a large error.

!!! example "Example 15.12"
    Comparing well-conditioned and ill-conditioned systems.

    **Well-conditioned system**: $A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}$, $\kappa_2(A) = 2$. A perturbation $\|\delta\mathbf{b}\|/\|\mathbf{b}\| = 10^{-6}$ leads to $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \leq 2 \times 10^{-6}$.

    **Ill-conditioned system**: $B = \begin{pmatrix} 1 & 1 \\ 1 & 1+10^{-10} \end{pmatrix}$, $\kappa_2(B) \approx 4 \times 10^{10}$. The same perturbation can lead to $\|\delta\mathbf{x}\|/\|\mathbf{x}\| \leq 4 \times 10^4$, meaning the solution is completely unreliable.
