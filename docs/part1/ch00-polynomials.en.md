# Chapter 0  Polynomial Algebra

<div class="context-flow" markdown>

**Prerequisites**: Basic concepts of fields ($\mathbb{Q}, \mathbb{R}, \mathbb{C}$) · **Chapter arc**: Polynomial rings → Divisibility and GCD → Coprime polynomials → Irreducible polynomials → Unique factorization → Roots → Multiple roots and discriminant → Irreducibles over $\mathbb{R}$ and $\mathbb{C}$ → Irreducibles over $\mathbb{Q}$ → Polynomial interpolation
**One-line essence**: Polynomial algebra is the preparatory language of linear algebra — characteristic polynomials, minimal polynomials, and Jordan normal forms all rest on the divisibility and factorization theory of polynomials

</div>

Polynomials are among the most fundamental objects in algebra. In linear algebra, polynomials play a central role: the characteristic polynomial determines eigenvalues, the minimal polynomial captures the essential structure of a linear transformation, and the Cayley-Hamilton theorem ties matrices intimately to polynomials. This chapter systematically develops the algebraic theory of polynomials in one variable, laying a solid foundation for subsequent chapters.

---

## 0.1 Polynomial Rings

<div class="context-flow" markdown>

**Core question**: How do we rigorously define a polynomial? → Distinguish polynomials (formal expressions) from polynomial functions → The totality of polynomials forms a ring

</div>

### Basic definitions

Let $\mathbb{F}$ be a field (in this course, typically $\mathbb{F} = \mathbb{Q}, \mathbb{R}$, or $\mathbb{C}$).

!!! definition "Definition 0.1 (Polynomial)"
    Let $\mathbb{F}$ be a field and $x$ an indeterminate. An expression of the form

    $$
    f(x) = a_n x^n + a_{n-1} x^{n-1} + \cdots + a_1 x + a_0
    $$

    is called a **polynomial in one variable** over $\mathbb{F}$, where $a_0, a_1, \ldots, a_n \in \mathbb{F}$ are called **coefficients**. If $a_n \neq 0$, then $n$ is called the **degree** of $f(x)$, written $\deg f = n$, and $a_n$ is the **leading coefficient**. If $a_n = 1$, then $f(x)$ is called a **monic polynomial**. The degree of the zero polynomial $f(x) = 0$ is defined to be $-\infty$.

!!! note "Note"
    A polynomial is a **formal expression**, not a function. Two polynomials are equal if and only if all corresponding coefficients are equal. Over infinite fields (such as $\mathbb{Q}, \mathbb{R}, \mathbb{C}$), polynomials and polynomial functions can be identified; but over finite fields (such as $\mathbb{F}_p$), the polynomial $x^p - x$ is identically zero as a function, yet it is not the zero polynomial.

!!! definition "Definition 0.2 (Polynomial ring)"
    The set of all polynomials in one variable over $\mathbb{F}$, equipped with the usual addition and multiplication, is denoted $\mathbb{F}[x]$ and called the **polynomial ring** over $\mathbb{F}$.

### Operations on polynomials

Let $f(x) = \sum_{i=0}^n a_i x^i$ and $g(x) = \sum_{j=0}^m b_j x^j$.

**Addition**: $f(x) + g(x) = \sum_{k=0}^{\max(n,m)} (a_k + b_k) x^k$ (padding higher-order coefficients with zero as needed).

**Multiplication**: $f(x) \cdot g(x) = \sum_{k=0}^{n+m} c_k x^k$, where $c_k = \sum_{i+j=k} a_i b_j$.

!!! theorem "Theorem 0.1 (Degree formula)"
    Let $f(x), g(x) \in \mathbb{F}[x]$ be nonzero polynomials. Then

    $$
    \deg(f \cdot g) = \deg f + \deg g.
    $$

    In particular, $\mathbb{F}[x]$ is an **integral domain**: if $f(x) \cdot g(x) = 0$, then $f(x) = 0$ or $g(x) = 0$.

??? proof "Proof"
    Let $\deg f = n$, $\deg g = m$, with leading coefficients $a_n \neq 0$ and $b_m \neq 0$ respectively. The highest-degree term of $f \cdot g$ is $a_n b_m x^{n+m}$. Since $\mathbb{F}$ is a field (hence an integral domain), $a_n b_m \neq 0$, so $\deg(f \cdot g) = n + m$.

    If $f \cdot g = 0$, then $\deg(f \cdot g) = -\infty$. If both $f, g$ were nonzero, then $\deg(f \cdot g) = \deg f + \deg g \geq 0$, a contradiction. $\blacksquare$

!!! theorem "Theorem 0.2 (Division algorithm)"
    Let $f(x), g(x) \in \mathbb{F}[x]$ with $g(x) \neq 0$. Then there exist unique $q(x), r(x) \in \mathbb{F}[x]$ such that

    $$
    f(x) = q(x) \cdot g(x) + r(x),
    $$

    where $\deg r < \deg g$ (or $r = 0$). Here $q(x)$ is called the **quotient** and $r(x)$ the **remainder**.

??? proof "Proof"
    **Existence**: By induction on $\deg f$. If $\deg f < \deg g$, take $q = 0$ and $r = f$. Otherwise, let $\deg f = n \geq m = \deg g$, with leading coefficients $a_n$ and $b_m$ respectively. Set

    $$
    f_1(x) = f(x) - \frac{a_n}{b_m} x^{n-m} g(x).
    $$

    Then $\deg f_1 < n$. By the induction hypothesis, $f_1 = q_1 g + r$ with $\deg r < \deg g$. Hence

    $$
    f = \left(\frac{a_n}{b_m} x^{n-m} + q_1\right) g + r.
    $$

    **Uniqueness**: Suppose $f = q_1 g + r_1 = q_2 g + r_2$ with $\deg r_1, \deg r_2 < \deg g$. Then $(q_1 - q_2)g = r_2 - r_1$. If $q_1 \neq q_2$, the left side has degree $\geq \deg g$, while the right side has degree $< \deg g$, a contradiction. So $q_1 = q_2$ and hence $r_1 = r_2$. $\blacksquare$

!!! example "Example 0.1"
    In $\mathbb{R}[x]$, let $f(x) = 2x^4 + 3x^3 - x + 5$ and $g(x) = x^2 + 1$. Find the quotient and remainder of $f$ divided by $g$.

    Performing polynomial long division:

    $$
    2x^4 + 3x^3 - x + 5 = (2x^2 + 3x - 2)(x^2 + 1) + (-4x + 7).
    $$

    Thus $q(x) = 2x^2 + 3x - 2$ and $r(x) = -4x + 7$. Check: $\deg r = 1 < 2 = \deg g$.

---

## 0.2 Divisibility and GCD

<div class="context-flow" markdown>

**Core question**: Divisibility relations among polynomials → Existence and computation of the GCD → The Euclidean algorithm

</div>

!!! definition "Definition 0.3 (Divisibility)"
    Let $f(x), g(x) \in \mathbb{F}[x]$. If there exists $q(x) \in \mathbb{F}[x]$ such that $f(x) = q(x) \cdot g(x)$, we say $g(x)$ **divides** $f(x)$, written $g(x) \mid f(x)$. In this case $g(x)$ is called a **factor** (or **divisor**) of $f(x)$, and $f(x)$ is a **multiple** of $g(x)$.

!!! definition "Definition 0.4 (Associates)"
    If $f(x) \mid g(x)$ and $g(x) \mid f(x)$, then $f(x)$ and $g(x)$ are called **associates**, written $f(x) \sim g(x)$. In $\mathbb{F}[x]$, $f \sim g$ if and only if $f = cg$ for some $c \in \mathbb{F}^*$ (nonzero constant).

!!! definition "Definition 0.5 (Greatest common divisor)"
    Let $f(x), g(x) \in \mathbb{F}[x]$, not both zero. A polynomial $d(x) \in \mathbb{F}[x]$ is called a **greatest common divisor (GCD)** of $f(x)$ and $g(x)$ if:

    1. $d(x) \mid f(x)$ and $d(x) \mid g(x)$ ($d$ is a common divisor);
    2. if $c(x) \mid f(x)$ and $c(x) \mid g(x)$, then $c(x) \mid d(x)$ ($d$ is the greatest).

    The unique monic GCD is denoted $\gcd(f, g)$ or $(f, g)$.

!!! theorem "Theorem 0.3 (Existence of GCD and Bezout's identity)"
    Let $f(x), g(x) \in \mathbb{F}[x]$, not both zero. Then:

    1. $\gcd(f, g)$ exists and is unique (up to associates).
    2. There exist $u(x), v(x) \in \mathbb{F}[x]$ such that

    $$
    \gcd(f, g) = u(x) f(x) + v(x) g(x).
    $$

    This is called **Bezout's identity**.

??? proof "Proof"
    Consider the set $I = \{u(x)f(x) + v(x)g(x) : u, v \in \mathbb{F}[x]\}$. This is an ideal of $\mathbb{F}[x]$. Let $d(x)$ be a nonzero element of $I$ of minimal degree (normalized to be monic).

    For any $h(x) \in I$, by the division algorithm $h = qd + r$ with $\deg r < \deg d$. Since $h, d \in I$, we have $r = h - qd \in I$. By minimality of $\deg d$, we get $r = 0$, so $d \mid h$.

    In particular, $f, g \in I$, so $d \mid f$ and $d \mid g$. If $c \mid f$ and $c \mid g$, then from $d = uf + vg$ we get $c \mid d$. Hence $d = \gcd(f, g)$. $\blacksquare$

### The Euclidean Algorithm

!!! theorem "Theorem 0.4 (Euclidean algorithm)"
    Let $f(x), g(x) \in \mathbb{F}[x]$, $g \neq 0$. Repeatedly apply the division algorithm:

    $$
    \begin{aligned}
    f &= q_1 g + r_1, & \deg r_1 &< \deg g, \\
    g &= q_2 r_1 + r_2, & \deg r_2 &< \deg r_1, \\
    r_1 &= q_3 r_2 + r_3, & \deg r_3 &< \deg r_2, \\
    &\;\;\vdots \\
    r_{k-2} &= q_k r_{k-1} + r_k, & \deg r_k &< \deg r_{k-1}, \\
    r_{k-1} &= q_{k+1} r_k.
    \end{aligned}
    $$

    Then $\gcd(f, g) \sim r_k$ (the last nonzero remainder).

??? proof "Proof"
    From $f = q_1 g + r_1$, a polynomial $d$ divides both $f$ and $g$ if and only if it divides both $g$ and $r_1$. Thus $\gcd(f, g) = \gcd(g, r_1) = \gcd(r_1, r_2) = \cdots = \gcd(r_{k-1}, r_k) = r_k$ (normalized to monic).

    Since $\deg g > \deg r_1 > \deg r_2 > \cdots$ is a strictly decreasing sequence of nonneg integers, the algorithm terminates in finitely many steps. $\blacksquare$

!!! example "Example 0.2"
    Find $\gcd(x^4 - 1, x^3 - 1)$ in $\mathbb{Q}[x]$.

    $$
    x^4 - 1 = x \cdot (x^3 - 1) + (x - 1),
    $$

    $$
    x^3 - 1 = (x^2 + x + 1)(x - 1) + 0.
    $$

    Therefore $\gcd(x^4 - 1, x^3 - 1) = x - 1$.

!!! example "Example 0.3"
    Find $\gcd(f, g)$ and the Bezout coefficients, where $f = x^3 + x + 1$ and $g = x^2 + x$ in $\mathbb{Q}[x]$.

    $$
    x^3 + x + 1 = (x - 1)(x^2 + x) + (2x + 1),
    $$

    $$
    x^2 + x = \left(\tfrac{1}{2}x + \tfrac{1}{4}\right)(2x + 1) + \tfrac{3}{4}.
    $$

    The remainder $\frac{3}{4}$ is a nonzero constant, so $\gcd(f, g) = 1$. Back-substitution yields Bezout coefficients $u(x), v(x)$ satisfying $u f + v g = 1$.

---

## 0.3 Coprime Polynomials

<div class="context-flow" markdown>

**Core question**: What does $\gcd(f, g) = 1$ mean? → Equivalent conditions for coprimality → Key properties of coprime polynomials

</div>

!!! definition "Definition 0.6 (Coprime)"
    If $\gcd(f(x), g(x)) = 1$, then $f(x)$ and $g(x)$ are called **coprime** (or **relatively prime**).

!!! theorem "Theorem 0.5 (Equivalent condition for coprimality)"
    $f(x)$ and $g(x)$ are coprime if and only if there exist $u(x), v(x) \in \mathbb{F}[x]$ such that

    $$
    u(x) f(x) + v(x) g(x) = 1.
    $$

??? proof "Proof"
    **Necessity**: If $\gcd(f, g) = 1$, this follows directly from Bezout's identity.

    **Sufficiency**: If $uf + vg = 1$, let $d = \gcd(f, g)$. Since $d \mid f$ and $d \mid g$, we have $d \mid (uf + vg) = 1$, so $d$ is a constant and $\gcd(f, g) = 1$. $\blacksquare$

!!! theorem "Theorem 0.6 (Multiplicative properties of coprimality)"
    Let $f(x), g(x), h(x) \in \mathbb{F}[x]$.

    1. If $f \mid gh$ and $\gcd(f, g) = 1$, then $f \mid h$.
    2. If $f \mid h$ and $g \mid h$, and $\gcd(f, g) = 1$, then $fg \mid h$.

??? proof "Proof"
    1. Since $\gcd(f, g) = 1$, there exist $u, v$ with $uf + vg = 1$. Multiplying by $h$: $ufh + vgh = h$. Since $f \mid ufh$ and $f \mid gh$ (so $f \mid vgh$), we get $f \mid h$.

    2. Since $f \mid h$, write $h = fk$. From $g \mid h = fk$ and $\gcd(f, g) = 1$, part (1) gives $g \mid k$, say $k = gl$. Then $h = fgl$, so $fg \mid h$. $\blacksquare$

!!! theorem "Theorem 0.7 (Pairwise coprime polynomials)"
    Let $f_1, f_2, \ldots, f_k \in \mathbb{F}[x]$ be pairwise coprime, and suppose each $f_i \mid h$. Then $f_1 f_2 \cdots f_k \mid h$.

??? proof "Proof"
    By induction on $k$. The case $k = 2$ is Theorem 0.6 (2). Assuming the result for $k - 1$, we have $f_1 \cdots f_{k-1} \mid h$. From $\gcd(f_i, f_k) = 1$ for $i = 1, \ldots, k-1$, one can show $\gcd(f_1 \cdots f_{k-1}, f_k) = 1$. Then Theorem 0.6 (2) gives $f_1 \cdots f_k \mid h$. $\blacksquare$

!!! example "Example 0.4"
    Let $f(x) = x^2 - 1 = (x-1)(x+1)$ and $g(x) = x^2 + x = x(x+1)$. Then

    $$
    \gcd(f, g) = x + 1 \neq 1,
    $$

    so $f$ and $g$ are not coprime. However, $x - 1$ and $x$ are coprime (their GCD is $1$).

---

## 0.4 Irreducible Polynomials

<div class="context-flow" markdown>

**Core question**: Which polynomials cannot be factored further? → Irreducible polynomials — the "primes" of the polynomial ring → Fundamental properties

</div>

!!! definition "Definition 0.7 (Irreducible polynomial)"
    Let $p(x) \in \mathbb{F}[x]$ have degree $\geq 1$. If $p(x)$ cannot be written as the product of two polynomials each of degree less than $\deg p$, then $p(x)$ is called an **irreducible polynomial** over $\mathbb{F}$. Otherwise $p(x)$ is called **reducible**.

    Equivalently, $p(x)$ is irreducible if and only if $p(x) = f(x)g(x)$ implies $f(x)$ or $g(x)$ is a constant.

!!! note "Note"
    Irreducibility depends on the base field $\mathbb{F}$. For example, $x^2 + 1$ is irreducible over $\mathbb{R}$, but reducible over $\mathbb{C}$ since $x^2 + 1 = (x + i)(x - i)$.

!!! theorem "Theorem 0.8 (Property of irreducible polynomials)"
    Let $p(x) \in \mathbb{F}[x]$ be irreducible and $f(x) \in \mathbb{F}[x]$. Then either $p \mid f$ or $\gcd(p, f) = 1$.

??? proof "Proof"
    Let $d = \gcd(p, f)$. Then $d \mid p$, i.e., $p = d \cdot q$. Since $p$ is irreducible, either $d$ or $q$ is a constant. If $d$ is a constant, then $\gcd(p, f) = 1$; if $q$ is a constant, then $d \sim p$, so $p \mid f$. $\blacksquare$

!!! theorem "Theorem 0.9 (Primality of irreducible polynomials)"
    Let $p(x) \in \mathbb{F}[x]$ be irreducible. If $p \mid f_1 f_2 \cdots f_k$, then $p \mid f_i$ for some $i$.

??? proof "Proof"
    By induction on $k$. The case $k = 1$ is trivial. For $k \geq 2$, suppose $p \mid f_1 \cdots f_k$. If $p \mid f_k$, we are done. Otherwise $\gcd(p, f_k) = 1$ (Theorem 0.8), so Theorem 0.6 (1) gives $p \mid f_1 \cdots f_{k-1}$, and the induction hypothesis applies. $\blacksquare$

!!! example "Example 0.5"
    In $\mathbb{R}[x]$, $x^2 + 1$ is irreducible. If $(x^2+1) \mid f(x)g(x)$, then $(x^2+1) \mid f(x)$ or $(x^2+1) \mid g(x)$.

    In $\mathbb{Q}[x]$, $x^2 - 2$ is irreducible (if it were reducible, it would have a rational root $\pm\sqrt{2}$, but $\sqrt{2} \notin \mathbb{Q}$).

---

## 0.5 Factorization of Polynomials

<div class="context-flow" markdown>

**Core question**: Can every polynomial be uniquely factored into irreducibles? → Analogy with the fundamental theorem of arithmetic → Unique factorization theorem

</div>

!!! theorem "Theorem 0.10 (Unique factorization theorem)"
    Every polynomial $f(x) \in \mathbb{F}[x]$ of degree $\geq 1$ can be factored as

    $$
    f(x) = c \cdot p_1(x)^{e_1} p_2(x)^{e_2} \cdots p_s(x)^{e_s},
    $$

    where $c \in \mathbb{F}^*$, $p_1, p_2, \ldots, p_s$ are pairwise non-associate monic irreducible polynomials, and $e_1, e_2, \ldots, e_s$ are positive integers. This factorization is unique up to the ordering of the factors.

??? proof "Proof"
    **Existence**: By induction on $\deg f$. If $\deg f = 1$, then $f$ is itself irreducible. For $\deg f \geq 2$: if $f$ is irreducible, the factorization is just $f$ itself. If $f$ is reducible, then $f = g \cdot h$ with $1 \leq \deg g, \deg h < \deg f$. By the induction hypothesis, $g$ and $h$ each have irreducible factorizations, and combining them gives a factorization of $f$.

    **Uniqueness**: Suppose $f = c \cdot p_1^{e_1} \cdots p_s^{e_s} = c' \cdot q_1^{d_1} \cdots q_t^{d_t}$. Since $p_1 \mid f = c' \cdot q_1^{d_1} \cdots q_t^{d_t}$ and $p_1$ is irreducible (hence prime by Theorem 0.9), $p_1 \mid q_j$ for some $j$. Since $q_j$ is also irreducible, $p_1 \sim q_j$. Pairing $p_1$ with $q_j$ and cancelling, we proceed by induction to obtain $s = t$ and, after reordering, $p_i \sim q_i$ and $e_i = d_i$. $\blacksquare$

!!! example "Example 0.6"
    Factor $f(x) = x^4 - 1$ in $\mathbb{Q}[x]$.

    $$
    x^4 - 1 = (x^2 - 1)(x^2 + 1) = (x - 1)(x + 1)(x^2 + 1).
    $$

    In $\mathbb{Q}[x]$, the factors $x - 1$, $x + 1$, and $x^2 + 1$ are all irreducible, giving the unique factorization.

    In $\mathbb{C}[x]$, $x^2 + 1 = (x - i)(x + i)$, so

    $$
    x^4 - 1 = (x - 1)(x + 1)(x - i)(x + i).
    $$

!!! example "Example 0.7"
    Factor $f(x) = x^6 - 1$ in $\mathbb{R}[x]$.

    $$
    x^6 - 1 = (x^3 - 1)(x^3 + 1) = (x-1)(x^2+x+1)(x+1)(x^2-x+1).
    $$

    The discriminants of $x^2 + x + 1$ and $x^2 - x + 1$ are both $-3 < 0$, so they are irreducible over $\mathbb{R}$.

---

## 0.6 Roots of Polynomials

<div class="context-flow" markdown>

**Core question**: When does a polynomial have a root? → Remainder theorem and factor theorem → Upper bound on the number of roots → Roots and factorization

</div>

!!! definition "Definition 0.8 (Root)"
    Let $f(x) \in \mathbb{F}[x]$ and $\alpha \in \mathbb{F}$. If $f(\alpha) = 0$, then $\alpha$ is called a **root** (or **zero**) of $f(x)$.

!!! theorem "Theorem 0.11 (Remainder theorem)"
    Let $f(x) \in \mathbb{F}[x]$ and $\alpha \in \mathbb{F}$. The remainder when $f(x)$ is divided by $(x - \alpha)$ equals $f(\alpha)$.

??? proof "Proof"
    By the division algorithm, $f(x) = q(x)(x - \alpha) + r$, where $r \in \mathbb{F}$ (since $\deg r < \deg(x - \alpha) = 1$, so $r$ is a constant). Setting $x = \alpha$: $f(\alpha) = q(\alpha) \cdot 0 + r = r$. $\blacksquare$

!!! theorem "Theorem 0.12 (Factor theorem)"
    $\alpha$ is a root of $f(x)$ if and only if $(x - \alpha) \mid f(x)$.

??? proof "Proof"
    By the remainder theorem, $f(x) = q(x)(x - \alpha) + f(\alpha)$. So $f(\alpha) = 0$ if and only if $f(x) = q(x)(x - \alpha)$, i.e., $(x - \alpha) \mid f(x)$. $\blacksquare$

!!! definition "Definition 0.9 (Multiplicity of a root)"
    Let $\alpha$ be a root of $f(x)$. If $(x - \alpha)^k \mid f(x)$ but $(x - \alpha)^{k+1} \nmid f(x)$, then $\alpha$ is called a **root of multiplicity $k$**. When $k = 1$, it is a **simple root**; when $k \geq 2$, a **multiple root**.

!!! theorem "Theorem 0.13 (Number of roots)"
    A polynomial of degree $n$ in $\mathbb{F}[x]$ has at most $n$ roots (counted with multiplicity).

??? proof "Proof"
    By induction on $n$. For $n = 0$, $f$ is a nonzero constant with no roots. For $n \geq 1$, suppose $f$ has a root $\alpha_1$ of multiplicity $k_1$. Then $f(x) = (x - \alpha_1)^{k_1} g(x)$ where $g(\alpha_1) \neq 0$ and $\deg g = n - k_1$. Every other root of $f$ must be a root of $g$ (since $(x - \alpha_1)^{k_1}$ is nonzero at $\alpha \neq \alpha_1$). By the induction hypothesis, $g$ has at most $n - k_1$ roots (with multiplicity), so $f$ has at most $k_1 + (n - k_1) = n$ roots. $\blacksquare$

!!! example "Example 0.8"
    $f(x) = x^3 - 3x + 2 = (x-1)^2(x+2)$.

    The roots are $x = 1$ (multiplicity $2$) and $x = -2$ (simple root), for a total of $3$ roots counted with multiplicity, equaling $\deg f = 3$.

!!! example "Example 0.9"
    Use the factor theorem to verify a special case of the Vandermonde determinant formula.

    Let $f(x) = \det \begin{pmatrix} 1 & 1 & 1 \\ a & b & x \\ a^2 & b^2 & x^2 \end{pmatrix}$. Expanding along the third column, $f(x)$ is a polynomial of degree $2$ in $x$.

    $f(a) = 0$ (first and third columns equal), $f(b) = 0$ (second and third columns equal). So $(x - a)(x - b) \mid f(x)$, and since $\deg f = 2$, we have $f(x) = c(x - a)(x - b)$. Comparing the coefficient of $x^2$ gives $c = a - b$ (details omitted), hence

    $$
    f(x) = (a - b)(x - a)(x - b).
    $$

---

## 0.7 Multiple Roots and the Discriminant

<div class="context-flow" markdown>

**Core question**: How to detect whether a polynomial has multiple roots? → Formal derivative → Criterion for multiple roots → The discriminant

</div>

!!! definition "Definition 0.10 (Formal derivative)"
    Let $f(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_1 x + a_0 \in \mathbb{F}[x]$. The **formal derivative** of $f(x)$ is defined as

    $$
    f'(x) = n a_n x^{n-1} + (n-1)a_{n-1}x^{n-2} + \cdots + a_1.
    $$

    This is a purely algebraic definition, independent of limits. The formal derivative satisfies the usual rules: $(f + g)' = f' + g'$ and $(fg)' = f'g + fg'$.

!!! theorem "Theorem 0.14 (Multiple root criterion)"
    Let $f(x) \in \mathbb{F}[x]$, $\deg f \geq 1$, and $\alpha \in \overline{\mathbb{F}}$ (the algebraic closure of $\mathbb{F}$). Then $\alpha$ is a root of multiplicity $k \geq 2$ if and only if $f(\alpha) = 0$ and $f'(\alpha) = 0$.

    More precisely, $\alpha$ is a root of multiplicity $k$ if and only if $f(\alpha) = f'(\alpha) = \cdots = f^{(k-1)}(\alpha) = 0$ and $f^{(k)}(\alpha) \neq 0$.

??? proof "Proof"
    Write $f(x) = (x - \alpha)^k g(x)$ with $g(\alpha) \neq 0$ and $k \geq 1$. Then

    $$
    f'(x) = k(x - \alpha)^{k-1} g(x) + (x - \alpha)^k g'(x) = (x - \alpha)^{k-1}[k g(x) + (x - \alpha)g'(x)].
    $$

    If $k \geq 2$: $f'(\alpha) = 0$, i.e., $(x - \alpha) \mid f'(x)$.

    If $k = 1$: $f'(\alpha) = 1 \cdot g(\alpha) \neq 0$.

    Conversely, if $f(\alpha) = 0$ and $f'(\alpha) = 0$, then $\alpha$ is a root of $f$ (multiplicity $\geq 1$), and the above computation shows the multiplicity is $\geq 2$. $\blacksquare$

!!! theorem "Theorem 0.15 (Criterion for no multiple roots)"
    $f(x)$ has no multiple roots (in $\overline{\mathbb{F}}$) if and only if $\gcd(f, f') = 1$.

??? proof "Proof"
    $f$ has a multiple root $\alpha$ ($k \geq 2$) $\Leftrightarrow$ $(x - \alpha) \mid f$ and $(x - \alpha) \mid f'$ $\Leftrightarrow$ $(x - \alpha) \mid \gcd(f, f')$ $\Leftrightarrow$ $\gcd(f, f') \neq 1$. $\blacksquare$

!!! definition "Definition 0.11 (Discriminant)"
    Let $f(x) = a_n x^n + \cdots + a_0 \in \mathbb{F}[x]$ with roots $\alpha_1, \ldots, \alpha_n$ (in the algebraic closure). The **discriminant** of $f(x)$ is defined as

    $$
    \Delta(f) = a_n^{2n-2} \prod_{i < j} (\alpha_i - \alpha_j)^2.
    $$

    $f$ has a multiple root if and only if $\Delta(f) = 0$.

For the quadratic $f(x) = ax^2 + bx + c$:

$$
\Delta = b^2 - 4ac.
$$

!!! example "Example 0.10"
    Determine whether $f(x) = x^3 - 3x + 2$ has multiple roots.

    $f'(x) = 3x^2 - 3 = 3(x-1)(x+1)$.

    $\gcd(f, f') = \gcd(x^3 - 3x + 2,\, 3x^2 - 3)$. By the Euclidean algorithm:

    $$
    x^3 - 3x + 2 = \tfrac{1}{3}x \cdot (3x^2 - 3) + (-2x + 2),
    $$

    $$
    3x^2 - 3 = \left(-\tfrac{3}{2}x - \tfrac{3}{2}\right)(-2x + 2) + 0.
    $$

    $\gcd(f, f') = x - 1 \neq 1$, so $f$ has a multiple root at $x = 1$. Check: $f(x) = (x-1)^2(x+2)$.

!!! example "Example 0.11"
    For the depressed cubic $f(x) = x^3 + px + q$, the discriminant is

    $$
    \Delta = -4p^3 - 27q^2.
    $$

    - $\Delta > 0$: three distinct real roots.
    - $\Delta = 0$: a multiple root.
    - $\Delta < 0$: one real root and two conjugate complex roots.

---

## 0.8 Irreducible Polynomials over $\mathbb{R}$ and $\mathbb{C}$

<div class="context-flow" markdown>

**Core question**: Which polynomials are irreducible over $\mathbb{C}$ and $\mathbb{R}$? → The fundamental theorem of algebra → Over $\mathbb{C}$, only degree-1 polynomials are irreducible → Over $\mathbb{R}$, only degree-1 and certain degree-2 polynomials are irreducible

</div>

!!! theorem "Theorem 0.16 (Fundamental theorem of algebra)"
    Every polynomial in $\mathbb{C}[x]$ of degree $\geq 1$ has at least one complex root.

!!! note "Note"
    The proof of the fundamental theorem of algebra requires tools from analysis (e.g., Liouville's theorem, the maximum modulus principle from complex analysis, or the fundamental group argument from topology) and lies outside the scope of pure algebra. We accept its conclusion here.

!!! theorem "Theorem 0.17 (Irreducible polynomials over $\\mathbb{C}$)"
    The irreducible polynomials in $\mathbb{C}[x]$ are precisely the polynomials of degree $1$. Consequently, every $f(x) \in \mathbb{C}[x]$ with $\deg f = n \geq 1$ factors as

    $$
    f(x) = a_n(x - \alpha_1)(x - \alpha_2) \cdots (x - \alpha_n),
    $$

    where $\alpha_1, \ldots, \alpha_n \in \mathbb{C}$ are all the roots of $f$ (counted with multiplicity) and $a_n$ is the leading coefficient.

??? proof "Proof"
    Suppose $p(x) \in \mathbb{C}[x]$ is irreducible with $\deg p \geq 2$. By the fundamental theorem of algebra, $p$ has a root $\alpha \in \mathbb{C}$, so $(x - \alpha) \mid p$. This contradicts the irreducibility of $p$. Hence every irreducible polynomial has degree exactly $1$.

    The factorization follows from the unique factorization theorem and the fact that every irreducible factor is linear. $\blacksquare$

!!! theorem "Theorem 0.18 (Conjugate root theorem)"
    Let $f(x) \in \mathbb{R}[x]$. If $\alpha = a + bi \in \mathbb{C}$ ($b \neq 0$) is a root of $f$, then $\overline{\alpha} = a - bi$ is also a root of $f$, with the same multiplicity.

??? proof "Proof"
    From $f(\alpha) = 0$, take complex conjugates: $\overline{f(\alpha)} = f(\overline{\alpha}) = 0$ (using the fact that the coefficients of $f$ are real, so $\overline{a_k \alpha^k} = a_k \overline{\alpha}^k$).

    For the multiplicity, if $f(x) = (x - \alpha)^k g(x)$ with $g(\alpha) \neq 0$, conjugation gives $f(x) = (x - \overline{\alpha})^k \overline{g}(x)$ where $\overline{g}$ denotes the polynomial with conjugated coefficients. Since the coefficients of $f$ are real, $\overline{g}$ also has real coefficients, and $\overline{g}(\overline{\alpha}) = \overline{g(\alpha)} \neq 0$. $\blacksquare$

!!! theorem "Theorem 0.19 (Irreducible polynomials over $\\mathbb{R}$)"
    The irreducible polynomials in $\mathbb{R}[x]$ are precisely:

    1. Linear polynomials $ax + b$ ($a \neq 0$);
    2. Quadratic polynomials $ax^2 + bx + c$ ($a \neq 0$) with negative discriminant $b^2 - 4ac < 0$.

??? proof "Proof"
    Linear polynomials are clearly irreducible. A quadratic with negative discriminant has no real roots, hence no linear factors, and is therefore irreducible.

    Suppose $p(x) \in \mathbb{R}[x]$ is irreducible with $\deg p \geq 3$. Then $p$ has a root $\alpha$ in $\mathbb{C}$. If $\alpha \in \mathbb{R}$, then $(x - \alpha) \mid p$, contradicting irreducibility. If $\alpha = a + bi$ ($b \neq 0$), then $\overline{\alpha}$ is also a root, and $(x - \alpha)(x - \overline{\alpha}) = x^2 - 2ax + a^2 + b^2 \in \mathbb{R}[x]$ divides $p$. Since $\deg p \geq 3$, this gives a proper factorization, a contradiction. Hence $\deg p \leq 2$.

    For $\deg p = 2$: if the discriminant is $\geq 0$, then $p$ has a real root and is reducible. $\blacksquare$

!!! example "Example 0.12"
    Factor $f(x) = x^4 + 4$ in $\mathbb{R}[x]$.

    $f$ has no real roots ($x^4 + 4 > 0$), but $f$ is nevertheless reducible:

    $$
    x^4 + 4 = x^4 + 4x^2 + 4 - 4x^2 = (x^2 + 2)^2 - (2x)^2 = (x^2 + 2x + 2)(x^2 - 2x + 2).
    $$

    The discriminants of $x^2 + 2x + 2$ and $x^2 - 2x + 2$ are both $4 - 8 = -4 < 0$, so each is irreducible.

---

## 0.9 Irreducible Polynomials over $\mathbb{Q}$

<div class="context-flow" markdown>

**Core question**: Irreducibility over $\mathbb{Q}$ is much harder to determine than over $\mathbb{R}$ or $\mathbb{C}$ → Rational root test → Gauss's lemma → Eisenstein's criterion

</div>

!!! theorem "Theorem 0.20 (Rational root theorem)"
    Let $f(x) = a_n x^n + \cdots + a_1 x + a_0 \in \mathbb{Z}[x]$, $a_n \neq 0$. If $\frac{p}{q}$ (with $\gcd(p,q) = 1$ and $q > 0$) is a rational root of $f$, then $p \mid a_0$ and $q \mid a_n$.

??? proof "Proof"
    From $f(p/q) = 0$, multiply both sides by $q^n$:

    $$
    a_n p^n + a_{n-1} p^{n-1} q + \cdots + a_1 p q^{n-1} + a_0 q^n = 0.
    $$

    We have $a_n p^n = -q(a_{n-1}p^{n-1} + \cdots + a_0 q^{n-1})$, so $q \mid a_n p^n$. Since $\gcd(p,q) = 1$, it follows that $q \mid a_n$. Similarly, $p \mid a_0$. $\blacksquare$

!!! definition "Definition 0.12 (Primitive polynomial)"
    An integer-coefficient polynomial $f(x) = a_n x^n + \cdots + a_0 \in \mathbb{Z}[x]$ is called **primitive** if $\gcd(a_n, a_{n-1}, \ldots, a_0) = 1$.

!!! theorem "Theorem 0.21 (Gauss's lemma)"
    The product of two primitive polynomials is primitive.

??? proof "Proof"
    Let $f = \sum a_i x^i$ and $g = \sum b_j x^j$ be primitive, and let $h = fg = \sum c_k x^k$. Suppose for contradiction that $h$ is not primitive. Then there exists a prime $p$ such that $p \mid c_k$ for all $k$.

    Since $f$ is primitive, let $a_r$ be the first coefficient not divisible by $p$ ($p \mid a_0, \ldots, a_{r-1}$, $p \nmid a_r$). Similarly, let $b_s$ be the first coefficient of $g$ not divisible by $p$. Consider $c_{r+s}$:

    $$
    c_{r+s} = \sum_{i+j=r+s} a_i b_j = a_r b_s + \sum_{\substack{i+j=r+s \\ i \neq r}} a_i b_j.
    $$

    When $i < r$, $p \mid a_i$; when $i > r$, $j = r+s-i < s$, so $p \mid b_j$. Thus $p \mid (c_{r+s} - a_r b_s)$, and since $p \mid c_{r+s}$, we get $p \mid a_r b_s$. Since $p$ is prime, $p \mid a_r$ or $p \mid b_s$, a contradiction. $\blacksquare$

!!! theorem "Theorem 0.22 (Gauss's theorem: reducibility over $\\mathbb{Q}$ implies reducibility over $\\mathbb{Z}$)"
    Let $f(x) \in \mathbb{Z}[x]$ be a primitive polynomial. If $f$ is reducible in $\mathbb{Q}[x]$, then $f$ can be written as the product of two integer-coefficient polynomials of lower degree.

??? proof "Proof"
    Suppose $f = gh$ with $g, h \in \mathbb{Q}[x]$, $\deg g, \deg h \geq 1$. Clearing denominators, write $\frac{a}{b} f = g_0 h_0$ with $g_0, h_0 \in \mathbb{Z}[x]$. Extract the content: $g_0 = d_1 g_1$, $h_0 = d_2 h_1$ with $g_1, h_1$ primitive. Then $\frac{a}{b} f = d_1 d_2 g_1 h_1$. By Gauss's lemma, $g_1 h_1$ is primitive, and since $f$ is also primitive, comparing contents gives $\frac{a}{b} = \pm d_1 d_2$, hence $f = \pm g_1 h_1$. $\blacksquare$

!!! theorem "Theorem 0.23 (Eisenstein's criterion)"
    Let $f(x) = a_n x^n + a_{n-1}x^{n-1} + \cdots + a_0 \in \mathbb{Z}[x]$, $n \geq 1$. If there exists a prime $p$ such that:

    1. $p \nmid a_n$;
    2. $p \mid a_{n-1}, p \mid a_{n-2}, \ldots, p \mid a_0$;
    3. $p^2 \nmid a_0$,

    then $f(x)$ is irreducible in $\mathbb{Q}[x]$.

??? proof "Proof"
    Suppose for contradiction that $f$ is reducible in $\mathbb{Q}[x]$. By Gauss's theorem, $f = gh$ with $g, h \in \mathbb{Z}[x]$, $\deg g = r \geq 1$, $\deg h = s \geq 1$, $r + s = n$. Write $g = b_r x^r + \cdots + b_0$ and $h = c_s x^s + \cdots + c_0$.

    From $a_0 = b_0 c_0$ and conditions (2)(3): $p \mid b_0 c_0$ but $p^2 \nmid b_0 c_0$, so $p$ divides exactly one of $b_0, c_0$. WLOG, $p \mid b_0$ and $p \nmid c_0$.

    From $a_n = b_r c_s$ and condition (1): $p \nmid b_r$. Let $b_k$ be the first coefficient of $g$ not divisible by $p$ ($1 \leq k \leq r < n$). Consider $a_k = b_k c_0 + b_{k-1}c_1 + \cdots + b_0 c_k$. Since $p \mid a_k$ (condition (2), as $k < n$) and $p \mid b_0, \ldots, b_{k-1}$, we get $p \mid b_k c_0$, hence $p \mid b_k$ (since $p \nmid c_0$), a contradiction. $\blacksquare$

!!! example "Example 0.13"
    Show that $f(x) = x^4 + 3x^3 + 9x + 3$ is irreducible over $\mathbb{Q}$.

    Take $p = 3$: $3 \nmid 1$ (leading coefficient), $3 \mid 3, 0, 9, 3$ (remaining coefficients), and $9 \nmid 3$ ($p^2 = 9 \nmid a_0 = 3$). By Eisenstein's criterion, $f$ is irreducible over $\mathbb{Q}$.

!!! example "Example 0.14"
    Show that the $p$-th cyclotomic polynomial $\Phi_p(x) = x^{p-1} + x^{p-2} + \cdots + x + 1$ (for $p$ prime) is irreducible over $\mathbb{Q}$.

    Substitute $y = x - 1$, i.e., $x = y + 1$:

    $$
    \Phi_p(x) = \frac{x^p - 1}{x - 1}, \quad \Phi_p(y + 1) = \frac{(y+1)^p - 1}{y} = \sum_{k=1}^{p} \binom{p}{k} y^{k-1} = y^{p-1} + \binom{p}{1} y^{p-2} + \cdots + \binom{p}{p-1}.
    $$

    Since $p \mid \binom{p}{k}$ for $1 \leq k \leq p-1$ and $p^2 \nmid \binom{p}{1} = p$, Eisenstein's criterion (with the prime $p$) shows that $\Phi_p(y+1)$ is irreducible over $\mathbb{Q}$, hence so is $\Phi_p(x)$.

---

## 0.10 Polynomial Interpolation

<div class="context-flow" markdown>

**Core question**: Given $n + 1$ points, can we uniquely determine a polynomial of degree $\leq n$ passing through them? → Lagrange interpolation → Newton interpolation → Applications

</div>

!!! theorem "Theorem 0.24 (Existence and uniqueness of polynomial interpolation)"
    Let $x_0, x_1, \ldots, x_n \in \mathbb{F}$ be pairwise distinct, and let $y_0, y_1, \ldots, y_n \in \mathbb{F}$ be given values. Then there exists a unique polynomial $p(x) \in \mathbb{F}[x]$ with $\deg p \leq n$ such that

    $$
    p(x_i) = y_i, \quad i = 0, 1, \ldots, n.
    $$

??? proof "Proof"
    **Uniqueness**: If $p$ and $q$ both satisfy the conditions with $\deg p, \deg q \leq n$, then $p - q$ has degree $\leq n$ and vanishes at the $n + 1$ distinct points $x_0, \ldots, x_n$. By Theorem 0.13, $p - q = 0$.

    **Existence**: Given by the Lagrange interpolation formula below. $\blacksquare$

!!! definition "Definition 0.13 (Lagrange basis polynomials)"
    For given distinct nodes $x_0, x_1, \ldots, x_n \in \mathbb{F}$, the **Lagrange basis polynomials** are defined by

    $$
    \ell_i(x) = \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}, \quad i = 0, 1, \ldots, n.
    $$

    They satisfy $\ell_i(x_j) = \delta_{ij}$ (the Kronecker delta).

!!! theorem "Theorem 0.25 (Lagrange interpolation formula)"
    The unique polynomial satisfying $p(x_i) = y_i$ ($i = 0, 1, \ldots, n$) is

    $$
    p(x) = \sum_{i=0}^{n} y_i \ell_i(x) = \sum_{i=0}^{n} y_i \prod_{\substack{j=0 \\ j \neq i}}^{n} \frac{x - x_j}{x_i - x_j}.
    $$

??? proof "Proof"
    From $\ell_i(x_k) = \delta_{ik}$, we directly verify $p(x_k) = \sum_i y_i \delta_{ik} = y_k$. That $\deg p \leq n$ is clear from the construction. $\blacksquare$

!!! example "Example 0.15"
    Find the quadratic interpolating polynomial through the points $(0, 1)$, $(1, 3)$, $(2, 7)$.

    Nodes: $x_0 = 0, x_1 = 1, x_2 = 2$; values: $y_0 = 1, y_1 = 3, y_2 = 7$.

    $$
    \ell_0(x) = \frac{(x-1)(x-2)}{(0-1)(0-2)} = \frac{(x-1)(x-2)}{2},
    $$

    $$
    \ell_1(x) = \frac{(x-0)(x-2)}{(1-0)(1-2)} = -x(x-2),
    $$

    $$
    \ell_2(x) = \frac{(x-0)(x-1)}{(2-0)(2-1)} = \frac{x(x-1)}{2}.
    $$

    $$
    p(x) = 1 \cdot \frac{(x-1)(x-2)}{2} + 3 \cdot [-x(x-2)] + 7 \cdot \frac{x(x-1)}{2}.
    $$

    Expanding and simplifying:

    $$
    p(x) = \frac{x^2 - 3x + 2}{2} - 3x^2 + 6x + \frac{7x^2 - 7x}{2} = \frac{2x^2 + 2x + 2}{2} = x^2 + x + 1.
    $$

    Check: $p(0) = 1$, $p(1) = 3$, $p(2) = 7$.

### Newton Interpolation

!!! definition "Definition 0.14 (Divided differences)"
    Let $x_0, x_1, \ldots, x_n$ be pairwise distinct. **Divided differences** are defined recursively:

    - Zeroth-order: $f[x_i] = f(x_i)$.
    - First-order: $f[x_i, x_{i+1}] = \dfrac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i}$.
    - $k$-th order: $f[x_i, x_{i+1}, \ldots, x_{i+k}] = \dfrac{f[x_{i+1}, \ldots, x_{i+k}] - f[x_i, \ldots, x_{i+k-1}]}{x_{i+k} - x_i}$.

!!! theorem "Theorem 0.26 (Newton interpolation formula)"
    The interpolating polynomial can be written in Newton form:

    $$
    p(x) = f[x_0] + f[x_0, x_1](x - x_0) + f[x_0, x_1, x_2](x - x_0)(x - x_1) + \cdots + f[x_0, \ldots, x_n]\prod_{i=0}^{n-1}(x - x_i).
    $$

??? proof "Proof"
    Let $p_k(x)$ be the interpolating polynomial through the first $k + 1$ points. Then $p_0(x) = f[x_0]$, and

    $$
    p_k(x) = p_{k-1}(x) + c_k \prod_{i=0}^{k-1}(x - x_i),
    $$

    where $c_k$ is determined by the condition $p_k(x_k) = y_k$. By induction, $c_k = f[x_0, x_1, \ldots, x_k]$. $\blacksquare$

!!! example "Example 0.16"
    Redo Example 0.15 using Newton interpolation (nodes $(0,1), (1,3), (2,7)$).

    Divided difference table:

    | $x_i$ | $f[x_i]$ | $f[x_i, x_{i+1}]$ | $f[x_i, x_{i+1}, x_{i+2}]$ |
    |---|---|---|---|
    | $0$ | $1$ | | |
    | | | $\frac{3-1}{1-0} = 2$ | |
    | $1$ | $3$ | | $\frac{4-2}{2-0} = 1$ |
    | | | $\frac{7-3}{2-1} = 4$ | |
    | $2$ | $7$ | | |

    $$
    p(x) = 1 + 2(x - 0) + 1 \cdot (x - 0)(x - 1) = 1 + 2x + x^2 - x = x^2 + x + 1.
    $$

    This agrees with the Lagrange result.

!!! note "Note"
    The advantage of Newton interpolation is its **recursive nature**: when a new node is added, only one additional term is needed, without recomputing everything. Moreover, divided differences are symmetric: $f[x_0, x_1, \ldots, x_k]$ does not depend on the ordering of the nodes.
