# Chapter 9  Quadratic Forms

<div class="context-flow" markdown>

**Prerequisites**: Ch8 Spectral theorem for symmetric matrices

**Chapter arc**: $\mathbf{x}^TA\mathbf{x}$ → Completing the square / Orthogonal diagonalization to canonical form → **Law of inertia** (signature invariance) → Positive definiteness criteria → Geometry (ellipsoids/hyperboloids)

**Further connections**：Quadratic forms pervade optimization (quadratic programming), statistics ($\chi^2$ distribution), and differential geometry (Riemannian metric); the antisymmetric case of bilinear forms yields symplectic forms — symplectic geometry is the mathematical language of classical mechanics (Hamilton's equations) and quantum field theory; Hermitian forms provide the mathematical structure for observables in quantum mechanics

</div>

A quadratic form is the algebraic theory of homogeneous quadratic polynomials, which has deep connections with symmetric matrices and inner product spaces. The study of quadratic forms is not only an important part of linear algebra but also has broad applications in differential geometry, optimization theory, statistics, and physics. This chapter systematically studies the definition, simplification methods, the law of inertia, and positive definiteness criteria for quadratic forms.

---

## 9.1 Definition of Quadratic Forms

<div class="context-flow" markdown>

**Symmetric matrix** $A$ $\leftrightarrow$ quadratic form $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ is a one-to-one correspondence; cross terms $x_ix_j$ are split as $a_{ij} = a_{ji}$

</div>

!!! definition "Definition 9.1 (Quadratic form)"
    Let $\mathbb{F} = \mathbb{R}$ (or $\mathbb{C}$). A **quadratic form** in $n$ variables $x_1, x_2, \ldots, x_n$ is a homogeneous quadratic polynomial of the form:

    $$Q(x_1, x_2, \ldots, x_n) = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i x_j$$

    where $a_{ij} \in \mathbb{F}$. In vector notation, let $\mathbf{x} = (x_1, \ldots, x_n)^T$, then

    $$Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$$

    where $A = (a_{ij})_{n \times n}$.

!!! definition "Definition 9.2 (Matrix of a quadratic form)"
    For the quadratic form $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$, we can always assume the matrix $A$ is **symmetric**. Indeed, for any matrix $A$, $\mathbf{x}^T A \mathbf{x} = \mathbf{x}^T \left(\frac{A + A^T}{2}\right) \mathbf{x}$, and $\frac{A + A^T}{2}$ is symmetric.

    The symmetric matrix $A$ is called the **matrix** of the quadratic form $Q$, and $\operatorname{rank}(A)$ is called the **rank** of $Q$.

!!! theorem "Theorem 9.1 (One-to-one correspondence between quadratic forms and symmetric matrices)"
    There is a one-to-one correspondence between real quadratic forms $Q(\mathbf{x})$ in $n$ variables and $n \times n$ real symmetric matrices $A$: $Q(\mathbf{x}) = \mathbf{x}^T A \mathbf{x}$.

??? proof "Proof"
    Given a symmetric matrix $A = (a_{ij})$, $\mathbf{x}^T A \mathbf{x} = \sum_{i,j} a_{ij}x_ix_j$ is a quadratic form.

    Conversely, given a quadratic form $Q(\mathbf{x}) = \sum_{i \leq j} c_{ij}x_ix_j$ (where $c_{ii}$ is the coefficient of $x_i^2$ and $c_{ij}$ ($i < j$) is the coefficient of $x_ix_j$), define the symmetric matrix $A$ with entries $a_{ii} = c_{ii}$, $a_{ij} = a_{ji} = c_{ij}/2$ ($i < j$). Then $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$.

    Uniqueness: If $\mathbf{x}^TA\mathbf{x} = \mathbf{x}^TB\mathbf{x}$ for all $\mathbf{x}$, with $A, B$ both symmetric, then $\mathbf{x}^T(A-B)\mathbf{x} = 0$ for all $\mathbf{x}$. Taking $\mathbf{x} = \mathbf{e}_i$ gives $a_{ii} = b_{ii}$; taking $\mathbf{x} = \mathbf{e}_i + \mathbf{e}_j$ gives $a_{ij} + a_{ji} = b_{ij} + b_{ji}$, and by symmetry $a_{ij} = b_{ij}$. $\blacksquare$

!!! example "Example 9.1"
    The symmetric matrix of the quadratic form $Q(x_1, x_2, x_3) = 2x_1^2 + 3x_2^2 - x_3^2 + 4x_1x_2 - 6x_1x_3 + 2x_2x_3$ is

    $$A = \begin{pmatrix} 2 & 2 & -3 \\ 2 & 3 & 1 \\ -3 & 1 & -1 \end{pmatrix}$$

    Note that the cross term $4x_1x_2$ is split as $a_{12} = a_{21} = 2$.

!!! example "Example 9.2"
    The quadratic form $Q(x_1, x_2) = x_1^2 + x_2^2$ on $\mathbb{R}^2$ corresponds to the matrix $A = I_2$, geometrically representing circles centered at the origin. The quadratic form $Q(x_1, x_2) = x_1^2 - x_2^2$ corresponds to the matrix $A = \operatorname{diag}(1, -1)$, geometrically representing a hyperbola.

---

## 9.2 Canonical Form of Quadratic Forms

<div class="context-flow" markdown>

Eliminate cross terms → Only $d_i y_i^2$ remains → **Completing the square** (Lagrange) is a constructive tool; any quadratic form can be reduced to canonical form

</div>

!!! definition "Definition 9.3 (Canonical form)"
    If a quadratic form $Q(\mathbf{x})$ contains only squared terms (no cross terms), i.e.,

    $$Q(\mathbf{x}) = d_1 x_1^2 + d_2 x_2^2 + \cdots + d_n x_n^2$$

    then $Q$ is said to be in **canonical form** (or diagonal form), with matrix $\operatorname{diag}(d_1, \ldots, d_n)$.

!!! definition "Definition 9.4 (Nonsingular linear substitution)"
    Let $\mathbf{x} = C\mathbf{y}$, where $C$ is an $n \times n$ invertible matrix. This is called a **nonsingular linear substitution**. Under this substitution,

    $$Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} = (C\mathbf{y})^TA(C\mathbf{y}) = \mathbf{y}^T(C^TAC)\mathbf{y}$$

    The matrix of the new quadratic form is $B = C^TAC$.

### Completing the Square

!!! theorem "Theorem 9.2 (Reduction to canonical form by completing the square)"
    Any real quadratic form $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ ($A$ symmetric) can be reduced to canonical form by a nonsingular linear substitution.

??? proof "Proof"
    By completing the square (Lagrange's method). Two cases:

    **Case 1:** If some $a_{ii} \neq 0$ (say $a_{11} \neq 0$), complete the square in all terms containing $x_1$:

    $$Q = a_{11}\left(x_1 + \frac{a_{12}}{a_{11}}x_2 + \cdots + \frac{a_{1n}}{a_{11}}x_n\right)^2 + Q_1(x_2, \ldots, x_n)$$

    Let $y_1 = x_1 + \frac{a_{12}}{a_{11}}x_2 + \cdots + \frac{a_{1n}}{a_{11}}x_n$, $y_i = x_i$ ($i \geq 2$); this is a nonsingular linear substitution. Apply the method recursively to $Q_1$.

    **Case 2:** If all $a_{ii} = 0$ but some $a_{ij} \neq 0$ ($i \neq j$). Let $x_i = y_i + y_j$, $x_j = y_i - y_j$, and $x_k = y_k$ for all other $k$. Then $2a_{ij}x_ix_j = 2a_{ij}(y_i^2 - y_j^2)$, producing squared terms, reducing to Case 1.

    By induction, after finitely many steps, $Q$ is reduced to canonical form. $\blacksquare$

!!! example "Example 9.3"
    Reduce to canonical form by completing the square: $Q(x_1, x_2, x_3) = 2x_1x_2 + 2x_1x_3 - 6x_2x_3$.

    All squared-term coefficients are zero (Case 2). Let $x_1 = y_1 + y_2$, $x_2 = y_1 - y_2$, $x_3 = y_3$:

    $$Q = 2(y_1+y_2)(y_1-y_2) + 2(y_1+y_2)y_3 - 6(y_1-y_2)y_3$$

    $$= 2y_1^2 - 2y_2^2 + 2y_1y_3 + 2y_2y_3 - 6y_1y_3 + 6y_2y_3$$

    $$= 2y_1^2 - 4y_1y_3 - 2y_2^2 + 8y_2y_3$$

    Completing the square in $y_1$: $2(y_1^2 - 2y_1y_3) = 2(y_1 - y_3)^2 - 2y_3^2$.

    Completing the square in $y_2$: $-2(y_2^2 - 4y_2y_3) = -2(y_2 - 2y_3)^2 + 8y_3^2$.

    $$Q = 2(y_1 - y_3)^2 - 2(y_2 - 2y_3)^2 + 6y_3^2$$

    Let $z_1 = y_1 - y_3$, $z_2 = y_2 - 2y_3$, $z_3 = y_3$, yielding the canonical form $Q = 2z_1^2 - 2z_2^2 + 6z_3^2$.

---

## 9.3 Orthogonal Diagonalization to Canonical Form

<div class="context-flow" markdown>

The substitution matrix from completing the square is not unique → Use the Ch8 **spectral theorem** $A = Q\Lambda Q^T$ for an orthogonal substitution → Canonical form coefficients = eigenvalues, transformation = isometry

</div>

The canonical form obtained by completing the square depends on the order of completion, and the transformation matrix is not unique. The orthogonal diagonalization method (using the spectral theorem) provides the most "natural" approach to canonical form reduction.

!!! theorem "Theorem 9.3 (Orthogonal diagonalization method)"
    Let $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$, where $A$ is an $n \times n$ real symmetric matrix. By the spectral theorem, there exists an orthogonal matrix $P$ such that

    $$P^TAP = \Lambda = \operatorname{diag}(\lambda_1, \lambda_2, \ldots, \lambda_n)$$

    where $\lambda_1, \ldots, \lambda_n$ are the eigenvalues of $A$. Let $\mathbf{x} = P\mathbf{y}$ (orthogonal substitution), then

    $$Q = \mathbf{y}^T\Lambda\mathbf{y} = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

??? proof "Proof"
    By the spectral theorem for real symmetric matrices (Theorem 8.15), there exists an orthogonal matrix $P$ (whose columns are orthonormal eigenvectors of $A$) such that $P^TAP = \Lambda$. Under the substitution $\mathbf{x} = P\mathbf{y}$:

    $$Q(\mathbf{x}) = (P\mathbf{y})^T A (P\mathbf{y}) = \mathbf{y}^T(P^TAP)\mathbf{y} = \mathbf{y}^T\Lambda\mathbf{y} = \sum_{i=1}^n \lambda_i y_i^2$$

    $\blacksquare$

!!! note "Note"
    The advantage of orthogonal substitution is that it preserves vector lengths and angles (it is an isometric transformation), making it particularly meaningful in geometric applications. The coefficients in the canonical form obtained by orthogonal diagonalization are exactly the eigenvalues, which is very natural from a theoretical standpoint.

!!! example "Example 9.4"
    Reduce to canonical form by orthogonal substitution: $Q(x_1, x_2) = 5x_1^2 + 4x_1x_2 + 8x_2^2$.

    Symmetric matrix $A = \begin{pmatrix} 5 & 2 \\ 2 & 8 \end{pmatrix}$.

    Characteristic polynomial: $\det(A - \lambda I) = (5-\lambda)(8-\lambda) - 4 = \lambda^2 - 13\lambda + 36 = (\lambda - 4)(\lambda - 9)$.

    $\lambda_1 = 4$: $(A - 4I)\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\begin{pmatrix} 1 & 2 \\ 2 & 4 \end{pmatrix}\mathbf{x} = \mathbf{0}$, $\mathbf{v}_1 = \frac{1}{\sqrt{5}}(-2, 1)^T$.

    $\lambda_2 = 9$: $(A - 9I)\mathbf{x} = \mathbf{0}$ $\Rightarrow$ $\begin{pmatrix} -4 & 2 \\ 2 & -1 \end{pmatrix}\mathbf{x} = \mathbf{0}$, $\mathbf{v}_2 = \frac{1}{\sqrt{5}}(1, 2)^T$.

    Orthogonal matrix $P = \frac{1}{\sqrt{5}}\begin{pmatrix} -2 & 1 \\ 1 & 2 \end{pmatrix}$, let $\mathbf{x} = P\mathbf{y}$, giving $Q = 4y_1^2 + 9y_2^2$.

---

## 9.4 Law of Inertia

<div class="context-flow" markdown>

The coefficients in the canonical form may vary, but the **number of positive coefficients $p$ and negative coefficients $q$ are invariant** — Sylvester's law of inertia is the core invariant of quadratic form theory

</div>

!!! definition "Definition 9.5 (Index of inertia)"
    Suppose the real quadratic form $Q(\mathbf{x})$ is reduced to the canonical form

    $$Q = d_1 y_1^2 + d_2 y_2^2 + \cdots + d_r y_r^2$$

    by a nonsingular linear substitution, where $d_i \neq 0$ ($i = 1, \ldots, r$) and $r = \operatorname{rank}(A)$. Let $p$ be the number of positive coefficients and $q = r - p$ the number of negative coefficients. Then

    - $p$ is called the **positive index of inertia**;
    - $q$ is called the **negative index of inertia**;
    - $(p, q)$ is called the **signature** of the quadratic form.

<div class="context-flow" markdown>

**Insight**: The core of the proof is a **dimension argument** — $V_1 \cap V_2 \neq \{0\}$ (since $\dim V_1 + \dim V_2 > n$) leads to a contradiction; this technique reappears in the Eckart-Young theorem in Ch11

</div>

!!! theorem "Theorem 9.4 (Sylvester's law of inertia)"
    The number of positive coefficients $p$ and negative coefficients $q$ in the canonical form of a real quadratic form are invariant, independent of the nonsingular linear substitution used. That is, $p$ and $q$ are determined solely by the quadratic form itself.

??? proof "Proof"
    Suppose the quadratic form $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ is reduced to canonical form by two different nonsingular linear substitutions $\mathbf{x} = C_1\mathbf{y}$ and $\mathbf{x} = C_2\mathbf{z}$:

    $$Q = y_1^2 + \cdots + y_p^2 - y_{p+1}^2 - \cdots - y_r^2$$

    $$Q = z_1^2 + \cdots + z_s^2 - z_{s+1}^2 - \cdots - z_r^2$$

    (Without loss of generality, the coefficients can be normalized to $\pm 1$.)

    By contradiction, assume $p \neq s$, say $p > s$.

    Let $\mathbf{y} = C_1^{-1}\mathbf{x}$, $\mathbf{z} = C_2^{-1}\mathbf{x}$. Consider the two subspaces:

    - $V_1 = \{\mathbf{x} \in \mathbb{R}^n : y_{p+1} = \cdots = y_n = 0\}$, $\dim V_1 = p$;
    - $V_2 = \{\mathbf{x} \in \mathbb{R}^n : z_1 = \cdots = z_s = 0\}$, $\dim V_2 = n - s$.

    Since $\dim V_1 + \dim V_2 = p + (n-s) > n$ (because $p > s$), we have $V_1 \cap V_2 \neq \{\mathbf{0}\}$.

    Let $\mathbf{0} \neq \mathbf{x}_0 \in V_1 \cap V_2$.

    - In $V_1$: $Q(\mathbf{x}_0) = y_1^2 + \cdots + y_p^2 > 0$ (since $\mathbf{x}_0 \neq \mathbf{0}$ implies at least one $y_i \neq 0$, $i \leq p$);
    - In $V_2$: $Q(\mathbf{x}_0) = -z_{s+1}^2 - \cdots - z_r^2 \leq 0$.

    Contradiction! Therefore $p = s$. $\blacksquare$

!!! corollary "Corollary 9.1"
    Two real quadratic forms are equivalent (i.e., can be transformed into each other by a nonsingular linear substitution) if and only if they have the same rank and the same positive index of inertia (or equivalently, the same signature).

!!! proposition "Proposition 9.1"
    The positive index of inertia of a real symmetric matrix $A$ equals the number of positive eigenvalues (counting multiplicity), and the negative index of inertia equals the number of negative eigenvalues (counting multiplicity).

??? proof "Proof"
    By orthogonal diagonalization, $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ can be reduced to $\lambda_1 y_1^2 + \cdots + \lambda_n y_n^2$ by an orthogonal substitution, where $\lambda_i$ are the eigenvalues of $A$. By the law of inertia, the positive index of inertia equals the number of positive eigenvalues, and the negative index of inertia equals the number of negative eigenvalues. $\blacksquare$

!!! example "Example 9.5"
    The matrix of the quadratic form $Q = x_1^2 + 4x_1x_2 + 4x_2^2 + 2x_3^2$ is

    $$A = \begin{pmatrix} 1 & 2 & 0 \\ 2 & 4 & 0 \\ 0 & 0 & 2 \end{pmatrix}$$

    Characteristic polynomial: $\det(A - \lambda I) = (2-\lambda)[(1-\lambda)(4-\lambda) - 4] = (2-\lambda)(\lambda^2 - 5\lambda) = -\lambda(2-\lambda)(\lambda - 5)$.

    Eigenvalues: $\lambda_1 = 0, \lambda_2 = 2, \lambda_3 = 5$. Positive index of inertia $p = 2$, negative index of inertia $q = 0$, rank $r = 2$.

---

## 9.5 Congruence Transformations and Congruent Matrices

<div class="context-flow" markdown>

**Similarity** $P^{-1}AP$ (preserves eigenvalues) vs **Congruence** $C^TAC$ (preserves signature) — Congruence is the natural equivalence relation for quadratic forms; symmetry and rank are preserved but eigenvalues may change

</div>

!!! definition "Definition 9.6 (Congruence)"
    Let $A, B$ be $n \times n$ real matrices. If there exists an invertible matrix $C$ such that

    $$B = C^TAC$$

    then $A$ and $B$ are said to be **congruent**, denoted $A \simeq B$. The map $\mathbf{x} \mapsto C\mathbf{x}$ is called a **congruence transformation**.

!!! proposition "Proposition 9.2 (Congruence is an equivalence relation)"
    Matrix congruence is an equivalence relation:

    1. **Reflexivity**: $A \simeq A$ (take $C = I$);
    2. **Symmetry**: If $A \simeq B$, then $B \simeq A$;
    3. **Transitivity**: If $A \simeq B$ and $B \simeq D$, then $A \simeq D$.

??? proof "Proof"
    (2) If $B = C^TAC$, then $A = (C^{-1})^T B C^{-1} = (C^{-1})^T B (C^{-1})$, so $B \simeq A$.

    (3) If $B = C_1^TAC_1$ and $D = C_2^TBC_2$, then $D = C_2^T(C_1^TAC_1)C_2 = (C_1C_2)^T A (C_1C_2)$. $\blacksquare$

!!! theorem "Theorem 9.5 (Congruence canonical form)"
    Any $n \times n$ real symmetric matrix $A$ (of rank $r$) is congruent to

    $$\begin{pmatrix} I_p & & \\ & -I_q & \\ & & O_{n-r} \end{pmatrix}$$

    where $p$ is the positive index of inertia and $q = r - p$ is the negative index of inertia. This canonical form is unique by the law of inertia.

??? proof "Proof"
    By completing the square (Theorem 9.2), there exists an invertible matrix $C_1$ such that $C_1^TAC_1 = \operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0)$, where $d_i \neq 0$. Without loss of generality, assume $d_1, \ldots, d_p > 0$ and $d_{p+1}, \ldots, d_r < 0$. Let

    $$C_2 = \operatorname{diag}\left(\frac{1}{\sqrt{|d_1|}}, \ldots, \frac{1}{\sqrt{|d_r|}}, 1, \ldots, 1\right)$$

    Then $C_2^T(\operatorname{diag}(d_1, \ldots, d_r, 0, \ldots, 0))C_2 = \operatorname{diag}(1, \ldots, 1, -1, \ldots, -1, 0, \ldots, 0)$.

    Taking $C = C_1C_2$ gives the result. $\blacksquare$

!!! note "Note"
    Congruence preserves symmetry and rank: if $A$ is symmetric and $B = C^TAC$, then $B$ is also symmetric and $\operatorname{rank}(B) = \operatorname{rank}(A)$. However, congruence does **not** preserve eigenvalues. By contrast, similarity preserves eigenvalues but does not necessarily preserve symmetry.

!!! example "Example 9.6"
    The matrix $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$ has eigenvalues $3$ and $-1$, so the positive index of inertia is $p = 1$ and the negative index of inertia is $q = 1$. $A$ is congruent to $\begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$.

    Verification: Taking $C = \begin{pmatrix} 1 & -1 \\ 0 & 1 \end{pmatrix}$ (corresponding to the completed square $(x_1+2x_2)^2 - 3x_2^2$ with rescaling), one can explicitly compute the congruence transformation matrix.

---

## 9.6 Positive Definite Quadratic Forms and Positive Definite Matrices

<div class="context-flow" markdown>

Signature $(n,0)$ $\leftrightarrow$ all eigenvalues $> 0$ $\leftrightarrow$ $A = C^TC$ $\leftrightarrow$ all leading principal minors positive → Positive definite matrices lead directly to Ch10 **Cholesky decomposition** $A = LL^T$

</div>

Positive definiteness is one of the most important properties of quadratic forms and symmetric matrices, playing a central role in optimization, statistics, differential equations, and other fields.

!!! definition "Definition 9.7 (Definiteness classification)"
    Let $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ be a real quadratic form with $A$ an $n \times n$ real symmetric matrix. We say $Q$ (or $A$) is:

    - **Positive definite**: if $Q(\mathbf{x}) > 0$ for all $\mathbf{x} \neq \mathbf{0}$;
    - **Positive semidefinite**: if $Q(\mathbf{x}) \geq 0$ for all $\mathbf{x}$;
    - **Negative definite**: if $Q(\mathbf{x}) < 0$ for all $\mathbf{x} \neq \mathbf{0}$;
    - **Negative semidefinite**: if $Q(\mathbf{x}) \leq 0$ for all $\mathbf{x}$;
    - **Indefinite**: if $Q$ takes both positive and negative values.

<div class="context-flow" markdown>

**Insight**: The five equivalent conditions unify the algebraic (eigenvalues), geometric ($A = C^TC$), and combinatorial (leading principal minors) perspectives

</div>

!!! theorem "Theorem 9.6 (Equivalent conditions for positive definiteness)"
    Let $A$ be an $n \times n$ real symmetric matrix. The following conditions are equivalent:

    1. $A$ is positive definite;
    2. All eigenvalues $\lambda_1, \ldots, \lambda_n$ of $A$ are positive;
    3. The positive index of inertia is $p = n$ (i.e., the canonical form of $Q$ is $y_1^2 + \cdots + y_n^2$);
    4. There exists an invertible matrix $C$ such that $A = C^TC$;
    5. All **leading principal minors** of $A$ are positive.

??? proof "Proof"
    **(1)$\Leftrightarrow$(2):** If $A$ is positive definite and $A\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v} \neq \mathbf{0}$, then $\lambda\|\mathbf{v}\|^2 = \mathbf{v}^TA\mathbf{v} > 0$, so $\lambda > 0$. Conversely, if all eigenvalues are positive, $Q(\mathbf{x}) = \mathbf{y}^T\Lambda\mathbf{y} = \sum \lambda_i y_i^2 > 0$ ($\mathbf{x} \neq \mathbf{0}$).

    **(2)$\Leftrightarrow$(3):** The positive index of inertia equals the number of positive eigenvalues (Proposition 9.1).

    **(1)$\Rightarrow$(4):** By orthogonal diagonalization $A = P\Lambda P^T$, $\Lambda = \operatorname{diag}(\lambda_1, \ldots, \lambda_n)$, $\lambda_i > 0$. Let $C = \Lambda^{1/2}P^T$ (where $\Lambda^{1/2} = \operatorname{diag}(\sqrt{\lambda_1}, \ldots, \sqrt{\lambda_n})$), then $C^TC = P\Lambda^{1/2}\Lambda^{1/2}P^T = P\Lambda P^T = A$.

    **(4)$\Rightarrow$(1):** $\mathbf{x}^TA\mathbf{x} = \mathbf{x}^TC^TC\mathbf{x} = \|C\mathbf{x}\|^2 \geq 0$, with equality if and only if $C\mathbf{x} = \mathbf{0}$, i.e., $\mathbf{x} = \mathbf{0}$ (since $C$ is invertible).

    **(1)$\Leftrightarrow$(5):** This is the Sylvester criterion. Let $A_k$ be the $k$-th leading principal submatrix of $A$ and $\Delta_k = \det(A_k)$.

    If $A$ is positive definite, then $A_k$ is also positive definite (for $\mathbf{x} = (x_1, \ldots, x_k, 0, \ldots, 0)^T$, $\mathbf{x}^TA\mathbf{x} = \mathbf{y}^TA_k\mathbf{y}$ where $\mathbf{y} = (x_1, \ldots, x_k)^T$). $A_k$ positive definite $\Rightarrow$ all eigenvalues positive $\Rightarrow$ $\Delta_k = \prod \lambda_i^{(k)} > 0$.

    Conversely, prove by induction on $n$. For $n=1$, $\Delta_1 = a_{11} > 0$ is positive definite. Assume the result holds for $n-1$. By $\Delta_1, \ldots, \Delta_{n-1} > 0$, $A_{n-1}$ is positive definite. Using the Schur complement, one can show $A$ is positive definite. $\blacksquare$

!!! theorem "Theorem 9.7 (Equivalent conditions for positive semidefiniteness)"
    Let $A$ be an $n \times n$ real symmetric matrix. The following conditions are equivalent:

    1. $A$ is positive semidefinite;
    2. All eigenvalues $\lambda_i \geq 0$;
    3. There exists a matrix $B$ (not necessarily invertible) such that $A = B^TB$;
    4. All **principal minors** (not just leading principal minors) of $A$ are nonnegative.

??? proof "Proof"
    The proofs of (1)$\Leftrightarrow$(2) and (1)$\Leftrightarrow$(3) are similar to the positive definite case.

    (4) Necessity: $A$ positive semidefinite $\Rightarrow$ every principal submatrix is also positive semidefinite $\Rightarrow$ principal minors (= products of eigenvalues) $\geq 0$. $\blacksquare$

!!! note "Note"
    To test for negative definiteness: $A$ is negative definite if and only if $-A$ is positive definite, equivalently $(-1)^k\Delta_k > 0$ ($k = 1, \ldots, n$), i.e., odd-order leading principal minors are negative and even-order ones are positive.

!!! example "Example 9.7"
    Determine the definiteness of $A = \begin{pmatrix} 2 & -1 & 0 \\ -1 & 2 & -1 \\ 0 & -1 & 2 \end{pmatrix}$.

    $\Delta_1 = 2 > 0$, $\Delta_2 = \det\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix} = 3 > 0$, $\Delta_3 = \det(A) = 2(4-1) - (-1)(-2) = 6 - 2 = 4 > 0$.

    All leading principal minors are positive, so $A$ is positive definite.

!!! example "Example 9.8"
    Determine the definiteness of $A = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$.

    $\Delta_1 = 1 > 0$, $\Delta_2 = 1 - 4 = -3 < 0$.

    $A$ is not positive definite. In fact, the eigenvalues of $A$ are $3$ and $-1$, so $A$ is indefinite.

---

## 9.7 Geometric Meaning of Quadratic Forms

<div class="context-flow" markdown>

$\mathbf{x}^TA\mathbf{x} = c$ defines a quadric surface → Orthogonal substitution along **eigenvector directions** (principal axes) eliminates cross terms → Surface type is determined by the signature $(p,q)$

</div>

Quadratic forms geometrically describe quadratic curves and surfaces.

!!! definition "Definition 9.8 (Quadric surface)"
    The set defined by the equation $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x} = c$ ($c$ a constant) in $\mathbb{R}^n$ is called a **quadric surface**. In $\mathbb{R}^2$ it is a **conic section**.

!!! theorem "Theorem 9.8 (Principal axis theorem)"
    Let $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ be a real quadratic form with $A$ an $n \times n$ real symmetric matrix. Through the orthogonal substitution $\mathbf{x} = P\mathbf{y}$ ($P$ has columns that are orthonormal eigenvectors of $A$), the quadratic form reduces to canonical form

    $$Q = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \cdots + \lambda_n y_n^2$$

    The directions of the new coordinate axes $y_1, y_2, \ldots, y_n$ (i.e., the columns of $P$) are called the **principal axes** of the quadric surface.

??? proof "Proof"
    This follows directly from the spectral theorem and orthogonal substitution. The orthogonal substitution amounts to rotating the coordinate system so that the new axes align with the eigenvector directions. In the new coordinates, the quadratic form has no cross terms, giving the simplest form of the quadric surface equation. $\blacksquare$

!!! example "Example 9.9"
    Classify the conic $5x_1^2 + 4x_1x_2 + 8x_2^2 = 36$ in $\mathbb{R}^2$.

    From Example 9.4, after orthogonal substitution we get $4y_1^2 + 9y_2^2 = 36$, i.e., $\dfrac{y_1^2}{9} + \dfrac{y_2^2}{4} = 1$. This is an **ellipse** with principal axes along the $y_1$ and $y_2$ axes, semi-major axis $a = 3$ and semi-minor axis $b = 2$.

    The principal axis directions are the eigenvectors of $A$: $\mathbf{v}_1 = \frac{1}{\sqrt{5}}(-2, 1)^T$ (for $\lambda = 4$) and $\mathbf{v}_2 = \frac{1}{\sqrt{5}}(1, 2)^T$ (for $\lambda = 9$).

!!! example "Example 9.10"
    Standard classification of quadric surfaces in $\mathbb{R}^3$ (with $Q = \lambda_1 y_1^2 + \lambda_2 y_2^2 + \lambda_3 y_3^2 = 1$):

    - **Ellipsoid**: $\lambda_1, \lambda_2, \lambda_3 > 0$, e.g., $\frac{y_1^2}{a^2} + \frac{y_2^2}{b^2} + \frac{y_3^2}{c^2} = 1$;
    - **Hyperboloid of one sheet**: Two positive, one negative, e.g., $\frac{y_1^2}{a^2} + \frac{y_2^2}{b^2} - \frac{y_3^2}{c^2} = 1$;
    - **Hyperboloid of two sheets**: One positive, two negative, e.g., $\frac{y_1^2}{a^2} - \frac{y_2^2}{b^2} - \frac{y_3^2}{c^2} = 1$;
    - If $Q = 0$ (homogeneous case), it is a **quadric cone**.

    The type of quadric surface is determined by the **signature** $(p, q)$ of the quadratic form.

---

## 9.8 Bilinear Forms

<div class="context-flow" markdown>

**Prerequisites**: Quadratic form $Q(\mathbf{x}) = \mathbf{x}^TA\mathbf{x}$ is a "one-variable function" → Bilinear form $f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^TA\mathbf{y}$ is a "two-variable function" → Through the **polarization identity** a quadratic form recovers the bilinear form → Symmetric/antisymmetric bilinear forms lead to orthogonal and symplectic groups

</div>

A quadratic form is the "diagonal" value $Q(\mathbf{x}) = f(\mathbf{x}, \mathbf{x})$. To fully understand quadratic forms, one must first understand the more general bilinear forms. Bilinear forms are the most fundamental "scalar functions of two vectors" in linear algebra — inner products, determinants, and area forms are all special cases.

!!! definition "Definition 9.9 (Bilinear form)"
    Let $V$ be an $n$-dimensional vector space over a field $\mathbb{F}$. A map $f: V \times V \to \mathbb{F}$ is called a **bilinear form** on $V$ if $f$ is linear in each variable:

    - In the first variable: $f(\alpha \mathbf{x}_1 + \beta \mathbf{x}_2, \mathbf{y}) = \alpha f(\mathbf{x}_1, \mathbf{y}) + \beta f(\mathbf{x}_2, \mathbf{y})$
    - In the second variable: $f(\mathbf{x}, \alpha \mathbf{y}_1 + \beta \mathbf{y}_2) = \alpha f(\mathbf{x}, \mathbf{y}_1) + \beta f(\mathbf{x}, \mathbf{y}_2)$

    for all $\mathbf{x}, \mathbf{x}_1, \mathbf{x}_2, \mathbf{y}, \mathbf{y}_1, \mathbf{y}_2 \in V$ and $\alpha, \beta \in \mathbb{F}$.

!!! definition "Definition 9.10 (Matrix of a bilinear form)"
    Let $\mathcal{B} = \{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ be a basis of $V$. The **matrix** (also called the **Gram matrix**) of the bilinear form $f$ with respect to $\mathcal{B}$ is defined as

    $$
    A = (a_{ij})_{n \times n}, \quad a_{ij} = f(\mathbf{e}_i, \mathbf{e}_j)
    $$

    If $\mathbf{x} = \sum x_i \mathbf{e}_i$ and $\mathbf{y} = \sum y_j \mathbf{e}_j$, then

    $$
    f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^T A \mathbf{y} = \sum_{i=1}^n \sum_{j=1}^n a_{ij} x_i y_j
    $$

!!! theorem "Theorem 9.9 (Change of basis for bilinear forms)"
    Let $f$ have matrix $A$ with respect to basis $\mathcal{B}$ and matrix $A'$ with respect to basis $\mathcal{B}'$. If $C$ is the change-of-basis matrix from $\mathcal{B}$ to $\mathcal{B}'$, then

    $$
    A' = C^T A C
    $$

??? proof "Proof"
    Let $\mathcal{B}' = \{\mathbf{e}_1', \ldots, \mathbf{e}_n'\}$ with $(\mathbf{e}_1', \ldots, \mathbf{e}_n') = (\mathbf{e}_1, \ldots, \mathbf{e}_n)C$. Then

    $$
    a_{ij}' = f(\mathbf{e}_i', \mathbf{e}_j') = f\!\left(\sum_k c_{ki}\mathbf{e}_k, \sum_l c_{lj}\mathbf{e}_l\right) = \sum_{k,l} c_{ki} a_{kl} c_{lj}
    $$

    In matrix form, this is $A' = C^TAC$. $\blacksquare$

!!! definition "Definition 9.11 (Rank and nondegeneracy)"
    The **rank** of a bilinear form $f$ is defined as the rank of its Gram matrix (by Theorem 9.9, this is independent of the choice of basis).

    If $\operatorname{rank}(f) = n$ (i.e., the Gram matrix is invertible), then $f$ is called **nondegenerate**. Equivalently, $f$ is nondegenerate if and only if: whenever $f(\mathbf{x}, \mathbf{y}) = 0$ for all $\mathbf{y} \in V$, it follows that $\mathbf{x} = \mathbf{0}$.

!!! theorem "Theorem 9.10 (Polarization identity)"
    Suppose $\operatorname{char}(\mathbb{F}) \neq 2$. A symmetric bilinear form $f$ and its associated quadratic form $Q(\mathbf{x}) = f(\mathbf{x}, \mathbf{x})$ mutually determine each other through the **polarization identity**:

    $$
    f(\mathbf{x}, \mathbf{y}) = \frac{1}{2}\bigl[Q(\mathbf{x} + \mathbf{y}) - Q(\mathbf{x}) - Q(\mathbf{y})\bigr]
    $$

??? proof "Proof"
    Expand $Q(\mathbf{x}+\mathbf{y}) = f(\mathbf{x}+\mathbf{y}, \mathbf{x}+\mathbf{y})$:

    $$
    Q(\mathbf{x}+\mathbf{y}) = f(\mathbf{x}, \mathbf{x}) + f(\mathbf{x}, \mathbf{y}) + f(\mathbf{y}, \mathbf{x}) + f(\mathbf{y}, \mathbf{y})
    $$

    By symmetry $f(\mathbf{x}, \mathbf{y}) = f(\mathbf{y}, \mathbf{x})$, so $Q(\mathbf{x}+\mathbf{y}) = Q(\mathbf{x}) + 2f(\mathbf{x}, \mathbf{y}) + Q(\mathbf{y})$. Solving for $f(\mathbf{x}, \mathbf{y})$ gives the result. $\blacksquare$

!!! definition "Definition 9.12 (Symmetric and antisymmetric bilinear forms)"
    A bilinear form $f$ is called

    - **symmetric**: if $f(\mathbf{x}, \mathbf{y}) = f(\mathbf{y}, \mathbf{x})$ for all $\mathbf{x}, \mathbf{y}$;
    - **antisymmetric** (skew-symmetric): if $f(\mathbf{x}, \mathbf{y}) = -f(\mathbf{y}, \mathbf{x})$ for all $\mathbf{x}, \mathbf{y}$.

    The matrix of a symmetric bilinear form satisfies $A = A^T$; the matrix of an antisymmetric bilinear form satisfies $A = -A^T$ (with zero diagonal entries).

!!! theorem "Theorem 9.11 (Decomposition of bilinear forms)"
    Suppose $\operatorname{char}(\mathbb{F}) \neq 2$. Any bilinear form $f$ can be uniquely decomposed into the sum of a symmetric part and an antisymmetric part:

    $$
    f = f_s + f_a, \quad f_s(\mathbf{x},\mathbf{y}) = \frac{f(\mathbf{x},\mathbf{y}) + f(\mathbf{y},\mathbf{x})}{2}, \quad f_a(\mathbf{x},\mathbf{y}) = \frac{f(\mathbf{x},\mathbf{y}) - f(\mathbf{y},\mathbf{x})}{2}
    $$

??? proof "Proof"
    Direct verification shows $f_s$ is symmetric, $f_a$ is antisymmetric, and $f = f_s + f_a$. Uniqueness: if $f = g_s + g_a$ ($g_s$ symmetric, $g_a$ antisymmetric), then $f(\mathbf{x},\mathbf{y}) + f(\mathbf{y},\mathbf{x}) = 2g_s(\mathbf{x},\mathbf{y})$, so $g_s = f_s$, and thus $g_a = f_a$. $\blacksquare$

!!! theorem "Theorem 9.12 (Canonical form of symmetric bilinear forms)"
    Suppose $\operatorname{char}(\mathbb{F}) \neq 2$ and $f$ is a symmetric bilinear form on $V$. Then there exists a basis $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ of $V$ such that

    $$
    f(\mathbf{e}_i, \mathbf{e}_j) = 0 \quad (i \neq j)
    $$

    That is, the matrix of $f$ with respect to this basis is diagonal.

??? proof "Proof"
    By induction on $n$. The case $n = 1$ is trivial. Assume the result holds for $(n-1)$-dimensional spaces.

    If $f \equiv 0$, any basis works. Otherwise, there exists $\mathbf{v}$ with $f(\mathbf{v}, \mathbf{v}) \neq 0$ (if $f(\mathbf{v}, \mathbf{v}) = 0$ for all $\mathbf{v}$, the polarization identity gives $f \equiv 0$, a contradiction). Let $\mathbf{e}_1 = \mathbf{v}$ and $d = f(\mathbf{e}_1, \mathbf{e}_1) \neq 0$.

    Let $W = \{\mathbf{w} \in V : f(\mathbf{e}_1, \mathbf{w}) = 0\}$. Then $V = \operatorname{span}\{\mathbf{e}_1\} \oplus W$ (for any $\mathbf{x} \in V$, let $\mathbf{w} = \mathbf{x} - \frac{f(\mathbf{e}_1, \mathbf{x})}{d}\mathbf{e}_1$, then $f(\mathbf{e}_1, \mathbf{w}) = 0$).

    By the induction hypothesis, $f|_W$ can be diagonalized with respect to some basis. Combining $\mathbf{e}_1$ with the basis of $W$ gives the result. $\blacksquare$

!!! definition "Definition 9.13 (Orthogonal complement)"
    Let $f$ be a bilinear form on $V$ and $S \subseteq V$. The **left orthogonal complement** and **right orthogonal complement** of $S$ with respect to $f$ are

    $$
    S^{\perp_L} = \{\mathbf{x} \in V : f(\mathbf{x}, \mathbf{s}) = 0, \forall \mathbf{s} \in S\}, \quad S^{\perp_R} = \{\mathbf{y} \in V : f(\mathbf{s}, \mathbf{y}) = 0, \forall \mathbf{s} \in S\}
    $$

    If $f$ is symmetric, then $S^{\perp_L} = S^{\perp_R}$, denoted simply $S^\perp$. The **radical** of $V$ is $\operatorname{rad}(f) = V^\perp$. The form $f$ is nondegenerate if and only if $\operatorname{rad}(f) = \{\mathbf{0}\}$.

!!! example "Example 9.11"
    The bilinear form $f(\mathbf{x}, \mathbf{y}) = x_1y_1 + x_1y_2 + x_2y_1 + 3x_2y_2 - x_3y_3$ on $\mathbb{R}^3$ has Gram matrix

    $$
    A = \begin{pmatrix} 1 & 1 & 0 \\ 1 & 3 & 0 \\ 0 & 0 & -1 \end{pmatrix}
    $$

    $\det(A) = (3-1)(-1) = -2 \neq 0$, so $f$ is nondegenerate. Since $A$ is symmetric, $f$ is a symmetric bilinear form.

!!! example "Example 9.12"
    The antisymmetric bilinear form $f(\mathbf{x}, \mathbf{y}) = x_1y_2 - x_2y_1$ on $\mathbb{R}^3$ has matrix

    $$
    A = \begin{pmatrix} 0 & 1 & 0 \\ -1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
    $$

    $\operatorname{rank}(A) = 2$, $\operatorname{rad}(f) = \operatorname{span}\{\mathbf{e}_3\}$, so $f$ is degenerate.

!!! example "Example 9.13"
    Verifying the polarization identity. Let $Q(x_1,x_2) = 2x_1^2 + 3x_1x_2 + x_2^2$, with symmetric matrix $A = \begin{pmatrix} 2 & 3/2 \\ 3/2 & 1 \end{pmatrix}$.

    Take $\mathbf{x} = (1,0)^T$, $\mathbf{y} = (0,1)^T$:

    $$
    Q(\mathbf{x}+\mathbf{y}) = Q(1,1) = 2 + 3 + 1 = 6, \quad Q(\mathbf{x}) = 2, \quad Q(\mathbf{y}) = 1
    $$

    $$
    f(\mathbf{x},\mathbf{y}) = \frac{6-2-1}{2} = \frac{3}{2} = a_{12} \checkmark
    $$

!!! theorem "Theorem 9.13 (Canonical form of antisymmetric bilinear forms)"
    Let $f$ be an antisymmetric bilinear form on a finite-dimensional vector space $V$. Then there exists a basis $\{\mathbf{e}_1, \mathbf{f}_1, \ldots, \mathbf{e}_m, \mathbf{f}_m, \mathbf{g}_1, \ldots, \mathbf{g}_k\}$ such that

    $$
    f(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}, \quad f(\mathbf{e}_i, \mathbf{e}_j) = f(\mathbf{f}_i, \mathbf{f}_j) = 0
    $$

    and $f(\mathbf{g}_s, \cdot) \equiv 0$. In particular, the rank of an antisymmetric bilinear form is always even, equal to $2m$.

??? proof "Proof"
    If $f \equiv 0$, the result is trivial. Otherwise, there exist $\mathbf{x}, \mathbf{y}$ with $f(\mathbf{x}, \mathbf{y}) = c \neq 0$. Let $\mathbf{e}_1 = \mathbf{x}$ and $\mathbf{f}_1 = \mathbf{y}/c$, so that $f(\mathbf{e}_1, \mathbf{f}_1) = 1$.

    Let $W = \{\mathbf{w} : f(\mathbf{e}_1, \mathbf{w}) = 0 \text{ and } f(\mathbf{f}_1, \mathbf{w}) = 0\}$. One can verify that $V = \operatorname{span}\{\mathbf{e}_1, \mathbf{f}_1\} \oplus W$ (for any $\mathbf{v}$, let $\mathbf{w} = \mathbf{v} - f(\mathbf{v}, \mathbf{f}_1)\mathbf{e}_1 + f(\mathbf{v}, \mathbf{e}_1)\mathbf{f}_1$, and check $\mathbf{w} \in W$).

    Apply the procedure recursively to $f|_W$. $\blacksquare$

---

## 9.9 Symplectic Spaces

<div class="context-flow" markdown>

Nondegenerate antisymmetric bilinear form = **symplectic form** → Symplectic spaces are necessarily even-dimensional → **Darboux theorem**: all symplectic spaces of the same dimension are equivalent → Symplectic geometry is the mathematical language of Hamiltonian mechanics

</div>

A symplectic space is a vector space equipped with a nondegenerate antisymmetric bilinear form. It plays a foundational role in classical mechanics (Hamiltonian systems), quantum mechanics, and modern differential geometry.

!!! definition "Definition 9.14 (Symplectic form and symplectic space)"
    Let $V$ be a finite-dimensional vector space over a field $\mathbb{F}$ ($\operatorname{char}(\mathbb{F}) \neq 2$). A **symplectic form** on $V$ is a nondegenerate antisymmetric bilinear form $\omega: V \times V \to \mathbb{F}$. The pair $(V, \omega)$ is called a **symplectic space**.

!!! theorem "Theorem 9.14 (Dimension of a symplectic space)"
    The dimension of a symplectic space is necessarily even.

??? proof "Proof"
    Let $(V, \omega)$ be a symplectic space with $\dim V = n$. The matrix $A$ of $\omega$ with respect to any basis satisfies $A^T = -A$ (antisymmetric), so

    $$
    \det(A) = \det(A^T) = \det(-A) = (-1)^n \det(A)
    $$

    Since $\omega$ is nondegenerate, $\det(A) \neq 0$, so $(-1)^n = 1$, meaning $n$ is even. $\blacksquare$

!!! definition "Definition 9.15 (Symplectic basis)"
    Let $(V, \omega)$ be a $2n$-dimensional symplectic space. A basis $\{\mathbf{e}_1, \ldots, \mathbf{e}_n, \mathbf{f}_1, \ldots, \mathbf{f}_n\}$ of $V$ is called a **symplectic basis** (or Darboux basis) if

    $$
    \omega(\mathbf{e}_i, \mathbf{e}_j) = 0, \quad \omega(\mathbf{f}_i, \mathbf{f}_j) = 0, \quad \omega(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}
    $$

    With respect to a symplectic basis, the matrix of $\omega$ is the standard symplectic matrix

    $$
    J_{2n} = \begin{pmatrix} O_n & I_n \\ -I_n & O_n \end{pmatrix}
    $$

!!! theorem "Theorem 9.15 (Darboux theorem — canonical form of symplectic spaces)"
    Every symplectic space $(V, \omega)$ admits a symplectic basis. Consequently, all $2n$-dimensional symplectic spaces (over the same field) are isomorphic to each other.

??? proof "Proof"
    This is a direct corollary of Theorem 9.13 in the nondegenerate case ($k = 0$). Since $\omega$ is nondegenerate, the rank equals $\dim V = 2m$, and Theorem 9.13 provides a basis $\{\mathbf{e}_1, \mathbf{f}_1, \ldots, \mathbf{e}_m, \mathbf{f}_m\}$ satisfying $\omega(\mathbf{e}_i, \mathbf{f}_j) = \delta_{ij}$ and $\omega(\mathbf{e}_i, \mathbf{e}_j) = \omega(\mathbf{f}_i, \mathbf{f}_j) = 0$. This is precisely a symplectic basis. $\blacksquare$

!!! definition "Definition 9.16 (Symplectic matrix and symplectic group)"
    A $2n \times 2n$ real matrix $M$ is called a **symplectic matrix** if

    $$
    M^T J_{2n} M = J_{2n}
    $$

    The set of all $2n \times 2n$ symplectic matrices forms a group under matrix multiplication, called the **symplectic group**, denoted $\operatorname{Sp}(2n, \mathbb{F})$.

!!! theorem "Theorem 9.16 (Determinant of a symplectic matrix)"
    The determinant of a symplectic matrix equals $1$.

??? proof "Proof"
    From $M^T J M = J$, taking determinants gives $\det(M)^2 \det(J) = \det(J)$. Since $\det(J) = 1$ (this can be computed directly or deduced from $J^2 = -I$ giving $\det(J)^2 = 1$), we get $\det(M)^2 = 1$, so $\det(M) = \pm 1$.

    To show $\det(M) = 1$, note that the symplectic group $\operatorname{Sp}(2n, \mathbb{R})$ is connected (one can continuously deform $M$ to $I_{2n}$), and the determinant is a continuous function with $\det(I_{2n}) = 1$, so $\det(M) = 1$.

    An alternative algebraic proof: write $M$ in block form $M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$. From $M^TJM = J$ one obtains $A^TD - C^TB = I$, and the sign of $\det(M)$ can be determined to be $+1$ via a Pfaffian argument. $\blacksquare$

!!! definition "Definition 9.17 (Lagrangian subspace)"
    Let $(V, \omega)$ be a $2n$-dimensional symplectic space. A subspace $L \subseteq V$ is called a **Lagrangian subspace** if $L = L^\perp$ (with respect to $\omega$), i.e.,

    $$
    \omega(\mathbf{x}, \mathbf{y}) = 0, \quad \forall \mathbf{x}, \mathbf{y} \in L
    $$

    and $\dim L = n$ (the maximum dimension of an isotropic subspace).

!!! theorem "Theorem 9.17 (Dimension bound for isotropic subspaces)"
    Let $(V, \omega)$ be a $2n$-dimensional symplectic space and $L$ an isotropic subspace (i.e., $\omega|_{L \times L} = 0$). Then $\dim L \leq n$, with equality if and only if $L$ is a Lagrangian subspace.

??? proof "Proof"
    Let $L^\perp = \{\mathbf{v} \in V : \omega(\mathbf{v}, \mathbf{l}) = 0, \forall \mathbf{l} \in L\}$. Since $\omega$ is nondegenerate, the map $V \to V^*$, $\mathbf{v} \mapsto \omega(\mathbf{v}, \cdot)$ is an isomorphism, so $\dim L^\perp = 2n - \dim L$.

    The isotropic condition means $L \subseteq L^\perp$, so $\dim L \leq \dim L^\perp = 2n - \dim L$, giving $\dim L \leq n$.

    Equality holds if and only if $L = L^\perp$, i.e., $L$ is a Lagrangian subspace. $\blacksquare$

!!! example "Example 9.14"
    The standard symplectic space $(\mathbb{R}^4, \omega)$ with symplectic form defined by the matrix $J_4 = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{pmatrix}$.

    The standard symplectic basis is $\mathbf{e}_1 = (1,0,0,0)^T$, $\mathbf{e}_2 = (0,1,0,0)^T$, $\mathbf{f}_1 = (0,0,1,0)^T$, $\mathbf{f}_2 = (0,0,0,1)^T$.

    $L_1 = \operatorname{span}\{\mathbf{e}_1, \mathbf{e}_2\}$ is a Lagrangian subspace: $\omega(\mathbf{e}_1, \mathbf{e}_2) = 0$ and $\dim L_1 = 2 = n$.

    $L_2 = \operatorname{span}\{\mathbf{f}_1, \mathbf{f}_2\}$ is also a Lagrangian subspace.

!!! example "Example 9.15"
    An example of a symplectic matrix. The matrix $M = \begin{pmatrix} a & 0 & 0 & b \\ 0 & d & -c & 0 \\ 0 & c & d & 0 \\ -b & 0 & 0 & a \end{pmatrix}$ (where $a^2 + b^2 = 1$ and $c^2 + d^2 = 1$) is a $4 \times 4$ symplectic matrix. The condition $M^TJ_4M = J_4$ can be verified by direct block computation. $\det(M) = (a^2+b^2)(c^2+d^2) = 1$.

!!! example "Example 9.16"
    Application to Hamiltonian mechanics. The phase space $\mathbb{R}^{2n}$ of classical mechanics uses generalized coordinates $(q_1, \ldots, q_n)$ and generalized momenta $(p_1, \ldots, p_n)$. The symplectic form is

    $$
    \omega = \sum_{i=1}^n dp_i \wedge dq_i
    $$

    Hamilton's equations $\dot{q}_i = \frac{\partial H}{\partial p_i}$, $\dot{p}_i = -\frac{\partial H}{\partial q_i}$ can be compactly written as

    $$
    \dot{\mathbf{z}} = J_{2n} \nabla H(\mathbf{z})
    $$

    where $\mathbf{z} = (q_1, \ldots, q_n, p_1, \ldots, p_n)^T$. Canonical transformations (variable substitutions that preserve the form of Hamilton's equations) correspond precisely to symplectic matrices.

---

## 9.10 Hermitian Forms

<div class="context-flow" markdown>

**Real** $\to$ **Complex**: inner products change from bilinear to **sesquilinear** → Hermitian form $H(\mathbf{x}) = \mathbf{x}^*A\mathbf{x}$ ($A = A^*$ Hermitian matrix) → **Law of inertia** generalized to the complex field → Signature remains a complete invariant

</div>

When the base field extends from $\mathbb{R}$ to $\mathbb{C}$, the natural generalization of symmetric bilinear forms is the Hermitian form. Hermitian forms play a fundamental role in quantum mechanics (observables correspond to Hermitian operators) and unitary geometry.

!!! definition "Definition 9.18 (Sesquilinear form)"
    Let $V$ be a complex vector space. A map $f: V \times V \to \mathbb{C}$ is called a **sesquilinear form** if

    - It is conjugate-linear in the first variable: $f(\alpha\mathbf{x}_1 + \beta\mathbf{x}_2, \mathbf{y}) = \bar{\alpha}f(\mathbf{x}_1, \mathbf{y}) + \bar{\beta}f(\mathbf{x}_2, \mathbf{y})$
    - It is linear in the second variable: $f(\mathbf{x}, \alpha\mathbf{y}_1 + \beta\mathbf{y}_2) = \alpha f(\mathbf{x}, \mathbf{y}_1) + \beta f(\mathbf{x}, \mathbf{y}_2)$

    (Here we follow the physics convention where the first variable is conjugated. Some mathematical references use the opposite convention.)

!!! definition "Definition 9.19 (Hermitian form)"
    A sesquilinear form $f$ is called a **Hermitian form** if it satisfies conjugate symmetry:

    $$
    f(\mathbf{x}, \mathbf{y}) = \overline{f(\mathbf{y}, \mathbf{x})}, \quad \forall \mathbf{x}, \mathbf{y} \in V
    $$

    In particular, $f(\mathbf{x}, \mathbf{x}) = \overline{f(\mathbf{x}, \mathbf{x})}$, so $f(\mathbf{x}, \mathbf{x}) \in \mathbb{R}$.

!!! definition "Definition 9.20 (Hermitian matrix representation)"
    Let $\mathcal{B} = \{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ be a basis of $V$, and let $A = (a_{ij})$ with $a_{ij} = f(\mathbf{e}_i, \mathbf{e}_j)$ be the matrix of the Hermitian form $f$. Then

    $$
    f(\mathbf{x}, \mathbf{y}) = \mathbf{x}^* A \mathbf{y} = \sum_{i,j} \bar{x}_i a_{ij} y_j
    $$

    where $\mathbf{x}^* = \bar{\mathbf{x}}^T$ is the conjugate transpose. The matrix $A$ satisfies $A^* = A$ (Hermitian matrix).

!!! theorem "Theorem 9.18 (Change of basis for Hermitian forms)"
    Let the Hermitian form $f$ have matrix $A$ with respect to basis $\mathcal{B}$ and matrix $A'$ with respect to basis $\mathcal{B}'$, with change-of-basis matrix $C$. Then

    $$
    A' = C^* A C
    $$

    That is, $A'$ and $A$ are **$^*$-congruent** (conjugate congruent).

??? proof "Proof"
    Let $\mathcal{B}' = \{\mathbf{e}_1', \ldots, \mathbf{e}_n'\}$ with $(\mathbf{e}_1', \ldots, \mathbf{e}_n') = (\mathbf{e}_1, \ldots, \mathbf{e}_n)C$. Then

    $$
    a_{ij}' = f(\mathbf{e}_i', \mathbf{e}_j') = f\!\left(\sum_k c_{ki}\mathbf{e}_k, \sum_l c_{lj}\mathbf{e}_l\right) = \sum_{k,l} \bar{c}_{ki} a_{kl} c_{lj}
    $$

    That is, $A' = C^*AC$. $\blacksquare$

!!! theorem "Theorem 9.19 (Canonical form of Hermitian forms)"
    Let $f$ be a Hermitian form on a complex vector space $V$. There exists a basis of $V$ with respect to which the matrix of $f$ is a diagonal matrix $\operatorname{diag}(d_1, \ldots, d_n)$ with $d_i \in \mathbb{R}$.

    Furthermore, after appropriate rescaling, this can be reduced to

    $$
    \operatorname{diag}(\underbrace{1, \ldots, 1}_{p}, \underbrace{-1, \ldots, -1}_{q}, \underbrace{0, \ldots, 0}_{n-r})
    $$

??? proof "Proof"
    The proof is analogous to the completing-the-square method for real symmetric bilinear forms. If there exists $\mathbf{v}$ with $f(\mathbf{v}, \mathbf{v}) \neq 0$ (this value is real), let $\mathbf{e}_1 = \mathbf{v}$, take the orthogonal complement $W = \{\mathbf{w} : f(\mathbf{e}_1, \mathbf{w}) = 0\}$, and diagonalize $f|_W$ by induction.

    If $f(\mathbf{v}, \mathbf{v}) = 0$ for all $\mathbf{v}$, the Hermitian version of the polarization identity

    $$
    f(\mathbf{x}, \mathbf{y}) = \frac{1}{4}\bigl[f(\mathbf{x}+\mathbf{y}, \mathbf{x}+\mathbf{y}) - f(\mathbf{x}-\mathbf{y}, \mathbf{x}-\mathbf{y}) + if(\mathbf{x}+i\mathbf{y}, \mathbf{x}+i\mathbf{y}) - if(\mathbf{x}-i\mathbf{y}, \mathbf{x}-i\mathbf{y})\bigr]
    $$

    shows $f \equiv 0$. Thus the induction proceeds.

    Rescaling: if $d_k > 0$, let $\mathbf{e}_k' = \mathbf{e}_k/\sqrt{d_k}$; if $d_k < 0$, let $\mathbf{e}_k' = \mathbf{e}_k/\sqrt{|d_k|}$. $\blacksquare$

!!! theorem "Theorem 9.20 (Law of inertia for Hermitian forms)"
    The number of positive terms $p$ and negative terms $q$ in the canonical form of a Hermitian form are invariants, independent of the choice of basis. The pair $(p, q)$ is called the **signature** of the Hermitian form.

??? proof "Proof"
    The proof is entirely analogous to Sylvester's law of inertia for real quadratic forms (Theorem 9.4). Assuming two canonical forms have different numbers of positive terms $p > s$, construct subspaces $V_1$ (dimension $p$, on which $f > 0$) and $V_2$ (dimension $n-s$, on which $f \leq 0$). The dimension argument $p + (n-s) > n$ gives $V_1 \cap V_2 \neq \{\mathbf{0}\}$, leading to a contradiction. $\blacksquare$

!!! theorem "Theorem 9.21 (Spectral theorem for Hermitian matrices)"
    A Hermitian matrix $A$ ($A^* = A$) has all real eigenvalues, and $A$ can be unitarily diagonalized: there exists a unitary matrix $U$ such that

    $$
    U^* A U = \operatorname{diag}(\lambda_1, \ldots, \lambda_n), \quad \lambda_i \in \mathbb{R}
    $$

??? proof "Proof"
    **Eigenvalues are real**: If $A\mathbf{v} = \lambda\mathbf{v}$, $\mathbf{v} \neq \mathbf{0}$, then $\lambda \mathbf{v}^*\mathbf{v} = \mathbf{v}^*A\mathbf{v} = (A\mathbf{v})^*\mathbf{v} = \bar{\lambda}\mathbf{v}^*\mathbf{v}$ (using $A^* = A$), so $\lambda = \bar{\lambda}$, meaning $\lambda \in \mathbb{R}$.

    **Unitary diagonalization**: Eigenvectors for distinct eigenvalues are orthogonal (the proof is the same as in the real case). Apply Gram-Schmidt orthogonalization within each eigenspace and combine to form the unitary matrix $U$. $\blacksquare$

!!! definition "Definition 9.21 (Positive definite Hermitian form)"
    A Hermitian form $f$ is called **positive definite** if $f(\mathbf{x}, \mathbf{x}) > 0$ for all $\mathbf{x} \neq \mathbf{0}$. The equivalent conditions are analogous to the real case:

    - All eigenvalues of the Hermitian matrix $A$ are positive;
    - There exists an invertible matrix $C$ such that $A = C^*C$;
    - All leading principal minors of $A$ are positive.

!!! example "Example 9.17"
    The Hermitian form $f(\mathbf{x}, \mathbf{y}) = 2\bar{x}_1y_1 + (1+i)\bar{x}_1y_2 + (1-i)\bar{x}_2y_1 + 3\bar{x}_2y_2$ has Hermitian matrix

    $$
    A = \begin{pmatrix} 2 & 1+i \\ 1-i & 3 \end{pmatrix}
    $$

    Verification that $A^* = A$: $\overline{(1+i)} = 1-i = a_{21}$ $\checkmark$. $\Delta_1 = 2 > 0$, $\Delta_2 = 6 - |1+i|^2 = 6 - 2 = 4 > 0$, so $A$ is positive definite.

!!! example "Example 9.18"
    Find the eigenvalues and unitary diagonalization of the Hermitian matrix $A = \begin{pmatrix} 1 & i \\ -i & 1 \end{pmatrix}$.

    Characteristic polynomial: $(1-\lambda)^2 - i(-i) = (1-\lambda)^2 - 1 = \lambda^2 - 2\lambda = \lambda(\lambda - 2)$.

    $\lambda_1 = 0$: $\mathbf{v}_1 = \frac{1}{\sqrt{2}}(i, 1)^T$.

    $\lambda_2 = 2$: $\mathbf{v}_2 = \frac{1}{\sqrt{2}}(-i, 1)^T$.

    Unitary matrix $U = \frac{1}{\sqrt{2}}\begin{pmatrix} i & -i \\ 1 & 1 \end{pmatrix}$, $U^*AU = \begin{pmatrix} 0 & 0 \\ 0 & 2 \end{pmatrix}$.

    Signature $(p, q) = (1, 0)$, rank $r = 1$. $A$ is positive semidefinite but not positive definite.

!!! example "Example 9.19"
    Anti-Hermitian (skew-Hermitian) forms. If a sesquilinear form $f$ satisfies $f(\mathbf{x}, \mathbf{y}) = -\overline{f(\mathbf{y}, \mathbf{x})}$, then $f(\mathbf{x}, \mathbf{x})$ is purely imaginary. The corresponding matrix $A$ satisfies $A^* = -A$ (anti-Hermitian matrix), and all its eigenvalues are purely imaginary.

    For example, $A = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$ is an anti-Hermitian matrix (also a real antisymmetric matrix) with eigenvalues $\pm i$.
