# Chapter 13B: λ-matrices and Rational Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Polynomial Algebra (Ch00) · Matrix Algebra (Ch02) · Jordan Canonical Form (Ch12) · Quotient Spaces (Ch13A)

**Chapter Outline**: From Number Fields to Polynomial Rings → Definition of $\lambda$-matrices and Elementary Operations → Smith Normal Form & Uniqueness Theorem → Determinantal Divisors and Invariant Factors → Elementary Divisors → Necessary and Sufficient Conditions for Matrix Similarity (Equivalence of Characteristic Matrices) → Construction of Companion Matrices → The Rational Canonical Form (RCF) Theorem → Relationship with Jordan Canonical Form → Applications: Matrix Analysis over General Fields, Canonical Forms in Control Systems

**Extension**: The Rational Canonical Form solves the problem of matrix canonical forms over general fields (e.g., $\mathbb{Q}$) without requiring the field to be algebraically closed (e.g., $\mathbb{C}$); it is a concrete application of the structure theorem for modules over a Principal Ideal Domain (PID).

</div>

While the Jordan Canonical Form is theoretically elegant, its construction relies on the existence of eigenvalues within the field. When working over the rational field $\mathbb{Q}$, the characteristic polynomial may not be factorable. The **Rational Canonical Form** (RCF) overcomes this limitation by providing a universal matrix standard form applicable to any field through the decomposition of polynomial rings. This chapter uses the elementary transformations of $\lambda$-matrices to reveal this profound algebraic construction.

---

## 13B.1 λ-matrices and Smith Normal Form

!!! definition "Definition 13B.1 (λ-matrix)"
    A matrix whose entries are polynomials in $\lambda$ is called a **$\lambda$-matrix**.
    The characteristic matrix $\lambda I - A$ is the most quintessential example.

!!! theorem "Theorem 13B.1 (Smith Normal Form)"
    Every $n \times n$ $\lambda$-matrix $A(\lambda)$ can be transformed via elementary operations into a unique diagonal form:
    $$S(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0)$$
    where each $d_i(\lambda)$ is a monic polynomial such that $d_i(\lambda) \mid d_{i+1}(\lambda)$.
    These polynomials are called the **Invariant Factors** of $A(\lambda)$.

---

## 13B.2 Similarity and Elementary Divisors

!!! theorem "Theorem 13B.2 (Similarity Criterion)"
    Two $n \times n$ matrices $A$ and $B$ are similar iff their characteristic matrices $\lambda I - A$ and $\lambda I - B$ have the same Smith Normal Form (i.e., identical invariant factors).

---

## 13B.3 Companion Matrices and Rational Canonical Form

!!! definition "Definition 13B.2 (Companion Matrix $C(p)$)"
    For a monic polynomial $p(\lambda) = \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_0$, the **companion matrix** is defined as:
    $$C(p) = \begin{pmatrix} 0 & 0 & \cdots & 0 & -a_0 \\ 1 & 0 & \cdots & 0 & -a_1 \\ 0 & 1 & \cdots & 0 & -a_2 \\ \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & \cdots & 1 & -a_{k-1} \end{pmatrix}$$
    **Property**: The characteristic and minimal polynomials of $C(p)$ are both $p(\lambda)$.

!!! theorem "Theorem 13B.3 (Rational Canonical Form)"
    Every square matrix $A$ is similar to a block diagonal matrix where each block is the companion matrix of an invariant factor:
    $$R = \operatorname{diag}(C(d_1), C(d_2), \ldots, C(d_k))$$
    This is known as the **Rational Canonical Form** of $A$. Note that the degree of $d_i$ increases with $i$.

---

## Exercises

**1. [Smith Form] Calculate the invariant factors of $\begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    1. **Determinantal Divisors**:
       - $D_1(\lambda) = \gcd(\lambda, 1, 0, \lambda) = 1$.
       - $D_2(\lambda) = \det = \lambda^2$.
    2. **Invariant Factors**:
       - $d_1(\lambda) = D_1 = 1$.
       - $d_2(\lambda) = D_2 / D_1 = \lambda^2$.
    **Conclusion**: The invariant factors are $1, \lambda^2$.

**2. [Similarity] If $A$ and $B$ have the same characteristic polynomial $\lambda^2$, are they necessarily similar?**

??? success "Solution"
    **Conclusion**: Not necessarily.
    **Analysis**:
    - Identical characteristic polynomials only mean the product of the invariant factors is the same.
    - Let $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$; its invariant factors are $1, \lambda^2$.
    - Let $B = \begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$; its invariant factors are $\lambda, \lambda$.
    - Since the sequences of invariant factors differ, they are not similar. This reflects the difference in Jordan block structures.

**3. [Companion] Write the companion matrix for $p(\lambda) = \lambda^2 - 3\lambda + 2$.**

??? success "Solution"
    **Construction:**
    By definition, $a_1 = -3, a_0 = 2$.
    The companion matrix is $\begin{pmatrix} 0 & -a_0 \\ 1 & -a_1 \end{pmatrix} = \begin{pmatrix} 0 & -2 \\ 1 & 3 \end{pmatrix}$.
    Verification: The eigenvalues are 1 and 2, which are the roots of the polynomial.

**4. [Minimal Polynomial] In the Rational Canonical Form, which block corresponds to the minimal polynomial?**

??? success "Solution"
    **Conclusion:**
    The companion matrix block corresponding to the last non-trivial invariant factor $d_k(\lambda)$.
    **Reasoning**: In the Smith form, $d_i \mid d_{i+1}$, so $d_k$ contains the highest powers of all elementary divisors, which is the definition of the minimal polynomial.

**5. [Divisibility] Prove $d_i(\lambda) \mid d_{i+1}(\lambda)$.**

??? success "Solution"
    **Logic:**
    1. Invariant factors are defined as $d_i = D_i / D_{i-1}$, where $D_i$ is the GCD of all $i \times i$ minors.
    2. Any $(i+1) \times (i+1)$ minor is a linear combination of $i \times i$ minors.
    3. Thus, $D_i$ must divide $D_{i+1}$.
    4. Advanced algebraic lemmas prove the quotient sequence also satisfies this divisibility chain.

**6. [Contrast] What is the primary difference between Jordan and Rational Canonical Forms?**

??? success "Solution"
    **Key Differences:**
    1. **Field Dependency**: Jordan form requires the characteristic polynomial to factor completely within the field (usually $\mathbb{C}$); RCF exists and is unique over any field (e.g., $\mathbb{Q}, \mathbb{R}$).
    2. **Decomposition Depth**: Jordan form decomposes the space into irreducible linear factor powers (eigenspaces); RCF decomposes into companion blocks of irreducible polynomials (cyclic spaces).

**7. [Elementary Divisors] If the invariant factors are $1, (\lambda-1)(\lambda-2)$, what are the elementary divisors?**

??? success "Solution"
    **Calculation:**
    Elementary divisors are the powers of irreducible factors obtained by factoring the invariant factors over an algebraically closed field.
    Factoring $(\lambda-1)(\lambda-2)$ yields $(\lambda-1)$ and $(\lambda-2)$.
    **Conclusion**: The elementary divisors are $(\lambda-1)$ and $(\lambda-2)$.

**8. [Rank] In the Smith form of $\lambda I - A$, why is the number of non-zero diagonal entries $r$ always $n$?**

??? success "Solution"
    **Conclusion:**
    $r = n$.
    **Reasoning**: The determinant of the characteristic matrix $\lambda I - A$ is a polynomial of degree $n$, which is not identically zero. Therefore, the $\lambda$-matrix is of full rank, and must have $n$ non-zero diagonal entries in its Smith form.

**9. [Calculation] Find the invariant factors of $J_2(\lambda_0)$.**

??? success "Solution"
    **Analysis:**
    $A = \begin{pmatrix} \lambda_0 & 1 \\ 0 & \lambda_0 \end{pmatrix} \implies \lambda I - A = \begin{pmatrix} \lambda-\lambda_0 & -1 \\ 0 & \lambda-\lambda_0 \end{pmatrix}$.
    1. $D_1 = \gcd(\lambda-\lambda_0, -1, 0, \lambda-\lambda_0) = 1$.
    2. $D_2 = (\lambda-\lambda_0)^2$.
    **Conclusion**: The invariant factors are $1, (\lambda-\lambda_0)^2$.

**10. [Application] Why is RCF important in computational algebra?**

??? success "Solution"
    **Significance:**
    1. **Avoiding Approximations**: Finding eigenvalues involves root-finding, which introduces precision loss. Computing the RCF only requires exact polynomial arithmetic.
    2. **Field Universality**: For matrices where exact eigenvalues cannot be found (e.g., over $\mathbb{Q}$), the RCF provides the only unique, exact way to describe the matrix structure.

## Chapter Summary

The Rational Canonical Form provides a universal characterization of matrix similarity classes over any field:

1.  **Field Independence**: By utilizing companion blocks of irreducible polynomials, RCF eliminates the dependency on the complex field, becoming the core of operator theory in abstract algebra.
2.  **Polynomial Logic**: The theory of invariant factors reveals the deep module structure behind the characteristic matrix $\lambda I - A$, establishing the ultimate algebraic criterion for similarity.
3.  **Computational Precision**: Compared to the instability of eigenvalue calculation, the construction of Smith forms based on elementary operations provides a robust algorithmic foundation for exact algebra software.
