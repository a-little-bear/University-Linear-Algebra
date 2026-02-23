# Chapter 13B: λ-matrices and Rational Canonical Form

<div class="context-flow" markdown>

**Prerequisites**: Polynomial Algebra (Ch00) · Matrix Algebra (Ch02) · Jordan Canonical Form (Ch12)

**Chapter Outline**: Definition of $\lambda$-matrices → Elementary Operations and Equivalence → Smith Normal Form → Invariant Factors → Elementary Divisors → Necessary and Sufficient Conditions for Similarity → Companion Matrices → Rational Canonical Form (RCF) → Structural Comparison (Jordan vs. Rational)

**Extension**: The Rational Canonical Form solves the problem of matrix canonical forms over general fields (e.g., $\mathbb{Q}$) without requiring the field to be algebraically closed (e.g., $\mathbb{C}$); it is a concrete application of the structure theorem for modules over a PID.

</div>

While the Jordan Canonical Form is theoretically perfect, its construction depends on the existence of eigenvalues within the underlying field. Working over the field of rational numbers $\mathbb{Q}$, the characteristic polynomial may not factor into linear terms. The Rational Canonical Form (RCF) overcomes this limitation by providing a universal standard form applicable over any field through the decomposition of polynomial rings.

---

## 13B.1 λ-matrices and Smith Normal Form

!!! definition "Definition 13B.1 (λ-matrix)"
    A matrix whose entries are polynomials in $\lambda$ is called a **$\lambda$-matrix** (or a matrix polynomial).

!!! theorem "Theorem 13B.1 (Smith Normal Form)"
    Every $n \times n$ $\lambda$-matrix $A(\lambda)$ can be transformed via elementary operations into a unique diagonal form:
    $$S(\lambda) = \operatorname{diag}(d_1(\lambda), d_2(\lambda), \ldots, d_r(\lambda), 0, \ldots, 0)$$
    where each $d_i(\lambda)$ is a monic polynomial such that $d_i(\lambda) \mid d_{i+1}(\lambda)$. These polynomials are called the **invariant factors** of $A(\lambda)$.

---

## 13B.2 Conditions for Similarity

!!! theorem "Theorem 13B.2 (Similarity Theorem)"
    Two $n \times n$ matrices $A$ and $B$ are similar if and only if their characteristic matrices $\lambda I - A$ and $\lambda I - B$ have the same Smith Normal Form (i.e., identical invariant factors).

---

## 13B.3 Companion Matrices and Rational Canonical Form

!!! definition "Definition 13B.2 (Companion Matrix $C(p)$)"
    For a monic polynomial $p(\lambda) = \lambda^k + a_{k-1}\lambda^{k-1} + \cdots + a_0$, the **companion matrix** is:
    $$C(p) = \begin{pmatrix} 0 & 0 & \cdots & -a_0 \\ 1 & 0 & \cdots & -a_1 \\ 0 & 1 & \cdots & -a_2 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & -a_{k-1} \end{pmatrix}$$
    **Property**: The characteristic and minimal polynomials of $C(p)$ are both $p(\lambda)$.

!!! theorem "Theorem 13B.3 (Rational Canonical Form)"
    Every square matrix $A$ is similar to a block diagonal matrix whose blocks are the companion matrices of its invariant factors:
    $$R = \operatorname{diag}(C(d_1), C(d_2), \ldots, C(d_k))$$
    This is known as the **Rational Canonical Form** of $A$.

---

## Exercises

1. **[Smith] Calculate the invariant factors of $\begin{pmatrix} \lambda & 1 \\ 0 & \lambda \end{pmatrix}$.**

   ??? success "Solution"
       The determinant is $\lambda^2$. The GCD of all $1 \times 1$ minors is $D_1(\lambda) = \gcd(\lambda, 1, 0, \lambda) = 1$.
       Thus $d_1 = D_1 = 1$ and $d_2 = D_2/D_1 = \lambda^2$. The invariant factors are $1, \lambda^2$.

2. **[Similarity] If $A$ and $B$ have the same characteristic polynomial $\lambda^2$, are they necessarily similar?**

   ??? success "Solution"
       No. Consider $\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ and $\begin{pmatrix} 0 & 0 \\ 0 & 0 \end{pmatrix}$. The first has invariant factors $1, \lambda^2$, while the second has $\lambda, \lambda$.

3. **[Companion] Write the companion matrix for $\lambda^2 - 3\lambda + 2$.**

   ??? success "Solution"
       $C(p) = \begin{pmatrix} 0 & -2 \\ 1 & 3 \end{pmatrix}$.

4. **[Minimal] In the Rational Canonical Form, which block corresponds to the minimal polynomial?**

   ??? success "Solution"
       The block corresponding to the last non-trivial invariant factor $d_k(\lambda)$.

5. **[Divisibility] Prove $d_i(\lambda) \mid d_{i+1}(\lambda)$.**

   ??? success "Solution"
       This is a structural property of the Smith Normal Form, derived from the definition of the GCD of $k$-th order minors.

6. **[Comparison] What is the main difference between JCF and RCF?**

   ??? success "Solution"
       JCF decomposes polynomials into powers of linear factors (elementary divisors), requiring an algebraically closed field. RCF uses invariant factors directly and works over any field.

7. **[Elementary Divisors] If the invariant factors are $1, (\lambda-1)(\lambda-2)$, what are the elementary divisors?**

   ??? success "Solution"
       $\lambda-1$ and $\lambda-2$.

8. **[Rank] In the Smith form of $\lambda I - A$, why is the number of non-zero diagonal entries always $n$?**

   ??? success "Solution"
       Because $\lambda I - A$ is always non-singular (its determinant is a degree $n$ polynomial).

9. **[Calculation] Find the invariant factors of $J_2(\lambda_0)$.**

   ??? success "Solution"
       $1, (\lambda-\lambda_0)^2$.

10. **[Application] Why is RCF important in computational algebra?**

   ??? success "Solution"
        It avoids the process of finding roots (which usually involves numerical approximation) and only requires elementary operations on matrices (exact algebraic computation).

## Chapter Summary

The Rational Canonical Form provides a universal characterization of matrix similarity classes over any field:

1.  **Field Independence**: By utilizing companion blocks of irreducible polynomials, RCF eliminates the dependency on the complex field, becoming the core of operator theory in abstract algebra.
2.  **Polynomial Logic**: The theory of invariant factors reveals the deep module structure behind the characteristic matrix $\lambda I - A$, establishing the ultimate algebraic criterion for similarity.
3.  **Computational Precision**: Compared to the instability of eigenvalue calculation, the construction of Smith forms based on elementary operations provides a robust algorithmic foundation for exact algebra software.
