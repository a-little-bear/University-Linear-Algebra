# Chapter 00: Polynomials and Fields

<div class="context-flow" markdown>

**Prerequisites**: Basic Algebra · Set Theory · Induction

**Chapter Outline**: Rings and Fields → Polynomial Rings $F[x]$ → Division Algorithm → Euclidean Algorithm and GCD → Irreducibility and Factorization → Fundamental Theorem of Algebra → Roots and Multiplicities → Vieta's Formulas → Interpolation → Minimal Polynomials

**Extension**: Polynomials are the precursors to linear operators; the spectral properties of a matrix are entirely encoded in its characteristic and minimal polynomials.

</div>

Linear algebra is built on the foundation of polynomial rings over fields. A matrix is essentially an operator whose behavior is governed by polynomial identities. Understanding the roots, divisibility, and factorization of polynomials is crucial for developing the theory of eigenvalues, canonical forms, and matrix functions. This chapter establishes the algebraic toolkit required for the transition from scalar equations to operator equations.

---

## 00.1 Polynomial Rings and Divisibility

!!! definition "Definition 00.1 (Field)"
    A **field** $F$ is a set with two operations (addition and multiplication) satisfying commutativity, associativity, distributivity, and the existence of identities and inverses. Common fields include $\mathbb{Q}, \mathbb{R}, \mathbb{C}$ and finite fields $\mathbb{F}_q$.

!!! theorem "Theorem 00.1 (Division Algorithm)"
    For any $f(x), g(x) \in F[x]$ with $g(x) \neq 0$, there exist unique $q(x), r(x) \in F[x]$ such that $f(x) = g(x)q(x) + r(x)$ and $\deg r < \deg g$.

---

## Exercises

1. **[Fundamentals] Determine the GCD of $x^3 - 1$ and $x^2 - 1$ over $\mathbb{Q}$.**
   ??? success "Solution"
       Using the Euclidean algorithm: $(x^3 - 1) = x(x^2 - 1) + (x - 1)$. Then $(x^2 - 1) = (x + 1)(x - 1) + 0$. The $\gcd$ is $x-1$.

2. **[Irreducibility] Is $x^2 + 1$ irreducible over $\mathbb{R}$? What about over $\mathbb{C}$ and $\mathbb{F}_2$?**
   ??? success "Solution"
       - Over $\mathbb{R}$: Yes, no real roots.
       - Over $\mathbb{C}$: No, factors as $(x-i)(x+i)$.
       - Over $\mathbb{F}_2$: No, factors as $(x+1)^2$ since $1+1=0$.

3. **[Roots] State the Fundamental Theorem of Algebra.**
   ??? success "Solution"
       Every non-constant polynomial in $\mathbb{C}[x]$ has at least one complex root. Consequently, every $n$-th degree polynomial splits into $n$ linear factors over $\mathbb{C}$.

4. **[Vieta] Relate the coefficients of $x^2 + ax + b = 0$ to its roots $\lambda_1, \lambda_2$.**
   ??? success "Solution"
       $\lambda_1 + \lambda_2 = -a$ and $\lambda_1 \lambda_2 = b$. In linear algebra, these correspond to the trace and determinant of a $2 \times 2$ matrix.

5. **[Interpolation] Find the Lagrange interpolating polynomial passing through $(0, 1), (1, 2), (2, 5)$.**
   ??? success "Solution"
       The polynomial is $p(x) = x^2 + 1$. Interpolation theory ensures a unique polynomial of degree $< n$ exists for $n$ distinct points.

6. **[Minimal Polynomial] Define the minimal polynomial of a matrix $A$.**
   ??? success "Solution"
       The unique monic polynomial $m(x)$ of lowest degree such that $m(A) = 0$. It must divide every other polynomial $p(x)$ for which $p(A) = 0$.

7. **[Multiplicity] Distinguish between algebraic and geometric multiplicity of a root.**
   ??? success "Solution"
       Algebraic multiplicity is the power of $(x-\lambda)$ in the characteristic polynomial. Geometric multiplicity is the dimension of the corresponding null space (eigenspace).

8. **[Derivative] Show that $\lambda$ is a multiple root of $f(x)$ iff $f(\lambda)=0$ and $f'(\lambda)=0$.**
   ??? success "Solution"
       If $f(x) = (x-\lambda)^k g(x)$ with $k \ge 2$, then $f'(x) = k(x-\lambda)^{k-1} g(x) + (x-\lambda)^k g'(x)$. Evaluation at $\lambda$ gives 0 for both $f$ and $f'$.

9. **[Characteristic] Why do we define the characteristic polynomial as $\det(\lambda I - A)$?**
   ??? success "Solution"
       Roots of this polynomial are the values for which $(\lambda I - A)$ is singular, meaning there exists a non-zero vector $v$ such that $Av = \lambda v$.

10. **[Euclidean Domain] Why is $F[x]$ called a Euclidean Domain?**
    ??? success "Solution"
        Because it possesses a division algorithm with a degree function that satisfies the Euclidean requirements, allowing for the construction of GCDs and unique factorization.

## Chapter Summary

This chapter establishes polynomials as the governing equations of linear algebra:

1. **Algebraic Rigor**: Defined fields and polynomial rings as the environment for matrix coefficients.
2. **Structural Tools**: Developed the division and Euclidean algorithms for manipulating operator identities.
3. **Spectral Encoding**: Linked the factorization of polynomials to the decomposition of vector spaces into invariant subspaces.
4. **Analytic Link**: Demonstrated through interpolation and root theory how scalar properties translate to matrix-valued functions.
