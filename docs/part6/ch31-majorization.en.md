# Chapter 31: Majorization and Schur-Horn Theorem

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Convexity (Ch64A) · Doubly Stochastic Matrices (Ch64A)

**Chapter Outline**: Definition of Majorization → Partial Order on $\mathbb{R}^n$ → Doubly Stochastic Matrices and Majorization (Hardy-Littlewood-Pólya) → Schur-Horn Theorem → Ky Fan's Maximum Principle → Lidskii's Theorem → Log-majorization → Applications in Probability and Physics

**Extension**: Majorization is the mathematical language of inequality and uncertainty (entropy) comparison; it connects the diagonal entries of a matrix to its eigenvalues.

</div>

Majorization provides a precise mathematical way to say that one vector is "more spread out" than another. It defines a partial order on vectors in $\mathbb{R}^n$ that is preserved by all doubly stochastic transformations. In matrix theory, majorization is the link between the entries of a matrix and its spectrum, as epitomized by the Schur-Horn theorem.

---

## 31.1 Definitions and Characterizations

!!! definition "Definition 31.1 (Majorization)"
    Let $x, y \in \mathbb{R}^n$. We say $x$ is **majorized** by $y$ (denoted $x \prec y$) if:
    $$\sum_{i=1}^k x_i^\downarrow \le \sum_{i=1}^k y_i^\downarrow, \quad k=1, \dots, n-1, \quad \text{and } \sum_{i=1}^n x_i = \sum_{i=1}^n y_i$$
    where $x_i^\downarrow$ are the components of $x$ sorted in non-increasing order.

!!! theorem "Theorem 31.1 (Hardy-Littlewood-Pólya)"
    $x \prec y$ if and only if $x = Dy$ for some doubly stochastic matrix $D$.

---

## Exercises

1. **[Fundamentals] Verify that $(1/2, 1/2) \prec (1, 0)$.**
   ??? success "Solution"
       Sums: $1/2 \le 1$ and $1/2+1/2 = 1+0 = 1$. The conditions hold. Geometrically, $(1/2, 1/2)$ is the average of $(1, 0)$ and $(0, 1)$, hence less "spread out."

2. **[Convexity] Prove that $x \prec y$ implies $\sum f(x_i) \le \sum f(y_i)$ for any convex function $f$.**
   ??? success "Solution"
       Since $x = Dy$, each $x_i$ is a convex combination of $y_j$. By Jensen's inequality, $f(x_i) \le \sum_j d_{ij} f(y_j)$. Summing over $i$ and using $\sum_i d_{ij} = 1$ yields the result.

3. **[Schur-Horn] State the Schur-Horn theorem.**
   ??? success "Solution"
       For a real symmetric matrix $A$, the vector of its diagonal entries $d$ is majorized by the vector of its eigenvalues $\lambda$: $d \prec \lambda$. Conversely, given any $d \prec \lambda$, there exists a symmetric matrix with these diagonals and eigenvalues.

4. **[Permutahedron] Describe the geometric set $\{x : x \prec y\}$.**
   ??? success "Solution"
       It is the **permutahedron** of $y$: the convex hull of all $n!$ permutations of the vector $y$. This is a direct consequence of Birkhoff's theorem on doubly stochastic matrices.

5. **[Ky Fan] How does Ky Fan's maximum principle relate to majorization?**
   ??? success "Solution"
       Ky Fan's principle states that $\sum_{i=1}^k \lambda_i(A) = \max \operatorname{tr}(U^T AU)$ where $U$ has $k$ orthonormal columns. This provides the variational upper bound needed to prove the majorization of diagonal entries.

6. **[Entropy] Prove that if $p \prec q$ are probability distributions, then the Shannon entropy $H(p) \ge H(q)$.**
   ??? success "Solution"
       The function $f(t) = t \log t$ is convex. $p \prec q \implies \sum p_i \log p_i \le \sum q_i \log q_i$, so $-H(p) \le -H(q)$, which means $H(p) \ge H(q)$. Majorization corresponds to a move toward the uniform distribution (maximum entropy).

7. **[Log-majorization] Define log-majorization $x \prec_{\log} y$ for positive vectors.**
   ??? success "Solution"
       $x \prec_{\log} y$ if $\prod_{i=1}^k x_i^\downarrow \le \prod_{i=1}^k y_i^\downarrow$ for all $k$, with equality at $k=n$. This is equivalent to $\log x \prec \log y$.

8. **[Singular Values] Relate majorization to the singular values of $A$ and $B$ when considering the product $AB$.**
   ??? success "Solution"
       The singular values satisfy $\sigma(AB) \prec_{\log} \sigma(A) \circ \sigma(B)$ (the entry-wise product of sorted singular values).

9. **[Lidskii's Theorem] State Lidskii's theorem regarding the eigenvalues of a sum of symmetric matrices.**
   ??? success "Solution"
       $\lambda(A+B) - \lambda(A) \prec \lambda(B)$. This provides a majorization-based constraint on how the spectrum shifts under additive perturbation.

10. **[Quantum] Explain the significance of majorization in quantum state transformation.**
    ??? success "Solution"
        A quantum state $\rho$ can be transformed into $\sigma$ using local operations and classical communication (LOCC) if and only if the vector of eigenvalues of $\rho$ is majorized by that of $\sigma$ ($spec(\rho) \prec spec(\sigma)$).

## Chapter Summary

This chapter formalizes the comparison of "spread" in vectors and its matrix implications:

1. **Ordering Diversity**: Defined majorization as a partial order characterizing the concentration of components.
2. **Matrix Diagonals**: Utilized the Schur-Horn theorem to establish that eigenvalues are the most "extreme" possible diagonal entries.
3. **Variational Principles**: Linked majorization to Ky Fan's trace maximization, providing a calculus for sums of eigenvalues.
4. **Information Metrics**: Demonstrated that majorization is the categorical tool for comparing entropy and uncertainty in probability and quantum theory.
