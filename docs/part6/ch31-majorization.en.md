# Chapter 31: Majorization and Doubly Stochastic Matrices

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Convex Sets (Ch25) · Matrix Inequalities (Ch18)

**Chapter Outline**: Geometric Definition of Majorization → Hardy-Littlewood-Pólya Theorem → Doubly Stochastic Matrices → Birkhoff-von Neumann Theorem (Convex Hull of Permutations) → Schur-Horn Theorem (Link between Spectra and Diagonals) → Schur-Convexity → Majorization of Entropy → Applications in Quantum Information and Inequality Metrics

**Extension**: Majorization is the unified mathematical framework for studying the "evenness" of distributions; it is central to pure state transformations in quantum mechanics and the Gini index in statistics.

</div>

In linear algebra, we frequently need to compare how "spread out" or "disordered" two vectors are. **Majorization theory** provides a robust mathematical framework for this comparison. It not only links the eigenvalues of a matrix to its diagonal entries (Schur-Horn Theorem) but also connects these relationships to convex geometry via **Doubly Stochastic Matrices**. This chapter reveals the deep structural laws governing the distribution of numerical data within matrices.

---

## 31.1 Definition of Majorization ($\prec$)

!!! definition "Definition 31.1 (Majorization)"
    Let $x, y \in \mathbb{R}^n$. Sort their components in non-increasing order as $x_{[1]} \ge x_{[2]} \ge \cdots \ge x_{[n]}$.
    We say $y$ **majorizes** $x$ (written $x \prec y$) if:
    1.  For $k=1, \ldots, n-1$, $\sum_{i=1}^k x_{[i]} \le \sum_{i=1}^k y_{[i]}$.
    2.  The total sums are equal: $\sum_{i=1}^n x_i = \sum_{i=1}^n y_i$.

!!! intuition "Intuition"
    $x \prec y$ means $x$ is "more even" or "less concentrated" than $y$. For example, for any non-negative vector $x$ summing to 1, $(1/n, \ldots, 1/n) \prec x \prec (1, 0, \ldots, 0)$.

---

## 31.2 Doubly Stochastic Matrices and Birkhoff's Theorem

!!! definition "Definition 31.2 (Doubly Stochastic Matrix)"
    A matrix $P \in M_n(\mathbb{R})$ is **doubly stochastic** if its entries are non-negative and every row and column sums to 1.

!!! theorem "Theorem 31.1 (Birkhoff-von Neumann Theorem)"
    The set of all $n \times n$ doubly stochastic matrices $\Omega_n$ is the **convex hull** of the set of all $n \times n$ **permutation matrices**.
    Thus, any doubly stochastic matrix can be written as $P = \sum \alpha_i P_{\sigma_i}$ where $\alpha_i \ge 0$ and $\sum \alpha_i = 1$.

!!! theorem "Theorem 31.2 (Hardy-Littlewood-Pólya)"
    For $x, y \in \mathbb{R}^n$, $x \prec y$ if and only if there exists a doubly stochastic matrix $P$ such that $x = Py$.

---

## 31.3 The Schur-Horn Theorem

!!! theorem "Theorem 31.3 (Schur-Horn Theorem)"
    Let $A$ be an $n \times n$ Hermitian matrix. Let $\mathbf{d} = (a_{11}, \ldots, a_{nn})$ be the vector of its diagonal entries, and $\boldsymbol{\lambda} = (\lambda_1, \ldots, \lambda_n)$ be the vector of its eigenvalues. Then:
    $$\mathbf{d} \prec \boldsymbol{\lambda}$$
    Conversely, if $\mathbf{d} \prec \boldsymbol{\lambda}$, there exists a Hermitian matrix with diagonal $\mathbf{d}$ and spectrum $\boldsymbol{\lambda}$.

---

## 31.4 Schur-Convexity

!!! definition "Definition 31.3 (Schur-Convex Function)"
    A function $\phi: \mathbb{R}^n \to \mathbb{R}$ is **Schur-convex** if $x \prec y \Rightarrow \phi(x) \le \phi(y)$.
    **Examples**: $\phi(x) = \sum x_i^2$ and $\phi(x) = \sum x_i \log x_i$ are Schur-convex.

---

## Exercises

1.  **[Fundamentals] Determine if $(2, 2, 2) \prec (3, 2, 1)$ holds.**
    ??? success "Solution"
        Yes. Both sum to 6. Partial sums: $2 < 3$ and $2+2=4 \le 3+2=5$. The conditions are met.

2.  **[Birkhoff] Decompose $\begin{pmatrix} 0.5 & 0.5 \\ 0.5 & 0.5 \end{pmatrix}$ into a combination of permutation matrices.**
    ??? success "Solution"
        $0.5 \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + 0.5 \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.

3.  **[Spectrum] Prove the largest eigenvalue of a doubly stochastic matrix is 1.**
    ??? success "Solution"
        Since $P \mathbf{1} = \mathbf{1}$, 1 is an eigenvalue. By Perron-Frobenius theory, the spectral radius of a non-negative matrix with row sums of 1 is 1.

4.  **[Schur-Horn] If a symmetric matrix has diagonal entries $(5, 5)$, can its eigenvalues be $(10, 0)$?**
    ??? success "Solution"
        Yes, because $(5, 5) \prec (10, 0)$ (sums match and $5 \le 10$).

5.  **[Entropy] Prove Shannon entropy $H(p) = -\sum p_i \log p_i$ is Schur-concave.**
    ??? success "Solution"
        Since $x \log x$ is convex, $\sum p_i \log p_i$ is Schur-convex. Negating it makes it Schur-concave. This implies that more uniform distributions (lower majorization) have higher entropy.

6.  **[Calculation] Find the vector $x$ from $\begin{pmatrix} 1/3 & 2/3 \\ 2/3 & 1/3 \end{pmatrix} \begin{pmatrix} 3 \\ 0 \end{pmatrix}$ and verify $x \prec (3, 0)^T$.**
    ??? success "Solution"
        $x = (1, 2)^T$. Sorted: $(2, 1) \prec (3, 0)$ because $2 \le 3$ and sums match.

7.  **[Extrema] Among all non-negative vectors summing to 1, which is the smallest under the majorization order?**
    ??? success "Solution"
        The uniform distribution $(1/n, \ldots, 1/n)$.

8.  **[Trace] Prove: If $A \ge 0$ is doubly stochastic, then $\operatorname{tr}(A) \le n$.**
    ??? success "Solution"
        Each diagonal entry $a_{ii} \le 1$ (since row sums are 1 and entries are non-negative), so the sum is $\le n$.

9.  **[Convexity] Is $\phi(x) = \max(x_i)$ Schur-convex?**
    ??? success "Solution"
        Yes. $x \prec y$ means $x$ is smoother than $y$, so its maximum value cannot exceed the maximum value of $y$.

10.  **[Application] Why is majorization relevant in wealth inequality?**
    ??? success "Solution"
        It quantifies the "inequality" between datasets. Using Lorenz curves and majorization, one can rigorously determine if the wealth gap in a population is expanding or shrinking.

## Chapter Summary

Majorization theory establishes the ordering of distributional shapes:

1.  **Metric of Uniformity**: Provided a rigorous mathematical characterization of "disorder" and "averaging," proving that the uniform distribution is the "floor" (minimum element) of all distributions.
2.  **Bridge between Operators and Geometry**: Theorems by Birkhoff and Schur-Horn showed how discrete structures (permutations) support continuous convex spaces and the essential constraints between a matrix's trace and its spectrum.
3.  **Information Foundations**: Through Schur-convexity, majorization became the natural language for processing entropy, energy dissipation, and quantum decoherence.
