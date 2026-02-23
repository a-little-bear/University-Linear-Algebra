# Chapter 18: Matrix Inequalities

<div class="context-flow" markdown>

**Prerequisites**: Positive Definiteness (Ch16) · Matrix Norms (Ch15) · Trace (Ch2) · Majorization (Ch31)

**Chapter Outline**: Lowner Partial Order ($A \preceq B$) → Cauchy-Schwarz for Matrices → Trace Inequalities (Hölder, Cauchy) → Minkowski Determinant Inequality → Kantorovich Inequality → Hadamard's Inequality → Weyl's Inequality for Eigenvalues → Golden-Thompson Inequality

**Extension**: Matrix inequalities are the bounding logic of optimization; they provide the limits on information capacity, energy dissipation, and computational error.

</div>

Matrix inequalities generalize scalar inequalities to linear operators. The foundation is the **Lowner partial order**, which compares matrices based on the positive definiteness of their difference. This chapter explores profound results like **Hadamard's Inequality** (bounding the determinant by the product of column lengths) and the **Golden-Thompson Inequality** (bounding the trace of matrix exponentials). These bounds are essential for information theory, control systems, and the analysis of non-linear operators.

---

## 18.1 The Lowner Order and Basic Bounds

!!! definition "Definition 18.1 (Lowner Partial Order)"
    Let $A, B$ be symmetric matrices. We say $A \preceq B$ if $B - A$ is positive semi-definite ($B - A \succeq 0$). This defines a partial order on $S_n$.

!!! theorem "Theorem 18.1 (Hadamard's Inequality)"
    For any $n \times n$ matrix $A = (a_1, \dots, a_n)$:
    $$|\det A| \le \prod_{i=1}^n \|a_i\|_2$$
    Volume is maximized when the column vectors are orthogonal.

---

## Exercises

1. **[Fundamentals] If $A \preceq B$, prove that $\operatorname{tr}(A) \le \operatorname{tr}(B)$.**
   ??? success "Solution"
       $B - A \succeq 0 \implies \operatorname{tr}(B - A) \ge 0$. By linearity of the trace, $\operatorname{tr}(B) - \operatorname{tr}(A) \ge 0$, so $\operatorname{tr}(A) \le \operatorname{tr}(B)$. Trace is a monotonic function in the Lowner order.

2. **[Monotonicity] Is $A \preceq B$ imply $A^2 \preceq B^2$?**
   ??? success "Solution"
       No. The map $X \mapsto X^2$ is **not** operator monotone. However, $X \mapsto \sqrt{X}$ and $X \mapsto -X^{-1}$ are operator monotone (Löwner's Theorem).

3. **[Hadamard] Calculate the Hadamard bound for $\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       $\|a_1\| = 1$, $\|a_2\| = \sqrt{2}$. Bound is $1 \cdot \sqrt{2} \approx 1.414$. Actual determinant is 1.

4. **[Cauchy-Schwarz] State the trace version of the Cauchy-Schwarz inequality.**
   ??? success "Solution"
       $|\operatorname{tr}(A^* B)|^2 \le \operatorname{tr}(A^* A) \operatorname{tr}(B^* B)$. This is exactly the Cauchy-Schwarz inequality for the Hilbert-Schmidt (Frobenius) inner product.

5. **[Weyl] What does Weyl's Inequality say about the eigenvalues of $A+B$?**
   ??? success "Solution"
       $\lambda_i(A+B) \le \lambda_j(A) + \lambda_k(B)$ for $i = j+k-n$. It bounds the spectrum of a sum by the sums of individual spectra.

6. **[Kantorovich] What does the Kantorovich inequality bound?**
   ??? success "Solution"
       The ratio $\frac{(x^T Ax)(x^T A^{-1}x)}{(x^T x)^2}$. It provides an upper bound on how much the "spread" of eigenvalues can amplify a quadratic form, often used in the analysis of the steepest descent method.

7. **[Golden-Thompson] State the Golden-Thompson inequality.**
   ??? success "Solution"
       $\operatorname{tr}(e^{A+B}) \le \operatorname{tr}(e^A e^B)$. This is a fundamental result in statistical mechanics, though note that the operator inequality $e^{A+B} \preceq e^A e^B$ is false.

8. **[Determinant] Prove Minkowski's determinant inequality for $A, B \succeq 0$.**
   ??? success "Solution"
       $(\det(A+B))^{1/n} \ge (\det A)^{1/n} + (\det B)^{1/n}$. The $n$-th root of the volume is a concave function on the PSD cone.

9. **[Fischer] What is Fischer's Inequality for partitioned matrices?**
   ??? success "Solution"
       For a PSD matrix $M = \begin{pmatrix} A & B \\ B^* & D \end{pmatrix}$, $\det M \le \det A \cdot \det D$. This generalizes Hadamard's inequality to blocks.

10. **[Entropy] How does matrix inequality relate to the Von Neumann entropy $S(\rho) = -\operatorname{tr}(\rho \log \rho)$?**
    ??? success "Solution"
        Properties like the subadditivity of entropy follow from the operator concavity of $f(x) = -x \log x$ (Lieb's Theorem). Matrix inequalities provide the bounds for information processing.

## Chapter Summary

This chapter establishes the relational calculus of linear operators:

1. **Partial Ordering**: Defined the Lowner order as the standard for comparing the "energy" of matrices.
2. **Volume Constraints**: Utilized Hadamard and Minkowski inequalities to establish the limits of geometric scaling.
3. **Spectral Bounds**: Formulated Weyl's inequalities to track the migration of eigenvalues under operator sums.
4. **Analytic Limits**: Leveraged trace and exponential inequalities to provide the foundations for information and statistical theory.
