# Chapter 18: Matrix Inequalities

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Eigenvalues (Ch6) · Norms (Ch15) · Matrix Analysis (Ch14)

**Chapter Outline**: Löwner Partial Order $(\succeq)$ → Inequalities for PD Matrices (Weyl, Ky Fan) → Singular Value Inequalities → Norm Inequalities (Hadamard, Cauchy-Schwarz) → Trace Inequalities (Golden-Thompson) → Variational Characterization of Eigenvalue Sums → Operator Monotone and Operator Convex Functions

**Extension**: Matrix inequalities are the analytical foundation for convex optimization, information theory, and quantum computing, generalizing scalar order to complex operator manifolds.

</div>

Matrix inequalities study the partial order relations between square matrices based on semi-definiteness. Compared to scalar inequalities, the existence of off-diagonal entries requires the use of spectral decomposition or variational principles to establish these relations.

---

## 18.1 Löwner Order and Spectral Inequalities

!!! definition "Definition 18.1 (Löwner Partial Order)"
    Let $A, B$ be Hermitian matrices. We say $A$ is greater than or equal to $B$ in the Löwner order (denoted $A \succeq B$) if $A - B$ is a positive semi-definite matrix.

!!! theorem "Theorem 18.3 (Weyl's Inequalities)"
    Let $\lambda_i(A)$ be the eigenvalues in descending order. For Hermitian matrices $A, B$:
    $$\lambda_{i+j-1}(A+B) \le \lambda_i(A) + \lambda_j(B)$$

---

## Exercises

1. **[Basic Property] Prove: If $A \succeq B \succeq 0$, then for any matrix $C$, $C^* A C \succeq C^* B C$.**
   ??? success "Solution"
       Consider the quadratic form: $x^* (C^* A C - C^* B C) x = (Cx)^* (A-B) (Cx)$.
       Let $y = Cx$. Since $A-B \succeq 0$, then $y^* (A-B) y \ge 0$.
       Since this holds for all $x$, $C^* A C \succeq C^* B C$.

2. **[Spectral Monotonicity] If $A \succeq B$, prove $\lambda_i(A) \ge \lambda_i(B)$ for all $i$.**
   ??? success "Solution"
       This follows from the Courant-Fischer variational characterization (minimax principle). Since $x^T Ax \ge x^T Bx$ for any subspace, the maximum/minimum values over those subspaces must satisfy the same order.

3. **[Inversion] Let $A \succeq B \succ 0$. Prove $B^{-1} \succeq A^{-1}$.**
   ??? success "Solution"
       Use congruence transformation. $A \succeq B \implies B^{-1/2} A B^{-1/2} \succeq I$.
       Let $X = B^{-1/2} A B^{-1/2}$, then $X \succeq I \implies X^{-1} \preceq I$.
       Thus $(B^{-1/2} A B^{-1/2})^{-1} \preceq I \implies B^{1/2} A^{-1} B^{1/2} \preceq I$.
       Multiply by $B^{-1/2}$ on both sides to get $A^{-1} \preceq B^{-1}$.

4. **[Hadamard] State Hadamard's Inequality and give its geometric interpretation.**
   ??? success "Solution"
       $\det(A) \le \prod a_{ii}$ (for PD matrices). Geometrically, this means the volume of the parallelotope formed by column vectors is maximized when the vectors are orthogonal (diagonal matrix).

5. **[Trace Inequality] Is $\operatorname{tr}(AB) \le \operatorname{tr}(A) \operatorname{tr}(B)$ always true for $A, B \succeq 0$?**
   ??? success "Solution"
       No. Counterexample: $A = B = I_{2 \times 2}$. $\operatorname{tr}(I^2) = 2$, but $\operatorname{tr}(I)\operatorname{tr}(I) = 4$. While $2 \le 4$ holds, this is not a general rule. A correct upper bound is usually $\operatorname{tr}(AB) \le \lambda_{\max}(A) \operatorname{tr}(B)$.

6. **[Ky Fan] State the Ky Fan $k$-norm sum inequality.**
   ??? success "Solution"
       $\sum_{i=1}^k \lambda_i(A+B) \le \sum_{i=1}^k \lambda_i(A) + \sum_{i=1}^k \lambda_i(B)$. This reflects the subadditivity of the sum of the largest eigenvalues.

7. **[Golden-Thompson] State the Golden-Thompson Inequality.**
   ??? success "Solution"
       $\operatorname{tr}(e^{A+B}) \le \operatorname{tr}(e^A e^B)$ for all Hermitian matrices $A, B$. This measures free energy bounds in statistical mechanics.

8. **[Operator Monotone] Give an example of a function that is increasing on scalars but not in the operator order.**
   ??? success "Solution"
       $f(t) = t^2$. Although $a > b > 0 \implies a^2 > b^2$, there exist $A \succeq B \succeq 0$ such that $A^2 \nsucceq B^2$. Only specific "operator monotone" functions like $\sqrt{t}, \log t, 1/t$ preserve the Löwner order.

9. **[Fiedler] State Fiedler's inequality regarding eigenvalues of sums of symmetric matrices.**
   ??? success "Solution"
       $\lambda(A+B)$ is constrained by $\lambda(A) + \lambda(B)$ in the sense of Majorization.

10. **[Application] How are matrix inequalities used in Compressed Sensing?**
    ??? success "Solution"
        Used to prove the Restricted Isometry Property (RIP). By bounding the fluctuation of singular values of the sensing matrix, it ensures that the energy of sparse signals is preserved after dimensionality reduction, allowing for exact reconstruction.

## Chapter Summary

Matrix inequalities are the "soft constraint" theory of linear algebra:

1. **Order Generalization**: Löwner order extends linear comparison from the real line to the PSD cone.
2. **Variational Essence**: Behind every spectral inequality lies the shadow of an extremum optimization problem.
3. **Analytical Rigidity**: Inequalities establish structural evolution boundaries for matrices under addition, multiplication, and functional mapping.
