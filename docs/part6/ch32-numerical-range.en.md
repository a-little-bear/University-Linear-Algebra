# Chapter 32: Numerical Range and Numerical Radius

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch08) · Eigenvalues (Ch06) · Convex Sets (Ch25)

**Chapter Outline**: Definition of Numerical Range (Field of Values) → Toeplitz-Hausdorff Convexity Theorem → Spectral Inclusion Properties → Numerical Radius $w(A)$ → Numerical Radius Inequalities vs. Spectral Norm → The Ellipticity Theorem for 2x2 Matrices → Numerical Range of Normal Matrices → Applications in Operator Stability

**Extension**: The numerical range is the "continuous expansion" of eigenvalues in the complex plane; it plays a vital role in analyzing transient growth of non-normal operators, convergence of Krylov subspace methods, and the range of quantum observables.

</div>

If eigenvalues $\sigma(A)$ are the "discrete fingerprints" of a matrix in the complex plane, then the **Numerical Range** $W(A)$ is its "continuous shadow." The numerical range captures all possible values of the quadratic form on the unit sphere. For normal matrices, this shadow is merely the convex hull of the eigenvalues; however, for non-normal matrices, the numerical range exhibits complex geometric forms that reveal the underlying instability of the operator.

---

## 32.1 Definition and Convexity

!!! definition "Definition 32.1 (Numerical Range / Field of Values)"
    The **Numerical Range** of a square matrix $A \in M_n(\mathbb{C})$ is the set of all values taken by the quadratic form on the unit sphere:
    $$W(A) = \{ \mathbf{x}^* A \mathbf{x} : \mathbf{x} \in \mathbb{C}^n, \|\mathbf{x}\| = 1 \}$$

!!! theorem "Theorem 32.1 (Toeplitz-Hausdorff Theorem)"
    For any matrix $A$, its numerical range $W(A)$ is always a **convex set** in the complex plane.
    **Significance**: Despite the unit sphere being a curved manifold, the mapping of the quadratic form into the complex plane exhibits perfect convexity.

---

## 32.2 Spectral Inclusion and Normality

!!! theorem "Theorem 32.2 (Spectral Inclusion)"
    All eigenvalues of a matrix $A$ are contained within its numerical range: $\sigma(A) \subseteq W(A)$.
    Furthermore, $W(A)$ contains the convex hull of the eigenvalues, $\operatorname{conv}(\sigma(A))$.

!!! theorem "Theorem 32.3 (Numerical Range of Normal Matrices)"
    $A$ is a normal matrix ($AA^* = A^*A$) $\iff$ $W(A) = \operatorname{conv}(\sigma(A))$.
    For non-normal matrices, $W(A)$ is typically strictly larger than the convex hull of its spectrum.

---

## 32.3 Numerical Radius $w(A)$

!!! definition "Definition 32.2 (Numerical Radius)"
    The maximum modulus of the elements in the numerical range is called the **Numerical Radius**:
    $$w(A) = \sup \{ |z| : z \in W(A) \ \}$$

!!! theorem "Theorem 32.4 (Numerical Radius Inequalities)"
    The numerical radius and the spectral norm $\|A\|_2$ satisfy the following two-sided inequality:
    $$\frac{1}{2} \|A\|_2 \le w(A) \le \|A\|_2$$
    This shows that $w(A)$ is equivalent to the operator norm, but $w(A)$ often provides a more refined characterization when studying the convergence of powers.

---

## Exercises

1.  **[Basic] Calculate the numerical range of the identity matrix $I$.**
    ??? success "Solution"
        $x^* I x = x^* x = 1$ for all $\|x\|=1$. Thus $W(I) = \{1\}$, a single point in the complex plane.

2.  **[Normal] Find the numerical range of $A = \operatorname{diag}(1, i)$.**
    ??? success "Solution"
        Since $A$ is normal, $W(A)$ is the convex hull of eigenvalues 1 and $i$, which is the line segment connecting them.

3.  **[Convexity] Give an example showing that $\sigma(A)$ itself is not necessarily convex.**
    ??? success "Solution"
        $A = \operatorname{diag}(1, -1)$. The spectrum is $\{1, -1\}$, which does not contain 0. However, $W(A) = [-1, 1]$ contains 0.

4.  **[Radius] Find the numerical range of $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**
    ??? success "Solution"
        By the Ellipticity Theorem for 2x2 matrices, its numerical range is a closed disk centered at the origin with radius 0.5.

5.  **[Property] Prove $W(A+cI) = W(A) + c$ for any scalar $c$.**
    ??? success "Solution"
        $x^*(A+cI)x = x^*Ax + c(x^*x) = x^*Ax + c$.

6.  **[Self-adjoint] Prove: If $A$ is Hermitian, then $W(A)$ is the real interval $[\lambda_{\min}, \lambda_{\max}]$.**
    ??? success "Solution"
        The quadratic form $x^*Ax$ of a Hermitian matrix is always real, and its range is bounded by the minimum and maximum eigenvalues (Rayleigh Quotient property).

7.  **[Norm] If $\|A\|_2 = 10$, what is the minimum possible value for the numerical radius $w(A)$?**
    ??? success "Solution"
        From the inequality $w(A) \ge \frac{1}{2}\|A\|_2$, we have $w(A) \ge 5$.

8.  **[Trace] How is the centroid of the numerical range related to the trace?**
    ??? success "Solution"
        For an $n \times n$ matrix, the point $\frac{1}{n} \operatorname{tr}(A)$ is always in $W(A)$ (it is the average of the quadratic form values for any orthonormal basis).

9.  **[Unitary] Prove $W(U^* A U) = W(A)$ for any unitary matrix $U$.**
    ??? success "Solution"
        $(Ux)^* A (Ux) = y^* A y$ where $y=Ux$. Since $U$ is unitary, it maps the unit sphere onto itself, preserving the set of values.

10. **[Application] Why is the numerical range more reliable than eigenvalues for the stability analysis of non-normal operators?**
    ??? success "Solution"
        Eigenvalues only determine long-term (asymptotic) stability, while the size of the numerical range limits the short-term (transient) growth of the norm of the solution. Even if all eigenvalues are in the left half-plane, if $W(A)$ extends into the right half-plane, the system may experience violent transient oscillations.

## Chapter Summary

The numerical range is the geometric extension of matrix spectral theory:

1.  **Convexity of Form**: The Toeplitz-Hausdorff theorem reveals a deep consistency in the action of linear operators—regardless of the operator's complexity, its "projection" is always a convex region.
2.  **Divide of Normality**: The degree to which the numerical range coincides with the convex hull of eigenvalues is a direct geometric measure of a matrix's "normality" (whether it possesses an orthonormal basis of eigenvectors).
3.  **Analytic Bounds**: The numerical radius establishes a new measure of operator magnitude, providing tighter upper bounds for transient analysis and iterative convergence rates than the spectral radius.
