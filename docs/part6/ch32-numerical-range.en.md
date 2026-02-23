# Chapter 32: Numerical Range

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch6) · Inner Product (Ch8) · Convexity (Ch64A) · Hermitian Matrices (Ch16)

**Chapter Outline**: Definition of Numerical Range (Field of Values) → Toeplitz-Hausdorff Theorem (Convexity) → Spectral Containment → Numerical Radius → Properties of $W(A)$ → Boundary and Elliptical Range Theorem → Relation to Norms and Stability

**Extension**: The numerical range is used to bound the spectrum of non-normal operators and to analyze the stability of numerical methods for PDEs.

</div>

The **numerical range** (also known as the field of values) of a matrix $A$ is the set of all possible expected values of the operator when acting on unit vectors. While the spectrum $\sigma(A)$ is a discrete set of points, the numerical range $W(A)$ is a continuous, convex subset of the complex plane that "envelopes" the eigenvalues. It provides a more robust description of an operator's magnitude and stability than the spectrum alone, especially for non-normal matrices.

---

## 32.1 Definitions and Convexity

!!! definition "Definition 32.1 (Numerical Range)"
    The numerical range of an $n \times n$ complex matrix $A$ is the set:
    $$W(A) = \{ x^* A x : x \in \mathbb{C}^n, \|x\| = 1 \}$$

!!! theorem "Theorem 32.1 (Toeplitz-Hausdorff Theorem)"
    The numerical range $W(A)$ is a **convex** subset of the complex plane.

---

## Exercises

1. **[Fundamentals] Determine $W(A)$ for a Hermitian matrix $A$.**
   ??? success "Solution"
       For a Hermitian matrix, $x^* A x$ is always real. $W(A)$ is the closed interval on the real axis $[\lambda_{\min}(A), \lambda_{\max}(A)]$.

2. **[Spectral Containment] Prove that the spectrum $\sigma(A)$ is always contained within the numerical range $W(A)$.**
   ??? success "Solution"
       If $\lambda \in \sigma(A)$, there exists a unit eigenvector $v$ such that $Av = \lambda v$. Then $v^* A v = v^* (\lambda v) = \lambda \|v\|^2 = \lambda$. Thus $\lambda \in W(A)$.

3. **[Numerical Radius] Define the numerical radius $r(A)$ and its relationship to the spectral norm.**
   ??? success "Solution"
       $r(A) = \max \{ |z| : z \in W(A) \}$. It satisfies $\frac{1}{2}\|A\| \le r(A) \le \|A\|$. For normal matrices, $r(A) = \|A\| = \rho(A)$.

4. **[Normal Matrices] For a normal matrix $A$, show that $W(A)$ is the convex hull of its eigenvalues.**
   ??? success "Solution"
       If $A$ is normal, it is unitarily diagonalizable. $x^* A x = \sum \lambda_i |y_i|^2$ where $y = U^* x$. Since $\sum |y_i|^2 = 1$, $x^* A x$ is a convex combination of the eigenvalues.

5. **[Ellipse Theorem] Compute $W(A)$ for $A = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$.**
   ??? success "Solution"
       For this $2 \times 2$ matrix, $W(A)$ is a circular disk centered at the origin with radius $1/2$. More generally, for $2 \times 2$ matrices, $W(A)$ is always an ellipse with the eigenvalues as foci.

6. **[Shift and Scale] How do $W(A+cI)$ and $W(cA)$ relate to $W(A)$?**
   ??? success "Solution"
       $W(A+cI) = W(A) + c$ (translation) and $W(cA) = c W(A)$ (scaling and rotation).

7. **[Subadditivity] Is $W(A+B)$ equal to $W(A) + W(B)$?**
   ??? success "Solution"
       No, in general $W(A+B) \subseteq W(A) + W(B)$. Equality holds if $A$ and $B$ are scalars or commute in specific ways.

8. **[Hermitian Parts] Relate $W(A)$ to the numerical ranges of its Hermitian and skew-Hermitian parts.**
   ??? success "Solution"
       Let $H = (A+A^*)/2$ and $K = (A-A^*)/2$. Then $\operatorname{Re}(W(A)) = W(H)$ and $\operatorname{Im}(W(A)) = W(K/i)$.

9. **[Stability] If $W(A)$ lies entirely in the left half-plane, what can you say about the stability of $\dot{x} = Ax$?**
   ??? success "Solution"
       The system is dissipative and stable. Since $\sigma(A) \subseteq W(A)$, all eigenvalues have negative real parts. Furthermore, the norm $\|e^{At}\|$ is strictly decreasing.

10. **[Boundary] Describe the significance of the "sharp points" on the boundary of $W(A)$.**
    ??? success "Solution"
        Any sharp point (corner) of the boundary of $W(A)$ must be an eigenvalue of $A$. This is because the numerical range of a $2 \times 2$ sub-block is an ellipse, which is smooth; corners can only occur when multiple eigenspaces intersect at the boundary.

## Chapter Summary

This chapter explores the continuous spectral envelope of an operator:

1. **Convex Geometry**: Established the Toeplitz-Hausdorff theorem as the definitive property of the field of values.
2. **Spectral Bound**: Positioned the numerical range as a convex container for the spectrum, particularly useful for non-normal operators.
3. **Radial Metrics**: Defined the numerical radius as a stable alternative to the spectral radius for non-commuting sums.
4. **Norm Relations**: Linked the boundary of $W(A)$ to the growth of the matrix exponential and the numerical stability of linear systems.
