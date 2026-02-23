# Chapter 37B: Vandermonde, Cauchy Matrices and Displacement Structure

<div class="context-flow" markdown>

**Prerequisites**: Matrix Operations (Ch2) · Determinants (Ch3) · Decompositions (Ch10) · Numerical Linear Algebra (Ch22) · Toeplitz/Hankel (Ch37A)

**Chapter Outline**: Vandermonde Matrices & Interpolation → Cauchy Matrices & Displacement Rank 1 → Resultant & Krylov Matrices → Sylvester/Stein Displacement Operators → Calculation of Displacement Rank → Generalized Schur Algorithm → Gohberg-Semencul Formula → Hierarchical Matrices

**Extension**: Displacement structure theory (Kailath-Kung-Morf, 1979) unifies the fast algorithm frameworks for all classical structured matrix classes.

</div>

While Ch 37A focused on Toeplitz and Hankel matrices, this chapter explores Vandermonde and Cauchy matrices, characterized by power and kernel structures. We introduce the **Displacement Structure** framework, which reveals that while these matrices aren't low-rank, they become low-rank under specific "displacement operators." This insight allows for $O(n^2)$ or even $O(n \log^2 n)$ solvers.

---

## 37B.1 Core Concepts

!!! definition "Definition 37B.1 (Vandermonde Matrix)"
    Given nodes $x_1, \ldots, x_n \in \mathbb{C}$, the Vandermonde matrix $V$ is defined by $V_{ij} = x_i^{j-1}$. It is non-singular if and only if all $x_i$ are distinct.

!!! definition "Definition 37B.8 (Sylvester Displacement)"
    The displacement of matrix $A$ with respect to operators $(F, F')$ is $
abla_{F,F'}(A) = FA - AF'$. The **displacement rank** is the rank of this resulting matrix.

---

## Exercises

1. **[Vandermonde] Compute the determinant of $V(1, 2, 3)$.**

   ??? success "Solution"
       $\det V = (2-1)(3-1)(3-2) = 1 \cdot 2 \cdot 1 = 2$.

2. **[Condition Number] Why are Vandermonde matrices with equidistant nodes on $[0, 1]$ considered ill-conditioned?**

   ??? success "Solution"
       The condition number grows exponentially with $n$. High-degree polynomials on equidistant nodes lead to the Runge phenomenon, making the interpolation coefficients extremely sensitive to small perturbations in data.

3. **[Cauchy] Show that the displacement rank of a Cauchy matrix $C = (1/(x_i - y_j))$ is 1 with respect to diagonal operators.**

   ??? success "Solution"
       Let $D_x = \operatorname{diag}(x_i)$ and $D_y = \operatorname{diag}(y_j)$. Then $(D_x C - C D_y)_{ij} = x_i \frac{1}{x_i - y_j} - \frac{1}{x_i - y_j} y_j = \frac{x_i - y_j}{x_i - y_j} = 1$. The result is the all-ones matrix, which has rank 1.

4. **[Toeplitz] Prove that the inverse of a Toeplitz matrix is not necessarily Toeplitz, but it is "Toeplitz-like."**

   ??? success "Solution"
       While the symmetry of entries is lost, the displacement rank is preserved. If $T$ has displacement rank 2, then $T^{-1}$ also has displacement rank 2 (with respect to shifted operators), allowing for fast matrix-vector products.

5. **[Krylov] Relate the Krylov matrix $K(A, b)$ to the Vandermonde matrix when $A$ is diagonal.**

   ??? success "Solution"
       If $A = \operatorname{diag}(x_i)$, then $(A^k b)_i = x_i^k b_i$. Thus $K(A, b) = \operatorname{diag}(b) V(x)$, which is a row-scaled Vandermonde matrix.

6. **[Complexity] Compare the complexity of solving a general $n 	imes n$ system versus a Toeplitz system.**

   ??? success "Solution"
       General: $O(n^3)$ (Gaussian elimination). Toeplitz: $O(n^2)$ (Levinson-Durbin) or $O(n \log^2 n)$ (Superfast solvers).

7. **[Gohberg-Semencul] What is the significance of the Gohberg-Semencul formula?**

   ??? success "Solution"
       It expresses the inverse of a Toeplitz matrix as a difference of products of triangular Toeplitz matrices. This allows applying the inverse to a vector in $O(n \log n)$ time using FFT.

8. **[Resultant] When is the Sylvester resultant matrix $\operatorname{Syl}(f, g)$ singular?**

   ??? success "Solution"
       It is singular if and only if $f(x)$ and $g(x)$ share a common root, i.e., $\gcd(f, g) 
eq 1$.

9. **[H-Matrices] Define the core idea behind Hierarchical ($\mathcal{H}$) matrices.**

   ??? success "Solution"
       While the whole matrix may be full rank, its off-diagonal blocks representing "far-field" interactions are low-rank. This allows for $O(n \log n)$ arithmetic.

10. **[Schur Algorithm] How does the Generalized Schur Algorithm use displacement generators?**

   ??? success "Solution"
        It performs Gaussian elimination directly on the low-rank generators $(G, H)$ instead of the $n^2$ matrix entries, reducing the work to $O(rn^2)$ where $r$ is the displacement rank.

## Chapter Summary

Displacement structure is the unifying theory of fast linear algebra:

1. **Low-rank Shadow**: Structured matrices are characterized by having a low-rank image under displacement operators.
2. **Inheritance**: This structure is preserved under inversion and Schur complementation.
3. **Algorithmic Efficiency**: By operating on generators, we break the $O(n^3)$ barrier for wide classes of engineering problems.
