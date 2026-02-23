# Chapter 47B: Fréchet Derivatives and Matrix Functions

<div class="context-flow" markdown>

**Prerequisites**: Matrix Functions (Ch13) · Matrix Calculus (Ch47A) · Norms (Ch15) · Linear Operators

**Chapter Outline**: Fréchet Derivative Definition → Condition Number of a Matrix Function → Kronecker Form of the Fréchet Derivative → The Daleckii-Krein Theorem (Divided Differences) → Estimating the Derivative via SVD → Complex Step Method → Applications in Stability and Sensitivity Analysis

**Extension**: Fréchet derivatives quantify the sensitivity of a matrix function (like $\exp(A)$ or $\log(A)$) to perturbations in the input matrix $A$.

</div>

While matrix calculus (Ch47A) often deals with scalar functions of matrices, **Fréchet derivatives** describe how a matrix function $f(A)$ itself changes when the input matrix $A$ is perturbed. For a non-linear function $f$, the Fréchet derivative $L_f(A, E)$ is a linear operator that maps the perturbation $E$ to the first-order change in $f(A)$. This is essential for computing the **condition number** of matrix functions and for analyzing the numerical stability of matrix algorithms.

---

## 47B.1 The Fréchet Derivative

!!! definition "Definition 47B.1 (Fréchet Derivative)"
    The Fréchet derivative of $f$ at $A$ is the linear operator $L_f(A, \cdot)$ such that:
    $$f(A + E) = f(A) + L_f(A, E) + o(\|E\|)$$
    The value $L_f(A, E)$ is the directional derivative of $f$ at $A$ in the direction $E$.

!!! theorem "Theorem 47B.1 (Daleckii-Krein Theorem)"
    If $A = Q \Lambda Q^*$ is diagonalizable, the entries of the Fréchet derivative in the eigenvector basis are:
    $$(Q^* L_f(A, E) Q)_{ij} = f[\lambda_i, \lambda_j] (Q^* E Q)_{ij}$$
    where $f[\lambda_i, \lambda_j]$ is the **divided difference**:
    $$f[\lambda_i, \lambda_j] = \begin{cases} \frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j} & 	ext{if } \lambda_i 
eq \lambda_j \ f'(\lambda_i) & 	ext{if } \lambda_i = \lambda_j \end{cases}$$

---

## Exercises

1. **[Fundamentals] Compute the Fréchet derivative of $f(A) = A^2$.**
   ??? success "Solution"
       $(A+E)^2 = A^2 + AE + EA + E^2$. The linear part in $E$ is $L_f(A, E) = AE + EA$.

2. **[Inverse] Compute the Fréchet derivative of $f(A) = A^{-1}$.**
   ??? success "Solution"
       From $d(A^{-1}) = -A^{-1}(dA)A^{-1}$, we have $L_f(A, E) = -A^{-1} E A^{-1}$.

3. **[Condition Number] Define the relative condition number $\kappa_f(A)$ of a matrix function.**
   ??? success "Solution"
       $\kappa_f(A) = \frac{\|L_f(A)\| \|A\|}{\|f(A)\|}$, where $\|L_f(A)\|$ is the operator norm of the Fréchet derivative. This measures how much a small relative change in $A$ is magnified in $f(A)$.

4. **[Commutativity] Does $L_f(A, E) = f'(A) E$?**
   ??? success "Solution"
       Only if $A$ and $E$ commute. In general, because $A$ and $E$ do not commute, the derivative is a more complex operator (like $AE + EA$ for $A^2$) rather than a simple multiplication.

5. **[Daleckii-Krein] Use the Daleckii-Krein theorem to find the Fréchet derivative of $f(A) = e^A$ for a diagonal matrix $A = \operatorname{diag}(\lambda_1, \lambda_2)$.**
   ??? success "Solution"
       $L_f(A, E) = \begin{pmatrix} e^{\lambda_1} E_{11} & \frac{e^{\lambda_1} - e^{\lambda_2}}{\lambda_1 - \lambda_2} E_{12} \ \frac{e^{\lambda_2} - e^{\lambda_1}}{\lambda_2 - \lambda_1} E_{21} & e^{\lambda_2} E_{22} \end{pmatrix}$.

6. **[Kronecker Form] Express $L_f(A, E)$ in the form $K_f(A) \operatorname{vec}(E)$.**
   ??? success "Solution"
       The operator $L_f(A, \cdot)$ can be represented by an $n^2 	imes n^2$ matrix $K_f(A)$. For $f(A) = A^2$, $K_f(A) = I \otimes A + A^T \otimes I$.

7. **[Exponential] What is the Fréchet derivative of the matrix exponential at the origin $A=0$?**
   ??? success "Solution"
       $L_{\exp}(0, E) = E$. Near the origin, the exponential map is the identity to first order.

8. **[Complex Step] Describe the complex step method for approximating $L_f(A, E)$.**
   ??? success "Solution"
       $L_f(A, E) \approx \operatorname{Im}(f(A + i h E) / h)$ for very small $h$. This method avoids the subtractive cancellation errors of standard finite differences.

9. **[Composition] State the chain rule for Fréchet derivatives.**
   ??? success "Solution"
       $L_{g \circ f}(A, E) = L_g(f(A), L_f(A, E))$. The derivative of the composition is the composition of the derivatives.

10. **[Symmetry] If $A$ and $E$ are Hermitian, is $L_f(A, E)$ always Hermitian?**
    ??? success "Solution"
        Yes, provided $f$ maps the real line to the real line. This is a consequence of the Daleckii-Krein formula and the symmetry of divided differences.

## Chapter Summary

This chapter explores the sensitivity and perturbation theory of matrix functions:

1. **Operator Calculus**: Defined the Fréchet derivative as the linear operator capturing the first-order response of matrix functions to perturbations.
2. **Spectral Sensitivity**: Utilized the Daleckii-Krein theorem to link the derivative to divided differences of eigenvalues.
3. **Stability Quantification**: Established the condition number as the definitive metric for the numerical reliability of matrix function evaluation.
4. **Computational Tools**: Outlined Kronecker representations and complex-step methods for the practical estimation of matrix derivatives.
