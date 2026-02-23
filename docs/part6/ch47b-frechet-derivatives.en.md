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
    $$f[\lambda_i, \lambda_j] = \begin{cases} \frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j} & \text{if } \lambda_i \neq \lambda_j \\ f'(\lambda_i) & \text{if } \lambda_i = \lambda_j \end{cases}$$

---

## Exercises

1. **[Fundamentals] Compute the Fréchet derivative of $f(A) = A^2$.**
   ??? success "Solution"
       Expand $(A+E)^2 = A^2 + AE + EA + E^2$. The part linear in $E$ is $L_f(A, E) = AE + EA$. This operator acts on $E$ by adding $A$ from both sides.

2. **[Inverse] Compute the Fréchet derivative of $f(A) = A^{-1}$.**
   ??? success "Solution"
       Using $d(A^{-1}) = -A^{-1}(dA)A^{-1}$, the Fréchet derivative is the operator $L_f(A, E) = -A^{-1} E A^{-1}$. This captures how the inverse shifts when the matrix is perturbed.

3. **[Condition Number] Define the relative condition number $\kappa_f(A)$ of a matrix function.**
   ??? success "Solution"
       $\kappa_f(A) = \frac{\|L_f(A)\| \|A\|}{\|f(A)\|}$, where $\|L_f(A)\|$ is the induced operator norm. It measures the sensitivity of the output function value to relative errors in the input matrix.

4. **[Commutativity] Under what condition does $L_f(A, E) = f'(A) E$?**
   ??? success "Solution"
       This holds if and only if $A$ and $E$ commute ($AE = EA$). Non-commutativity is the primary reason why matrix derivatives are operators rather than simple multipliers.

5. **[Daleckii-Krein] Use the Daleckii-Krein theorem to find $L_f(A, E)$ for $f(A) = e^A$ when $A = \operatorname{diag}(\lambda_1, \lambda_2)$.**
   ??? success "Solution"
       The entries of $L_f(A, E)$ are $L_{ij} = \frac{e^{\lambda_i} - e^{\lambda_j}}{\lambda_i - \lambda_j} E_{ij}$ for $i \neq j$, and $L_{ii} = e^{\lambda_i} E_{ii}$ for $i=j$. This formula links the sensitivity of the exponential to the spread of the eigenvalues.

6. **[Kronecker Form] Express the operator $L_f(A, E) = AE + EA$ in its $n^2 \times n^2$ Kronecker matrix form.**
   ??? success "Solution"
       The matrix representation is $K_f(A) = I \otimes A + A^T \otimes I$. Acting on $\operatorname{vec}(E)$ with this matrix yields $\operatorname{vec}(L_f(A, E))$.

7. **[Exponential Identity] What is the Fréchet derivative of the matrix exponential at the origin $A=0$?**
   ??? success "Solution"
       $L_{\exp}(0, E) = E$. To first order, the exponential map is the identity near zero.

8. **[Complex Step] Describe the advantage of the complex step method for estimating $L_f(A, E)$.**
   ??? success "Solution"
       $L_f(A, E) \approx \operatorname{Im}(f(A + i h E) / h)$. Because it doesn't involve subtracting two close numbers, it is immune to the subtractive cancellation that plagues standard finite difference methods.

9. **[Composition] State the chain rule for matrix Fréchet derivatives.**
   ??? success "Solution"
       For $h(A) = g(f(A))$, the derivative is the composition of the linear operators: $L_h(A, E) = L_g(f(A), L_f(A, E))$.

10. **[Symmetry] Prove that if $A$ and $E$ are Hermitian, then $L_f(A, E)$ is also Hermitian (assuming $f$ is real-valued on $\mathbb{R}$).**
    ??? success "Solution"
        From Daleckii-Krein, $(Q^* L Q)_{ij} = f[\lambda_i, \lambda_j] (Q^* E Q)_{ij}$. Since $f[\lambda_i, \lambda_j] = f[\lambda_j, \lambda_i]$ and $Q^*EQ$ is Hermitian, their product is Hermitian. Thus $L$ is Hermitian.

## Chapter Summary

This chapter explores the sensitivity and perturbation theory of matrix functions:

1. **Operator Calculus**: Defined the Fréchet derivative as the linear operator capturing the first-order response of matrix functions to perturbations.
2. **Spectral Sensitivity**: Utilized the Daleckii-Krein theorem to link the derivative to divided differences of eigenvalues.
3. **Stability Quantification**: Established the condition number as the definitive metric for the numerical reliability of matrix function evaluation.
4. **Computational Tools**: Outlined Kronecker representations and complex-step methods for the practical estimation of matrix derivatives.
