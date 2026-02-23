# Chapter 47B: Fréchet Derivatives and High-order Theory

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Matrix Functions (Ch13) · Matrix Calculus Basics (Ch47A) · Norms (Ch15)

**Chapter Outline**: Definition of Fréchet Derivatives (Linear Approximation) → Relation to Gâteaux Derivatives → Derivatives of Matrix Functions (Daleckii-Krein Theorem) → Integral Representations → Derivatives of Eigenvalues and Eigenvectors → Derivatives of Singular Value Decomposition (SVD) → Higher-order Derivatives and Taylor Expansions → Applications: Sensitivity Analysis and Newton’s Method for Nonlinear Matrix Equations

**Extension**: The Fréchet derivative elevates matrix calculus to the level of operator analysis; it is the most rigorous tool for studying tangent spaces on matrix manifolds (Ch24) and quantifying the impact of perturbations in numerical stability analysis.

</div>

In elementary matrix calculus, we focus on the arrangement of numerical values. In high-order theory, we view a matrix function $f(A)$ as a mapping between operator spaces. The **Fréchet Derivative** provides a coordinate-independent definition of a derivative in a pure operator sense. It not only handles infinitesimal changes in matrices but also reveals the fine evolution of non-linear invariants such as eigenvalues and singular values under perturbation.

---

## 47B.1 Definition of the Fréchet Derivative

!!! definition "Definition 47B.1 (Fréchet Derivative)"
    Let $f: M_n 	o M_n$ be a mapping on matrix space. If there exists a linear operator $L_f(A, \cdot)$ such that:
    $$f(A+E) = f(A) + L_f(A, E) + o(\|E\|)$$
    then $L_f(A, E)$ is called the **Fréchet derivative** of $f$ at $A$ in the direction $E$.

!!! note "Relation to Gâteaux Derivatives"
    The Gâteaux derivative is the directional derivative along $E$: $G_f(A, E) = \frac{d}{dt} f(A+tE) |_{t=0}$.
    If the Fréchet derivative exists, the two are equal. The Fréchet derivative requires the approximation to be uniform across all directions.

---

## 47B.2 Derivatives of Matrix Functions and Daleckii-Krein

!!! theorem "Theorem 47B.1 (Daleckii-Krein Formula)"
    Let $A = V \Lambda V^{-1}$ be diagonalizable. The derivative $L = L_f(A, E)$ of the matrix function $f(A)$, expressed in the basis of eigenvectors, is:
    $$(V^{-1} L V)_{ij} = \frac{f(\lambda_i) - f(\lambda_j)}{\lambda_i - \lambda_j} (V^{-1} E V)_{ij}$$
    where the quotient is taken as $f'(\lambda_i)$ if $\lambda_i = \lambda_j$.
    **Significance**: This formula reveals that the variation of a matrix function is closely linked to the divided differences of the original matrix's eigenvalues (Löwner Matrix, Ch46).

---

## 47B.3 Derivatives of Eigenvalues and Singular Values

!!! theorem "Theorem 47B.2 (Derivative of Eigenvalues)"
    If $\lambda$ is a simple eigenvalue of $A$, and $x, y$ are its corresponding right and left eigenvectors (normalized such that $y^* x = 1$), then:
    $$\frac{d\lambda}{dA}(E) = y^* E x$$
    **Physical Meaning**: The change in an eigenvalue is the "projection" of the perturbation onto the eigen-direction.

---

## Exercises


****
??? success "Solution"
     $(A+E)^2 - A^2 = AE + EA + E^2 \approx AE + EA$. Thus $L = AE + EA$.


****
??? success "Solution"
     $(A+E)^{-1} - A^{-1} = (A(I+A^{-1}E))^{-1} - A^{-1} \approx (I - A^{-1}E)A^{-1} - A^{-1} = -A^{-1}EA^{-1}$.


****
??? success "Solution"
     The divided difference matrix is $\begin{pmatrix} 3 & 7 \ 7 & 12 \end{pmatrix}$. Thus $L = \begin{pmatrix} 3e_{11} & 7e_{12} \ 7e_{21} & 12e_{22} \end{pmatrix}$.


****
??? success "Solution"
     $x = e_1, y = e_1$. $\Delta \lambda \approx e_1^* E e_1 = E_{11} = 0.1$.


****
??? success "Solution"
     By Jacobi's formula, the first-order term is $\operatorname{tr}((
abla \det) E) = \operatorname{tr}(\det(A) A^{-T} E) = \det(A) \operatorname{tr}(A^{-1}E)$.


****
??? success "Solution"
     If $A = U\Sigma V^*$, then $\dot{\sigma} = u^* \dot{A} v$.


****
??? success "Solution"
     $L_f(A, E) = \int_0^1 e^{A(1-s)} E e^{As} ds$.


****
??? success "Solution"
     $L^{(2)}(A, E, H) = EH + HE$.


****
??? success "Solution"
     No. Since the trace is linear, its derivative is always $E \mapsto \operatorname{tr}(E)$.

****
??? success "Solution"
    ## Chapter Summary

Fréchet derivatives establish the rigorous logic of operator variation:


****: They flatten complex matrix functions locally, allowing us to describe non-linear evolution using the language of linear operators.

****: The Daleckii-Krein theorem reveals how the clustering of eigenvalues directly amplifies calculation errors via the divided difference effect—a theoretical early warning for numerical analysis.

****: Derivative formulas for eigenvalues and singular values provide precise analytical paths for understanding how complex systems (e.g., structural mechanics, neural networks) respond to fine parameter tuning.
