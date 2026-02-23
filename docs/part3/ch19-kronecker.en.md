# Chapter 19: Kronecker Product and Vec Operator

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Matrix Equations (Ch20) · Multilinear Algebra & Tensors (Ch21)

**Chapter Outline**: From Standard Multiplication to Tensor Products → Definition of the Kronecker Product and its Block Structure → Key Algebraic Properties (Associativity, Transpose, Inverse) → The Mixed-Product Property → Spectral Properties: Eigenvalues and Trace of Kronecker Products → The Vec Operator and its Linearity → The Fundamental Identity: $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$ → The Kronecker Sum ($A \oplus B$) → Applications: Vectorized Solutions to Linear Matrix Equations (Sylvester, Lyapunov) and Modeling of Multivariate Systems

**Extension**: The Kronecker product is the concrete realization of the tensor product in the category of matrices; by performing "dimension multiplication," it constructs an algebraic space describing the coupling of multiple physical quantities, serving as an essential tool for studying Quantum Entanglement (Ch28) and high-dimensional signal processing.

</div>

When dealing with complex equations involving the interaction of multiple matrices (such as $AX + XB = C$), traditional matrix arithmetic often fails to provide a direct closed-form solution. The **Kronecker Product** and the **Vec Operator** provide a powerful toolkit for transforming "matrix equations" into standard "vector equations." This strategy of "dimensional reduction" allows us to solve high-dimensional operator interactions using classic linear systems theory.

---

## 19.1 Definition and Properties of the Kronecker Product

!!! definition "Definition 19.1 (Kronecker Product)"
    Let $A$ be an $m \times n$ matrix and $B$ be a $p \times q$ matrix. Their **Kronecker product** $A \otimes B$ is an $mp \times nq$ block matrix:
    $$A \otimes B = \begin{pmatrix} a_{11}B & a_{12}B & \cdots & a_{1n}B \\ a_{21}B & a_{22}B & \cdots & a_{2n}B \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1}B & a_{m2}B & \cdots & a_{mn}B \end{pmatrix}$$

!!! theorem "Theorem 19.1 (Mixed-Product Property)"
    If the matrix dimensions are compatible, then:
    $$(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$$
    **Significance**: This property allows us to deconstruct a complex composite transformation into two independent transformations in lower-dimensional spaces.

---

## 19.2 Spectral Properties and Trace

!!! theorem "Theorem 19.2 (Eigenvalues and Trace)"
    1.  **Eigenvalues**: If $A$ has eigenvalues $\{\lambda_i\}$ and $B$ has eigenvalues $\{\mu_j\}$, then $A \otimes B$ has eigenvalues $\{\lambda_i \mu_j\}$ for all possible pairs.
    2.  **Trace**: $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$.
    3.  **Determinant**: $\det(A \otimes B) = (\det A)^p (\det B)^m$ (where $B$ is $p \times p$ and $A$ is $m \times m$).

---

## 19.3 The Vec Operator and Vectorization Identity

!!! definition "Definition 19.2 (Vec Operator)"
    For an $m \times n$ matrix $X$, $\operatorname{vec}(X)$ is the $mn \times 1$ vector obtained by stacking the columns of $X$ one below the other in order.

!!! theorem "Theorem 19.3 (Vectorization Identity)"
    For matrices $A, X, B$ of appropriate dimensions:
    $$\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$$
    **Significance**: This is the "golden key" for solving matrix equations. It successfully strips the unknown matrix $X$ out of its surroundings, transforming the problem into a standard linear system.

---

## Exercises

**1. [Calculation] Compute $\begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} \otimes \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$.**

??? success "Solution"
    **Steps:**
    By definition, multiply each element of the left matrix by the right matrix:
    $$A \otimes B = \begin{pmatrix} 1 \cdot B & 0 \cdot B \\ 0 \cdot B & 2 \cdot B \end{pmatrix} = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 2 \\ 0 & 0 & 2 & 0 \end{pmatrix}$$

**2. [Eigenvalues] If $A$ has eigenvalues 1 and 2, and $B$ has eigenvalues 3 and 4, find all eigenvalues of $A \otimes B$.**

??? success "Solution"
    **Theorem Application:**
    The eigenvalues of a Kronecker product are the products of every possible pair of eigenvalues from the factors.
    1. $1 \cdot 3 = 3$
    2. $1 \cdot 4 = 4$
    3. $2 \cdot 3 = 6$
    4. $2 \cdot 4 = 8$
    **Conclusion**: The eigenvalues are $\{3, 4, 6, 8\}$.

**3. [Kronecker Sum] Find the eigenvalues of $A \oplus B = A \otimes I + I \otimes B$ for the matrices in the previous problem.**

??? success "Solution"
    **Theorem Application:**
    The eigenvalues of a Kronecker sum are the **sums** of the eigenvalue pairs.
    1. $1 + 3 = 4$
    2. $1 + 4 = 5$
    3. $2 + 3 = 5$
    4. $2 + 4 = 6$
    **Conclusion**: The eigenvalues are $\{4, 5, 5, 6\}$. This is often used in spectral analysis for solving $AX + XB = C$.

**4. [Vectorization] Transform the matrix equation $AX = B$ into the standard form $My = f$ (where $y = \operatorname{vec}(X)$).**

??? success "Solution"
    **Derivation:**
    1. The equation can be written as $AXI = B$ (where $I$ is an identity matrix of appropriate size).
    2. Apply the identity $\operatorname{vec}(AXB) = (B^T \otimes A) \operatorname{vec}(X)$.
    3. Set $B$ (in the formula) $= I$, $X$ $= X$, and $A$ $= A$.
    **Conclusion**: $(I \otimes A) \operatorname{vec}(X) = \operatorname{vec}(B)$.

**5. [Trace] Prove $\operatorname{tr}(A \otimes B) = \operatorname{tr}(B \otimes A)$.**

??? success "Solution"
    **Proof:**
    1. $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A)\operatorname{tr}(B)$.
    2. Since scalar multiplication is commutative, $\operatorname{tr}(A)\operatorname{tr}(B) = \operatorname{tr}(B)\operatorname{tr}(A)$.
    3. $\operatorname{tr}(B \otimes A) = \operatorname{tr}(B)\operatorname{tr}(A)$.
    **Conclusion**: While $A \otimes B \neq B \otimes A$ in general, their traces are identical.

**6. [Inverse] If $A$ and $B$ are both invertible, find $(A \otimes B)^{-1}$.**

??? success "Solution"
    **Using Mixed-Product Property:**
    1. Propose $A^{-1} \otimes B^{-1}$ as the inverse.
    2. Multiply: $(A \otimes B)(A^{-1} \otimes B^{-1}) = (A A^{-1}) \otimes (B B^{-1})$.
    3. $= I \otimes I = I$.
    **Conclusion**: $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$.

**7. [Rank] Prove $\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$.**

??? success "Solution"
    **Proof Strategy:**
    1. Use SVD: Let $A = U_1 \Sigma_1 V_1^*$ and $B = U_2 \Sigma_2 V_2^*$.
    2. Then $A \otimes B = (U_1 \otimes U_2) (\Sigma_1 \otimes \Sigma_2) (V_1 \otimes V_2)^*$.
    3. Since $U_1 \otimes U_2$ and $V_1 \otimes V_2$ remain unitary, this forms the SVD of $A \otimes B$.
    4. The singular value matrix is $\Sigma_1 \otimes \Sigma_2$, and the number of non-zero entries is clearly $\operatorname{rank}(A) \cdot \operatorname{rank}(B)$.

**8. [Vec Operation] For $X = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$, write $\operatorname{vec}(X)$.**

??? success "Solution"
    **Steps:**
    1. Extract first column: $(a, c)^T$.
    2. Extract second column: $(b, d)^T$.
    3. Stack vertically.
    **Conclusion**: $\operatorname{vec}(X) = (a, c, b, d)^T$. Note: Stacking is column-wise, not row-wise.

**9. [Lyapunov] Vectorize the Lyapunov equation $AX + XA^T = Q$.**

??? success "Solution"
    **Steps:**
    1. Vectorize each term: $\operatorname{vec}(AXI) + \operatorname{vec}(IXA^T) = \operatorname{vec}(Q)$.
    2. Apply formula: $(I \otimes A) \operatorname{vec}(X) + (A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(Q)$.
    3. Factor out $\operatorname{vec}(X)$: $(I \otimes A + A \otimes I) \operatorname{vec}(X) = \operatorname{vec}(Q)$.
    **Conclusion**: $(A \oplus A) \operatorname{vec}(X) = \operatorname{vec}(Q)$.

**10. [Application] Why are joint states of two particles in quantum mechanics represented by a tensor product (Kronecker product)?**

??? success "Solution"
    **Physical Logic:**
    1. Each particle's state is described by a vector space.
    2. The degrees of freedom of a joint system are the combination of the subsystems' degrees of freedom.
    3. If system 1 has $n$ basis states and system 2 has $m$ basis states, the joint system has $n \times m$ basis states.
    **Algebraic Mapping**: The Kronecker product constructs exactly such an $nm$-dimensional space through a "complete permutation" of basis pairs, perfectly characterizing both independence and **entanglement** (states that cannot be factored into $v_1 \otimes v_2$).

## Chapter Summary

The Kronecker product and Vec operator provide a scheme for "dimension elevation and reduction" in matrix algebra:

1.  **Dimension Multiplication**: The Kronecker product integrates the actions of two independent operators into a single massive composite operator via a tiling-and-nesting approach—the only language for multi-body interactions.
2.  **Operator Deconstruction**: The Vec operator eliminates the two-dimensional topology of a matrix, reducing it to its most basic vector form in linear space, thereby unleashing the full power of classical linear solvers.
3.  **Equation Unification**: The vectorization identity bridges matrix equation theory and numerical linear algebra, proving that all linear matrix equations are essentially the same linear system viewed under different bases.
