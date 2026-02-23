# Chapter 52: Linear Algebra over Finite Fields

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Polynomial Algebra (Ch00) · Abstract Algebra Basics

**Chapter Outline**: From Continuous to Discrete Fields → Definition and Construction of Finite Fields $GF(q)$ → Characteristic of a Field → Vector Spaces over Finite Fields → Counting Principles for Dimensions and Basis Vectors → Matrix Rank and RREF → Order of the General Linear Group $GL(n, q)$ → Applications: Error-Correcting Codes (Linear Codes, Hamming Codes), Cryptography (Field Arithmetic in AES), and Design Theory

**Extension**: Linear algebra over finite fields is the backbone of modern digital communication; it proves that the concept of linear space does not depend on the magnitude of values, but only on the completeness of algebraic structure. It is the ultimate mathematical language for processing binary and hexadecimal data.

</div>

In previous chapters, we typically assumed scalars came from the real field $\mathbb{R}$ or the complex field $\mathbb{C}$. However, in computer science and communication engineering, scalars are often taken from a finite set. **Linear Algebra over Finite Fields** (Galois Fields) deals with "discrete" linear spaces. Although intuition differs (e.g., $1+1=0$ in binary arithmetic), core concepts like bases, rank, and eigenvalues remain perfectly valid. This chapter introduces the algebraic system underlying error correction and cryptography.

---

## 52.1 Properties of Finite Fields $GF(q)$

!!! definition "Definition 52.1 (Finite Field)"
    A field containing a finite number of elements is a **Finite Field**. Its order $q$ must be a power of a prime $p$, denoted $GF(p^k)$.
    - **$GF(p)$**: The ring of integers modulo $p$, $\mathbb{Z}/p\mathbb{Z}$.
    - **Characteristic**: The smallest positive integer $p$ such that $p \cdot 1 = 0$. For $GF(p^k)$, the characteristic is $p$.

---

## 52.2 Dimensions and Counting Principles

!!! theorem "Theorem 52.1 (Size of Vector Spaces)"
    Let $V$ be an $n$-dimensional vector space over $GF(q)$. Then $V$ contains exactly **$q^n$** distinct vectors.

!!! technique "Counting: Order of $GL(n, q)$"
    The order of the group of all invertible $n \times n$ matrices over $GF(q)$ is:
    $$|GL(n, q)| = (q^n - 1)(q^n - q)(q^n - q^2) \cdots (q^n - q^{n-1})$$
    This is derived from the number of ways to choose linearly independent column vectors.

---

## 52.3 Application: Linear Codes

!!! technique "Application: Error-Correcting Codes"
    An $(n, k)$ **Linear Code** is a $k$-dimensional subspace of $GF(q)^n$.
    - **Generator Matrix $G$**: A matrix whose rows form a basis for the subspace.
    - **Parity-Check Matrix $H$**: A matrix whose rows form a basis for the dual space (null space).
    - **Decoding**: Error detection and localization are achieved via the syndrome $\mathbf{s} = H\mathbf{r}^T$ through matrix multiplication.

---

## Exercises

**1. [Basics] In $GF(2)$, calculate $(1, 0, 1) + (1, 1, 0)$.**

??? success "Solution"
    **Steps:**
    1. Addition in $GF(2)$ corresponds to the XOR operation: $1+1=0, 1+0=1, 0+0=0$.
    2. Component 1: $1+1 = 0$.
    3. Component 2: $0+1 = 1$.
    4. Component 3: $1+0 = 1$.
    **Conclusion**: The result is $(0, 1, 1)$.

**2. [Basis Counting] In a 2D space over $GF(3)$, how many non-zero vectors are there? How many total bases?**

??? success "Solution"
    **Calculation:**
    1. Total vectors: $3^2 = 9$.
    2. Non-zero vectors: $9 - 1 = 8$.
    3. **Counting Bases**:
       - Choice for the first vector: 8 (any non-zero vector).
       - Choice for the second: Must not be a multiple of the first. Multiples are $\{0, v, 2v\}$, so 3 vectors are excluded. Remaining choices: $9 - 3 = 6$.
       - Total ordered bases: $8 \times 6 = 48$.
       - Total unordered bases: $48 / 2! = 24$.

**3. [Matrix Rank] Determine the rank of $\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ over $GF(2)$.**

??? success "Solution"
    **Elimination:**
    1. $R_2 \gets R_2 - R_1$.
    2. Calculate: $1-1=0$. (Note: $-1 = 1$ in $GF(2)$).
    3. The matrix becomes $\begin{pmatrix} 1 & 1 \\ 0 & 0 \end{pmatrix}$.
    **Conclusion**: Rank is 1. This is consistent with results over $\mathbb{R}$.

**4. [Eigenvalues] Find the eigenvalues of $\begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ over $GF(2)$.**

??? success "Solution"
    **Characteristic Equation:**
    1. $p(\lambda) = \det(\lambda I - A) = \det \begin{pmatrix} \lambda & 1 \\ 1 & \lambda \end{pmatrix} = \lambda^2 - 1$.
    2. In $GF(2)$, $-1 = 1$, so $p(\lambda) = \lambda^2 + 1$.
    3. Since $(\lambda+1)^2 = \lambda^2 + 2\lambda + 1 = \lambda^2 + 1$ over $GF(2)$.
    **Conclusion**: The eigenvalue is $\lambda = 1$ with multiplicity 2.

**5. [Linear Code] Given generator $G = \begin{pmatrix} 1 & 0 & 1 \\ 0 & 1 & 1 \end{pmatrix}$ over $GF(2)$. What codewords are in the code?**

??? success "Solution"
    **Enumerate Combinations:**
    1. $0G_1 + 0G_2 = (0, 0, 0)$
    2. $1G_1 + 0G_2 = (1, 0, 1)$
    3. $0G_1 + 1G_2 = (0, 1, 1)$
    4. $1G_1 + 1G_2 = (1, 1, 0)$
    **Conclusion**: The code is $\{(0,0,0), (1,0,1), (0,1,1), (1,1,0)\}$.

**6. [Parity-Check] Find the parity-check matrix $H$ for the code in the previous problem.**

??? success "Solution"
    **Finding the Null Space:**
    1. Codeword $(x_1, x_2, x_3)$ must satisfy $x_1+x_3=0$ and $x_2+x_3=0$.
    2. This implies $x_1=x_3$ and $x_2=x_3$.
    3. Setting $x_3=1$ gives $(1, 1, 1)$.
    **Conclusion**: $H = \begin{pmatrix} 1 & 1 & 1 \end{pmatrix}$.

**7. [Inverse] Calculate the inverse of 2 in $GF(5)$.**

??? success "Solution"
    **Solve $2x \equiv 1 \pmod 5$:**
    1. $2 \times 3 = 6 \equiv 1 \pmod 5$.
    **Conclusion**: $2^{-1} = 3$. This is used frequently in field matrix inversion.

**8. [Property] Does the concept of "length" exist in $GF(p)$ vector spaces?**

??? success "Solution"
    **Analysis:**
    Traditional Euclidean length (root of squares) is not well-defined. While one can define $x^T x$, it fails the positivity axiom (e.g., in $GF(2)$, $(1,1)$ has "length squared" $1+1=0$). These are known as **Pseudo-Euclidean geometries**, focusing on orthogonality rather than magnitude.

**9. [Groups] Compute the order of $GL(2, 2)$.**

??? success "Solution"
    **Formula:**
    $|GL(2, 2)| = (2^2-1)(2^2-2) = (3)(2) = 6$.
    **Note**: These 6 matrices are the 3 permutation matrices and their linear combinations.

**10. [Application] Briefly state the use of linear algebra in the AES encryption algorithm.**

??? success "Solution"
    The **MixColumns** step in AES is a fixed linear transformation over $GF(2^8)$. It multiplies each column of the state matrix by a Maximum Distance Separable (MDS) matrix, ensuring rapid diffusion of information so that small plaintext changes lead to massive ciphertext fluctuations.

## Chapter Summary

Linear algebra over finite fields is the algebraic soul of the digital world:

1.  **Consistency of Axioms**: It proves that the logic of linear spaces remains perfectly intact even when continuity and magnitude are abandoned, showcasing the purity of mathematical truth.
2.  **Power of Counting**: Finiteness allow us to precisely quantify space sizes and group orders, providing benchmarks for search space estimates and security strength in information theory.
3.  **Bedrock of Error Correction**: By modeling communication sequences as subspaces over finite fields, linear algebra not only enables information transfer but also empowers data with self-healing capabilities, supporting the physical layer of the modern internet.
