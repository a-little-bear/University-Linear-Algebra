# Chapter 52: Linear Algebra over Finite Fields

<div class="context-flow" markdown>

**Prerequisites**: Vector Spaces (Ch04) · Matrix Algebra (Ch02) · Linear Coding Basics (Ch30)

**Chapter Outline**: Definition of Finite Fields (Galois Fields) → Characteristic and Construction → Vector Spaces over $\mathbb{F}_q$ → Subspace Counting and Gaussian Binomial Coefficients → Order of the General Linear Group $GL(n, q)$ → Algebraic Construction of Linear Codes → Generator and Parity Check Matrices → Syndrome Decoding → Applications: AES Encryption, Error-Correcting Codes (Reed-Solomon), and Network Coding

**Extension**: Linear algebra over finite fields is the intersection of discrete mathematics and continuous algebra; it is the underlying logic of the modern digital world, powering everything from 5G communications to banking encryption systems.

</div>

Linear algebra is typically taught over the fields of real or complex numbers. However, in computer science and communication theory, the most important domains are **Finite Fields** (also known as Galois Fields). In these fields, there is no "infinity" or "continuity," yet the logic of linear systems, subspaces, and bases remains perfectly intact. This chapter demonstrates how this discrete linear algebra serves as the ultimate tool for error correction and cryptography.

---

## 52.1 Construction of Finite Fields $\mathbb{F}_q$

!!! definition "Definition 52.1 (Finite Field)"
    A field containing a finite number of elements is a **finite field**.
    - The number of elements $q$ must be a power of a prime $p$, i.e., $q = p^n$.
    - **Characteristic**: The characteristic of a finite field is the smallest positive integer $p$ such that $1$ summed $p$ times equals $0$.

!!! example "Example 52.1"
    $\mathbb{F}_2$ is the simplest finite field, containing $\{0, 1\}$. Addition is equivalent to XOR, and multiplication is equivalent to logical AND.

---

## 52.2 Counting and Gaussian Binomials

!!! theorem "Theorem 52.1 (Counting Subspaces)"
    In an $n$-dimensional vector space $\mathbb{F}_q^n$, the number of $k$-dimensional subspaces is given by the **Gaussian Binomial Coefficient** (or $q$-binomial coefficient):
    $$\begin{bmatrix} n \\ k \end{bmatrix}_q = \frac{(q^n-1)(q^{n-1}-1)\cdots(q^{n-k+1}-1)}{(q^k-1)(q^{k-1}-1)\cdots(q-1)}$$
    **Contrast**: As $q \to 1$, this formula approaches the classical binomial coefficient $\binom{n}{k}$.

!!! theorem "Theorem 52.2 (Order of $GL(n, q)$)"
    The number of elements in the group of all invertible $n \times n$ matrices over $\mathbb{F}_q$ is:
    $$|GL(n, q)| = \prod_{i=0}^{n-1} (q^n - q^i)$$

---

## 52.3 Algebraic Structure of Linear Codes

!!! definition "Definition 52.2 (Linear Code)"
    An $[n, k]$ **linear code** over $\mathbb{F}_q$ is a $k$-dimensional subspace of $\mathbb{F}_q^n$.
    - **Generator Matrix $G$**: A matrix whose rows span the subspace.
    - **Parity Check Matrix $H$**: A matrix whose kernel (nullspace) is the subspace. It satisfies $GH^T = 0$.

---

## Exercises

1.  **[Basics] In $\mathbb{F}_2^3$, how many non-zero vectors are there?**
    ??? success "Solution"
        $2^3 - 1 = 7$.

2.  **[Arithmetic] Calculate $\begin{pmatrix} 1 & 1 \\ 1 & 0 \end{pmatrix} + \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix}$ over $\mathbb{F}_2$.**
    ??? success "Solution"
        $\begin{pmatrix} 0 & 1 \\ 0 & 1 \end{pmatrix}$. Note that $1+1=0$ in characteristic 2.

3.  **[GL Group] Calculate $|GL(2, 2)|$.**
    ??? success "Solution"
        $(2^2 - 1)(2^2 - 2) = 3 \times 2 = 6$. These 6 matrices form the symmetric group $S_3$.

4.  **[Counting] How many 1-dimensional subspaces are there in $\mathbb{F}_2^3$?**
    ??? success "Solution"
        $\begin{bmatrix} 3 \\ 1 \end{bmatrix}_2 = \frac{2^3-1}{2^1-1} = 7$. This corresponds to the points in the Fano plane $PG(2, 2)$.

5.  **[Dimension] If a linear code has generator matrix $G = [1 \ 1 \ 0]$, what is its dimension $k$?**
    ??? success "Solution"
        $k=1$ (the number of rows).

6.  **[Weight] What is the Hamming distance between $(1, 0, 1)$ and $(0, 1, 1)$ over $\mathbb{F}_2$?**
    ??? success "Solution"
        The number of positions where they differ is 2.

7.  **[Check Matrix] If an $[n, k]$ code has a $k \times n$ generator matrix, what are the dimensions of its parity check matrix $H$?**
    ??? success "Solution"
        $(n-k) \times n$. Its rows correspond to the annihilators of the subspace.

8.  **[Freshman's Dream] Prove $(a+b)^p = a^p + b^p$ in a field of characteristic $p$.**
    ??? success "Solution"
        Using the binomial expansion, all intermediate terms $\binom{p}{k}$ are divisible by $p$, which is 0 in the field.

9.  **[Application] Why are finite fields important in Linear Feedback Shift Registers (LFSR)?**
    ??? success "Solution"
        The state evolution of an LFSR is equivalent to powers of a companion matrix acting on a vector space over a finite field, determining the period of pseudo-random sequences.

10. **[Cryptography] Upon which finite field is the S-Box of the AES algorithm based?**
    ??? success "Solution"
        $\mathbb{F}_{2^8} = \mathbb{F}_2[x] / (x^8+x^4+x^3+x+1)$.

## Chapter Summary

Linear algebra over finite fields is the universal language of the discrete world:

1.  **Unity of Logic**: It proves that core concepts like basis, rank, and projection do not depend on "continuity" and remain robust on discrete points.
2.  **Counting Power**: Gaussian binomial coefficients perfectly merge combinatorial counting with linear subspaces, establishing mathematical criteria for describing complex discrete structures like coding spaces.
3.  **Real-world Foundation**: By transforming algebraic operations into bitwise logic, linear algebra over finite fields achieves the ultimate unity of mathematical beauty and engineering efficiency, forming the technical lifeline of modern information security and reliable transmission.
