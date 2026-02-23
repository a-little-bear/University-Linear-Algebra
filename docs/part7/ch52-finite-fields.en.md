# Chapter 52: Linear Algebra over Finite Fields

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch2) · Groups, Rings, and Fields · Vector Spaces (Ch4)

**Chapter Outline**: Galois Fields $\mathbb{F}_q$ → Vector Spaces over Finite Fields → Subspace Counting (q-Binomial Coefficients) → Linear Codes → Hamming Code → Reed-Solomon Codes → Cyclic Codes → Linear Feedback Shift Registers (LFSR) → Application in Cryptography

**Extension**: Linear algebra over finite fields is the foundation of digital communication (error correction) and modern cryptography (AES, Elliptic Curve Cryptography).

</div>

Linear algebra can be defined over any field, including **finite fields** (Galois fields $\mathbb{F}_q$). While many properties (like Gaussian elimination and rank) remain identical to the real case, finite fields introduce new phenomena: there is no concept of "ordering," and the number of vectors and subspaces is finite. These properties are leveraged in coding theory to detect and correct errors in data transmission.

---

## 52.1 Vector Spaces and Counting

!!! definition "Definition 52.1 (q-Binomial Coefficient)"
    The number of $k$-dimensional subspaces of an $n$-dimensional vector space over $\mathbb{F}_q$ is given by the **Gaussian binomial coefficient**:
    $$\binom{n}{k}_q = \frac{(q^n - 1)(q^n - q) \dots (q^n - q^{k-1})}{(q^k - 1)(q^k - q) \dots (q^k - q^{k-1})}$$

!!! theorem "Theorem 52.1 (Singleton Bound)"
    For any $(n, k, d)$ linear code (length $n$, dimension $k$, minimum distance $d$), $d \le n - k + 1$. Codes achieving this bound are called **MDS** (Maximum Distance Separable) codes.

---

## Exercises

1. **[Counting] How many vectors are there in the space $\mathbb{F}_2^3$?**
   ??? success "Solution"
       The size is $q^n = 2^3 = 8$. These are all 3-bit binary strings.

2. **[Subspaces] Calculate the number of 1-dimensional subspaces (lines) in $\mathbb{F}_3^2$.**
   ??? success "Solution"
       $\binom{2}{1}_3 = \frac{3^2 - 1}{3^1 - 1} = \frac{8}{2} = 4$. Geometrically, this corresponds to the 4 possible slopes in the affine plane over $\mathbb{F}_3$.

3. **[Linear Codes] Define the generator matrix $G$ and parity-check matrix $H$ of a linear code.**
   ??? success "Solution"
       $G$ is a $k 	imes n$ matrix whose rows span the code $C$. $H$ is an $(n-k) 	imes n$ matrix such that $v \in C \iff Hv^T = 0$. $H$ defines the parity-check constraints of the subspace.

4. **[Minimum Distance] Show that the minimum distance $d$ of a code is equal to the minimum number of linearly dependent columns of the parity-check matrix $H$.**
   ??? success "Solution"
       $Hv^T = 0$ is a linear combination of the columns of $H$. If a codeword $v$ has weight $w$, then $w$ columns of $H$ are dependent. The smallest such $w$ is the minimum distance $d$.

5. **[Hamming Code] Analyze the $(7, 4)$ Hamming code parity-check matrix $H$. What is its minimum distance?**
   ??? success "Solution"
       The columns of $H$ are all $2^3 - 1 = 7$ non-zero binary vectors of length 3. Since no two columns are identical (linearly dependent in $\mathbb{F}_2$) but there exist sets of three dependent columns, $d = 3$. It can correct 1 error.

6. **[Characteristic] What happens to the identity $v + v = 0$ in $\mathbb{F}_2^n$?**
   ??? success "Solution"
       In characteristic 2, $1 + 1 = 0$, so $v + v = 0$ for all $v$. This means a vector is its own additive inverse, and there is no distinction between addition and subtraction.

7. **[LFSR] Describe how a Linear Feedback Shift Register implements a linear recurrence over $\mathbb{F}_2$.**
   ??? success "Solution"
       An LFSR uses a state vector $x_k$. The next state is $x_{k+1} = A x_k$, where $A$ is a companion matrix of the feedback polynomial. If the polynomial is primitive, the LFSR generates a sequence of maximum period $2^n - 1$.

8. **[MDS Codes] Why are Reed-Solomon codes considered MDS?**
   ??? success "Solution"
       Reed-Solomon codes are based on evaluating polynomials. Since a polynomial of degree $k-1$ has at most $k-1$ roots, any $k$ evaluation points uniquely determine the polynomial. This ensures the minimum distance $d = n - k + 1$, attaining the Singleton bound.

9. **[Orthogonality] Define the dual code $C^\perp$. Can a code be its own dual?**
   ??? success "Solution"
       $C^\perp = \{u : u \cdot v = 0, \forall v \in C\}$. A code is **self-dual** if $C = C^\perp$. For example, the $(8, 4)$ extended Hamming code is self-dual.

10. **[Cryptography] How does the AES (Advanced Encryption Standard) use matrix inversion over $\mathbb{F}_{2^8}$?**
    ??? success "Solution"
        In the SubBytes step, AES performs a multiplicative inversion in the Galois field $\mathbb{F}_{2^8}$, followed by an affine transformation. The algebraic complexity of inversion over finite fields provides resistance against linear cryptanalysis.

## Chapter Summary

This chapter explores linear algebra in discrete finite settings:

1. **Subspace Combinatorics**: Used q-binomial coefficients to count configurations in finite vector spaces.
2. **Error Correction**: Formulated the theory of linear block codes using generator and parity-check matrices.
3. **Algebraic Codes**: Examined the MDS property of polynomial-based codes like Reed-Solomon.
4. **Digital Logic**: Linked matrix recurrences to the design of LFSRs and cryptographic primitives.
