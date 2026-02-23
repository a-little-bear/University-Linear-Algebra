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
       The number of vectors is $q^n = 2^3 = 8$. These are all possible 3-bit binary strings.

2. **[Subspaces] Calculate the number of 1-dimensional subspaces (lines through the origin) in $\mathbb{F}_3^2$.**
   ??? success "Solution"
       Using the q-binomial coefficient: $\binom{2}{1}_3 = \frac{3^2 - 1}{3^1 - 1} = \frac{8}{2} = 4$. Geometrically, these correspond to the 4 possible slopes in the affine plane over $\mathbb{F}_3$.

3. **[Linear Codes] Define the generator matrix $G$ and parity-check matrix $H$ of a linear code.**
   ??? success "Solution"
       $G$ is a $k \times n$ matrix whose rows form a basis for the code $C$. $H$ is an $(n-k) \times n$ matrix such that $v \in C \iff Hv^T = 0$. $H$ defines the subspace as the intersection of several hyperplanes.

4. **[Minimum Distance] Show that the minimum distance $d$ of a linear code is equal to the minimum number of linearly dependent columns of its parity-check matrix $H$.**
   ??? success "Solution"
       The equation $Hv^T = 0$ is a linear combination of the columns of $H$. If a codeword $v$ has Hamming weight $w$, then $w$ columns of $H$ must be linearly dependent. The smallest such $w$ for a non-zero $v$ is the minimum distance $d$.

5. **[Hamming Code] Analyze the $(7, 4)$ Hamming code parity-check matrix $H$. What is its minimum distance?**
   ??? success "Solution"
       The columns of $H$ are all $2^3 - 1 = 7$ non-zero binary vectors of length 3. Since no two columns are identical (they would be dependent in $\mathbb{F}_2$) but there exist sets of three columns that sum to zero, the minimum distance is $d = 3$. This code can correct any single error.

6. **[Characteristic] Explain why subtraction is identical to addition in $\mathbb{F}_2^n$.**
   ??? success "Solution"
       In a field of characteristic 2, $1 + 1 = 0 \pmod 2$, which implies $1 = -1$. Thus $v + v = 0$ for any vector $v$, making addition and subtraction the same operation.

7. **[LFSR] Describe how a Linear Feedback Shift Register implements a linear recurrence over $\mathbb{F}_2$.**
   ??? success "Solution"
       An LFSR uses a state vector $x_k$ and a transition matrix $A$ (often a companion matrix). The next state is $x_{k+1} = A x_k \pmod 2$. If the feedback polynomial is primitive, the sequence of states has the maximum possible period of $2^n - 1$.

8. **[MDS Codes] Why are Reed-Solomon codes considered MDS?**
   ??? success "Solution"
       Reed-Solomon codes are constructed by evaluating polynomials of degree $< k$. Since a non-zero polynomial of degree $k-1$ has at most $k-1$ roots, any $k$ evaluation points uniquely determine the polynomial. This property ensures $d = n - k + 1$, attaining the Singleton bound.

9. **[Orthogonality] Define the dual code $C^\perp$ and the concept of a self-dual code.**
   ??? success "Solution"
       The dual code $C^\perp = \{u : u \cdot v = 0, \forall v \in C\}$. A code is **self-dual** if $C = C^\perp$. Self-dual codes are fundamental in the study of lattices and sphere packings.

10. **[Cryptography] How is matrix inversion over finite fields utilized in the AES (Advanced Encryption Standard)?**
    ??? success "Solution"
        In the SubBytes transformation, AES computes the multiplicative inverse of each byte in the Galois field $\mathbb{F}_{2^8}$. The non-linearity of the inversion operation over a finite field provides strong resistance against linear and differential cryptanalysis.

## Chapter Summary

This chapter explores linear algebra in discrete finite settings:

1. **Subspace Combinatorics**: Used q-binomial coefficients to count configurations in finite vector spaces.
2. **Error Correction**: Formulated the theory of linear block codes using generator and parity-check matrices.
3. **Algebraic Codes**: Examined the MDS property of polynomial-based codes like Reed-Solomon.
4. **Digital Logic**: Linked matrix recurrences to the design of LFSRs and cryptographic primitives.
