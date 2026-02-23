# Chapter 28: Linear Algebra in Quantum Information

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch08) · Eigenvalues (Ch06) · Unitary Matrices (Ch55A) · Tensor Products (Ch21)

**Chapter Outline**: Algebraic Axioms of the Quantum World → States as Unit-norm Vectors (Bra-ket Notation) → The Superposition Principle → Evolution: Unitary Transformations ($U$) → Measurement: Eigenvalues and Projections of Hermitian Operators ($M$) → Composite Systems: Tensor Products and Entanglement → Pure vs. Mixed States: The Density Matrix → Matrix Representation of Quantum Gates (Hadamard, CNOT) → Applications: Qubits, Quantum Teleportation, and Shor’s Algorithm

**Extension**: Quantum information is the most avant-garde and pure application field of linear algebra; it equates physical reality entirely with operator operations in complex Hilbert spaces. It proves that the essence of the microscopic world is high-dimensional linear superposition—the only mathematical ticket to understanding the next computing revolution.

</div>

In classical computers, information is stored as 0 or 1. In quantum computers, information exists as complex vectors. **Quantum Information** is built entirely upon the axioms of complex linear algebra. By utilizing **Superposition** (linear combination) and **Entanglement** (irreducibility of tensor products), quantum systems demonstrate computational power beyond classical logic. This chapter introduces how to use the language of linear algebra to describe this microscopic world of probability and waves.

---

## 28.1 Quantum States and Bra-ket Notation

!!! definition "Definition 28.1 (Qubit)"
    A qubit is a unit vector in $\mathbb{C}^2$:
    $$|\psi\rangle = \alpha |0\rangle + \beta |1\rangle, \quad |\alpha|^2 + |\beta|^2 = 1$$
    - $|0\rangle = (1, 0)^T$ and $|1\rangle = (0, 1)^T$ form the computational basis.
    - $\langle \phi | \psi \rangle$ denotes the complex inner product.

---

## 28.2 Evolution and Measurement

!!! theorem "Theorem 28.1 (Quantum Evolution)"
    The evolution of a closed quantum system is described by a **Unitary Matrix** $U$: $|\psi'\rangle = U |\psi\rangle$.
    **Physical Meaning**: Unitarity ensures that the sum of probabilities remains exactly 1.

!!! technique "Axiom: Measurement"
    Measurement of a physical observable $M$ (a Hermitian matrix) yields one of its eigenvalues $\lambda_i$. Upon measurement, the state collapses into the corresponding eigenspace.

---

## 28.3 Entanglement and Tensor Products

!!! theorem "Theorem 28.2 (Entanglement Criterion)"
    For a composite system $H_A \otimes H_B$, a state $|\Psi\rangle$ is **Entangled** if it cannot be decomposed as a tensor product of individual states $|\psi_A\rangle \otimes |\psi_B\rangle$.
    **Algebraic Essence**: Entangled states correspond to tensors with rank greater than 1.

---

## Exercises

**1. [Basics] Determine if the qubit $|\psi\rangle = \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$ is normalized.**

??? success "Solution"
    **Verification:**
    Compute the squared norm: $\| \psi \|^2 = |1/\sqrt{2}|^2 + |1/\sqrt{2}|^2 = 1/2 + 1/2 = 1$.
    **Conclusion**: Yes. This state represents a uniform superposition of 0 and 1.

**2. [Calculation] Find the result of the Hadamard gate $H = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ acting on $|0\rangle$.**

??? success "Solution"
    **Matrix Multiplication:**
    $H \begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \end{pmatrix}$.
    **Conclusion**: The Hadamard gate transforms the basis state $|0\rangle$ into the superposition state $\frac{|0\rangle + |1\rangle}{\sqrt{2}}$.

**3. [Measurement] For the state $|\psi\rangle = \frac{\sqrt{3}}{2} |0\rangle + \frac{1}{2} |1\rangle$, what is the probability of measuring 1 in the computational basis?**

??? success "Solution"
    **Born Rule:**
    The probability $P(1) = |\langle 1 | \psi \rangle|^2$.
    The inner product is $1 \cdot (1/2) = 1/2$.
    **Conclusion**: $P(1) = (1/2)^2 = 1/4$ or 25%.

**4. [Entanglement] Is the Bell state $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$ entangled?**

??? success "Solution"
    **Analysis:**
    1. Write as a coefficient matrix: $M = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$.
    2. Compute the rank: $\operatorname{rank}(M) = 2$.
    3. Separable (non-entangled) states must have coefficient matrices of rank 1.
    **Conclusion**: Since rank > 1, this is a **maximally entangled state**.

**5. [Density Matrix] Calculate the density matrix $\rho$ for a pure state $|\psi\rangle$.**

??? success "Solution"
    **Formula:**
    $\rho = |\psi\rangle \langle \psi|$ (the outer product).
    **Properties**: Pure state density matrices satisfy $\rho^2 = \rho$ and $\operatorname{tr}(\rho) = 1$.

**6. [Unitarity] Verify the Pauli-X gate $X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$ is unitary.**

??? success "Solution"
    **Calculation:**
    $X^* X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix} = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$.
    **Conclusion**: It is unitary. In quantum logic, it acts as the NOT gate.

**7. [Hermiticity] Why must measurement operators be Hermitian?**

??? success "Solution"
    **Physical Reason:**
    Measurement outcomes (observables) must be **real numbers**. In linear algebra, the class of operators that guarantees all eigenvalues are real is the class of **Hermitian operators**.

**8. [Composition] Express the composite basis $|0\rangle \otimes |1\rangle$ as a vector.**

??? success "Solution"
    **Using Kronecker Product:**
    $\begin{pmatrix} 1 \\ 0 \end{pmatrix} \otimes \begin{pmatrix} 0 \\ 1 \end{pmatrix} = \begin{pmatrix} 0 \\ 1 \\ 0 \\ 0 \end{pmatrix}$.
    This is the second standard basis vector in 4D space, often written as $|01\rangle$.

**9. [Uncertainty] Briefly explain why $[A, B] \neq 0$ leads to measurement incompatibility.**

??? success "Solution"
    **Algebraic Reason:**
    If operators do not commute, they do not share a common basis of eigenvectors (Ch63A). Measuring $A$ collapses the state to an eigenspace of $A$ that is typically not an eigenspace of $B$. Therefore, a subsequent measurement of $B$ introduces new randomness, making it impossible for both quantities to have sharp values simultaneously.

**10. [Application] What is the "Bloch Sphere" representation of a qubit?**

??? success "Solution"
    A qubit is a 2D complex vector of norm 1, and overall phase is physically irrelevant. Its state can be mapped uniquely to a point on a 3D unit sphere: $|\psi\rangle = \cos(\theta/2)|0\rangle + e^{i\phi}\sin(\theta/2)|1\rangle$. This proves that the algebraic structure of a single qubit is isomorphic to geometric rotations of a 3D sphere.

## Chapter Summary

Linear algebra is the only syntax for describing quantum reality:

1.  **Algebraic Essence of Superposition**: Parallelism in the quantum world stems from linear combinations in vector spaces; through basis changes, we realize flexible transitions between different physical attributes.
2.  **Tensor Nature of Entanglement**: Entangled states reveal powerful correlations between irreducible elements in tensor spaces, proving that the information of a whole can exceed the simple sum of its parts.
3.  **Conservation via Unitarity**: Unitary matrix evolution establishes the logical reversibility and probability conservation of quantum computing, serving as the highest governing principle for designing quantum algorithms and error-correction protocols.
