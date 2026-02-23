# Chapter 28: Linear Algebra in Quantum Information

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces & Unitary Transforms (Ch07-08) · Tensor Products (Ch21) · Positive Definite Matrices (Ch16) · SVD (Ch11) · Matrix Inequalities (Ch18)

**Chapter Outline**: Quantum States (Unit Vectors in Hilbert Spaces) → Dirac Notation → Quantum Gates (Unitary Matrices) → Tensor Products & Composite Systems → Entanglement (Bell States & Schmidt Decomposition) → Measurement (Projection & POVM) → Density Matrices & Mixed States → Quantum Channels (Kraus Operators) → Quantum Error Correction → Matrix Inequalities (von Neumann Entropy)

**Extension**: Quantum information science is entirely built upon the foundation of linear algebra. States are vectors, evolution is unitary, measurement is projection, and entanglement is a non-separable tensor structure.

</div>

Quantum information science is arguably the most direct and profound application of linear algebra in modern physics. In this realm, the bits of classical computing ($0$ and $1$) are replaced by **qubits**, which are unit vectors in a complex Hilbert space. All quantum operations are unitary transformations (rotations in complex space), and the "spooky action at a distance" known as entanglement is simply a property of the tensor product of vector spaces. This chapter explores how matrix theory provides the unique language for the quantum world.

---

## 28.1 Quantum States and Hilbert Space

<div class="context-flow" markdown>

**Basic Correspondence**: Classical $\{0, 1\}$ → Qubits $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle \in \mathbb{C}^2$ → Dirac Notation → Pure vs. Mixed States.

</div>

!!! definition "Definition 28.1 (Quantum State)"
    A **quantum state** is a unit vector $|\psi\rangle$ in a Hilbert space $\mathcal{H}$, satisfying $\langle\psi|\psi\rangle = 1$. In **Dirac notation**:
    - **Ket** $|\psi\rangle$ represents a column vector.
    - **Bra** $\langle\psi|$ represents the conjugate transpose (row vector).
    - **Inner product** $\langle\phi|\psi\rangle$ is a complex scalar.
    - **Outer product** $|\psi\rangle\langle\phi|$ is a matrix (linear operator).

!!! definition "Definition 28.2 (Density Matrix)"
    A matrix $\rho$ is a valid **density matrix** representing a quantum state if:
    1.  $\rho \succeq 0$ (positive semi-definite).
    2.  $\operatorname{tr}(\rho) = 1$ (trace is 1).
    - A state is **pure** iff $\rho^2 = \rho$, meaning $\rho = |\psi\rangle\langle\psi|$. Otherwise, it is a **mixed state**.

---

## 28.2 Quantum Gates and Unitary Evolution

<div class="context-flow" markdown>

**Principle**: The evolution of a closed quantum system is described by a unitary matrix $U$ ($U^\dagger U = I$), which preserves probabilities.

</div>

!!! definition "Definition 28.3 (Standard Gates)"
    - **Pauli Gates**: $X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$, $Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}$, $Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$.
    - **Hadamard Gate**: $H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ (creates superposition).
    - **CNOT Gate**: A two-qubit gate $|0\rangle\langle 0| \otimes I + |1\rangle\langle 1| \otimes X$ that flips the second qubit if the first is $|1\rangle$.

---

## 28.3 Entanglement and Schmidt Decomposition

<div class="context-flow" markdown>

**The Essence of Quantum**: A composite state is entangled if it cannot be written as a product $|\psi\rangle_A \otimes |\psi\rangle_B$.

</div>

!!! definition "Definition 28.4 (Bell States)"
    The four **Bell states** are maximally entangled states in $\mathbb{C}^2 \otimes \mathbb{C}^2$:
    $$|\Phi^{\pm}\rangle = \frac{|00\rangle \pm |11\rangle}{\sqrt{2}}, \quad |\Psi^{\pm}\rangle = \frac{|01\rangle \pm |10\rangle}{\sqrt{2}}$$

!!! theorem "Theorem 28.1 (Schmidt Decomposition)"
    Any pure state $|\psi\rangle$ in a composite system $\mathcal{H}_A \otimes \mathcal{H}_B$ can be written as:
    $$|\psi\rangle = \sum_{i=1}^r \sqrt{\lambda_i} |u_i\rangle_A \otimes |v_i\rangle_B$$
    where $r$ is the **Schmidt rank**. The state is entangled iff $r > 1$. This is the quantum version of the Singular Value Decomposition (SVD).

---

## 28.4 Measurement and POVMs

!!! theorem "Theorem 28.2 (Born Rule)"
    For a measurement defined by a set of projectors $\{P_m\}$, the probability of outcome $m$ is:
    $$p(m) = \langle\psi|P_m|\psi\rangle = \operatorname{tr}(P_m \rho)$$
    After measurement, the state "collapses" to $\frac{P_m |\psi\rangle}{\sqrt{p(m)}}$.

---

## 28.5 Entropy and Information

!!! definition "Definition 28.5 (von Neumann Entropy)"
    The entropy of a quantum state $\rho$ is $S(\rho) = -\operatorname{tr}(\rho \log \rho)$. It quantifies the lack of information or the "mixedness" of a state. For a pure state, $S(\rho) = 0$.

---

## Exercises

1.  **[Pure/Mixed] Is $\rho = \frac{1}{2} \begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix}$ a pure state?**
    ??? success "Solution"
        $\rho^2 = \frac{1}{4} \begin{pmatrix} 2 & 2 \\ 2 & 2 \end{pmatrix} = \rho$. Yes, it is a pure state ($|+\rangle\langle+|$).

2.  **[Unitary] Prove that the Pauli $X$ gate is both Hermitian and Unitary.**
    ??? success "Solution"
        $X = X^\dagger = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$. $X^\dagger X = X^2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = I$. Thus, it is both.

3.  **[Superposition] Calculate $H|0\rangle$.**
    ??? success "Solution"
        $\frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 \\ 0 \end{pmatrix} = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 \\ 1 \end{pmatrix} = \frac{|0\rangle + |1\rangle}{\sqrt{2}} = |+\rangle$.

4.  **[Entanglement] Use Schmidt rank to show that $|00\rangle$ is not entangled.**
    ??? success "Solution"
        The coefficient matrix is $C = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$. Its rank is 1. Thus, the Schmidt rank is 1, so it is a product state (not entangled).

5.  **[Bell] Show that for $|\Phi^+\rangle$, the local state of qubit A is maximally mixed.**
    ??? success "Solution"
        The partial trace $\rho_A = \operatorname{tr}_B(|\Phi^+\rangle\langle\Phi^+|) = \frac{1}{2} (|0\rangle\langle 0| + |1\rangle\langle 1|) = \frac{1}{2} I$. This is the maximally mixed state.

6.  **[Measurement] If we measure $|\psi\rangle = \frac{\sqrt{3}}{2} |0\rangle + \frac{1}{2} |1\rangle$ in the computational basis, what is the probability of getting $|1\rangle$?**
    ??? success "Solution"
        $p(1) = |\langle 1 | \psi \rangle|^2 = |1/2|^2 = 1/4$.

7.  **[Gate] Find the matrix for $H \otimes H$.**
    ??? success "Solution"
        $\frac{1}{2} \begin{pmatrix} 1 & 1 & 1 & 1 \\ 1 & -1 & 1 & -1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & -1 & 1 \end{pmatrix}$.

8.  **[Entropy] Calculate $S(\rho)$ for $\rho = \frac{1}{2} I$ (single qubit).**
    ??? success "Solution"
        Eigenvalues are $1/2, 1/2$. $S(\rho) = -(1/2 \log 1/2 + 1/2 \log 1/2) = \log 2 = 1$ bit.

9.  **[Trace] Prove $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \operatorname{tr}(B)$.**
    ??? success "Solution"
        $\sum_{i,j} (A \otimes B)_{ii,jj} = \sum_i A_{ii} \sum_j B_{jj} = \operatorname{tr}(A)\operatorname{tr}(B)$.

****

??? success "Solution"
    

## Chapter Summary

Linear algebra is the "operating system" of the quantum universe:

1.  **State Logic**: Replaced discrete bits with complex vectors, using Hilbert spaces to accommodate the principle of superposition.
2.  **Symmetry & Evolution**: Identified unitary matrices as the only valid operators for physical change, ensuring that probability is always conserved.
3.  **Compositional Complexity**: Used tensor products to explain how the state space of a system grows exponentially, providing the mathematical basis for quantum advantage.
4.  **Information Geometry**: Leveraged density matrices and entropy to bridge the gap between microscopic pure states and macroscopic statistical observations.
