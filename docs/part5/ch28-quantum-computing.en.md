# Chapter 28: Linear Algebra in Quantum Computing

<div class="context-flow" markdown>

**Prerequisites**: Complex Vector Spaces (Ch4) · Inner Product (Ch8) · Unitary Matrices (Ch7) · Kronecker Product (Ch19) · Probability

**Chapter Outline**: Qubits as Unit Vectors in $\mathbb{C}^2$ → Superposition and State Space → Quantum Gates as Unitary Operators → Measurement and Projections → Entanglement and Tensor Products → The Pauli Matrices → Density Matrices and Mixed States → Bloch Sphere Representation → Quantum Algorithms Overview (Shor, Grover)

**Extension**: Quantum computing is pure linear algebra performed on complex Hilbert spaces; it replaces deterministic state transitions with unitary evolutions.

</div>

Quantum computing is perhaps the most direct physical application of linear algebra. In this domain, a **qubit** is a unit vector in a 2D complex space, and **quantum gates** are unitary matrices that rotate these vectors while preserving their length. The phenomenon of **entanglement** is described by the Kronecker product of vector spaces, creating states that cannot be decomposed into individual qubits. This chapter translates the mysterious laws of quantum mechanics into the clear language of complex matrices and operators.

---

## 28.1 Qubits and Unitary Evolution

!!! definition "Definition 28.1 (Qubit State)"
    A qubit is a state vector $|\psi\rangle = \alpha |0\rangle + \beta |1\rangle \in \mathbb{C}^2$, where $|0\rangle = (1, 0)^T$ and $|1\rangle = (0, 1)^T$. Conservation of probability requires $|\alpha|^2 + |\beta|^2 = 1$.

!!! theorem "Theorem 28.1 (Unitary Gates)"
    Every reversible quantum operation is a **unitary matrix** $U$ ($U^* U = I$). Unitary matrices preserve the inner product and thus keep state vectors normalized.

---

## Exercises

1. **[Fundamentals] Is the Hadamard gate $H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ unitary?**
   ??? success "Solution"
       $H^* H = \frac{1}{2} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} = \frac{1}{2} \begin{pmatrix} 2 & 0 \\ 0 & 2 \end{pmatrix} = I$. Yes, it is unitary.

2. **[Measurement] If $|\psi\rangle = \frac{\sqrt{3}}{2} |0\rangle + \frac{1}{2} |1\rangle$, what is the probability of measuring $|1\rangle$?**
   ??? success "Solution"
       $P(1) = |\langle 1 | \psi \rangle|^2 = |1/2|^2 = 1/4 = 25\%$.

3. **[Pauli Matrices] Write the Pauli-X, Y, and Z matrices.**
   ??? success "Solution"
       $X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$. These are the generators of rotations on the Bloch sphere.

4. **[Superposition] What state results from applying $H$ to $|0\rangle$?**
   ??? success "Solution"
       $H|0\rangle = \frac{1}{\sqrt{2}} (|0\rangle + |1\rangle)$. This is an equal superposition state.

5. **[Entanglement] Show that the Bell state $|\Phi^+\rangle = \frac{1}{\sqrt{2}} (|00\rangle + |11\rangle)$ cannot be written as a product state $|a\rangle \otimes |b\rangle$.**
   ??? success "Solution"
       Let $|a\rangle = (\alpha_0, \alpha_1)$ and $|b\rangle = (\beta_0, \beta_1)$. Then $|a\rangle \otimes |b\rangle = (\alpha_0 \beta_0, \alpha_0 \beta_1, \alpha_1 \beta_0, \alpha_1 \beta_1)$. For $|00\rangle + |11\rangle$, we need $\alpha_0 \beta_1 = 0$ and $\alpha_1 \beta_0 = 0$. But this would force either $\alpha_0 \beta_0 = 0$ or $\alpha_1 \beta_1 = 0$, making it impossible to have both terms. Thus, the state is entangled.

6. **[Bloch Sphere] Relate any qubit state to a point on a sphere.**
   ??? success "Solution"
       $|\psi\rangle = \cos(\theta/2)|0\rangle + e^{i\phi}\sin(\theta/2)|1\rangle$. $(\theta, \phi)$ are spherical coordinates. Pure states lie on the surface; mixed states (density matrices) lie inside.

7. **[Density Matrix] Define the density matrix $\rho$ for a pure state $|\psi\rangle$.**
   ??? success "Solution"
       $\rho = |\psi\rangle \langle \psi|$. It is a rank-1, positive semi-definite matrix with $\operatorname{tr}(\rho) = 1$.

8. **[CNOT] Write the $4 \times 4$ matrix for the CNOT gate.**
   ??? success "Solution"
       $CNOT = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix}$. It flips the second qubit iff the first is $|1\rangle$.

9. **[Projection] Express measurement in basis $\{|v_i\rangle\}$ as a set of operators.**
   ??? success "Solution"
       The projection operators are $P_i = |v_i\rangle \langle v_i|$. The probability of outcome $i$ is $\langle \psi | P_i | \psi \rangle$.

10. **[No-Cloning] Why does the linearity of transformations forbid copying an unknown quantum state?**
    ??? success "Solution"
        A cloning operator $U$ would need to satisfy $U(|\psi\rangle|0\rangle) = |\psi\rangle|\psi\rangle$ for all $\psi$. But if $\psi = a\phi + b\xi$, linearity requires $U(\psi|0\rangle) = a\phi\phi + b\xi\xi$, which is not $(a\phi+b\xi)(a\phi+b\xi)$.

## Chapter Summary

This chapter explores the linear algebraic language of the subatomic world:

1. **State Formalism**: Defined qubits as unit vectors in $\mathbb{C}^2$, replacing binary bits with continuous superposition.
2. **Unitary Logic**: Established unitary matrices as the reversible gates of quantum computing.
3. **Tensor Synergy**: Utilized Kronecker products to describe entangled systems and multi-qubit registers.
4. **Projective Measurement**: Linked the extraction of information to the collapse of state vectors onto eigenspaces.
