# Chapter 28: Applications of Linear Algebra in Quantum Information

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Spaces (Ch8) · Kronecker Product (Ch19) · Normal Operators (Ch8) · SVD (Ch11)

**Chapter Outline**: Quantum States and Hilbert Space → Operators and Observables (Hermitian Operators) → Spectral Decomposition → Measurement and Projection Operators → Superposition and Tensor Products → Quantum Entanglement and Schmidt Decomposition → Density Matrices and Trace-Class Operators → Applications (Quantum Gates, Quantum Search Algorithms, Bell Inequalities)

**Extension**: Quantum mechanics is linear algebra over complex inner product spaces; every axiom of quantum info is essentially a physical description of linear operator properties.

</div>

Quantum information science extends classical logic bits (0 or 1) into quantum bits (unit vectors in state space). In this microscopic world, all physical processes—from state evolution to info retrieval—strictly follow the theory of linear operators on complex linear spaces.

---

## 28.1 Core Concepts and Axioms

!!! definition "Definition 28.1 (Quantum State and Operator)"
    - The state of an isolated quantum system is represented by a unit vector $|\psiangle$ in a Hilbert space.
    - Physical observables are represented by **Hermitian operators** $H$ (satisfying $H = H^*$).

!!! theorem "Theorem 28.3 (Schmidt Decomposition)"
    Any state $|\Psiangle$ of a composite system $AB$ can be decomposed as $|\Psiangle = \sum \sqrt{\lambda_i} |u_iangle |v_iangle$, where $\lambda_i$ are eigenvalues of the reduced density matrix.

---

## Exercises

1. **[Fundamentals] Prove the condition for a superposition $|\psiangle = a|0angle + b|1angle$ to be a unit vector.**
   ??? success "Solution"
       $\langle \psi | \psi angle = 1$.
       Expanding gives $(a^* \langle 0| + b^* \langle 1|)(a|0angle + b|1angle) = |a|^2 \langle 0|0angle + |b|^2 \langle 1|1angle = |a|^2 + |b|^2 = 1$. This is the probability normalization condition.

2. **[Hermitian] If $H$ represents energy (Hamiltonian), prove its eigenvalues are real.**
   ??? success "Solution"
       This is a fundamental property of Hermitian operators. Let $H|vangle = E|vangle$, then $\langle v|H|vangle = E \langle v|vangle$.
       Taking the adjoint gives $\langle v|H^*|vangle = E^* \langle v|vangle$.
       Since $H^* = H$, we have $E = E^*$, so $E$ is real. This ensures physical measurements (energy) are observable real numbers.

3. **[Quantum Gates] Verify the Hadamard gate $H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix}$ is a unitary matrix.**
   ??? success "Solution"
       Calculate $HH^*$: $\frac{1}{2} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} \begin{pmatrix} 1 & 1 \ 1 & -1 \end{pmatrix} = \frac{1}{2} \begin{pmatrix} 2 & 0 \ 0 & 2 \end{pmatrix} = I$.
       Unitarity ensures quantum evolution preserves the length of probability amplitudes, i.e., no info is lost.

4. **[Tensor Product] Compute the composite state $|+angle \otimes |-angle$, where $|+angle = \frac{1}{\sqrt{2}}(|0angle+|1angle)$ and $|-angle = \frac{1}{\sqrt{2}}(|0angle-|1angle)$.**
   ??? success "Solution"
       $|+angle \otimes |-angle = \frac{1}{2} \begin{pmatrix} 1 \ 1 \end{pmatrix} \otimes \begin{pmatrix} 1 \ -1 \end{pmatrix} = \frac{1}{2} [1, -1, 1, -1]^T$.
       Expanded as $|00angle - |01angle + |10angle - |11angle$.

5. **[Entanglement] Determine if the Bell state $|\Phi^+angle = \frac{1}{\sqrt{2}}(|00angle + |11angle)$ is separable.**
   ??? success "Solution"
       Not separable. If it were, then $|\Phi^+angle = (a|0angle+b|1angle) \otimes (c|0angle+d|1angle) = ac|00angle + ad|01angle + bc|10angle + bd|11angle$.
       Equating coefficients gives $ad=0$ and $bc=0$. This would imply $ac=0$ or $bd=0$, leading to a contradiction. Hence it is an entangled state.

6. **[Measurement] Given state $|\psiangle = \frac{\sqrt{3}}{2}|0angle + \frac{1}{2}|1angle$ and measurement basis $\{|0angle, |1angle\}$, find the probability of obtaining 0.**
   ??? success "Solution"
       $P(0) = |\langle 0|\psiangle|^2 = | \frac{\sqrt{3}}{2} \langle 0|0angle + \frac{1}{2} \langle 0|1angle |^2 = (\frac{\sqrt{3}}{2})^2 = 3/4$.

7. **[Density Matrix] What idempotency does the density matrix $ho = |\psiangle\langle \psi|$ of a pure state satisfy?**
   ??? success "Solution"
       $ho^2 = (|\psiangle\langle \psi|)(|\psiangle\langle \psi|) = |\psiangle (\langle \psi|\psiangle) \langle \psi| = |\psiangle \langle \psi| = ho$.
       The density matrix of a pure state is an idempotent projection operator with $\operatorname{tr}(ho) = 1$.

8. **[Partial Trace] Why does the Partial Trace operation correspond to "ignoring a subsystem"?**
   ??? success "Solution"
       The partial trace $ho_A = \operatorname{tr}_B(ho_{AB})$ is a linear map from the composite density matrix to the subsystem's. Algebraically, it eliminates all info from B by summing over its basis, leaving only the marginal probability distribution for A.

9. **[Pauli Matrices] Calculate the commutator $[X, Z]$ of Pauli-X and Pauli-Z.**
   ??? success "Solution"
       $X = \begin{pmatrix} 0 & 1 \ 1 & 0 \end{pmatrix}, Z = \begin{pmatrix} 1 & 0 \ 0 & -1 \end{pmatrix}$.
       $XZ = \begin{pmatrix} 0 & -1 \ 1 & 0 \end{pmatrix}, ZX = \begin{pmatrix} 0 & 1 \ -1 & 0 \end{pmatrix}$.
       $[X, Z] = XZ - ZX = \begin{pmatrix} 0 & -2 \ 2 & 0 \end{pmatrix} = -2iY$. Non-zero commutators reflect quantum uncertainty.

10. **[Schmidt Rank] What does the number of non-zero coefficients in a Schmidt decomposition represent?**
    ??? success "Solution"
        The number of non-zero coefficients is the **Schmidt Rank**. If Schmidt Rank $> 1$, the composite system is in an entangled state. It is the most basic algebraic metric for quantifying entanglement.

## Chapter Summary

Quantum information is the leading-edge application of linear algebra:

1. **Superposition and Tensors**: Qubits occupy space via linear superposition; composite systems expand explosively via tensor products.
2. **Operators as Logic**: All gate operations in quantum computing are unitary rotations in complex space.
3. **Eigenvalues as Results**: All physical measurements are essentially extracting eigenvalues and their probabilities from the eigenspaces defined by projection operators.
