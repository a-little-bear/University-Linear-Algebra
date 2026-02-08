# Chapter 28  Linear Algebra in Quantum Information

<div class="context-flow" markdown>

**Prerequisites**: Inner product spaces and unitary transformations (Ch7-8) · tensor products (Ch21) · positive definite matrices (Ch16) · SVD (Ch11) · matrix inequalities (Ch18) · Kronecker product (Ch19) · **Arc**: Quantum states (unit vectors in Hilbert spaces) → quantum gates (unitary matrices) → tensor products (multi-body systems) → entanglement (Bell states / Schmidt decomposition) → measurement (projective / POVM) → density matrices and quantum channels (Kraus operators) → quantum error correction (stabilizer codes) → matrix inequalities (von Neumann entropy / strong subadditivity)
**Essence**: The mathematical language of quantum mechanics is linear algebra — states are vectors, evolution is unitary matrices, observables are self-adjoint operators, entanglement is nontrivial structure in tensor products, and information measures are inequalities on matrix functions

</div>

Quantum information science is built entirely on the foundation of linear algebra. Quantum states are vectors in Hilbert spaces, quantum evolution is unitary transformation, quantum measurement is projection, quantum entanglement is non-separable structure in tensor product spaces, and quantum channels are completely positive trace-preserving maps. This chapter proceeds from quantum states and Hilbert spaces, systematically demonstrating the central role of linear algebra tools in quantum computing, quantum communication, and quantum error correction, concluding with matrix inequalities in quantum information.

---

## 28.1 Quantum States and Hilbert Spaces

<div class="context-flow" markdown>

**Basic correspondence**: Classical bit $\{0, 1\}$ → qubit $|\psi\rangle = \alpha|0\rangle + \beta|1\rangle \in \mathbb{C}^2$ → Dirac notation → pure vs mixed states
**Link**: Direct application of inner product spaces from Ch8

</div>

The qubit is the fundamental unit of quantum information, mathematically described as a unit vector in a two-dimensional complex Hilbert space. Dirac notation provides a concise and powerful linear algebra language for quantum mechanics.

!!! definition "Definition 28.1 (Quantum state)"
    A **quantum state** is a unit vector $|\psi\rangle$ in a Hilbert space $\mathcal{H}$, satisfying $\langle\psi|\psi\rangle = 1$. In Dirac notation:

    - A **ket** $|\psi\rangle$ represents a column vector.
    - A **bra** $\langle\psi|$ represents the conjugate transpose of $|\psi\rangle$, i.e., a row vector.
    - The **inner product** $\langle\phi|\psi\rangle$ is a complex number.
    - The **outer product** $|\psi\rangle\langle\phi|$ is a matrix (linear operator).

    A global phase $e^{i\theta}|\psi\rangle$ represents the same physical state as $|\psi\rangle$, so the physical state space is actually the projective Hilbert space.

!!! definition "Definition 28.2 (Qubit)"
    The state of a **qubit** is a unit vector in $\mathbb{C}^2$:

    $$
    |\psi\rangle = \alpha|0\rangle + \beta|1\rangle, \quad |\alpha|^2 + |\beta|^2 = 1,
    $$

    where $|0\rangle = \begin{pmatrix} 1 \\ 0 \end{pmatrix}$, $|1\rangle = \begin{pmatrix} 0 \\ 1 \end{pmatrix}$ form the **computational basis**, and $\alpha, \beta \in \mathbb{C}$. The quantities $|\alpha|^2$ and $|\beta|^2$ are the probabilities of measuring $|0\rangle$ and $|1\rangle$, respectively.

    The joint state space of $n$ qubits is the tensor product $\mathcal{H} = (\mathbb{C}^2)^{\otimes n} \cong \mathbb{C}^{2^n}$, whose dimension grows exponentially with the number of qubits.

!!! theorem "Theorem 28.1 (Density matrix characterization)"
    An operator $\rho$ is a valid density matrix (representing a quantum state) if and only if it satisfies two conditions:

    1. $\rho \succeq 0$ (positive semidefinite).
    2. $\operatorname{tr}(\rho) = 1$.

    Furthermore:

    - $\rho$ represents a **pure state** $\iff$ $\rho^2 = \rho$ (idempotent) $\iff$ $\operatorname{rank}(\rho) = 1$ $\iff$ $\operatorname{tr}(\rho^2) = 1$.
    - $\rho$ represents a **mixed state** $\iff$ $\operatorname{tr}(\rho^2) < 1$, in which case $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ ($p_i > 0$, $\sum p_i = 1$).

??? proof "Proof"
    **Necessity**: A pure state $\rho = |\psi\rangle\langle\psi|$ is clearly positive semidefinite (for any $|\phi\rangle$, $\langle\phi|\rho|\phi\rangle = |\langle\phi|\psi\rangle|^2 \ge 0$), and $\operatorname{tr}(\rho) = \langle\psi|\psi\rangle = 1$. A mixed state $\rho = \sum_i p_i |\psi_i\rangle\langle\psi_i|$ is a convex combination of positive semidefinite matrices, hence positive semidefinite, and $\operatorname{tr}(\rho) = \sum_i p_i = 1$.

    **Sufficiency**: If $\rho \succeq 0$ and $\operatorname{tr}(\rho) = 1$, then $\rho$ has a spectral decomposition $\rho = \sum_i \lambda_i |e_i\rangle\langle e_i|$, where $\lambda_i \ge 0$ (positive semidefinite) and $\sum_i \lambda_i = 1$ (trace one). Thus the $\lambda_i$ form a probability distribution, and $\rho$ is indeed a probabilistic mixture of pure states.

    **Pure state criterion**: If $\rho = |\psi\rangle\langle\psi|$, then $\rho^2 = |\psi\rangle\langle\psi|\psi\rangle\langle\psi| = |\psi\rangle\langle\psi| = \rho$. Conversely, suppose $\operatorname{tr}(\rho^2) = 1$. The eigenvalues $\lambda_i \ge 0$ satisfy $\sum \lambda_i = 1$ and $\sum \lambda_i^2 = 1$. By the Cauchy-Schwarz inequality $\sum \lambda_i^2 \le (\max \lambda_i)\sum \lambda_i = \max \lambda_i$, so $\max \lambda_i \ge 1$. Since $\lambda_i \le \sum \lambda_i = 1$, exactly one eigenvalue equals $1$ and the rest are $0$. $\blacksquare$

!!! example "Example 28.1"
    **The Bloch sphere representation.** A single-qubit pure state (up to global phase) can be parameterized as

    $$
    |\psi\rangle = \cos\frac{\theta}{2}|0\rangle + e^{i\phi}\sin\frac{\theta}{2}|1\rangle, \quad \theta \in [0, \pi],\, \phi \in [0, 2\pi).
    $$

    The corresponding density matrix is $\rho = |\psi\rangle\langle\psi| = \frac{1}{2}(I + \mathbf{r} \cdot \boldsymbol{\sigma})$, where the Bloch vector $\mathbf{r} = (\sin\theta\cos\phi, \sin\theta\sin\phi, \cos\theta)$, $|\mathbf{r}| = 1$.

    Mixed states correspond to $|\mathbf{r}| < 1$ (interior of the sphere). The maximally mixed state $\rho = I/2$ corresponds to $\mathbf{r} = \mathbf{0}$ (center).

    | State | Bloch vector | Position |
    |:---:|:---:|:---:|
    | $\|0\rangle$ | $(0, 0, 1)$ | North pole |
    | $\|1\rangle$ | $(0, 0, -1)$ | South pole |
    | $\|+\rangle = \frac{1}{\sqrt{2}}(\|0\rangle + \|1\rangle)$ | $(1, 0, 0)$ | Positive $x$-axis |
    | $\|-\rangle = \frac{1}{\sqrt{2}}(\|0\rangle - \|1\rangle)$ | $(-1, 0, 0)$ | Negative $x$-axis |

---

## 28.2 Quantum Gates and Unitary Transformations

<div class="context-flow" markdown>

**Core principle**: Quantum evolution = unitary transformation $U^\dagger U = I$ → preserves inner products (probability conservation) → single-qubit gates (Pauli, Hadamard, S, T) → multi-qubit gates (CNOT, Toffoli) → universal gate sets
**Link**: Direct application of unitary matrix theory from Ch7

</div>

The evolution of a closed quantum system is described by unitary matrices. Quantum gates are the fundamental building blocks of quantum computation.

!!! definition "Definition 28.3 (Quantum gate)"
    A **quantum gate** is a unitary matrix $U$ ($U^\dagger U = UU^\dagger = I$) acting on qubits.

    **Single-qubit Pauli gates**:

    $$
    X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix},\  Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix},\  Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.
    $$

    **Hadamard gate**: $H = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$, satisfying $H|0\rangle = |+\rangle$, $H|1\rangle = |-\rangle$.

    **Phase gates**: $S = \begin{pmatrix} 1 & 0 \\ 0 & i \end{pmatrix}$, $T = \begin{pmatrix} 1 & 0 \\ 0 & e^{i\pi/4} \end{pmatrix}$. Note that $S = T^2$ and $Z = S^2$.

!!! definition "Definition 28.4 (Multi-qubit gates)"
    The **CNOT (controlled-NOT) gate** is the most important two-qubit gate:

    $$
    \text{CNOT} = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix} = |0\rangle\langle 0| \otimes I + |1\rangle\langle 1| \otimes X.
    $$

    CNOT flips the second qubit (target) conditioned on the first qubit (control) being $|1\rangle$.

    The **Toffoli gate** (CCNOT) is a three-qubit gate: it flips the target qubit if and only if both control qubits are $|1\rangle$. The Toffoli gate is universal for classical reversible computation.

!!! theorem "Theorem 28.2 (Universality of gate sets)"
    The gate set $\{H, T, \text{CNOT}\}$ is **universal**: for any $n$-qubit unitary matrix $U \in U(2^n)$ and any $\varepsilon > 0$, there exists a finite quantum circuit $\tilde{U}$ composed of $H$, $T$, and $\text{CNOT}$ gates such that

    $$
    \|U - \tilde{U}\| < \varepsilon.
    $$

    More precisely, the required gate count is $O(4^n \cdot \operatorname{poly}(\log(1/\varepsilon)))$.

??? proof "Proof"
    The proof proceeds in three steps.

    **Step 1**: Any $n$-qubit unitary can be decomposed as a product of two-qubit gates. Via the quantum analogue of QR decomposition (Givens rotations), any unitary in $U(2^n)$ can be written as a product of $O(4^n)$ two-qubit unitaries.

    **Step 2**: Any two-qubit gate can be implemented using CNOT and single-qubit gates. This is achieved through the KAK decomposition (Cartan decomposition): any element of $SU(4)$ can be written as $(A_1 \otimes A_2) \cdot \exp(i(a X \otimes X + b Y \otimes Y + c Z \otimes Z)) \cdot (A_3 \otimes A_4)$, where $A_i \in SU(2)$. The entangling part in the middle can be realized with at most $3$ CNOT gates.

    **Step 3**: $\{H, T\}$ is dense in $SU(2)$ (Solovay-Kitaev theorem). The group generated by $H$ and $T$ is dense in $SU(2)$, and any $SU(2)$ element can be approximated to precision $\varepsilon$ using $O(\log^c(1/\varepsilon))$ gates from $\{H, T\}$, where $c \approx 3.97$.

    Combining the three steps, the total gate count is $O(4^n \cdot \operatorname{poly}(\log(1/\varepsilon)))$. $\blacksquare$

!!! example "Example 28.2"
    **Bell state preparation circuit.** Using the Hadamard and CNOT gates to prepare a Bell state from $|00\rangle$:

    $$
    |00\rangle \xrightarrow{H \otimes I} \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)|0\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |10\rangle) \xrightarrow{\text{CNOT}} \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle) = |\Phi^+\rangle.
    $$

    In matrix form:

    $$
    \text{CNOT} \cdot (H \otimes I) = \begin{pmatrix}1&0&0&0\\0&1&0&0\\0&0&0&1\\0&0&1&0\end{pmatrix} \cdot \frac{1}{\sqrt{2}}\begin{pmatrix}1&0&1&0\\0&1&0&1\\1&0&-1&0\\0&1&0&-1\end{pmatrix} = \frac{1}{\sqrt{2}}\begin{pmatrix}1&0&1&0\\0&1&0&1\\0&1&0&-1\\1&0&-1&0\end{pmatrix}.
    $$

    Acting on $|00\rangle = (1,0,0,0)^T$ yields $\frac{1}{\sqrt{2}}(1,0,0,1)^T = |\Phi^+\rangle$.

---

## 28.3 Tensor Products and Multi-Qubit Systems

<div class="context-flow" markdown>

**Core**: $\mathcal{H}_A \otimes \mathcal{H}_B$ = state space of the composite system → Kronecker product gives the matrix representation → exponential dimension growth = source of quantum advantage
**Link**: Direct application of Kronecker products from Ch19 and tensor algebra from Ch21

</div>

The mathematical description of multi-body quantum systems relies on tensor products. The state space of a composite system is the tensor product of the state spaces of its subsystems, leading to exponential dimensional growth.

!!! definition "Definition 28.5 (Composite quantum system)"
    The state space of a **composite quantum system** consisting of subsystem $A$ (state space $\mathcal{H}_A$, $\dim = d_A$) and $B$ (state space $\mathcal{H}_B$, $\dim = d_B$) is

    $$
    \mathcal{H}_{AB} = \mathcal{H}_A \otimes \mathcal{H}_B, \quad \dim \mathcal{H}_{AB} = d_A \cdot d_B.
    $$

    If $\{|i\rangle_A\}$ and $\{|j\rangle_B\}$ are orthonormal bases for $\mathcal{H}_A$ and $\mathcal{H}_B$ respectively, then $\{|i\rangle_A \otimes |j\rangle_B\}$ is an orthonormal basis for $\mathcal{H}_{AB}$.

    For operators, the joint action of an operator $O_A$ on $A$ and $O_B$ on $B$ is given by the Kronecker product $O_A \otimes O_B$.

!!! theorem "Theorem 28.3 (Properties of tensor products)"
    Let $A \in \mathbb{C}^{m \times m}$, $B \in \mathbb{C}^{n \times n}$. The Kronecker product $A \otimes B \in \mathbb{C}^{mn \times mn}$ satisfies:

    1. **Mixed product property**: $(A \otimes B)(C \otimes D) = AC \otimes BD$.
    2. **Spectral property**: If $A$ has eigenvalues $\{\alpha_i\}$ and $B$ has eigenvalues $\{\beta_j\}$, then $A \otimes B$ has eigenvalues $\{\alpha_i \beta_j\}$.
    3. **Trace**: $\operatorname{tr}(A \otimes B) = \operatorname{tr}(A) \cdot \operatorname{tr}(B)$.
    4. **Determinant**: $\det(A \otimes B) = (\det A)^n (\det B)^m$.
    5. **Unitarity preservation**: If $A$ and $B$ are both unitary, then $A \otimes B$ is also unitary.

??? proof "Proof"
    (1): The $(i,j)$-block of $(A \otimes B)(C \otimes D)$ is $\sum_k A_{ik}C_{kj} \cdot BD = (AC)_{ij} BD$, which is the $(i,j)$-block of $AC \otimes BD$.

    (2): If $A\mathbf{u} = \alpha\mathbf{u}$ and $B\mathbf{v} = \beta\mathbf{v}$, then $(A \otimes B)(\mathbf{u} \otimes \mathbf{v}) = A\mathbf{u} \otimes B\mathbf{v} = \alpha\beta(\mathbf{u} \otimes \mathbf{v})$.

    (3): $\operatorname{tr}(A \otimes B) = \sum_{i,j}(A \otimes B)_{(i,j),(i,j)} = \sum_i A_{ii} \sum_j B_{jj} = \operatorname{tr}(A)\operatorname{tr}(B)$.

    (5): $(A \otimes B)^\dagger(A \otimes B) = (A^\dagger \otimes B^\dagger)(A \otimes B) = A^\dagger A \otimes B^\dagger B = I \otimes I = I$. $\blacksquare$

!!! example "Example 28.3"
    **Matrix representation of a two-qubit system.** Consider the two-qubit state $|\psi\rangle = \frac{1}{\sqrt{2}}|00\rangle + \frac{1}{2}|01\rangle + \frac{1}{2}|10\rangle$.

    Vector representation: $|\psi\rangle = \frac{1}{\sqrt{2}}(1, 0)^T \otimes (1, 0)^T + \frac{1}{2}(1, 0)^T \otimes (0, 1)^T + \frac{1}{2}(0, 1)^T \otimes (1, 0)^T = (\frac{1}{\sqrt{2}}, \frac{1}{2}, \frac{1}{2}, 0)^T$.

    Normalization check: $\frac{1}{2} + \frac{1}{4} + \frac{1}{4} + 0 = 1$.

    Applying $X \otimes Z$ to this state:

    $$
    (X \otimes Z)|\psi\rangle = \begin{pmatrix}0&0&1&0\\0&0&0&-1\\1&0&0&0\\0&-1&0&0\end{pmatrix}\begin{pmatrix}1/\sqrt{2}\\1/2\\1/2\\0\end{pmatrix} = \begin{pmatrix}1/2\\0\\1/\sqrt{2}\\-1/2\end{pmatrix}.
    $$

---

## 28.4 Entanglement and Bell States

<div class="context-flow" markdown>

**Core concept**: Product state $|\alpha\rangle \otimes |\beta\rangle$ vs entangled state (non-separable) → Bell states (maximally entangled) → Schmidt decomposition = quantum version of SVD → Schmidt rank determines entanglement
**Link**: Core application of SVD from Ch11 and tensor products from Ch21

</div>

Quantum entanglement is the most essential distinction between quantum and classical information, with its mathematical essence being non-separable structure in tensor product spaces.

!!! definition "Definition 28.6 (Entangled and separable states)"
    A bipartite pure state $|\psi\rangle_{AB} \in \mathcal{H}_A \otimes \mathcal{H}_B$ is called **separable** if there exist $|\alpha\rangle \in \mathcal{H}_A$ and $|\beta\rangle \in \mathcal{H}_B$ such that

    $$
    |\psi\rangle_{AB} = |\alpha\rangle \otimes |\beta\rangle.
    $$

    Otherwise it is called an **entangled state**.

    For mixed states, $\rho_{AB}$ is separable if and only if there exists a decomposition $\rho_{AB} = \sum_i p_i \rho_A^{(i)} \otimes \rho_B^{(i)}$ ($p_i > 0$, $\sum p_i = 1$).

!!! definition "Definition 28.7 (Bell states)"
    The four **Bell states** form an orthonormal basis for $\mathbb{C}^2 \otimes \mathbb{C}^2$:

    $$
    |\Phi^{\pm}\rangle = \frac{1}{\sqrt{2}}(|00\rangle \pm |11\rangle), \quad |\Psi^{\pm}\rangle = \frac{1}{\sqrt{2}}(|01\rangle \pm |10\rangle).
    $$

    Bell states are **maximally entangled**: each subsystem's reduced density matrix is $\rho_A = \rho_B = \frac{I}{2}$ (maximally mixed state), and the von Neumann entropy achieves its maximum value $S = \log 2 = 1$ bit.

!!! theorem "Theorem 28.4 (Schmidt decomposition theorem)"
    Any bipartite pure state $|\psi\rangle \in \mathcal{H}_A \otimes \mathcal{H}_B$ ($\dim \mathcal{H}_A = d_A$, $\dim \mathcal{H}_B = d_B$) can be written as

    $$
    |\psi\rangle = \sum_{i=1}^{r} \sqrt{\lambda_i}\, |u_i\rangle_A \otimes |v_i\rangle_B,
    $$

    where $r \le \min(d_A, d_B)$ is the **Schmidt rank**, $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_r > 0$ are the **Schmidt coefficients** ($\sum_{i=1}^r \lambda_i = 1$), and $\{|u_i\rangle\}$, $\{|v_i\rangle\}$ are orthonormal sets in $\mathcal{H}_A$ and $\mathcal{H}_B$, respectively. The Schmidt coefficients are uniquely determined; when they are all distinct, the vectors are also unique up to phase factors.

    $|\psi\rangle$ is a product state $\iff$ Schmidt rank $r = 1$.

??? proof "Proof"
    Arrange the coefficients of $|\psi\rangle = \sum_{j,k} c_{jk} |j\rangle_A |k\rangle_B$ into a matrix $C \in \mathbb{C}^{d_A \times d_B}$, where $C_{jk} = c_{jk}$. Perform the singular value decomposition (SVD, Ch11): $C = U\Sigma V^\dagger$, where $U \in \mathbb{C}^{d_A \times d_A}$ and $V \in \mathbb{C}^{d_B \times d_B}$ are unitary, and $\Sigma$ is diagonal with entries $\sigma_i \ge 0$.

    Then $c_{jk} = \sum_{i=1}^{r} \sigma_i U_{ji} V_{ki}^*$. Define

    $$
    |u_i\rangle_A = \sum_j U_{ji} |j\rangle_A, \quad |v_i\rangle_B = \sum_k V_{ki}^* |k\rangle_B.
    $$

    The unitarity of $U$ and $V$ ensures that $\{|u_i\rangle\}$ and $\{|v_i\rangle\}$ are orthonormal. Normalization gives $\sum_i \sigma_i^2 = \|C\|_F^2 = \langle\psi|\psi\rangle = 1$. Setting $\lambda_i = \sigma_i^2$:

    $$
    |\psi\rangle = \sum_{i=1}^r \sigma_i |u_i\rangle_A |v_i\rangle_B = \sum_{i=1}^r \sqrt{\lambda_i} |u_i\rangle_A |v_i\rangle_B.
    $$

    Uniqueness follows from the uniqueness of singular values in the SVD. $\blacksquare$

!!! example "Example 28.4"
    **Schmidt decomposition and entanglement detection for Bell states.**

    $|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$. The coefficient matrix is $C = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} = \frac{1}{\sqrt{2}}I$. SVD: $U = V = I$, $\sigma_1 = \sigma_2 = \frac{1}{\sqrt{2}}$. Schmidt rank $r = 2 > 1$, so $|\Phi^+\rangle$ is entangled. The Schmidt coefficients $\lambda_1 = \lambda_2 = \frac{1}{2}$ (uniform), confirming maximal entanglement.

    Comparison: $|\psi\rangle = |01\rangle$. Coefficient matrix $C = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$. SVD gives $\sigma_1 = 1$, Schmidt rank $r = 1$, so $|01\rangle = |0\rangle \otimes |1\rangle$ is a product state (separable).

---

## 28.5 Quantum Measurement and POVM

<div class="context-flow" markdown>

**Projective measurement**: Orthogonal decomposition $I = \sum_m P_m$ → Born rule $p(m) = \langle\psi|P_m|\psi\rangle$ → post-measurement state collapse → POVM generalizes to non-orthogonal measurements → Naimark dilation theorem
**Link**: Direct application of orthogonal projections and spectral decomposition from Ch7

</div>

Quantum measurement combines the theory of projection operators from linear algebra with probability theory.

!!! definition "Definition 28.8 (Projective measurement and the Born rule)"
    A **projective measurement** is defined by a set of orthogonal projectors $\{P_m\}$ satisfying

    $$
    P_m P_{m'} = \delta_{mm'} P_m, \quad \sum_m P_m = I, \quad P_m = P_m^\dagger.
    $$

    **Born rule**: For a state $|\psi\rangle$, the probability of measurement outcome $m$ is

    $$
    p(m) = \langle\psi|P_m|\psi\rangle = \operatorname{tr}(P_m |\psi\rangle\langle\psi|).
    $$

    After measurement, the state collapses to $\frac{P_m|\psi\rangle}{\sqrt{p(m)}}$. For a density matrix $\rho$, $p(m) = \operatorname{tr}(P_m \rho)$.

!!! definition "Definition 28.9 (POVM)"
    A **Positive Operator-Valued Measure (POVM)** is a generalized measurement framework. A POVM is defined by a set of **positive operators** $\{E_m\}$ satisfying

    $$
    E_m \succeq 0, \quad \sum_m E_m = I.
    $$

    The probability of outcome $m$ is $p(m) = \operatorname{tr}(E_m \rho)$. POVM elements $E_m$ need not be projectors, so there can be more measurement outcomes than the dimension of the space. A POVM describes information extraction but does not specify the post-measurement state.

!!! theorem "Theorem 28.5 (Naimark dilation theorem)"
    Any POVM $\{E_m\}_{m=1}^M$ acting on $\mathcal{H}$ can be realized as a projective measurement on an extended space $\mathcal{H} \otimes \mathcal{H}_{\text{aux}}$. That is, there exist an auxiliary space $\mathcal{H}_{\text{aux}}$, an initial auxiliary state $|0\rangle_{\text{aux}}$, and orthogonal projectors $\{\Pi_m\}$ on the extended space such that

    $$
    E_m = \langle 0|_{\text{aux}} \Pi_m |0\rangle_{\text{aux}}, \quad \forall m.
    $$

    Equivalently, $\operatorname{tr}(E_m \rho) = \operatorname{tr}(\Pi_m (\rho \otimes |0\rangle\langle 0|))$.

??? proof "Proof"
    Construct $\mathcal{H}_{\text{aux}} = \mathbb{C}^M$. For each POVM element $E_m$, since $E_m \succeq 0$, there exists $A_m$ such that $E_m = A_m^\dagger A_m$ (take $A_m = \sqrt{E_m}$). Define the isometry $V: \mathcal{H} \to \mathcal{H} \otimes \mathcal{H}_{\text{aux}}$:

    $$
    V|\psi\rangle = \sum_{m=1}^M (A_m|\psi\rangle) \otimes |m\rangle.
    $$

    Verify that $V$ is an isometry: $\langle\psi|V^\dagger V|\psi\rangle = \sum_m \langle\psi|A_m^\dagger A_m|\psi\rangle = \langle\psi|\sum_m E_m|\psi\rangle = \langle\psi|\psi\rangle$.

    Define $\Pi_m = I_{\mathcal{H}} \otimes |m\rangle\langle m|$ (projector on the extended space). Then

    $$
    \langle\psi|V^\dagger \Pi_m V|\psi\rangle = \langle\psi|A_m^\dagger A_m|\psi\rangle = \langle\psi|E_m|\psi\rangle = p(m).
    $$

    Thus POVM measurement is equivalent to embedding into an extended space followed by a projective measurement. $\blacksquare$

!!! example "Example 28.5"
    **Three-element POVM for distinguishing non-orthogonal states.** Consider receiving $|0\rangle$ or $|+\rangle = \frac{1}{\sqrt{2}}(|0\rangle + |1\rangle)$ with equal probability $1/2$. Since $\langle 0|+\rangle = \frac{1}{\sqrt{2}} \ne 0$, the two states are non-orthogonal and cannot be perfectly distinguished by a projective measurement.

    Construct a three-element POVM for unambiguous state discrimination (USD) $\{E_0, E_+, E_?\}$:

    - $E_0 = \frac{1}{1 + 1/\sqrt{2}} |1\rangle\langle 1|$: when triggered, the state is definitely $|+\rangle$ (since $\langle 1|0\rangle = 0$, so $\operatorname{tr}(E_0 |0\rangle\langle 0|) = 0$).
    - $E_+ = \frac{1}{1 + 1/\sqrt{2}} |-\rangle\langle -|$: when triggered, the state is definitely $|0\rangle$ (since $\langle -|+\rangle = 0$).
    - $E_? = I - E_0 - E_+$: inconclusive result.

    The optimal unambiguous discrimination success probability is $p_{\text{succ}} = 1 - |\langle 0|+\rangle| = 1 - \frac{1}{\sqrt{2}} \approx 0.293$.

---

## 28.6 Density Matrices and Quantum Channels

<div class="context-flow" markdown>

**Open systems**: Mixed state $\rho$ → partial trace (reduced density matrix) → quantum channel $\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger$ (Kraus representation) → CPTP maps → Choi-Kraus representation theorem
**Link**: Extension of positive definite matrices / PSD cone from Ch16 and tensor products from Ch21

</div>

Quantum channels describe the evolution of open quantum systems, and their mathematical essence is completely positive trace-preserving (CPTP) maps.

!!! definition "Definition 28.10 (Partial trace and reduced density matrix)"
    For a composite system state $\rho_{AB} \in \mathcal{B}(\mathcal{H}_A \otimes \mathcal{H}_B)$, the **reduced density matrix** of subsystem $A$ is

    $$
    \rho_A = \operatorname{tr}_B(\rho_{AB}) = \sum_j (I_A \otimes \langle j|_B) \rho_{AB} (I_A \otimes |j\rangle_B),
    $$

    where $\{|j\rangle_B\}$ is any orthonormal basis for $\mathcal{H}_B$. The partial trace is the unique operation satisfying: for all operators $O_A$ on $A$,

    $$
    \operatorname{tr}(O_A \rho_A) = \operatorname{tr}((O_A \otimes I_B)\rho_{AB}).
    $$

!!! definition "Definition 28.11 (Quantum channel and Kraus representation)"
    A **quantum channel** is a completely positive trace-preserving (CPTP) map $\mathcal{E}: \mathcal{B}(\mathcal{H}_{\text{in}}) \to \mathcal{B}(\mathcal{H}_{\text{out}})$. Its **Kraus representation** is

    $$
    \mathcal{E}(\rho) = \sum_{k=1}^{r} K_k \rho K_k^\dagger, \quad \sum_{k=1}^{r} K_k^\dagger K_k = I.
    $$

    The operators $\{K_k\}$ are called **Kraus operators**. The condition $\sum_k K_k^\dagger K_k = I$ ensures trace preservation ($\operatorname{tr}(\mathcal{E}(\rho)) = \operatorname{tr}(\rho) = 1$).

    Common quantum channels include:

    - **Depolarizing channel**: $\mathcal{E}(\rho) = (1-p)\rho + \frac{p}{d}I$, replacing the state with the maximally mixed state with probability $p$.
    - **Amplitude damping channel**: $K_0 = \begin{pmatrix} 1 & 0 \\ 0 & \sqrt{1-\gamma} \end{pmatrix}$, $K_1 = \begin{pmatrix} 0 & \sqrt{\gamma} \\ 0 & 0 \end{pmatrix}$, modeling spontaneous emission.

!!! theorem "Theorem 28.6 (Choi-Kraus representation theorem)"
    A linear map $\mathcal{E}: \mathcal{B}(\mathcal{H}_{\text{in}}) \to \mathcal{B}(\mathcal{H}_{\text{out}})$ is completely positive if and only if there exist operators $\{K_k\}_{k=1}^r$ such that

    $$
    \mathcal{E}(\rho) = \sum_{k=1}^r K_k \rho K_k^\dagger.
    $$

    Equivalently, $\mathcal{E}$ is completely positive if and only if its **Choi matrix**

    $$
    J(\mathcal{E}) = \sum_{i,j} |i\rangle\langle j| \otimes \mathcal{E}(|i\rangle\langle j|) \in \mathcal{B}(\mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}})
    $$

    is positive semidefinite. Additionally, $\mathcal{E}$ is trace-preserving if and only if $\operatorname{tr}_{\text{out}}(J(\mathcal{E})) = I_{\text{in}}$.

??? proof "Proof"
    **Kraus $\Rightarrow$ CP**: For any reference space $\mathcal{H}_R$ and $\rho_{RA} \succeq 0$,

    $$
    (\operatorname{id}_R \otimes \mathcal{E})(\rho_{RA}) = \sum_k (I_R \otimes K_k) \rho_{RA} (I_R \otimes K_k)^\dagger \succeq 0.
    $$

    Each term $X\rho X^\dagger \succeq 0$ (preserves positive semidefiniteness), so the sum is also positive semidefinite.

    **CP $\Rightarrow$ Kraus**: Suppose $\mathcal{E}$ is completely positive. Then the Choi matrix $J(\mathcal{E}) \succeq 0$. Perform a spectral decomposition $J(\mathcal{E}) = \sum_k |\phi_k\rangle\langle\phi_k|$, and reshape each $|\phi_k\rangle \in \mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}}$ into an operator $K_k: \mathcal{H}_{\text{in}} \to \mathcal{H}_{\text{out}}$ (via the natural isomorphism $\mathcal{H}_{\text{in}} \otimes \mathcal{H}_{\text{out}} \cong \operatorname{Hom}(\mathcal{H}_{\text{in}}, \mathcal{H}_{\text{out}})$). One can verify that $\mathcal{E}(\rho) = \sum_k K_k \rho K_k^\dagger$.

    The trace-preservation condition $\sum_k K_k^\dagger K_k = I$ is equivalent to $\operatorname{tr}_{\text{out}}(J(\mathcal{E})) = I_{\text{in}}$, which can be verified by direct computation. $\blacksquare$

!!! example "Example 28.6"
    **The dephasing channel.** The single-qubit dephasing channel $\mathcal{E}(\rho) = (1-p)\rho + pZ\rho Z$ has Kraus operators $K_0 = \sqrt{1-p}\, I$, $K_1 = \sqrt{p}\, Z$.

    Acting on a general state $\rho = \begin{pmatrix} a & b \\ b^* & 1-a \end{pmatrix}$:

    $$
    \mathcal{E}(\rho) = (1-p)\begin{pmatrix} a & b \\ b^* & 1-a \end{pmatrix} + p\begin{pmatrix} a & -b \\ -b^* & 1-a \end{pmatrix} = \begin{pmatrix} a & (1-2p)b \\ (1-2p)b^* & 1-a \end{pmatrix}.
    $$

    The diagonal elements are unchanged, while the off-diagonal elements decay by a factor of $(1-2p)$. When $p = 1/2$, the off-diagonal elements vanish and coherence is completely lost. The Choi matrix is

    $$
    J(\mathcal{E}) = \begin{pmatrix} 1 & 0 & 0 & 1-2p \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 1-2p & 0 & 0 & 1 \end{pmatrix},
    $$

    with eigenvalues $\{1 + (1-2p), 1 - (1-2p), 0, 0\} = \{2-2p, 2p, 0, 0\} \ge 0$, confirming complete positivity.

---

## 28.7 Quantum Error Correction

<div class="context-flow" markdown>

**Core idea**: Logical qubits encoded into subspace of physical qubits → errors = linear operators → Knill-Laflamme conditions = orthogonality conditions on the code space → stabilizer codes → CSS codes
**Link**: Core application of subspaces from Ch4 and projections from Ch7

</div>

Quantum error correction is the foundation of fault-tolerant quantum computation.

!!! definition "Definition 28.12 (Quantum error-correcting code)"
    An $[[n, k, d]]$ **quantum error-correcting code** encodes $k$ logical qubits into $n$ physical qubits. The code space $\mathcal{C}$ is a $2^k$-dimensional subspace of $(\mathbb{C}^2)^{\otimes n}$. The code distance $d$ is the maximum error weight that can be detected plus $1$; the code can correct $\lfloor(d-1)/2\rfloor$ arbitrary single-qubit errors.

!!! theorem "Theorem 28.7 (Knill-Laflamme quantum error correction conditions)"
    A quantum code $\mathcal{C}$ (with code space projector $P$) can correct an error set $\{E_a\}$ if and only if there exists a Hermitian matrix $(\alpha_{ab})$ such that

    $$
    P E_a^\dagger E_b P = \alpha_{ab} P, \quad \forall a, b.
    $$

    Equivalently, for any orthonormal basis $\{|\psi_i\rangle\}$ of the code space,

    $$
    \langle\psi_i|E_a^\dagger E_b|\psi_j\rangle = \alpha_{ab}\delta_{ij}, \quad \forall a, b, i, j.
    $$

    Intuitive meaning: different errors must either map the code space to orthogonal subspaces (distinguishable), or act proportionally on the code space (equivalent errors).

??? proof "Proof"
    **Sufficiency**: The Knill-Laflamme conditions guarantee that $\{E_a P\}$ maps the code space into distinguishable subspaces. Unitarily diagonalize the Hermitian matrix $\alpha = (\alpha_{ab})$ as $\alpha = W D W^\dagger$, and define "canonical errors" $\tilde{E}_c = \sum_a W_{ac} E_a$. Then

    $$
    P \tilde{E}_c^\dagger \tilde{E}_{c'} P = d_c \delta_{cc'} P,
    $$

    meaning that canonical errors map the code space to mutually orthogonal subspaces. Error correction consists of projecting onto the error subspace (identifying which error occurred) + inverse rotation (restoring the original state).

    **Necessity**: If error correction is possible, there exists a recovery operation $\mathcal{R}$ such that $\mathcal{R}(\sum_a E_a \rho E_a^\dagger) = \rho$ for all $\rho$ supported on the code space. Letting $\{R_l\}$ be the Kraus operators of $\mathcal{R}$, the recovery condition implies that $\langle\psi_i|E_a^\dagger E_b|\psi_j\rangle$ is independent of $i, j$ (when $i = j$) and zero when $i \ne j$. $\blacksquare$

!!! definition "Definition 28.13 (Stabilizer codes and CSS codes)"
    A **stabilizer code** is defined by an abelian subgroup $\mathcal{S}$ (the stabilizer group) of the $n$-qubit Pauli group. The code space is

    $$
    \mathcal{C} = \{|\psi\rangle : S|\psi\rangle = |\psi\rangle,\, \forall S \in \mathcal{S}\}.
    $$

    If $\mathcal{S}$ is generated by $n - k$ independent generators, then $\dim \mathcal{C} = 2^k$.

    **CSS codes** (Calderbank-Shor-Steane) are a special class of stabilizer codes constructed from two classical linear codes $C_1 \supset C_2$, such that $X$ errors and $Z$ errors can be corrected independently.

!!! example "Example 28.7"
    **Shor's 9-qubit code and the Steane 7-qubit code.**

    **Shor's $[[9,1,3]]$ code** encodes $1$ logical qubit into $9$ physical qubits:

    $$
    |0_L\rangle = \frac{1}{2\sqrt{2}}(|000\rangle + |111\rangle)^{\otimes 3}, \quad |1_L\rangle = \frac{1}{2\sqrt{2}}(|000\rangle - |111\rangle)^{\otimes 3}.
    $$

    This code can correct arbitrary single-qubit errors. Structure: the outer 3-qubit repetition code corrects bit-flip ($X$) errors, and the inner 3-qubit phase code corrects phase-flip ($Z$) errors.

    **The Steane $[[7,1,3]]$ code** is a CSS code based on the classical $[7,4,3]$ Hamming code. Since the Hamming code is self-dual ($C_2 = C_1^\perp \subset C_1$), the Steane code requires only $7$ physical qubits to correct arbitrary single-qubit errors, making it more compact than Shor's code.

    Verifying the Knill-Laflamme conditions: for the Steane code, for any two Pauli errors $E_a, E_b$ of weight $\le 1$, $E_a^\dagger E_b$ has weight $\le 2 < d = 3$, so $P E_a^\dagger E_b P = \alpha_{ab} P$.

!!! example "Example 28.8"
    **Structure of stabilizer codes.** The $[[5,1,3]]$ code is the smallest quantum code that can correct arbitrary single-qubit errors. Its stabilizer group is generated by the following $4 = n - k$ generators:

    $$
    g_1 = XZZXI, \quad g_2 = IXZZX, \quad g_3 = XIXZZ, \quad g_4 = ZXIXZ.
    $$

    The code space consists of all states satisfying $g_i|\psi\rangle = |\psi\rangle$ ($i = 1, 2, 3, 4$), with dimension $2^1 = 2$.

    Code distance $d = 3$: any Pauli operator of weight $\le 2$ either belongs to the stabilizer group or anticommutes with some generator, so it can be detected. The $5$-qubit code saturates the quantum Singleton bound $n - k \ge 2(d-1)$, i.e., $4 \ge 4$.

---

## 28.8 Matrix Inequalities in Quantum Information

<div class="context-flow" markdown>

**Core tools**: Von Neumann entropy $S(\rho) = -\operatorname{tr}(\rho \log \rho)$ → quantum relative entropy → strong subadditivity → data processing inequality
**Link**: Quantum information extensions of matrix inequalities from Ch18

</div>

The core inequalities in quantum information theory can all be formulated as inequalities on matrix functions.

!!! definition "Definition 28.14 (Von Neumann entropy)"
    The **von Neumann entropy** of a density matrix $\rho$ is defined as

    $$
    S(\rho) = -\operatorname{tr}(\rho \log \rho) = -\sum_i \lambda_i \log \lambda_i,
    $$

    where $\{\lambda_i\}$ are the eigenvalues of $\rho$, with the convention $0 \log 0 = 0$.

    - $S(\rho) = 0$ if and only if $\rho$ is a pure state.
    - $S(\rho) \le \log d$ ($d = \dim \mathcal{H}$), with equality if and only if $\rho = I/d$ (maximally mixed state).

!!! definition "Definition 28.15 (Quantum relative entropy)"
    The **quantum relative entropy** between two density matrices $\rho$ and $\sigma$ is

    $$
    S(\rho \| \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma),
    $$

    which is finite when $\ker(\sigma) \cap \operatorname{supp}(\rho) = \{0\}$, and defined to be $+\infty$ otherwise.

!!! theorem "Theorem 28.8 (Strong subadditivity)"
    For any density matrix $\rho_{ABC}$ of a tripartite system $ABC$, the von Neumann entropy satisfies **strong subadditivity** (SSA):

    $$
    S(\rho_{ABC}) + S(\rho_B) \le S(\rho_{AB}) + S(\rho_{BC}),
    $$

    where $\rho_{AB} = \operatorname{tr}_C(\rho_{ABC})$, etc., are reduced density matrices.

    Equivalent forms:

    - **Conditional entropy decrease**: $S(A|BC) \le S(A|B)$, where $S(A|B) = S(\rho_{AB}) - S(\rho_B)$.
    - **Non-negativity of conditional mutual information**: $I(A:C|B) = S(A|B) - S(A|BC) \ge 0$.

??? proof "Proof"
    The standard proof of strong subadditivity uses the Lieb concavity theorem. Define the map $f(X) = \operatorname{tr}(\exp(\log M + \log X))$ ($M$ a fixed positive definite matrix); Lieb proved that $f$ is concave.

    An alternative proof is based on the monotonicity of quantum relative entropy (data processing inequality). The quantum relative entropy is monotonically decreasing under quantum channels (CPTP maps):

    $$
    S(\mathcal{E}(\rho) \| \mathcal{E}(\sigma)) \le S(\rho \| \sigma).
    $$

    Taking $\mathcal{E} = \operatorname{tr}_C$ (the partial trace is a CPTP map), $\rho = \rho_{ABC}$, $\sigma = I_A \otimes \rho_{BC} / d_A$, and simplifying yields strong subadditivity.

    Specifically, $S(\rho_{ABC} \| I_A/d_A \otimes \rho_{BC}) = \log d_A - S(A|BC)$ and $S(\rho_{AB} \| I_A/d_A \otimes \rho_B) = \log d_A - S(A|B)$. Monotonicity under partial trace gives $\log d_A - S(A|BC) \ge \log d_A - S(A|B)$, i.e., $S(A|B) \ge S(A|BC)$, which is strong subadditivity. $\blacksquare$

!!! theorem "Theorem 28.9 (Klein's inequality and non-negativity of quantum relative entropy)"
    For any density matrices $\rho$ and $\sigma$,

    $$
    S(\rho \| \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma) \ge 0,
    $$

    with equality if and only if $\rho = \sigma$.

??? proof "Proof"
    This is a corollary of Klein's inequality. Klein's inequality states that for a convex function $f$ and self-adjoint operators $A, B$:

    $$
    \operatorname{tr}(f(A) - f(B) - f'(B)(A - B)) \ge 0.
    $$

    Take $f(t) = t \log t$ (a convex function), $A = \rho$, $B = \sigma$, $f'(t) = \log t + 1$. Then

    $$
    \operatorname{tr}(\rho \log \rho - \sigma \log \sigma - (\log \sigma + I)(\rho - \sigma)) \ge 0.
    $$

    Expanding and using $\operatorname{tr}(\rho) = \operatorname{tr}(\sigma) = 1$:

    $$
    \operatorname{tr}(\rho \log \rho - \sigma \log \sigma - \rho \log \sigma - \rho + \sigma \log \sigma + \sigma) = \operatorname{tr}(\rho \log \rho - \rho \log \sigma) = S(\rho \| \sigma) \ge 0.
    $$

    The equality condition follows from the strict convexity of $f$: $\rho = \sigma$. $\blacksquare$

!!! example "Example 28.9"
    **Verifying strong subadditivity: GHZ state.** The three-qubit GHZ state $|\text{GHZ}\rangle = \frac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$.

    Reduced density matrices:

    - $\rho_{ABC} = |\text{GHZ}\rangle\langle\text{GHZ}|$ (pure state), $S(\rho_{ABC}) = 0$.
    - $\rho_{AB} = \frac{1}{2}(|00\rangle\langle 00| + |11\rangle\langle 11|)$, $S(\rho_{AB}) = \log 2 = 1$.
    - $\rho_{BC} = \frac{1}{2}(|00\rangle\langle 00| + |11\rangle\langle 11|)$, $S(\rho_{BC}) = 1$.
    - $\rho_B = \frac{I}{2}$, $S(\rho_B) = 1$.

    SSA verification: $S(\rho_{ABC}) + S(\rho_B) = 0 + 1 = 1 \le S(\rho_{AB}) + S(\rho_{BC}) = 1 + 1 = 2$. The inequality holds.

    Conditional mutual information $I(A:C|B) = S(A|B) - S(A|BC) = (S(\rho_{AB}) - S(\rho_B)) - (S(\rho_{ABC}) - S(\rho_{BC})) = (1-1) - (0-1) = 1 \ge 0$.

!!! note "Note"
    The linear algebra tools in this chapter permeate every area of quantum information: Hilbert spaces provide the framework for state spaces (Ch8), unitary matrices describe quantum evolution (Ch7), SVD yields the Schmidt decomposition (Ch11), Kronecker products describe multi-body systems (Ch19), positive definite matrix theory characterizes density matrices and POVMs (Ch16), and matrix inequalities constrain the limits of information processing (Ch18). Quantum information science is one of the deepest and most active application areas of linear algebra.
