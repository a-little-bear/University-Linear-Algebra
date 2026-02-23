# Chapter 23: Random Matrices

<div class="context-flow" markdown>

**Prerequisites**: Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16) · Probability Theory · Statistics

**Chapter Outline**: Introduction to Random Matrix Theory (RMT) → Gaussian Ensembles (GOE, GUE, GSE) → Wigner’s Semicircle Law → Wishart Matrices and Sample Covariance → Marchenko-Pastur Law → Tracy-Widom Distribution & Extreme Eigenvalues → Level Repulsion → Applications: High-Dimensional Statistics, Quantum Physics, Deep Learning, and Compressed Sensing

**Extension**: RMT is the "Linear Algebra of Chaos"; it describes the universal laws governing the spectrum of large-scale systems where individual entries are unknown but follow statistical distributions.

</div>

In traditional linear algebra, matrices are deterministic. However, in modern big data, quantum physics, and neural network analysis, we often deal with matrices where the entries are random variables. **Random Matrix Theory** (RMT) studies the behavior of eigenvalues and eigenvectors of such matrices as their dimension $n \to \infty$. The most remarkable discovery of RMT is **Universality**: the statistical properties of the spectrum often depend only on the symmetry of the matrix, not on the specific distribution of its entries.

---

## 23.1 Gaussian Ensembles and Symmetry

<div class="context-flow" markdown>

**Motivation**: To study random matrices, we classify them by their symmetry groups. This leads to the three classical "Gaussian Ensembles."

</div>

!!! definition "Definition 23.1 (Gaussian Orthogonal Ensemble - GOE)"
    The **GOE** consists of real symmetric matrices $A = (A + A^T)/2$ where entries are i.i.d. Gaussian. It is invariant under orthogonal transformations ($A \mapsto Q^T A Q$).

!!! definition "Definition 23.2 (Gaussian Unitary Ensemble - GUE)"
    The **GUE** consists of complex Hermitian matrices invariant under unitary transformations. It describes systems without time-reversal symmetry.

!!! definition "Definition 23.3 (Gaussian Symplectic Ensemble - GSE)"
    The **GSE** consists of quaternionic Hermitian matrices invariant under symplectic transformations.

---

## 23.2 Wigner's Semicircle Law

<div class="context-flow" markdown>

**The Bulk of the Spectrum**: What is the overall shape of the eigenvalue distribution for a large random symmetric matrix?

</div>

!!! theorem "Theorem 23.1 (Wigner's Semicircle Law)"
    Let $A$ be an $n \times n$ symmetric matrix where $a_{ij}$ are i.i.d. with mean 0 and variance $\sigma^2$. As $n \to \infty$, the empirical spectral distribution of $A/\sqrt{n}$ converges to the **semicircle distribution**:
    $$\rho(x) = \begin{cases} \frac{1}{2\pi \sigma^2} \sqrt{4\sigma^2 - x^2} & \text{if } |x| \le 2\sigma \\ 0 & \text{otherwise} \end{cases}$$
    This shows that eigenvalues are not uniformly distributed but are more concentrated near the origin.

---

## 23.3 Wishart Matrices and Marchenko-Pastur

<div class="context-flow" markdown>

**Statistics Perspective**: In data science, we care about the eigenvalues of the **sample covariance matrix** $S = \frac{1}{n} XX^T$. These are called Wishart matrices.

</div>

!!! theorem "Theorem 23.2 (Marchenko-Pastur Law)"
    Let $X$ be an $n \times p$ matrix with i.i.d. entries (mean 0, variance $\sigma^2$). Let $n, p \to \infty$ such that the aspect ratio $p/n \to \gamma \in (0, \infty)$. The distribution of eigenvalues of $S = \frac{1}{n} XX^T$ converges to:
    $$\rho_{\gamma}(x) = \frac{1}{2\pi \sigma^2 \gamma x} \sqrt{(b-x)(x-a)}$$
    where $a = \sigma^2(1-\sqrt{\gamma})^2$ and $b = \sigma^2(1+\sqrt{\gamma})^2$.
    *Note: If $\gamma > 1$, there is also a point mass of size $1 - 1/\gamma$ at $x = 0$.*

---

## 23.4 Extreme Eigenvalues and Tracy-Widom

!!! theorem "Theorem 23.3 (Tracy-Widom Distribution)"
    The fluctuations of the largest eigenvalue $\lambda_{\max}$ of a GUE matrix do not follow a Gaussian distribution. Instead, they follow the **Tracy-Widom distribution**, which has a significantly fatter tail on the left. This distribution appears universally in finance, growth models, and the length of the longest increasing subsequence.

---

## 23.5 Level Repulsion

!!! technique "Repulsion: Eigenvalues Hate Each Other"
    In a random matrix, eigenvalues are "repelled" from one another. The probability of finding two eigenvalues very close together ($s \approx 0$) is nearly zero. This is in stark contrast to independent random variables (Poisson process) where clustering is common. This "Level Repulsion" is a hallmark of quantum chaos.

---

## Exercises

1.  **[Semicircle] Calculate the support of the semicircle distribution for $A/\sqrt{n}$ if the entries have variance $\sigma^2 = 1$.**
    ??? success "Solution"
        The support is $[-2\sigma, 2\sigma] = [-2, 2]$.

2.  **[M-P Law] If we have 1000 samples ($n=1000$) and 500 features ($p=500$), what is the aspect ratio $\gamma$ and the range of eigenvalues for the sample covariance matrix (assuming $\sigma=1$)?**
    ??? success "Solution"
        $\gamma = p/n = 0.5$.
        $a = (1 - \sqrt{0.5})^2 \approx (1 - 0.707)^2 \approx 0.086$.
        $b = (1 + \sqrt{0.5})^2 \approx (1 + 0.707)^2 \approx 2.91$.
        The spectrum is spread across $[0.086, 2.91]$.

3.  **[Spikes] What happens to the M-P law if the true covariance matrix is not $I$ but has one very large eigenvalue (a "spike")?**
    ??? success "Solution"
        If the spike is large enough (above the BBP threshold), one eigenvalue will "pop out" of the bulk distribution $[a, b]$ and become visible. This is the basis for signal detection in noise.

4.  **[Universality] Does the Semicircle Law hold if the entries are not Gaussian (e.g., Uniform or Bernoulli)?**
    ??? success "Solution"
        Yes. As long as the entries are independent and have finite second moments, the semicircle law holds. This is the **Universality** of RMT.

5.  **[Trace] Use the Semicircle Law to estimate $\frac{1}{n} \operatorname{tr}(A^2)$ for large $n$ (assuming $\sigma=1$).**
    ??? success "Solution"
        The average of $\lambda^2$ is $\int_{-2}^2 x^2 \rho(x) dx$. For the semicircle distribution, the second moment is $\sigma^2 = 1$. Thus, $\operatorname{tr}(A^2) \approx n$.

6.  **[Invariance] Why is the GUE called "Unitary"?**
    ??? success "Solution"
        Because its probability density $P(A) \propto \exp(-\operatorname{tr}(A^2))$ is invariant under the change of basis $A \mapsto UAU^*$, where $U$ is a unitary matrix.

7.  **[RMT in AI] How is RMT used in deep learning?**
    ??? success "Solution"
        It is used to analyze the initialization of weights. If the eigenvalues of the weight matrices are outside the "safe" range, gradients will either explode or vanish. RMT helps design "Orthogonal Initialization" to keep the spectrum stable.

8.  **[Wishart] Is a Wishart matrix always positive semi-definite?**
    ??? success "Solution"
        Yes. Since $S = \frac{1}{n} XX^T$, the quadratic form $v^T S v = \frac{1}{n} \|X^T v\|^2 \ge 0$.

9.  **[Compressed Sensing] How does RMT relate to the Restricted Isometry Property (RIP)?**
    ??? success "Solution"
        RMT proves that random Gaussian or Bernoulli matrices satisfy RIP with high probability, meaning they act like an isometry on sparse vectors.

10. **[Tracy-Widom] Compare Tracy-Widom to the Normal distribution.**

   ??? success "Solution"
        Tracy-Widom is asymmetric and has a much steeper decay on the right and a slower decay on the left. It describes the "edge" of the spectrum, whereas Normal describes the sum of many independent events.

## Chapter Summary

Random Matrix Theory reveals the order hidden within high-dimensional chaos:

1.  **Emergent Shapes**: Demonstrated that while individual eigenvalues are random, their collective density follows rigid, predictable shapes like the Semicircle or Marchenko-Pastur laws.
2.  **Universal Laws**: Highlighted that spectral properties of large systems are often independent of the fine details of their components, depending only on their symmetry.
3.  **Interaction Dynamics**: Identified level repulsion as the key difference between random matrices and simple random sequences, linking linear algebra to quantum physics.
4.  **Practical Diagnostic**: Established RMT as a tool for distinguishing between true signals and random noise in high-dimensional data analysis.
