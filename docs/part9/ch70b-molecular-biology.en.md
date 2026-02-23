# Chapter 70B: Molecular Biology and Genomics

<div class="context-flow" markdown>

**Prerequisites**: Probability Theory · Matrix Exponential (Ch13) · Markov Chains (Ch71) · Vector Spaces (Ch04)

**Chapter Outline**: Digital Encoding of Genetic Information → Vector Representation of DNA Bases $\{A, C, G, T\}$ → Substitution Probability Matrices → Core Model: The Jukes-Cantor (JC69) Model and its Algebraic Structure → The Kimura Two-Parameter Model → Algebraic Derivation of Evolutionary Distance via Matrix Logarithms → Linear Scoring Functions in Sequence Alignment → Applications: Phylogenetic Tree Construction, Paternity Testing, and Clustering of Gene Expression Profiles

**Extension**: Genomics is the "informational reconstruction" of linear algebra; it simplifies the history of life's evolution into the cumulative product of mutation matrices over millions of years. It proves that differences between species can be quantified as the distance between operators in a state space—the algebraic engine for understanding molecular evolution and precision medicine.

</div>

The code of life is a sequence woven from four bases: A, C, G, and T. At the molecular level, evolution is essentially the stochastic substitution of these bases over time. **Linear Algebra** provides the framework to describe this probabilistic evolution via **Substitution Matrices**. Using matrix exponentials and logarithms, we can infer common ancestors from millions of years ago based on current sequence discrepancies. This chapter introduces the algebraic theory at the heart of computational biology and evolutionary genomics.

---

## 70B.1 Base States and Substitution Matrices

!!! definition "Definition 70B.1 (State Space)"
    The state of a DNA site is represented by 4D basis vectors:
    $|A\rangle = (1,0,0,0)^T, |C\rangle = (0,1,0,0)^T, |G\rangle = (0,0,1,0)^T, |T\rangle = (0,0,0,1)^T$.
    A substitution matrix $M(t)$ has entries $m_{ij}$ representing the probability that base $j$ changes to base $i$ over time $t$.

---

## 70B.2 The Jukes-Cantor (JC69) Model

!!! definition "Definition 70B.2 (JC69 Model)"
    Assuming equal mutation rates among all bases, the rate matrix $Q$ has a symmetric structure:
    $$Q = \begin{pmatrix} -3\alpha & \alpha & \alpha & \alpha \\ \alpha & -3\alpha & \alpha & \alpha \\ \alpha & \alpha & -3\alpha & \alpha \\ \alpha & \alpha & \alpha & -3\alpha \end{pmatrix}$$
    The evolution matrix is $M(t) = e^{Qt}$.

---

## 70B.3 Algebraic Calculation of Evolutionary Distance

!!! theorem "Theorem 70B.1 (Distance Formula)"
    If the observed proportion of mismatched sites between two sequences is $p$, the evolutionary distance $d$ under the JC69 model is:
    $$d = -\frac{3}{4} \ln(1 - \frac{4}{3}p)$$
    **Significance**: This logarithmic formula corrects for underestimation caused by multiple mutations returning a site to its original state (back-mutations).

---

## Exercises

**1. [Basics] Prove that each row of the JC69 rate matrix $Q$ sums to 0.**

??? success "Solution"
    **Proof:**
    1. Each row contains one $-3\alpha$ and three $\alpha$s.
    2. Total sum $S = -3\alpha + \alpha + \alpha + \alpha = 0$.
    **Physical Meaning**: This ensures the rows of the transition matrix $M(t) = e^{Qt}$ sum to 1 (probability conservation), consistent with the definition of a stochastic matrix.

**2. [Calculation] Find the multiplicity of the non-zero eigenvalue $-4\alpha$ for the matrix $Q$ in the JC69 model.**

??? success "Solution"
    **Analysis:**
    1. $Q$ can be written as $4\alpha (\frac{1}{4}J - I)$, where $J$ is the all-ones matrix.
    2. The eigenvalues of $J$ are $\{4, 0, 0, 0\}$.
    3. The eigenvalues of $Q$ are $4\alpha(\frac{1}{4} \cdot 4 - 1) = 0$ and three instances of $4\alpha(\frac{1}{4} \cdot 0 - 1) = -4\alpha$.
    **Conclusion**: The multiplicity is 3.

**3. [Calculation] Find the diagonal entry $P_{AA}(t)$ of the JC69 transition matrix $M(t)$.**

??? success "Solution"
    **Using the Matrix Exponential:**
    1. Diagonalizing $Q$ yields eigenvalues 0 and $-4\alpha$.
    2. The exponential map gives $e^0=1$ and $e^{-4\alpha t}$.
    3. The combination yields $P_{AA}(t) = \frac{1}{4} + \frac{3}{4}e^{-4\alpha t}$.
    **Conclusion**: As $t \to \infty$, the probability of remaining base A approaches $1/4$ (uniform distribution).

**4. [Distance] For two DNA sequences of 100 bp with 10 differences, calculate the JC distance.**

??? success "Solution"
    **Steps:**
    1. Proportion of differences $p = 10/100 = 0.1$.
    2. Formula: $d = -0.75 \ln(1 - 1.333 \cdot 0.1) = -0.75 \ln(0.8667)$.
    3. $\ln(0.8667) \approx -0.143$.
    4. $d \approx -0.75 \cdot (-0.143) \approx 0.107$.
    **Conclusion**: The estimated actual number of mutations is 10.7 per 100 sites (correcting for hidden mutations).

**5. [Kimura] What is the Kimura Two-Parameter (K2P) model, and how does its matrix structure differ?**

??? success "Solution"
    **Difference:**
    K2P distinguishes between **Transitions** (A-G or C-T) and **Transversions**. In the rate matrix, the single rate $\alpha$ is split into $\alpha$ and $\beta$. This reflects the biological reality that substitutions between bases with similar chemical structures are more frequent.

**6. [Property] Prove that the eigenvalues $\lambda$ of a substitution matrix $M(t)$ satisfy $|\lambda| \le 1$.**

??? success "Solution"
    **Reasoning:**
    $M(t)$ is a stochastic matrix (rows sum to 1 and entries are non-negative). By spectral radius properties (Ch17), the maximum eigenvalue of a stochastic matrix is 1, and all others lie within the unit circle.

**7. [Phylogeny] Briefly state the role of linear algebra in building phylogenetic trees.**

??? success "Solution"
    By computing a distance matrix $D$ between species, algorithms like **Neighbor-Joining** use linear transformations to recursively reduce the matrix and group the closest nodes, reconstructing the branching structure of life's evolution.

**8. [Limits] If two sequences are completely random ($p=0.75$), what happens to the JC distance?**

??? success "Solution"
    **Calculation:**
    $1 - \frac{4}{3}(0.75) = 1 - 1 = 0$.
    $\ln(0) \to -\infty$.
    **Significance**: When differences reach 75% (the limit for 4 random bases), all information of a common ancestor is lost. The distance tends to infinity, and algebraic reconstruction becomes impossible.

**9. [PCA] Why is PCA useful for gene expression data?**

??? success "Solution"
    Gene expression matrices have thousands of dimensions. PCA projects this data onto principal components like "tissue type" or "disease state." By finding the eigenvectors of the covariance matrix, linear algebra automatically identifies clusters of core genes responsible for sample variance.

**10. [Application] What is matrix analysis of "Codon Usage Bias"?**

??? success "Solution"
    Different species favor different triplets (codons) to encode the same amino acid. By constructing a $64 \times N$ matrix of codon frequencies and performing Correspondence Analysis (a generalized SVD), researchers can reveal differences in translation efficiency and evolutionary positioning between species.

## Chapter Summary

Linear algebra is the "molecular clock" decoding the history of life:

1.  **Algebraization of Probability**: Substitution matrices transform microscopic random changes into macroscopic matrix exponential evolutions, establishing dynamical standards for evolutionary processes.
2.  **Topological Correction of Distance**: The matrix logarithm formula proves that visible variations are just the tip of the iceberg, providing the only reliable algebraic correction for true genetic distance.
3.  **Dimension Reduction of Information**: From sequence alignment to gene clustering, linear algebra extracts critical structural modes from massive biological datasets, supporting modern computational biology and precision medicine.
