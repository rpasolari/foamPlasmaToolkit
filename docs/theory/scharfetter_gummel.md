# Scharfetter-Gummel Flux Scheme

## 1. Introduction & Motivation

In low-temperature plasma and semiconductor simulations, the transport of charged species (electrons and ions) is governed by the **Drift-Diffusion** approximation. The flux $\mathbf{J}$ of a species with number density $n$ ($1/m^3$), mobility $\mu$, and diffusivity $D$ is given by:

$$
\mathbf{J} = \underbrace{q \mu \mathbf{E} n}_{\text{Drift}} - \underbrace{q D \nabla n}_{\text{Diffusion}}
$$

Where $\mathbf{E} = -\nabla V$ is the electric field calculated from the potential $V$.

### The Numerical Challenge
This equation represents a convection-diffusion problem. A key dimensionless parameter governing stability is the **local Peclet number** ($Pe$), which compares the strength of drift (convection) to diffusion across a grid cell of width $\Delta x$:

$$
Pe = \frac{\text{Drift Velocity} \cdot \Delta x}{\text{Diffusion}} = \frac{\mu E \Delta x}{D}
$$

* **Low Field ($Pe \ll 1$):** Diffusion dominates. Standard **Central Difference** schemes work well.
* **High Field ($Pe \gg 1$):** Drift dominates. Standard Central Difference schemes become numerically unstable, producing **spurious oscillations** (wiggles) and, critically, **negative densities**.

In plasma simulations, negative density is unphysical and causes immediate solver failure. While **Upwind** schemes are stable, they are only first-order accurate and introduce excessive artificial diffusion (smearing).

### The Solution
The **Scharfetter-Gummel (SG)** scheme provides a robust solution. Instead of assuming a linear profile for density (like Central Difference), it assumes an **exponential profile** that matches the physics of the drift-diffusion equation exactly. This guarantees:
1.  **Unconditional Stability:** Stable for any electric field strength.
2.  **Positivity Preservation:** The density $n$ never becomes negative.
3.  **Correct Limiting Behavior:** Automatically transitions between Central Difference (low field) and Upwind (high field).

---

## 2. Derivation: The Analytical 1D Solution

The core idea of the SG scheme is to solve the drift-diffusion equation **analytically** in 1D along the edge connecting two computational nodes.

### Assumptions
Consider two neighboring cell centers, $P$ and $N$, separated by a distance $d$. We make two key assumptions for the interval $[0, d]$:

1.  The **Flux** $J$ is constant along the edge.
2.  The **Electric Field** $E$ is constant along the edge (implying the potential $V$ varies linearly).

### The Differential Equation
In 1D, the steady-state flux equation is given by the sum of drift and diffusion:

$$
J = q \mu n E - q D \frac{dn}{dx}
$$

Recalling that $E = -\frac{dV}{dx}$ and assuming a constant field, the equation is rearranged for the derivative of density:

$$
\frac{dn}{dx} - \frac{\mu}{D} \frac{dV}{dx} n = \frac{J}{qD}
$$

Defining the inverse thermal voltage scale $\lambda = \frac{1}{V_T} \frac{dV}{dx} = \frac{\mu}{D} \frac{dV}{dx}$ (using the Einstein relation), the ordinary differential equation becomes:

$$
\frac{dn}{dx} - \lambda n = \frac{J}{qD}
$$

### Analytical Solution
Multiplying the equation by the integrating factor $e^{-\lambda x}$ yields:

$$
e^{-\lambda x} \frac{dn}{dx} - \lambda e^{-\lambda x} n = \frac{J}{qD} e^{-\lambda x}
$$

Which can be written as a total derivative:

$$
\frac{d}{dx} \left( n e^{-\lambda x} \right) = \frac{J}{qD} e^{-\lambda x}
$$

Integrating over the interval $[0, d]$ from node $P$ to node $N$, and rearranging for the flux J:

$$
J = \frac{q D \lambda}{1 - e^{\lambda d}} \left( n_N - n_P e^{\lambda d} \right)
$$

### The Scharfetter-Gummel Form
It is convenient to express this in terms of the local Peclet number (or normalized potential drop), $Pe = \lambda d$. Substituting $\lambda = Pe/d$ gives:

$$
J = \frac{qD}{d} \frac{Pe}{1 - e^{Pe}} \left( n_N  - n_P e^{Pe} \right)
$$

To ensure numerical stability, this is expressed using the **Bernoulli function**, defined as:

$$
B(x) = \frac{x}{e^x - 1}
$$

Using the $B(-x) = - \frac{x e^x}{1 - e^x}$, the flux equation simplifies to the final Scharfetter-Gummel form:

$$
J = \frac{qD}{d} \left[ B(-Pe) n_P - B(Pe) n_N \right]
$$

This formulation guarantees that the flux is calculated using an exponential profile, providing unconditional stability regardless of the electric field magnitude.

Here is the text for Section 3, detailing the Boundary Condition implementation.

It follows your structure: explaining why standard schemes fail, introducing the surface area normalization, deriving the weights, and proving the validity for fixedValue and fixedGradient.

Markdown

## 3. Boundary Condition Implementation in OpenFOAM

### 3.1 The Implementation Challenge
In standard OpenFOAM schemes (like `fvm::laplacian`), boundary conditions are handled by passing linear coefficients ($\alpha$ and $\beta$) to the discretization operator. The operator assumes the flux is linearly dependent on the gradient.

However, the Scharfetter-Gummel scheme relies on **non-linear weighting** via the Bernoulli function, which depends on the local electric field ($Pe$). Standard operators cannot "see" this field dependence inside the boundary condition. Therefore, we must manually calculate the matrix coefficients to ensure the exponential profile is preserved right up to the wall.

### 3.2 The Discretized Face Flux
To implement this in a Finite Volume framework, we must include the face surface area magnitude $|\mathbf{S}_f|$. We also remove the elementary charge $q$ (assuming the equation is solved for number density flux).

The total particle flux $F_f$ leaving cell $P$ through boundary face $b$ is:

$$
F_f = \frac{D |\mathbf{S}_f|}{d} \left[ B(-Pe) n_P - B(Pe) n_b \right]
$$

We define the **Scharfetter-Gummel weights** for the internal cell ($w_P$) and the boundary face ($w_b$):

$$
w_P = \frac{D |\mathbf{S}_f|}{d} B(-Pe)
$$
$$
w_b = \frac{D |\mathbf{S}_f|}{d} B(Pe)
$$

Thus, the flux simplifies to:

$$
F_f = w_P n_P - w_b n_b
$$

### 3.3 The "Universal" BC Formulation
OpenFOAM boundary conditions (derived from `fvPatchField`) generalize the relationship between the boundary value $n_b$ and the cell center value $n_P$ using two coefficients, `valueInternalCoeffs` ($\alpha$) and `valueBoundaryCoeffs` ($\beta$):

$$
n_b = \alpha n_P + \beta
$$

By substituting this relation into our SG flux equation, we can derive the exact contributions to the linear system matrix ($Ax = b$).

$$
F_f = w_P n_P - w_b (\alpha n_P + \beta)
$$

Grouping the terms by the unknown $n_P$:

$$
F_f = \underbrace{(w_P - w_b \alpha)}_{\text{Diagonal Contribution}} n_P - \underbrace{w_b \beta}_{\text{Source Contribution}}
$$

In OpenFOAM code, these map directly to the matrix coefficients:
* `intCoeffs` (Diagonal): $w_P - w_b \alpha$
* `bdryCoeffs` (Source): $w_b \beta$

---

### 3.4 Verification: Standard Cases

To prove this implementation is correct, we apply it to the two most common boundary conditions.

#### Case A: Fixed Value (Dirichlet)
* **Physics:** The density at the wall is fixed to a specific value $n_{ref}$.
* **OpenFOAM Coeffs:** $\alpha = 0$, $\beta = n_{ref}$.

**Substitution:**
$$
\text{intCoeffs} = w_P - (w_b \cdot 0) = w_P
$$
$$
\text{bdryCoeffs} = w_b \cdot n_{ref}
$$

**Resulting Flux:**
$$
F_f = w_P n_P - w_b n_{ref}
$$
$$
F_f = \frac{D |\mathbf{S}_f|}{d} \left[ B(-Pe) n_P - B(Pe) n_{ref} \right]
$$
> **Verdict:** Matches the analytical SG formula exactly. The diagonal term stabilizes the matrix, and the fixed value acts as a source/sink weighted by the field.

#### Case B: Fixed Gradient (Neumann)
* **Physics:** The gradient is fixed such that $\frac{n_b - n_P}{d} = g$. Rearranging for the boundary value: $n_b = 1 \cdot n_P + g d$.
* **OpenFOAM Coeffs:** $\alpha = 1$, $\beta = g d$.

**Substitution:**
$$
\text{intCoeffs} = w_P - (w_b \cdot 1) = w_P - w_b
$$
$$
\text{bdryCoeffs} = w_b (g d)
$$

**Resulting Flux:**
$$
F_f = (w_P - w_b) n_P - w_b (g d)
$$

Using the Bernoulli identity $B(-x) - B(x) = x$, the diagonal term simplifies to the Drift Flux:
$$
w_P - w_b \propto B(-Pe) - B(Pe) = Pe
$$
$$
F_f = \underbrace{(v_{drift} |\mathbf{S}_f|) n_P}_{\text{Pure Drift Outflow}} - \underbrace{D |\mathbf{S}_f| B(Pe) g}_{\text{Gradient Source}}
$$

> **Verdict:** Correct. The boundary condition correctly calculates that the drift flux pulls carriers out of the cell (diagonal term), while the fixed gradient adds a specific diffusive flux contribution (source term).

### 3.5 The "Zero Gradient" Trap: Diffusion vs. Total Flux

It is critical to distinguish between the mathematical definition of **Gradient** and the physical definition of **Flux**. In standard CFD, `zeroGradient` often implies a "wall" or "symmetry" plane where no flow occurs. **In Drift-Diffusion, this is incorrect.**

The total flux $\mathbf{J}$ consists of two components:

$$
\mathbf{J} = \underbrace{q \mu \mathbf{E} n}_{\text{Drift}} - \underbrace{q D \nabla n}_{\text{Diffusion}}
$$

The standard OpenFOAM boundary conditions `fixedGradient` and `zeroGradient` operate **only** on the $\nabla n$ term (Diffusion). They do not control the Electric Field ($\mathbf{E}$).

#### The `zeroGradient` Case
If you apply `zeroGradient`, you enforce $\nabla n = 0$. The diffusion term vanishes, but the drift term remains:

$$
\mathbf{J}_{boundary} = q \mu \mathbf{E} n + 0
$$

**Physical Consequence:** If there is a non-zero Electric Field at the boundary, the standard `zeroGradient` condition acts as an **Open Boundary** (Outflow/Inflow) where carriers are dragged through the wall by the electric field. It is **not** an insulator.

#### The `fixedGradient` Case
If you apply `fixedGradient` (e.g., $g$), you enforce $\nabla n = g$.

$$
\mathbf{J}_{boundary} = q \mu \mathbf{E} n - q D g
$$

**Physical Consequence:** This adds a constant **Diffusive Source**, but the total current entering the domain is still modulated by the local density $n$ and the electric field. This is physically valid for current injection, but complex to control if the drift velocity is unknown.

## 4. References

1.  **The Original Derivation**
    * Scharfetter, D. L., & Gummel, H. K. (1969). "Large-signal analysis of a silicon Read diode oscillator." *IEEE Transactions on Electron Devices*, 16(1), 64-77.
    * *Note:* The foundational paper introducing the scheme for diode simulation.
