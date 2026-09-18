# Symmetric Relaxation: From Mathematical to Energy Version — Complete Logic and Formula Specification (English Translation)

> The energy-version symmetric relaxation algorithm was produced by **Kai Xu** with the assistance of the **Hy4 large language model (LLM)**, converting the mathematical-version symmetric relaxation algorithm (Xu et al. 2026, arXiv:2606.18604) into an energy version.
> This document fixes only the algorithm logic, formulas, and parameter definitions. It contains no validation, testing, or numerical-experiment content.
> Implementation is entered only after the logic is confirmed.

***

## 0. Positioning and Overview

**Mathematical version**: a set of geometric rules — fit $n_c$ optimal rays → take the triangle centroid as the target point → move; boundary points slide along the boundary; use "cancel the move if it increases $\sum\alpha^2$" as a hard angle gate.

**Energy version**: find a scalar energy $E$ that makes the above rules a special case of "gradient descent on $E$ (over-damped motion)". The original algorithm is upgraded from a "rule set" to a "numerical realization of a variational principle".

One-line summary:
$E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}},\qquad \mathbf F_v=-\nabla_{\mathbf r_v}E,\qquad \Delta\mathbf r_v=\mu\,\mathbf F_v$

***

## 1. Topology and Notation Conventions

| Symbol                     | Meaning                                                                                                                                                                                                                      |
| -------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Interior vertex            | shared by 3 cells (degree 3); its 3 incident edges connect 3 neighbors                                                                                                                                                       |
| Boundary (marginal) vertex | shared by 2 cells (degree 2), lying on the free boundary and sliding along the boundary edge. It has **3 incident edges**: 2 **boundary edges** (each belonging to only 1 cell) + 1 **internal edge** (belonging to 2 cells) |
| Interior cell              | all vertices are interior vertices                                                                                                                                                                                           |
| Boundary (marginal) cell   | contains 2 boundary vertices + on average about 3 interior vertices                                                                                                                                                          |

- The vertices of cell $c$ are ordered **counter-clockwise (CCW)** as $\mathbf r_{c,0},\dots,\mathbf r_{c,n_c-1}$, with $n_c$ the vertex count.
- **Each cell holds its own copy of the vertex coordinates**; the same physical vertex coincides across cells and must be clustered by key into unique vertices.
- $J$ denotes the **+90° rotation** operator: $J(x,y)=(-y,x)$.
- Angles are taken as **lifted / unwrapped** angles, i.e. continuously unfolded along the cyclic order, without modulo $2\pi$.

***

## 2. Energy Definition

### 2.1 Total Energy

$E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$

- **Dimensionless** (depends only on angles, see §7), hence carries no physical force/energy dimensions.
- It does not contain line tension $\Lambda\sum l_e$ (not part of symmetric relaxation; mentioned only as an extension in the Discussion).

### 2.2 Center-Angle Term $E_{\mathrm{cen}}$ (Cell Level)

For cell $c$:

- Centroid (**vertex centroid**, see §8.1):
  $\mathbf O_c=\frac{1}{n_c}\sum_{k=0}^{n_c-1}\mathbf r_{c,k}$
- Radial vector and direction angle of vertex $i$:
  $\mathbf u_{c,i}=\mathbf r_{c,i}-\mathbf O_c,\qquad \rho_{c,i}=|\mathbf u_{c,i}|,\qquad \hat{\mathbf u}_{c,i}=\frac{\mathbf u_{c,i}}{\rho_{c,i}},\qquad \psi^{\uparrow}_{c,i}=\arg(\mathbf u_{c,i})$

> where $\psi^{\uparrow}_{c,i}$ is the **lifted angle** (continuously unfolded along the cyclic order, no modulo $2\pi$, see §1).

- **Best-fit phase** (arithmetic mean, closed form):
  $\varphi_c=\frac{1}{n_c}\sum_{k=0}^{n_c-1}\left(\psi^{\uparrow}_{c,k}-\frac{2\pi k}{n_c}\right)$
- **Center-angle deviation**:
  $\boxed{\ \delta_{c,i}=\psi^{\uparrow}_{c,i}-\varphi_c-\frac{2\pi i}{n_c}\ }$
- **Center-angle energy** (**bare sum**, without $\kappa_c$):
  $\boxed{\ E_{\mathrm{cen}}=\sum_c\sum_{i=0}^{n_c-1}\delta_{c,i}^2\ }$

> **Important convention**: $E_{\mathrm{cen}}$ and $E_{\mathrm{ang}}$ are always defined as **bare sums without stiffness**; stiffness appears only in the total energy $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$ (§2.1). Only then does it stay consistent with the gradient weights $2\kappa_c\delta$, $2\kappa_a(\alpha-\tau)$ of §3. If $\kappa$ were written into the definitions of $E_{\mathrm{cen}}$/$E_{\mathrm{ang}}$, $E$ would contain $\kappa^2$ and all of §3's formulas would break.

**Identity**: $\displaystyle\sum_{i}\delta_{c,i}\equiv 0$ (follows directly from $\varphi_c$ being the arithmetic mean). This identity is the key to the later simplicity of the gradient and the pure-azimuthal nature of the force.

**Meaning**: $\delta_{c,i}$ measures how much vertex $i$'s actual direction deviates from "the best arrangement of $n_c$ equally-spaced directions". $E_{\mathrm{cen}}=0\iff$ every adjacent pair of vertex directions spans an angle $2\pi/n_c$.

#### Relationship between "Best-Fit Phase" and "Optimal Rays"

The full meaning of the mathematical version's "find optimal rays" is: from $\mathbf O_c$, draw $n_c$ **rigid** rays (fixed adjacent angle $2\pi/n_c$); the whole fan is allowed only an overall rotation by an angle $\varphi$ so as to best fit the actual $n_c$ vertex directions. Minimizing $\sum_i(\psi^{\uparrow}_{c,i}-\varphi-2\pi i/n_c)^2$ gives exactly the $\varphi_c$ above.

> **Rigidity is the essence of the center-angle term**: if each ray could rotate on its own to aim at its vertex, the residual would be identically zero and the very concept of "center angle" would vanish. The only degree of freedom is the overall phase.

Therefore "computing $\varphi_c$" and "finding optimal rays" **are the solution of the same least-squares problem**, but only the first of its three outputs is retained:

| Layer                                                                  | What the mathematical version uses it for | Energy version                                   |
| ---------------------------------------------------------------------- | ----------------------------------------- | ------------------------------------------------ |
| ① Phase $\varphi_c$                                                    | Sets the ray orientation                  | **Kept**: as a reference for measuring deviation |
| ② Full geometry of the rays ($n_c$ half-lines)                         | Drawn out                                 | **Discarded**: no rays are constructed           |
| ③ Rays intersect into a triangle; its centroid is the **target point** | Vertices move directly there              | **Discarded**: no target point is generated      |

That is: **the energy version inherits the solution but discards the geometric construction** — it only asks "how far is the actual direction off", not "where is the ideal position".

#### Algebraic Equivalence with the Mathematical Version's Least Squares

The mathematical version's ray fit is usually written in ascending-angle form: rank the vertex directions $\theta_1<\theta_2<\cdots<\theta_{n_c}$, set the optimal-ray directions to $\delta_i=\delta_1+2\pi(i-1)/n_c$, and choose

$$\delta_1=\frac{1}{n_c}\Bigl(\sum_{j=1}^{n_c}\theta_j-\pi(n_c-1)\Bigr)$$

to minimize $\sum_{j=1}^{n_c}(\theta_j-\delta_j)^2$. Differentiating the objective w.r.t. $\delta_1$ and setting to zero gives $\delta_1=\frac{1}{n_c}\sum_j\bigl(\theta_j-2\pi(j-1)/n_c\bigr)$, which is **algebraically identical** to the code's $\varphi_c=\frac{1}{n_c}\sum_k(\psi^{\uparrow}_{c,k}-2\pi k/n_c)$ (the index shift $i=k$ vs. $j-1$ is immaterial once the vertices are ordered cyclically and lifted, see §1).

Consequently, the quantity the mathematical version minimizes **once** — $\sum_{j=1}^{n_c}(\theta_j-\delta_j)^2$ — is **exactly** the single-cell center-angle energy $E_{\mathrm{cen},c}=\sum_i\delta_{c,i}^2$. In other words: **$E_{\mathrm{cen}}$ *is* the mathematical version's ray-fit residual, evaluated every iteration.** The mathematical version uses that residual only to orient the rays and then removes it by construction (the triangle-target step); the energy version promotes the residual itself to the objective and reduces it by gradient descent.

**Per (cell, vertex), the deviation is measured against *that cell's own* single ray.** An interior vertex touches three rays (one from each of its three incident cells) and a marginal vertex touches two; but each ray serves only as the angular reference of *its* cell's term, and the rays are **never combined into a triangle or a target point** — that assembly is precisely the geometric construction discarded above.

#### The Role of $\varphi_c$: Gauge Fixing, Not a Target Configuration

$E_{\mathrm{cen}}$ measures "how far the vertex directions deviate from a rigid equiangular fan", but the fan's **overall orientation is irrelevant** — rotating the whole cell should not change its degree of symmetry. The role of $\varphi_c$ is precisely to absorb this redundant orientational degree of freedom.

- If one forced $\varphi_c\equiv0$, a **perfectly regular** hexagon rotated by only $30^\circ$ would yield a huge $E_{\mathrm{cen}}$, purely because the template did not rotate along.
- Thus $\varphi_c$ is not an optional decoration but the **guarantee of the energy's rotation invariance**: $E_{\mathrm{cen}}$ is a measure of symmetry only after the overall rotation is modded out.
- Contrast: **the mathematical version treats** **$\varphi_c$** **as a "target" (what to become); the energy version treats** **$\varphi_c$** **as a "ruler" (the reference for measuring deviation).** Same number, completely different semantics.

#### Constrains Direction Only, Not Radius

- A ray is a **half-line** that gives only a direction, not how far along it a vertex stands. The mathematical version must therefore additionally use "rays intersect into a triangle, take its centroid" to supply a definite **position** — that step fixes both direction and radius at once.
- The energy version's $E_{\mathrm{cen}}$ depends only on the direction angle $\psi^{\uparrow}$ and **not at all on** **$\rho$**. In the force (§3.1) this shows up as the center-angle force always pointing along $J\hat{\mathbf u}$, perpendicular to the radial direction → **pure azimuthal force, doing no radial work**.
- Corollary: **the cell size (radius) is not determined by** **$E_{\mathrm{cen}}$** but indirectly by $E_{\mathrm{ang}}$ and initial conditions. This is two sides of the same fact as §7.3's "uniform expansion is a zero mode, absolute size is not determined by the energy".

### 2.3 Interior-Angle Term $E_{\mathrm{ang}}$ (Vertex Level)

For every **interior angle** (each cell has one interior angle at each of its vertices):

- The interior angle at vertex $w$ is spanned by two edges: $\mathbf a=\mathbf r_p-\mathbf r_w$ (pointing to the CCW predecessor $p$), $\mathbf b=\mathbf r_q-\mathbf r_w$ (pointing to the CCW successor $q$).
- Interior angle
  $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$
  where the lifted angles satisfy $\psi^{\uparrow}_a\in\bigl(\psi^{\uparrow}_b,\ \psi^{\uparrow}_b+2\pi\bigr)$. **The direction must not be reversed**: it is the angle swept counterclockwise from "the outgoing-edge direction $\mathbf b$" to "the reversed incoming edge $\mathbf a$" — i.e. the angle on the **inner** side of the CCW polygon.
  Equivalent form: $\alpha=\pi-(\theta_{\mathrm{out}}-\theta_{\mathrm{in}})$ ($\theta$ being the lifted azimuth of the CCW directed edge), with $\alpha\in(0,\pi)$ for convex vertices.

> **Common mistake**: writing $\alpha=\psi^{\uparrow}_b-\psi^{\uparrow}_a$. This gives $-\alpha$ and flips the sign of all gradients in §3.2 (and the edge form of §3.2.2 with it).

**Target angle** $\tau$ depends on the vertex degree:

| Vertex type     | Degree | $\tau$             | Source                                                                         |
| --------------- | ------ | ------------------ | ------------------------------------------------------------------------------ |
| Interior vertex | 3      | $2\pi/3=120^\circ$ | foam three-film force balance (Plateau), **derived**                           |
| Boundary vertex | 2      | $\pi/2=90^\circ$   | free-boundary orthogonality condition (boundary analogue of 120°), **derived** |

$\boxed{\ E_{\mathrm{ang}}=\sum_{\text{corner}}(\alpha-\tau)^2\ }$

> The sum runs over **every interior angle** (one per cell-vertex). Again a **bare sum**; stiffness $\kappa_a$ appears only in the total energy.

**Key equivalence with the mathematical version** (interior vertex): the three interior angles sum to $2\pi$, hence
$\sum_{k=1}^{3}(\alpha_k-\tfrac{2\pi}{3})^2=\sum_k\alpha_k^2-\frac{4\pi}{3}\underbrace{\sum_k\alpha_k}_{=2\pi}+3\left(\frac{2\pi}{3}\right)^2=\sum_k\alpha_k^2+\text{const}$
i.e. **minimizing** **$E_{\mathrm{ang}}$** **≡ minimizing** **$\sum\alpha^2$**, and $\sum\alpha^2$ is exactly the quantity guarded by the mathematical version's hard angle constraint.

**Decomposition at a boundary vertex**: when the two boundary interior angles $\beta_1,\beta_2$ take $90^\circ$,
$(\beta_1-\tfrac{\pi}{2})^2+(\beta_2-\tfrac{\pi}{2})^2=\tfrac12(\beta_1-\beta_2)^2+\tfrac12(\beta_1+\beta_2-\pi)^2$

- $\tfrac12(\beta_1-\beta_2)^2$: the two films are symmetric about the boundary (boundary-orthogonality term).
- $\tfrac12(\beta_1+\beta_2-\pi)^2$: the discrete turning angle at the boundary (free-boundary bending term).
- The third angle (between the two internal edges $=2\pi-(\beta_1+\beta_2)$) has **no penalty term** and is left to geometry.

***

## 3. Complete Force Formulas

Vertex force = the negative gradient of the total energy w\.r.t. position:
$\boxed{\ \mathbf F_v=-\nabla_{\mathbf r_v}E=\underbrace{\left(-\kappa_c\sum_{c\ni v}\nabla_{\mathbf r_v}E_{{\mathrm{cen}},c}\right)}_{\text{centre-angle force}}+\underbrace{\left(-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}\right)}_{\text{interior-angle force}}\ }$

> The two terms are the **center-angle force** (§3.1) and the **interior-angle force** (§3.2).

### 3.1 Center-Angle Force (**Core formula — mind the parentheses**)

Define cell $c$'s **cell-wide shared average term** (identical for all vertices, no index $i$):
$\boxed{\ M_c=\frac{1}{n_c}\sum_{j=0}^{n_c-1}\delta_{c,j}\frac{J\hat{\mathbf u}_{c,j}}{\rho_{c,j}}\ }$

Then the center-angle force on vertex $i$ in cell $c$ is:
$\boxed{\ \mathbf F^{(c)}_{v}=-2\kappa_c\left(\delta_{c,i}\frac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)\ }$

- **Direct term** $\delta_{c,i}J\hat{\mathbf u}_{c,i}/\rho_{c,i}$: depends only on the vertex itself; direction perpendicular to radial (azimuthal force), magnitude $\propto$ the vertex's angular deviation, $\propto 1/\rho$.
- **Shared term** $M_c$: identical for **every vertex** of cell $c$ → pulls all the cell's vertices together in the same direction. This is the variational expression of **"the center angle binds all the cell's vertices together"**.

**Vertex membership count** (how many force contributions to sum):

- Interior vertex → sum of center-angle forces from its **3** cells;
- Boundary vertex → sum of center-angle forces from its **2** cells (each cell still binds its 2 boundary + \~3 interior vertices).

**Force zero-sum property** (direct consequence of the vertex centroid + $\sum\delta\equiv0$):
$\sum_{i=0}^{n_c-1}\mathbf F^{(c)}_{v_i}=-2\kappa_c\left(\sum_i\delta_{c,i}\frac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-n_cM_c\right)=\mathbf 0$
i.e. the center-angle force is a **purely azimuthal force**: it sums to zero over a cell, doing **no work under uniform expansion**.

### 3.2 Interior-Angle Force (Including Endpoint Contributions — Easy to Miss)

For a single interior angle $\alpha$ (at vertex $w$, arms pointing to $p,q$, $\mathbf a=\mathbf r_p-\mathbf r_w$, $\mathbf b=\mathbf r_q-\mathbf r_w$), write the weight
$w_\alpha=2\kappa_a(\alpha-\tau)$

The **gradient contribution** of this interior angle to the three relevant vertices:

**Derivation** ($J$ is +90° rotation): $\nabla_{\mathbf r_w}\arg(\mathbf r_p-\mathbf r_w)=-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$, $\nabla_{\mathbf r_p}\arg(\mathbf r_p-\mathbf r_w)=+\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$, and analogously for $\mathbf b$; substituting $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$ (§2.3) gives the table below.

| Vertex | Role                               | Gradient contribution $\nabla(\alpha)$                                                             | Contribution to $\kappa_a\nabla E_{\mathrm{ang}}$                                                                       |
| ------ | ---------------------------------- | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| $w$    | apex                               | $-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}+\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$ | $w_\alpha\left(-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}+\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}\right)$ |
| $p$    | forearm endpoint (CCW predecessor) | $+\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$                                                 | $+w_\alpha\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$                                                              |
| $q$    | rear-arm endpoint (CCW successor)  | $-\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$                                                 | $-w_\alpha\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$                                                              |

**Force = the negative of the above contributions**, i.e. $\mathbf F^{\mathrm{ang}}=-\kappa_a\nabla E_{\mathrm{ang}}$.

- The three terms sum to zero → the interior-angle force is an **internal force**, producing no net translation.
- **Key point**: moving a vertex $v$ changes not only the angle "above its own head" but also the angles in which it participates as an **arm endpoint** at its 3 neighbors. Hence when assembling $\mathbf F_v$ one must **accumulate both** the "angle where $v$ is the apex" and the "angle where $v$ is an arm endpoint" contributions, otherwise the gradient is incomplete.

**Number of interior angles per vertex (fixed, independent of network size)**:

- Interior vertex (degree 3, 3 **internal edges**): 3 as apex (one per cell) + 6 as arm endpoint (3 edges, each belonging to 2 cells) → **9 terms total**.
- Boundary vertex (degree 2, 2 **boundary edges** + 1 **internal edge**, see §1): 2 as apex + 4 as arm endpoint (internal edge belongs to 2 cells → 2 terms; the 2 boundary edges each belong to 1 cell → 1 term each) → **6 terms total**.

### 3.2.1 Scope: One Hop Only, **No Cascading Propagation**

$\nabla_{\mathbf r_v}E_{\mathrm{ang}}$ is defined by "treat $\mathbf r_v$ as the variable, all else as constant, and differentiate" — therefore it **accumulates only terms that explicitly contain** **$\mathbf r_v$**. A single interior angle depends on only 3 positions (apex + two arm endpoints), so the angles involving $\mathbf r_v$ are **only the fixed few terms in the one-hop neighborhood listed above**.

- "Moving $v$ changes neighbors' angles" is an **outcome**, not a **differentiation term**; it is already recorded, to first order and exactly, by these 9 (or 6) terms.
- Farther indirect effects (neighbors of neighbors …) do not appear in the expression for $\partial E/\partial\mathbf r_v$; they are the responsibility of the **next iteration** — that is the job of iteration, not of a single-step gradient.
- Hence **no cascade, no** **$O(V^2)$** **or exponential blow-up**: $O(1)$ terms per vertex, $O(V)$ globally.

### 3.2.2 Edge (Pairwise) Equivalent Form — **Recommended Implementation**

After merging the three rows of the table above, the "apex contribution" and "arm-endpoint contribution" on the same edge exactly **subtract**. Define the scalar weight of each **corner** (corner = a (vertex, cell) pair)
$w_{v,c}=2\kappa_a(\alpha_{v,c}-\tau_v)$

and denote the **interior-angle gradient with stiffness included** by
$\mathbf G_v\equiv\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}\qquad(\mathbf F^{\mathrm{ang}}_v=-\mathbf G_v)$

**Traversal scope**: iterate cell by cell along **CCW** over **all** edges of the cell —

- the **2 boundary edges** of marginal cells participate too (angles at boundary vertices are weighted as usual);
- an **internal edge is traversed once by each of its two cells**; the two traversals have different weights ($\alpha$ belongs to different cells), **both must be accumulated**, no de-duplication.

For each directed edge $X\to Y$ (length $l$, unit vector $\hat{\mathbf e}_{X\to Y}$), set
$\boxed{\ \mathbf g_c(X\!\to\!Y)=\bigl(w_{X,c}-w_{Y,c}\bigr)\frac{J\hat{\mathbf e}_{X\to Y}}{l}\ }$

then
$\mathbf G_X\;+{=}\;\mathbf g_c(X\!\to\!Y)$
$\mathbf G_Y\;-{=}\;\mathbf g_c(X\!\to\!Y)$

**Equivalence derivation** (assume CCW order $\dots\to X\to Y\to\dots$, the edge is $\mathbf e$; note $\mathbf b_X=\mathbf e$, $\mathbf a_Y=-\mathbf e$):

- **Vertex** **$X$**: ① as apex of $\alpha_X$, take the $\mathbf b$ part of the apex row → $+w_X\dfrac{J\hat{\mathbf e}}{l}$; ② as the forearm endpoint of $\alpha_Y$ ($X$ is the CCW predecessor of $Y$, $\mathbf a_Y=-\mathbf e$) → $+w_Y\dfrac{J(-\hat{\mathbf e})}{l}$. Total $(w_X-w_Y)\dfrac{J\hat{\mathbf e}}{l}=\mathbf g$.
- **Vertex** **$Y$**: ① as apex of $\alpha_Y$, take the $\mathbf a$ part of the apex row → $-w_Y\dfrac{J(-\hat{\mathbf e})}{l}=+w_Y\dfrac{J\hat{\mathbf e}}{l}$; ② as the rear-arm endpoint of $\alpha_X$ ($\mathbf b_X=\mathbf e$) → $-w_X\dfrac{J\hat{\mathbf e}}{l}$. Total $(w_Y-w_X)\dfrac{J\hat{\mathbf e}}{l}=-\mathbf g$.

Thus the contributions of all three roles (apex / forearm / rear-arm) in the §3.2 table are **all included, with no omission and no duplication**.

> **Sign self-check**: when $w_X=w_Y$ (the two angles have equal deviation) the edge contribution is zero — as intuition demands: two equally-weighted adjacent angles, pushing the common endpoint along the edge normal produces no first-order energy change.

**Advantages**:

1. **Two lines of code, one pass**; no "forearm / rear-arm" branch logic, less error-prone than accumulating apex-by-apex.
2. **Force zero-sum holds automatically**: each edge contributes a pair $+\mathbf g$ and $-\mathbf g$ → $\sum_v\mathbf F_v=\mathbf 0$ holds explicitly, no spurious bulk translation (§7.3).
3. **Same-order complexity**: the traversal volume is all (cell, edge) incidences $=\sum_c n_c=\sum_v\deg(v)\approx 3V$ ($\deg$ per §1 is the number of cells the vertex belongs to), **same order** as "apex-only" — the scatter writes 2 vectors per step instead of 1, a constant factor of about 2. The real cost driver is computing $\alpha$ itself (\~$3V$ inverse trig calls), which both versions pay, so the difference between scatter and accumulation is largely masked.

### 3.2.3 The Arm-Endpoint Term **Must Not Be Omitted** (Not a Precision Trade-off)

If only the apex term is kept, one obtains a **different operator**, with clear costs:

1. **Variational principle breaks**: $\mathbf F_v\neq-\nabla_{\mathbf r_v}E$ → what the iteration minimizes is **not** the written $E$, and the core claim "the original algorithm is a special case of gradient descent on $E$" fails.
2. **Spurious bulk drift**: $E_{\mathrm{ang}}$ depends only on relative positions; the full gradient strictly satisfies $\sum_v\mathbf F_v=\mathbf 0$; after truncation this no longer holds, the tissue drifts/rotates from nothing, and this motion is misread as relaxation progress.
3. **Line search stalls**: the truncated direction is not a descent direction → energy does not drop → step is forced to halve down to zero → stagnation. Trying to save effort makes it less stable.

**It cannot be dropped on magnitude grounds either**: near equilibrium $w\propto(\alpha-\tau)$ is a **first-order** small quantity, the same order as the apex term — not a negligible higher-order correction.

**If one insists on simplifying**, the only option is to degrade to a **block-local rule** (each vertex moves only by the angle above its own head, a Jacobi-style local projection). It runs, but loses $E$ as a Lyapunov function, loses the monotone-convergence guarantee, and cannot be called a variational model. **This spec does not adopt it.**

### 3.3 Assembly

$\mathbf F_v=\mathbf F_v^{\mathrm{cen}}+\mathbf F_v^{\mathrm{ang}}$

- $\displaystyle\mathbf F_v^{\mathrm{cen}}=\sum_{c\ni v}\mathbf F^{(c)}_{v}$ — **center-angle force** (§3.1, $\mathbf F^{(c)}$ already includes $\kappa_c$).
- $\displaystyle\mathbf F_v^{\mathrm{ang}}=-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}=-\mathbf G_v$ — **interior-angle force** (§3.2, $E_{\mathrm{ang}}$ is a bare sum, stiffness $\kappa_a$ is multiplied in here; $w_{v,c}$ already contains $\kappa_a$, so the §3.2.2 accumulation yields $\mathbf G_v$ rather than $\nabla E_{\mathrm{ang}}$).

**The interior-angle force is assembled in edge form**: first compute $\alpha_{v,c}$ and the scalar weight $w_{v,c}$ per corner, then traverse each cell's boundary CCW and scatter $\mathbf g_c(X\!\to\!Y)$ pairwise to the two endpoints per §3.2.2. This disposes of all three contributions (apex + arm endpoints) in one pass, with no need to enumerate roles separately.

### 3.4 Vertex-Level Aggregation: From Scalar Energy to Resultant Force

**Conceptual bridge.** §2 defines two *scalar* energy terms; §3.1–§3.2 give the *force* of each. The present section explains how, at a single vertex $v$, the several scalar deviations it carries are turned into one resultant force direction — and why this is *not* the geometric version's triangle-centroid step.

**A scalar deviation has no direction.** Whether $\delta_{c,i}^{2}$ (center-angle) or $(\alpha-\tau)^{2}$ (interior-angle), the deviation is a single number. The force direction does not come from the number itself; it comes from the *spatial gradient* of that number with respect to the vertex position:

$$\mathbf F_{v}=-\nabla_{\mathbf r_v}\big(\text{scalar terms involving }\mathbf r_v\big)
=-\frac{\partial(\text{deviation}^{2})}{\partial(\text{deviation})}\cdot\frac{\partial(\text{deviation})}{\partial\mathbf r_v}
=-2(\text{deviation})\cdot\Big(\frac{\partial\text{deviation}}{\partial\mathbf r_v}\Big).$$

The first factor (magnitude & sign) is the scalar; the second factor (a geometric vector) is the direction.

**Center-angle contribution — tangential, magnitude $\propto|\delta|/\rho$.** For a cell $c$ with local index $i$,

$$\psi=\operatorname{azimuth}(\mathbf r_v-\mathbf O_c),\qquad
\hat{\mathbf u}=\frac{\mathbf r_v-\mathbf O_c}{\rho},\qquad \rho=|\mathbf r_v-\mathbf O_c|.$$

The gradient of the azimuth is a *tangential* vector:

$$\nabla_{\mathbf r_v}\psi=\frac{1}{\rho}\,J\hat{\mathbf u},\qquad J=\text{+90° rotation}.$$

Hence the per-cell center-angle force on $v$ (ignoring the cell-shared term, §3.1) is

$$\mathbf F^{(c)}_{v}=-2\kappa_c\,\delta_{c,i}\,\frac{1}{\rho}\,J\hat{\mathbf u}\;-\;\mathbf F^{\mathrm{shared}}_c.$$

Three consequences, all visible from this formula:
- **Direction is tangential** ($\perp$ the centroid line), never radial toward $\mathbf O_c$.
- **Magnitude $\propto|\delta_{c,i}|$**: larger deviation $\rightarrow$ stronger force; sign flip $\rightarrow$ reversed force.
- **Magnitude $\propto 1/\rho$**: farther from the centroid, the azimuth is less sensitive to translation, so the force is weaker.
- **Radial force = 0**: sliding $v$ along the centroid line does not change $\psi$, so $E_{\mathrm{cen}}$ does not drive cell *size* — it is a pure shape force.

(The cell-shared term $-2\kappa_c\mathbf{shared}_c$ is a cell-constant, identical for all vertices of the cell; it enforces the translation zero-mode and cell rigidity, but does not change the per-force *tangential* character above.)

**Interior-angle contribution — along the angle-opening direction.** For each corner (one per cell-vertex) at $v$ with arms $\mathbf a,\mathbf b$ to the previous/next neighbours,

$$\mathbf g_v=\frac{J\hat{\mathbf a}}{|\mathbf a|}-\frac{J\hat{\mathbf b}}{|\mathbf b|},
\qquad \mathbf F^{\mathrm{corner}}_v=-2\kappa_a\,(\alpha-\tau_v)\,\mathbf g_v.$$

The direction $\mathbf g_v$ points along the angle bisector (opening/closing the angle); the magnitude $\propto|\alpha-\tau_v|$. The same corner also couples to its two neighbouring vertices (§3.2.2), so each interior-angle contribution is scattered to three vertices.

**Resultant = vector sum, no triangle, no target point.** For a degree-3 interior vertex shared by cells $c_1,c_2,c_3$ with local indices $i_1,i_2,i_3$ and three interior angles $\alpha_1,\alpha_2,\alpha_3$,

$$\mathbf F_v=-\sum_{k=1}^{3}\Big[\,2\kappa_c\big(\delta_{c_k,i_k}\tfrac{1}{\rho_k}J\hat{\mathbf u}_k-\mathbf{shared}_{c_k}\big)\Big]
\;-\;\sum_{k=1}^{3}\Big[\,2\kappa_a(\alpha_k-\tfrac{2\pi}{3})\,\mathbf g^{(k)}_v\Big]
\;+\;\text{(neighbour-coupling terms)}.$$

All center-angle and interior-angle gradient vectors are accumulated onto the shared vertex by `np.add.at` (energy.py). **The three rays / three deviations are never assembled into a triangle or an absolute target position** — that geometric construction is precisely what the energy version discards (§2.2). Instead, the three scalar energies are differentiated and their gradient vectors are added.

**Contrast with the mathematical version.** The mathematical version fits the $n_c$ optimal rays, then for each vertex takes the *three* rays from its three cells, forms the *triangle* of three rays, and moves the vertex to the triangle's centroid — a single absolute target point. The energy version replaces this with: each ray yields one scalar residual $\delta^{2}$, each forms one gradient vector, and the three vectors are summed into a resultant force; the vertex moves along $-\mathbf F_v$ (then line-searched). Because the two versions optimize *different objectives* (an explicit target point vs. the net-force-zero of a single scalar $E$), their optima do **not** coincide exactly — which is exactly why the measured $\lambda=1$ result only *approximately* reproduces the mathematical version.

***

## 4. Dynamics and Motion

### 4.1 Over-Damped Motion

$\zeta\frac{d\mathbf r_v}{dt}=\mathbf F_v\quad\Longrightarrow\quad \mathbf v_v=\frac{1}{\zeta}\mathbf F_v$

Forward Euler discretization (time step $\Delta t$):
$\boxed{\ \mathbf r_v\leftarrow\mathbf r_v+\text{step}\cdot\mathbf F_v\ },\qquad \text{step}=\mu L_0^{2}\;(\text{L0 normalization, default}),\quad \mu\equiv\frac{\Delta t}{\zeta}$

- $\mu$ = **mobility** (step-size constant), the **only** time-scale parameter that needs to be given.
- The step uses **L0 normalization**: let $L_0$ be the median edge length of the tissue; then $\text{step}=\mu L_0^{2}$ and the per-round displacement $\approx\mu L_0$ — $\mu$ becomes a scale-invariant dimensionless step, needing no adjustment across networks of different sizes (on by default).

### 4.2 Direction of Motion

- **Interior vertex**: free 2D, along $\mathbf F_v$, $\Delta\mathbf r_v=\mu\,\mathbf F_v$.
- **Boundary vertex**: can only slide along the boundary. It has **2 boundary edges** (§1); their tangents differ and form a kink at the vertex, so **the tangent is not uniquely determined** and must be chosen by "steepest descent":
  $\hat{\mathbf t}_v=\arg\max_{\hat{\mathbf t}\in\{\hat{\mathbf t}_1,\ \hat{\mathbf t}_2\}}\ \bigl\lvert\mathbf F_v\cdot\hat{\mathbf t}\bigr\rvert,\qquad \Delta\mathbf r_v=\mu\,(\mathbf F_v\cdot\hat{\mathbf t}_v)\,\hat{\mathbf t}_v$
  where $\hat{\mathbf t}_1,\hat{\mathbf t}_2$ are the unit directions of the vertex's two boundary edges (pointing to their respective outer neighbors).
  - Meaning: project the total force onto the boundary edge **most parallel** to it and slide along that edge;
  - if both edge projections are near zero (force nearly perpendicular to the boundary), the vertex **does not move** this round;
  - this is a **dynamic constraint** (restricting the displacement direction), not an energy term — the energy is computed as usual per §2, and only the force is projected when applied.

### 4.3 Step-Size Adaptation (Line Search, Serial Per-Vertex)

The **serial descending-order** update (§5) is used: free vertices are sorted by $|\mathbf F_v|$ in **descending** order and moved one by one (Gauss–Seidel, later vertices use already-updated coordinates). **Each vertex does its own line search**: starting from that vertex, the displacement first takes the $\text{step}$ of §4.1; if it causes

1. total energy to rise ($E_{\mathrm{new}}>E_{\mathrm{old}}$, comparing only the local energy of the cells incident to that vertex), or
2. some **cell incident to that vertex** that "was convex at the start of this round" to become non-convex

then the step is halved (up to 60 times) and the single-vertex step is recomputed until the acceptance condition holds; if all halvings fail, that vertex **does not move** this round.

The convexity guard is **localized to the vertex being moved** — a fragile cell constrains only the vertex it belongs to, and unlike the synchronous scheme it does not veto the whole global step because one cell flips non-convex (under the synchronous scheme a globally rejected step freezes large random networks). Cells that are already non-convex at the start are **exempt** (they are not required to become convex; the angle energy will naturally pull them back).

### 4.4 Velocity Decay

The force magnitude naturally shrinks as equilibrium is approached ($\mathbf F\propto$ deviation) → displacement automatically shrinks → over-damped relaxation, no extra cooling schedule needed.

***

## 5. Algorithm Flow (Single Round)

1. **Build vertex index**: cluster the duplicate coordinates held by cells into unique vertices by key; record each vertex's degree (number of cells sharing it) and its references in each cell; for boundary vertices **record the directions of their 2 boundary edges** $\hat{\mathbf t}_1,\hat{\mathbf t}_2$ (for the §4.2 tangent choice).
2. **Compute center-angle quantities per cell**: $\mathbf O_c$ (vertex centroid) → lifted direction angle $\psi^{\uparrow}$ → phase $\varphi_c$ → deviation $\delta_{c,i}$; also compute the shared term $M_c$.
3. **Accumulate center-angle force**: per §3.1, directly accumulate each cell's $\mathbf F^{(c)}_v$ (which already includes $\kappa_c$ and is a **force**, not a gradient) onto its vertices.
4. **Compute angle quantities and weights per interior angle**: for each corner (vertex, cell) compute $\alpha_{v,c}=\psi^{\uparrow}_a-\psi^{\uparrow}_b$, take $\tau_v$ by vertex degree, and compute the scalar weight $w_{v,c}=2\kappa_a(\alpha_{v,c}-\tau_v)$. Then **traverse all edges of each cell CCW** (including marginal cells' boundary edges; an internal edge is traversed once by each of its two cells and both are accumulated), and per §3.2.2 scatter $\mathbf g_c(X\!\to\!Y)=(w_{X,c}-w_{Y,c})J\hat{\mathbf e}_{X\to Y}/l$ pairwise to the two endpoints — apex and arm-endpoint terms disposed of in one pass.
5. **Assemble total force** $\mathbf F_v=\mathbf F^{\mathrm{cen}}_v+\mathbf F^{\mathrm{ang}}_v$ (§3.3).
6. **Serial update in descending** **$|\mathbf F_v|$** **order**: sort free vertices by force magnitude $|\mathbf F_v|$ in **descending** order; move each vertex along its own force $\Delta\mathbf r_v=\text{step}\cdot\mathbf F_v$ (boundary vertices take the §4.2 tangent component), using Gauss–Seidel (later vertices computed from already-updated coordinates), with **each vertex doing its own line search** (§4.3).
7. **Refresh geometry**: write the new coordinates back into each cell's coordinate copy.
   > Note: the recomputed **area** is for display/export only; **the energy model does not use area**; the **centroid** is computed directly from vertex coordinates in step 2 of the next round, so no caching is needed here.
8. Repeat 1–7 until the energy change falls below tolerance or the maximum number of rounds is reached.

***

## 6. Correspondence with the Mathematical Version

| Mathematical version (geometric rules)                                     | Energy version (variational interpretation) |
| -------------------------------------------------------------------------- | ------------------------------------------- |
| Least-squares fit of $n_c$ optimal rays with angle $2\pi/n_c$              | The same least-squares problem, but **only its solution $\varphi_c$ is taken** (closed form: arithmetic mean) as the measurement reference; the geometric construction of rays is discarded (§2.2) |
| Three optimal rays form a triangle; its centroid is the target point       | **No target point is generated**; instead a gradient-descent step on $\mathbf r_v$ replaces that closed-form step. Consequence: the radius is no longer directly constrained by $E_{\mathrm{cen}}$ (§2.2) |
| Cancel the move if it increases $\sum\alpha^2$                             | The hard-constraint limit of $E_{\mathrm{ang}}$ ($\lambda\to\infty$) |
| Boundary point moves toward the midpoint along the edge of the small angle | Steepest descent on the **full boundary interior-angle energy** $\tfrac12(\beta_1-\beta_2)^2+\tfrac12(\beta_1+\beta_2-\pi)^2$ on the 1-D boundary manifold, subject to the tangent-projection constraint (§4.2, taking the larger projection among the two boundary edges) |
| Move serially in **descending** annealing-distance order                   | **Serial descending order** (Gauss–Seidel serial update in descending $|\mathbf F_v|$ order): vertices move one by one in descending force magnitude, each doing its own line search with the convexity guard localized to its incident cells. This aligns semantically with the geometric version's "descending annealing-distance" order; unlike the synchronous scheme — where one cell flipping non-convex vetoes the whole step globally and freezes large random networks — the serial scheme does not. |

**Overall relationship**: the mathematical version = alternating minimization of $E_{\mathrm{cen}}$ (first fix the phase in closed form, then fix the target point analytically) + treating $E_{\mathrm{ang}}$ as a hard constraint; the energy version = **serial-descending (Gauss–Seidel) soft-weighted minimization** of $E_{\mathrm{cen}}+\lambda E_{\mathrm{ang}}$. **The two belong to the same family; no claim of full equivalence is made.**

### 6.1 The Interior-Angle Constraint "On/Off" → Continuous Knob $\lambda$

$\lambda\equiv\frac{\kappa_a}{\kappa_c}$

| Mathematical version switch | Energy version $\lambda$                                                                                                  |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| Angle constraint "off"      | $\lambda=0$ ($\kappa_a=0$, angle force off)                                                                               |
| Angle constraint "on"       | $\lambda\to\infty$ (any move increasing the angle deviation makes the energy blow up and is rejected by gradient descent) |
| ——                          | **Intermediate** **$\lambda$**: a new regime the mathematical version cannot reach                                        |

**$\lambda$** **sweep**: fix $\kappa_c=1$ and let $\kappa_a$ (i.e. $\lambda$) scan over a sequence of values (e.g. $0,\,0.2,\,0.5,\,1,\,2,\,5,\,10,\,50$), relaxing to steady state at each $\lambda$ and recording statistics. The two endpoints should respectively reproduce the mathematical version's "off / on".

***

## 7. Intrinsic Properties and Limits of the Model (Boundaries Must Be Stated)

1. **Energy is dimensionless and strictly scale-invariant**: $E_{\mathrm{cen}},E_{\mathrm{ang}}$ depend only on angles → $E(\lambda\mathbf r)=E(\mathbf r)$.
2. **Force decays with scale**: $|\mathbf F|\sim\kappa\,\delta/\rho$, i.e. the force magnitude is inversely proportional to the tissue's length scale ($\propto 1/\text{length}$).
3. **Three zero modes**: the energy is insensitive to the three global transformations of translation, rotation, and uniform expansion, and has no restoring force for any of them.

| Transformation                                    | Basis for invariance                                                                               | Equivalent identity                                                                           |
| ------------------------------------------------- | -------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| Translation $\mathbf r_v\to\mathbf r_v+\mathbf c$ | $E$ depends only on **relative** positions                                                         | $\sum_v\mathbf F_v=\mathbf 0$ (§3.1 cell-wise sum zero; §3.2 three-vertex-per-angle sum zero) |
| Rotation $\mathbf r_v\to R\mathbf r_v$            | $\varphi_c$ absorbs the overall rotation (§2.2); interior angles are themselves rotation-invariant | $\sum_v\mathbf r_v\times\mathbf F_v=0$ (no net torque)                                        |
| Expansion $\mathbf r_v\to\lambda\mathbf r_v$      | the **scale invariance** of §7.1                                                                   | $\sum_v\mathbf r_v\cdot\mathbf F_v=0$ (expansion does no work)                                |

**Proof of the expansion zero mode**: differentiate $E(\lambda\mathbf r)=E(\mathbf r)$ w\.r.t. $\lambda$ at $\lambda=1$:
$0=\frac{d}{d\lambda}E(\lambda\mathbf r)\Big|_{\lambda=1}=\sum_v\nabla_{\mathbf r_v}E\cdot\mathbf r_v=-\sum_v\mathbf F_v\cdot\mathbf r_v$

> ⚠️ **Easy to get wrong**: "force sum is zero" proves the **translation** zero mode and **cannot** be used to argue the expansion zero mode — the latter must use the virial-type identity above, $\sum_v\mathbf F_v\cdot\mathbf r_v=0$.

- Corollary: **the absolute size of the tissue is not determined by this energy** but by initial conditions and boundaries. (In the classic vertex model the size is fixed by the balance of $\Lambda$ and $\Gamma$; this model has no line tension, hence no such mechanism.)
- Corollary: **there is no intrinsic length scale**; the model fixes shape, not size.
- Corollary: combined with §7.2, a larger tissue has weaker force → **absolute displacement is slower**. Under the mobility-$\mu$ formulation this is a genuine property of the model, not a numerical artifact.

1. **No line tension → no coarsening**: a pure symmetry energy only relaxes to a frustrated minimum; it cannot drive T1 topological transitions or von Neumann–Mullins coarsening. This is a limitation inherent to symmetric relaxation.
2. **Effect of centroid choice** (see §8.1): under the vertex centroid the zero-sum holds and the gradient formula is exact; under the area centroid it no longer holds strictly.

***

## 8. Parameter Tables

### 8.1 Model and Physical Parameters

| Symbol                             | Name                       | Definition / value                                                                               | Meaning                                                                                                              |
| ---------------------------------- | -------------------------- | ------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------- |
| $E$                                | Total energy               | $\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$                                              | The scalar objective to minimize (dimensionless)                                                                     |
| $E_{\mathrm{cen}}$                 | Center-angle energy        | $\sum\delta^2$                                                                                   | Cell level: equal spacing of directions                                                                              |
| $E_{\mathrm{ang}}$                 | Interior-angle energy      | $\sum(\alpha-\tau)^2$                                                                            | Vertex level: angular regularity                                                                                     |
| $\kappa_c$                         | Center-angle stiffness     | $>0$, usually normalized to 1                                                                    | Weight of the center-angle term                                                                                      |
| $\kappa_a$                         | Interior-angle stiffness   | $\ge 0$                                                                                          | Weight of the interior-angle term                                                                                    |
| $\lambda$                          | Relative angle weight      | $\kappa_a/\kappa_c$                                                                              | **Continuous knob**; 0 = constraint off, $\infty$ = constraint on                                                    |
| $\mathbf O_c$                      | Cell center                | see below                                                                                        | Pole for measuring center angles                                                                                     |
| $n_c$                              | Cell vertex count          | integer                                                                                          | Determines the equal-spacing angle $2\pi/n_c$                                                                        |
| $\psi^{\uparrow}_{c,i}$            | Lifted direction angle     | $\arg(\mathbf r_{c,i}-\mathbf O_c)$, cyclically unfolded                                         | Basic observable of the center angle                                                                                 |
| $\varphi_c$                        | Best-fit phase             | $\frac1{n_c}\sum(\psi^\uparrow-2\pi k/n_c)$                                                      | **Measurement reference / gauge fixing** (absorbs overall rotation), **not** a target configuration; arithmetic mean |
| $\delta_{c,i}$                     | Center-angle deviation     | $\psi^\uparrow_{c,i}-\varphi_c-2\pi i/n_c$                                                       | Basic quantity of $E_{\mathrm{cen}}$; $\sum_i\delta\equiv0$                                                          |
| $\rho_{c,i},\hat{\mathbf u}_{c,i}$ | Radius, radial unit vector | $\lvert\mathbf u\rvert,\ \mathbf u/\lvert\mathbf u\rvert$                                        | Force decays as $1/\rho$                                                                                             |
| $J$                                | +90° rotation              | $J(x,y)=(-y,x)$                                                                                  | turns radial into tangential (azimuthal force)                                                                       |
| $M_c$                              | Cell shared-average term   | $\frac1{n_c}\sum\delta_j J\hat{\mathbf u}_j/\rho_j$                                              | **Vertex-binding term**; identical for all vertices of the cell                                                      |
| $\tau$                             | Target interior angle      | $2\pi/3$ interior; $\pi/2$ boundary                                                              | Boundary condition derived from foam force balance                                                                   |
| $\alpha,\beta$                     | Interior angle             | $\psi^{\uparrow}_a-\psi^{\uparrow}_b$, **order must not be reversed** (§2.3)                     | Basic quantity of $E_{\mathrm{ang}}$; $\mathbf a$ = arm to predecessor, $\mathbf b$ = arm to successor               |
| $\mathbf F_v$                      | Vertex force               | $-\nabla_{\mathbf r_v}E$                                                                         | The quantity driving vertex motion                                                                                   |
| $\mu$                              | Mobility                   | $\Delta t/\zeta>0$                                                                               | **Step-size constant** (the only time-scale parameter)                                                               |
| $\zeta$                            | Damping coefficient        | material friction (dimensionless in this model, $\zeta$ has no independent physical calibration) | merged with $\Delta t$ into $\mu$                                                                                    |
| $\hat{\mathbf t}$                  | Boundary tangent           | unit direction of the boundary edge associated with the boundary vertex                          | a boundary vertex can only slide along this direction                                                                |

### 8.2 On the Choice of Centroid (Decided: Vertex Centroid)

| Scheme                       | Gradient                                                         | Force zero-sum (pure azimuthal, expansion does no work) | vs. original mathematical version                               |
| ---------------------------- | ---------------------------------------------------------------- | ------------------------------------------------------- | --------------------------------------------------------------- |
| **Vertex centroid (chosen)** | formula §3.1 holds exactly                                       | holds strictly                                          | the original used area centroid; the model uses vertex centroid |
| Area centroid                | formula not exact (needs extra centroid-motion correction terms) | only approximate                                        | consistent with the original                                    |

**Reason for choice**: exact gradient + the clean property "center-angle force is purely azimuthal and uniform expansion does no work"; the difference from the original is only an $O(\sigma/R)$ geometric fluctuation.

### 8.3 Algorithm Control Parameters (To Be Determined)

| Parameter                      | Meaning                                            | TBD |
| ------------------------------ | -------------------------------------------------- | --- |
| $\mu$                          | initial mobility / whether to schedule over rounds | TBD |
| Convergence tolerance          | energy-change threshold                            | TBD |
| Max rounds                     | upper bound on single-round iterations             | TBD |
| Max line-search halvings       | upper bound on step protection                     | TBD |
| $\lambda$ sweep value sequence | list of $\kappa_a$ for the sweep                   | TBD |

***

## 9. List of Settled Decisions

1. The energy contains only two terms: $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$, **no line tension** **$\Lambda$**.
2. Centroid is the **vertex centroid**.
3. $\varphi_c$ takes the **arithmetic mean** (not circular-mean).
4. Interior target angle $2\pi/3$, boundary target angle $\pi/2$, both **derived from foam force balance**, not empirical constants.
5. Motion = **over-damped gradient descent** (forward Euler), $\Delta\mathbf r=\text{step}\cdot\mathbf F$; **introduce** **$L_0$** **median-edge-length normalization**: $\text{step}=\mu L_0^{2}$, per-round displacement $\approx\mu L_0$ (on by default; can be turned off via `PhysicalAnnealerParams.normalize_L0=False`).
6. Vertices use **descending-** **$|\mathbf F_v|$** **Gauss–Seidel serial update** (on by default); the convexity guard is localized to incident cells. The synchronous Jacobi update remains an optional mode (set `serial_descending=False`) but easily freezes large random networks, so it is not recommended.
7. Boundary vertices **slide along the boundary tangent** (tangent projection).
8. Line search: **under the serial scheme each vertex does its own**; the convexity guard is localized to the cells incident to that vertex; cells already non-convex at the start are exempt (when `serial_descending=False` this degrades to a global line search).
9. The interior-angle constraint is upgraded from a binary switch to a **continuous knob** **$\lambda=\kappa_a/\kappa_c$**.
10. $\varphi_c$ and "optimal rays" share the same solution, but **only the phase is taken**; its role is **gauge fixing** (absorbing overall rotation, guaranteeing the energy's rotation invariance), **not** a target configuration. Neither the ray geometry nor a target point is constructed.
11. $E_{\mathrm{cen}}$ **constrains direction only, not radius** (the center-angle force is purely azimuthal); cell size is determined by $E_{\mathrm{ang}}$ and initial conditions.
12. The interior-angle gradient **must include the arm-endpoint term** (§3.2); its scope is the **one-hop neighborhood**, a fixed 9 terms (interior) / 6 terms (boundary), **no cascade, no** **$O(V^2)$**.
13. The interior-angle force is **implemented in edge (pairwise) form** (§3.2.2): traverse cell-by-cell CCW, scatter each edge by $(w_{X,c}-w_{Y,c})J\hat{\mathbf e}/l$ pairwise to the two endpoints. Same-order complexity as the apex version, and $\sum_v\mathbf F_v=\mathbf 0$ holds explicitly.
14. The interior angle is defined as $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$ (**swept CCW from the outgoing-edge direction to the reversed incoming edge**). Both the §3.2 gradient table and the §3.2.2 edge form take this as the sole basis.
15. $E_{\mathrm{cen}}$ and $E_{\mathrm{ang}}$ are always **bare sums without stiffness**; $\kappa_c,\kappa_a$ appear only once, in the total energy $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$.
16. A boundary vertex has **2 boundary edges**; the sliding direction is the tangent of the edge with the larger $\lvert\mathbf F\cdot\hat{\mathbf t}\rvert$ among the two (§4.2).
17. The model has **three zero modes**: translation, rotation, expansion. The expansion zero mode is given by the virial identity $\sum_v\mathbf F_v\cdot\mathbf r_v=0$ and **cannot** be argued from "force sum is zero" (which only proves the translation zero mode).
18. $\lambda$ **defaults to 1**; supports a **dynamic schedule** (homotopy / continuation method): once enabled, $\lambda$ decays linearly from `lam_start` to `lam_end` (default $0\to1$), recomputed each round — first pushing vertices toward the 120° manifold with a larger $\lambda$ as warm-up, then dropping to 1 for fine convergence. The static-$\lambda$ mode leaves the schedule unchecked.

***

## 10. Revision History

- **Correction**: in the center-angle force formula, $\delta_{c,i}$ multiplies **only the first term**, not the whole parenthesis.
  Correct: $\mathbf F^{(c)}=-2\kappa_c\left(\delta_{c,i}\dfrac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)$
  Wrong (deprecated): $\mathbf F^{(c)}=-2\kappa_c\,\delta_{c,i}\left(\dfrac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)$
  Basis: strictly derived from the vertex centroid $+\sum_i\delta_{c,i}\equiv0$; the former satisfies the cell force zero-sum, the latter does not.
- **Correction**: the step removed the median edge length $L_0$, switching to the pure-mobility form $\Delta\mathbf r_v=\mu\mathbf F_v$. *(Historical note — later superseded by the L0-normalization decision in §9 #5.)*
- **Addition**: the interior-angle gradient includes the **arm-endpoint contribution** (§3.2); apex-only is incomplete.
- **Clarification (§2.2)**: $\varphi_c$ and "optimal rays" have the **same solution but only the phase is taken**; its identity is **gauge fixing**, not a target configuration — it absorbs the overall rotation, guaranteeing $E_{\mathrm{cen}}$ measures shape asymmetry, not orientation. Fixing $\varphi_c\equiv0$ would assign spurious energy to a perfectly rotated regular polygon.
- **Clarification (§2.2)**: $E_{\mathrm{cen}}$ **constrains direction only, not radius**. The mathematical version's "rays + triangle centroid" fixes both direction and position together; after discarding the target point, the center-angle force is always a pure azimuthal force perpendicular to the radial direction, doing no radial work, so cell size is not determined by $E_{\mathrm{cen}}$. This is the equivalent statement of §7.3's "uniform expansion is a zero mode".
- **Clarification (§3.2.1)**: the interior-angle gradient has **no cascading propagation**. The gradient is first-order and accumulates only terms that explicitly contain $\mathbf r_v$; "moving $v$ changes neighbors' angles" is an outcome, already recorded exactly to first order by the fixed 9 (interior) / 6 (boundary) terms in the one-hop neighborhood, with farther effects handled by iteration. Complexity is $O(1)$ per vertex, $O(V)$ globally.
- **Addition (§3.2.2)**: the **edge (pairwise) equivalent form** of the interior-angle gradient $\mathbf g_c(X\!\to\!Y)=(w_{X,c}-w_{Y,c})J\hat{\mathbf e}_{X\to Y}/l$, obtained by one CCW traversal per cell, with apex and arm-endpoint contributions merged automatically. Plus: complexity argument and the automatic $\sum_v\mathbf F_v=\mathbf 0$.
- **Addition (§3.2.3)**: the three reasons the arm-endpoint term **must not be omitted** (variational principle fails / spurious bulk drift / line search stalls) and the magnitude argument "same order, not a higher-order correction"; the degraded scheme (block-local rule) is explicitly not adopted.
- **Correction (sign error, severe)**: the interior-angle definition changed from $\alpha=\psi^{\uparrow}_b-\psi^{\uparrow}_a$ to $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$. The old form gives $-\alpha$ (the opposite of the interior angle), flipping all three rows of the §3.2 gradient table in sign and contradicting the §3.2.2 edge form (the two differ by an overall sign). Basis: under CCW vertex order $\mathbf a$=arm to predecessor, $\mathbf b$=arm to successor, $\psi_a=\theta_{\mathrm{in}}+\pi$, $\psi_b=\theta_{\mathrm{out}}$, hence $\alpha=\pi-(\theta_{\mathrm{out}}-\theta_{\mathrm{in}})=\psi_a-\psi_b$. Hand check on a square: $\nabla_{\mathbf r_i}\alpha_i=\dfrac{J\hat{\mathbf e}_{\mathrm{in}}}{l_{\mathrm{in}}}+\dfrac{J\hat{\mathbf e}_{\mathrm{out}}}{l_{\mathrm{out}}}$, consistent with the direct computation that the interior angle increases when the vertex moves along $(-1,1)$. Consequence: signs in the §3.2 three-row gradient table, the §2.3 interior-angle definition, and the §8.1 parameter table were all corrected. The §3.2.2 edge-form **conclusion is unchanged** (verified independently and by sign self-check: zero contribution when $w_X=w_Y$).
- **Correction (stiffness double-counting)**: the boxed formulas in §2.2/§2.3 originally wrote $E_{\mathrm{cen}},E_{\mathrm{ang}}$ as $\kappa\sum(\cdot)^2$, so stacking onto the total energy $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$ double-counted stiffness (producing $\kappa^2$). Now unified as **bare sums** $E_{\mathrm{cen}}=\sum\delta^2$, $E_{\mathrm{ang}}=\sum(\alpha-\tau)^2$, with stiffness appearing only once in the total energy, consistent with the gradient weights $2\kappa_c\delta$, $2\kappa_a(\alpha-\tau)$ of §3.
- **Correction (§3.3 missing stiffness)**: $\mathbf F^{\mathrm{ang}}_v$ was written as $-\nabla_{\mathbf r_v}E_{\mathrm{ang}}$, missing $\kappa_a$; now corrected to $-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}=-\mathbf G_v$, with the note that §3.2.2 accumulates $\mathbf G_v$ ($w$ already contains $\kappa_a$).
- **Correction (§4.2 tangent undefined)**: a boundary vertex actually has **2 boundary edges** (§1 added: 2 boundary edges + 1 internal edge), and the tangent is not unique at a kink. Now the rule is to take the larger $\lvert\mathbf F\cdot\hat{\mathbf t}\rvert$, and it is stated explicitly that this is a **dynamic constraint, not an energy term**.
- **Correction (§7.3 zero-mode argument error)**: the original used "force sum is zero" to argue the expansion zero mode, but that identity proves the **translation** zero mode. Now split into the three zero modes of translation / rotation / expansion; the expansion zero mode is proven instead by the virial identity $\sum_v\mathbf F_v\cdot\mathbf r_v=0$ (obtained by differentiating the scale invariance w\.r.t. $\lambda$); the rotation zero mode is newly added.
- **Correction (§6 correspondence table)**: the boundary row originally read "steepest descent on $(\beta_1-\beta_2)^2$", missing the boundary-bending term $\tfrac12(\beta_1+\beta_2-\pi)^2$; now corrected to the full boundary interior-angle energy. Also added a "descending-annealing-distance serial ↔ update order" row (that row's first draft read "synchronous update", later corrected in this round to serial Gauss–Seidel update in descending $|\mathbf F_v|$ order).
- **Clarification (§3.2.2 traversal scope)**: an internal edge is traversed once by **each of its two cells and both must be accumulated** (different weights, no de-duplication); marginal cells' boundary edges also participate in the traversal. Plus added the step-by-step equivalence derivation and the sign self-check.
- **Correction (§3.2 boundary-vertex term count note)**: the count of 6 is correct, but the old explanation "2 edges each belonging to 2 cells" was wrong; in fact the internal edge belongs to 2 cells (2 terms) + the 2 boundary edges each belong to 1 cell (1 term each) = 4 terms.
- **Clarification (§5 flow)**: step 1 records the **2** boundary-edge directions of boundary vertices; step 3 accumulates **force** not gradient; step 7's recomputed **area does not enter the energy** (display/export only), and the centroid is computed directly from vertex coordinates in step 2 of the next round.
- **This-round revision (aligned with code / introduction)**: ① the update order was changed from "synchronous update" to **Gauss–Seidel serial update in descending** **$|\mathbf F_v|$** **order** (default), with the convexity guard localized to incident cells (§4.3, §5 step 6, §6 correspondence table, §9 decision #6); ② **L0 median-edge-length normalization** was introduced ($\text{step}=\mu L_0^{2}$, displacement $\approx\mu L_0$, on by default, `normalize_L0=False` to turn off), superseding the old "no L0" wording (§4.1, §9 decision #5); ③ the line search was changed from "whole-step global veto" to **per-vertex, with localized convexity guard** (§4.3, §9 decision #8); ④ $\lambda$ **default = 1** and a **dynamic schedule** $0\to1$ (homotopy, recomputed each round) were added to §9 decision #18; ⑤ the synchronous Jacobi update is retained as an optional mode (`serial_descending=False`).

***

# 对称退火 · 数学版 → 能量版 完整逻辑与公式规范

> 能量版对称退火算法由 Kai Xu 借助 Hy4 大语言模型（LLM）将数学版对称退火算法（Xu et al. 2026, arXiv:2606.18604）转换为能量版得到。
> 本文只确定算法逻辑、公式与参数定义。不含任何验证、测试、数值实验内容。
> 逻辑确认后再进入实现层。

***

## 0. 定位与总纲

**数学版**：一组几何规则——拟合 $n_c$ 条最优射线 → 取三角形质心为目标点 → 移动；边缘点沿边界滑动；用「移动若使 $\sum\alpha^2$ 增大则取消」做角度硬把关。

**能量版**：找一个标量能量 $E$，使上述规则成为「对 $E$ 做梯度下降（过阻尼运动）」的特例。原算法从「规则集」升级为「变分原理的数值实现」。

总纲一句话：
$E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}},\qquad \mathbf F_v=-\nabla_{\mathbf r_v}E,\qquad \Delta\mathbf r_v=\mu\,\mathbf F_v$

***

## 1. 拓扑与记号约定

| 记号   | 含义                                                                                   |
| ---- | ------------------------------------------------------------------------------------ |
| 内部顶点 | 被 3 个细胞共享（度 3），3 条关联边连 3 个邻居                                                         |
| 边缘顶点 | 被 2 个细胞共享（度 2），贴自由边界，沿边界边滑动。有 **3 条关联边**：2 条**边界边**（各只属于 1 个细胞）+ 1 条**内边**（属于 2 个细胞） |
| 内部细胞 | 全部顶点都是内部顶点                                                                           |
| 边缘细胞 | 含 2 个边缘顶点 + 平均约 3 个内部顶点                                                              |

- 细胞 $c$ 的顶点按\*\*逆时针（CCW）\*\*记为 $\mathbf r_{c,0},\dots,\mathbf r_{c,n_c-1}$，$n_c$ 为顶点数。
- **每个细胞各自持有一份顶点坐标副本**；同一物理顶点在不同细胞中坐标重合，需按键值聚类成唯一顶点。
- $J$ 表示 **+90° 旋转**算子：$J(x,y)=(-y,x)$。
- 角度取**提升角**（lifted / unwrapped），即按循环序连续展开，不取模 $2\pi$。

***

## 2. 能量定义

### 2.1 总能量

$E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$

- **无量纲**（只依赖角度，见 §7），因此没有物理的力/能量量纲。
- 不含线张力 $\Lambda\sum l_e$（不属于对称松弛，仅在 Discussion 中作为扩展提及）。

### 2.2 中心角项 $E_{\mathrm{cen}}$（细胞层级）

对细胞 $c$：

- 质心（**顶点质心**，见 §8.1）：
  $\mathbf O_c=\frac{1}{n_c}\sum_{k=0}^{n_c-1}\mathbf r_{c,k}$
- 顶点 $i$ 的径向量与方向角：
  $\mathbf u_{c,i}=\mathbf r_{c,i}-\mathbf O_c,\qquad \rho_{c,i}=|\mathbf u_{c,i}|,\qquad \hat{\mathbf u}_{c,i}=\frac{\mathbf u_{c,i}}{\rho_{c,i}},\qquad \psi^{\uparrow}_{c,i}=\arg(\mathbf u_{c,i})$

> 其中 $\psi^{\uparrow}_{c,i}$ 为**提升角**（按循环序连续展开，不取模 $2\pi$，见 §1）。

- **最佳拟合相位**（算术平均，闭式）：
  $\varphi_c=\frac{1}{n_c}\sum_{k=0}^{n_c-1}\left(\psi^{\uparrow}_{c,k}-\frac{2\pi k}{n_c}\right)$
- **中心角偏差**：
  $\boxed{\ \delta_{c,i}=\psi^{\uparrow}_{c,i}-\varphi_c-\frac{2\pi i}{n_c}\ }$
- **中心角能量**（**裸和**，不含 $\kappa_c$）：
  $\boxed{\ E_{\mathrm{cen}}=\sum_c\sum_{i=0}^{n_c-1}\delta_{c,i}^2\ }$

> **重要约定**：$E_{\mathrm{cen}}$、$E_{\mathrm{ang}}$ 一律定义为**不含刚度的裸和**；刚度只出现在总能量 $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$（§2.1）中。这样才能与 §3 的梯度权重 $2\kappa_c\delta$、$2\kappa_a(\alpha-\tau)$ 一致。若把 $\kappa$ 写进 $E_{\mathrm{cen}}$/$E_{\mathrm{ang}}$ 的定义，$E$ 里就会出现 $\kappa^2$，§3 全部公式作废。

**恒等式**：$\displaystyle\sum_{i}\delta_{c,i}\equiv 0$（由 $\varphi_c$ 是算术平均直接推出）。这条恒等式是后面梯度简洁、力纯方位的关键。

**含义**：$\delta_{c,i}$ 度量「顶点 $i$ 的实际方向」相对「$n_c$ 个等间隔方向的最佳排布」的偏离。$E_{\mathrm{cen}}=0\iff$ 各相邻顶点方向夹角都等于 $2\pi/n_c$。

#### 「最佳拟合相位」与「最优射线」的关系

数学版「找最优射线」的完整含义是：从 $\mathbf O_c$ 出发作 $n_c$ 条**刚性**射线（相邻夹角固定为 $2\pi/n_c$），整个扇形只允许**整体旋转**一个角度 $\varphi$，使其最贴合实际的 $n_c$ 个顶点方向。最小化 $\sum_i(\psi^{\uparrow}_{c,i}-\varphi-2\pi i/n_c)^2$ 得到的正是上式 $\varphi_c$。

> **刚性是中心角项的本质**：若每条射线可各自转动去对准自己的顶点，残差恒为零，「中心角」这个概念即不存在。唯一的自由度就是整体相位。

因此「算 $\varphi_c$」与「找最优射线」**是同一个最小二乘问题的解**，但只包含三层产出中的第一层：

| 层                      | 数学版用它做什么 | 能量版                |
| ---------------------- | -------- | ------------------ |
| ① 相位 $\varphi_c$       | 定射线的朝向   | **保留**：作为测量偏差的参考基准 |
| ② 射线的完整几何（$n_c$ 条半直线）  | 画出来      | **丢弃**：不构造任何射线     |
| ③ 射线相交成三角形、取质心作**目标点** | 顶点直接移过去  | **丢弃**：不生成目标点      |

即：**能量版继承解，丢弃几何构造**——只问「实际方向偏了多少」，不问「理想位置在哪」。

#### 与数学版最小二乘的具体代数等价

数学版射线拟合常写作升序角度形式：把顶点方向按 $\theta_1<\theta_2<\cdots<\theta_{n_c}$ 排序，令最优射线方向为 $\delta_i=\delta_1+2\pi(i-1)/n_c$，并取

$$\delta_1=\frac{1}{n_c}\Bigl(\sum_{j=1}^{n_c}\theta_j-\pi(n_c-1)\Bigr)$$

以最小化 $\sum_{j=1}^{n_c}(\theta_j-\delta_j)^2$。对该目标关于 $\delta_1$ 求导并令为零，即 $\delta_1=\frac{1}{n_c}\sum_j\bigl(\theta_j-2\pi(j-1)/n_c\bigr)$——这与代码中的 $\varphi_c=\frac{1}{n_c}\sum_k(\psi^{\uparrow}_{c,k}-2\pi k/n_c)$ **代数恒等**（顶点按循环序排列并提升后，指标 $i=k$ 与 $j-1$ 之差无关紧要，见 §1）。

于是，数学版**一次性**最小化的那个量 $\sum_{j=1}^{n_c}(\theta_j-\delta_j)^2$，**正是**单细胞中心角能量 $E_{\mathrm{cen},c}=\sum_i\delta_{c,i}^2$。换言之：**$E_{\mathrm{cen}}$ 就是数学版射线拟合的残差，在每一步迭代中持续求值。** 数学版只拿残差来定射线朝向、再用作图（三角形目标点）把它消除；能量版把残差本身立为目标，用梯度下降持续压低它。

**对每个（细胞, 顶点），偏差只对「该细胞自己的那一条射线」量。** 内顶点会碰到 3 条射线（来自其 3 个所属细胞）、边缘顶点碰到 2 条；但每条射线只充当**其所在细胞那一项**的量角基准，射线之间**从不拼成三角形或目标点**——那个组合正是上文丢弃的几何构造。

#### $\varphi_c$ 的角色：规范固定（gauge fixing），而非目标构型

$E_{\mathrm{cen}}$ 度量的是「顶点方向相对一个刚性等角扇形的偏差」，而该扇形的**整体朝向是无关的**——细胞整体旋转不应改变它的对称程度。$\varphi_c$ 的作用正是吸收这个多余的朝向自由度。

- 若强行取 $\varphi_c\equiv0$，则一个**完美规则**的六边形只要整体旋转 $30^\circ$，就会算出巨大的 $E_{\mathrm{cen}}$，纯粹因为模板没跟着转。
- 故 $\varphi_c$ 不是可选装饰，而是**能量旋转不变性的保证**：$E_{\mathrm{cen}}$ 只在模掉整体旋转之后才是对称性的度量。
- 对比：**数学版视** **$\varphi_c$** **为「目标」（要变成的样子）；能量版视** **$\varphi_c$** **为「尺子」（量偏差用的基准）。** 同一个数，语义完全不同。

#### 只约束方向，不约束半径

- 射线是**半直线**，只给方向、不给顶点沿射线站多远。数学版因此需额外用「射线相交成三角形、取质心」补一个确定的**位置**——那一步同时定下了方向与半径。
- 能量版的 $E_{\mathrm{cen}}$ 只依赖方向角 $\psi^{\uparrow}$，**完全不依赖** **$\rho$**。反映在力（§3.1）上即中心角力恒沿 $J\hat{\mathbf u}$ 方向、垂直于径向 → **纯方位力，不做径向功**。
- 推论：**细胞大小（半径）不由** **$E_{\mathrm{cen}}$** **决定**，而由 $E_{\mathrm{ang}}$ 与初始条件间接确定。这与 §7.3「均匀膨胀是零模、绝对尺寸不由能量决定」是同一事实的两个侧面。

### 2.3 内角项 $E_{\mathrm{ang}}$（顶点层级）

对每一个**内角**（每个细胞的每个顶点处各有一个内角）：

- 顶点 $w$ 处的内角由两条边张成：$\mathbf a=\mathbf r_p-\mathbf r_w$（指向逆时针前驱 $p$），$\mathbf b=\mathbf r_q-\mathbf r_w$（指向逆时针后继 $q$）。
- 内角
  $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$
  其中提升角取值满足 $\psi^{\uparrow}_a\in\bigl(\psi^{\uparrow}_b,\ \psi^{\uparrow}_b+2\pi\bigr)$。**方向不可写反**：它是从「出边方向 $\mathbf b$」逆时针扫到「入边反向 $\mathbf a$」的角，即 CCW 多边形**内侧**的夹角。
  等价写法：$\alpha=\pi-(\theta_{\mathrm{out}}-\theta_{\mathrm{in}})$（$\theta$ 为 CCW 有向边的提升方位角），凸顶点下 $\alpha\in(0,\pi)$。

> **常见错误**：写成 $\alpha=\psi^{\uparrow}_b-\psi^{\uparrow}_a$。这给出的是 $-\alpha$，会使 §3.2 全部梯度符号翻转（§3.2.2 的边形式也随之翻转）。

**目标角** $\tau$ 按顶点度数定：

| 顶点类型 | 度数 | $\tau$             | 来源                          |
| ---- | -- | ------------------ | --------------------------- |
| 内部顶点 | 3  | $2\pi/3=120^\circ$ | 泡沫三膜力平衡（Plateau），**导出**     |
| 边缘顶点 | 2  | $\pi/2=90^\circ$   | 自由边界正交条件（120° 的边界类比），**导出** |

$\boxed{\ E_{\mathrm{ang}}=\sum_{\text{corner}}(\alpha-\tau)^2\ }$

> 求和遍历**每一个内角**（每个细胞的每个顶点处各一个内角）。同样为**裸和**，刚度 $\kappa_a$ 只在总能量中出现。

**与数学版的关键等价**（内部顶点）：三个内角和为 $2\pi$，故
$\sum_{k=1}^{3}(\alpha_k-\tfrac{2\pi}{3})^2=\sum_k\alpha_k^2-\frac{4\pi}{3}\underbrace{\sum_k\alpha_k}_{=2\pi}+3\left(\frac{2\pi}{3}\right)^2=\sum_k\alpha_k^2+\text{const}$
即 **最小化** **$E_{\mathrm{ang}}$** **≡ 最小化** **$\sum\alpha^2$**，而 $\sum\alpha^2$ 正是数学版角度硬约束所守的量。

**边缘顶点的分解**：两个边界内角 $\beta_1,\beta_2$ 对 $90^\circ$ 取值时
$(\beta_1-\tfrac{\pi}{2})^2+(\beta_2-\tfrac{\pi}{2})^2=\tfrac12(\beta_1-\beta_2)^2+\tfrac12(\beta_1+\beta_2-\pi)^2$

- $\tfrac12(\beta_1-\beta_2)^2$：两条膜关于边界对称（边界正交项）。
- $\tfrac12(\beta_1+\beta_2-\pi)^2$：边界离散转角（自由边界弯曲项）。
- 第三个角（两内边夹角 $=2\pi-(\beta_1+\beta_2)$）**不设罚项**，由几何自行决定。

***

## 3. 力的完整公式

顶点力 = 总能量对位置的负梯度：
$\boxed{\ \mathbf F_v=-\nabla_{\mathbf r_v}E=\underbrace{\left(-\kappa_c\sum_{c\ni v}\nabla_{\mathbf r_v}E_{{\mathrm{cen}},c}\right)}_{\text{centre-angle force}}+\underbrace{\left(-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}\right)}_{\text{interior-angle force}}\ }$

> 右端两项即**中心角力**（§3.1）与**内角力**（§3.2）。

### 3.1 中心角力（**核心公式，注意括号位置**）

定义细胞 $c$ 的**全细胞共享平均项**（对所有顶点相同，不含下标 $i$）：
$\boxed{\ M_c=\frac{1}{n_c}\sum_{j=0}^{n_c-1}\delta_{c,j}\frac{J\hat{\mathbf u}_{c,j}}{\rho_{c,j}}\ }$

则细胞 $c$ 中第 $i$ 个顶点所受的中心角力为：
$\boxed{\ \mathbf F^{(c)}_{v}=-2\kappa_c\left(\delta_{c,i}\frac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)\ }$

- **直接项** $\delta_{c,i}J\hat{\mathbf u}_{c,i}/\rho_{c,i}$：只关该顶点自己；方向垂直于径向（方位力），大小 $\propto$ 该顶点的角偏差、$\propto 1/\rho$。
- **共享项** $M_c$：对细胞 $c$ 的**每一个顶点都相同** → 把全细胞顶点朝同一方向整体拉。这就是\*\*「中心角把细胞所有顶点绑定」的变分表达\*\*。

**顶点归属计数**（力的求和份数）：

- 内部顶点 → 所属 **3 个**细胞的中心角力之和；
- 边缘顶点 → 所属 **2 个**细胞的中心角力之和（每个细胞仍绑定其 2 边缘 + 约 3 内部顶点）。

**力的零和性质**（顶点质心 + $\sum\delta\equiv0$ 的直接后果）：
$\sum_{i=0}^{n_c-1}\mathbf F^{(c)}_{v_i}=-2\kappa_c\left(\sum_i\delta_{c,i}\frac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-n_cM_c\right)=\mathbf 0$
即中心角力是**纯方位力**：对一个细胞求和为零，**均匀膨胀不做功**。

### 3.2 内角力（含端点贡献，容易漏）

对单个内角 $\alpha$（顶点 $w$ 处，两臂指向 $p,q$，$\mathbf a=\mathbf r_p-\mathbf r_w$，$\mathbf b=\mathbf r_q-\mathbf r_w$），记权重
$w_\alpha=2\kappa_a(\alpha-\tau)$

该内角对三个相关顶点的**梯度贡献**：

**推导**（$J$ 为 +90° 旋转）：$\nabla_{\mathbf r_w}\arg(\mathbf r_p-\mathbf r_w)=-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$、$\nabla_{\mathbf r_p}\arg(\mathbf r_p-\mathbf r_w)=+\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$，对 $\mathbf b$ 同理；再代入 $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$（§2.3）即得下表。

| 顶点  | 角色           | 梯度贡献 $\nabla(\alpha)$                                                                              | 对 $\kappa_a\nabla E_{\mathrm{ang}}$ 的贡献                                                                                 |
| --- | ------------ | -------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| $w$ | 角顶           | $-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}+\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$ | $w_\alpha\left(-\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}+\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}\right)$ |
| $p$ | 前臂端点（CCW 前驱） | $+\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$                                                 | $+w_\alpha\dfrac{J\hat{\mathbf a}}{\lvert\mathbf a\rvert}$                                                              |
| $q$ | 后臂端点（CCW 后继） | $-\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$                                                 | $-w_\alpha\dfrac{J\hat{\mathbf b}}{\lvert\mathbf b\rvert}$                                                              |

**力 = 上述贡献取负**，即 $\mathbf F^{\mathrm{ang}}=-\kappa_a\nabla E_{\mathrm{ang}}$。

- 三项之和为零 → 内角力是**内力**，不产生净平移。
- **要点**：动一个顶点 $v$，不只改变它**自己头上**的角，还改变它作为**臂端点**参与的、位于其 3 个邻居处的角。因此装配 $\mathbf F_v$ 时，必须把「$v$ 为角顶的角」与「$v$ 为臂端点的角」两类贡献**全部累加**，否则梯度不完整。

**每个顶点涉及的内角数（固定，与组织规模无关）**：

- 内部顶点（度 3，3 条**内边**）：作角顶 3 个（每细胞一个）+ 作臂端点 6 个（3 条边，每条边属于 2 个细胞）→ **共 9 项**。
- 边缘顶点（度 2，2 条**边界边** + 1 条**内边**，见 §1）：作角顶 2 个 + 作臂端点 4 个（内边属于 2 个细胞，贡献 2 项；2 条边界边各只属于 1 个细胞，各贡献 1 项）→ **共 6 项**。

### 3.2.1 作用范围：一跳即止，**不存在连锁传递**

$\nabla_{\mathbf r_v}E_{\mathrm{ang}}$ 的定义是「把 $\mathbf r_v$ 当变量、其余全当常数求导」，因此**只累加公式里显式含有** **$\mathbf r_v$** **的项**。单个内角只依赖 3 个位置（角顶 + 两臂端点），故涉及 $\mathbf r_v$ 的角**只有上面列出的那一跳邻域内的固定几项**。

- 「动 $v$ 会改变邻居的角」是**结果**，不是**求导项**；它已被这 9 项（或 6 项）**一阶地、精确地**记录。
- 更远的间接效应（邻居的邻居 …）不出现在 $\partial E/\partial\mathbf r_v$ 的表达式中，由**下一轮迭代**承担——那是迭代的职责，不是单步梯度的职责。
- 因此**没有级联传播，没有** **$O(V^2)$** **或指数爆炸**：每顶点 $O(1)$ 项，全局 $O(V)$。

### 3.2.2 边（成对）等价形式 —— **推荐的实现写法**

把上表三行合并后，同一条边上的「角顶贡献」与「臂端点贡献」权重恰好**相减**。定义每个**角**（corner，即 (顶点, 细胞) 对）的标量权重
$w_{v,c}=2\kappa_a(\alpha_{v,c}-\tau_v)$

并记**含刚度的内角梯度**为
$\mathbf G_v\equiv\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}\qquad(\mathbf F^{\mathrm{ang}}_v=-\mathbf G_v)$

**遍历范围**：逐细胞沿 **CCW** 遍历该细胞的**全部**边——

- 边缘细胞的 **2 条边界边**同样参与（边界顶点处的角照常计权）；
- **内边被它所属的两个细胞各遍历一次**，两次权重不同（$\alpha$ 属不同细胞），**两次都要累加**，不可去重。

对每条有向边 $X\to Y$（长度 $l$，单位矢 $\hat{\mathbf e}_{X\to Y}$），令
$\boxed{\ \mathbf g_c(X\!\to\!Y)=\bigl(w_{X,c}-w_{Y,c}\bigr)\frac{J\hat{\mathbf e}_{X\to Y}}{l}\ }$

则
$\mathbf G_X\;+{=}\;\mathbf g_c(X\!\to\!Y)$
$\mathbf G_Y\;-{=}\;\mathbf g_c(X\!\to\!Y)$

**等价性推导**（设 CCW 序为 $\dots\to X\to Y\to\dots$，该边即 $\mathbf e$；注意 $\mathbf b_X=\mathbf e$、$\mathbf a_Y=-\mathbf e$）：

- **顶点** **$X$**：① 作 $\alpha_X$ 的角顶，取表中角顶行的 $\mathbf b$ 部分 → $+w_X\dfrac{J\hat{\mathbf e}}{l}$；② 作 $\alpha_Y$ 的前臂端点（$X$ 是 $Y$ 的 CCW 前驱，$\mathbf a_Y=-\mathbf e$）→ $+w_Y\dfrac{J(-\hat{\mathbf e})}{l}$。合计 $(w_X-w_Y)\dfrac{J\hat{\mathbf e}}{l}=\mathbf g$。
- **顶点** **$Y$**：① 作 $\alpha_Y$ 的角顶，取角顶行的 $\mathbf a$ 部分 → $-w_Y\dfrac{J(-\hat{\mathbf e})}{l}=+w_Y\dfrac{J\hat{\mathbf e}}{l}$；② 作 $\alpha_X$ 的后臂端点（$\mathbf b_X=\mathbf e$）→ $-w_X\dfrac{J\hat{\mathbf e}}{l}$。合计 $(w_Y-w_X)\dfrac{J\hat{\mathbf e}}{l}=-\mathbf g$。

即 §3.2 表中三种角色（角顶 / 前臂 / 后臂）的贡献**全部包含在内，无遗漏、无重复**。

> **符号自检**：$w_X=w_Y$（两角偏差相同）时该边贡献为零——符合直觉：等权的相邻两角，把公共端点沿边法向推不产生一阶能量变化。

**优点**：

1. **两行代码、一次遍历**；无「前臂 / 后臂」分支判断，比按角顶逐项累加更不易写错。
2. **力零和自动成立**：每条边贡献 $+\mathbf g$ 与 $-\mathbf g$ 一对 → $\sum_v\mathbf F_v=\mathbf 0$ 显式成立，无虚假整体平动（§7.3）。
3. **复杂度同阶**：遍历量为所有 (细胞, 边) 关联数 $=\sum_c n_c=\sum_v\deg(v)\approx 3V$（$\deg$ 按 §1 记为该顶点所属细胞数），与「只算角顶」**同阶**——散射步每次写 2 个向量而非 1 个，常数约 2 倍。真正的成本大头是算 $\alpha$ 本身（约 $3V$ 次反三角函数），两版都要付，故散射方式的差别基本被掩盖。

### 3.2.3 臂端点项**不可省略**（不是精度取舍）

若只保留角顶项，得到的是一个**不同的算子**，代价明确：

1. **变分原理失效**：$\mathbf F_v\neq-\nabla_{\mathbf r_v}E$ → 迭代最小化的**不是**写出的那个 $E$，「原算法是 $E$ 的梯度下降特例」这一核心论断不成立。
2. **虚假整体漂移**：$E_{\mathrm{ang}}$ 只依赖相对位置，完整梯度严格满足 $\sum_v\mathbf F_v=\mathbf 0$；截断后不再满足，组织会凭空平动/转动，且这部分位移会被误读为弛豫进度。
3. **线搜索卡死**：截断方向不是下降方向 → 能量不降 → 步长被迫折半至零 → 停滞。想省事反而更不稳。

**量级上也不可丢**：平衡点附近 $w\propto(\alpha-\tau)$ 为**一阶**小量，臂端点项与角顶项同阶，不是可忽略的高阶修正。

**若坚持要简化**，唯一选项是退化为**分块局部规则**（每顶点只按自己头顶的角移动，Jacobi 式局部投影）。它能跑，但失去 $E$ 作为 Lyapunov 函数、失去单调收敛保证、不能称为变分模型。**本规范不采用。**

### 3.3 装配

$\mathbf F_v=\mathbf F_v^{\mathrm{cen}}+\mathbf F_v^{\mathrm{ang}}$

- $\displaystyle\mathbf F_v^{\mathrm{cen}}=\sum_{c\ni v}\mathbf F^{(c)}_{v}$ —— **中心角力**（§3.1，$\mathbf F^{(c)}$ 已含 $\kappa_c$）。
- $\displaystyle\mathbf F_v^{\mathrm{ang}}=-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}=-\mathbf G_v$ —— **内角力**（§3.2，$E_{\mathrm{ang}}$ 为裸和，刚度 $\kappa_a$ 在此乘入；$w_{v,c}$ 已含 $\kappa_a$，故 §3.2.2 累加得到的是 $\mathbf G_v$ 而非 $\nabla E_{\mathrm{ang}}$）。

**内角力采用边形式装配**：先逐角算 $\alpha_{v,c}$ 与标量权重 $w_{v,c}$，再逐细胞沿 CCW 遍历其边界，按 §3.2.2 把 $\mathbf g_c(X\!\to\!Y)$ 成对散到两端点。这样「角顶 + 臂端点」三类贡献一次到位，无需单独枚举角色。

### 3.4 顶点层力聚合：从标量能量到合力

**概念桥梁。** §2 定义两个**标量**能量项；§3.1–§3.2 给出各项的**力**。本节说明：在一个顶点 $v$ 处，它所携带的若干标量偏差，如何汇成**一个合力方向**——以及为何这**不是**几何版「三角形质心」那一步。

**标量偏差本身没有方向。** 无论 $\delta_{c,i}^{2}$（中心角）还是 $(\alpha-\tau)^{2}$（内角），偏差只是一个数。力的方向不来自这个数本身，而来自这个数**对顶点位置的空间梯度**：

$$\mathbf F_{v}=-\nabla_{\mathbf r_v}\big(\text{含 }\mathbf r_v\text{ 的标量项}\big)
=-\frac{\partial(\text{偏差}^{2})}{\partial(\text{偏差})}\cdot\frac{\partial(\text{偏差})}{\partial\mathbf r_v}
=-2(\text{偏差})\cdot\Big(\frac{\partial\text{偏差}}{\partial\mathbf r_v}\Big).$$

第一因子（大小与符号）是标量；第二因子（一个几何向量）才是方向。

**中心角贡献 —— 切向，大小 $\propto|\delta|/\rho$。** 对局部索引为 $i$ 的细胞 $c$，

$$\psi=\operatorname{azimuth}(\mathbf r_v-\mathbf O_c),\qquad
\hat{\mathbf u}=\frac{\mathbf r_v-\mathbf O_c}{\rho},\qquad \rho=|\mathbf r_v-\mathbf O_c|.$$

方位角的梯度是一个**切向向量**：

$$\nabla_{\mathbf r_v}\psi=\frac{1}{\rho}\,J\hat{\mathbf u},\qquad J=\text{+90° 旋转}.$$

故细胞 $c$ 作用在 $v$ 上的中心角力（暂不计细胞共享项，见 §3.1）为

$$\mathbf F^{(c)}_{v}=-2\kappa_c\,\delta_{c,i}\,\frac{1}{\rho}\,J\hat{\mathbf u}\;-\;\mathbf F^{\mathrm{shared}}_c.$$

由式可见四点：
- **方向为切向**（⊥ 质心连线），而非指向质心 $\mathbf O_c$。
- **大小 $\propto|\delta_{c,i}|$**：偏差越大力越强；偏差变号力反向。
- **大小 $\propto 1/\rho$**：$v$ 离质心越远，方位角对平移越不敏感，力越弱。
- **径向力 = 0**：沿质心连线滑动 $v$ 不改变 $\psi$，故 $E_{\mathrm{cen}}$ 不驱动细胞**大小**——它是纯形状力。

（细胞共享项 $-2\kappa_c\mathbf{shared}_c$ 是细胞级常量，对细胞内各顶点相同；它负责平移零模与细胞刚性，但不改变上述每力**切向**的本质。）

**内角贡献 —— 沿张角方向。** 对 $v$ 处每个角（每个 (顶点, 细胞) 对），设其连向前后邻居的两条臂为 $\mathbf a,\mathbf b$，

$$\mathbf g_v=\frac{J\hat{\mathbf a}}{|\mathbf a|}-\frac{J\hat{\mathbf b}}{|\mathbf b|},
\qquad \mathbf F^{\mathrm{corner}}_v=-2\kappa_a\,(\alpha-\tau_v)\,\mathbf g_v.$$

方向 $\mathbf g_v$ 沿角平分线（张开/收紧该角）；大小 $\propto|\alpha-\tau_v|$。同一角还通过 §3.2.2 与它的两个相邻顶点耦合，故每个内角贡献被散到三个顶点。

**合力 = 矢量求和，无三角形、无目标点。** 对由细胞 $c_1,c_2,c_3$ 共享、局部索引分别为 $i_1,i_2,i_3$、带三个内角 $\alpha_1,\alpha_2,\alpha_3$ 的度 3 内部顶点，

$$\mathbf F_v=-\sum_{k=1}^{3}\Big[\,2\kappa_c\big(\delta_{c_k,i_k}\tfrac{1}{\rho_k}J\hat{\mathbf u}_k-\mathbf{shared}_{c_k}\big)\Big]
\;-\;\sum_{k=1}^{3}\Big[\,2\kappa_a(\alpha_k-\tfrac{2\pi}{3})\,\mathbf g^{(k)}_v\Big]
\;+\;\text{（邻居耦合项）}.$$

所有中心角与内角梯度向量均通过 `np.add.at` 累加到共享顶点上（energy.py）。**三条射线 / 三个偏差从不拼成三角形或绝对目标位置**——那个几何构造正是能量版丢弃的部分（§2.2）。取而代之的是：三个标量能量各自求导，其梯度向量矢量求和。

**与数学版对照。** 数学版拟合 $n_c$ 条最优射线后，对每个顶点取其 3 个所属细胞的 3 条射线、围成三角形、把顶点移到三角形质心——一个绝对目标点。能量版用以下方式取代：每条射线给出一个标量残差 $\delta^{2}$，每个残差形成一个梯度向量，三个向量求和成合力；顶点沿 $-\mathbf F_v$ 移动（再经线搜索）。由于两版优化的**目标不同**（显式目标点 vs. 单一标量 $E$ 的合力为零），其最优解**不会完全重合**——这也正是实测「$\lambda=1$ 大致复现数学版、但不完全一样」的原因。

***

## 4. 动力学与移动

### 4.1 过阻尼运动

$\zeta\frac{d\mathbf r_v}{dt}=\mathbf F_v\quad\Longrightarrow\quad \mathbf v_v=\frac{1}{\zeta}\mathbf F_v$

前向 Euler 离散（时间步 $\Delta t$）：
$\boxed{\ \mathbf r_v\leftarrow\mathbf r_v+\text{step}\cdot\mathbf F_v\ },\qquad \text{step}=\mu L_0^{2}\;(\text{L0 归一化，默认}),\quad \mu\equiv\frac{\Delta t}{\zeta}$

- $\mu$ = **迁移率**（步长常数），是**唯一**需要给定的时间尺度参数。
- 步长采用 **L0 归一化**：令 $L_0$ 为组织中位数边长，则 $\text{step}=\mu L_0^{2}$，每轮位移 $\approx\mu L_0$——$\mu$ 成为尺度无关的无量纲步长，跨不同规模网络无需调整（默认开启）。

### 4.2 移动方向

- **内部顶点**：自由 2D，沿 $\mathbf F_v$ 方向，$\Delta\mathbf r_v=\mu\,\mathbf F_v$。
- **边缘顶点**：只能沿边界滑动。它有 **2 条边界边**（§1），两条边切向不同、在顶点处形成折角，**切向不是唯一确定的**，须按「下降最多」选取：
  $\hat{\mathbf t}_v=\arg\max_{\hat{\mathbf t}\in\{\hat{\mathbf t}_1,\ \hat{\mathbf t}_2\}}\ \bigl\lvert\mathbf F_v\cdot\hat{\mathbf t}\bigr\rvert,\qquad \Delta\mathbf r_v=\mu\,(\mathbf F_v\cdot\hat{\mathbf t}_v)\,\hat{\mathbf t}_v$
  其中 $\hat{\mathbf t}_1,\hat{\mathbf t}_2$ 为该顶点两条边界边的单位方向（指向各自的外侧邻居）。
  - 含义：把总力投影到与它**最平行**的那条边界边上，沿该边滑动；
  - 若两条边的投影都近于零（力几乎垂直于边界），该顶点本轮**不动**；
  - 这是**动力学约束**（限制位移方向），不是能量项——能量照常按 §2 计算，只在施力时投影。

### 4.3 步长自适应（线搜索，串行逐顶点）

采用**串行降序**更新（§5）：自由顶点按 $|\mathbf F_v|$ **降序**排序后逐个移动（Gauss–Seidel，后续顶点用已更新坐标）。**每个顶点各自做线搜索**：从该顶点出发的位移先按 §4.1 的 $\text{step}$ 走，若使

1. 总能量上升（$E_{\mathrm{new}}>E_{\mathrm{old}}$，仅比较该顶点入射细胞的局部能量），或
2. 该顶点**入射的某一个**「本轮开始时是凸的」细胞变成非凸

则将步长折半（最多 60 次）重算该顶点的单步，直到满足接受条件；全部折半仍失败则该顶点本轮**不动**。

凸性守卫被**局部化到正在移动的那个顶点**——某一脆弱细胞只约束它自己所在的顶点，不会像同步方案那样因一个细胞翻非凸而否决整网整步（同步方案下整步被整体拒绝会令大随机网冻结）。起始即为非凸的细胞**豁免**（不要求它变凸，角度能量会自然把它拉回来）。

### 4.4 速度衰减

力的大小随逼近平衡自然减小（$\mathbf F\propto$ 偏差）→ 位移自动变小 → 过阻尼弛豫，无需额外降温调度。

***

## 5. 算法流程（单轮）

1. **建顶点索引**：把各细胞持有的重复坐标按键值聚类成唯一顶点；记录每个顶点的度（共享它的细胞数）、在各细胞中的引用；对边缘顶点**记录其 2 条边界边的方向** $\hat{\mathbf t}_1,\hat{\mathbf t}_2$（供 §4.2 选切向）。
2. **逐细胞算中心角量**：$\mathbf O_c$（顶点质心）→ 提升方向角 $\psi^{\uparrow}$ → 相位 $\varphi_c$ → 偏差 $\delta_{c,i}$；同时算共享项 $M_c$。
3. **累加中心角力**：按 §3.1 直接把每个细胞的 $\mathbf F^{(c)}_v$（已含 $\kappa_c$，是**力**不是梯度）累加到其各顶点。
4. **逐内角算角度量与权重**：算每个角（顶点, 细胞）的 $\alpha_{v,c}=\psi^{\uparrow}_a-\psi^{\uparrow}_b$、按顶点度数取 $\tau_v$、算标量权重 $w_{v,c}=2\kappa_a(\alpha_{v,c}-\tau_v)$。
   再**逐细胞沿 CCW 遍历其全部边**（含边缘细胞的边界边；内边被两细胞各遍历一次、都要累加），按 §3.2.2 把 $\mathbf g_c(X\!\to\!Y)=(w_{X,c}-w_{Y,c})J\hat{\mathbf e}_{X\to Y}/l$ 成对散到两端点——角顶项与臂端点项一次到位。
5. **装配总力** $\mathbf F_v=\mathbf F^{\mathrm{cen}}_v+\mathbf F^{\mathrm{ang}}_v$（§3.3）。
6. **按** **$|\mathbf F_v|$** **降序串行更新**：将自由顶点按受力大小 $|\mathbf F_v|$ **降序**排序；逐个顶点沿自身力 $\Delta\mathbf r_v=\text{step}\cdot\mathbf F_v$（边缘顶点取 §4.2 切向分量）移动，采用 Gauss–Seidel（后续顶点基于已更新的坐标计算），**每个顶点各自做线搜索**（§4.3）。
7. **刷新几何**：把新坐标写回各细胞的坐标副本。
   > 注：重算的**面积**仅供显示/导出，**能量模型不使用面积**；**质心**在下一轮步骤 2 由顶点坐标直接算出，无需在此缓存。
8. 重复 1–7，直到能量变化小于容差或达到最大轮数。

***

## 6. 与数学版的对应关系

| 数学版（几何规则）                         | 能量版（变分解释） |
| --------------------------------- | -------- |
| 最小二乘拟合 $n_c$ 条夹角 $2\pi/n_c$ 的最优射线 | 同一个最小二乘问题，**只取其解 $\varphi_c$**（闭式：算术平均）作为测量基准；射线的几何构造被丢弃（§2.2） |
| 三条最优射线构成三角形，取质心为目标点               | **不生成目标点**；改由对 $\mathbf r_v$ 的梯度下降替代该闭式步。连带后果：半径不再被 $E_{\mathrm{cen}}$ 直接约束（§2.2） |
| 移动若使 $\sum\alpha^2$ 增大则取消         | $E_{\mathrm{ang}}$ 的硬约束极限（$\lambda\to\infty$） |
| 边缘点沿小角对应边向中点移动                    | 在边界一维流形上对**完整的边缘内角能量** $\tfrac12(\beta_1-\beta_2)^2+\tfrac12(\beta_1+\beta_2-\pi)^2$ 做最速下降，并受切向投影约束（§4.2，两条边界边中取投影最大者） |
| 按退火距离**降序**、逐个串行移动                | **串行降序**（按 $|\mathbf F_v|$ 降序的 Gauss–Seidel 串行更新）：顶点按受力大小降序逐个移动，每个顶点各自做线搜索、凸性守卫局部化到其入射细胞。这与几何版「退火距离降序」语义对齐；与同步方案不同——同步时任一起始凸细胞被翻成非凸会整体否决整步、令大随机网冻结，串行则不会。 |

**总关系**：数学版 = 对 $E_{\mathrm{cen}}$ 交替最小化（先闭式定相位、再解析定目标点）+ 把 $E_{\mathrm{ang}}$ 当硬约束；能量版 = 对 $E_{\mathrm{cen}}+\lambda E_{\mathrm{ang}}$ **串行降序（Gauss–Seidel）软加权最小化**。**二者属同一族，不声称完全等价。**

### 6.1 内角约束的「开/关」→ 连续旋钮 $\lambda$

$\lambda\equiv\frac{\kappa_a}{\kappa_c}$

| 数学版开关   | 能量版 $\lambda$                               |
| ------- | ------------------------------------------- |
| 角度约束「关」 | $\lambda=0$（$\kappa_a=0$，角度力关闭）             |
| 角度约束「开」 | $\lambda\to\infty$（任何增大角偏差的移动使能量暴涨，被梯度下降拒绝） |
| ——      | **中间** **$\lambda$**：数学版到不了的新区间             |

**$\lambda$** **扫描**：固定 $\kappa_c=1$，让 $\kappa_a$（即 $\lambda$）沿一串取值（如 $0,\,0.2,\,0.5,\,1,\,2,\,5,\,10,\,50$）扫描，每个 $\lambda$ 弛豫到稳态并记录统计量。两个端点应分别复现数学版的「关 / 开」。

***

## 7. 模型的内在性质与边界（必须写明的限定）

1. **能量无量纲、严格尺度不变**：$E_{\mathrm{cen}},E_{\mathrm{ang}}$ 只依赖角度 → $E(\lambda\mathbf r)=E(\mathbf r)$。
2. **力随尺度衰减**：$|\mathbf F|\sim\kappa\,\delta/\rho$，即力的大小与组织长度尺度成反比（$\propto 1/\text{length}$）。
3. **三个零模**：能量对平移、旋转、均匀膨胀三类整体变换**均不敏感**，均无恢复力。

| 变换                                       | 不变性依据                             | 等价恒等式                                                      |
| ---------------------------------------- | --------------------------------- | ---------------------------------------------------------- |
| 平移 $\mathbf r_v\to\mathbf r_v+\mathbf c$ | $E$ 只依赖**相对**位置                   | $\sum_v\mathbf F_v=\mathbf 0$（§3.1 逐细胞求和为零；§3.2 逐角三顶点求和为零） |
| 旋转 $\mathbf r_v\to R\mathbf r_v$         | $\varphi_c$ 吸收整体旋转（§2.2）；内角本身旋转不变 | $\sum_v\mathbf r_v\times\mathbf F_v=0$（无净力矩）               |
| 膨胀 $\mathbf r_v\to\lambda\mathbf r_v$    | §7.1 的**尺度不变性**                   | $\sum_v\mathbf r_v\cdot\mathbf F_v=0$（膨胀不做功）               |

**膨胀零模的证明**：对 $E(\lambda\mathbf r)=E(\mathbf r)$ 在 $\lambda=1$ 处求导，
$0=\frac{d}{d\lambda}E(\lambda\mathbf r)\Big|_{\lambda=1}=\sum_v\nabla_{\mathbf r_v}E\cdot\mathbf r_v=-\sum_v\mathbf F_v\cdot\mathbf r_v$

> ⚠️ **易错**：「力求和为零」证明的是**平移**零模，**不能**用来论证膨胀零模——后者必须用上面的位力（virial）型恒等式 $\sum_v\mathbf F_v\cdot\mathbf r_v=0$。

- 推论：**组织的绝对尺寸不由本能量决定**，由初始条件与边界决定。（经典 vertex model 的尺寸由 $\Lambda$ 与 $\Gamma$ 的平衡定死；本模型无线张力，故无此机制。）
- 推论：**无内在长度尺度**，模型只定形状不定大小。
- 推论：结合 §7.2，大组织力更弱 → **绝对位移更慢**。用迁移率 $\mu$ 的写法下这是模型的真实性质，非数值缺陷。

1. **不含线张力 → 不粗化**：纯对称能量只会弛豫到一个受挫极小，不会驱动 T1 拓扑跃迁与 von Neumann–Mullins 粗化。这是对称松弛自身的限定。
2. **质心选择的影响**（见 §8.1）：顶点质心下方为零和、梯度公式精确；面积质心下不再严格成立。

***

## 8. 参数总表

### 8.1 模型与物理参数

| 符号                                 | 名称       | 定义 / 取值                                                   | 意义                                                        |
| ---------------------------------- | -------- | --------------------------------------------------------- | --------------------------------------------------------- |
| $E$                                | 总能量      | $\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$       | 被最小化的标量目标（无量纲）                                            |
| $E_{\mathrm{cen}}$                 | 中心角能量    | $\sum\delta^2$                                            | 细胞层级：方向等间隔性                                               |
| $E_{\mathrm{ang}}$                 | 内角能量     | $\sum(\alpha-\tau)^2$                                     | 顶点层级：角度规整性                                                |
| $\kappa_c$                         | 中心角刚度    | $>0$，常归一为 1                                               | 中心角项的权重                                                   |
| $\kappa_a$                         | 内角刚度     | $\ge 0$                                                   | 内角项的权重                                                    |
| $\lambda$                          | 相对角权重    | $\kappa_a/\kappa_c$                                       | **连续旋钮**；0 =约束关，$\infty$=约束开                              |
| $\mathbf O_c$                      | 细胞中心     | 见下                                                        | 测中心角的极点                                                   |
| $n_c$                              | 细胞顶点数    | 整数                                                        | 决定等间隔角 $2\pi/n_c$                                         |
| $\psi^{\uparrow}_{c,i}$            | 提升方向角    | $\arg(\mathbf r_{c,i}-\mathbf O_c)$，循环展开                  | 中心角的基本观测量                                                 |
| $\varphi_c$                        | 最佳拟合相位   | $\frac1{n_c}\sum(\psi^\uparrow-2\pi k/n_c)$               | **测量基准 / 规范固定**（吸收整体旋转），**非**目标构型；算术平均                    |
| $\delta_{c,i}$                     | 中心角偏差    | $\psi^\uparrow_{c,i}-\varphi_c-2\pi i/n_c$                | $E_{\mathrm{cen}}$ 的基本量；$\sum_i\delta\equiv0$             |
| $\rho_{c,i},\hat{\mathbf u}_{c,i}$ | 半径、径向单位矢 | $\lvert\mathbf u\rvert,\ \mathbf u/\lvert\mathbf u\rvert$ | 力按 $1/\rho$ 衰减                                            |
| $J$                                | +90° 旋转  | $J(x,y)=(-y,x)$                                           | 把径向转成切向（方位力）                                              |
| $M_c$                              | 细胞共享平均项  | $\frac1{n_c}\sum\delta_j J\hat{\mathbf u}_j/\rho_j$       | **顶点绑定项**；细胞全顶点相同                                         |
| $\tau$                             | 目标内角     | 内部 $2\pi/3$；边缘 $\pi/2$                                    | 泡沫力平衡导出的边界条件                                              |
| $\alpha,\beta$                     | 内角       | $\psi^{\uparrow}_a-\psi^{\uparrow}_b$，**顺序不可颠倒**（§2.3）    | $E_{\mathrm{ang}}$ 的基本量；$\mathbf a$=指前驱臂，$\mathbf b$=指后继臂 |
| $\mathbf F_v$                      | 顶点力      | $-\nabla_{\mathbf r_v}E$                                  | 驱动顶点运动的量                                                  |
| $\mu$                              | 迁移率      | $\Delta t/\zeta>0$                                        | **步长常数**（唯一时间尺度参数）                                        |
| $\zeta$                            | 阻尼系数     | 材料摩擦（本模型无量纲，$\zeta$ 无独立物理标定）                              | 与 $\Delta t$ 合并为 $\mu$                                    |
| $\hat{\mathbf t}$                  | 边界切向     | 边缘顶点关联边界边的单位方向                                            | 边缘顶点只能沿此方向滑动                                              |

### 8.2 关于质心的选择（已定：顶点质心）

| 方案           | 梯度                | 力零和（纯方位、膨胀不做功） | 与原数学版            |
| ------------ | ----------------- | -------------- | ---------------- |
| **顶点质心（选定）** | 公式 §3.1 精确成立      | 严格成立           | 原版用面积质心，模型改用顶点质心 |
| 面积质心         | 公式不严格（需额外质心运动修正项） | 仅为近似           | 与原版一致            |

**选定理由**：梯度精确 +「中心角力是纯方位力、均匀膨胀不做功」这条性质干净成立；与原版的差异只是 $O(\sigma/R)$ 的几何涨落。

### 8.3 算法控制参数（待定）

| 参数               | 含义                | 待定 |
| ---------------- | ----------------- | -- |
| $\mu$            | 迁移率初值 / 是否随轮次调度   | 待定 |
| 收敛容差             | 能量变化阈值            | 待定 |
| 最大轮数             | 单轮迭代上限            | 待定 |
| 线搜索最大折半次数        | 步长保护上限            | 待定 |
| $\lambda$ 扫描取值序列 | 扫描用 $\kappa_a$ 列表 | 待定 |

***

## 9. 已定决策清单

1. 能量只含两项：$E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$，**无线张力** **$\Lambda$**。
2. 质心取**顶点质心**。
3. $\varphi_c$ 取**算术平均**（非圆量平均）。
4. 内部目标角 $2\pi/3$、边缘目标角 $\pi/2$，二者均为**泡沫力平衡导出**，非经验常数。
5. 移动 = **过阻尼梯度下降**（前向 Euler），$\Delta\mathbf r=\text{step}\cdot\mathbf F$；**引入** **$L_0$** **中位边长归一化**：$\text{step}=\mu L_0^{2}$，每轮位移 $\approx\mu L_0$（默认开启，可在 `PhysicalAnnealerParams.normalize_L0=False` 时关）。
6. 顶点按 $|\mathbf F_v|$ **降序的 Gauss–Seidel 串行更新**（默认开启）；凸性守卫局部化到入射细胞。同步 Jacobi 更新仍为可选模式（设 `serial_descending=False`），但易致大随机网冻结，不推荐。
7. 边缘顶点**沿边界切向滑动**（切向投影）。
8. 线搜索：**串行方案下每个顶点各自做**；凸性守卫局部化到该顶点入射的细胞；起始即非凸的细胞豁免（设 `serial_descending=False` 时退化为全局线搜索）。
9. 内角约束从二元开关升级为**连续旋钮** **$\lambda=\kappa_a/\kappa_c$**。
10. $\varphi_c$ 与「最优射线」同解，但**只取相位**；其作用为**规范固定**（吸收整体旋转、保证能量旋转不变性），**不是**目标构型。射线几何与目标点均不构造。
11. $E_{\mathrm{cen}}$ **只约束方向、不约束半径**（中心角力为纯方位力）；细胞大小由 $E_{\mathrm{ang}}$ 与初始条件决定。
12. 内角梯度**必须含臂端点项**（§3.2）；其作用范围为**一跳邻域**，固定 9 项（内部）/ 6 项（边缘），**无连锁传递、无** **$O(V^2)$**。
13. 内角力**采用边（成对）形式实现**（§3.2.2）：逐细胞 CCW 遍历，每条边按 $(w_{X,c}-w_{Y,c})J\hat{\mathbf e}/l$ 成对散到两端点。复杂度与角顶版同阶，且 $\sum_v\mathbf F_v=\mathbf 0$ 显式成立。
14. 内角定义为 $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$（**从出边方向逆时针扫到入边反向**）。§3.2 梯度表与 §3.2.2 边形式均以此为唯一基准。
15. $E_{\mathrm{cen}}$、$E_{\mathrm{ang}}$ 一律为**不含刚度的裸和**；$\kappa_c,\kappa_a$ 只在总能量 $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$ 中出现一次。
16. 边缘顶点有 **2 条边界边**；滑动方向取两条边中 $\lvert\mathbf F\cdot\hat{\mathbf t}\rvert$ 较大者的切向（§4.2）。
17. 模型有**三个零模**：平移、旋转、膨胀。膨胀零模由位力恒等式 $\sum_v\mathbf F_v\cdot\mathbf r_v=0$ 给出，**不能**用「力求和为零」论证（那只证明平移零模）。
18. $\lambda$ **默认 = 1**；支持**动态调度**（同伦/连续法）：勾选后 $\lambda$ 从 `lam_start` 线性衰减到 `lam_end`（默认 $0\to1$），每轮重算——先以较大 $\lambda$ 把顶点推向 120° 流形预热，再降到 1 精细收敛。静态 $\lambda$ 模式不勾选调度。

***

## 10. 修订记录

- 修正：中心角力公式中 $\delta_{c,i}$ **只乘第一项**，不乘整个括号。
  正确：$\mathbf F^{(c)}=-2\kappa_c\left(\delta_{c,i}\dfrac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)$
  错误写法（已废弃）：$\mathbf F^{(c)}=-2\kappa_c\,\delta_{c,i}\left(\dfrac{J\hat{\mathbf u}_{c,i}}{\rho_{c,i}}-M_c\right)$
  依据：由顶点质心 $+\sum_i\delta_{c,i}\equiv0$ 严格推出；前者满足细胞力零和，后者不满足。
- 修正：步长去掉中位边长 $L_0$，改为纯迁移率形式 $\Delta\mathbf r_v=\mu\mathbf F_v$。
- 补充：内角梯度包含**臂端点贡献**（§3.2），只算角顶不完整。
- 澄清（§2.2）：$\varphi_c$ 与「最优射线」**同解但只取相位**，身份是\*\*规范固定（gauge fixing）\*\*而非目标构型——它吸收整体旋转自由度，保证 $E_{\mathrm{cen}}$ 只度量形状不对称、不度量朝向。若固定 $\varphi_c\equiv0$，整体旋转一个完美多边形会产生虚假能量。
- 澄清（§2.2）：$E_{\mathrm{cen}}$ **只约束方向、不约束半径**。数学版的「射线+三角形质心」同时定了方向与位置；能量版丢弃目标点后，中心角力恒为垂直于径向的纯方位力，不做径向功，故细胞大小不由 $E_{\mathrm{cen}}$ 决定。此为 §7.3「均匀膨胀是零模」的等价表述。
- 澄清（§3.2.1）：内角梯度**不存在连锁传递**。梯度是一阶的，只累加显式含 $\mathbf r_v$ 的项；「动 $v$ 会改邻居的角」是结果，已被一跳邻域内固定 9 项（内部）/ 6 项（边缘）精确记录，更远效应由迭代承担。复杂度每顶点 $O(1)$、全局 $O(V)$。
- 新增（§3.2.2）：内角梯度的**边（成对）等价形式** $\mathbf g_c(X\!\to\!Y)=(w_{X,c}-w_{Y,c})J\hat{\mathbf e}_{X\to Y}/l$，逐细胞 CCW 遍历一次即可，角顶与臂端点贡献自动合并。附：复杂度论证与 $\sum_v\mathbf F_v=\mathbf 0$ 自动成立。
- 新增（§3.2.3）：臂端点项**不可省略**的三条理由（变分原理失效 / 虚假整体漂移 / 线搜索卡死）与「同阶小量、非高阶修正」的量级论证；退化方案（分块局部规则）明确不采用。
- **修正（符号错误，严重）**：内角定义由 $\alpha=\psi^{\uparrow}_b-\psi^{\uparrow}_a$ 改为 $\alpha=\psi^{\uparrow}_a-\psi^{\uparrow}_b$。
  原写法给出的是 $-\alpha$（内角的相反数），导致 §3.2 三行梯度表**全部反号**，与 §3.2.2 的边形式自相矛盾（两者差一个整体负号）。
  依据：CCW 顶点序下 $\mathbf a$=指向前驱的臂、$\mathbf b$=指向后继的臂，$\psi_a=\theta_{\mathrm{in}}+\pi$、$\psi_b=\theta_{\mathrm{out}}$，故 $\alpha=\pi-(\theta_{\mathrm{out}}-\theta_{\mathrm{in}})=\psi_a-\psi_b$。正方形的手算校验：$\nabla_{\mathbf r_i}\alpha_i=\dfrac{J\hat{\mathbf e}_{\mathrm{in}}}{l_{\mathrm{in}}}+\dfrac{J\hat{\mathbf e}_{\mathrm{out}}}{l_{\mathrm{out}}}$，与顶点沿 $(-1,1)$ 方向移动时内角增大的直接计算一致。
  连带修正：§3.2 梯度表三行符号、§2.3 内角定义、§8.1 参数表内角定义。
  §3.2.2 边形式的**结论不变**（经独立推导与符号自检：$w_X=w_Y$ 时该边贡献为零）。
- **修正（刚度重复计入）**：§2.2/§2.3 盒式公式中的 $E_{\mathrm{cen}}$、$E_{\mathrm{ang}}$ 原写作 $\kappa\sum(\cdot)^2$，与 §2.1 总能量 $E=\kappa_cE_{\mathrm{cen}}+\kappa_aE_{\mathrm{ang}}$ 叠加后刚度被乘两次（出现 $\kappa^2$）。现统一改为**裸和** $E_{\mathrm{cen}}=\sum\delta^2$、$E_{\mathrm{ang}}=\sum(\alpha-\tau)^2$，刚度只在总能量中出现一次，与 §3 的梯度权重 $2\kappa_c\delta$、$2\kappa_a(\alpha-\tau)$ 一致。
- **修正（§3.3 漏刚度）**：$\mathbf F^{\mathrm{ang}}_v$ 原写作 $-\nabla_{\mathbf r_v}E_{\mathrm{ang}}$，缺 $\kappa_a$；现改为 $-\kappa_a\nabla_{\mathbf r_v}E_{\mathrm{ang}}=-\mathbf G_v$，并说明 §3.2.2 累加得到的是 $\mathbf G_v$（$w$ 已含 $\kappa_a$）。
- **修正（§4.2 切向未定义）**：边缘顶点实有 **2 条**边界边（§1 已补充：2 边界边 + 1 内边），切向在折角处不唯一。现规定取 $\lvert\mathbf F\cdot\hat{\mathbf t}\rvert$ 较大者，并明确这是**动力学约束而非能量项**。
- **修正（§7.3 零模论证错误）**：原文用「力求和为零」论证膨胀零模，但该恒等式证明的是**平移**零模。现拆分为平移 / 旋转 / 膨胀三个零模，膨胀零模改用位力恒等式 $\sum_v\mathbf F_v\cdot\mathbf r_v=0$（由尺度不变性对 $\lambda$ 求导得到）证明；旋转零模为新增。
- **修正（§6 对应表）**：边缘行原写「对 $(\beta_1-\beta_2)^2$ 最速下降」，漏掉边界弯曲项 $\tfrac12(\beta_1+\beta_2-\pi)^2$；现补为完整的边缘内角能量。另新增「退火距离降序串行 ↔ 更新顺序」一行（该行列最初稿曾写作「同步更新」，后于本轮修订修正为按 $|\mathbf F_v|$ 降序的串行 Gauss–Seidel 更新）。
- **澄清（§3.2.2 遍历范围）**：内边被所属两个细胞**各遍历一次且都要累加**（权重不同，不可去重）；边缘细胞的边界边同样参与遍历。另补充等价性的逐步推导与符号自检。
- **修正（§3.2 边缘顶点项数说明）**：项数 6 正确，但原解释「2 条边各属于 2 个细胞」有误；实为内边属 2 细胞（2 项）+ 2 条边界边各属 1 细胞（各 1 项）= 4 项。
- **澄清（§5 流程）**：步骤 1 记录边缘顶点的 **2 条**边界边方向；步骤 3 累加的是**力**而非梯度；步骤 7 的重算**面积不参与能量**（仅供显示/导出），质心由下一轮步骤 2 直接算出。
- **本轮修订（与代码 / introduction 对齐）**：① 移动顺序由「同步更新」改为**按** **$|\mathbf F_v|$** **降序的 Gauss–Seidel 串行更新**（默认），凸性守卫局部化到入射细胞（§4.3、§5 步骤 6、§6 对应表、§9 决策 #6）；② 引入 **L0 中位边长归一化**（$\text{step}=\mu L_0^{2}$，位移 $\approx\mu L_0$，默认开启，`normalize_L0=False` 可关），撤销原「不引入 L0」的写法（§4.1、§9 决策 #5）；③ 线搜索由「整步全局否决」改为**逐顶点各自做**、凸性守卫局部化（§4.3、§9 决策 #8）；④ 补充 $\lambda$ **默认 = 1** 与**动态调度** $0\to1$（同伦法，每轮重算）写入 §9 决策 #18；⑤ 同步 Jacobi 保留为可选模式（`serial_descending=False`）。

