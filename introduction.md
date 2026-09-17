# Energy-Version Symmetric Relaxation for 2D Cellular Network

## English

**Energy-Version Symmetric Relaxation for 2D Cellular Network — A Brief Introduction**

The energy-version symmetric relaxation algorithm was produced by **Kai Xu** with the assistance of the **Hy4 large language model (LLM)**, converting the mathematical-version symmetric relaxation algorithm (Xu et al. 2026, arXiv:2606.18604) into an energy version; it relaxes a 2D cellular network toward physical equilibrium by minimizing a single scalar energy

$$E = E_{\text{cen}} + \lambda\,E_{\text{ang}},$$

where $E_{\text{cen}}$ is the center-angle (cell-barycenter) term and $E_{\text{ang}}$ the vertex interior-angle term, with $k_{\text{center}}$ fixed to 1 so that $\lambda = k_{\text{angle}}/k_{\text{center}}$ is the sole dimensionless knob linking the method to the geometric version.

The two terms encode symmetry at two distinct scales. The **center-angle term** $E_{\text{cen}}$ enforces *cell-level* symmetry: each cell's vertices are driven to be symmetrically arranged about the cell's own center (barycenter), so minimizing it makes every cell a self-symmetric — in the ideal limit, a regular — polygon. The **interior-angle term** $E_{\text{ang}}$ enforces *vertex-level* symmetry at every vertex: at each interior (degree-3) vertex shared by three cells, the three angles are driven toward the equal 120° configuration (the Plateau/foam rule), i.e. a symmetric three-fold junction; at each boundary (degree-2 marginal) vertex shared by two cells, the two interior angles are driven toward 90° — the foam free-boundary (Young) condition where the internal edge meets the boundary orthogonally — and the vertex is constrained to slide only along the boundary tangent rather than move freely. $E_{\text{ang}}$ is exactly the energy analogue of the geometric version's interior-angle-square-sum guard.

Each free vertex is updated by over-damped gradient descent, $r_v \leftarrow r_v + \mu\,F_v$ with $F_v = -\nabla_v E$. The step is **L0-normalized** (effective step $\mu L_0^2$, per-round displacement $\approx \mu L_0$), making $\mu$ a scale-invariant fraction of the median edge length that works unchanged across network sizes. Vertices are swept in **descending order of $|F_v|$** (serial-descending Gauss–Seidel), each with its own line search and a convexity guard localized to the cells it touches — this stops one fragile cell from vetoing the whole update, letting large random cellular networks relax instead of freezing.

The weight $\lambda$ tunes the balance between the two symmetries: $\lambda=0$ keeps only the center-angle symmetry, $\lambda=1$ balances the two terms equally, and large $\lambda$ (10–100) makes the angle-term symmetry dominant, approaching the hard interior-angle constraint in the $\lambda\to\infty$ limit. In a sweep over $\lambda\in\{0,1,10,100\}$ on three initial cellular networks — slightly disordered ($k=0.2$), strongly disordered ($k=1$) hexagonal-grid Voronoi, and uniform-random Voronoi — with 400 seed points × 10 repeats, static $\lambda$, serial-descending update and $\mu=0.1$, the algorithm at **$\lambda=1$ relaxes best: the evolution of cell area, edge length, and interior angles closely matches that of the geometric version, which applies center-angle symmetry together with an interior-angle-square-sum guard**.

## 中文

**能量版二维细胞网络对称退火方法 —— 简介**

能量版对称退火算法由 Kai Xu 借助 Hy4 大语言模型（LLM）将数学版对称退火算法（Xu et al. 2026, arXiv:2606.18604）转换为能量版得到；它通过最小化一个标量能量，将二维细胞网络退火到物理平衡态：

$$E = E_{\text{cen}} + \lambda\,E_{\text{ang}},$$

其中 $E_{\text{cen}}$ 为中心角（细胞重心）项，$E_{\text{ang}}$ 为顶点内角项；$k_{\text{center}}$ 固定为 1，因此 $\lambda = k_{\text{angle}}/k_{\text{center}}$ 是连接该方法与几何版的唯一无量纲旋钮。

两项分别在两个尺度上编码对称性。**中心角项** $E_{\text{cen}}$ 强制*细胞级*对称：驱动每个细胞的顶点围绕其自身中心（重心）对称排布，最小化它使每个细胞成为自对称——理想极限下为正则——多边形。**顶点内角项** $E_{\text{ang}}$ 在*每个顶点*上强制*顶点级*对称：在每个由三个细胞共享的内部（degree‑3）顶点处，三个内角被驱动趋向相等的 120° 构型（Plateau/泡沫规则），即对称的三叉汇点；在每个由两细胞共享的边界（degree‑2 边缘）顶点处，两个内角被驱动趋向 90°——即泡沫自由边界（Young）条件，内部边与边界正交——且该顶点被约束只能沿边界切向滑动，而非自由移动。$E_{\text{ang}}$ 正是几何版"内角平方和守卫"的能量对应物。

每个自由顶点通过过阻尼梯度下降更新：$r_v \leftarrow r_v + \mu\,F_v$，其中 $F_v = -\nabla_v E$。步长采用 **L0 归一化**（有效步长 $\mu L_0^2$，每轮位移 $\approx \mu L_0$），使 $\mu$ 成为中位数边长的一个尺度无关的分数，在不同规模的网络上均无需调整。顶点按 $|F_v|$ **降序**扫描（串行降序 Gauss–Seidel），每个顶点各自做线搜索，且凸性守卫局部化到其所属的细胞——这避免单个脆弱细胞否决整步更新，使大型随机细胞网络能够退火而非冻结。

权重 $\lambda$ 调节两种对称性之间的平衡：$\lambda=0$ 时仅保留中心角对称；$\lambda=1$ 时两项等权；$\lambda$ 较大（10–100）时内角对称主导，在 $\lambda\to\infty$ 极限下趋近硬内角约束。在 $\lambda\in\{0,1,10,100\}$ 的取值扫描中，对三类初始细胞网络——轻微扰动（$k=0.2$）、强扰动（$k=1$）的六边形网格 Voronoi，以及均匀随机 Voronoi——各取 400 个种子点 × 10 次重复，采用静态 $\lambda$、串行降序更新、$\mu=0.1$，算法在 **$\lambda=1$ 时退火效果最好：细胞面积、边长与内角的演化，与几何版（应用中心角对称 + 内角平方和守卫）的结果高度相似**。
