# Cell Relaxation Tool — 细胞退火工具

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.14+](https://img.shields.io/badge/python-3.14+-blue.svg)](https://www.python.org/downloads/)

A Python implementation of the **symmetric relaxation (energy version) algorithm** for entire 2D cellular networks. 

Based on:

- Xu K., Weng L., Wang Z., Lian Y., Huang B. (2026). *A symmetric relaxation method for entire two-dimensional cellular networks and its implications*. (arXiv: XuSR20260616)
- Xu K. (2021). *A geometry-based relaxation algorithm for equilibrating a trivalent polygonal network in two dimensions and its implications*. Philosophical Magazine, 101(14), 1632-1653.

***

## Algorithm Introduction (Energy Version)

> This tool implements the **energy-version** symmetric relaxation algorithm. The brief introduction below is given in both English and Chinese.

### English

**Energy-Version Symmetric Relaxation for 2D Cellular Network — A Brief Introduction**

The energy-version symmetric relaxation algorithm was produced by **Kai Xu** with the assistance of the **Hy4 large language model (LLM)**, converting the mathematical-version symmetric relaxation algorithm (Xu et al. 2026, arXiv:2606.18604) into an energy version; it relaxes a 2D cellular network toward physical equilibrium by minimizing a single scalar energy

$$E = E_{\text{cen}} + \lambda\,E_{\text{ang}},$$

where $E_{\text{cen}}$ is the center-angle (cell-barycenter) term and $E_{\text{ang}}$ the vertex interior-angle term, with $k_{\text{center}}$ fixed to 1 so that $\lambda = k_{\text{angle}}/k_{\text{center}}$ is the sole dimensionless knob linking the method to the geometric version.

The two terms encode symmetry at two distinct scales. The **center-angle term** $E_{\text{cen}}$ enforces *cell-level* symmetry: each cell's vertices are driven to be symmetrically arranged about the cell's own center (barycenter), so minimizing it makes every cell a self-symmetric — in the ideal limit, a regular — polygon. The **interior-angle term** $E_{\text{ang}}$ enforces *vertex-level* symmetry at every vertex: at each interior (degree-3) vertex shared by three cells, the three angles are driven toward the equal 120° configuration (the Plateau/foam rule), i.e. a symmetric three-fold junction; at each boundary (degree-2 marginal) vertex shared by two cells, the two interior angles are driven toward 90° — the foam free-boundary (Young) condition where the internal edge meets the boundary orthogonally — and the vertex is constrained to slide only along the boundary tangent rather than move freely. $E_{\text{ang}}$ is exactly the energy analogue of the geometric version's interior-angle-square-sum guard.

Each free vertex is updated by over-damped gradient descent, $r_v \leftarrow r_v + \mu\,F_v$ with $F_v = -\nabla_v E$. The step is **L0-normalized** (effective step $\mu L_0^2$, per-round displacement $\approx \mu L_0$), making $\mu$ a scale-invariant fraction of the median edge length that works unchanged across network sizes. Vertices are swept in **descending order of $|F_v|$** (serial-descending Gauss–Seidel), each with its own line search and a convexity guard localized to the cells it touches — this stops one fragile cell from vetoing the whole update, letting large random cellular networks relax instead of freezing.

The weight $\lambda$ tunes the balance between the two symmetries: $\lambda=0$ keeps only the center-angle symmetry, $\lambda=1$ balances the two terms equally, and large $\lambda$ (10–100) makes the angle-term symmetry dominant, approaching the hard interior-angle constraint in the $\lambda\to\infty$ limit.

### 中文

**能量版二维细胞网络对称退火方法 —— 简介**

能量版对称退火算法由 Kai Xu 借助 Hy4 大语言模型（LLM）将数学版对称退火算法（Xu et al. 2026, arXiv:2606.18604）转换为能量版得到；它通过最小化一个标量能量，将二维细胞网络退火到物理平衡态：

$$E = E_{\text{cen}} + \lambda\,E_{\text{ang}},$$

其中 $E_{\text{cen}}$ 为中心角（细胞重心）项，$E_{\text{ang}}$ 为顶点内角项；$k_{\text{center}}$ 固定为 1，因此 $\lambda = k_{\text{angle}}/k_{\text{center}}$ 是连接该方法与几何版的唯一无量纲旋钮。

两项分别在两个尺度上编码对称性。**中心角项** $E_{\text{cen}}$ 强制*细胞级*对称：驱动每个细胞的顶点围绕其自身中心（重心）对称排布，最小化它使每个细胞成为自对称——理想极限下为正则——多边形。**顶点内角项** $E_{\text{ang}}$ 在*每个顶点*上强制*顶点级*对称：在每个由三个细胞共享的内部（degree‑3）顶点处，三个内角被驱动趋向相等的 120° 构型（Plateau/泡沫规则），即对称的三叉汇点；在每个由两细胞共享的边界（degree‑2 边缘）顶点处，两个内角被驱动趋向 90°——即泡沫自由边界（Young）条件，内部边与边界正交——且该顶点被约束只能沿边界切向滑动，而非自由移动。$E_{\text{ang}}$ 正是几何版"内角平方和守卫"的能量对应物。

每个自由顶点通过过阻尼梯度下降更新：$r_v \leftarrow r_v + \mu\,F_v$，其中 $F_v = -\nabla_v E$。步长采用 **L0 归一化**（有效步长 $\mu L_0^2$，每轮位移 $\approx \mu L_0$），使 $\mu$ 成为中位数边长的一个尺度无关的分数，在不同规模的网络上均无需调整。顶点按 $|F_v|$ **降序**扫描（串行降序 Gauss–Seidel），每个顶点各自做线搜索，且凸性守卫局部化到其所属的细胞——这避免单个脆弱细胞否决整步更新，使大型随机细胞网络能够退火而非冻结。

权重 $\lambda$ 调节两种对称性之间的平衡：$\lambda=0$ 时仅保留中心角对称；$\lambda=1$ 时两项等权；$\lambda$ 较大（10–100）时内角对称主导，在 $\lambda\to\infty$ 极限下趋近硬内角约束。

***

## Features

### Voronoi Network Generation

Generate trimmed Voronoi networks:

- **Hexagon random-disordered Voronoi (n×n)** — Randomly perturb the coordinates of seed points in a regular hexagonal grid using the rhombus-based sampling method, with perturbation range confined within the hexagon, controlled by parameter *k* (0 = regular hexagons, 1 = high disorder).
- **Uniform Random Voronoi** — Based on uniformly random points.

### Symmetric Relaxation (Annealing)

Equilibrate a 2D cellular network toward physical equilibrium by minimizing a single scalar energy $E = E_{\text{cen}} + \lambda\,E_{\text{ang}}$:

- **Energy model** — The center-angle term $E_{\text{cen}}$ enforces *cell-level* symmetry (driving each cell toward self-symmetry / a regular polygon); the vertex interior-angle term $E_{\text{ang}}$ enforces *vertex-level* symmetry (120° at inner vertices, 90° orthogonal-to-boundary at marginal vertices).
- **Gradient descent** — Each free vertex is updated by over-damped gradient descent with an L0-normalized step (displacement ≈ μ·L0), scale-invariant; vertices are swept in descending order of $|F_v|$ (serial-descending Gauss–Seidel), each with its own line search and a local convexity guard so that one fragile cell cannot veto the whole update.
- **Configurable** — Relaxation factor μ (step size), center-angle energy weight λ (best at λ=1), marginal-vertex participation toggle, and dynamic λ scheduling (0→1).

### Ellipse Fitting

- **R-LMG** (via `conicfit`) — High-precision geometric fitting, preferred for triangles and quadrilaterals.
- **Algebraic Least Squares** (numpy SVD) — Fast fitting, preferred for polygons with 5+ sides.
- The area ratio of the ellipse to the n-gon must be between 1 and 3. If the preferred method fails, fit using algebraic least squares on real points + mirror points (2n total); if that fails, fit using R-LMG on real points + mirror points (2n total); if that still fails, fall back to the minimum bounding rectangle (MBR) on real points + mirror points (2n total).

### Data Export

Export results to Excel (`.xlsx`):

- **Ellipse & Cell perimeter and area** — Cell ID, Ellipse Centroid, Major Semi-axis, Minor Semi-axis, Inclination Angle of the Major Semi-axis, Sum of Neighbor Edge Numbers, Cell Edge Number, Cell Perimeter, Cell Area, Original Cell Layer, Current Cell Layer (layer definition: outermost boundary cells = layer 1, incrementing inward).
- **Edge Angle Data** — Cell ID, Edge Number, Sum of Neighbor Edge Numbers, Interior Angle, Adjacent Edge 1, Adjacent Edge 2, Original Cell Layer, Current Cell Layer (one row per vertex).
- **ME/MA Statistics** — Cell ID, Margin Edge Length (ME), Margin Angle 1 (MA1), Margin Angle 2 (MA2), Original Cell Layer, Current Cell Layer (boundary cells only, one row per cell).

### Visualization

- Real-time cell network display with color-coded layers.
- Ellipse overlay for fitted cells.
- Image export to EPS (vector) and PNG formats.

### Energy Observation Window

- Plots the energy curve in real time (accumulated per "Execute Annealing"), with annealing step on the x-axis and energy on the y-axis, breaking down into $E_{\text{total}}$, $E_{\text{cen}}$, and $\lambda E_{\text{ang}}$.
- Shows convergence status: Stable / Converging / Not decreasing.
- Supports configuring λ (center-angle energy weight) and its dynamic scheduling (0→1).

***

## Requirements

- **Python 3.14+** (managed by [uv](https://docs.astral.sh/uv/); the setup script installs Python 3.14.6 automatically)
- Python packages: see [requirements.txt](requirements.txt)
- **R** language environment (for R-LMG ellipse fitting via `rpy2`), optional but recommended
- R packages: `conicfit`, `sp`, `shotGroups`

***

## Installation

> Prerequisite: install [uv](https://docs.astral.sh/uv/getting-started/installation/) first (a fast Python package manager; works on Windows / macOS / Linux).

### Quick Start (one-command setup)

**Windows:**

```bash
setup.bat
```

**Linux / macOS:**

```bash
chmod +x setup.sh
./setup.sh
```

The script will automatically:

1. Check uv, install Python via `uv python install 3.14.6`, and create a `.venv`.
2. Install all Python packages via `uv pip install -r requirements.txt`.
3. Detect the bundled `R_Dist/` portable directory and configure it.
4. Fall back to system R if not found.
5. **Automatically download and install R** if still not found (Windows/macOS).
6. Install required R packages (`conicfit`, `sp`, `shotGroups`).

### Manual Setup

```bash
# 1. Install Python 3.14.6 and create a virtual environment with uv
uv python install 3.14.6
uv venv .venv --python 3.14.6
uv pip install -r requirements.txt

# 2. Install R and required R packages (optional, for R-LMG fitting)
#    Install R from https://www.r-project.org/, then:
#    Rscript -e "install.packages(c('conicfit', 'sp', 'shotGroups'))"

# 3. Run
uv run python energy_annealing_main.py
```

### Build EXE (PyInstaller)

Run [build_exe.bat](build_exe.bat) directly, or invoke manually:

```bash
uv pip install pyinstaller
uv run pyinstaller --onefile --windowed --name "CellAnnealing" ^
    --add-data "annealing;annealing" ^
    --add-data "cell;cell" ^
    --add-data "utillib;utillib" ^
    --add-data "initVoronoi.py;." ^
    --hidden-import "scipy.spatial" ^
    --hidden-import "scipy.optimize" ^
    --hidden-import "openpyxl" ^
    --hidden-import "pyenvelope" ^
    energy_annealing_main.py
# Output: dist/CellAnnealing.exe
```

> After building, place `R_Dist/` next to the EXE to enable R-LMG fitting.

***

## Project Structure

```
only_annealing/
├── energy_annealing_main.py      # Main entry point (Tkinter GUI)
├── requirements.txt            # Python dependencies
├── setup.bat                   # Windows one-click setup script (uv)
├── setup.sh                    # Linux/macOS one-click setup script (uv)
├── build_exe.bat               # PyInstaller build script (uv)
├── .gitignore
├── .gitattributes
├── LICENSE                     # MIT License
├── README_EN.md                # English documentation
├── README_CN.md                # Chinese documentation
├── ANNEALING_FLOW_CN.md        # Annealing algorithm flow documentation (Chinese)
├── PARAMS_CN_EN.md             # Parameter Chinese-English reference
│
├── initVoronoi.py              # Voronoi network initialization
│
├── annealing/                  # Relaxation (annealing) algorithms
│   ├── AnnealingGUI.py         # Annealer core (parameterized)
│   └── annealerUtil.py         # Annealing utility functions
│
├── cell/                       # Core cell data & statistics
│   └── CellData.py             # Cell data management
│
├── scripts/                    # Auxiliary scripts (data collection/testing)
│   ├── collect_data.py         # Data collection script
│   ├── collect_iteration_data.py  # Iteration data collection
│   ├── test_area.py            # Area calculation unit tests
│   └── DATA_COLLECTION_PLAN.md # Batch data collection plan
│
└── utillib/                    # Utility libraries
    ├── mylib.py                # Core data structures (Cell, Point, Line, etc.)
    ├── fittinglib.py           # Ellipse fitting (R-LMG + least squares)
    ├── exportUtils.py          # Excel export utilities
    ├── layerMarker.py          # Cell layer marking
    └── i18n.py                 # Internationalization (Chinese/English)
```

***

## R Environment

The `setup.bat` / `setup.sh` scripts handle R automatically with the following priority:

| Priority | Source | Detection Method |
|---|---|---|
| 1 (best) | Bundled `R_Dist/` | `R_Dist/bin/R` exists locally |
| 2 | System R | `R_HOME` environment variable or `R` in PATH |
| 3 (fallback) | Auto-download | Script downloads and installs R automatically |

If no R environment is found, ellipse fitting falls back to pure Python (numpy/scipy least squares), which remains fully functional but may produce slightly different results for triangles and quadrilaterals.

> **Note:** `R_Dist/` is excluded from git (see `.gitignore`). Run `setup.bat` after cloning to configure R, or copy `R_Dist/` from a complete distribution.

***

## Usage

> **Note:** When running the annealing program on a regular computer, the number of cells is recommended not to exceed 400, otherwise it will be very slow or crash; after closing the graphical visualization interface, the cell limit can be moderately increased.

```bash
uv run python energy_annealing_main.py
```

The GUI provides:

1. **Voronoi Initialization** — Configure network type and parameters, then generate.
2. **Annealing** — Set relaxation factor μ, center-angle energy weight λ (and its dynamic scheduling), and marginal-vertex participation, then run single/multiple iterations; the energy observation window accumulates the energy curve in real time.
3. **Energy Observation Window** — Display the energy curve and convergence status (Stable / Converging / Not decreasing) in real time.
4. **Ellipse Fitting** — Fit ellipses to all cells and overlay on display.
5. **Data Export** — Save geometric and topological parameters of all cells to three Excel files.
6. **Image Export** — Save current view as EPS or PNG via the matplotlib toolbar.

***

## Citation

If you use this code in academic work, please cite:

```
Xu K., Weng L., Wang Z., Lian Y., Huang B. (2026). A symmetric relaxation method for entire two-dimensional cellular networks and its implications. arXiv: XuSR20260616.
```

```
Xu K. (2021). A geometry-based relaxation algorithm for equilibrating a trivalent polygonal network in two dimensions and its implications. Philosophical Magazine, 101(14), 1632-1653.
```

***

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

***

## Contributors

- **Kai Xu (许凯)** 
- Email: <kaixu@jmu.edu.cn> / <kxu2013@gmail.com>
- ORCID: [0000-0002-1341-1525](https://orcid.org/0000-0002-1341-1525)
- Affiliation: Fisheries College, Jimei University, Xiamen, China

For questions, bug reports, or collaboration inquiries, please contact the corresponding author.
