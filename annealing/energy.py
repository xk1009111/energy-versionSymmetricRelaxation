# -*- coding: utf-8 -*-
"""
对称松弛 · 能量版（变分）—— 由数学版（几何规则）退火修改而来。

本模块从数学版细胞退火（原 ``annealerUtil.py`` / ``math_annealer.py`` 的几何
规则实现，已删除）修改而来：把「拟合最优射线 → 三角形质心目标点 → 移动」的
规则集，升级为对单个标量能量做（过阻尼）梯度下降的变分原理。全部公式与符号
以 ``SymRelax_Energy_Model_Spec.md``（完整规范）与
``SymRelax_Energy_Force_Scheme.md``（力计算方案）为准。能量定义（规范 §2，
E_cen / E_ang 均为不含刚度的裸和，κ 只在总能量中出现一次，决策 #15）：

    E = k_center * E_cen  +  k_angle * E_ang,   lambda = k_angle / k_center

    E_cen = sum_c sum_i  delta_{c,i}^2                    (规范 §2.2，中心角项)
        delta_{c,i} = psi^up_{c,i} - phi_c - 2*pi*i/n_c
        phi_c       = mean_i(psi^up_{c,i} - 2*pi*i/n_c)   算术平均；规范固定
                      (gauge fixing)：测量基准而非目标构型（决策 #10）

    E_ang = sum_{corners} (alpha - tau)^2                  (规范 §2.3，内角项)
        alpha = psi^up_a - psi^up_b   a=前驱臂, b=后继臂（顺序不可颠倒）
        tau = 2*pi/3  内部顶点（度 3，泡沫三膜力平衡 / Plateau 导出）
            = pi/2    边缘顶点（度 2，自由边界正交条件导出）

与数学版的对应（规范 §6）
--------------------------
    数学版「最小二乘拟合 n_c 条最优射线」
        <=> 同一最小二乘问题只取解 phi_c 作测量基准；射线几何不构造
    数学版「三条射线三角形取质心为目标点」
        <=> 不生成目标点，由梯度下降替代（半径不再被 E_cen 约束）
    数学版「移动若使 sum(alpha^2) 增大则取消」
        <=> E_ang 在 lambda -> +inf 的极限
    数学版「边缘点沿小角对应边向中点移动」
        <=> tau = pi/2 的完整边缘内角能量在边界切向上的最速下降
           （切向投影由 PhysicalAnnealer 施加，本模块只算能量与梯度）
    数学版「按退火距离降序串行移动」
        <=> 全体顶点同步更新（决策 #6）

关键正确性要点（同规范修订记录）
--------------------------------
1.  phi_c 取提升角（lifted / unwrapped，规范 §1）的算术平均，非圆量平均——
    只有算术平均给出 sum_i delta_i = 0，梯度公式（规范 §3.1）才精确成立。
    数值验证：算术平均梯度误差 ~1e-10，圆量平均误差 ~O(1)。

2.  质心 = 顶点质心（已定 A，规范 §8.2）：梯度公式精确，且中心角力为纯方位力
    （逐细胞求和为零 → 均匀膨胀不做功）。数学版用面积质心，二者差 O(sigma/R)。
    'area' 模式仅保留用于与原方法做数值对照。

3.  两项均为无量纲角度泛函 → 能量严格尺度不变（规范 §7.1），力 |F| ~
    kappa*delta/rho 随尺度衰减（规范 §7.2）；平移/旋转/膨胀三个零模（规范 §7.3）。

梯度公式（解析；中心差分校验 ~1e-10；规范 §3）
------------------------------------------------
记 J 为 +90° 旋转，J(x, y) = (-y, x)。内角按规范记号 alpha = psi^up_a - psi^up_b
（a = 前驱臂 r_prev - r_w，b = 后继臂 r_next - r_w）。注意本模块代码变量 a/b
与规范相互置换（代码 a_hat=后继臂、b_hat=前驱臂），梯度相应镜像、数值一致：

    grad_{r_w} alpha = -J a_hat / |a| + J b_hat / |b|     (角顶 w)
    grad_{r_p} alpha = +J a_hat / |a|                     (前臂端点 p)
    grad_{r_q} alpha = -J b_hat / |b|                     (后臂端点 q)

中心角项（u_i = (r_i - O) / rho_i）：

    d psi_{c,i} / d r_{c,k} = (delta_ik - 1/n_c) * J u_i / rho_i
    grad_{r_k} E_cen = 2 k_center * [ delta_k * J u_k / rho_k
                                    - (1/n_c) * sum_i delta_i * J u_i / rho_i ]

  括号内第二项即细胞共享平均项 M_c（规范 §3.1）：对细胞内每个顶点相同——
  「中心角把细胞所有顶点绑定」的变分表达。

数据模型
----------
``cells`` 为 ``utillib.mylib.Cell`` 列表。每个细胞各自持有一份顶点坐标副本
（``cell.points``，[x, y] 列表；规范 §1），共享边在两细胞中以重合坐标出现。
:func:`build_vertex_index` 将重合坐标聚类为唯一全局顶点，并把逐细胞梯度累加
回对应全局顶点。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

# NOTE: this module deliberately does NOT import utillib.mylib / fittinglib, so
# the physical energy core stays independent of the (heavy, R/matplotlib-backed)
# ellipse-fitting dependency.  ``cells`` are only used structurally (each cell
# exposes ``cell.points``, a list of [x, y]); the type hints below are lazy
# thanks to ``from __future__ import annotations`` and never need ``Cell``.

TWO_PI = 2.0 * np.pi
_EPS = 1e-12


# --------------------------------------------------------------------------
# parameters
# --------------------------------------------------------------------------

@dataclass
class SymParams:
    """Coefficients for the two symmetric (angle-based) energy terms.

    The total energy is E = k_center * E_cen + k_angle * E_ang with the
    dimensionless ratio lambda = k_angle / k_center controlling the balance
    between cell-level shape regularity and vertex-level mechanical balance.

    Attributes:
        k_center: weight of the centre-angle term (cell-level shape regularity).
        k_angle: weight of the vertex-angle term (vertex-level balance).
        target_deg3: target sector angle at a degree-3 (tri-valent) vertex.
            Default 2*pi/3 is the Plateau / equal-tension force-balance value.
        target_deg2: target sector angle at a degree-2 (marginal) vertex.
            Default pi/2 is the foam free-boundary value (orthogonality of the
            internal edge to the boundary).
        centroid_mode: 'vertex' (default, barycentre -- decided A) gives the
            exact gradient formula and makes the central-angle force a pure
            azimuthal force (sum over a cell is zero, so uniform scaling does
            no work).  The geometric algorithm measures central angles from the
            *area* centroid (Cell.center_point, shoelace); 'area' is kept only
            for numerical equivalence checks against the original method.  The
            two differ by O(sigma/R) geometric fluctuations.
        NOTE: degree-1 vertices (patch tips) are excluded from E_ang and frozen
            by the annealer; the network is trivalent, so no degree >= 4
            targets exist.
    """
    k_center: float = 1.0
    k_angle: float = 1.0
    target_deg3: float = TWO_PI / 3.0
    target_deg2: float = np.pi / 2.0
    centroid_mode: str = "vertex"

    @property
    def active(self) -> bool:
        return (self.k_center != 0.0) or (self.k_angle != 0.0)

    @property
    def lam(self) -> Optional[float]:
        """Dimensionless ratio k_angle / k_center, or None if k_center == 0."""
        if self.k_center == 0.0:
            return None
        return self.k_angle / self.k_center


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def J(v: np.ndarray) -> np.ndarray:
    """Rotate by +90 degrees (counter-clockwise): J(x, y) = (-y, x)."""
    out = np.empty_like(v)
    out[..., 0] = -v[..., 1]
    out[..., 1] = v[..., 0]
    return out


def _as_ccw(verts: np.ndarray) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Return (verts_ccw, perm) with verts_ccw = verts[perm].

    ``perm`` is None when ``verts`` was already CCW.  The continuous lift used by
    :func:`centre_angle_terms` assumes a CCW traversal; feeding it a clockwise
    polygon silently winds n-1 times instead of once, which destroys the energy.
    """
    if verts.ndim != 2 or verts.shape[0] < 3:
        return verts, None
    x = verts[:, 0]
    y = verts[:, 1]
    signed = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
    if signed < 0:
        perm = np.arange(verts.shape[0] - 1, -1, -1)
        return verts[perm].copy(), perm
    return verts, None


def _unpermute(g: np.ndarray, perm: Optional[np.ndarray]) -> np.ndarray:
    """Map a gradient taken w.r.t. verts[perm] back to one w.r.t. verts."""
    if perm is None:
        return g
    out = np.empty_like(g)
    out[perm] = g
    return out


def _key_of(p) -> Tuple[float, float]:
    """Coordinate cluster key (9-decimal rounding, matches geometric code)."""
    return (round(float(p[0]), 9), round(float(p[1]), 9))


def _centroid(verts: np.ndarray, mode: str) -> np.ndarray:
    """Area centroid (shoelace) by default; vertex barycentre if mode=='vertex'."""
    if mode == "vertex":
        return verts.mean(axis=0)
    x = verts[:, 0]
    y = verts[:, 1]
    A = 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
    if abs(A) < _EPS:
        return verts.mean(axis=0)
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    cx = np.sum((x + np.roll(x, -1)) * cross) / (6.0 * A)
    cy = np.sum((y + np.roll(y, -1)) * cross) / (6.0 * A)
    return np.array([cx, cy], dtype=float)


# --------------------------------------------------------------------------
# centre-angle term (cell level)
# --------------------------------------------------------------------------

def centre_angle_terms(verts: np.ndarray,
                       centroid_mode: str = "vertex"
                       ) -> Tuple[float, np.ndarray, float, np.ndarray, np.ndarray]:
    """Return (E_raw, delta, phi_scalar, u, rho) for one cell.

    E_raw = sum_i delta_i^2 (multiply by k_center for the energy).  delta sums to
    zero by construction.  The cell is normalised to CCW first.
    """
    verts, _perm = _as_ccw(np.asarray(verts, dtype=float))
    n = verts.shape[0]
    if n < 3:
        return 0.0, np.zeros(n), 0.0, np.zeros((n, 2)), np.zeros(n)

    O = _centroid(verts, centroid_mode)
    d = verts - O
    rho = np.linalg.norm(d, axis=1)
    rho_safe = np.maximum(rho, _EPS)
    u = d / rho_safe[:, None]

    psi0 = np.arctan2(d[:, 1], d[:, 0])
    gaps = (np.roll(psi0, -1) - psi0) % TWO_PI
    psi_up = psi0[0] + np.concatenate(([0.0], np.cumsum(gaps[:-1])))

    i = np.arange(n)
    z = psi_up - TWO_PI * i / n

    # NOTE: arithmetic mean, not circular mean -- see module docstring note 1.
    phi = float(np.mean(z))
    delta = z - phi

    E_raw = float(np.sum(delta ** 2))
    return E_raw, delta, phi, u, rho_safe


def centre_angle_energy_and_grad(verts: np.ndarray,
                                 k_center: float,
                                 centroid_mode: str = "vertex"
                                 ) -> Tuple[float, np.ndarray]:
    """Centre-angle energy and its gradient w.r.t. each vertex of the cell."""
    verts_ccw, perm = _as_ccw(np.asarray(verts, dtype=float))
    if verts_ccw.shape[0] < 3 or k_center == 0.0:
        return 0.0, np.zeros_like(verts)

    E_raw, delta, _phi, u, rho = centre_angle_terms(verts_ccw, centroid_mode)
    n = verts_ccw.shape[0]

    Ju_over_rho = J(u) / rho[:, None]                       # (n, 2)
    shared = np.mean(delta[:, None] * Ju_over_rho, axis=0)  # (2,) centroid reaction

    grad = 2.0 * k_center * (delta[:, None] * Ju_over_rho - shared[None, :])
    return k_center * E_raw, _unpermute(grad, perm)


# --------------------------------------------------------------------------
# vertex-angle term (vertex level)
# --------------------------------------------------------------------------

def interior_angle_terms(verts: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                                     np.ndarray, np.ndarray]:
    """Interior angles of a CCW polygon and the pieces needed for gradients.

    规范 §2.3 记号：alpha = psi^up_a - psi^up_b，a = 前驱臂 (r_{i-1} - r_i)、
    b = 后继臂 (r_{i+1} - r_i)，顺序不可颠倒。本模块代码变量与规范相互置换：
    这里 a_i = r_{i+1} - r_i（后继/next 臂）、b_i = r_{i-1} - r_i（前驱/prev 臂），
    alpha 仍为同一内角（数值一致），梯度按置换后记号书写（见模块 docstring）。

    Returns (alpha, a_hat, b_hat, l_a, l_b) where, at local vertex i,
    alpha_i is the CCW angle from a_i to b_i with
        a_i = r_{i+1} - r_i   (arm towards the NEXT vertex, 规范的"后继臂")
        b_i = r_{i-1} - r_i   (arm towards the PREVIOUS vertex, 规范的"前驱臂")
    """
    n = verts.shape[0]
    nxt = np.roll(verts, -1, axis=0)
    prv = np.roll(verts, 1, axis=0)

    a = nxt - verts
    b = prv - verts
    l_a = np.maximum(np.linalg.norm(a, axis=1), _EPS)
    l_b = np.maximum(np.linalg.norm(b, axis=1), _EPS)
    a_hat = a / l_a[:, None]
    b_hat = b / l_b[:, None]

    alpha = np.arctan2(a_hat[:, 0] * b_hat[:, 1] - a_hat[:, 1] * b_hat[:, 0],
                       a_hat[:, 0] * b_hat[:, 0] + a_hat[:, 1] * b_hat[:, 1]) % TWO_PI
    return alpha, a_hat, b_hat, l_a, l_b


def vertex_targets(degree: np.ndarray,
                   sp: SymParams) -> Tuple[np.ndarray, np.ndarray]:
    """Target sector angle and inclusion mask at each global vertex.

    Returns (tau, mask) with mask == 1 where the vertex takes part in E_ang and
    0 where it is excluded.  Degree-1 vertices (patch tips, a single incident
    cell) are excluded from E_ang and frozen by the annealer; the network is
    trivalent, so no degree >= 4 targets exist.
    """
    degree = np.asarray(degree, dtype=int)
    tau = np.full(degree.shape[0], sp.target_deg3, dtype=float)
    mask = np.ones(degree.shape[0], dtype=float)

    tau[degree == 2] = sp.target_deg2
    mask[degree <= 1] = 0.0
    return tau, mask


# --------------------------------------------------------------------------
# global vertex index (clusters duplicated per-cell coordinates)
# --------------------------------------------------------------------------

def build_vertex_index(cells: List[Cell],
                        tol: float = 1e-9
                        ) -> Dict:
    """Cluster the duplicated per-cell vertex coordinates into global vertices.

    Returns a dict with keys:
        'coords'        : (G, 2) float array of unique vertex positions
        'refs'          : list (length G) of lists of (cell, local_index)
        'degree'        : (G,) int array, number of cells sharing each vertex
        'boundary_edges' : list (length G) of lists of unit tangent vectors of the
                           boundary edges incident to a marginal vertex (empty for
                           interior / frozen vertices)
        'gmap'          : dict mapping coordinate key -> global id (for fast lookup)

    The cluster key uses a 1e-9 rounding so vertices coinciding up to float noise
    are merged.  ``tol`` is accepted for API symmetry but the rounding precision
    is fixed at 9 decimals (consistent with the geometric code's ``points_equal``).
    """
    gmap: Dict[Tuple[float, float], int] = {}
    coords: List[List[float]] = []
    refs: List[List[Tuple[Cell, int]]] = []

    for cell in cells:
        for li, p in enumerate(cell.points):
            key = _key_of(p)
            gid = gmap.get(key)
            if gid is None:
                gid = len(coords)
                gmap[key] = gid
                coords.append([float(p[0]), float(p[1])])
                refs.append([])
            refs[gid].append((cell, li))

    G = len(coords)
    gcoords = np.asarray(coords, dtype=float)
    degree = np.array([len(r) for r in refs], dtype=int)

    # For every global vertex, the cells/edges it participates in, with the
    # global ids of the two polygon neighbours within each cell.
    neighbour_info: List[List[Tuple[Cell, int, int, int]]] = [[] for _ in range(G)]
    for cell in cells:
        local_keys = [_key_of(p) for p in cell.points]
        gid_of = [gmap[k] for k in local_keys]
        n_loc = len(cell.points)
        for li in range(n_loc):
            gid = gid_of[li]
            prev_gid = gid_of[(li - 1) % n_loc]
            next_gid = gid_of[(li + 1) % n_loc]
            neighbour_info[gid].append((cell, li, prev_gid, next_gid))

    boundary_edges: List[List[np.ndarray]] = [[] for _ in range(G)]
    for gid in range(G):
        if degree[gid] != 2:
            continue  # only marginal vertices slide along a boundary edge
        # A marginal (degree-2) vertex is shared by exactly two cells.  Within
        # each cell it has two polygon neighbours; the neighbour shared by BOTH
        # cells is the internal vertex O, and the other two neighbours (one per
        # cell) lie on the boundary rim.  The boundary edges are gid -> those two
        # rim neighbours, regardless of whether the rim vertices themselves are
        # degree 2 (closed rim) or degree 1 (open tip).
        from collections import defaultdict
        by_cell = defaultdict(list)
        for (_cell, _li, prev_gid, next_gid) in neighbour_info[gid]:
            by_cell[id(_cell)].append((prev_gid, next_gid))
        if len(by_cell) < 2:
            continue
        cell_nb_sets = []
        for cell_pairs in by_cell.values():
            s = set()
            for (prev_gid, next_gid) in cell_pairs:
                s.add(prev_gid)
                s.add(next_gid)
            cell_nb_sets.append(s)
        shared = cell_nb_sets[0] & cell_nb_sets[1]
        o_gid = next(iter(shared)) if shared else None
        seen_edges = set()
        for (_cell, _li, prev_gid, next_gid) in neighbour_info[gid]:
            for nb in (prev_gid, next_gid):
                if nb == o_gid:
                    continue
                if nb in seen_edges:
                    continue
                seen_edges.add(nb)
                t = gcoords[nb] - gcoords[gid]
                nrm = np.linalg.norm(t)
                if nrm > _EPS:
                    boundary_edges[gid].append(t / nrm)

    return {
        'coords': gcoords,
        'refs': refs,
        'degree': degree,
        'boundary_edges': boundary_edges,
        'neighbour_info': neighbour_info,
        'gmap': gmap,
    }


# --------------------------------------------------------------------------
# tissue-level energy and gradient
# --------------------------------------------------------------------------

@dataclass
class SymEnergyReport:
    """Breakdown returned by :func:`total_energy` / :func:`total_energy_gradient`."""
    total: float = 0.0
    centre: float = 0.0
    angle: float = 0.0
    centre_per_cell: Dict = None
    angle_per_cell: Dict = None

    def __post_init__(self):
        if self.centre_per_cell is None:
            self.centre_per_cell = {}
        if self.angle_per_cell is None:
            self.angle_per_cell = {}


def _cell_gids(cell: Cell, verts_ccw: np.ndarray, perm, vidx: Dict) -> List[int]:
    """Global ids aligned with ``verts_ccw`` order (handles CCW reversal)."""
    raw = np.asarray(cell.points, dtype=float)
    orig = [_key_of(raw[li]) for li in range(raw.shape[0])]
    gmap = vidx['gmap']
    if perm is None:
        return [gmap[k] for k in orig]
    return [gmap[orig[p]] for p in perm]


def total_energy_gradient(cells: List[Cell],
                           sp: SymParams,
                           vidx: Optional[Dict] = None
                           ) -> Tuple[SymEnergyReport, np.ndarray]:
    """Symmetric energy and its gradient over the global (clustered) vertices.

    Returns (report, grad) with grad of shape (G, 2) matching vidx['coords'].
    The gradient is dE/d r_v; the mechanical force is its negative.
    """
    if vidx is None:
        vidx = build_vertex_index(cells)
    G = vidx['coords'].shape[0]
    grad = np.zeros((G, 2), dtype=float)
    rep = SymEnergyReport()
    if G == 0 or not sp.active:
        return rep, grad

    degree = vidx['degree']
    tau, mask = vertex_targets(degree, sp)

    for cell in cells:
        raw = np.asarray(cell.points, dtype=float)
        verts_ccw, perm = _as_ccw(raw)
        if verts_ccw.shape[0] < 3:
            continue
        n_loc = verts_ccw.shape[0]
        gids = _cell_gids(cell, verts_ccw, perm, vidx)
        gids_arr = np.asarray(gids, dtype=int)

        # ---- centre-angle term (per cell) ----------------------------------
        if sp.k_center != 0.0:
            E_c, g_c = centre_angle_energy_and_grad(verts_ccw, sp.k_center, sp.centroid_mode)
            rep.centre += E_c
            rep.centre_per_cell[id(cell)] = rep.centre_per_cell.get(id(cell), 0.0) + E_c
            np.add.at(grad, gids_arr, g_c)

        # ---- vertex-angle term (per cell) ----------------------------------
        if sp.k_angle != 0.0:
            alpha, a_hat, b_hat, l_a, l_b = interior_angle_terms(verts_ccw)
            Ja_over_la = J(a_hat) / l_a[:, None]      # arm to NEXT
            Jb_over_lb = J(b_hat) / l_b[:, None]      # arm to PREV
            g_v = Ja_over_la - Jb_over_lb             # d alpha_i / d r_i
            nxt_arr = np.roll(gids_arr, -1)
            prv_arr = np.roll(gids_arr, 1)
            for ll in range(n_loc):
                gg = gids_arr[ll]
                m = mask[gg]
                if m <= 0.0:
                    continue
                dev = alpha[ll] - tau[gg]
                w = 2.0 * sp.k_angle * dev
                np.add.at(grad, gg, w * g_v[ll])
                np.add.at(grad, nxt_arr[ll], w * (-Ja_over_la[ll]))
                np.add.at(grad, prv_arr[ll], w * (Jb_over_lb[ll]))
                rep.angle += sp.k_angle * dev * dev

    rep.total = rep.centre + rep.angle
    return rep, grad


def total_energy(cells: List[Cell],
                 sp: SymParams,
                 vidx: Optional[Dict] = None
                 ) -> SymEnergyReport:
    """Total symmetric energy and a per-cell breakdown (delegates to gradient)."""
    rep, _g = total_energy_gradient(cells, sp, vidx)
    return rep
