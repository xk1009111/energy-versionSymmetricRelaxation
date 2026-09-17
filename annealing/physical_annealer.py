# -*- coding: utf-8 -*-
"""
Physical (energy-based) annealer —— 由数学版（几何规则）退火修改而来。

用单个标量对称能量上的过阻尼梯度下降（前向 Euler），替代数学版细胞退火
（原 ``annealerUtil.py`` 的逐顶点几何移动规则：最优射线 → 三角形质心目标点 →
移动，已删除）。全部公式与符号以 ``SymRelax_Energy_Model_Spec.md``（完整规范）
与 ``SymRelax_Energy_Force_Scheme.md``（力计算方案）为准。能量（规范 §2，
裸和，决策 #15）：

    E_sym = kappa_c * E_cen + kappa_a * E_ang

implemented in :mod:`annealing.energy`.  Every relaxation round:

    1.  build the global vertex index (cluster duplicated per-cell coordinates);
    2.  compute the analytic gradient dE/d r_v of the whole network
        (force F_v = -dE/d r_v: centre-angle force = direct term + shared
        term M_c per incident cell, spec §3.1; angle force = deviation x edge
        gradient, spec §3.2);
    3.  explicit Euler step (前向 Euler，规范 §4.1).  Default: applied to ALL
        vertices simultaneously (同步更新 / Jacobi，决策 #6).  If
        ``serial_descending=True`` the round instead does a serial Gauss-Seidel
        sweep: free vertices are ordered by descending force magnitude |F_v|
        (energy-high-first, the physical counterpart of the geometric version's
        descending anneal distance) and moved one at a time along their own
        force, each with its own line search + local convexity guard:
            r_v <- r_v + eta * F_v ,   eta = mobility (mu = dt / zeta)
        (mu = 迁移率，唯一时间尺度参数；默认按中位边长 L0 归一化
        （reversing 决策 #5）：有效步 = mu·L0²，位移 ≈ mu·L0，
        与网络绝对尺度无关):
          * interior vertices (degree == 3): free 2-D move;
          * marginal vertices (degree == 2): the force is projected onto the
            boundary tangent and the vertex slides *along* the edge (规范
            §4.2：两条边界边中取 |F·t| 较大者的切向；动力学约束而非能量项，
            对应数学版「边缘点沿小角对应边移动」，决策 #7);
          * degree-1 (patch-tip) vertices: frozen;
    4.  step-size line search (线搜索，规范 §4.3 / 决策 #8): if the step
        raises the total energy or breaks a cell that was convex at round
        start ("起始即凸"), eta is halved and recomputed -- monotone
        convergence; cells starting non-convex are exempt (起始即非凸的细胞豁免).
        Synchronous mode (default): the whole simultaneous step is recomputed
        from the round-start positions.  Serial mode (serial_descending=True):
        each vertex has its OWN step halved independently and the convexity
        guard is localised to the cells it touches, so one fragile cell cannot
        veto the whole network (see :meth:`_one_round_serial`);
    5.  flush (recompute each cell's centre / area; 面积仅供显示/导出，
        能量模型不使用面积，规范 §5 步骤 7).

This is the physical counterpart of ``Annealer`` and shares the same
``annealing(cells)`` entry point so the GUI can switch between the two modes.

NOTE on scale (规范 §7.2-7.3): the model energy E is strictly scale-invariant
(only angles matter), so the force scales as |F| ~ kappa * delta / rho.  With
normalize_L0=True (default, reversing spec decision #5) the effective step is
mu * L0**2, so each round moves a vertex by ~ mu * L0 -- a dimensionless
fraction of the median edge length, identical on networks of any absolute
scale.  The GUI "步长 μ" knob is therefore scale-robust: pick mu once and it
behaves the same on small and large tessellations.  Set normalize_L0=False to
recover the plain absolute-step behaviour (mu = dt/zeta, scale-dependent).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from annealing import energy as E


# Local reimplementation of the convex-polygon test from annealerUtil, so this
# module does not pull in the (heavy, R/matplotlib-backed) fittinglib chain via
# utillib.mylib.  points: list of [x, y] in polygon order.
def _is_polygon_convex(points, eps=1e-9):
    pts = [tuple(p) for p in points]
    n = len(pts)
    if n < 3:
        return False
    area = 0.0
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    orient = 1 if area >= 0 else -1
    prev_sign = 0
    for i in range(n):
        p0 = pts[i]
        p1 = pts[(i + 1) % n]
        p2 = pts[(i + 2) % n]
        v1 = (p1[0] - p0[0], p1[1] - p0[1])
        v2 = (p2[0] - p1[0], p2[1] - p1[1])
        cross = v1[0] * v2[1] - v1[1] * v2[0]
        if abs(cross) < eps:
            continue
        sign = 1 if cross > 0 else -1
        if prev_sign == 0:
            prev_sign = sign
        elif sign != prev_sign:
            return False
    if prev_sign == 0 or prev_sign != orient:
        return False
    return True


def _median_edge_length(cells, eps=1e-12):
    """Median length of all cell edges (L0, used for scale normalisation).

    Edges are taken as consecutive vertex pairs within each cell (including the
    closing edge).  Shared edges are double-counted, which does not bias the
    median.  Falls back to 1.0 on a degenerate (edgeless) tessellation.
    """
    lengths = []
    for cell in cells:
        pts = cell.points
        n = len(pts)
        for i in range(n):
            p0 = np.asarray(pts[i], dtype=float)
            p1 = np.asarray(pts[(i + 1) % n], dtype=float)
            d = float(np.linalg.norm(p1 - p0))
            if d > eps:
                lengths.append(d)
    if not lengths:
        return 1.0
    return float(np.median(np.asarray(lengths)))


@dataclass
class PhysicalAnnealerParams:
    """Parameters for :class:`PhysicalAnnealer`.

    Attributes:
        mobility: the over-damped step constant mu = dt / zeta (spec §4.1).
            One round moves each free vertex by mobility * |F_v|; the line
            search may shrink it.  This is the *single* time-scale parameter.
            With normalize_L0=True (default, reversing spec decision #5) the
            effective step is mobility * L0**2 (L0 = median edge length), so
            mobility becomes a dimensionless fraction and each vertex moves
            ~ mobility * L0 per round -- identical convergence on any network
            scale.  Set normalize_L0=False for the old absolute-step behaviour.
        k_center: weight of the centre-angle (cell-level) term.  Fixed to 1.0 by
            the energy-model convention; lambda = k_angle / k_center is the
            single dimensionless input.
        k_angle: weight of the vertex-angle (vertex-level) term (the lambda
            input; k_center = 1).
        centroid_mode: passed through to :class:`E.SymParams` ('vertex' =
            vertex barycentre, decided A).
        convex_guard: if True (default), the line search also requires that
            cells which were convex at round start stay convex (cells starting
            non-convex are exempt so the angle energy can convexify them).
        freeze_marginal: if True, marginal vertices are frozen (degree-2 excluded
            from the move, like the geometric default marginal_point_judge=0).
            The physical model slides them along the edge, so default False.
        serial_descending: if True (default), every free vertex is moved in a
            serial Gauss-Seidel sweep ordered by descending force magnitude
            |F_v| (energy-high-first -- the physical counterpart of the
            geometric version's descending anneal distance), each with its own
            line search and a convexity guard localised to the cells it touches.
            This stops one fragile cell from vetoing the whole network (决策 #6
            的死穴) and usually converges faster; it trades the round-level
            monotone-energy guarantee for overall faster convergence.  If False,
            all free vertices are moved SIMULTANEOUSLY (Jacobi, 决策 #6).
    """
    mobility: float = 0.1
    k_center: float = 1.0
    k_angle: float = 1.0
    centroid_mode: str = "vertex"
    convex_guard: bool = True
    freeze_marginal: bool = False
    normalize_L0: bool = True
    serial_descending: bool = True

    @property
    def lam(self) -> float:
        """Dimensionless ratio lambda = k_angle / k_center."""
        return self.k_angle / self.k_center

    def to_sym_params(self) -> E.SymParams:
        return E.SymParams(
            k_center=self.k_center,
            k_angle=self.k_angle,
            centroid_mode=self.centroid_mode,
        )


class PhysicalAnnealer:
    """Over-damped gradient-descent annealer on the symmetric energy."""

    def __init__(self, params: Optional[PhysicalAnnealerParams] = None):
        self.params = params or PhysicalAnnealerParams()
        self.inner_points: int = 0
        self.marginal_points: int = 0
        self.frozen_points: int = 0

    # ------------------------------------------------------------------
    def one_round(self,
                  cells: List[Cell],
                  sp: Optional[E.SymParams] = None,
                  eta: Optional[float] = None
                  ) -> Tuple[E.SymEnergyReport, int]:
        """Perform one relaxation round.

        Synchronous mode (default, ``serial_descending=False``): explicit Euler,
        all free vertices moved simultaneously (Jacobi, 决策 #6) with a single
        global line search.  Serial mode (``serial_descending=True``): free
        vertices are ordered by descending force magnitude |F_v| and moved one
        at a time (Gauss-Seidel) along their own force, each with its own line
        search + local convexity guard -- see :meth:`_one_round_serial`.

        Returns (energy_report, n_moved).  Mutates ``cells`` in place.

        Per SymRelax_Energy_Model_Spec.md sec. 3-4: all free vertices move
        SIMULTANEOUSLY along the force F_v = -dE/d r_v with the over-damped
        step r_v <- r_v + step * F_v, step = eta * L0**2 (L0 = median edge
        length) when normalize_L0=True, else step = eta (mu = dt/zeta).
        The speed decays automatically as the force vanishes at equilibrium.

            * interior vertices (degree == 3): free 2-D move eta * F_v;
            * marginal vertices (degree == 2): the force is projected onto the
              boundary tangent, eta * (F_v . t) t -- the vertex only slides
              *along* the boundary edge;
            * degree-1 (patch-tip) vertices: frozen.

        Step-size line search (spec §4.3): if the step raises the total
        energy or breaks a cell that was convex at round start ("起始即凸"),
        eta is halved.  Synchronous mode recomputes the WHOLE step; serial mode
        halves each vertex's own step.  Cells that start non-convex are exempt,
        so the angle energy can still convexify them.
        """
        if sp is None:
            sp = self.params.to_sym_params()
        vidx = E.build_vertex_index(cells)
        coords = vidx['coords']
        degree = vidx['degree']
        refs = vidx['refs']
        boundary_edges = vidx['boundary_edges']
        G = coords.shape[0]
        if G == 0:
            return E.SymEnergyReport(), 0

        # eta = mobility (mu = dt / zeta).  Scale-invariant normalisation
        # (reversing spec decision #5): the effective step is eta * L0**2 where
        # L0 is the median edge length, so each round moves a vertex by
        # ~ eta * L0 (eta is now a dimensionless fraction, independent of the
        # absolute network scale).  Set params.normalize_L0 = False to recover
        # the old absolute-step behaviour.
        if eta is None:
            eta = self.params.mobility
        L0 = _median_edge_length(cells)
        step = eta * (L0 * L0) if self.params.normalize_L0 else eta

        rep, grad = E.total_energy_gradient(cells, sp, vidx)
        E0 = rep.total
        base = coords.copy()

        # Convexity of each cell at the START of the round.  The line search
        # only forbids turning a *convex* cell into a non-convex one; a cell
        # that is already non-convex (e.g. a dented boundary cell from the
        # initial tessellation) is exempt, and the angle energy naturally
        # convexifies it.
        convex0 = {id(cell): _is_polygon_convex(cell.points) for cell in cells}

        free = (degree > 1)
        if self.params.freeze_marginal:
            free = free & (degree != 2)

        # Serial descending sweep (energy-high-first, Gauss-Seidel) -- see the
        # class docstring / previous turn.  Keeps the synchronous Jacobi path
        # exactly as-is when the flag is off.
        if self.params.serial_descending:
            return self._one_round_serial(
                cells, sp, step, vidx, grad, base, convex0, free, degree)

        def _apply(cand):
            for gid in range(G):
                for (cell, li) in refs[gid]:
                    cell.points[li] = [float(cand[gid, 0]), float(cand[gid, 1])]

        def _flush():
            for cell in cells:
                cell.setVertex()
                cell.setArea()

        def _convex_ok():
            if not self.params.convex_guard:
                return True
            for cell in cells:
                if convex0.get(id(cell), True) and not _is_polygon_convex(cell.points):
                    return False
            return True

        # line search: halve eta and recompute the whole simultaneous step
        # until the energy does not rise AND starting-convex cells stay convex
        new_rep = None
        accepted = False
        for _try in range(60):
            # explicit Euler candidate: r_v <- r_v + step * F_v (F_v = -grad),
            # applied to ALL free vertices simultaneously
            cand = base.copy()
            if free.any():
                for gid in range(G):
                    if not free[gid]:
                        continue  # frozen patch tip / marginal
                    g = grad[gid]
                    if degree[gid] == 2 and boundary_edges[gid]:
                        # 边缘顶点切向投影（规范 §4.2 / 决策 #16）：t 取两条
                        # 边界边中 |F·t| 较大者的单位方向，沿边界边滑动
                        t = None
                        best = -1.0
                        for e in boundary_edges[gid]:
                            align = float(np.dot(g, e))
                            if abs(align) > best:
                                best = abs(align)
                                t = e
                        cand[gid] = base[gid] - step * float(np.dot(g, t)) * t
                    else:
                        cand[gid] = base[gid] - step * g
            _apply(cand)

            if not _convex_ok():
                _apply(base)
                step *= 0.5
                continue

            _flush()
            r_try = E.total_energy(cells, sp)
            if r_try.total <= E0 + 1e-9:
                new_rep = r_try
                accepted = True
                break

            # energy rose: restore the start positions and halve the step
            _apply(base)
            step *= 0.5

        if not accepted:
            _apply(base)
            _flush()
            new_rep = rep

        # count vertices that actually moved from their starting position
        moved = 0
        n_interior = 0
        n_marginal = 0
        if accepted:
            for gid in range(G):
                if not free[gid]:
                    continue
                (cell0, li0) = refs[gid][0]
                cur = np.asarray(cell0.points[li0], dtype=float)
                if np.linalg.norm(cur - base[gid]) > 1e-12:
                    moved += 1
                    if degree[gid] == 2:
                        n_marginal += 1
                    else:
                        n_interior += 1

        self.inner_points = n_interior
        self.marginal_points = n_marginal
        self.frozen_points = int((degree <= 1).sum())
        # report the FINAL energy of this round so a caller comparing successive
        # rounds sees the energy trajectory.
        return new_rep, moved

    # ------------------------------------------------------------------
    def _one_round_serial(self,
                          cells: List[Cell],
                          sp: E.SymParams,
                          step: float,
                          vidx: dict,
                          grad: np.ndarray,
                          base: np.ndarray,
                          convex0: dict,
                          free: np.ndarray,
                          degree: np.ndarray
                          ) -> Tuple[E.SymEnergyReport, int]:
        """Serial Gauss-Seidel sweep, free vertices ordered by descending force
        magnitude |F_v| (energy-high-first).

        Physical counterpart of the geometric annealer's descending anneal-
        distance sweep.  Each free vertex is moved one at a time along its own
        force F_v = -grad[v] (tangent-projected for marginals), with its own
        line search and a convexity guard restricted to the cells it touches --
        so a single fragile cell can no longer veto the whole network the way
        the synchronous Jacobi step can (决策 #6 的死穴).  Returns
        (energy_report, n_moved) and sets the self.*_points counters.
        """
        refs = vidx['refs']
        G = base.shape[0]
        boundary_edges = vidx['boundary_edges']

        # Per-cell ordered global-id list (local-index aligned), frozen at round
        # start so it stays valid after vertices move (the gmap keys would
        # otherwise drift once a vertex is written to a new coordinate).
        cell_gids_ordered: Dict = {}
        for gid in range(G):
            for (cell, li) in refs[gid]:
                cid = id(cell)
                lst = cell_gids_ordered.get(cid)
                if lst is None:
                    lst = [None] * len(cell.points)
                    cell_gids_ordered[cid] = lst
                lst[li] = gid

        tau, mask = E.vertex_targets(degree, sp)

        def _apply_vertex(gid, pos):
            for (cell, li) in refs[gid]:
                cell.points[li] = [float(pos[0]), float(pos[1])]

        def _flush_vertex(gid):
            for (cell, li) in refs[gid]:
                cell.setVertex()
                cell.setArea()

        def _cell_energy(cell):
            raw = np.asarray(cell.points, dtype=float)
            verts, perm = E._as_ccw(raw)
            n = verts.shape[0]
            if n < 3:
                return 0.0
            gid_list = cell_gids_ordered[id(cell)]
            gids_ccw = [gid_list[perm[i]] for i in range(n)] if perm is not None else gid_list
            e = 0.0
            if sp.k_center != 0.0:
                Ec, _ = E.centre_angle_energy_and_grad(verts, sp.k_center, sp.centroid_mode)
                e += Ec
            if sp.k_angle != 0.0:
                alpha, a_hat, b_hat, l_a, l_b = E.interior_angle_terms(verts)
                for ll in range(n):
                    gg = gids_ccw[ll]
                    if mask[gg] <= 0.0:
                        continue
                    dev = alpha[ll] - tau[gg]
                    e += sp.k_angle * dev * dev
            return e

        def _local_convex_ok(gid):
            if not self.params.convex_guard:
                return True
            for (cell, li) in refs[gid]:
                if convex0.get(id(cell), True) and not _is_polygon_convex(cell.points):
                    return False
            return True

        def _local_energy(gid):
            e = 0.0
            for (cell, li) in refs[gid]:
                e += _cell_energy(cell)
            return e

        moved = 0
        n_interior = 0
        n_marginal = 0
        if free.any():
            force_mag = np.linalg.norm(grad, axis=1)
            order = np.argsort(-force_mag)  # descending |F_v|
            for gid in order:
                if not free[gid]:
                    continue
                g = grad[gid]
                # displacement direction per unit step (F_v = -grad)
                if degree[gid] == 2 and boundary_edges[gid]:
                    t = None
                    best = -1.0
                    for e in boundary_edges[gid]:
                        align = float(np.dot(g, e))
                        if abs(align) > best:
                            best = abs(align)
                            t = e
                    d = -float(np.dot(g, t)) * t
                else:
                    d = -g
                e_before = _local_energy(gid)
                accepted_v = False
                st = step
                for _try in range(60):
                    cand = base[gid] + st * d
                    _apply_vertex(gid, cand)
                    if not _local_convex_ok(gid):
                        _apply_vertex(gid, base[gid])
                        st *= 0.5
                        continue
                    e_after = _local_energy(gid)
                    if e_after <= e_before + 1e-9:
                        _flush_vertex(gid)
                        base[gid] = cand.copy()
                        accepted_v = True
                        moved += 1
                        if degree[gid] == 2:
                            n_marginal += 1
                        else:
                            n_interior += 1
                        break
                    _apply_vertex(gid, base[gid])
                    st *= 0.5
                # if rejected on all 60 halvings, vertex stays at base[gid]

        self.inner_points = n_interior
        self.marginal_points = n_marginal
        self.frozen_points = int((degree <= 1).sum())
        new_rep = E.total_energy(cells, sp)
        return new_rep, moved

    # convenience: relax until convergence / cap --------------------------
    def relax(self,
              cells: List[Cell],
              max_rounds: int = 2000,
              tol: float = 1e-9,
              record_every: int = 0
              ) -> List[float]:
        """Relax until the energy change per round drops below ``tol`` or
        ``max_rounds`` is reached.  Returns the list of total energies per round
        (every round, or every ``record_every`` rounds if > 0)."""
        sp = self.params.to_sym_params()
        history: List[float] = []
        prev = None
        for r in range(max_rounds):
            rep, n = self.one_round(cells, sp)
            if record_every and (r % record_every == 0):
                history.append(rep.total)
            elif not record_every:
                history.append(rep.total)
            # convergence must be judged on rounds that actually moved: a
            # controller-rejected round reports an unchanged energy (n == 0)
            # and must not be mistaken for convergence.
            if n > 0 and prev is not None and abs(rep.total - prev) < tol:
                break
            prev = rep.total
        return history
