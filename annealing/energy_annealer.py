# -*- coding: utf-8 -*-
"""
能量版（变分）细胞退火器 —— 由数学版（几何规则）退火修改而来。

按《对称退火 · 数学版 → 能量版 完整逻辑与公式规范》
（``SymRelax_Energy_Model_Spec.md``，下称「规范」）与
《对称松弛 · 能量版力计算方案》（``SymRelax_Energy_Force_Scheme.md``）实现：
对对称能量（规范 §2，裸和，决策 #15）

    E = kappa_c * E_cen + kappa_a * E_ang,   lambda = kappa_a / kappa_c

做过阻尼梯度下降（规范 §4.1），前向 Euler 同步更新所有顶点（决策 #6；
kappa_c 固定为 1，lambda 为唯一无量纲输入，规范 §6.1 连续旋钮）。
边缘顶点沿边界切向滑动（规范 §4.2，决策 #7），度数 1（贴片尖端）顶点冻结。

与数学版（原几何规则实现，已删除）保持相同的调用接口（类名、init(params)、
annealing(cells)、统计属性 inner_points / marginal_points），GUI 无需区分。

入口参数（params 字典）：
    - annealingRate   : 迁移率 μ = Δt/ζ（步长，规范 §4.1）。能量版默认按中位
                        边长 L0 归一化（reversing 决策 #5）：有效步 = μ·L0²，
                        位移 ≈ μ·L0，μ 为无量纲分数、与网络绝对尺度无关，换
                        不同尺度网络无需重调；设 PhysicalAnnealerParams
                        .normalize_L0=False 可还原纯绝对步长（旧行为）。
    - marginal_point_judge : 边缘顶点是否参与退火 (True/False 或 1/0)
    - inner_angle_sq_guard : 内角平方和守卫。勾选时默认 lambda=10（模拟原算法
                             角度硬约束在 lambda→∞ 的极限，规范 §6.1）；不勾选
                             时 lambda=1（软加权）。λ 扫描时可用下方 lam 显式覆盖。
    - lam            : 显式输入 lambda (kappa_a/kappa_c)；给出时优先于守卫默认。

注：能量模型本身（annealing.energy / annealing.physical_annealer）严格按规范
实现；步长默认按 L0 归一化（有效步 = μ·L0²，位移 ≈ μ·L0），使 μ 成为无量纲
分数、跨尺度一致。置 normalize_L0=False 可还原纯绝对步长。
"""

from __future__ import annotations

from annealing import energy as E
from annealing.physical_annealer import PhysicalAnnealer, PhysicalAnnealerParams


class Annealer:
    """细胞退火器（能量版）。"""

    annealingRate = 0.0
    marginal_point_judge = False
    k_center = 1.0
    lam = 1.0
    inner_points = 0
    marginal_points = 0
    _serial_descending = True

    def __init__(self, params=None):
        """初始化退火器（能量版）。:param params: 参数字典，见模块 docstring。"""
        if params is None:
            params = {}
        self.init(params)

    # 细胞退火器参数初始化（参数化版本）
    def init(self, params=None):
        """参数化初始化方法。:param params: 参数字典，见模块 docstring。"""
        if params is None:
            params = {}

        print("↓↓===细胞退火器参数设置（能量版：E = E_cen + λ·E_ang）===↓↓")
        self.annealingRate = float(params.get('annealingRate', 0.5))  # 默认0.5

        inner_angle_sq_guard = params.get('inner_angle_sq_guard', 1)
        if isinstance(inner_angle_sq_guard, bool):
            guard = inner_angle_sq_guard
        elif isinstance(inner_angle_sq_guard, (int, str)):
            guard = int(inner_angle_sq_guard) == 1
        else:
            guard = True

        marginal_point_judge = params.get('marginal_point_judge', 0)  # 默认不参与
        if isinstance(marginal_point_judge, bool):
            self.marginal_point_judge = marginal_point_judge
        elif isinstance(marginal_point_judge, (int, str)):
            self.marginal_point_judge = int(marginal_point_judge) == 1
        else:
            self.marginal_point_judge = False

        # λ = kappa_a / kappa_c；kappa_c 固定 1。显式输入优先，否则由守卫开关给默认。
        lam = params.get('lam', None)
        if lam is None:
            lam = 10.0 if guard else 1.0
        self.lam = float(lam)
        self.k_center = 1.0
        self._freeze_marginal = not self.marginal_point_judge

        # 串行降序退火（按 |F_v| 降序的 Gauss-Seidel 串行更新；默认开）
        sd = params.get('serial_descending', 1)
        if isinstance(sd, bool):
            self._serial_descending = sd
        elif isinstance(sd, (int, str)):
            self._serial_descending = int(sd) == 1
        else:
            self._serial_descending = False

        print(f"mobility μ(=annealingRate)={self.annealingRate}, lambda(=kappa_a/kappa_c)={self.lam}, "
              f"kappa_c={self.k_center}, marginal_slide={self.marginal_point_judge}")
        print("----------------------------")

    def annealing(self, cells):
        """一轮退火：对所有自由顶点做一次显式 Euler 同步梯度下降（能量版）。"""
        sp = E.SymParams(k_center=self.k_center, k_angle=self.lam,
                         centroid_mode="vertex")  # 已定 A：顶点质心
        # 迁移率 μ = Δt/ζ；默认 normalize_L0=True：有效步 = μ·L0²（L0=中位边长），
        # 每轮位移 ≈ μ·L0，μ 为无量纲分数、跨尺度一致（reversing 决策 #5）。
        # 设 normalize_L0=False 可还原纯绝对步长（μ = Δt/ζ，与尺度耦合）。
        mu = self.annealingRate

        pa = PhysicalAnnealer(PhysicalAnnealerParams(
            mobility=mu,
            k_center=self.k_center,
            k_angle=self.lam,
            centroid_mode="vertex",
            convex_guard=True,
            freeze_marginal=self._freeze_marginal,
            serial_descending=self._serial_descending,
        ))
        pa.one_round(cells, sp, eta=mu)
        self.inner_points = int(pa.inner_points)
        self.marginal_points = int(pa.marginal_points)
