# -*- coding: utf-8 -*-
"""退火器统一入口（能量版）—— 由数学版（几何规则）退火修改而来。

历史说明：本文件曾是「数学版（几何规则）/ 能量版（变分）」的双模式分发器。
数学版（math_annealer.py 及其依赖 annealerUtil.py）已删除，现仅保留能量版，
``Annealer`` 直接再导出 ``energy_annealer.Annealer``。

能量版公式与符号以 ``SymRelax_Energy_Model_Spec.md``（完整规范）与
``SymRelax_Energy_Force_Scheme.md``（力计算方案）为准：对对称能量
E = κ_c·E_cen + κ_a·E_ang（裸和，κ_c=1、λ=κ_a/κ_c 为唯一无量纲输入）做
过阻尼梯度下降（前向 Euler，规范 §4.1 / 决策 #6）。默认按中位边长 L0 归一化
（reversing 决策 #5）：有效步 = μ·L0²，每轮位移 ≈ μ·L0，μ 为无量纲分数、
跨尺度一致。全体顶点同步更新（决策 #6）；边缘顶点沿边界切向滑动（规范 §4.2），
度数 1（贴片尖端）顶点冻结。

GUI 通过 ``from annealing.AnnealingGUI import Annealer`` 使用，接口不变：
``init(params)`` / ``annealing(cells)`` / 统计属性 ``inner_points`` /
``marginal_points`` / ``marginal_point_judge``。

params 支持键（详见 energy_annealer 模块 docstring）：
    - annealingRate         : 迁移率 μ（默认 L0 归一化下为无量纲分数）
    - marginal_point_judge  : 边缘顶点是否参与退火
    - inner_angle_sq_guard  : 守卫开关（未显式给 λ 时决定 λ 默认值；
                              勾选时 λ=10 模拟数学版硬约束 λ→∞ 极限）
    - lam                   : 显式 λ = κ_a/κ_c，给出时优先于守卫默认值
"""

from annealing.energy_annealer import Annealer

__all__ = ["Annealer"]
