# -*- coding: utf-8 -*-
import os
import sys
import tkinter.messagebox as msg

# =========================================================================
# 步骤 1：定义环境配置与弹窗提示函数
# =========================================================================
def _check_system_r():
    """
    检测系统是否已安装 R
    返回: True (系统 R 可用), False (不可用)
    """
    # 检查 R_HOME 环境变量
    r_home = os.environ.get('R_HOME', '')
    if r_home and os.path.exists(r_home):
        return True

    # 检查 PATH 中是否有 R
    r_paths = os.environ.get('PATH', '').split(os.pathsep)
    for p in r_paths:
        r_exe = os.path.join(p, 'R.exe')
        if os.path.exists(r_exe):
            return True

    # 尝试通过 where 命令查找 (Windows)
    try:
        import subprocess
        result = subprocess.run(['where', 'R'], capture_output=True, text=True, timeout=3)
        if result.returncode == 0 and result.stdout.strip():
            return True
    except Exception:
        pass

    return False


def _use_portable_r(base_path):
    """配置并使用便携版 R_Dist"""
    portable_r_home = os.path.join(base_path, 'R_Dist')
    portable_r_bin_root = os.path.join(portable_r_home, 'bin')
    portable_r_bin = os.path.join(portable_r_bin_root, 'x64')

    os.environ['R_HOME'] = portable_r_home
    os.environ['R_USER'] = portable_r_home
    os.environ['PATH'] = portable_r_bin + os.pathsep + portable_r_bin_root + os.pathsep + os.environ.get('PATH', '')

    if hasattr(os, 'add_dll_directory'):
        try:
            os.add_dll_directory(portable_r_bin)
            os.add_dll_directory(portable_r_bin_root)
        except Exception:
            pass

    return portable_r_home


def setup_portable_r():
    """
    配置 R 环境，优先级：
    1. 用户系统已安装的 R（通过 R_HOME / PATH 检测）
    2. 便携版 R_Dist（位于本程序目录下）
    3. 都不可用 → 降级为纯 Python 拟合（使用 scipy/numpy）
    """
    if getattr(sys, 'frozen', False):
        # exe 运行时：R_Dist 在 exe 同目录
        base_path = os.path.dirname(sys.executable)
    else:
        # 源码运行时：R_Dist 在 fittinglib.py 的上上级目录
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    status_title = "R 语言环境检测"
    _is_portable = False

    # —— 优先级 1：检测系统 R ——
    if _check_system_r():
        r_home = os.environ.get('R_HOME', '系统 PATH 中')
        msg.showinfo(status_title,
            f"✅ 检测到系统中已安装 R (System R detected)\n\n"
            f"R_HOME: {r_home}\n\n"
            "将使用系统 R 进行椭圆拟合。(Will use system R for ellipse fitting.)")
        return

    # —— 优先级 2：检测便携版 R_Dist ——
    portable_r_home_candidate = os.path.join(base_path, 'R_Dist')
    if os.path.exists(portable_r_home_candidate):
        r_path = _use_portable_r(base_path)
        _is_portable = True
        msg.showinfo(status_title,
            f"✅ 使用便携版 R 环境 (Using portable R environment)\n\n"
            f"路径: {r_path}\n\n"
            "说明：未检测到系统 R，已自动挂载内置 R。(Note: System R not detected, built-in R auto-mounted.)")
        return

    # —— 均不可用 ——
    system_r_home = os.environ.get('R_HOME', '未检测到')
    msg.showwarning(status_title,
        f"⚠️ 未找到可用的 R 环境 (No available R environment found)\n\n"
        f"系统 R_HOME: {system_r_home}\n\n"
        "椭圆拟合将降级为纯 Python 方式（最小二乘 / scipy LMG），仍可正常使用，但精度可能略低于 R-LMG 算法。\n"
        "(Ellipse fitting will fall back to pure Python mode (Least Squares / scipy LMG). Still usable, but accuracy may be slightly lower than R-LMG algorithm.)\n\n"
        "如需恢复 R 拟合，请：(To restore R fitting, please:)\n"
        "• 安装 R 并配置 R_HOME 环境变量，或 (• Install R and configure R_HOME environment variable, or)\n"
        f"• 将 R_Dist 文件夹放在: {base_path} (• Place R_Dist folder at: {base_path})")

# =========================================================================
# 步骤 2：立即执行配置 (必须在 import rpy2 之前！)
# =========================================================================
setup_portable_r()

# =========================================================================
# 步骤 3：设置语言环境 (解决中文乱码)
# =========================================================================
os.environ["LC_ALL"] = "C"
# 必须在设置完环境变量后，再 import rpy2
# [升级适配 Python 3.14] rpy2 为可选依赖：装得上则用 R-LMG 高精度拟合；
# 装不上（如目标 Python 暂无对应 wheel，或 R 运行时缺失）则降级为
# 纯 Python/scipy 拟合，不阻塞主程序启动。fitting_call_R_conicfit
# 在不可用时直接返回 None，由 fitting() 的多阶段策略走代数法兜底。
_RPY2_AVAILABLE = False
robjects = None
importr = None
numpy2ri = None
conversion = None
try:
    import rpy2.robjects as robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import numpy2ri
    # 激活 numpy 到 R 矩阵的自动转换
    from rpy2.robjects import conversion
    _RPY2_AVAILABLE = True
except Exception as _rpy2_err:
    print(f"[fittinglib] rpy2 不可用，椭圆拟合将使用纯 Python 方式 (scipy/最小二乘)。原因: {_rpy2_err}")
import math
import numpy as np
from pyenvelope import get_minimum_bounding_rectangle # [新增] MBR支持


def re_ellipse_fitting(points, eps=1e-12):
    """
    使用 SVD 伪逆 手搓最小二乘（数值稳定）
    等价于 numpy.linalg.lstsq 的数学本质

    模型:  x^2 + Bxy + Cy^2 + Dx + Ey + F = 0
    即:    A @ can_shu ≈ b
    """

    # --- 1. 数据准备 ---
    pts = np.array(points, dtype=float)
    x = pts[:, 0]
    y = pts[:, 1]

    # --- 2. 构建设计矩阵 A (N, 5) ---
    # 每一行: [xy, y^2, x, y, 1]
    A = np.column_stack([
        x * y,          # B
        y ** 2,         # C
        x,              # D
        y,              # E
        np.ones_like(x) # F
    ])

    # --- 3. 构建目标向量 b (N, 1) ---
    b = -(x ** 2)
    b = b.reshape(-1, 1)

    # --- 4. SVD 分解 A = U Σ V^T ---
    # full_matrices=False → 经济型 SVD，更高效
    U, S, Vt = np.linalg.svd(A, full_matrices=False)

    # --- 5. 构造 Σ^+（奇异值的伪逆）---
    # 对非常小的奇异值进行截断，防止数值爆炸
    S_inv = np.zeros_like(S)
    for i in range(len(S)):
        if S[i] > eps:
            S_inv[i] = 1.0 / S[i]
        else:
            S_inv[i] = 0.0

    # 也可以一行写完（等价）：
    # S_inv = np.where(S > eps, 1.0 / S, 0.0)

    # --- 6. 计算 A^+ = V Σ^+ U^T ---
    A_pinv = Vt.T @ np.diag(S_inv) @ U.T   # (5, N)

    # --- 7. 最小二乘解 can_shu = A^+ b ---
    can_shu = A_pinv @ b                   # (5, 1)

    return can_shu

#---------------------------------------------------------
# 代数→几何转换缓存（_algebraic_to_geometric_internal 使用）
_geo_cache = {}

#正规方程求解: (A^T A)^-1 A^T b
def _algebraic_to_geometric_internal(can_shu):
    """
    使用标准特征值方法将代数参数转换为几何参数
    参数: can_shu - [B, C, D, E, F] (A=1)
    返回: dict {'h', 'k', 'a', 'b', 'theta'}
    """
    can_shu = np.asarray(can_shu).flatten()
    key = tuple(can_shu.tolist())
    
    if key in _geo_cache:
        return _geo_cache[key]
    
    A, B, C, D, E, F = 1.0, can_shu[0], can_shu[1], can_shu[2], can_shu[3], can_shu[4]
    
    discriminant = B * B - 4 * A * C
    if discriminant >= 0:
        raise ValueError("方程不是椭圆")
    
    Q = np.array([[A, B / 2.0], [B / 2.0, C]], dtype=float)
    rhs = -np.array([D, E], dtype=float)
    
    try:
        h, k = np.linalg.solve(2 * Q, rhs)
    except np.linalg.LinAlgError:
        raise ValueError("无法求解中心，矩阵奇异")
    
    F_prime = A * h * h + B * h * k + C * k * k + D * h + E * k + F
    K = -F_prime
    
    eigvals, eigvecs = np.linalg.eigh(Q)
    lambda1, lambda2 = eigvals[0], eigvals[1]
    
    a = np.sqrt(K / lambda1)
    b = np.sqrt(K / lambda2)
    
    if a < b:
        a, b = b, a
        long_axis_vec = eigvecs[:, 1]
    else:
        long_axis_vec = eigvecs[:, 0]
    
    theta = np.arctan2(long_axis_vec[1], long_axis_vec[0])
    
    if theta > np.pi / 2.0:
        theta -= np.pi
    elif theta < -np.pi / 2.0:
        theta += np.pi
    
    result = {'h': h, 'k': k, 'a': a, 'b': b, 'theta': theta}
    _geo_cache[key] = result
    
    return result


#-----------------------辅助函数begin-test-----------------------------------------------

def mirror_points_180(points, centroid):
    """
    围绕形心镜像旋转180度，生成2n个点
    参数:
        points: 原始顶点集合
        centroid: 形心坐标（包含x, y属性）
    返回: 镜像后的2n个点（原始点 + 镜像点）
    """
    mirrored = []
    for p in points:
        mx = 2 * centroid.x - p[0]
        my = 2 * centroid.y - p[1]
        mirrored.append([mx, my])
    
    return points + mirrored


def algebraic_to_geometric(can_shu):
    """
    将代数参数 can_shu 转换为几何参数数组
    返回: [cx, cy, a, b, theta]（与R返回格式一致）
    """
    r = _algebraic_to_geometric_internal(can_shu)
    return np.array([r['h'], r['k'], r['a'], r['b'], r['theta']])


def calculate_teacher_pargini(points, init_points=None):
    """
    [新增/修改]  ParGini 初始值计算逻辑
    参数：
        points: 真实顶点集合（用于计算重心）
        init_points: 可选，用于计算MBR的点集合。如果为None，使用points计算MBR
    参数顺序：
    1. 多边形重心 X (并非矩形中心)
    2. 多边形重心 Y
    3. 矩形长 * 0.5 (半长轴)
    4. 矩形宽 * 0.5 (半短轴)
    5. 矩形长轴与 X 轴夹角
    """
    pts_np = np.array(points)

    # --- A. 计算多边形重心 (使用真实顶点) ---
    poly_centroid = np.mean(pts_np, axis=0)
    center_x = poly_centroid[0]
    center_y = poly_centroid[1]

    # --- B. 计算最小外接矩形 (MBR) ---
    # 如果提供了init_points（如镜像点），使用init_points计算MBR
    mbr_input = init_points if init_points is not None else points
    mbr_points = get_minimum_bounding_rectangle(mbr_input)
    mbr = np.array(mbr_points)

    # 取前三个点计算边向量 (p0->p1, p1->p2)
    p0, p1, p2 = mbr[0], mbr[1], mbr[2]
    vec1 = p1 - p0
    vec2 = p2 - p1
    len1 = np.linalg.norm(vec1)
    len2 = np.linalg.norm(vec2)

    # --- C. 区分长宽与角度 ---
    if len1 >= len2:
        rect_len = len1
        rect_width = len2
        # 计算长边向量的角度 (atan2 返回弧度)
        angle = math.atan2(vec1[1], vec1[0])
    else:
        rect_len = len2
        rect_width = len1
        angle = math.atan2(vec2[1], vec2[0])

    # --- D. 组装参数 ---
    # LMG 算法需要的 a, b 为半轴长，所以乘 0.5
    init_a = rect_len * 0.5
    init_b = rect_width * 0.5

    # 返回列表 [Cx, Cy, a, b, theta]
    return [center_x, center_y, init_a, init_b, angle]


def check_ls_quality(can_shu, points, poly_area):
    """
    [修改] 按照新规则判断拟合质量，仅支持几何参数格式
    参数:
        can_shu: 几何参数 [cx, cy, a, b, theta]
    返回: True (合格), False (异常)
    """
    if can_shu is None:
        return False

    try:
        can_shu_flat = np.asarray(can_shu).flatten()
        if len(can_shu_flat) != 5:
            return False

        a = float(can_shu_flat[2])
        b = float(can_shu_flat[3])

        if math.isnan(a) or math.isnan(b) or a <= 0 or b <= 0:
            return False

        ellipse_area = math.pi * a * b
        if poly_area > 0:
            ratio = ellipse_area / poly_area
            if ratio < 1.0 or ratio > 3.0:
                return False

        return True

    except Exception as e:
        print(f"校验过程出错: {e}")
        return False

def fitting_call_R_conicfit(points, init_points=None):
    """
    [修改后] 调用 R 语言 fit.ellipseLMG，使用老师指定的 ParGini 初值
    参数：
        points: 真实顶点集合（用于拟合）
        init_points: 可选，用于计算MBR的点集合（如镜像点）。如果为None，使用points计算MBR
    返回：
        几何参数数组 [cx, cy, a, b, theta]，R失败时返回 None（不回退）
    """
    if not _RPY2_AVAILABLE:
        return None
    try:
        pts_np = np.array(points)

        # 使用 init_points（如镜像点）计算MBR，真实顶点计算重心
        pargini_list = calculate_teacher_pargini(points, init_points)

        # 构造 R 对象
        pargini_vec = robjects.FloatVector(pargini_list)
        par_gini_r = robjects.r.matrix(pargini_vec, ncol=1)

        # 导入包
        conicfit = importr('conicfit')

        # 使用局部转换器调用 R
        with conversion.localconverter(robjects.default_converter + numpy2ri.converter):
            res_r = conicfit.fit_ellipseLMG(pts_np, par_gini_r, 1e-5)

            geo_params_r = res_r[0]
            geo_params = np.array(geo_params_r).flatten()

        # 直接返回几何参数，不再转换为代数参数
        return geo_params

    except Exception as e:
        print(f"R语言接口调用失败: {e}")
        return None

def fitting(points, center_point, area):
    """
    [重写 V13.0] 三阶段拟合策略，返回几何参数
    策略：
    1. n >= 5：
       - round1: 代数法，基于真实点（n个）
       - round2: 代数法，基于真实点+镜像点（2n个）
       - round3: R LMG，基于真实点+镜像点（2n个），用2n点的MBR做初值
    2. n < 5：
       - round1: R LMG，基于真实点（n个），用真实点MBR做初值
       - round2: 代数法，基于真实点+镜像点（2n个）
       - round3: R LMG，基于真实点+镜像点（2n个），用2n点的MBR做初值
    3. 面积比合格范围：1.0 ~ 3.0
    4. 每轮失败直接进下一轮，不回退
    5. 三轮全失败：用真实点重心(cx,cy) + 2n点MBR(半长轴/半短轴/倾角)作为椭圆保底
    6. 返回值：统一为几何参数 [cx, cy, a, b, theta]
    """
    n_sides = len(points)
    final_geo = None

    mirrored_points = mirror_points_180(points[:], center_point)

    if n_sides >= 5:
        # --- n >= 5：代数法真实点 → 代数法2n点 → R LMG 2n点(2n初值) ---

        # Round 1: 代数法，真实点
        try:
            can_shu = re_ellipse_fitting(points)
            geo_round1 = algebraic_to_geometric(can_shu)
        except Exception as e:
            print(f"round1 代数拟合失败: {e}")
            geo_round1 = None

        if geo_round1 is not None and check_ls_quality(geo_round1, points, area):
            final_geo = geo_round1
        else:
            # Round 2: 代数法，2n点
            try:
                can_shu = re_ellipse_fitting(mirrored_points)
                geo_round2 = algebraic_to_geometric(can_shu)
            except Exception as e:
                print(f"round2 代数拟合失败: {e}")
                geo_round2 = None

            if geo_round2 is not None and check_ls_quality(geo_round2, points, area):
                final_geo = geo_round2
            else:
                # Round 3: R LMG，2n点，2n点MBR初值
                try:
                    geo_round3 = fitting_call_R_conicfit(mirrored_points, init_points=mirrored_points)
                except Exception as e:
                    print(f"round3 R LMG拟合失败: {e}")
                    geo_round3 = None

                if geo_round3 is not None and check_ls_quality(geo_round3, points, area):
                    final_geo = geo_round3
                # Round 3 失败或不合格，直接进保底

    else:
        # --- n < 5：R LMG真实点 → 代数法2n点 → R LMG 2n点(2n初值) ---

        # Round 1: R LMG，真实点，真实点MBR初值
        try:
            geo_round1 = fitting_call_R_conicfit(points)
        except Exception as e:
            print(f"round1 R LMG拟合失败: {e}")
            geo_round1 = None

        if geo_round1 is not None and check_ls_quality(geo_round1, points, area):
            final_geo = geo_round1
        else:
            # Round 2: 代数法，2n点
            try:
                can_shu = re_ellipse_fitting(mirrored_points)
                geo_round2 = algebraic_to_geometric(can_shu)
            except Exception as e:
                print(f"round2 代数拟合失败: {e}")
                geo_round2 = None

            if geo_round2 is not None and check_ls_quality(geo_round2, points, area):
                final_geo = geo_round2
            else:
                # Round 3: R LMG，2n点，2n点MBR初值
                try:
                    geo_round3 = fitting_call_R_conicfit(mirrored_points, init_points=mirrored_points)
                except Exception as e:
                    print(f"round3 R LMG拟合失败: {e}")
                    geo_round3 = None

                if geo_round3 is not None and check_ls_quality(geo_round3, points, area):
                    final_geo = geo_round3
                # Round 3 失败或不合格，直接进保底

    # 最终保底：用真实点重心(cx,cy) + 2n点MBR(半长轴/半短轴/倾角)作为椭圆
    # MBR 为最终兜底，不再用 try-except 吞异常返回 None（否则下游 make_final_ellipse(None) 报更难定位的错）
    if final_geo is None:
        # calculate_teacher_pargini 返回 [cx, cy, a, b, theta]
        # points 用于计算重心，mirrored_points(2n点) 用于计算 MBR
        final_geo = calculate_teacher_pargini(points, mirrored_points)
        print("三轮拟合全不合格，使用 MBR 保底")

    return final_geo
