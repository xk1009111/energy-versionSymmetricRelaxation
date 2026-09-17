#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
batch_lam_net.py —— λ × 初始网络类型 批量退火数据生成（能量版）

数据文件 = 软件「输出数据」按钮输出的三种 xlsx（与 GUI 导出完全同构，
经 exportUtils.create 生成）：
  - {label}_{state}_ellipse.xlsx    每细胞 1 行（椭圆拟合 + n/P/A/mn/层数）
  - {label}_{state}_edgeAngle.xlsx  每顶角 1 行（内角 + 夹边长 + 层数）
  - {label}_{state}_ME_MA.xlsx      每边缘细胞 1 行（ME/MA + 层数）

实验矩阵（2026-09-13 方案）：
  - λ ∈ {0, 1, 10}（显式 k_angle，k_center=1）
  - 初始网络（目标 ~400 细胞）：
      grid02  : 扰动维诺图 getCells(k=0.2, param2=0, n=20) → 400 细胞
      grid1   : 扰动维诺图 getCells(k=1.0, param2=0, n=20) → ~397 细胞
      random  : 均匀随机 getCellsRandom(n_seeds=400) → ~334 细胞（GUI 输入习惯）
  - μ = 0.1（normalize_L0=True：有效步 = μ·L0²，每轮位移 ≈ μ·L0，无量纲）
  - 每组 10 seeds（rep=0..9），共 3×3×10 = 90 次运行
  - 每次固定 100 轮；同 rep 跨 λ 共享同一初始网络（配对比较 Δ）
  - 退火器配置与 GUI 完全一致：边缘顶点参与（freeze_marginal=False）、
    凸性守卫开、串行降序扫描（serial_descending=True）、L0 归一化开

输出布局（每 run 一个子文件夹，初始态 + 终态各 3 个 xlsx）：
  data/{label}/  其中 label = {net}_L{λ}_rep{r}，如 grid02_L0_rep0

实现要点（沿用 pilot_lambda_mu.py 框架）：
  - tkinter 桩（无界面环境导入 utillib.mylib 所需）
  - 导出链与 GUI「输出数据」按钮逐行同构：
        cellData.list_line_of_cell() → 逐细胞 like_ellipse() →
        exportUtils.create("export", cells, lineOfCell, currentTimes) →
        move 重命名（GUI: {base}_ellipse/_edgeAngle/_ME_MA.xlsx）
  - 并行安全：exportUtils 输出文件名用秒级时间戳，10 进程并行必然撞车，
    改用全局唯一 run-id 作 currentTimes（init=run_id*2, final=run_id*2+1）
  - 退火梯度只依赖顶点坐标，循环内不做全量 flush（初/终各一次）
  - 按 rep 多进程并行（Windows spawn，模块顶层完成 chdir/sys.path）
"""
import os
import sys
import time
import copy
import random
import shutil
import warnings
import traceback
import multiprocessing as mp

import numpy as np

warnings.filterwarnings("ignore")
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())


def _install_tkinter_stub():
    """无界面（无显示器）环境里 tkinter 会卡死；批量脚本不需要 GUI，
    用无操作桩替换 tkinter.messagebox，仅让 fittinglib 的导入期
    showinfo/showwarning 变成空操作。不影响任何数值计算。"""
    import types
    _tk = types.ModuleType("tkinter")
    _mb = types.ModuleType("tkinter.messagebox")
    _noop = lambda *a, **k: None
    for fn in ("showinfo", "showwarning", "showerror",
               "askyesno", "askokcancel", "askquestion"):
        setattr(_mb, fn, _noop)
    _tk.messagebox = _mb
    sys.modules.setdefault("tkinter", _tk)
    sys.modules.setdefault("tkinter.messagebox", _mb)


_install_tkinter_stub()

import initVoronoi as Viro
from cell.CellData import CellData
import annealing.energy as E
import utillib.exportUtils as exportUtils
import utillib.layerMarker as layerMarker
from annealing.physical_annealer import PhysicalAnnealer, PhysicalAnnealerParams


# ========================== 实验配置 ==========================
N_GRID = 20            # 扰动网格边长（k=0.2 → 400 细胞）
N_RANDOM_SEEDS = 400   # 均匀随机种子点数（预处理后 ~334 细胞）
NETWORKS = [           # (network, k)；k 仅 grid 模式使用
    ("grid", 0.2),
    ("grid", 1.0),
    ("random", None),
]
LAMBDAS = [100.0]
MU = 0.1
ITERS = 100            # 每次运行固定退火轮数
REPS = 10              # 每组种子数
MARGINAL = 1           # 边缘顶点参与退火（freeze_marginal=False）
PARALLEL = True
WORKERS = 10           # 18 核机器：10 个 rep 各占一进程（每 rep 内 9 次运行串行）

OUT_DIR = os.path.join(os.getcwd(), "data")
os.makedirs(OUT_DIR, exist_ok=True)

# 能量测量统一用裸和（k_center=1, k_angle=1），仅用于控制台进度日志
SP_REPORT = dict(k_center=1.0, k_angle=1.0, centroid_mode="vertex")


# ========================== 工具函数 ==========================
def netkey(network, k):
    return {"grid": f"grid{str(k).replace('.', '')}"}.get(network, network)


def make_initial(network, k, seed):
    """用固定种子生成初始网络（同 rep 跨 λ 共享；grid 两 k 值共享同一
    扰动随机流——rand_in_unit_hex 的抽样与 k 无关，k 仅缩放偏移量）。"""
    random.seed(seed)
    np.random.seed(seed % (2 ** 32))
    if network == "random":
        cells = Viro.getCellsRandom(N_RANDOM_SEEDS)
    else:
        cells = Viro.getCells(k, 0, N_GRID)
    for c in cells:
        for i in range(len(c.points)):
            pt = c.points[i]
            c.points[i] = pt.tolist() if hasattr(pt, "tolist") else list(pt)
    return cells


def export_state(cells, cd, label, state, times):
    """复刻 GUI「输出数据」按钮的导出链，输出到 data/{label}/。

    与 energy_annealing_main._export_data 逐行对应：
        逐细胞 like_ellipse() → list_line_of_cell() →
        exportUtils.create("export", cells, lineOfCell, currentTimes) →
        move 重命名为 {label}_{state}_{ellipse|edgeAngle|ME_MA}.xlsx
    times 必须全局唯一（并行防撞车），由调用方用 run-id 生成。
    """
    out_dir = os.path.join(OUT_DIR, label)
    os.makedirs(out_dir, exist_ok=True)

    # --- GUI 导出准备（_export_data 的 try 块） ---
    cd.list_line_of_cell()
    for c in cells:
        c.like_ellipse()

    # --- exportUtils.create（与 GUI 同参：excelName="export"） ---
    if not exportUtils.create("export", cells, cd.lineOfCell, currentTimes=times):
        raise RuntimeError("exportUtils.create returned False")

    # --- move 重命名（与 GUI 同构，文件名按实验设定） ---
    pairs = [
        (f"ellipse_{times}.xlsx", f"{label}_{state}_ellipse.xlsx"),
        (f"edgeAngle_{times}.xlsx", f"{label}_{state}_edgeAngle.xlsx"),
        (f"export_ME_MA_{times}.xlsx", f"{label}_{state}_ME_MA.xlsx"),
    ]
    for src, dst in pairs:
        if not os.path.exists(src):
            raise FileNotFoundError(f"导出产物缺失: {src}")
        shutil.move(src, os.path.join(out_dir, dst))


def cleanup_stray(times):
    """导出中途失败时，清理工作目录下可能残留的产物文件。"""
    for name in (f"ellipse_{times}.xlsx", f"edgeAngle_{times}.xlsx",
                 f"export_ME_MA_{times}.xlsx"):
        try:
            if os.path.exists(name):
                os.remove(name)
        except OSError:
            pass


# ========================== 单次运行 ==========================
def run_single(init_cells, network, k, lam, rep, run_id, seed):
    """一次运行：初态导出 → ITERS 轮退火 → 终态导出。

    run_id ∈ [0, 90)：全局唯一，currentTimes 用 run_id*2（初态）/ *2+1（终态）。
    返回 (label, E_init, E_final, elapsed_sec)（仅供控制台日志）。
    """
    nk = netkey(network, k)
    label = f"{nk}_L{lam:g}_rep{rep}"
    t0 = time.time()

    cells = copy.deepcopy(init_cells)
    cd = CellData(cells)
    cd.flush()

    # GUI 导出前必调：layer_mark2 从边缘细胞（层 1）BFS 扩散，
    # 给所有非边缘细胞赋予正确的层号（2, 3, ...）。
    # 若漏掉，非边缘细胞保持层 0（flush → setting_layer 只标边缘细胞为层 1）。
    # 对应 GUI 第 857 行（初始态）与第 1199 行（终态）。
    layerMarker.layer_mark2(cells, seed)

    sp_run = E.SymParams(k_center=1.0, k_angle=lam, centroid_mode="vertex")
    sp_rep = E.SymParams(**SP_REPORT)

    # ---- 初态导出 ----
    t_init = run_id * 2
    try:
        export_state(cells, cd, label, "init", t_init)
    except Exception:
        cleanup_stray(t_init)
        raise

    e0 = E.total_energy(cells, sp_rep)

    # ---- 退火（配置与 GUI Annealer 完全一致）----
    pa = PhysicalAnnealer(PhysicalAnnealerParams(
        mobility=MU,
        k_center=1.0,
        k_angle=lam,
        centroid_mode="vertex",
        convex_guard=True,
        freeze_marginal=not bool(MARGINAL),
        normalize_L0=True,       # GUI 有效默认：有效步 = μ·L0²
        serial_descending=True,  # GUI 有效默认：按 |F| 降序串行扫描
    ))
    for _ in range(ITERS):
        pa.one_round(cells, sp_run, eta=MU)

    # ---- 终态导出 ----
    cd.flush()
    # 退火后拓扑/几何改变，层号可能变化；GUI 在终态导出前再次调 layer_mark2
    # （对应第 1199 行），确保边缘细胞判定与 BFS 扩散基于退火后的网络。
    layerMarker.layer_mark2(cells, seed)
    e1 = E.total_energy(cells, sp_rep)
    t_fin = run_id * 2 + 1
    try:
        export_state(cells, cd, label, "final", t_fin)
    except Exception:
        cleanup_stray(t_fin)
        raise

    return label, float(e0.total), float(e1.total), round(time.time() - t0, 2)


# ========================== 单 rep 流程 ==========================
def run_rep(rep):
    """跑完某个 rep 下全部 (网络 × λ) 组合，共 9 次。"""
    seed = 1000 + rep * 100
    total = len(NETWORKS) * len(LAMBDAS)
    done = 0
    for net_i, (network, k) in enumerate(NETWORKS):
        init_cells = make_initial(network, k, seed)
        for lam_i, lam in enumerate(LAMBDAS):
            done += 1
            run_id = rep * (len(NETWORKS) * len(LAMBDAS)) + net_i * len(LAMBDAS) + lam_i
            nk = netkey(network, k)
            label = f"{nk}_L{lam:g}_rep{rep}"
            print(f"[rep{rep} {done}/{total}] {label} "
                  f"(n_init={len(init_cells)})", flush=True)
            try:
                lb, e0, e1, dt = run_single(init_cells, network, k, lam, rep, run_id, seed)
                print(f"[rep{rep}] {lb} done "
                      f"E_total {e0:.3f}->{e1:.3f} ({dt}s)", flush=True)
            except Exception as e:
                traceback.print_exc()
                print(f"[rep{rep}] {label} FAILED: {e}", flush=True)
    return rep


# ========================== 主流程 ==========================
def run():
    t_start = time.time()
    if PARALLEL and REPS > 1:
        workers = min(WORKERS, REPS)
        with mp.Pool(workers) as pool:
            for _ in pool.imap_unordered(run_rep, range(REPS)):
                print(f"[main] 一个 rep 完成 "
                      f"(累计耗时 {time.time() - t_start:.0f}s)", flush=True)
    else:
        for rep in range(REPS):
            run_rep(rep)

    # 清点产物：90 个子文件夹 × 6 个 xlsx
    n_dirs = 0
    n_files = 0
    for name in sorted(os.listdir(OUT_DIR)):
        sub = os.path.join(OUT_DIR, name)
        if os.path.isdir(sub):
            n_dirs += 1
            n_files += sum(1 for f in os.listdir(sub) if f.endswith(".xlsx"))
    print(f"\n完成：{n_dirs} 个 run 文件夹 / {n_files} 个 xlsx"
          f"（期望 {REPS * len(NETWORKS) * len(LAMBDAS)} 文件夹 / "
          f"{REPS * len(NETWORKS) * len(LAMBDAS) * 6} 文件），"
          f"总耗时 {time.time() - t_start:.0f}s", flush=True)
    print(f"输出目录：{OUT_DIR}", flush=True)


if __name__ == "__main__":
    run()
