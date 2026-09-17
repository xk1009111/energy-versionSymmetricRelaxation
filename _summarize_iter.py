# -*- coding: utf-8 -*-
"""
_summarize_iter.py — 批量退火数据汇总（纯 openpyxl + numpy）。

数据结构: data_fixed/{net}/L{lam}/run{r}/{iter}_{phase}_{type}.xlsx
  net ∈ {grid02, grid10, random}, lam ∈ {0, 1, 10}, r ∈ 1..10, iter ∈ {0, 1}

输出: 3 个 {net}_iteration_汇总.xlsx + 1 个 iteration_C1_C2_汇总对比.xlsx
每组汇总内含 3 个 λ (L0/L1/L10) 的 init/final/Δ 并列对比。

分类规则:
  - mn 汇总三分类: MC=层1, NC-MC=层2, IC=层≥3
  - 其他汇总二分类: MC=层1, IC=层≥2
内角保持度数。
相同随机数下 init 态拓扑/几何量与 k、λ 无关，统计时合并 4λ×全部 runs；final 态按 λ 分列。
"""
import os, sys, re
import numpy as np
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

# ---------- 样式常量 ----------
THIN = Side(style='thin', color='000000')
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER = Alignment(horizontal='center', vertical='center', wrap_text=True)
HEADER_FONT = Font(bold=True, size=11)
TITLE_FONT = Font(bold=True, size=12)
HEADER_FILL = PatternFill(start_color='B8CCE4', end_color='B8CCE4', fill_type='solid')
GROUP_FILL = PatternFill(start_color='D5E8D4', end_color='D5E8D4', fill_type='solid')
DELTA_FILL = PatternFill(start_color='F8CECC', end_color='F8CECC', fill_type='solid')
LENGTH_LABELS = ['<0.1', '0.1-0.3', '0.3-0.5', '0.5-0.7', '0.7-0.9', '0.9-1.1', '1.1-1.3', '1.3-1.5', '1.5-1.7', '1.7-1.9', '>1.9']
ANGLE_LABELS = ['<30°', '30-50°', '50-70°', '70-90°', '90-110°', '110-130°', '130-150°', '150-170°', '170-180°']
LAMS = ['L0', 'L1', 'L10', 'L100']
ITERATIONS = [0, 1]
PHASE_MAP = {0: 'initial', 1: 'final'}
DELTA_PAIRS = [('final-initial', 1, 0)]

def style_header(cell):
    cell.font = HEADER_FONT; cell.fill = HEADER_FILL
    cell.alignment = CENTER; cell.border = BORDER

def length_bin(v):
    bins = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, float('inf')]
    for i, b in enumerate(bins):
        if v <= b: return i
    return len(bins) - 1

def angle_bin(v):
    bins = [30, 50, 70, 90, 110, 130, 150, 170, 180]
    for i, b in enumerate(bins):
        if v <= b: return i
    return len(bins) - 1

# ---------- 发现 ----------
def discover_nets(base):
    return sorted([d for d in os.listdir(base)
                   if os.path.isdir(os.path.join(base, d))
                   and any(os.path.isdir(os.path.join(base, d, l)) for l in LAMS)])

# ---------- 读取 ----------
def read_xlsx_rows(path):
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    if not rows: return []
    header = list(rows[0])
    return [dict(zip(header, row)) for row in rows[1:]]

def collect(base, net, lam, runs, iter_num, phase, typ):
    frames = []
    for run in runs:
        path = os.path.join(base, net, lam, f'run{run}', f'{iter_num}_{phase}_{typ}.xlsx')
        if not os.path.exists(path): continue
        for row in read_xlsx_rows(path):
            row['run'] = run
            frames.append(row)
    return frames

# ---------- 分类 ----------
def mc_ic(rows):
    """二分类: MC=层1, IC=层≥2"""
    mc = [r for r in rows if int(r.get('当前细胞层', 0)) == 1]
    ic = [r for r in rows if int(r.get('当前细胞层', 0)) >= 2]
    return mc, ic

def three_way(rows):
    """三分类: MC=层1, NC-MC=层2, IC=层≥3（用于 mn 汇总）"""
    mc = [r for r in rows if int(r.get('当前细胞层', 0)) == 1]
    ncmc = [r for r in rows if int(r.get('当前细胞层', 0)) == 2]
    ic3 = [r for r in rows if int(r.get('当前细胞层', 0)) >= 3]
    return mc, ncmc, ic3

def mean_std(vals):
    vals = [float(v) for v in vals if v is not None and not (isinstance(v, float) and np.isnan(v))]
    if len(vals) == 0: return None, None
    return float(np.mean(vals)), float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0

def group_by_n(rows, n_col='细胞边数'):
    d = {}
    for r in rows:
        n = int(r.get(n_col, 0))
        d.setdefault(n, []).append(r)
    return d

# ---------- 均值|标准差交替写入 ----------
def _write_cell(ws, r, c, val, fmt=None):
    cell = ws.cell(row=r, column=c, value=val)
    cell.border = BORDER
    if fmt and isinstance(val, (int, float)): cell.number_format = fmt
    return cell

def _write_mean_std_alt(ws, start_r, title, groups_dict, col_key, fmt='0.0000'):
    """均值|标准差交替列序。groups_dict 键为组名，值为 rows 列表。"""
    ws.cell(row=start_r, column=1, value=title).font = TITLE_FONT
    r = start_r + 1
    headers = ['组']
    for gname in groups_dict.keys():
        headers += [f'{gname}_均值', f'{gname}_标准差']
    for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j + 1, value=h))
    r += 1
    for gname, rows in groups_dict.items():
        m, s = mean_std([row.get(col_key) for row in rows])
        _write_cell(ws, r, 1, gname)
        col = 2
        for gn2, rows2 in groups_dict.items():
            m2, s2 = mean_std([row.get(col_key) for row in rows2])
            _write_cell(ws, r, col, m2, fmt)
            _write_cell(ws, r, col + 1, s2, fmt)
            col += 2
        break  # 只写第一行（因为每行一个组，列遍历所有组）
    # 上面逻辑有问题：每行是一个"行组"（MC/IC等），列是所有 gname 的均值|标准差
    # 重写：
    ws.delete_rows(start_r + 1, 1)
    r = start_r + 1
    row_labels = list(groups_dict.keys())
    for label in row_labels:
        _write_cell(ws, r, 1, label)
        col = 2
        for gname, rows in groups_dict.items():
            m, s = mean_std([row.get(col_key) for row in rows])
            _write_cell(ws, r, col, m, fmt)
            _write_cell(ws, r, col + 1, s, fmt)
            col += 2
        r += 1
    return r + 1

def _write_dist_alt(ws, start_r, title, groups_dict, col_key, bin_fn, labels, is_angle=False):
    """分布统计，均值|标准差交替。"""
    ws.cell(row=start_r, column=1, value=title).font = TITLE_FONT
    r = start_r + 1
    n_bins = len(labels)
    for gname, rows in groups_dict.items():
        ws.cell(row=r, column=1, value=gname).font = Font(bold=True)
        r += 1
        headers = ['统计量'] + labels + ['总数'] + ['均值'] + ['标准差']
        for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j+1, value=h))
        r += 1
        vals_raw = [row.get(col_key) for row in rows if row.get(col_key) is not None]
        if is_angle:
            vals = [float(v) * 180.0 / np.pi if isinstance(v, float) and v < 6.5 else float(v)
                    for v in vals_raw]
        else:
            vals = [float(v) for v in vals_raw]
        bins = [0] * n_bins
        for v in vals:
            idx = bin_fn(v)
            if 0 <= idx < n_bins: bins[idx] += 1
        total = len(vals)
        for row_label, data in [('数量', bins)]:
            _write_cell(ws, r, 1, row_label)
            for j, v in enumerate(data): _write_cell(ws, r, 2 + j, v)
            _write_cell(ws, r, 2 + n_bins, total)
            m, s = mean_std(vals)
            _write_cell(ws, r, 2 + n_bins + 1, m, '0.0000')
            _write_cell(ws, r, 2 + n_bins + 2, s, '0.0000')
            r += 1
        # 频率行（均值|标准差交替不适用于频率，频率直接放数量右边）
        freq = [f'{b} ({b/total:.1%})' if total else '0' for b in bins]
        _write_cell(ws, r, 1, '频率')
        for j, v in enumerate(freq): _write_cell(ws, r, 2 + j, v)
        _write_cell(ws, r, 2 + n_bins, total)
        r += 1
        r += 1
    return r + 1

def _write_delta_block_c(ws, start_r, title, all_pairs_data):
    ws.cell(row=start_r, column=1, value=title).font = TITLE_FONT
    r = start_r + 1
    delta_metrics = ['ΔP', 'ΔA', 'ΔSI', 'ΔEavg']
    pair_names = [p[0] for p in all_pairs_data]
    for (mc_label, rows_list) in [('MC', [p[1] for p in all_pairs_data]),
                                   ('IC', [p[2] for p in all_pairs_data])]:
        for metric in delta_metrics:
            ws.cell(row=r, column=1, value=f'{mc_label} - {metric}').font = Font(bold=True)
            r += 1
            all_n = set()
            for rows in rows_list:
                for row in rows:
                    n = row.get('n')
                    if n is not None: all_n.add(int(n))
            all_n = sorted(all_n)
            headers = ['n值']
            for pn in pair_names: headers += [f'{pn}_均值', f'{pn}_标准差']
            for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j+1, value=h))
            r += 1
            for nv in all_n:
                _write_cell(ws, r, 1, nv)
                col = 2
                for pi, rows in enumerate(rows_list):
                    sub = [row for row in rows if int(row.get('n', -1)) == nv]
                    m, s = mean_std([row.get(metric) for row in sub])
                    _write_cell(ws, r, col, m, '0.0000')
                    _write_cell(ws, r, col + 1, s, '0.0000')
                    col += 2
                r += 1
            r += 1
    return r + 1

# ---------- Delta / C1C2 ----------
def compute_delta_rows(from_rows, to_rows):
    from_map = {(int(r['run']), int(r['细胞序号'])): r for r in from_rows}
    to_map = {(int(r['run']), int(r['细胞序号'])): r for r in to_rows}
    keys = set(from_map.keys()) & set(to_map.keys())
    out = []
    for key in keys:
        fr, tr = from_map[key], to_map[key]
        Pf = float(fr['细胞周长']); Pt = float(tr['细胞周长'])
        Af = float(fr['细胞面积']); At = float(tr['细胞面积'])
        nf = int(fr.get('细胞边数', fr.get('边数', 0)))
        layer = int(fr.get('当前细胞层', 0))
        SI_f = Pf / np.sqrt(max(Af, 1e-12))
        SI_t = Pt / np.sqrt(max(At, 1e-12))
        out.append({
            'run': key[0], 'No.': key[1], 'n': nf,
            'ΔP': Pt - Pf, 'ΔA': At - Af,
            'ΔSI': SI_t - SI_f, 'ΔEavg': (Pt - Pf) / max(nf, 1),
            '当前细胞层': layer, 'A_from': Af,
        })
    return out

def fit_c1_c2(delta_rows):
    results = {}
    for (label, n_offset, mask_fn) in [
        ('MC', 5, lambda r: int(r.get('当前细胞层', 0)) == 1),
        ('IC', 6, lambda r: int(r.get('当前细胞层', 0)) >= 2)
    ]:
        sub = [r for r in delta_rows if mask_fn(r)]
        if not sub:
            results[label] = {'C1': 0, 'C2': 0, 'RMSE': 0, 'n_cells': 0, 'An': {}, 'details': []}
            continue
        An = {}
        for r in sub:
            nv = int(r['n'])
            An.setdefault(nv, []).append(float(r['A_from']))
        An = {k: float(np.mean(v)) for k, v in An.items()}
        n_arr = np.array([int(r['n']) for r in sub], dtype=float)
        da_arr = np.array([float(r['ΔA']) for r in sub], dtype=float)
        a_arr = np.array([float(r['A_from']) for r in sub], dtype=float)
        dn = n_arr - n_offset
        mask = dn != 0
        c1 = float(np.sum(dn[mask] * da_arr[mask]) / np.sum(dn[mask] ** 2)) if mask.any() else 0.0
        correction = np.array([1 - a_arr[i] / An.get(int(n_arr[i]), 1.0) for i in range(len(sub))])
        best_c2, best_rmse = 0.0, float('inf')
        for c2_val in np.arange(0.0, 3.001, 0.01):
            pred = c1 * dn + c2_val * correction
            rmse = float(np.sqrt(np.mean((da_arr - pred) ** 2)))
            if rmse < best_rmse:
                best_rmse, best_c2 = rmse, float(c2_val)
        results[label] = {'C1': c1, 'C2': best_c2, 'RMSE': best_rmse,
                           'n_cells': len(sub), 'An': An, 'details': sub}
    return results

# ---------- 写入原始数据块 ----------
def write_raw_block(ws, titles, data_lists, headers):
    n = len(titles); cols_per = len(headers)
    for i in range(n):
        start_col = 1 + i * cols_per
        ws.merge_cells(start_row=1, start_column=start_col,
                       end_row=1, end_column=start_col + cols_per - 1)
        c = ws.cell(row=1, column=start_col, value=titles[i])
        c.font = HEADER_FONT; c.fill = DELTA_FILL; c.alignment = CENTER
        for cc in range(start_col, start_col + cols_per):
            ws.cell(row=1, column=cc).fill = DELTA_FILL
        for j, h in enumerate(headers):
            style_header(ws.cell(row=2, column=start_col + j, value=h))
        rows = data_lists[i]
        for ri, row_data in enumerate(rows):
            for j, h in enumerate(headers):
                val = row_data.get(h)
                c = ws.cell(row=3 + ri, column=start_col + j, value=val)
                c.border = BORDER
                if isinstance(val, float): c.number_format = '0.000000'
    ws.freeze_panes = 'A3'

def write_delta_detail_block(ws, all_pairs_data):
    delta_metrics = ['ΔP', 'ΔA', 'ΔSI', 'ΔEavg']
    pair_names = [p[0] for p in all_pairs_data]
    n_pairs = len(pair_names); cols_per_group = n_pairs * len(delta_metrics)
    col = 1
    for label, idx in [('MC', 1), ('IC', 2)]:
        merged = {}
        for pi, (pair_name, mc_r, ic_r) in enumerate(all_pairs_data):
            target = mc_r if label == 'MC' else ic_r
            for r in target:
                key = (r.get('run'), r.get('No.'))
                merged.setdefault(key, {}).update(r)
                for dm in delta_metrics:
                    merged[key].update({f'{dm}_{pair_name}': r.get(dm)})
        sorted_rows = [merged[k] for k in sorted(merged.keys())]
        start_col = col
        ws.merge_cells(start_row=1, start_column=start_col,
                       end_row=1, end_column=start_col + cols_per_group - 1)
        c = ws.cell(row=1, column=start_col, value=label)
        c.font = HEADER_FONT; c.fill = DELTA_FILL; c.alignment = CENTER
        for cc in range(start_col, start_col + cols_per_group):
            ws.cell(row=1, column=cc).fill = DELTA_FILL
        for mi, metric in enumerate(delta_metrics):
            metric_start = start_col + mi * n_pairs
            ws.merge_cells(start_row=2, start_column=metric_start,
                           end_row=2, end_column=metric_start + n_pairs - 1)
            c2 = ws.cell(row=2, column=metric_start, value=metric)
            c2.font = HEADER_FONT; c2.fill = DELTA_FILL; c2.alignment = CENTER
            for cc in range(metric_start, metric_start + n_pairs):
                ws.cell(row=2, column=cc).fill = DELTA_FILL
            for j, pn in enumerate(pair_names):
                style_header(ws.cell(row=3, column=metric_start + j, value=pn))
                ws.cell(row=3, column=metric_start + j).fill = DELTA_FILL
            col = metric_start + n_pairs
        for ri, row_data in enumerate(sorted_rows):
            for mi, metric in enumerate(delta_metrics):
                for j, pn in enumerate(pair_names):
                    key = f'{metric}_{pn}'
                    val = row_data.get(key)
                    c = ws.cell(row=4 + ri, column=start_col + mi * n_pairs + j,
                                value=float(val) if val is not None else None)
                    c.border = BORDER
                    if val is not None: c.number_format = '0.000000'
    ws.freeze_panes = 'A4'

# ---------- 主处理（单网络） ----------
def process_net(base, net, all_c1):
    out_path = os.path.join(base, f'{net}_iteration_汇总.xlsx')
    runs = list(range(1, 11))

    # 读 3 λ × 2 iter
    ell, ea, mm = {}, {}, {}
    for lam in LAMS:
        for it in ITERATIONS:
            phase = PHASE_MAP[it]
            ell.setdefault(lam, {})[it] = collect(base, net, lam, runs, it, phase, 'ellipse')
            ea.setdefault(lam, {})[it] = collect(base, net, lam, runs, it, phase, 'edgeAngle')
            mm.setdefault(lam, {})[it] = collect(base, net, lam, runs, it, phase, 'ME_MA')

    # MC/IC 二分类 + three_way 三分类
    mc_ell, ic_ell, tway = {}, {}, {}
    mc_ea, ic_ea = {}, {}
    for lam in LAMS:
        mc_ell[lam], ic_ell[lam] = {}, {}
        mc_ea[lam], ic_ea[lam] = {}, {}
        for it in ITERATIONS:
            mc_ell[lam][it], ic_ell[lam][it] = mc_ic(ell[lam][it])
            mc_ea[lam][it], ic_ea[lam][it] = mc_ic(ea[lam][it])
        tway[lam] = three_way(ell[lam][0])  # mn 用 init 态

    # init 态合并池 (相同随机数, 与 k/λ 无关)
    ell_init_pool = []
    for lam in LAMS:
        ell_init_pool += ell[lam][0]
    tway_pool = three_way(ell_init_pool)

    # Delta + C1C2（final - initial）
    all_pairs_data = []
    c1_c2_results = {}
    for lam in LAMS:
        delta = compute_delta_rows(ell[lam][0], ell[lam][1])
        mc_d, ic_d = mc_ic(delta)
        all_pairs_data.append((lam, mc_d, ic_d))
        c1_c2_results[lam] = fit_c1_c2(delta)

    wb = Workbook(); wb.remove(wb.active)

    # ---- ellipse 原始 ----
    ws = wb.create_sheet('ellipse')
    titles_raw = []; data_raw = []
    for lam in LAMS:
        for it in ITERATIONS:
            titles_raw += [f'{lam}_{PHASE_MAP[it]}_MC', f'{lam}_{PHASE_MAP[it]}_IC']
            data_raw += [mc_ell[lam][it], ic_ell[lam][it]]
    write_raw_block(ws, titles_raw, data_raw,
                    ['run', '细胞序号', '长半轴', '短半轴', '相邻细胞边数和', '细胞边数', '细胞周长', '细胞面积', '当前细胞层'])

    # ---- delta_detail ----
    ws = wb.create_sheet('delta_detail')
    write_delta_detail_block(ws, all_pairs_data)

    # ---- ellipse_统计 ----
    ws = wb.create_sheet('ellipse_统计')
    r = 1

    # Block A-1: mn 三分类（init 合并 4λ×runs, 行=MC/NC-MC/IC≥3）
    ws.cell(row=r, column=1, value='Block A-1: mn 三分类 (MC/NC-MC/IC≥3, init 合并4λ×runs)').font = TITLE_FONT
    r += 1
    for j, h in enumerate(['分类', '均值', '标准差']): style_header(ws.cell(row=r, column=j + 1, value=h))
    r += 1
    for cat_idx, cat_name in enumerate(['MC', 'NC-MC', 'IC≥3']):
        m, s = mean_std([row.get('相邻细胞边数和') for row in tway_pool[cat_idx]])
        _write_cell(ws, r, 1, cat_name)
        _write_cell(ws, r, 2, m, '0.0000')
        _write_cell(ws, r, 3, s, '0.0000')
        r += 1
    r += 1

    # Block A-2: n/P/A（init 合并 4λ×runs, final × 4λ, MC/IC）
    for col_key, short in [('细胞边数', 'n'), ('细胞周长', 'P'), ('细胞面积', 'A')]:
        ws.cell(row=r, column=1, value=f'Block A-2-{short}: {short} (init 合并4λ×runs, final×4λ, MC/IC)').font = TITLE_FONT
        r += 1
        headers_a2 = ['分类', 'init_均值', 'init_标准差']
        for lam in LAMS:
            headers_a2 += [f'{lam}_final_均值', f'{lam}_final_标准差']
        for j, h in enumerate(headers_a2): style_header(ws.cell(row=r, column=j + 1, value=h))
        r += 1
        for sub_label, mask_fn in [('MC', lambda rr: int(rr.get('当前细胞层', 0)) == 1),
                                    ('IC', lambda rr: int(rr.get('当前细胞层', 0)) >= 2)]:
            _write_cell(ws, r, 1, sub_label)
            col = 2
            m, s = mean_std([row.get(col_key) for row in ell_init_pool if mask_fn(row)])
            _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
            col += 2
            for lam in LAMS:
                rows = [row for row in ell[lam][1] if mask_fn(row)]
                m, s = mean_std([row.get(col_key) for row in rows])
                _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
                col += 2
            r += 1
        r += 1

    # Block B-0: n 边形分布（init 合并 4λ×runs, MC/IC）
    ws.cell(row=r, column=1, value='Block B-0: n 边形分布 (init 合并4λ×runs, MC/IC)').font = TITLE_FONT
    r += 1
    for sub_label, mask_fn in [('MC', lambda rr: int(rr.get('当前细胞层', 0)) == 1),
                                ('IC', lambda rr: int(rr.get('当前细胞层', 0)) >= 2)]:
        rows = [row for row in ell_init_pool if mask_fn(row)]
        ws.cell(row=r, column=1, value=sub_label).font = Font(bold=True)
        r += 1
        headers = ['n值', '数量', '频率']
        for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j+1, value=h))
        r += 1
        by_n = group_by_n(rows)
        total = sum(len(v) for v in by_n.values())
        for nv in sorted(by_n.keys()):
            cnt = len(by_n[nv])
            _write_cell(ws, r, 1, nv)
            _write_cell(ws, r, 2, cnt)
            _write_cell(ws, r, 3, f'{cnt/total:.1%}' if total else '0')
            r += 1
        r += 1

    # Block B: P/A 按 n 分组（4λ × 2态 × MC/IC）
    for col_key, short in [('细胞周长', 'P'), ('细胞面积', 'A')]:
        ws.cell(row=r, column=1, value=f'Block B-{short}: {short} 按 n 分组').font = TITLE_FONT
        r += 1
        for lam in LAMS:
            for it in ITERATIONS:
                for sub_label, rows in [('MC', mc_ell[lam][it]), ('IC', ic_ell[lam][it])]:
                    ws.cell(row=r, column=1, value=f'{lam}_{PHASE_MAP[it]}_{sub_label}').font = Font(bold=True)
                    r += 1
                    headers = ['n值', f'{short}_均值', f'{short}_标准差']
                    for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j+1, value=h))
                    r += 1
                    by_n = group_by_n(rows)
                    for nv in sorted(by_n.keys()):
                        vals = [float(row[col_key]) for row in by_n[nv] if row.get(col_key) is not None]
                        m, s = mean_std(vals)
                        _write_cell(ws, r, 1, nv)
                        _write_cell(ws, r, 2, m, '0.0000')
                        _write_cell(ws, r, 3, s, '0.0000')
                        r += 1
                    r += 1
            r += 1

    # Block C: Delta 按 n
    r = _write_delta_block_c(ws, r, 'Block C: Delta 按 n (4λ × MC/IC)', all_pairs_data)

    # Block D: C1/C2
    ws.cell(row=r, column=1, value='Block D: C1/C2 拟合 (final-initial)').font = TITLE_FONT
    r += 1
    headers = ['λ', '参数', 'MC', 'IC']
    for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j+1, value=h))
    r += 1
    for lam in LAMS:
        results = c1_c2_results[lam]
        for param_name, key in [('C1', 'C1'), ('C2', 'C2'), ('RMSE', 'RMSE'), ('细胞数', 'n_cells')]:
            _write_cell(ws, r, 1, lam)
            _write_cell(ws, r, 2, param_name)
            for ci, label in enumerate(['MC', 'IC']):
                val = results.get(label, {}).get(key, 0)
                fmt = '0.000000' if isinstance(val, float) else None
                _write_cell(ws, r, 3 + ci, int(val) if key == 'n_cells' else val, fmt)
            r += 1
    r += 1

    all_c1[net] = c1_c2_results

    # ---- edgeAngle ----
    ws = wb.create_sheet('edgeAngle')
    titles_ea, data_ea = [], []
    for lam in LAMS:
        for it in ITERATIONS:
            titles_ea += [f'{lam}_{PHASE_MAP[it]}_MC', f'{lam}_{PHASE_MAP[it]}_IC']
            data_ea += [mc_ea[lam][it], ic_ea[lam][it]]
    write_raw_block(ws, titles_ea, data_ea,
                    ['run', '细胞序号', '边数', '相邻细胞边数和', '内角', '夹边1', '当前细胞层'])

    ws = wb.create_sheet('edgeAngle_统计')
    r = 1
    ea_init_pool = []
    for lam in LAMS:
        ea_init_pool += ea[lam][0]
    for col_key, label in [('内角', '内角[度]'), ('夹边1', '夹边1')]:
        ws.cell(row=r, column=1, value=f'Block A: 整体 ({label}, init 合并4λ×runs)').font = TITLE_FONT
        r += 1
        headers = ['分类', 'init_均值', 'init_标准差']
        for lam in LAMS:
            headers += [f'{lam}_final_均值', f'{lam}_final_标准差']
        for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j + 1, value=h))
        r += 1
        for sub_label, mask_fn in [('MC', lambda rr: int(rr.get('当前细胞层', 0)) == 1),
                                    ('IC', lambda rr: int(rr.get('当前细胞层', 0)) >= 2)]:
            _write_cell(ws, r, 1, sub_label)
            col = 2
            m, s = mean_std([row.get(col_key) for row in ea_init_pool if mask_fn(row)])
            _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
            col += 2
            for lam in LAMS:
                rows = [row for row in ea[lam][1] if mask_fn(row)]
                m, s = mean_std([row.get(col_key) for row in rows])
                _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
                col += 2
            r += 1
        r += 1
    # 分布块
    groups_ea = {'init_MC': [row for row in ea_init_pool if int(row.get('当前细胞层', 0)) == 1],
                 'init_IC': [row for row in ea_init_pool if int(row.get('当前细胞层', 0)) >= 2]}
    for lam in LAMS:
        for sub_label, rows in [('MC', mc_ea[lam][1]), ('IC', ic_ea[lam][1])]:
            groups_ea[f'{lam}_final_{sub_label}'] = rows
    r = _write_dist_alt(ws, r, 'Block B: 边长分布 (夹边1)', groups_ea, '夹边1', length_bin, LENGTH_LABELS)
    r = _write_dist_alt(ws, r, 'Block C: 角度分布 (内角, 已为度)', groups_ea, '内角', angle_bin, ANGLE_LABELS)

    # ---- ME_MA ----
    ws = wb.create_sheet('ME_MA')
    mm_mc = {lam: {it: [r for r in mm[lam][it] if int(r.get('当前细胞层', 0)) == 1]
                     for it in ITERATIONS} for lam in LAMS}
    titles_mm, data_mm = [], []
    for lam in LAMS:
        for it in ITERATIONS:
            titles_mm += [f'{lam}_{PHASE_MAP[it]}_MC']
            data_mm += [mm_mc[lam][it]]
    write_raw_block(ws, titles_mm, data_mm,
                    ['run', '细胞序号', '边缘边长(ME)', '边缘角1(MA1)', '边缘角2(MA2)', '当前细胞层'])

    ws = wb.create_sheet('ME_MA_统计')
    r = 1
    mm_init_pool = []
    for lam in LAMS:
        mm_init_pool += mm_mc[lam][0]
    for col_key, label in [('边缘边长(ME)', 'ME'), ('边缘角1(MA1)', 'MA1'), ('边缘角2(MA2)', 'MA2')]:
        ws.cell(row=r, column=1, value=f'Block A: 整体 ({label}, init 合并4λ×runs)').font = TITLE_FONT
        r += 1
        headers = ['init_均值', 'init_标准差']
        for lam in LAMS:
            headers += [f'{lam}_final_均值', f'{lam}_final_标准差']
        for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j + 1, value=h))
        r += 1
        col = 1
        m, s = mean_std([row.get(col_key) for row in mm_init_pool])
        _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
        col += 2
        for lam in LAMS:
            rows = mm_mc[lam][1]
            m, s = mean_std([row.get(col_key) for row in rows])
            _write_cell(ws, r, col, m, '0.0000'); _write_cell(ws, r, col + 1, s, '0.0000')
            col += 2
        r += 1
    groups_mm = {'init_MC': mm_init_pool}
    for lam in LAMS:
        groups_mm[f'{lam}_final_MC'] = mm_mc[lam][1]
    r = _write_dist_alt(ws, r, 'Block B: ME 分布', groups_mm, '边缘边长(ME)', length_bin, LENGTH_LABELS)
    r = _write_dist_alt(ws, r, 'Block C: MA1 分布 (已为度)', groups_mm, '边缘角1(MA1)', angle_bin, ANGLE_LABELS)
    r = _write_dist_alt(ws, r, 'Block D: MA2 分布 (已为度)', groups_mm, '边缘角2(MA2)', angle_bin, ANGLE_LABELS)

    for ws in wb.worksheets:
        for cc in range(1, ws.max_column + 1):
            ws.column_dimensions[get_column_letter(cc)].width = 14

    wb.save(out_path)
    print(f'  Saved: {out_path}')

# ---------- 跨组对比 ----------
def generate_c1_summary(base, all_c1):
    out_path = os.path.join(base, 'iteration_C1_C2_汇总对比.xlsx')
    wb = Workbook(); wb.remove(wb.active)
    ws = wb.create_sheet('C1_C2汇总')
    headers = ['初始网络', 'λ', 'C1_MC', 'C1_IC', 'C2_MC', 'C2_IC',
               'RMSE_MC', 'RMSE_IC', '细胞数_MC', '细胞数_IC']
    for cc in range(1, len(headers) + 1):
        ws.column_dimensions[get_column_letter(cc)].width = 14
    r = 1
    for net in sorted(all_c1.keys()):
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=len(headers))
        c = ws.cell(row=r, column=1, value=f'===== {net} =====')
        c.font = TITLE_FONT; c.fill = GROUP_FILL; c.alignment = CENTER
        r += 1
        for j, h in enumerate(headers): style_header(ws.cell(row=r, column=j + 1, value=h))
        r += 1
        for lam in LAMS:
            results = all_c1[net][lam]
            _write_cell(ws, r, 1, net)
            _write_cell(ws, r, 2, lam)
            for ci, key in enumerate(['C1', 'C2', 'RMSE']):
                for li, label in enumerate(['MC', 'IC']):
                    val = results.get(label, {}).get(key, 0)
                    _write_cell(ws, r, 3 + ci * 2 + li, val, '0.000000')
            for li, label in enumerate(['MC', 'IC']):
                val = results.get(label, {}).get('n_cells', 0)
                _write_cell(ws, r, 9 + li, int(val))
            r += 1
        r += 1
    wb.save(out_path)
    print(f'Saved: {out_path}')

if __name__ == '__main__':
    base_dir = r'data_fixed'
    nets = discover_nets(base_dir)
    print(f'发现 {len(nets)} 个初始网络: {nets}')
    print('-' * 50)
    all_c1 = {}
    for net in sorted(nets):
        print(f'Processing: {net} ...')
        process_net(base_dir, net, all_c1)
        print(f'  Done: {net}')
    print('-' * 50)
    if all_c1: generate_c1_summary(base_dir, all_c1)
    print('-' * 50)
    print('ALL DONE!')
