# -*- coding: utf-8 -*-
"""
_convert_for_iteration.py — 把 data/ 下的批量输出转换为 iteration_summarize 可读结构。

转换映射:
  data/{net}_L{lam}_rep{r}/{net}_L{lam}_rep{r}_{init|final}_{type}.xlsx
    → data_fixed/{net}/L{lam}/run{r+1}/{0|1}_{initial|final}_{type}.xlsx

示例:
  data/grid02_L0_rep0/grid02_L0_rep0_init_ellipse.xlsx
    → data_fixed/grid02/L0/run1/0_initial_ellipse.xlsx

列名: 英文 → 中文（仅改表头，不改数值；内角保持度数不转换）
"""
import os
import re
import shutil
import openpyxl

SRC = r'data'
DST = r'data_fixed'

# 列名映射（英文 -> 中文）
COL_MAP = {
    # ellipse
    'Cell ID': '细胞序号',
    'Ellipse Centroid': '椭圆中心',
    'Major Semi-axis': '长半轴',
    'Minor Semi-axis': '短半轴',
    'Inclination Angle of the Major Semi-axis': '倾角',
    'Sum of Neighbor Edge Numbers': '相邻细胞边数和',
    'Cell Edge Number': '细胞边数',
    'Cell Perimeter': '细胞周长',
    'Cell Area': '细胞面积',
    'Original Cell Layer': '原始细胞层',
    'Current Cell Layer': '当前细胞层',
    # edgeAngle
    'Edge Number': '边数',
    'Interior Angle': '内角',
    'Adjacent Edge 1': '夹边1',
    'Adjacent Edge 2': '夹边2',
    # ME_MA
    'Marginal Edge Length (ME)': '边缘边长(ME)',
    'Marginal Angle 1 (MA1)': '边缘角1(MA1)',
    'Marginal Angle 2 (MA2)': '边缘角2(MA2)',
}

def convert_xlsx(src_path, dst_path):
    """读 xlsx，替换表头，写回新路径。"""
    wb = openpyxl.load_workbook(src_path)
    for ws in wb.worksheets:
        for cell in ws[1]:
            if cell.value in COL_MAP:
                cell.value = COL_MAP[cell.value]
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
    wb.save(dst_path)

def main():
    if os.path.exists(DST):
        shutil.rmtree(DST)
    os.makedirs(DST)

    rep_dirs = sorted([d for d in os.listdir(SRC) if os.path.isdir(os.path.join(SRC, d))])
    print(f'发现 {len(rep_dirs)} 个 rep 文件夹')

    n_files = 0
    for rep_dir in rep_dirs:
        m = re.match(r'^(.+)_L(\d+)_rep(\d+)$', rep_dir)
        if not m:
            print(f'  跳过: {rep_dir} (命名不匹配)')
            continue
        net, lam, rep = m.group(1), m.group(2), int(m.group(3))
        run = rep + 1

        src_sub = os.path.join(SRC, rep_dir)
        dst_sub = os.path.join(DST, net, f'L{lam}', f'run{run}')
        os.makedirs(dst_sub, exist_ok=True)

        for fname in os.listdir(src_sub):
            if not fname.endswith('.xlsx'):
                continue
            m2 = re.match(rf'^{re.escape(rep_dir)}_(init|final)_(ellipse|edgeAngle|ME_MA)\.xlsx$', fname)
            if not m2:
                continue
            phase_src, typ = m2.group(1), m2.group(2)
            iter_num = 0 if phase_src == 'init' else 1
            phase_dst = 'initial' if phase_src == 'init' else 'final'
            new_name = f'{iter_num}_{phase_dst}_{typ}.xlsx'
            convert_xlsx(os.path.join(src_sub, fname),
                        os.path.join(dst_sub, new_name))
            n_files += 1

    print(f'完成: {n_files} 个 xlsx 转换 → {DST}')
    # 输出目录结构摘要
    nets = sorted([d for d in os.listdir(DST) if os.path.isdir(os.path.join(DST, d))])
    for net in nets:
        lams = sorted([d for d in os.listdir(os.path.join(DST, net))
                      if os.path.isdir(os.path.join(DST, net, d))])
        print(f'  {net}/: {lams}')

if __name__ == '__main__':
    main()
