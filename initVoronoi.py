# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial import Voronoi
import random
import math
from utillib.mylib import Cell


def rand_in_unit_hex(scale=1.0):
    """
    菱形采样法：在单位正六边形内生成均匀随机点
    将正六边形分成3个菱形，随机选择一个菱形，用两个[0,1]的随机数生成点

    :param scale: 缩放因子，用于控制扰动强度（对应原来的param）
    :return: (dx, dy) 正六边形内的随机扰动偏移量
    """
    # 正六边形核心参数：边长 s = sqrt(3)/3
    s = math.sqrt(3) / 3

    # 单位正六边形（尖头垂直朝上）的3组相邻顶点向量（从原点指向顶点）
    vectors = [
        (0.0, 1.0),                    # 向量1：正上方
        (math.sqrt(3)/2, -1/2),        # 向量2：右下方
        (-math.sqrt(3)/2, -1/2)         # 向量3：左下方
    ]

    # 1. 随机选择1组相邻向量（3组中选1组）
    vec_idx = random.randrange(3)
    v1 = vectors[vec_idx]
    v2 = vectors[(vec_idx + 1) % 3]  # 第二个相邻向量（循环取余，保证相邻）

    # 2. 生成2个[0,1]区间的均匀随机数（菱形采样核心）
    u = random.random()
    v = random.random()

    # 3. 向量线性组合得到单位正六边形内的点
    unit_hex_dx = u * v1[0] + v * v2[0]
    unit_hex_dy = u * v1[1] + v * v2[1]

    # 4. 缩放至目标正六边形（边长s），并应用扰动强度scale
    dx = unit_hex_dx * s * scale
    dy = unit_hex_dy * s * scale

    return (dx, dy)


def new_vor(param, param2, n):  # 构造随机维诺图数据 Construct random venogram data
    """
    使用菱形采样法生成随机扰动的Voronoi图种子点
    替换原来的拒绝法，提高均匀度
    """
    arr_list = list()
    for x in range(n):
        for y in range(n+1):
            tmp_param = param
            # 边缘区域不进行扰动
            if x in range(0,param2) or x in range(n-param2,n) or y in range(0,param2+1) or y in range(n-param2, n+2) :
                tmp_param = 0

            # 计算正六边形网格的基础坐标
            point_x = x
            point_y = y
            if y % 2 == 0:
                point_x += 0.5

            point_y = 0.5 * point_y * np.sqrt(3)

            # 使用菱形采样法生成扰动（替换原来的拒绝法）
            if tmp_param > 0:
                tmp_x, tmp_y = rand_in_unit_hex(tmp_param)
            else:
                tmp_x, tmp_y = 0.0, 0.0

            point_x = point_x + tmp_x
            point_y = point_y + tmp_y
            p = (point_x, point_y)
            arr_list.append(p)
    return arr_list


def new_vor_random(n, scale=None):
    """
    生成均匀随机的Voronoi图种子点（类似R代码的deldir方式）

    :param n: 种子点数量
    :param scale: 坐标范围缩放因子。如果为None，则使用sqrt(n)；如果指定，则使用该值作为范围上限
    :return: 种子点列表
    """
    if scale is None:
        # 默认使用sqrt(n)作为坐标范围上限
        scale = math.sqrt(n)

    arr_list = []
    for i in range(n):
        # 生成[0, scale)范围内的随机点
        x = random.uniform(0, scale)
        y = random.uniform(0, scale)
        arr_list.append((x, y))

    return arr_list


def just(vor, n=None, x_max=None, y_max=None):
    """
    处理Voronoi图，提取有效的细胞区域

    :param vor: Voronoi对象
    :param n: 网格大小（用于规则网格的边界判断，如果为None则使用x_max和y_max）
    :param x_max: x坐标最大值（用于随机点生成的边界判断）
    :param y_max: y坐标最大值（用于随机点生成的边界判断）
    """
    regions = vor.regions  # 所有细胞的点位置 The point position of all the cells
    vertices = vor.vertices  # 细胞点坐标 Cell point coordinates
    cells = list()

    # 确定边界范围
    if n is not None:
        # 规则网格模式：使用n+1作为边界
        x_bound = n + 1
        y_bound = n + 1
    elif x_max is not None and y_max is not None:
        # 随机点模式：使用指定的边界
        x_bound = x_max
        y_bound = y_max
    else:
        # 如果没有指定边界，尝试从vertices中推断
        if len(vertices) > 0:
            x_bound = max(v[0] for v in vertices if len(v) > 0) + 1
            y_bound = max(v[1] for v in vertices if len(v) > 0) + 1
        else:
            x_bound = y_bound = 100  # 默认值

    for region in regions:
        flag = 0
        points = []
        for po in region:
            if po == -1:
                flag = 1
                break
            elif len(vertices[po]) == 0:
                break
            elif vertices[po][0] < -1 or vertices[po][0] > x_bound or vertices[po][1] < -1 or vertices[po][1] > y_bound:
                flag = 1
                break
            point = vertices[po]
            points.append(point)

        if len(points) < 3:
            continue
        if flag == 0 and points != []:
            c1 = Cell(points)
            if c1.area < 0:
                Cell.cell_no -= 1
                c1 = Cell(points[::-1])
            cells.append(c1)

    return cells


'''
    细胞集预处理，含有：
    1、去除边缘独立细胞
    2、细胞块边缘整齐化
'''
def pretreatment(cells):
    # num = 0
    points_times = {}
    for c in cells:
        for p in c.points:
            p_key = '{:.10f}-{:.10f}'.format(p[0], p[1])
            if p_key in points_times.keys():
                points_times[p_key] += 1
            else:
                points_times[p_key] = 1

    # 去除边缘独立细胞
    i = -1
    while True:
        i += 1
        if i >= len(cells):
            break
        flag_count = 0
        for p in cells[i].points:
            p_key = '{:.10f}-{:.10f}'.format(p[0], p[1])
            if points_times[p_key] > 1:
                flag_count += 1
        # print(flag_count)
        if flag_count < 3:
            # print(flag_count)
            # num+=1
            c = cells.pop(i)
            i -= 1
            # print(11111)
            for p in c.points:
                p_key = '{:.10f}-{:.10f}'.format(p[0], p[1])
                points_times[p_key] -= 1

    # 细胞块边缘整齐化
    for c in cells:
        p_i = -1
        while True:
            p_i += 1
            if p_i >= len(c.points):
                break
            p = [c.points[p_i][0], c.points[p_i][1]]
            p_key = '{:.10f}-{:.10f}'.format(p[0], p[1])
            if points_times[p_key] == 1:
                c.points.pop(p_i)
                p_i -= 1
        c.setArea()
        c.setVertex()
    # print("num", num)
    return cells


def getCells(param, param2, n):
    param_float = float(param)
    n_int = int(n) + 1
    param2_int = int(param2)
    arr = new_vor(param_float, param2_int, n_int)
    vor = Voronoi(arr)
    cells = just(vor, n=n_int)
    cells = pretreatment(cells)
    return cells


def getCellsRandom(n, scale=None):
    """
    生成均匀随机的Voronoi图（类似R代码的deldir方式）

    :param n: 种子点数量
    :param scale: 坐标范围缩放因子。如果为None，则使用sqrt(n)；如果指定，则使用该值作为范围上限
    :return: Cell对象列表
    """
    arr = new_vor_random(n, scale)
    # 根据随机点的范围计算边界
    xs = [p[0] for p in arr]
    ys = [p[1] for p in arr]
    x_max = max(xs) if xs else 1.0
    y_max = max(ys) if ys else 1.0

    vor = Voronoi(arr)
    cells = just(vor, x_max=x_max, y_max=y_max)
    cells = pretreatment(cells)
    return cells


def getCellsBySource(source_type, param=0, param2=0, n=10, scale=None):
    """
    统一的细胞生成入口，根据source_type选择生成方式

    :param source_type: 'grid'（随机扰动维诺图n×n）, 'random'（均匀随机维诺图）
    :param param: 扰动参数（grid模式）或种子点数（random模式）
    :param param2: 边缘层数（grid模式）
    :param n: 种子数（random模式）
    :param scale: 坐标范围（random模式）
    :return: Cell对象列表
    """
    if source_type == 'random':
        return getCellsRandom(n, scale)
    else:
        # 默认：随机扰动维诺图（grid模式）
        return getCells(param, param2, n)
