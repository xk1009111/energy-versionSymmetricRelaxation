# 细胞退火过程详细流程

> 本文档基于当前源码整理，描述从 GUI 触发到顶点移动的完整退火链路。
> 所有行号均与源文件当前状态一致，可点击跳转。

---

## 1. 概述

退火的核心目标：**把每个细胞顶点朝"理想位置"逐步移动，使细胞形状趋向规则（接近圆形 / 等角分布），同时保证拓扑与几何合法性（不破坏凸性、不恶化内角）。**

退火分两类顶点：

| 类型 | 定义 | 移动目标 |
|------|------|----------|
| **内部顶点 (inner)** | 被 **3 个**细胞共享的顶点 | 由 3 条预最优线交出的"拟合三角形"的重心 |
| **边缘顶点 (marginal)** | 只被 **2 个**细胞共享的顶点 | 使两个边缘角 ∠AVO、∠BVO 趋于相等的目标点（V 与较小角邻点的中点） |

两类顶点在 [`move_point`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1143) 中被统一收集到一个队列，按"退火距离"降序处理后逐个移动。

> **关于面积**：退火决策**不直接使用面积**——面积只是顶点移动的被动副产品。详见 [第 11 节](#11-面积被动结果与影响因素)。

---

## 2. 整体流程图

```
┌──────────────────────────────────────────────────────────────┐
│  GUI: _run_anneal  (energy_annealing_main.py:841)              │
│  读取参数 → 构造 Annealer → 循环 times 轮                     │
└───────────────────────────┬──────────────────────────────────┘
                            │ 每轮调用
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  Annealer.annealing(cells)  (AnnealingGUI.py:53)             │
│                                                              │
│  ① 为每个细胞构建连线：                                        │
│     actual_lines  → get_actual_lines           (Util:129)    │
│     delta         → get_best_rotate_delta_by_calculation (Util:189) │
│     pre_best_lines→ get_pre_best_lines         (Util:161)    │
│                                                              │
│  ② 提取交汇细胞块：                                           │
│     intersection_cell_blocks → get_intersection_cell_blocks (Util:242) │
│                                                              │
│  ③ 统一移动：                                                 │
│     move_point(blocks, rate, marginal_judge, cells) (Util:1143) │
│        返回 (annealing_count, {inner_points, marginal_points}) │
└───────────────────────────┬──────────────────────────────────┘
                            │
                            ▼
┌──────────────────────────────────────────────────────────────┐
│  move_point 内部（统一队列）                                   │
│                                                              │
│  Step1 收集顶点                                               │
│    ├─ marginal: get_all_marginal_points (Util:669)           │
│    │            calculate_marginal_annealing_distance (Util:946) → D │
│    └─ inner:    来自 intersection_cell_blocks                  │
│                 D = 当前点到三角形重心的距离                    │
│                                                              │
│  Step2 按 D 降序排序（一次排序，本轮不变）                      │
│                                                              │
│  Step3 依序遍历                                               │
│    ├─ marginal → get_marginal_move_point (Util:996)          │
│    └─ inner    → get_point_of_destination (Util:321)         │
│                 judge_if_annealing (Util:649)                 │
│                 通过 → 写回 3 个细胞 + setVertex               │
└──────────────────────────────────────────────────────────────┘
```

---

## 3. 驱动层：多轮循环

入口：[`_run_anneal`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/energy_annealing_main.py#L841)

1. **校验**：必须先生成细胞（`cellData` 非空）；校验退火速率、边缘开关、内角约束开关。
2. **构造退火器**：
   ```python
   anneal_params = {
       'annealingRate': float(rate_str),          # 单步退火比例 0~1
       'marginal_point_judge': int(edge_str),     # 0=不退火边缘, 1=退火边缘
       'inner_angle_sq_guard': int(self.anneal_inner_angle_guard.get())  # 0/1
   }
   annealer = Annealer(anneal_params)
   ```
3. **循环 `times` 轮**（`times` 来自"退火次数"输入框，默认 1），每轮：
   - `annealer.annealing(cellData.cells)` 执行一轮
   - `cellData.flush()` 刷新细胞数据（重算中心点 / 面积 / 拓扑）
   - 更新右侧统计面板（内部顶点 / 边缘顶点 / 总细胞数 / 当前轮次）
   - 重绘画布；多轮时每轮间隔 0.5s 可视化

---

## 4. 单轮退火：`Annealer.annealing`

位置：[`AnnealingGUI.py:53`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/AnnealingGUI.py#L53)

### 4.1 为每个细胞构建连线

对 `cells` 中每个细胞：

- `max_distance = get_distance_centerpoint_point(cell)` —— 细胞中心点到顶点的**最大距离**（作为射线的极径 r）。见 [`Util:106`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L106)
- `cell.actual_lines = get_actual_lines(cell, max_distance)` —— 从中心点指向每个顶点的实际射线（极坐标 θ, r）。见 [`Util:129`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L129)
- `delta = get_best_rotate_delta_by_calculation(cell)` —— **计算法**求使"实际角与理想均分角偏差平方和"最小的旋转增量。见 [`Util:189`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L189)
- `cell.pre_best_lines = get_pre_best_lines(cell, max_distance, delta)` —— 旋转 `delta` 后、把 2π/N 均匀分布的"预最优射线"对齐到第一个顶点。见 [`Util:161`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L161)

> 直觉：`pre_best_lines` 是"细胞形状完全规则时各顶点应在的方向"；`actual_lines` 是当前实际方向。两者越接近，细胞越规则。

### 4.2 提取交汇细胞块

`intersection_cell_blocks = get_intersection_cell_blocks(cells)` —— 见 [`Util:242`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L242)

逻辑：
1. 遍历所有细胞的所有顶点，统计每个点被几个细胞共享。
2. **被恰好 3 个细胞共享的点** → 构造一个 `CellBlock`，记录这 3 个细胞（`cell1/2/3`）及该点在各细胞中的索引（`index1/2/3`）。
3. 用每个细胞的 `pre_best_lines`（过中心点的射线）两两求交，得到 3 个交点 → 拼成 **拟合三角形** `cell_block.triangle`（3 个顶点）。

这个三角形就是内部顶点的"理想落点区域"，其重心即退火目标。

### 4.3 统一移动

```python
result = util.move_point(intersection_cell_blocks, self.annealingRate,
                         self.marginal_point_judge, cells)
```

返回 `(annealing_count, {'inner_points': int, 'marginal_points': int})`，`Annealer` 把统计写入 `self.inner_points / self.marginal_points` 供面板显示。

---

## 5. `move_point` 统一队列详解

位置：[`Util:1143`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1143)

### 5.1 Step 1 — 收集顶点到队列 `vertex_queue`

**边缘顶点**（仅当 `marginal_point_judge=True`）：
- `all_marginal_points = get_all_marginal_points(cells)`（[`Util:669`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L669)）—— 筛选只被 2 个细胞共享的顶点。
- 对每个 `point_v`：`D = calculate_marginal_annealing_distance(point_v, cells)`（[`Util:946`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L946)，[行 1175](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1175)）—— **不实际移动**，只算退火距离用于排序。
- 入队：`{'type':'marginal', 'point':point_v, 'distance':D}`

**内部顶点**：
- 遍历每个 `cb in intersection_cell_blocks`：
  - `point_g = cb.getTriCentreOfGravity()` —— 拟合三角形重心（`Point` 对象，[`mylib:218`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/utillib/mylib.py#L218)）。
  - `current_point = cb.cell1.points[cb.index1]` —— 该顶点当前坐标 `[x,y]`。
  - `D = ‖current_point - point_g‖`（[行 1187](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1187)）—— 当前点到重心的欧氏距离。
- 入队：`{'type':'inner', 'cb':cb, 'distance':D}`

> 注意：`calculate_marginal_annealing_distance` **不做角度阈值判断**，恒返回 V 到目标中点的正距离（用于排序）；阈值判断只发生在移动阶段（`get_marginal_move_point` 内）。

### 5.2 Step 2 — 按距离 D 降序排序

```python
vertex_queue.sort(key=lambda x: x['distance'], reverse=True)   # 行 1197
```

**一次排序，本轮不再重排**：距离大的（偏离最优最严重）先处理，保证最需要退火的顶点优先得到机会。

### 5.3 Step 3 — 依序遍历执行（[行 1202-1235](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1202)）

- **marginal 项**：**不重算排序距离 D**，直接调用 `get_marginal_move_point(point_v, annealing_rate, cells)`。该函数内部会用最新坐标重新找关键点、重算两角与阈值、做凸性检查；返回 `>0` 表示成功移动（`marginal_annealing_points += 1`），返回 `0` 表示跳过（角度差低于阈值或凸性失败）。
- **inner 项**：见下一节。

> **排序距离只算一次**：marginal 顶点的 D 仅用于排序，循环内不再重算；实际是否移动由 `get_marginal_move_point` 内部的角度阈值与凸性检查决定——返回 `>0` 表示移动成功，返回 `0` 表示跳过。

---

## 6. 内部顶点退火 (inner)

### 6.1 目标点：三角形重心

每次都用 `cb.getTriCentreOfGravity()` **重新**取重心（坐标可能已被前序移动更新）。重心 = 拟合三角形三个顶点坐标的平均。

### 6.2 计算移动目标

`move_point_result = get_point_of_destination(current_point, point_g, annealing_rate)` —— 见 [`Util:321`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L321)

公式：
```
move_point = current_point + annealing_rate * (point_g - current_point)
```
即从当前位置朝重心走 `annealing_rate` 比例（`annealing_rate=0.5` 走一半）；位移大小 = `annealing_rate × D`。

### 6.3 可行性判定：`judge_if_annealing`

`flag_index = judge_if_annealing(cb, move_point_result)` —— 见 [`Util:649`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L649)

按顺序检查，返回第一个不满足的原因码：

| flag_index | 含义 | 判定函数 | 在 move_point 中的处理 |
|:---:|------|------|------|
| **-1** | 当前顶点已在拟合三角形**内部**，已接近最优 | [`is_point_in_triangle`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L564)（向量法 u,v>0 且 u+v<1） | 不移动，`continue`（不单独计数） |
| **-3** | 移动后会**破坏凸性**（某细胞变非凸） | [`judge_by_intersection_cell_blocks`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L413) → [`judge_by_cell`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L434) → [`judge_by_change`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L466) | `judge_180_count++`，不移动 |
| **-4** | 移动后该顶点对应的 **3 个内角平方和会增大** | [`judge_sum_inner_angle2`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L523)（仅当 `USE_INNER_ANGLE_SQ_GUARD=True`） | 不移动，`continue`（不单独计数） |
| **0** | 全部通过，可以移动 | — | `move_flag=True`，`inner_annealing_points++` |

> 当前 `move_point` 只显式统计 `judge_180_count`（flag -3）；flag -1 与 -4 直接 `continue` 跳过、不单独计数。
> **内角平方和约束**是可选开关，由 GUI "退火设定"中的 `inner_angle_sq_guard` 控制，通过 [`set_annealing_options`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L643) 写入全局变量 `USE_INNER_ANGLE_SQ_GUARD`。默认开启。

### 6.4 应用移动

通过判定后，把同一个 `move_point_result` 写回三个细胞对应索引位置，并重算中心点（[`setVertex`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/utillib/mylib.py#L99)）：

```python
cb.cell1.points[cb.index1] = [move_point_result.x, move_point_result.y]
cb.cell2.points[cb.index2] = [move_point_result.x, move_point_result.y]
cb.cell3.points[cb.index3] = [move_point_result.x, move_point_result.y]
cb.cell1.setVertex(); cb.cell2.setVertex(); cb.cell3.setVertex()
```

> 因为三个细胞共享同一个顶点，必须**同步更新**，否则会产生拓扑裂缝。
> 注意 `setVertex()` 只重算**中心点**，**不**重算面积（面积在每轮 `flush` 时统一重算，见第 11 节）。

---

## 7. 边缘顶点退火 (marginal)

### 7.1 找出所有边缘顶点

[`get_all_marginal_points`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L669)：统计每个顶点被几个细胞共享，**恰好 2 个**的就是边缘顶点。

### 7.2 找关键点：V → A, B → 边缘细胞 → O

[`find_marginal_key_points_new`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L710) 对给定边缘顶点 V 还原其局部拓扑：

```
        A ─────── V ─────── B
         \ cell1  |  cell2 /
          \       |       /
           \      O      /
            ───────────
        (O 是 cell1、cell2 共有的内部邻点)
```

步骤：
1. 确认 V 是边缘顶点（被 2 个细胞共享）。
2. 在 V 的邻点中找**也是边缘顶点**的两个点 A、B（分别属于两个不同的 `layer==1` 边缘细胞 `marginal_cell1`、`marginal_cell2`）。
3. 找 O：在 `marginal_cell1` 中 V 的邻点里，找一个**同时也在 `marginal_cell2` 中且与 V 相邻**的点（即两细胞的共边内端点）。
4. 返回 `{point_v, point_a, point_b, point_o, marginal_cell1, marginal_cell2, idx_va, idx_vb, idx_oa, idx_ob}`，找不到返回 `None`。

### 7.3 计算退火距离（仅用于排序）

[`calculate_marginal_annealing_distance`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L946)：

1. 用 `find_marginal_key_points_new` 取关键点。
2. 用向量夹角（`acos`）算两个边缘角：`angle_AVO = ∠(A, V, O)`、`angle_BVO = ∠(B, V, O)`。
3. 确定目标点（向**较小角**方向移动，缩短对应边）：`angle_AVO < angle_BVO` → `target = midpoint(V, A)`，否则 → `target = midpoint(V, B)`。
4. 返回 `D = ‖V - target‖`（**不乘退火速率，不做阈值判断**，统一用于排序）。

> 阈值判断不在此函数，而在下面的 `get_marginal_move_point` 内。

### 7.4 执行边缘移动

[`get_marginal_move_point`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L996)：

1. 重新调用 `find_marginal_key_points_new` 取关键点（用最新坐标）。
2. 重新算 `angle_AVO`、`angle_BVO`。
3. **形状审查设阈值**（[行 1062-1067](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1062)）：
   - `marginal_cell1` 或 `marginal_cell2` 是三角形（3 边）→ 阈值 **60°**
   - 否则 → 阈值 **20°**
4. 若 `|angle_AVO - angle_BVO| < 阈值` → `return 0`（跳过，[行 1073-1076](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1073)）。
5. **生成候选点**（[行 1092-1093](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1092)）：
   ```python
   target_point = midpoint(V, 较小角对应的 A 或 B)
   candidate_V = V + annealing_rate * (target_point - V)
   ```
   位移大小 = `annealing_rate × D`。
6. **凸性检查**（保证移动后两个边缘细胞仍为凸多边形，[行 1098-1105](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1098)）：
   - `is_cell_convex_after_move(marginal_cell1, idx_va, candidate_V)`（[`Util:626`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L626)）
   - `is_cell_convex_after_move(marginal_cell2, idx_vb, candidate_V)`
   - 任一不满足 → 打印 `[ConvexCheck STOP]`，`return 0`。
   - 凸性判定由 [`is_polygon_convex`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L583) 实现：按有向面积定方向，检查所有相邻边叉积符号一致。
7. **应用移动**（[行 1110-1111](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1110)）：把 `candidate_V` 写回 `marginal_cell1.points[idx_va]` 与 `marginal_cell2.points[idx_b]`。
8. 打印前后角度对比，`return 1`。

> ⚠️ **marginal 写回不调用 `setVertex()`**（与 inner 不同）。边缘细胞的中心点 / 面积要等本轮 `flush` 时统一重算。

---

## 8. 关键判定函数说明

### 8.1 `judge_by_change` —— 同侧改变法

[`Util:466`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L466)

对受影响的 3 个内角，逐一检查"移动前后，细胞重心与中间点相对于对角线是否发生同侧/异侧翻转"：
- 把点代入直线一般式 `ax+by+c`，用两点的乘积判断同侧（>0）还是异侧（<0）。
- 若乘积 ≈ 0（点落到线上）→ 不可移动。
- 若移动前后**同侧/异侧状态发生改变** → 说明该角越过临界（细胞将变非凸）→ 返回 `True`（不可移动）。

`judge_by_cell` 对移动顶点的前、中、后三个角各调用一次 `judge_by_change`；`judge_by_intersection_cell_blocks` 对 cell1/cell2/cell3 各调用一次 `judge_by_cell`，任一返回 `True` 即判定破坏凸性（flag=-3）。

### 8.2 `judge_sum_inner_angle2` —— 内角平方和约束

[`Util:523`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L523)

计算移动前后该内部顶点对应的 3 个内角（分别属于三个细胞）的**平方和**：
- 移动前 `be_sia = angle1² + angle2² + angle3²`
- 移动后 `af_sia = af_angle1² + af_angle2² + af_angle3²`
- 若 `af_sia > be_sia` → 返回 `True`（内角分布变差，拒绝退火，flag=-4）。

> 数学背景：对固定和（内角和）的角集，平方和越小越接近等角；平方和增大意味着角度分布更不均匀。

### 8.3 `is_point_in_triangle` —— 向量法

[`Util:564`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L564)：用重心坐标 `u, v` 判断点是否在三角形内（`u>0 and v>0 and u+v<1`）。返回 `True` → flag=-1（已接近最优，不移动）。

---

## 9. 返回值与统计

`move_point` 末尾（[行 1238-1253](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L1238)）：

```python
actual_annealed_cell_blocks = marginal_annealing_points + inner_annealing_points
stats = {
    'marginal_points': marginal_annealing_points,   # 实际成功移动的边缘顶点数
    'inner_points':     inner_annealing_points,      # 实际成功移动的内部顶点数
}
return actual_annealed_cell_blocks, stats
```

> 与旧版不同：当前返回的计数就是"实际成功移动的 marginal + inner 之和"，**不再**用 `now_count - best_count - judge_180_count - ...` 的扣减公式；`best_count` / `judge_inner_angle_count` 等变量已不存在，仅保留 `judge_180_count` 用于调试输出。

`Annealer.annealing` 把 `stats` 写入 `self.inner_points / self.marginal_points`，GUI 面板以 `本轮移动数 / 总数` 形式展示。

---

## 10. 关键参数一览

| 参数 | 来源 | 含义 | 默认 |
|------|------|------|------|
| `annealingRate` | GUI 退火速率 | 单步朝目标点移动的比例 | 0.5 |
| `marginal_point_judge` | GUI 边缘开关(0/1) | 是否对边缘顶点退火 | 0（关闭） |
| `inner_angle_sq_guard` | GUI 内角约束(0/1) | 是否启用内角平方和约束(flag=-4) | 1（开启） |
| `times` | GUI 退火次数 | 单次点击执行的轮数 | 1 |
| 边缘角阈值 | 代码内 | 三角形→60°，非三角形→20° | — |

---

## 11. 面积：被动结果与影响因素

**核心结论：整个退火流程不把面积当作约束、目标或阈值——面积只是顶点移动的被动副产品。**

### 11.1 面积的计算与存储

- [`Cell.setArea()`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/utillib/mylib.py#L129)：**鞋带公式**，存到 `cell.area`。
- [`Cell.setVertex()`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/utillib/mylib.py#L99)：只算**中心点（重心）**，**不**算面积——两者是独立的方法。中心点被退火实际使用（构造 `actual_lines` / `pre_best_lines` 射线），面积不被决策使用。

### 11.2 退火决策里面积出现的唯一处

[`is_polygon_convex`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/annealing/annealerUtil.py#L583) 里算有向面积，**只用其符号**判定多边形绕向（CCW/CW）来统一叉积符号判凸性；**数值大小从不参与任何比较**。

### 11.3 每轮刷新

[`CellData.flush()`](file:///c:/Users/kaixu/Desktop/github/SRpython3.14/cell/CellData.py#L35)（形参是 `isListLineOfCell`，非 `isGrow`）→ `topo_grow()` → 对每个细胞调 `setVertex()` + `setArea()`，面积被动重算。但 `CellData` **没有面积聚合字段**（无总面积、ΔA 等）。一轮内移动顶点后 `cell.area` 是旧值，要等 `flush` 才更新；因决策不用面积，无影响。

### 11.4 GUI

退火参数只有 3 个（`annealingRate` / `marginal_point_judge` / `inner_angle_sq_guard`）；统计面板只显示内部顶点数 / 边缘顶点数 / 细胞总数 / 当前轮次，**无任何面积输入、开关或 ΔA 显示**。

### 11.5 间接影响面积的因素

既然面积无人直控，影响都来自位移幅度与是否被允许：

| 影响因素 | 作用机制 |
|---|---|
| `annealingRate` | 步长比例。越大→单轮位移越大→面积变化越剧烈（也越易被守卫拒绝） |
| `times`（退火轮数） | 累积位移。轮数越多→面积累计漂移越大 |
| 凸性守卫（flag -3） | 拒绝会使细胞变非凸的移动→变相抑制大幅面积变化 |
| 内角平方和守卫（flag -4） | 拒绝使角度分布变差的移动→驱动形状趋等角，间接让面积趋向"均衡态"（但非守恒） |
| `marginal_point_judge` | 是否退火边缘顶点。开→边缘细胞面积会变；关→边缘细胞面积冻结 |
| 目标点定义本身 | inner 朝重心、marginal 朝中点——目标定义决定面积漂移的"方向" |

---

## 附：调用关系总览

```
_run_anneal (energy_annealing_main.py:841)
└─ Annealer.annealing (AnnealingGUI.py:53)
   ├─ get_distance_centerpoint_point        (Util:106)
   ├─ get_actual_lines                      (Util:129)
   ├─ get_best_rotate_delta_by_calculation  (Util:189)
   ├─ get_pre_best_lines                    (Util:161)
   ├─ get_intersection_cell_blocks          (Util:242)
   │  ├─ get_point_index_in_cell            (Util:228)
   │  ├─ get_triangle_by_lines              (Util:395)
   │  │  └─ get_crossover_point             (Util:355)
   │  │     └─ get_slope_by_xy              (Util:341)
   │  └─ CellBlock.getTriCentreOfGravity    (mylib.py:218)
   └─ move_point                            (Util:1143)
      ├─ get_all_marginal_points            (Util:669)
      ├─ calculate_marginal_annealing_distance (Util:946)
      │  └─ find_marginal_key_points_new    (Util:710)
      ├─ get_marginal_move_point            (Util:996)
      │  ├─ find_marginal_key_points_new    (Util:710)
      │  ├─ angle_by_three_points          (mylib)
      │  └─ is_cell_convex_after_move       (Util:626)
      │     └─ is_polygon_convex            (Util:583)
      ├─ get_point_of_destination           (Util:321)
      ├─ judge_if_annealing                 (Util:649)
      │  ├─ is_point_in_triangle            (Util:564)
      │  ├─ judge_by_intersection_cell_blocks (Util:413)
      │  │  └─ judge_by_cell                (Util:434)
      │  │     └─ judge_by_change           (Util:466)
      │  └─ judge_sum_inner_angle2          (Util:523)
      │     └─ angle_by_three_points       (mylib)
      └─ set_annealing_options              (Util:643)   ← 初始化时调用
```
