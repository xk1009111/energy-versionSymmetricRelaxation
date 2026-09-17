# 中英文参数对照表

## 界面标题与标签

| 键名 | 中文 | 英文 |
|------|------|------|
| window_title | 细胞退火工具 | Cell Annealing (Relaxation) Tool |
| main_panel_title | 细胞退火参数配置 | Cell Annealing Configuration |
| font_label | 字体: | Font: |
| reset | 重置 | Reset |

## 维诺图初始化

| 键名 | 中文 | 英文 |
|------|------|------|
| voronoi_init | 维诺图初始化 | Voronoi Initialization |
| voronoi_type | 维诺图类型: | Voronoi Type: |
| voronoi_grid | 随机扰动 n×n | hexagon random-disordered Voronoi (n×n) |
| voronoi_random | 均匀随机 | Uniform Random Voronoi |
| voronoi_seed | n（种子值）: | n (Seed Value): |
| voronoi_param | k（扰动 0~1）: | k (Irregularity 0~1): |
| random_seed | 种子数（平方数）: | Number of Seeds (Square Number): |
| generate_voronoi | 生成维诺图并初始化 | Generate Voronoi & Initialize |

## 细胞退火器

| 键名 | 中文 | 英文 |
|------|------|------|
| cell_annealer | 细胞退火器 | Cell Annealer |
| anneal_rate | 退火因子 (γ): | Relaxation factor (γ): |
| anneal_times | 退火轮次: | Iterations: |
| anneal_edge | 边缘顶点参与退火: | Include Marginal Vertices in Annealing: |
| yes | 是 | Yes |
| no | 否 | No |
| inner_vertices_hint | 内部顶点：若移动后对应3个内角平方和增大则禁止退火 | Inner Vertices: Forbid annealing if sum of squares of 3 interior angles increases after movement |
| execute_annealing | 执行退火 | Execute Annealing |

## 辅助功能

| 键名 | 中文 | 英文 |
|------|------|------|
| auxiliary_lines | 辅助线 | Auxiliary Lines |
| ellipse_fitting | 椭圆拟合 | Ellipse Fitting |
| export_data | 输出数据 | Export Data |
| undo | 回退 | Undo |

## 统计信息

| 键名 | 中文 | 英文 |
|------|------|------|
| statistics | 统计信息 | Statistics |
| inner_vertices_stat | 内部顶点总数/上轮退火数 | Total Inner Vertices / Moved Vertices (Previous Iteration) |
| edge_vertices_stat | 边缘顶点总数/上轮退火数 | Total Marginal Vertices / Moved Vertices (Previous Iteration) |
| total_cells | 细胞总数 | Total Cells |

## 工具栏

| 键名 | 中文 | 英文 |
|------|------|------|
| toolbar_rect_zoom | 框选放大 | Zoom |
| toolbar_pan | 移动 | Pan |
| toolbar_reset_view | 重置视图 | Reset View |
| toolbar_save | 保存 | Save |
| font_reset | 重置 | Reset |

## 错误与提示

| 键名 | 中文 | 英文 |
|------|------|------|
| error_missing_params | 参数缺失 | Missing Parameters |
| error_fill_voronoi | 请填写完整的维诺图参数 | Please fill in complete Voronoi parameters |
| error | 错误 | Error |
| error_seed_valid | 种子值需为数字且大于3 | Seed value must be a number greater than 3 |
| error_irregularity_range | Irregularity must be between 0.0 and 1.0 | Irregularity must be between 0.0 and 1.0 |
| error_irregularity_number | Irregularity must be a number, current value: {k_val} | Irregularity must be a number, current value: {k_val} |
| error_fill_seeds | 填写种子数 | Input the number of seeds |
| confirm | 确认 | Confirm |
| confirm_overwrite | 检测到已有细胞图像，确定要生成新的初始图像并覆盖吗？ | Detected existing cell image, confirm to generate a new initial image and overwrite? |
| success | 成功 | Success |
| success_init | 系统初始化成功！ | System initialization successful! |
| error_init | 初始化失败：{error} | Initialization failed: {error} |
| warning | 提示 | Warning |
| warning_generate_first | 请先生成维诺图 | Please generate Voronoi first |
| error_param | 参数错误：{error} | Parameter error: {error} |
| complete | 完成 | Complete |
| complete_annealing | 退火完成。 | Annealing completed. |
| error_display_rays | 显示最优中心射线失败：{error} | Failed to display optimal rays: {error} |
| error_ellipse | 椭圆拟合失败：{error} | Ellipse fitting failed: {error} |

## 文件操作

| 键名 | 中文 | 英文 |
|------|------|------|
| save_data | 保存数据 | Save Data |
| excel_workbook | Excel 工作簿 | Excel Workbook |
| all_files | 所有文件 | All Files |
| error_export | 导出数据失败 | Export data failed |
| error_export_prepare | 准备导出数据失败：{error} | Failed to prepare export data: {error} |
| error_export_not_found | 未找到导出的临时文件 | Exported temporary file not found |
| success_export | 数据已导出。 | Data exported. |
| save_image | 保存图片 | Save Image |
| eps_vector | EPS 矢量图 | EPS Vector |
| png_image | PNG 图片 | PNG Image |
| error_save_file | 文件未成功生成：{file_path} | File not generated successfully: {file_path} |
| success_save | 保存成功 | Save Successful |
| success_save_path | 图片已保存到：{file_path} | Image saved to: {file_path} |
| error_save | 保存失败 | Save Failed |
| error_export_image | 导出失败：{error} | Export failed: {str_e} |
| error_save_permission | 保存失败：没有权限写入文件 {file_path}，请尝试保存到其他位置（如桌面） | Save failed: No permission to write to {file_path}, Please try saving to another location (e.g., Desktop) |

## 回退功能

| 键名 | 中文 | 英文 |
|------|------|------|
| error_no_snapshot | 没有可回退的快照 | No snapshot to undo |
| confirm_undo | 确定要回退到上次操作前的状态吗？ | Confirm to undo to the state before the last operation? |
| success_undo | 已回退。 | Undone. |
| error_undo | 回退失败。 | Undo failed. |

## 语言切换

| 键名 | 中文 | 英文 |
|------|------|------|
| lang_chinese | 中文 | 中文 |
| lang_english | English | English |

## Excel 表头 - 椭圆拟合数据

| 中文 | 英文 |
|------|------|
| 细胞序号 | Cell ID |
| 椭圆形心 | Ellipse Centroid |
| 长半轴 | Major Semi-axis |
| 短半轴 | Minor Semi-axis |
| 短半轴倾角 | Inclination Angle of the Minor Semi-axis |
| 相邻细胞边数和 | Sum of Neighbor Edge Numbers |
| 细胞边数 | Cell Edge Number |
| 细胞周长 | Cell Perimeter |
| 细胞面积 | Cell Area |
| 最初母细胞层 | Original Cell Layer |
| 当前细胞层 | Current Cell Layer |

## Excel 表头 - 边角数据

| 中文 | 英文 |
|------|------|
| 细胞序号 | Cell ID |
| 边数 | Edge Number |
| 相邻细胞边数和 | Sum of Neighbor Edge Numbers |
| 内角 | Interior Angle |
| 夹边1 | Adjacent Edge 1 |
| 夹边2 | Adjacent Edge 2 |
| 最初母细胞层 | Original Cell Layer |
| 当前细胞层 | Current Cell Layer |

## Excel 表头 - 边缘数据

| 中文 | 英文 |
|------|------|
| 细胞序号 | Cell ID |
| 边缘边长(ME) | Margin Edge Length (ME) |
| 边缘角1(MA1) | Margin Angle 1 (MA1) |
| 边缘角2(MA2) | Margin Angle 2 (MA2) |
| 最初母细胞层 | Original Cell Layer |
| 当前细胞层 | Current Cell Layer |

## Excel 表头 - 边缘数据详情

| 中文 | 英文 |
|------|------|
| 细胞序号 | Cell ID |
| 参数类型 | Parameter Type |
| 数值 | Value |
| 单位 | Unit |
| 位置索引(点/边) | Position Index (Point/Edge) |
| 细胞层数 | Cell Layer |

## Excel 表头 - 附加数据

| 中文 | 英文 |
|------|------|
| 是否分裂 | Divided |
| 子细胞面积比值 | Daughter Cell Area Ratio |

## Excel 表头 - 分类标题

| 中文 | 英文 |
|------|------|
| 拟合椭圆数据 | Fitted Ellipse Data |
| 细胞基本数据 | Cell Basic Data |
