# -*- coding: utf-8 -*-


'''
    周边探索法。
    Peripheral exploration method
    :param cells: 细胞集，Cell对象列表 Cell set, cell object list
    :param N: 细胞维诺图种子数  Seed number of cell Vinot map
'''
def layer_mark2(cells, N):
    """
    周边探索法（改进版：确保所有细胞都被标记，从第1层开始）。
    通过逐层向外扩散的方式标记所有细胞，避免出现"第0层"的未标记细胞。
    当一个细胞与多个不同层的细胞相邻时，取最小相邻层号+1。
    :param cells: 细胞集，Cell对象列表 Cell set, cell object list
    :param N: 细胞维诺图种子数（保留参数，不再用于限制循环次数） Seed number of cell Vinot map (reserved, no longer used for loop bound)
    """
    # 初始时，边缘细胞已由 setting_layer() 标记为第1层
    # 使用BFS逐层标记：第2层只标记与第1层相邻的，第3层只标记与第2层相邻的...
    
    # 获取第1层的细胞作为初始队列
    current_layer_cells = [c for c in cells if c.layer == 1]
    current_layer = 1
    
    while current_layer_cells:
        next_layer_cells = []
        next_layer = current_layer + 1
        
        for c in current_layer_cells:
            # 找到所有与当前层细胞相邻的未标记细胞
            for tc in cells:
                if tc.layer != 0:  # 已标记的跳过
                    continue
                
                # 判断是否相邻（共享至少一个点）
                is_neighbor = False
                for cp in c.points:
                    if cp in tc.points:
                        is_neighbor = True
                        break
                
                if is_neighbor:
                    tc.layer = next_layer
                    next_layer_cells.append(tc)
        
        current_layer_cells = next_layer_cells
        current_layer = next_layer


