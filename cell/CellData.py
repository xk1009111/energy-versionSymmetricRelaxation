class CellData:
    cells = []
    length = 0
    lineOfCell = {}

    def __init__(self, cells):
        self.init(cells)

    # 细胞集参数初始化
    def init(self, cells):
        self.cells = cells
        self.length = len(cells)

        for i in cells:
            i.setArea()

        self.list_line_of_cell()

    def topo_grow(self):
        for c in self.cells:
            c.setVertex()
            c.setArea()
        self.length = len(self.cells)

    def flush(self, isListLineOfCell=True):
        self.length = len(self.cells)
        self.topo_grow()

        if isListLineOfCell:
            self.list_line_of_cell()

        #刷新边缘细胞
        self.setting_layer()
        self.length = len(self.cells)

    """
        生成细胞边归属字典，用于记录每条线段归属于哪个细胞
        A cell edge attribution dictionary is generated to record which cell each line segment belongs to
        :param cells: 细胞集，Cell对象列表 Cell set, cell object list
        :return lineOfCell: 细胞边归属字典，记录每条线段归属于哪个细胞 Cell edge attribution dictionary, recording which cell each line segment belongs to
    """
    def list_line_of_cell(self):
        # 将初始图的每条边的归属细胞记录
        # Record the attribution cell of each edge of the initial diagram
        # 注意：这里存储的是单个细胞，不是列表
        # 对于共享边，后遍历的细胞会覆盖前面的
        # 但这不影响 setting_layer() 的判断，因为它只检查边是否存在
        lineOfCell = {}
        for j in range(len(self.cells)):
            cell = self.cells[j]
            for i in range(len(cell.points)):
                # 使用固定的顺序：(i-1, i)，确保与 setting_layer() 一致
                string = '{0}-{1}-{2}-{3}'.format(cell.points[i - 1][0], cell.points[i - 1][1], cell.points[i][0],
                                                  cell.points[i][1])
                # 存储边和对应的细胞（后遍历的会覆盖前面的，但不影响判断）
                lineOfCell[string] = cell
        self.lineOfCell = lineOfCell

    def setting_layer(self):
        # 初始化为正六边形分裂则不排除
        # Regular hexagon splitting is not ruled out
        # 重新构建边计数，统计每条边被多少个细胞共享
        # 边缘边只被一个细胞拥有，内部边被两个细胞共享
        # 关键：需要规范化边的表示，确保同一条边（无论方向）都使用相同的字符串

        def normalize_point(p):
            """将点转换为元组，确保类型一致"""
            if isinstance(p, (list, tuple)):
                return (float(p[0]), float(p[1]))
            elif hasattr(p, '__getitem__'):
                return (float(p[0]), float(p[1]))
            else:
                return (float(p.x), float(p.y))

        def normalize_edge(p1, p2):
            """规范化边：总是使用较小的点作为起点"""
            p1_tuple = normalize_point(p1)
            p2_tuple = normalize_point(p2)
            if p1_tuple < p2_tuple:
                return '{0}-{1}-{2}-{3}'.format(p1_tuple[0], p1_tuple[1], p2_tuple[0], p2_tuple[1])
            else:
                return '{0}-{1}-{2}-{3}'.format(p2_tuple[0], p2_tuple[1], p1_tuple[0], p1_tuple[1])

        edge_count = {}  # 统计每条边被多少个细胞拥有（使用规范化的边表示）
        for cell in self.cells:
            for i in range(len(cell.points)):
                p1 = cell.points[i - 1]
                p2 = cell.points[i]
                string = normalize_edge(p1, p2)
                edge_count[string] = edge_count.get(string, 0) + 1

        # 判断边缘细胞：如果细胞的某条边只被一个细胞拥有（计数为1），则该细胞是边缘细胞
        for cell in self.cells:
            flag = False
            # 若边没有相邻细胞，则此细胞在边缘
            # If there are no adjacent cells on the edge, the cell is on the edge
            for i in range(len(cell.points)):
                p1 = cell.points[i - 1]
                p2 = cell.points[i]
                string = normalize_edge(p1, p2)

                # 如果边只被一个细胞拥有（计数为1），说明是边缘边
                # 内部边被2个细胞共享，所以计数为2
                if edge_count.get(string, 0) == 1:
                    flag = True
                    break
            # 如果是边缘细胞，将周围细胞全部排除
            # If it is a marginal cell, exclude all surrounding cells
            if flag:
                cell.ok = False
                cell.layer = 1
