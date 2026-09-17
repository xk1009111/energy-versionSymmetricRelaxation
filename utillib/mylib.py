# -*- coding: utf-8 -*-
from decimal import Decimal
import math
from utillib import fittinglib


def get_distance_point_point(p1, p2):
    """
    计算两点之间的距离，支持多种格式

    Args:
        p1: 点1，可以是 Point 对象、列表 [x,y] 或元组 (x,y)
        p2: 点2，格式同 p1

    Returns:
        float: 两点间的距离
    """
    # 处理 p1
    if hasattr(p1, 'x') and hasattr(p1, 'y'):
        x1, y1 = p1.x, p1.y
    elif isinstance(p1, (list, tuple)) and len(p1) >= 2:
        x1, y1 = p1[0], p1[1]
    else:
        raise ValueError(f"不支持的坐标格式: {type(p1)}")

    # 处理 p2
    if hasattr(p2, 'x') and hasattr(p2, 'y'):
        x2, y2 = p2.x, p2.y
    elif isinstance(p2, (list, tuple)) and len(p2) >= 2:
        x2, y2 = p2[0], p2[1]
    else:
        raise ValueError(f"不支持的坐标格式: {type(p2)}")

    distance = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))
    return distance


def angle_by_three_points(p1, vertex, p2):
    """
    根据三个点计算几何夹角，vertex 为角顶点。

    返回值范围为 [0, π]（恒 ≤180°），即纯几何夹角，不区分凹角/优角。
    适用于细胞内角计算：退火凸性守卫已保证细胞恒为凸多边形，
    故每个顶点的真实内角均 ≤180°，无需优角分支。

    Args:
        p1: 角的一条边的端点，list/tuple [x, y]
        vertex: 角顶点，list/tuple [x, y]
        p2: 角的另一条边的端点，list/tuple [x, y]

    Returns:
        float: 夹角（弧度），范围 [0, π]
    """
    v1 = (p1[0] - vertex[0], p1[1] - vertex[1])
    v2 = (p2[0] - vertex[0], p2[1] - vertex[1])
    dot_product = v1[0] * v2[0] + v1[1] * v2[1]
    mag1 = math.hypot(*v1)
    mag2 = math.hypot(*v2)
    if math.isclose(mag1 * mag2, 0):
        return 0
    cos_val = max(-1.0, min(1.0, dot_product / (mag1 * mag2)))
    return math.acos(cos_val)


# 点类
# Point class
class Point:
    x = 0
    y = 0
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return str((self.x, self.y))

# 线段类（直角坐标系）
# Line segment class (rectangular coordinate system)
class Line:
    p1 = Point(0,0)
    p2 = Point(0,0)
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2

    def __str__(self):
        return str([self.p1, self.p2])


# 线段类（极坐标系）
# Line segment class (polar coordinate system)
class Line_in_Polar_Coordinate_System:
    cita = 0
    r = 1
    def __init__(self, cita, r):
        self.cita = cita
        self.r = r

    def __str__(self):
        return str([self.cita, self.r])

    def getX(self):
        # print(math.cos(self.cita))
        return self.r * math.cos(self.cita)

    def getY(self):
        return self.r * math.sin(self.cita)



# 细胞类
# Cell class
class Cell:
    no = 0
    cell_no = 0
    ok = True # 是否最外两层 Is it the outermost two layers
    layer = 0 # 层数，默认为0 The number of layers is 0 by default
    points = []
    center_point = Point(0, 0)
    vx = 0
    vy = 0
    area = 0
    actual_lines = []
    pre_best_lines = []

    # 拟合数据
    data = None

    def __init__(self,points):
        self.ok = True  # 是否不在边缘 Not on the edge
        self.points = points  # 逆时针存储点 Anticlockwise storage point
        self.setNo()
        self.setVertex()
        self.setArea()
        self.layer = 0



    def __str__(self):
        return str(self.center_point)

    # 设置细胞编号
    # Set cell number
    def setNo(self):
        Cell.cell_no+=1
        self.no=Cell.cell_no

    # 计算中心点
    # Computing center
    def setVertex(self):
        points=self.points
        if len(points) <= 2:
            return

        area = Decimal(0.0)
        x, y = Decimal(0.0), Decimal(0.0)
        for i in range(len(points)):
            # print('points[i]', points[i])
            # print('points[i-1]', points[i-1])
            lng = Decimal(points[i][0])
            lat = Decimal(points[i][1])
            nextlng = Decimal(points[i-1][0])
            nextlat = Decimal(points[i-1][1])

            tmp_area = (nextlng*lat - nextlat*lng)/Decimal(2.0)
            area += tmp_area
            x += tmp_area*(lng+nextlng)/Decimal(3.0)
            y += tmp_area*(lat+nextlat)/Decimal(3.0)
        # print('x, area', x, area)
        # print('self.vx, self.vy', self.vx, self.vy)
        if math.fabs(x * area) < 1e-10:
            return
        x = x/area
        y = y/area
        self.vx=float(x)
        self.vy=float(y)
        self.center_point=Point(float(x), float(y))
    # 设置面积
    # Set the area
    def setArea(self):
        s=0
        i=0
        points=self.points
        length=len(points)
        while i<length :
            s+=points[i-1][0]*points[i][1]-points[i][0]*points[i-1][1]
            i+=1
        s/=2
        self.area= s


    # 椭圆拟合部分
    # value 为阈值
    def like_ellipse(self):  # 椭圆拟合 Ellipse fitting
        points = self.points[:]  # 获取细胞点集合 Get cell point collection
        # Since the starting point and the ending point coincide in the cell class, the number of vertices existing
        # in the cell class should be the number of cell edges minus the length of the point set
        # fitting() (V13.0) 直接返回几何参数 [cx, cy, a, b, theta]
        geo = fittinglib.fitting(points, self.center_point, self.area)

        data = self.make_final_ellipse(geo)

        # print("拟合前后点集",self.points, points)
        self.data = data
        return data

    def make_final_ellipse(self, geo):
        # geo 为 fitting() 返回的几何参数 [cx, cy, a, b, theta]；
        # fitting() 已完成代数→几何转换，这里直接映射即可
        # geo is the geometric params [cx, cy, a, b, theta] returned by fitting();
        # fitting() has already done the algebraic→geometric conversion, so just
        # map them directly.
        ellipse_data = {}  # 声明存放椭圆参数的变量 Declare the variable that holds the ellipse parameters

        cx = float(geo[0])
        cy = float(geo[1])
        ellipse_data['cp'] = Point(cx, cy)  # 设置椭圆中心点 Set ellipse center point
        ellipse_data['a'] = float(geo[2])  # 设置椭圆长半轴数据 Set ellipse long half axis data
        ellipse_data['b'] = float(geo[3])  # 设置椭圆短半轴数据 Set ellipse minor axis data
        ellipse_data['angle'] = float(geo[4])  # 设置椭圆与水平轴倾斜角数据 Set the data of inclination angle between ellipse and horizontal axis

        return ellipse_data  # 返回所有椭圆数据信息 Returns all ellipse data information


# 退火细胞块类
# Annealed cell block class
class CellBlock:
    cell1 = object
    index1 = 0
    cell2 = object
    index2 = 0
    cell3 = object
    index3 = 0
    point = Point(0,0)
    triangle = [] # 拟合三角形点集，包含逆时针排列的三个点 The set of fitted triangle points contains three points arranged anticlockwise

    def __init__(self):
        self.cell1 = object
        self.index1 = 0
        self.cell2 = object
        self.index2 = 0
        self.cell3 = object
        self.index3 = 0
        self.point = Point(0,0)

    def setCell1(self, c, i):
        self.cell1 = c
        self.index1 = i

    def setCell2(self, c, i):
        self.cell2 = c
        self.index2 = i

    def setCell3(self, c, i):
        self.cell3 = c
        self.index3 = i

    def setPoint(self, p):
        self.point = p

    #计算三角形重心 Calculate the center of gravity of the triangle
    def getTriCentreOfGravity(self):
        if self.triangle[2] is None:
            raise ValueError("triangle[2] is None, cannot compute centre of gravity")
        x1 = self.triangle[0].x
        y1 = self.triangle[0].y
        x2 = self.triangle[1].x
        y2 = self.triangle[1].y
        x3 = self.triangle[2].x
        y3 = self.triangle[2].y

        xg = (x1+x2+x3) / 3 ;
        yg = (y1+y2+y3) / 3 ;
        #print('in:',xg,yg)
        return Point(xg, yg)

