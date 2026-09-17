# 细胞退火工具 — Cell Relaxation (Annealing) Tool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.14+](https://img.shields.io/badge/python-3.14+-blue.svg)](https://www.python.org/downloads/)

> **[English](README_EN.md)**

本项目是**二维细胞网络对称松弛（退火）算法**的 Python 实现。从完整版 CellSymRelax 项目中剥离了细胞增殖、拓扑变换等扩展功能，仅保留核心松弛功能：维诺图生成、顶点平衡、椭圆拟合和数据导出。

基于以下研究：

- Xu K., Weng L., Wang Z., Lian Y., Huang B. (2026). *A symmetric relaxation method for entire two-dimensional cellular networks and its implications*. (arXiv: XuSR20260616)
- Xu K. (2021). *A geometry-based relaxation algorithm for equilibrating a trivalent polygonal network in two dimensions and its implications*. Philosophical Magazine, 101(14), 1632-1653.

***

## 功能特性

### 维诺图生成

生成修剪的 Voronoi 网络：

- **随机扰动 (n×n)** — 随机扰动正六边形网格的种子点坐标，扰动（菱形采样法）范围限定在正六边形内，参数 *k* 控制不规则度（0=正六边形，1=高度无序）。
- **均匀随机** — 基于均匀随机点。

### 对称松弛（退火）

通过几何对称性平衡内部顶点和边缘顶点：

- **内部顶点** — 由关联细胞的中心角对称性驱动，可选内角对称限制（即内角倾向于120°）。
- **边缘顶点** — 沿边界边移动，目标边缘角度 90°；也受中心角对称影响。
- **可配置** — 退火因子（可看作顶点移动速度）、边缘顶点参与开关、内角平方和守卫。

### 椭圆拟合

- **R-LMG**（`conicfit` 包）— 高精度几何拟合，三角形/四边形的首选方法。
- **代数最小二乘**（numpy SVD）— 快速拟合，五边及以上的首选方法。
- 要求椭圆和n边形的面积比值在1和3之间，首选失败则用代数法基于真实点+镜像点（共2n）拟合；次选失败则用R-LMG基于真实点+镜像点（共2n）拟合；还失败则用基于真实点+镜像点（共2n）的保底最小外接矩形MBR。

### 数据导出

导出为 Excel（`.xlsx`）：

- **椭圆及细胞周长和面积** — 细胞序号、椭圆形心、长半轴、短半轴、长半轴倾角、相邻细胞边数和、细胞边数、细胞周长、细胞面积、最初母细胞层、当前细胞层（层数定义：最外层边缘细胞的层数为1，从外往内依次加1）。
- **边角数据** — 细胞序号、边数、相邻细胞边数和、内角、夹边1、夹边2、最初母细胞层、当前细胞层（每个顶点一行）。
- **ME/MA 统计** — 细胞序号、边缘边长(ME)、边缘角1(MA1)、边缘角2(MA2)、最初母细胞层、当前细胞层（仅边缘细胞，每个细胞一行）。

### 可视化

- 实时细胞网络显示，按层数着色。
- 辅助线叠加显示形心和顶点连线，以及每个顶点对应的三条射线组成的三角形。
- 椭圆叠加显示拟合结果。
- 图片导出为 EPS（矢量图）和 PNG。

***

## 环境要求

- **Python 3.14+**（由 [uv](https://docs.astral.sh/uv/) 托管，安装脚本会自动安装 Python 3.14.6）
- Python 依赖见 [requirements.txt](requirements.txt)
- **R** 语言环境（用于 R-LMG 椭圆拟合，通过 `rpy2` 调用），可选但推荐
- R 包：`conicfit`、`sp`、`shotGroups`

***

## 安装

> 前置：先安装 [uv](https://docs.astral.sh/uv/getting-started/installation/)（极快的 Python 包管理器，Windows / macOS / Linux 均支持）。

### 一键安装

**Windows：**

```bash
setup.bat
```

**Linux / macOS：**

```bash
chmod +x setup.sh
./setup.sh
```

脚本将自动：

1. 检查 uv，通过 `uv python install 3.14.6` 安装 Python 并创建 `.venv`。
2. 通过 `uv pip install -r requirements.txt` 安装所有 Python 包。
3. 检测 `R_Dist/` 便携版目录并配置。
4. 未找到则回退到系统 R。
5. 仍未找到则**自动下载安装 R**（Windows/macOS）。
6. 安装所需 R 包（`conicfit`、`sp`、`shotGroups`）。

### 手动安装

```bash
# 1. 用 uv 安装 Python 3.14.6 并创建虚拟环境
uv python install 3.14.6
uv venv .venv --python 3.14.6
uv pip install -r requirements.txt

# 2. 安装 R 及 R 包（可选，用于 R-LMG 拟合）
#    从 https://www.r-project.org/ 安装 R，然后：
#    Rscript -e "install.packages(c('conicfit', 'sp', 'shotGroups'))"

# 3. 运行
uv run python energy_annealing_main.py
```

### 打包为 EXE

可直接运行 [build_exe.bat](build_exe.bat)，或手动执行：

```bash
uv pip install pyinstaller
uv run pyinstaller --onefile --windowed --name "CellAnnealing" ^
    --add-data "annealing;annealing" ^
    --add-data "cell;cell" ^
    --add-data "utillib;utillib" ^
    --add-data "initVoronoi.py;." ^
    --hidden-import "scipy.spatial" ^
    --hidden-import "scipy.optimize" ^
    --hidden-import "openpyxl" ^
    --hidden-import "pyenvelope" ^
    energy_annealing_main.py
# 输出：dist/CellAnnealing.exe
```

> 打包后将 `R_Dist/` 放在 exe 同目录下即可启用 R-LMG 拟合。

***

## 文件结构

```
only_annealing/
├── energy_annealing_main.py      # 主入口（Tkinter GUI）
├── requirements.txt            # Python 依赖
├── setup.bat                   # Windows 一键安装脚本（uv）
├── setup.sh                    # Linux/macOS 一键安装脚本（uv）
├── build_exe.bat               # PyInstaller 打包脚本（uv）
├── .gitignore
├── .gitattributes
├── LICENSE                     # MIT 许可证
├── README_EN.md                # 英文文档
├── README_CN.md                # 中文文档
├── ANNEALING_FLOW_CN.md        # 退火算法流程说明（中文）
├── PARAMS_CN_EN.md             # 参数中英文对照表
│
├── initVoronoi.py              # 维诺图初始化
│
├── annealing/                  # 松弛（退火）算法
│   ├── AnnealingGUI.py         # 退火器核心（参数化）
│   └── annealerUtil.py         # 退火工具函数
│
├── cell/                       # 细胞数据与统计
│   └── CellData.py             # 细胞数据管理
│
├── scripts/                    # 辅助脚本（数据收集/测试）
│   ├── collect_data.py         # 数据收集脚本
│   ├── collect_iteration_data.py  # 迭代数据收集
│   ├── test_area.py            # 面积计算单元测试
│   └── DATA_COLLECTION_PLAN.md # 批量数据收集方案
│
└── utillib/                    # 工具库
    ├── mylib.py                # 数据结构（Cell、Point、Line 等）
    ├── fittinglib.py           # 椭圆拟合（R-LMG + 最小二乘）
    ├── exportUtils.py          # Excel 数据导出
    ├── layerMarker.py          # 细胞层数标记
    └── i18n.py                 # 中英文国际化
```

***

## R 环境配置

`setup.bat` / `setup.sh` 按以下优先级处理 R 环境：

| 优先级 | 来源 | 检测方式 |
|---|---|---|
| 1（最佳） | 便携版 `R_Dist/` | 本地存在 `R_Dist/bin/R` |
| 2 | 系统 R | `R_HOME` 环境变量或 `R` 在 PATH 中 |
| 3（保底） | 自动下载 | 脚本自动下载安装 R |

无 R 环境时，椭圆拟合自动降级为纯 Python（numpy/scipy 最小二乘），功能不受影响，但对三角形/四边形的拟合精度可能略有不同。

> **注意：** `R_Dist/` 已加入 `.gitignore`，不会上传到 Git。克隆后运行 `setup.bat` 即可自动配置，或从完整分发版复制 `R_Dist/`。

***

## 使用说明

> **注意：** 普通计算机运行退火程序时，细胞数量建议不超过 400，否则会很慢或死机；关闭图形可视化界面后，可适度提升细胞上限。

```bash
uv run python energy_annealing_main.py
```

GUI 提供以下功能：

1. **维诺图初始化** — 配置网络类型和参数，生成初始网络。
2. **退火** — 设置退火速率、边缘顶点是否参与、内角平方和守卫，执行单次或多次迭代。
3. **辅助线** — 显示最优射线，叠加显示形心和顶点连线，以及每个顶点对应的三条射线组成的三角形。
4. **椭圆拟合** — 对所有细胞进行椭圆拟合并叠加显示。
5. **数据导出** — 将多边形几何和拓扑参数数据保存为三个 Excel 文件。
6. **图片导出** — 通过 matplotlib 工具栏保存当前视图为 EPS 或 PNG。

***

## 引用

如果本项目用于学术研究，请引用：

```
Xu K., Weng L., Wang Z., Lian Y., Huang B. (2026). A symmetric relaxation method for entire two-dimensional cellular networks and its implications. arXiv: XuSR20260616.
```

```
Xu K. (2021). A geometry-based relaxation algorithm for equilibrating a trivalent polygonal network in two dimensions and its implications. Philosophical Magazine, 101(14), 1632-1653.
```

***

## 许可证

本项目采用 MIT 许可证 — 详见 [LICENSE](LICENSE) 文件。

***

## 贡献者

### 项目负责人 & 核心算法设计

- **许凯 (Kai Xu)** — 项目发起人，所有核心算法的设计师，代码架构规划与设计。在 Trae 协助下开发了中英文双语 GUI，并优化了代码。

### 核心技术贡献者

- **史国威 (Guowei Shi)** — 实现了初始基础代码和原始模拟退火原型。
- **翁力凡 (Lifan Weng)、王子涵 (Zihan Wang)、连钰洋 (Yuyang Lian)、黄斌 (Bin Huang)** — 进行了后续代码重构、计算性能优化、GUI功能改进和补充国际化支持。

***

## 联系方式

**许凯 (Kai Xu)**

- 邮箱：<kaixu@jmu.edu.cn> / <kxu2013@gmail.com>
- ORCID: [0000-0002-1341-1525](https://orcid.org/0000-0002-1341-1525)
- 单位：集美大学 水产学院

如有问题、Bug 报告或合作意向，请联系通讯作者。
