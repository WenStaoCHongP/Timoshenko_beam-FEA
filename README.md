# Timoshenko Beam Finite Element Analysis

![Example Visualization](docs/deformation_plot.png) *(示例：梁的变形和应力分布)*

## 📌 项目概述
基于Timoshenko梁理论的二维有限元分析程序，支持：
- 线性/二次四边形单元
- 静态载荷分析
- 应力和位移可视化
- 与理论解的误差计算

## ✨ 核心功能
- **材料模型**  
  支持各向同性材料（弹性模量E、泊松比ν）
- **单元类型**  
  `LINEAR`：4节点线性单元  
  `QUADRATIC`：8节点二次单元
- **求解器**  
  直接求解器(sparse)和迭代求解器(CG/MINRES)
- **后处理**  
  自动生成：  
  ✓ 变形网格图  
  ✓ 应力/位移云图  
  ✓ 误差分析报告

## 🛠️ 安装指南
### 依赖环境
```bash
# 必需库
pip install numpy scipy matplotlib
# 可选加速库
pip install numba
```

### 获取代码
```bash
git clone https://github.com/your_username/timoshenko-fem.git
cd timoshenko-fem
```

## 🚀 快速开始
### 基本分析
```python
from src.Main import run_analysis

# 默认参数运行（20x4网格，线性单元）
results = run_analysis()
```

### 自定义参数
```python
results = run_analysis(
    nx=30,              # x方向单元数
    ny=6,               # y方向单元数
    element_type='quadratic',  # 单元类型
    length=10.0,       # 梁长度(m)
    height=0.5,         # 梁高度(m) 
    E=210e9,           # 弹性模量(Pa)
    nu=0.3,            # 泊松比
    P=1000             # 分布载荷(N/m)
)
```

## 📊 典型输出
```text
=== 分析结果 ===
最大垂直位移: -3.21e-3 m  
最大弯曲应力(σ_xx): 582.4 Pa  
L2误差(位移): 1.45e-5
```

## 📝 文件结构
src/

├── Core.py # 核心FEM计算

├── SetMesh.py # 网格生成

├── Solver.py # 线性求解器

├── Visualization.py # 结果绘图

└── Main.py # 主入口
## 📜 许可证
MIT License - 详见 [LICENSE](LICENSE) 文件
