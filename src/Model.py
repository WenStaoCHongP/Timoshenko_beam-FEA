from enum import Enum
import numpy as np

class ElementType(Enum):
    LINEAR = 1    # 线性四边形单元
    QUADRATIC = 2 # 二次四边形单元

class Material:
    def __init__(self, E=1e9, nu=0.25):
        """初始化材料参数
        
        参数:
            E: 弹性模量(Pa) (必须>0)
            nu: 泊松比 (必须0<nu<0.5)
            
        异常:
            ValueError: 如果参数无效
        """
        if E <= 0:
            raise ValueError("弹性模量必须为正数")
        if not 0 < nu < 0.5:
            raise ValueError("泊松比必须在0和0.5之间")
            
        self.E = E       # 弹性模量(Pa)
        self.nu = nu     # 泊松比
        self.G = E/(2*(1+nu))  # 剪切模量
        # 平面应力问题的弹性矩阵
        self.D = self.E/(1-self.nu**2) * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1-self.nu)/2]
        ])
        
        # 验证材料矩阵
        if np.any(np.isnan(self.D)) or np.any(np.isinf(self.D)):
            raise ValueError("材料矩阵计算错误，包含NaN或Inf值")

class Geometry:
    def __init__(self, length=10.0, height=1.0):
        """初始化几何参数
        
        参数:
            length: 长度(m) (必须>0)
            height: 高度(m) (必须>0)
            
        异常:
            ValueError: 如果参数无效
        """
        if length <= 0 or height <= 0:
            raise ValueError("长度和高度必须为正数")
            
        self.L = length  # 梁长度(m)
        self.H = height  # 梁高度(m)
        
    @property
    def area(self):
        """计算横截面积"""
        if self.L <= 0 or self.H <= 0:
            raise ValueError("长度和高度必须为正数")
        return self.L * self.H
        
    @property
    def I(self):
        """计算截面惯性矩"""
        if self.H <= 0:
            raise ValueError("高度必须为正数")
        return self.H**3 / 12  # 截面惯性矩

class Load:
    def __init__(self, P=10.0):
        self.P = P       # 分布载荷(N/m)
        
class FEMResults:
    """存储有限元分析结果的容器类"""
    def __init__(self):
        self.nodes = None
        self.elements = None
        self.displacement = None
        self.stresses = None
        self.errors = None
        
    def validate_results(self):
        """验证FEM结果数据的有效性
        
        异常:
            ValueError: 如果结果数据无效
        """
        if self.displacement is not None:
            if np.any(np.isnan(self.displacement)) or np.any(np.isinf(self.displacement)):
                raise ValueError("位移结果包含NaN或Inf值")
                
        if self.stresses is not None:
            if np.any(np.isnan(self.stresses)) or np.any(np.isinf(self.stresses)):
                raise ValueError("应力结果包含NaN或Inf值")
                
        if self.errors is not None:
            if np.any(np.isnan(self.errors)) or np.any(np.isinf(self.errors)):
                raise ValueError("误差结果包含NaN或Inf值")