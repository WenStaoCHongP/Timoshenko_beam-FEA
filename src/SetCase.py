from Model import Material, Geometry, Load, ElementType
import numpy as np

class CaseSetup:
    def __init__(self, length=10.0, height=1.0, elastic_modulus=1e9, 
                 poisson_ratio=0.25, distributed_load=10.0, 
                 element_type=ElementType.QUADRATIC, debug_mode=False):
        """初始化案例参数
        
        参数:
            length: 梁长度(m) (必须>0)
            height: 梁高度(m) (必须>0)
            elastic_modulus: 弹性模量(Pa) (必须>0)
            poisson_ratio: 泊松比 (必须0<nu<0.5)
            distributed_load: 分布载荷(N/m) (必须>=0)
            element_type: 单元类型(ElementType枚举)
            
        异常:
            ValueError: 如果输入参数无效
        """
        # 验证输入参数
        if length <= 0 or height <= 0 or elastic_modulus <= 0:
            raise ValueError("长度、高度和弹性模量必须为正数")
        if not 0 < poisson_ratio < 0.5:
            raise ValueError("泊松比必须在0和0.5之间")
        if distributed_load < 0:
            raise ValueError("分布载荷不能为负")
            
        self.material = Material(E=elastic_modulus, nu=poisson_ratio)
        self.geometry = Geometry(length=length, height=height)
        self.load = Load(P=distributed_load)
        self.element_type = element_type
        
        # 验证材料矩阵
        if np.any(np.isnan(self.material.D)) or np.any(np.isinf(self.material.D)):
            raise ValueError("材料矩阵计算错误，包含NaN或Inf值")
        
    def print_parameters(self):
        """打印当前案例参数"""
        print("\n=== 输入参数 ===")
        print(f"梁长度(L): {self.geometry.L} m")
        print(f"梁高度(H): {self.geometry.H} m")
        print(f"弹性模量(E): {self.material.E} Pa")
        print(f"泊松比(nu): {self.material.nu}")
        print(f"分布载荷(P): {self.load.P} N/m")
        print(f"单元类型: {self.element_type.name}")
        print(f"剪切模量(G): {self.material.G:.2e} Pa")
        
    def update_parameters(self, **kwargs):
        """更新案例参数"""
        for key, value in kwargs.items():
            if hasattr(self.material, key):
                setattr(self.material, key, value)
                # 如果更新了E或nu，需要重新计算G
                if key in ['E', 'nu']:
                    self.material.G = self.material.E/(2*(1+self.material.nu))
            elif hasattr(self.geometry, key):
                setattr(self.geometry, key, value)
            elif hasattr(self.load, key):
                setattr(self.load, key, value)
            elif key == 'element_type':
                self.element_type = value
                
    def get_material_matrix(self):
        """获取平面应力问题的材料矩阵D"""
        return self.material.D