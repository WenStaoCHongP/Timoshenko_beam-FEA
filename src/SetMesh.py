import numpy as np
from Model import ElementType

class MeshGenerator:
    def __init__(self, case):
        """初始化网格生成器
        
        参数:
            case: CaseSetup对象，包含几何和材料参数
        """
        self.case = case
        
    def generate(self, nx=20, ny=4):
        """生成有限元网格
        
        参数:
            nx: x方向单元数(必须>=1)
            ny: y方向单元数(必须>=1)
            
        返回:
            nodes: 节点坐标数组 (n_nodes x 2)
            elements: 单元连接列表 (n_elements x 4或8)
            
        异常:
            ValueError: 如果输入无效或网格生成失败
        """
        # 详细验证输入参数
        if not isinstance(nx, int) or nx < 1:
            raise ValueError(f"无效的nx值: {nx} (必须为正整数)")
        if not isinstance(ny, int) or ny < 1:
            raise ValueError(f"无效的ny值: {ny} (必须为正整数)")
        if not hasattr(self.case, 'geometry'):
            raise ValueError("案例几何参数未初始化")
        if self.case.geometry.L <= 0 or self.case.geometry.H <= 0:
            raise ValueError("梁长度和高度必须为正数")
            
        # 检查单元数是否合理
        max_elements = 10000  # 最大允许单元数
        if nx * ny > max_elements:
            raise ValueError(f"单元数{nx*ny}超过最大限制{max_elements}")
            
        if self.case.element_type == ElementType.LINEAR:
            return self._generate_linear_mesh(nx, ny)
        else:
            return self._generate_quadratic_mesh(nx, ny)
            
    def _generate_linear_mesh(self, nx, ny):
        """生成线性四边形单元网格
        
        返回:
            nodes: 节点坐标数组 (n_nodes x 2)
            elements: 单元连接列表 (n_elements x 4)
            
        异常:
            ValueError: 如果网格生成失败或质量不合格
        """
        try:
            # 生成节点坐标
            x = np.linspace(0, self.case.geometry.L, nx+1)
            y = np.linspace(0, self.case.geometry.H, ny+1)
            xx, yy = np.meshgrid(x, y)
            nodes = np.column_stack([xx.ravel(), yy.ravel()])
            
            # 验证节点生成
            if len(nodes) != (nx+1)*(ny+1):
                raise ValueError(f"生成的节点数{len(nodes)}与预期{(nx+1)*(ny+1)}不符")
            if nodes.shape[1] != 2:
                raise ValueError("节点坐标维度不正确")
            
            # 生成单元连接
            elements = []
            for j in range(ny):
                for i in range(nx):
                    n0 = i + j*(nx+1)
                    n1 = n0 + 1
                    n2 = n1 + (nx+1)
                    n3 = n0 + (nx+1)
                    
                    # 验证节点索引
                    for node_idx in [n0, n1, n2, n3]:
                        if node_idx < 0 or node_idx >= len(nodes):
                            raise ValueError(f"无效节点索引: {node_idx} (总节点数: {len(nodes)})")
                    
                    elements.append([n0, n1, n2, n3])
            
            # 验证单元连接性
            if len(elements) != nx * ny:
                raise ValueError(f"生成的单元数{len(elements)}与预期{nx * ny}不符")
            
            # 基本网格质量检查
            for elem in elements:
                coords = nodes[elem]
                vec1 = coords[1] - coords[0]
                vec2 = coords[3] - coords[0]
                area = 0.5 * np.abs(np.cross(vec1, vec2))
                if area < 1e-10:
                    raise ValueError(f"单元{elem}面积过小 ({area:.2e}), 可能导致雅可比矩阵奇异")
            
            print(f"生成线性单元网格: {len(nodes)}节点, {len(elements)}单元")
            return nodes, np.array(elements)
            
        except Exception as e:
            raise ValueError(f"线性网格生成失败: {str(e)}")
        
    def _generate_quadratic_mesh(self, nx, ny):
        """生成二次四边形单元网格(优化版本)
        
        返回:
            nodes: 节点坐标数组 (n_nodes x 2)
            elements: 单元连接列表 (n_elements x 8)
            
        节点顺序:
            0-3: 角节点 (顺时针)
            4-7: 边中点 (底,右,顶,左)
            
        异常:
            ValueError: 如果网格生成失败或质量不合格
        """
        try:
            # 生成角节点坐标
            x = np.linspace(0, self.case.geometry.L, nx+1)
            y = np.linspace(0, self.case.geometry.H, ny+1)
            xx, yy = np.meshgrid(x, y)
            nodes = np.column_stack([xx.ravel(), yy.ravel()])
            corner_count = len(nodes)
            
            # 验证角节点生成
            if len(nodes) != (nx+1)*(ny+1):
                raise ValueError(f"生成的角节点数{len(nodes)}与预期{(nx+1)*(ny+1)}不符")
            
            # 生成边中点坐标 (优化存储方式)
            # 水平边中点 (底边和顶边)
            hor_midpoints = []
            for j in range(ny+1):
                for i in range(nx):
                    x_mid = (x[i] + x[i+1])/2
                    if not (x[i] <= x_mid <= x[i+1]):
                        raise ValueError(f"水平边中点x坐标{x_mid}超出范围[{x[i]}, {x[i+1]}]")
                    hor_midpoints.append([x_mid, y[j]])
            
            # 垂直边中点 (左边和右边)
            ver_midpoints = []
            for j in range(ny):
                for i in range(nx+1):
                    y_mid = (y[j] + y[j+1])/2
                    if not (y[j] <= y_mid <= y[j+1]):
                        raise ValueError(f"垂直边中点y坐标{y_mid}超出范围[{y[j]}, {y[j+1]}]")
                    ver_midpoints.append([x[i], y_mid])
            
            # 合并所有节点: 角节点 -> 水平边中点 -> 垂直边中点
            nodes = np.vstack([nodes, hor_midpoints, ver_midpoints])
            
            # 验证合并后的节点
            expected_nodes = corner_count + len(hor_midpoints) + len(ver_midpoints)
            if len(nodes) != expected_nodes:
                raise ValueError(f"合并后节点数{len(nodes)}与预期{expected_nodes}不符")
        
            # 生成单元连接 (优化连接顺序)
            elements = []
            for j in range(ny):
                for i in range(nx):
                    # 角节点索引
                    n0 = i + j*(nx+1)
                    n1 = n0 + 1
                    n2 = n1 + (nx+1)
                    n3 = n0 + (nx+1)
                    
                    # 边中点索引 (优化计算方式)
                    # 水平边中点索引 = corner_count + j*nx + i
                    # 垂直边中点索引 = corner_count + len(hor_midpoints) + j*(nx+1) + i
                    bottom_mid = corner_count + j*nx + i
                    right_mid = corner_count + len(hor_midpoints) + j*(nx+1) + i + 1
                    top_mid = corner_count + (j+1)*nx + i
                    left_mid = corner_count + len(hor_midpoints) + j*(nx+1) + i
                    
                    # 验证所有节点索引
                    all_nodes = [n0, n1, n2, n3, bottom_mid, right_mid, top_mid, left_mid]
                    for node_idx in all_nodes:
                        if node_idx < 0 or node_idx >= len(nodes):
                            raise ValueError(f"无效节点索引: {node_idx} (总节点数: {len(nodes)})")
                    
                    # 检查重复节点
                    if len(set(all_nodes)) != len(all_nodes):
                        raise ValueError(f"单元包含重复节点: {all_nodes}")
                    
                    # 标准8节点单元顺序: 4角点(顺时针) + 4边中点(底,右,顶,左)
                    elements.append(all_nodes)
            
            # 验证单元连接性
            if len(elements) != nx * ny:
                raise ValueError(f"生成的单元数{len(elements)}与预期{nx * ny}不符")
            
            # 网格质量检查
            for elem_idx, elem in enumerate(elements):
                coords = nodes[elem]
                
                # 检查角节点是否形成凸四边形
                vec1 = coords[1] - coords[0]
                vec2 = coords[3] - coords[0]
                area = 0.5 * np.abs(np.cross(vec1, vec2))
                if area < 1e-10:
                    raise ValueError(f"单元{elem_idx}角节点面积过小 ({area:.2e}), 可能导致雅可比矩阵奇异")
                
                # 检查边中点是否在边的中间
                for side in range(4):
                    start = coords[side]
                    end = coords[(side+1)%4]
                    mid = coords[4+side]
                    expected_mid = (start + end)/2
                    dist = np.linalg.norm(mid - expected_mid)
                    if dist > 1e-6:
                        raise ValueError(f"单元{elem_idx}边{side}中点偏离预期位置 (距离={dist:.2e})")
            
            print(f"生成优化二次单元网格: {len(nodes)}节点, {len(elements)}单元")
            print(f"节点分布: {corner_count}角节点, {len(hor_midpoints)}水平边中点, {len(ver_midpoints)}垂直边中点")
            
            return nodes, np.array(elements)
        except Exception as e:
            raise ValueError(f"二次网格生成失败: {str(e)}")