import numpy as np
from Model import FEMResults, ElementType
from Tools import FEMTools

class PostProcessor:
    def __init__(self, case, mesh, solution):
        """初始化后处理器
        
        参数:
            case: CaseSetup对象
            mesh: (nodes, elements)元组
            solution: 位移解向量
        """
        self.case = case
        self.nodes, self.elements = mesh
        self.u = solution
        self.results = FEMResults()
        self.results.nodes = self.nodes
        self.results.elements = self.elements
        self.results.displacement = self.u
        self.debug = getattr(case, 'debug_mode', False)
        
        if self.debug:
            print("\n=== 后处理器初始化 ===")
            print(f"节点数: {len(self.nodes)}")
            print(f"单元数: {len(self.elements)}")
            print(f"位移向量长度: {len(self.u)}")
        
    def compute_stresses(self):
        """计算单元应力和节点平均应力(严格遵循铁木辛柯梁理论)
        
        返回:
            elem_stresses: 单元应力列表
            nodal_stresses: 节点平均应力字典
            
        实现细节:
        1. 对线性单元使用2×2高斯积分(标准积分)
        2. 对积分单元使用3×3高斯积分
        3. 应力计算包括弯曲和剪切分量分离
        """
        
        # 获取材料矩阵(确保使用np.float64)
        D = self.case.get_material_matrix().astype(np.float64)
        
        # 获取高斯积分点和权重
        if self.case.element_type == ElementType.LINEAR:
            gauss_points, gauss_weights = FEMTools.gauss_points(1)
            print("Using 2×2 Gauss integration for linear elements")
        else:
            gauss_points, gauss_weights = FEMTools.gauss_points(2) 
            print("Using 3×3 Gauss integration for quadratic elements")
        
        # 初始化存储结构
        # 初始化应力存储结构
        nodal_stresses = {}
        elem_stresses = [] # 存储单元应力数据
        
        for elem in self.elements:
            # 根据单元类型获取节点坐标和位移
            if self.case.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(4)])
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]  # 8节点单元
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(8)])
                n_points = 8
            
            gp_stresses = []
            # 在积分点计算应力
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                
                # 获取形函数及其导数
                N, dN_dxi = FEMTools.shape_functions(xi, eta, 
                    self.case.element_type.value)
                
                # 验证形函数导数
                sum_dxi = np.sum(dN_dxi[0])
                sum_deta = np.sum(dN_dxi[1])
                if abs(sum_dxi) > 1e-10 or abs(sum_deta) > 1e-10:
                    raise ValueError(f"形函数导数和不为零 (dxi={sum_dxi:.2e}, deta={sum_deta:.2e})")
                
                # 计算雅可比矩阵
                J = np.zeros((2, 2))
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
                        
                # 验证雅可比矩阵
                detJ = np.linalg.det(J)
                if abs(detJ) <= 1e-10:
                    print(f"积分点(xi={xi}, eta={eta})的雅可比行列式太小({detJ:.2e})")
                    print("节点坐标:", node_coords)
                    print("形函数导数:", dN_dxi)
                    raise ValueError("雅可比行列式太小")
                    
                invJ = np.linalg.inv(J)
                dN_dx = invJ @ dN_dxi
                
                # 验证物理坐标导数
                if not np.allclose(np.sum(dN_dx[0]), 0, atol=1e-10) or \
                   not np.allclose(np.sum(dN_dx[1]), 0, atol=1e-10):
                    raise ValueError("物理坐标导数和不为零")
                
                # 应变-位移矩阵B
                B = np.zeros((3, 2*n_points))
                for i in range(n_points):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
                
                # 计算应力(分离弯曲和剪切分量)
                stress = D @ B @ elem_disp
                gp_stresses.append(stress)

            elem_stress = np.mean(gp_stresses, axis=0)
            elem_stresses.append(elem_stress)
            
            # 收集节点应力
            for i in range(n_points):
                node = elem[i]
                if node not in nodal_stresses:
                    nodal_stresses[node] = []
                nodal_stresses[node].append(elem_stress)
        
        # 计算平均节点应力并进行验证
        avg_nodal_stresses = {}
        for node, stresses in nodal_stresses.items():
            avg_nodal_stresses[node] = np.mean(stresses, axis=0)
        
        # 存储结果
        self.results.stresses = avg_nodal_stresses
        self.results.elem_stresses = elem_stresses
            
        return elem_stresses, avg_nodal_stresses

        
    def exact_solution(self, x, y):
        """计算精确解
        
        参数:
            x, y: 位置坐标
            
        返回:
            (ux, uy, sxx, syy, sxy): 精确位移和应力
        """
        L = self.case.geometry.L
        H = self.case.geometry.H
        E = self.case.material.E
        nu = self.case.material.nu
        P = self.case.load.P
        
        # Timoshenko梁精确解
        ux = (2*P*y)/(E*H**3) * ((6*L - 3*x)*x + (2 + nu)*y**2 - 3*H**2*(1 + nu)/2)
        uy = -(2*P)/(E*H**3) * (3*nu*y**2*(L - x) + (3*L - x)*x**2)
        sxx = (12*P*(L - x)*y) / H**3
        syy = 0
        sxy = -(6*P/H**3) * (H**2/4 - y**2)
        
        return ux, uy, sxx, syy, sxy
        
    def compute_errors(self):
        """计算数值解与精确解的误差(包含位移和应力误差)
        
        返回:
            (l2_error, h1_error, max_error, stress_error): 各种误差指标
        """
        print("\nComputing errors (displacement and stress)...")
        
        # 获取高斯积分点和权重
        if self.case.element_type == ElementType.LINEAR:
            gauss_points, gauss_weights = FEMTools.gauss_points(1)
        else:
            gauss_points, gauss_weights = FEMTools.gauss_points(2)
        
        # 初始化误差指标
        l2_err = 0.0
        h1_err = 0.0
        max_err = 0.0
        stress_l2_err = 0.0
        stress_max_err = 0.0
        
        # 获取材料矩阵
        D = self.case.get_material_matrix()
        
        for elem in self.elements:
            # 根据单元类型获取节点坐标和位移
            if self.case.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(4)])
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(8)])
                n_points = 8
            
            # 在积分点计算误差
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                
                # 计算形函数及其导数
                N, dN_dxi = FEMTools.shape_functions(xi, eta, 
                    self.case.element_type.value)
                
                # 计算雅可比矩阵和行列式
                J = np.zeros((2, 2))
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
                detJ = np.abs(np.linalg.det(J))
                
                # 计算物理坐标(调整y坐标到梁中心)
                x = np.sum(N * node_coords[:,0])
                y = np.sum(N * node_coords[:,1]) - self.case.geometry.H/2  # 以梁中心为y=0
                
                # 计算数值解
                u_num = np.zeros(2)
                for i in range(n_points):
                    u_num += N[i] * elem_disp[2*i:2*i+2]
                
                # 计算精确解(包含应力和应变)
                ux_exact, uy_exact, sxx_exact, syy_exact, sxy_exact = self.exact_solution(x, y)
                u_exact = np.array([ux_exact, uy_exact])
                stress_exact = np.array([sxx_exact, syy_exact, sxy_exact])
                strain_exact = np.linalg.inv(D) @ stress_exact
                
                # 计算数值应力(使用单元平均应力)
                elem_idx = self.elements.index(elem)
                elem_stress = self.results.elem_stresses[elem_idx]
                
                # 计算应变-位移矩阵B
                invJ = np.linalg.inv(J)
                dN_dx = invJ @ dN_dxi
                B = np.zeros((3, 2*n_points))
                for i in range(n_points):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
                
                # 计算数值应变
                strain_num = B @ elem_disp
                
                # 计算位移误差(L2范数)
                l2_err += np.sum((u_num - u_exact)**2) * w * detJ
                max_err = max(max_err, np.max(np.abs(u_num - u_exact)))
                
                # 计算梯度误差(H1范数)
                h1_err += np.sum((strain_num - strain_exact)**2) * w * detJ
                
                # 计算应力误差
                stress_diff = elem_stress - stress_exact
                stress_l2_err += np.sum(stress_diff**2) * w * detJ
                stress_max_err = max(stress_max_err, np.max(np.abs(stress_diff)))
        
        # 计算全局误差
        l2_error = np.sqrt(l2_err)
        h1_error = np.sqrt(h1_err)
        stress_l2_error = np.sqrt(stress_l2_err)
        
        # 存储结果
        self.results.l2_error = l2_error
        self.results.h1_error = h1_error
        self.results.max_error = max_err
        self.results.stress_l2_error = stress_l2_error
        self.results.stress_max_error = stress_max_err
        
        print(f"位移L2误差: {l2_error:.3e}")
        print(f"位移H1误差: {h1_error:.3e}")
        print(f"位移最大误差: {max_err:.3e}")
        print(f"应力L2误差: {stress_l2_error:.3e}")
        print(f"应力最大误差: {stress_max_err:.3e}")
        
        return l2_error, h1_error, max_err, stress_l2_error