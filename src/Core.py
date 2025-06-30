from scipy.sparse import lil_matrix, csr_matrix, dok_matrix
import numpy as np
from Model import FEMResults, ElementType
from Tools import FEMTools

class FEMCore:
    def __init__(self, case, mesh):
        """初始化有限元核心计算
        
        参数:
            case: CaseSetup对象，包含材料和几何参数
            mesh: (nodes, elements)元组，网格数据
        """
        # 首先验证形函数导数计算是否正确
        from Tools import FEMTools
        if not FEMTools.test_shape_function_derivatives():
            raise ValueError("形函数导数测试失败，请检查形函数实现")
            
        self.case = case
        self.nodes, self.elements = mesh
        self.results = FEMResults()
        self.K = None  # 全局刚度矩阵
        self.F = None  # 载荷向量
        self.debug = getattr(case, 'debug_mode', False)  # 获取调试模式设置
        
    def assemble_stiffness_matrix(self):
        """组装全局刚度矩阵
        
        返回:
            K: 全局刚度矩阵(scipy.sparse.lil_matrix)
            
        异常:
            ValueError: 如果网格数据无效或组装失败
        """
        print("\nAssembling stiffness matrix...")
        
        # 详细验证网格数据
        if self.nodes is None or self.elements is None:
            raise ValueError("网格数据未初始化")
        if len(self.nodes) == 0:
            raise ValueError("节点列表为空")
        if len(self.elements) == 0:
            raise ValueError("单元列表为空")
        if self.nodes.shape[1] != 2:
            raise ValueError("节点坐标维度不正确")
            
        # 在开始计算前验证形函数导数
        print("\n=== 验证形函数导数 ===")
        test_points = [(0,0), (1,0), (0,1), (-1,0), (0,-1)]
        for xi, eta in test_points:
            N, dN_dxi = FEMTools.shape_functions(xi, eta, self.case.element_type.value)
            print(f"\n测试点 (xi={xi}, eta={eta}):")
            print("形函数N:", N)
            print("dN/dxi:", dN_dxi[0])
            print("dN/deta:", dN_dxi[1])
            print("dN/dxi之和:", np.sum(dN_dxi[0]))
            print("dN/deta之和:", np.sum(dN_dxi[1]))
            
        n_nodes = len(self.nodes)
        dof_per_node = 2
        total_dofs = dof_per_node * n_nodes
        
        # 初始化稀疏矩阵
        self.K = dok_matrix((total_dofs, total_dofs))
        
        # 获取材料矩阵(使用更高精度)
        D = self.case.get_material_matrix().astype(np.float64)
        
        # 获取高斯积分点和权重(优化数值精度)
        if self.case.element_type == ElementType.LINEAR:
            gauss_points, gauss_weights = FEMTools.gauss_points(1)
            gauss_weights = np.array(gauss_weights, dtype=np.float64)
        else:
            gauss_points, gauss_weights = FEMTools.gauss_points(2)
            gauss_weights = np.array(gauss_weights, dtype=np.float64)
        
        # 组装单元刚度矩阵
        for elem_idx, elem in enumerate(self.elements):
            # 根据单元类型获取节点坐标
            if self.case.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]  # 8节点单元
                n_points = 8
                
            ke = np.zeros((2*n_points, 2*n_points))
            
            # 打印单元节点信息
            print(f"\n单元{elem_idx}节点坐标:")
            for i, node in enumerate(elem[:8]):
                print(f"节点{i}: {self.nodes[node]}")
            
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                print(f"\n积分点(xi={xi}, eta={eta})")
                
                # 获取形函数及其导数
                N, dN_dxi = FEMTools.shape_functions(xi, eta, 
                    self.case.element_type.value)
                
                # 详细验证形函数导数
                print(f"\n单元{elem_idx}在积分点(xi={xi}, eta={eta})")
                print("形函数N:", N)
                print("dN/dxi:", dN_dxi[0])
                print("dN/deta:", dN_dxi[1])
                print("dN/dxi之和:", np.sum(dN_dxi[0]))
                print("dN/deta之和:", np.sum(dN_dxi[1]))
                
                # 验证形函数导数性质
                sum_dxi = np.sum(dN_dxi[0])
                sum_deta = np.sum(dN_dxi[1])
                if abs(sum_dxi) > 1e-10 or abs(sum_deta) > 1e-10:
                    print(f"警告: 形函数导数和不为零! (dxi={sum_dxi:.2e}, deta={sum_deta:.2e})")
                    print("形函数导数:", dN_dxi)
                    raise ValueError("形函数导数和不为零")
                
                # 计算雅可比矩阵(优化数值稳定性)
                J = np.zeros((2, 2), dtype=np.float64)
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :].astype(np.float64) * 
                                       node_coords[:, j].astype(np.float64))
                        print(f"J[{i},{j}] = sum(dN_dxi[{i},:] * node_coords[:,{j}]) = {J[i,j]}")
                
                # 验证雅可比矩阵
                detJ = np.linalg.det(J)
                if abs(detJ) <= 1e-10:
                    print(f"\n警告: 单元{elem_idx}的雅可比行列式太小({detJ:.2e})")
                    print("节点坐标:", node_coords)
                    print("形函数导数:", dN_dxi)
                    print("雅可比矩阵:", J)
                    raise ValueError(f"单元{elem_idx}的雅可比行列式太小({detJ:.2e}) - 请检查网格质量")
                    
                # 检查雅可比矩阵条件数
                cond = np.linalg.cond(J)
                if cond > 1e6:
                    print(f"警告: 单元{elem_idx}的雅可比矩阵条件数过大({cond:.2e})")
                    print("雅可比矩阵:", J)
                    
                # 检查节点顺序
                if detJ < 0:
                    print(f"警告: 单元{elem_idx}的雅可比行列式为负({detJ:.2e}) - 可能节点顺序错误")
                    print("建议检查单元节点顺序是否为逆时针方向")
                    
                invJ = np.linalg.inv(J)
                dN_dx = invJ @ dN_dxi
                
                # 验证雅可比矩阵
                if abs(detJ) < 1e-10:
                    print(f"\n单元{elem_idx}在积分点(xi={xi}, eta={eta})的雅可比矩阵:")
                    print(f"J = [[{J[0,0]:.6e}, {J[0,1]:.6e}]")
                    print(f"     [{J[1,0]:.6e}, {J[1,1]:.6e}]]")
                    print(f"det(J) = {detJ:.6e}")
                    raise ValueError(f"单元{elem_idx}的雅可比行列式太小({detJ:.2e})")
                
                # 验证逆矩阵计算
                if np.any(np.isnan(invJ)) or np.any(np.isinf(invJ)):
                    raise ValueError(f"单元{elem_idx}的雅可比矩阵不可逆")
                
                # 应变-位移矩阵B
                B = np.zeros((3, 2*n_points))
                for i in range(n_points):
                    B[0, 2*i] = dN_dx[0, i] # ε_xx
                    B[1, 2*i+1] = dN_dx[1, i] # ε_yy 
                    B[2, 2*i] = dN_dx[1, i] # γ_xy
                    B[2, 2*i+1] = dN_dx[0, i] # γ_xy
                
                # 单元刚度矩阵贡献(优化数值稳定性)
                B = B.astype(np.float64)
                ke_part = B.T @ D @ B * detJ * w
                ke += ke_part.astype(np.float64)
                
                # 检查刚度矩阵元素是否合理
                max_ke = np.max(np.abs(ke_part))
                if max_ke > 1e20 or np.any(np.isnan(ke_part)):
                    raise ValueError(f"单元{elem_idx}刚度矩阵值异常 (max={max_ke:.2e})")
            
            # 组装到全局矩阵
            for i in range(n_points):
                for j in range(n_points):
                    for di in range(2):
                        for dj in range(2):
                            row = 2*elem[i] + di
                            col = 2*elem[j] + dj
                            self.K[row, col] += ke[2*i+di, 2*j+dj]
        
        print(f"刚度矩阵组装完成, 总自由度: {total_dofs}")
        return self.K
        
    def apply_boundary_conditions(self):
        """施加边界条件和载荷
        
        返回:
            K_ff: 自由自由度的刚度矩阵(csr_matrix)
            F_f: 自由自由度的载荷向量
            free_dofs: 自由自由度索引列表
            
        异常:
            ValueError: 如果输入无效或处理失败
        """
        print("\nApplying boundary conditions...")
        
        # 详细验证输入
        if self.K is None:
            raise ValueError("刚度矩阵未初始化")
        if not hasattr(self, 'nodes') or self.nodes is None:
            raise ValueError("节点数据未初始化")
        if self.K.shape[0] != 2 * len(self.nodes):
            raise ValueError(f"刚度矩阵维度{self.K.shape}与节点数{len(self.nodes)}不匹配")
            
        n_nodes = len(self.nodes)
        
        # 找到左端节点 (x≈0)
        left_nodes = np.where(np.isclose(self.nodes[:,0], 0.0, atol=1e-6))[0]
        
        # 计算梁中间高度
        min_y = np.min(self.nodes[:,1])
        max_y = np.max(self.nodes[:,1])
        mid_height = min_y + (max_y - min_y)/2
        
        # 应用边界条件
        fixed_dofs = []
        print("\n固定节点信息:")
        for node in left_nodes:
            y_pos = self.nodes[node,1]
            if np.isclose(y_pos, mid_height, atol=1e-6):
                # 中间节点固定u和v
                fixed_dofs.extend([2*node, 2*node+1])
                print(f"节点{node} (y={y_pos:.4f}): 固定u和v")
            else:
                # 其他节点仅限制u
                fixed_dofs.append(2*node)
                print(f"节点{node} (y={y_pos:.4f}): 仅固定u")
        
        # 去除重复并排序
        fixed_dofs = sorted(set(fixed_dofs))
        print(f"\n固定自由度总数: {len(fixed_dofs)}")
        print("固定自由度列表:", fixed_dofs)
        
        # 自由自由度
        free_dofs = [i for i in range(self.K.shape[0]) if i not in fixed_dofs]
        
        # 提取子矩阵并转换为CSR格式
        self.K_ff = self.K[free_dofs, :][:, free_dofs].tocsr()
        
        # 验证边界条件处理后的矩阵
        if self.K_ff.shape[0] == 0:
            raise ValueError("处理后刚度矩阵大小为0 - 检查边界条件设置")
        if np.any(np.isnan(self.K_ff.data)) or np.any(np.isinf(self.K_ff.data)):
            raise ValueError("处理后刚度矩阵包含NaN或Inf值")
            
        # 施加载荷(优化数值精度)
        self.F = np.zeros(self.K.shape[0], dtype=np.float64)
        right_nodes = np.where(np.isclose(
            self.nodes[:,0], self.case.geometry.L, atol=1e-6, rtol=1e-6))[0]
            
        if len(right_nodes) > 0:
            total_load = np.float64(self.case.load.P) * np.float64(self.case.geometry.H)
            node_load = np.float64(total_load / len(right_nodes))
            print(f"\n载荷信息:")
            print(f"总载荷 P*H = {total_load:.4f}")
            print(f"节点数: {len(right_nodes)}")
            print(f"每个节点载荷: {node_load:.4f}")
            
            for node in right_nodes:
                self.F[2*node+1] = -node_load  # 负号表示向下
                print(f"节点{node} (y={self.nodes[node,1]:.4f}): 施加垂直载荷 {-node_load:.4f}")
            
            applied_load = -np.sum(self.F[2*np.array(right_nodes)+1])
            print(f"实际施加总载荷: {applied_load:.4f} (应为 {total_load:.4f})")
            if not np.isclose(applied_load, total_load, rtol=1e-4):
                raise ValueError(f"施加载荷与总载荷不匹配 (差值={abs(applied_load-total_load):.2e})")
        else:
            print("警告: 未找到右端节点施加载荷")
            
        self.F_f = self.F[free_dofs]
        
        print(f"边界条件处理完成: 固定{len(fixed_dofs)}个自由度, 剩余{len(free_dofs)}个自由度")
        return self.K_ff, self.F_f, free_dofs
        
    def compute_stresses(self, u):
        """计算单元应力
        
        参数:
            u: 位移解向量
            
        返回:
            elem_stresses: 单元应力列表
            nodal_stresses: 节点平均应力字典
            
        异常:
            ValueError: 如果输入无效或计算失败
        """
        # 增强输入验证
        if u is None:
            raise ValueError("位移向量不能为None")
        if len(u) != 2 * len(self.nodes):
            raise ValueError(f"位移向量长度{len(u)}与节点数{2*len(self.nodes)}不匹配")
        
        max_disp = np.max(np.abs(u))
        print(f"最大位移: {max_disp:.4e}")
        
        # 获取材料矩阵并验证
        D = self.case.get_material_matrix()
        if D is None or D.shape != (3,3):
            print("无效的材料矩阵:", D)
            raise ValueError("无效的材料矩阵")
        print("\n材料矩阵验证通过")
        
        # 显示右端节点位移
        right_nodes = np.where(np.isclose(self.nodes[:,0], 
                            self.case.geometry.L, atol=1e-6))[0]
        if len(right_nodes) > 0:
            print("\n右端节点位移:")
            for node in right_nodes:
                print(f"节点{node}: u={u[2*node]:.4e}, v={u[2*node+1]:.4e}")
        
        print("\nComputing stresses...")

        # 获取材料矩阵
        D = self.case.get_material_matrix()
        
        # 获取高斯积分点和权重
        if self.case.element_type == ElementType.LINEAR:
            gauss_points, gauss_weights = FEMTools.gauss_points(1)
        else:
            gauss_points, gauss_weights = FEMTools.gauss_points(2)
        
        # 初始化应力存储结构
        nodal_stresses = {}
        elem_stresses = []
        
        for elem in self.elements:
            # 获取单元节点坐标和位移
            if self.case.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(4)])
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]
                elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(8)])
                n_points = 8
            
            # 在积分点计算应力
            weighted_stress = np.zeros(3)
            total_weight = 0
            gp_stresses = []
            
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                
                # 获取形函数及其导数
                N, dN_dxi = FEMTools.shape_functions(xi, eta, 
                    self.case.element_type.value)
                
                # 计算雅可比矩阵
                J = np.zeros((2, 2))
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
                        
                detJ = np.linalg.det(J)
                invJ = np.linalg.inv(J)
                dN_dx = invJ @ dN_dxi
                
                # 应变-位移矩阵B
                B = np.zeros((3, 2*n_points))
                for i in range(n_points):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
                
                # 计算应力
                stress = D @ B @ elem_disp
                gp_stresses.append(stress)
                weighted_stress += D @ B @ elem_disp * detJ * w
                total_weight += detJ * w
            
            # 平均积分点应力
            elem_stress = np.mean(gp_stresses, axis=0)
            elem_stresses.append(elem_stress)
            
            # 收集节点应力
            for i in range(n_points):
                node = elem[i]
                if node not in nodal_stresses:
                    nodal_stresses[node] = []
                nodal_stresses[node].append(elem_stress)
        
        # 计算并显示最大应力
        max_sxx = max([s[0] for s in elem_stresses])
        max_syy = max([s[1] for s in elem_stresses])
        max_txy = max([s[2] for s in elem_stresses])
        
        # 计算von Mises应力
        von_mises = [np.sqrt(s[0]**2 - s[0]*s[1] + s[1]**2 + 3*s[2]**2) 
                    for s in elem_stresses]
        max_vm = max(von_mises)
        
        print("\n=== 应力结果 ===")
        print(f"最大正应力 σ_xx: {max_sxx:.4e}")
        print(f"最大正应力 σ_yy: {max_syy:.4e}")
        print(f"最大剪应力 τ_xy: {max_txy:.4e}")
        print(f"最大von Mises应力: {max_vm:.4e}")
        
        # 显示右端节点应力
        if len(right_nodes) > 0:
            print("\n右端节点应力:")
            for node in right_nodes:
                if node in nodal_stresses:
                    s_avg = np.mean(nodal_stresses[node], axis=0)
                    s_vm = np.sqrt(s_avg[0]**2 - s_avg[0]*s_avg[1] + s_avg[1]**2 + 3*s_avg[2]**2)
                    print(f"节点{node}: σ_xx={s_avg[0]:.4e}, σ_yy={s_avg[1]:.4e}, "
                         f"τ_xy={s_avg[2]:.4e}, σ_vm={s_vm:.4e}")
        
        return elem_stresses, nodal_stresses
        