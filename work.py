"""
Timoshenko悬臂梁有限元分析程序 - 模块化版本
"""
from enum import Enum

class ElementType(Enum):
    """有限单元类型枚举"""
    LINEAR = 1    # 线性四边形单元
    QUADRATIC = 2 # 二次四边形单元

class TimoshenkoBeamSolver:
    """
    Timoshenko悬臂梁有限元分析求解器
    
    实现功能:
    - 支持线性和二次四边形单元
    - 自动网格生成
    - 精确解计算
    - 误差分析
    - 结果可视化
    
    参数:
        length (float): 梁长度(m)，默认10.0
        height (float): 梁高度(m)，默认1.0
        elastic_modulus (float): 弹性模量(Pa)，默认1e9
        poisson_ratio (float): 泊松比，默认0.25
        distributed_load (float): 分布载荷(N/m)，默认10.0
        element_type (ElementType): 单元类型，默认ElementType.QUADRATIC
    
    示例:
        >>> solver = TimoshenkoBeamSolver()
        >>> solver.solve(nx=20, ny=4)
        >>> solver.analyze_convergence()
    """
    def __init__(self, length=10.0, height=1.0, elastic_modulus=1e9, 
                 poisson_ratio=0.25, distributed_load=10.0, 
                 element_type=ElementType.QUADRATIC):
        self.L = length
        self.H = height
        self.E = elastic_modulus
        self.nu = poisson_ratio
        self.P = distributed_load
        self.element_type = element_type
        
        # 初始化其他属性
        self.nodes = None
        self.elements = None
        self.u = None  # 位移解
        self.stresses = None  # 应力结果
        
    def print_parameters(self):
        """打印输入参数"""
        print("\n=== 输入参数 ===")
        print(f"梁长度(L): {self.L} m")
        print(f"梁高度(H): {self.H} m")
        print(f"弹性模量(E): {self.E} Pa")
        print(f"泊松比(nu): {self.nu}")
        print(f"分布载荷(P): {self.P} N/m")
        print(f"单元类型: {self.element_type.name}")

    def generate_mesh(self, nx=20, ny=4):
        """
        生成有限元网格
        
        参数:
            nx: x方向单元数
            ny: y方向单元数
        """
        print("\nGenerating mesh...")
        x = np.linspace(0, self.L, nx+1)
        y = np.linspace(0, self.H, ny+1)
        xx, yy = np.meshgrid(x, y)
        self.nodes = np.column_stack([xx.ravel(), yy.ravel()])
        
        if self.element_type == ElementType.LINEAR:
            # 线性单元(4节点)
            self.elements = []
            for j in range(ny):
                for i in range(nx):
                    n0 = i + j*(nx+1)
                    n1 = n0 + 1
                    n2 = n1 + (nx+1)
                    n3 = n0 + (nx+1)
                    self.elements.append([n0, n1, n2, n3])
        else:
            # 二次单元(8节点)
            self.elements = []
            for j in range(ny):
                for i in range(nx):
                    # 4个角节点
                    n0 = i + j*(nx+1)
                    n1 = n0 + 1
                    n2 = n1 + (nx+1)
                    n3 = n0 + (nx+1)
                    
                    # 4个边中点节点
                    n4 = len(self.nodes) + i + j*(2*nx+1)
                    n5 = len(self.nodes) + (nx+ny+1) + i + j*(nx+1)
                    n6 = n4 + 1
                    n7 = n5 - (nx+1)
                    
                    self.elements.append([n0, n1, n2, n3, n4, n5, n6, n7])
            
            # 生成边中点节点
            edge_nodes = []
            # 水平边中点(每行nx个)
            for j in range(ny+1):
                for i in range(nx):
                    x_mid = (x[i] + x[i+1])/2
                    edge_nodes.append([x_mid, y[j]])
            # 垂直边中点(每列ny个)
            for j in range(ny):
                for i in range(nx+1):
                    y_mid = (y[j] + y[j+1])/2
                    edge_nodes.append([x[i], y_mid])
            
            # 合并所有节点
            self.nodes = np.vstack([self.nodes, np.array(edge_nodes)])
            
        self.elements = np.array(self.elements)
        print(f"生成网格: {len(self.nodes)}个节点, {len(self.elements)}个单元")

    def assemble_stiffness_matrix(self):
        """组装全局刚度矩阵"""
        print("\nAssembling stiffness matrix...")
        n_nodes = len(self.nodes)
        
        # 根据单元类型确定自由度
        if self.element_type == ElementType.LINEAR:
            dof_per_node = 2
            total_dofs = dof_per_node * n_nodes
        else:
            dof_per_node = 2
            total_dofs = dof_per_node * (n_nodes + len(self.elements)*4)
            
        self.K = lil_matrix((total_dofs, total_dofs))
        
        # 材料矩阵D (平面应力)
        self.D = self.E/(1-self.nu**2) * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1-self.nu)/2]
        ])
        
        # 高斯积分点
        if self.element_type == ElementType.LINEAR:
            a = 1.0/np.sqrt(3)
            gauss_points = [[-a, -a], [a, -a], [a, a], [-a, a]]
            gauss_weights = [1.0, 1.0, 1.0, 1.0]
        else:
            # 二次单元使用3x3高斯积分
            a = np.sqrt(0.6)
            w1 = 5/9
            w2 = 8/9
            gauss_points = []
            gauss_weights = []
            for xi in [-a, 0, a]:
                for eta in [-a, 0, a]:
                    gauss_points.append([xi, eta])
                    if xi == 0 or eta == 0:
                        gauss_weights.append(w1 * w2)
                    else:
                        gauss_weights.append(w1 * w1)
        
        # 组装单元刚度矩阵
        for elem in self.elements:
            if self.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]  # 8节点单元
                n_points = 8
                
            ke = np.zeros((2*n_points, 2*n_points))
            
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                
                # 计算形函数导数
                if self.element_type == ElementType.LINEAR:
                    dN_dxi = 0.25 * np.array([
                        [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                        [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                    ])
                else:
                    # 二次单元形函数导数
                    dN_dxi = np.array([
                        [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                         0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                         -xi*(1-eta), 0.5*(1-eta**2), 
                         -xi*(1+eta), -0.5*(1-eta**2)],
                         
                        [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                         0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                         -0.5*(1-xi**2), -eta*(1+xi),
                         0.5*(1-xi**2), -eta*(1-xi)]
                    ])
                
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
                
                # 单元刚度矩阵贡献
                ke += B.T @ self.D @ B * detJ * w
            
            # 组装到全局矩阵
            for i in range(n_points):
                for j in range(n_points):
                    for di in range(2):
                        for dj in range(2):
                            row = 2*elem[i] + di
                            col = 2*elem[j] + dj
                            self.K[row, col] += ke[2*i+di, 2*j+dj]
        
        print(f"刚度矩阵组装完成, 总自由度: {total_dofs}")

    def apply_boundary_conditions(self):
        """施加边界条件和载荷"""
        print("\nApplying boundary conditions...")
        n_nodes = len(self.nodes)
        
        # 找到左端节点
        left_nodes = np.where(np.isclose(self.nodes[:,0], 0.0))[0]
        
        # 计算梁中间高度
        min_y = np.min(self.nodes[:,1])
        max_y = np.max(self.nodes[:,1])
        mid_height = min_y + (max_y - min_y)/2
        
        # 应用边界条件
        fixed_dofs = []
        for node in left_nodes:
            y_pos = self.nodes[node,1]
            if np.isclose(y_pos, mid_height, atol=1e-6):
                # 中间节点固定u和v
                fixed_dofs.extend([2*node, 2*node+1])
            else:
                # 其他节点仅限制u
                fixed_dofs.append(2*node)
        
        # 去除重复并排序
        fixed_dofs = sorted(set(fixed_dofs))
        
        # 自由自由度
        free_dofs = [i for i in range(self.K.shape[0]) if i not in fixed_dofs]
        
        # 提取子矩阵
        self.K_ff = self.K[free_dofs, :][:, free_dofs].tocsr()
        
        # 施加载荷
        self.F = np.zeros(self.K.shape[0])
        right_nodes = np.where(np.isclose(self.nodes[:,0], self.L))[0]
        node_load = (self.P * self.H) / len(right_nodes)
        
        for node in right_nodes:
            self.F[2*node+1] = -node_load  # 负号表示向下
            
        self.F_f = self.F[free_dofs]
        
        print(f"边界条件处理完成: 固定{len(fixed_dofs)}个自由度, 剩余{len(free_dofs)}个自由度")

    def solve_system(self):
        """求解线性系统"""
        print("\nSolving system...")
        start_time = time.time()
        
        try:
            u_f = spsolve(self.K_ff, self.F_f)
            solve_time = time.time() - start_time
            
            # 组装完整位移向量
            self.u = np.zeros(self.K.shape[0])
            self.u[self.free_dofs] = u_f
            
            # 计算残差
            residual = np.linalg.norm(self.K_ff.dot(u_f) - self.F_f)
            
            print(f"系统求解完成, 耗时: {solve_time:.3f}秒")
            print(f"残差范数: {residual:.2e}")
            
        except Exception as e:
            print(f"求解失败: {str(e)}")
            print("尝试最小二乘解...")
            try:
                u_f = np.linalg.lstsq(self.K_ff.toarray(), self.F_f, rcond=None)[0]
                self.u = np.zeros(self.K.shape[0])
                self.u[self.free_dofs] = u_f
                print("使用最小二乘解完成")
            except Exception as e2:
                print(f"最小二乘解失败: {str(e2)}")
                self.u = np.zeros(self.K.shape[0])
                print("使用零位移作为默认值")

    def compute_stresses(self):
        """计算单元应力"""
        print("\nComputing stresses...")
        
        # 准备高斯积分点
        if self.element_type == ElementType.LINEAR:
            a = 1.0/np.sqrt(3)
            gauss_points = [[-a, -a], [a, -a], [a, a], [-a, a]]
            gauss_weights = [1.0, 1.0, 1.0, 1.0]
        else:
            a = np.sqrt(0.6)
            w1 = 5/9
            w2 = 8/9
            gauss_points = []
            gauss_weights = []
            for xi in [-a, 0, a]:
                for eta in [-a, 0, a]:
                    gauss_points.append([xi, eta])
                    if xi == 0 or eta == 0:
                        gauss_weights.append(w1 * w2)
                    else:
                        gauss_weights.append(w1 * w1)
        
        # 存储应力结果
        self.nodal_stresses = {}
        self.elem_stresses = []
        
        for elem in self.elements:
            if self.element_type == ElementType.LINEAR:
                node_coords = self.nodes[elem]
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(4)])
                n_points = 4
            else:
                node_coords = self.nodes[elem[:8]]
                elem_disp = np.concatenate([self.u[2*elem[i]:2*elem[i]+2] for i in range(8)])
                n_points = 8
            
            # 在积分点计算应力
            gp_stresses = []
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                
                # 计算形函数导数
                if self.element_type == ElementType.LINEAR:
                    dN_dxi = 0.25 * np.array([
                        [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                        [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                    ])
                else:
                    dN_dxi = np.array([
                        [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                         0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                         -xi*(1-eta), 0.5*(1-eta**2), 
                         -xi*(1+eta), -0.5*(1-eta**2)],
                         
                        [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                         0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                         -0.5*(1-xi**2), -eta*(1+xi),
                         0.5*(1-xi**2), -eta*(1-xi)]
                    ])
                
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
                stress = self.D @ B @ elem_disp
                gp_stresses.append(stress)
            
            # 平均积分点应力
            elem_stress = np.mean(gp_stresses, axis=0)
            self.elem_stresses.append(elem_stress)
            
            # 收集节点应力
            for i in range(n_points):
                node = elem[i]
                if node not in self.nodal_stresses:
                    self.nodal_stresses[node] = []
                self.nodal_stresses[node].append(elem_stress)
        
        # 计算平均节点应力
        self.avg_nodal_stresses = {}
        for node, stresses in self.nodal_stresses.items():
            self.avg_nodal_stresses[node] = np.mean(stresses, axis=0)
        
        print("应力计算完成")

    def exact_solution(self, x, y):
        """
        计算精确解
        返回: (ux, uy, sxx, syy, sxy)
        """
        ux = (2*self.P*y)/(self.E*self.H**3) * ((6*self.L - 3*x)*x + (2 + self.nu)*y**2 - 3*self.H**2*(1 + self.nu)/2)
        uy = -(2*self.P)/(self.E*self.H**3) * (3*self.nu*y**2*(self.L - x) + (3*self.L - x)*x**2)
        sxx = (12*self.P*(self.L - x)*y) / self.H**3
        syy = 0
        sxy = -(6*self.P/self.H**3) * (self.H**2/4 - y**2)
        return ux, uy, sxx, syy, sxy

    def compute_errors(self):
        """计算数值解与精确解的误差"""
        print("\nComputing errors...")
        
        # 准备高斯积分点
        if self.element_type == ElementType.LINEAR:
            a = 1.0/np.sqrt(3)
            gauss_points = [[-a, -a], [a, -a], [a, a], [-a, a]]
            gauss_weights = [1.0, 1.0, 1.0, 1.0]
        else:
            a = np.sqrt(0.6)
            w1 = 5/9
            w2 = 8/9
            gauss_points = []
            gauss_weights = []
            for xi in [-a, 0, a]:
                for eta in [-a, 0, a]:
                    gauss_points.append([xi, eta])
                    if xi == 0 or eta == 0:
                        gauss_weights.append(w1 * w2)
                    else:
                        gauss_weights.append(w1 * w1)
        
        # 初始化误差指标
        l2_err = 0.0
        h1_err = 0.0
        max_err = 0.0
        
        for elem in self.elements:
            if self.element_type == ElementType.LINEAR:
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
                
                # 计算形函数
                if self.element_type == ElementType.LINEAR:
                    N = 0.25 * np.array([
                        (1-xi)*(1-eta),
                        (1+xi)*(1-eta),
                        (1+xi)*(1+eta),
                        (1-xi)*(1+eta)
                    ])
                else:
                    N = np.array([
                        0.25*(1-xi)*(1-eta)*(-xi-eta-1),
                        0.25*(1+xi)*(1-eta)*(xi-eta-1),
                        0.25*(1+xi)*(1+eta)*(xi+eta-1),
                        0.25*(1-xi)*(1+eta)*(-xi+eta-1),
                        0.5*(1-xi**2)*(1-eta),
                        0.5*(1+xi)*(1-eta**2),
                        0.5*(1-xi**2)*(1+eta),
                        0.5*(1-xi)*(1-eta**2)
                    ])
                
                # 计算物理坐标
                x = np.sum(N * node_coords[:,0])
                y = np.sum(N * node_coords[:,1])
                
                # 数值解
                u_num = np.zeros(2)
                for i in range(n_points):
                    u_num += N[i] * elem_disp[2*i:2*i+2]
                
                # 精确解
                ux_exact, uy_exact, _, _, _ = self.exact_solution(x, y)
                u_exact = np.array([ux_exact, uy_exact])
                
                # 计算形函数导数
                if self.element_type == ElementType.LINEAR:
                    dN_dxi = 0.25 * np.array([
                        [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                        [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                    ])
                else:
                    dN_dxi = np.array([
                        [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                         0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                         -xi*(1-eta), 0.5*(1-eta**2), 
                         -xi*(1+eta), -0.5*(1-eta**2)],
                         
                        [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                         0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                         -0.5*(1-xi**2), -eta*(1+xi),
                         0.5*(1-xi**2), -eta*(1-xi)]
                    ])
                
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
                
                strain_num = B @ elem_disp
                
                # 精确应变
                sxx_exact, _, sxy_exact = self.exact_solution(x, y)[2:]
                strain_exact = np.array([
                    sxx_exact/self.E,
                    -self.nu*sxx_exact/self.E,
                    sxy_exact/(self.E/(2*(1+self.nu)))
                ])
                
                # 累加误差
                l2_err += np.sum((u_num - u_exact)**2) * detJ * w
                h1_err += np.sum((strain_num - strain_exact)**2) * detJ * w
                max_err = max(max_err, np.max(np.abs(u_num - u_exact)))
        
        # 最终误差指标
        self.l2_error = np.sqrt(l2_err)
        self.h1_error = np.sqrt(h1_err)
        self.max_error = max_err
        
        print(f"L2误差: {self.l2_error:.3e}")
        print(f"H1半范数误差: {self.h1_error:.3e}")
        print(f"最大位移误差: {self.max_error:.3e}")

    def visualize_results(self):
        """可视化结果"""
        print("\nVisualizing results...")
        
        # 准备应力数据
        sxx_values = np.array([self.avg_nodal_stresses.get(node, [0,0,0])[0] 
                             for node in range(len(self.nodes))])
        
        # 计算误差分布
        error_values = np.zeros(len(self.nodes))
        for i, (x, y) in enumerate(self.nodes):
            ux_num = self.u[2*i]
            uy_num = self.u[2*i+1]
            ux_exact, uy_exact, _, _, _ = self.exact_solution(x, y)
            error_values[i] = np.sqrt((ux_num-ux_exact)**2 + (uy_num-uy_exact)**2)
        
        # 创建绘图
        fig = plt.figure(figsize=(18, 12))
        
        # 子图1: 变形网格
        ax1 = plt.subplot(2, 2, 1)
        scale_factor = 50.0
        deformed_nodes = self.nodes + scale_factor * self.u.reshape(-1, 2)
        
        for elem in self.elements:
            # 原始网格
            x = self.nodes[elem, 0]
            y = self.nodes[elem, 1]
            ax1.fill(x, y, edgecolor='gray', fill=False, linestyle=':', linewidth=0.8)
            
            # 变形后网格
            xd = deformed_nodes[elem, 0]
            yd = deformed_nodes[elem, 1]
            ax1.fill(xd, yd, edgecolor='red', fill=False, linewidth=1.5)
        
        ax1.set_title(f'Mesh Deformation (Magnified x{scale_factor})')
        ax1.set_xlabel('x (m)')
        ax1.set_ylabel('y (m)')
        ax1.axis('equal')
        
        # 子图2: 应力云图
        ax2 = plt.subplot(2, 2, 2)
        levels = 20
        tcf = ax2.tricontourf(self.nodes[:,0], self.nodes[:,1], sxx_values, 
                             levels=levels, cmap='jet')
        plt.colorbar(tcf, ax=ax2).set_label(r'Stress $\sigma_{xx}$ (Pa)')
        ax2.set_title(r'Stress $\sigma_{xx}$ Distribution')
        ax2.set_xlabel('x (m)')
        ax2.set_ylabel('y (m)')
        ax2.axis('equal')
        
        # 子图3: 误差分布
        ax3 = plt.subplot(2, 2, 3)
        tcf_err = ax3.tricontourf(self.nodes[:,0], self.nodes[:,1], error_values, 
                                 levels=20, cmap='hot_r')
        plt.colorbar(tcf_err, ax=ax3).set_label('Displacement Error (m)')
        ax3.set_title('Displacement Error Distribution')
        ax3.set_xlabel('x (m)')
        ax3.set_ylabel('y (m)')
        ax3.axis('equal')
        
        # 子图4: 结果摘要
        ax4 = plt.subplot(2, 2, 4)
        ax4.text(0.1, 0.8, f"L2 Error: {self.l2_error:.3e}", fontsize=10)
        ax4.text(0.1, 0.7, f"H1 Error: {self.h1_error:.3e}", fontsize=10)
        ax4.text(0.1, 0.6, f"Max Error: {self.max_error:.3e}", fontsize=10)
        ax4.text(0.1, 0.5, f"Max Stress: {np.max(sxx_values):.2f} Pa", fontsize=10)
        ax4.text(0.1, 0.4, f"Min Stress: {np.min(sxx_values):.2f} Pa", fontsize=10)
        ax4.axis('off')
        ax4.set_title('Results Summary')
        
        plt.tight_layout()
        
        # 保存图形
        output_file = 'timoshenko_beam_results.png'
        fig.savefig(output_file, dpi=300)
        plt.close(fig)
        print(f"结果已保存到 {output_file}")

    def analyze_convergence(self, nx_list=[10, 20, 40], ny_list=[2, 4, 8]):
        """
        运行收敛性分析
        
        参数:
            nx_list: x方向单元数列表
            ny_list: y方向单元数列表
        """
        print("\n=== Running Convergence Analysis ===")
        
        results = []
        for nx, ny in zip(nx_list, ny_list):
            print(f"\nRunning case: nx={nx}, ny={ny}")
            self.solve(nx, ny)
            results.append({
                'nx': nx,
                'ny': ny,
                'h': self.L/nx,
                'l2_error': self.l2_error,
                'h1_error': self.h1_error,
                'max_error': self.max_error
            })
        
        # 绘制收敛曲线
        plt.figure(figsize=(12, 5))
        
        # L2误差收敛
        plt.subplot(1, 2, 1)
        h = [r['h'] for r in results]
        l2_err = [r['l2_error'] for r in results]
        plt.loglog(h, l2_err, 'o-', label='L2 Error')
        
        # 计算收敛速率
        rate = np.polyfit(np.log(h), np.log(l2_err), 1)[0]
        plt.title(f'L2 Error Convergence (rate={rate:.2f})')
        plt.xlabel('Mesh size h')
        plt.ylabel('Error')
        plt.grid(True, which="both", ls="-")
        plt.legend()
        
        # H1误差收敛
        plt.subplot(1, 2, 2)
        h1_err = [r['h1_error'] for r in results]
        plt.loglog(h, h1_err, 'o-', label='H1 Error')
        
        # 计算收敛速率
        rate = np.polyfit(np.log(h), np.log(h1_err), 1)[0]
        plt.title(f'H1 Error Convergence (rate={rate:.2f})')
        plt.xlabel('Mesh size h')
        plt.ylabel('Error')
        plt.grid(True, which="both", ls="-")
        plt.legend()
        
        plt.tight_layout()
        plt.savefig('convergence_analysis.png')
        plt.close()
        print("\n收敛性分析完成，结果已保存到 convergence_analysis.png")

    def solve(self, nx=20, ny=4, verbose=True):
        """
        运行有限元分析
        
        参数:
            nx: x方向单元数
            ny: y方向单元数
            verbose: 是否打印详细信息
        """
        if verbose:
            self.print_parameters()
            print(f"\nRunning simulation with nx={nx}, ny={ny}")
        
        self.generate_mesh(nx, ny)
        self.assemble_stiffness_matrix()
        self.apply_boundary_conditions()
        self.solve_system()
        self.compute_stresses()
        self.compute_errors()
        
        if verbose:
            self.visualize_results()
            print("\n分析完成")
            print(f"关键结果:")
            print(f"- 最大应力: {np.max([s[0] for s in self.avg_nodal_stresses.values()]):.2f} Pa")
            print(f"- 最小应力: {np.min([s[0] for s in self.avg_nodal_stresses.values()]):.2f} Pa")
            print(f"- L2误差: {self.l2_error:.3e}")
            print(f"- H1误差: {self.h1_error:.3e}")

if __name__ == "__main__":
    """
    Timoshenko梁有限元分析程序 - 主程序
    
    使用示例:
    1. 基本分析:
       solver = TimoshenkoBeamSolver()
       solver.solve(nx=20, ny=4)
    
    2. 使用线性单元:
       solver = TimoshenkoBeamSolver(element_type=ElementType.LINEAR)
       solver.solve(nx=40, ny=8)
    
    3. 收敛性分析:
       solver = TimoshenkoBeamSolver()
       solver.analyze_convergence(nx_list=[10, 20, 40], ny_list=[2, 4, 8])
    
    4. 自定义参数:
       solver = TimoshenkoBeamSolver(
           length=8.0, 
           height=0.8,
           elastic_modulus=2e9,
           poisson_ratio=0.3,
           distributed_load=15.0
       )
       solver.solve(nx=30, ny=6)
    """
    
    # 默认运行示例
    print("=== Timoshenko梁有限元分析 ===")
    print("1. 运行标准分析 (二次单元)")
    solver = TimoshenkoBeamSolver(element_type=ElementType.QUADRATIC)
    solver.solve(nx=20, ny=4)
    
    print("\n2. 运行收敛性分析")
    solver.analyze_convergence(nx_list=[10, 20, 40], ny_list=[2, 4, 8])
    
    print("\n分析全部完成!")

if __name__ == "__main__":
    # 创建求解器实例
    solver = TimoshenkoBeamSolver(element_type=ElementType.QUADRATIC)
    
    # 运行分析
    solver.solve(nx=20, ny=4)

import numpy as np
import time
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def exact_displacement(x, y, L, P, E, d, v):
    """
    Timoshenko梁的精确位移解
    参数:
        x, y: 位置坐标
        L: 梁长度
        P: 集中载荷
        E: 弹性模量
        d: 梁高度
        v: 泊松比
    返回:
        ux, uy: x和y方向的位移
    """
    ux = (2*P*y)/(E*d**3) * ((6*L - 3*x)*x + (2 + v)*y**2 - 3*d**2*(1 + v)/2)
    uy = -(2*P)/(E*d**3) * (3*v*y**2*(L - x) + (3*L - x)*x**2)
    return ux, uy

def exact_stress(x, y, L, P, d):
    """
    Timoshenko梁的精确应力解
    参数:
        x, y: 位置坐标
        L: 梁长度
        P: 集中载荷
        d: 梁高度
    返回:
        sxx, syy, sxy: 应力分量
    """
    sxx = (12*P*(L - x)*y) / d**3
    syy = 0
    sxy = -(6*P/d**3) * (d**2/4 - y**2)
    return sxx, syy, sxy

# 精确解函数
def exact_solution(x, L, P, E, I, G, A, k):
    """
    Timoshenko梁的精确解
    参数:
        x: 沿梁长度的位置坐标
        L: 梁长度
        P: 集中载荷
        E: 弹性模量
        I: 截面惯性矩
        G: 剪切模量
        A: 截面积
        k: 剪切修正系数
    返回:
        w: 挠度
        theta: 转角
    """
    # 这里需要根据具体问题设置精确解
    # 示例: 简支梁受集中载荷
    w = (P * x * (L**2 - x**2)) / (6 * E * I)
    theta = (P * (L**2 - 3 * x**2)) / (6 * E * I)
    return w, theta

def generate_mesh(L, H, nx, ny, element_order=1):
    """
    生成矩形区域的四边形网格
    参数:
        L: 长度
        H: 高度
        nx: x方向单元数
        ny: y方向单元数
        element_order: 单元阶次 (1=线性, 2=二次)
    返回:
        nodes: 节点坐标数组 (n_nodes x 2)
        elements: 单元连接列表 (n_elements x 4或8)
    """
    # 生成节点坐标
    x = np.linspace(0, L, nx+1)
    y = np.linspace(0, H, ny+1)
    xx, yy = np.meshgrid(x, y)
    nodes = np.column_stack([xx.ravel(), yy.ravel()])
    
    if element_order == 1:
        # 线性单元(4节点)
        elements = []
        for j in range(ny):
            for i in range(nx):
                n0 = i + j*(nx+1)
                n1 = n0 + 1
                n2 = n1 + (nx+1)
                n3 = n0 + (nx+1)
                elements.append([n0, n1, n2, n3])
    else:
        # 二次单元(8节点)
        # 首先生成角节点(与线性单元相同)
        elements = []
        for j in range(ny):
            for i in range(nx):
                # 4个角节点
                n0 = i + j*(nx+1)
                n1 = n0 + 1
                n2 = n1 + (nx+1)
                n3 = n0 + (nx+1)
                
                # 4个边中点节点
                # 水平边中点(右)
                n4 = len(nodes) + i + j*(2*nx+1)
                # 垂直边中点(上)
                n5 = len(nodes) + (nx+ny+1) + i + j*(nx+1)
                # 水平边中点(左)
                n6 = n4 + 1
                # 垂直边中点(下)
                n7 = n5 - (nx+1)
                
                elements.append([n0, n1, n2, n3, n4, n5, n6, n7])
        
        # 生成边中点节点
        edge_nodes = []
        # 水平边中点(每行nx个)
        for j in range(ny+1):
            for i in range(nx):
                x_mid = (x[i] + x[i+1])/2
                edge_nodes.append([x_mid, y[j]])
        # 垂直边中点(每列ny个)
        for j in range(ny):
            for i in range(nx+1):
                y_mid = (y[j] + y[j+1])/2
                edge_nodes.append([x[i], y_mid])
        
        # 合并所有节点
        nodes = np.vstack([nodes, np.array(edge_nodes)])
        
        # 验证单元节点数
        for elem in elements:
            if len(elem) != 8:
                raise ValueError(f"二次单元应有8个节点，实际得到{len(elem)}个节点")
    
    return nodes, np.array(elements)

# 材料参数
E = 1e9  # 弹性模量(Pa)
nu = 0.25  # 泊松比
P = 10  # 均布载荷(N/m)

# 几何参数
L = 10.0  # 梁长度(m)
H = 1.0   # 梁高度(m)

def run_simulation(nx=20, ny=4, element_order=1):
    """运行单次模拟并返回误差指标
    Args:
        nx: x方向单元数
        ny: y方向单元数
        element_order: 单元阶次 (1=线性, 2=二次)
    Returns:
        dict: 包含网格参数和误差指标的字典
    """
    print(f"\nRunning simulation with nx={nx}, ny={ny}, element_order={element_order}")
    print("Timoshenko Beam Finite Element Analysis - DEBUG MODE")
    
    # 输入参数
    L = 10.0    # 梁长度(m)
    H = 1.0     # 梁高度(m)
    E = 1e9     # 弹性模量(Pa)
    nu = 0.25   # 泊松比
    P = 10      # 分布载荷(N/m)
    
    # Debug: Print input parameters
    print("\n=== DEBUG: Input Parameters ===")
    print(f"Length (L): {L} m")
    print(f"Height (H): {H} m")
    print(f"Elastic modulus (E): {E} Pa")
    print(f"Poisson's ratio (nu): {nu}")
    print(f"Distributed load (P): {P} N/m")
    
    # 1. 网格生成
    print("\nGenerating mesh...")
    nodes, elements = generate_mesh(L, H, nx, ny, element_order)
    
    # 2. 组装刚度矩阵
    print("Assembling stiffness matrix...")
    n_nodes = len(nodes)
    if element_order == 1:
        dof_per_node = 2
        K = lil_matrix((dof_per_node*n_nodes, dof_per_node*n_nodes))
    else:
        # 二次单元有额外的边中点节点
        dof_per_node = 2
        K = lil_matrix((dof_per_node*(n_nodes + len(elements)*4), dof_per_node*(n_nodes + len(elements)*4)))
    
    # 材料矩阵D (平面应力)
    D = E/(1-nu**2) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1-nu)/2]
    ])
    
    if element_order == 1:
        # 线性单元(4节点)的高斯积分
        a = 1.0/np.sqrt(3)
        gauss_points = [
            [-a, -a],  # 第一象限点
            [a, -a],   # 第二象限点
            [a, a],    # 第三象限点
            [-a, a]    # 第四象限点
        ]
        gauss_weights = [1.0, 1.0, 1.0, 1.0]
        
        for elem in elements:
            ke = np.zeros((8, 8))  # 4节点单元，每个节点2个自由度
            node_coords = nodes[elem]
            
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                # 形函数在自然坐标下的导数
                dN_dxi = 0.25 * np.array([
                    [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                    [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                ])
                
                # 详细调试信息
                print("\n=== 线性单元计算调试 ===")
                print(f"单元节点编号: {elem}")
                print(f"节点坐标:\n{node_coords}")
                print(f"节点坐标形状: {node_coords.shape}")
                print(f"形函数导数矩阵:\n{dN_dxi}")
                print(f"形函数导数矩阵形状: {dN_dxi.shape}")
                
                # 确保节点坐标是4x2矩阵
                if node_coords.shape != (4, 2):
                    node_coords = node_coords.reshape(4, 2)
                    print(f"重整后的节点坐标形状: {node_coords.shape}")
                
                # 计算雅可比矩阵
                try:
                    J = np.zeros((2, 2))
                    for i in range(2):
                        for j in range(2):
                            J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
                    print(f"计算得到的雅可比矩阵:\n{J}")
                except Exception as e:
                    print("雅可比矩阵计算错误！")
                    print(f"错误详情: {str(e)}")
                    raise
                detJ = np.linalg.det(J)
                invJ = np.linalg.inv(J)
                
                # 形函数对物理坐标的导数
                dN_dx = invJ @ dN_dxi
                
                # 应变-位移矩阵B
                B = np.zeros((3, 8))
                for i in range(4):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
                
                # 单元刚度矩阵贡献
                ke += B.T @ D @ B * detJ * w
            
            # 将单元刚度矩阵组装到全局矩阵
            for i in range(4):
                for j in range(4):
                    for di in range(2):
                        for dj in range(2):
                            row = 2*elem[i] + di
                            col = 2*elem[j] + dj
                            K[row, col] += ke[2*i+di, 2*j+dj]
    else:
        # 二次单元(8节点)的高斯积分
        # 使用3x3高斯积分
        gauss_points = []
        gauss_weights = []
        a = np.sqrt(0.6)
        w1 = 5/9
        w2 = 8/9
        for xi in [-a, 0, a]:
            for eta in [-a, 0, a]:
                gauss_points.append([xi, eta])
                if xi == 0 or eta == 0:
                    gauss_weights.append(w1 * w2)
                else:
                    gauss_weights.append(w1 * w1)
        
        # 二次单元(8节点)刚度矩阵组装
        for elem in elements:
            ke = np.zeros((16, 16))  # 8节点单元，每个节点2个自由度
            node_coords = nodes[elem[:8]]  # 取前8个节点(角节点+边中点)
            
            for gp, w in zip(gauss_points, gauss_weights):
                xi, eta = gp
                # 8节点二次形函数在自然坐标下的导数
                # 形函数公式
                N = np.array([
                    0.25*(1-xi)*(1-eta)*(-xi-eta-1),  # N1
                    0.25*(1+xi)*(1-eta)*(xi-eta-1),   # N2
                    0.25*(1+xi)*(1+eta)*(xi+eta-1),   # N3
                    0.25*(1-xi)*(1+eta)*(-xi+eta-1),  # N4
                    0.5*(1-xi**2)*(1-eta),            # N5
                    0.5*(1+xi)*(1-eta**2),            # N6
                    0.5*(1-xi**2)*(1+eta),            # N7
                    0.5*(1-xi)*(1-eta**2)             # N8
                ])
                
                # 形函数导数 (确保形状为2x8)
                dN_dxi = np.array([
                    [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                     0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                     -xi*(1-eta), 0.5*(1-eta**2), 
                     -xi*(1+eta), -0.5*(1-eta**2)],  # dN1-8/dxi
                     
                    [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                     0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                     -0.5*(1-xi**2), -eta*(1+xi),
                     0.5*(1-xi**2), -eta*(1-xi)]     # dN1-8/deta
                ])
                
                # 详细调试信息
                print("\n=== 二次单元计算调试 ===")
                print(f"单元节点数: {len(elem)} (应为8)")
                print(f"节点坐标矩阵形状: {node_coords.shape} (应为8x2)")
                print(f"形函数导数矩阵形状: {dN_dxi.shape} (应为2x8)")
                
                # 确保节点坐标是8x2矩阵
                if node_coords.shape != (8, 2):
                    try:
                        node_coords = node_coords.reshape(8, 2)
                        print(f"重整后的节点坐标形状: {node_coords.shape}")
                    except ValueError:
                        raise ValueError(f"无法将形状{node_coords.shape}重整为8x2矩阵")
                
                # 确保形函数导数是2x8矩阵
                if dN_dxi.shape != (2, 8):
                    raise ValueError(f"形函数导数矩阵形状应为2x8，实际为{dN_dxi.shape}")
                
                # 安全的雅可比矩阵计算
                try:
                    J = np.zeros((2, 2))
                    for i in range(2):
                        for j in range(2):
                            J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
                    print(f"计算得到的雅可比矩阵:\n{J}")
                except Exception as e:
                    print("雅可比矩阵计算错误！")
                    print(f"节点坐标:\n{node_coords}")
                    print(f"形函数导数:\n{dN_dxi}")
                    raise
                detJ = np.linalg.det(J)
                invJ = np.linalg.inv(J)
                
                # 形函数对物理坐标的导数
                print(f"\n=== 矩阵维度调试 ===")
                print(f"invJ形状: {invJ.shape} (应为2x2)")
                print(f"dN_dxi形状: {dN_dxi.shape} (应为2x8)")
                
                # 正确的矩阵乘法: (2x2) @ (2x8) = (2x8)
                dN_dx = invJ @ dN_dxi
                
                print(f"dN_dx形状: {dN_dx.shape} (应为2x8)")
                
                # 应变-位移矩阵B
                B = np.zeros((3, 16))
                for i in range(8):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
                
                # 单元刚度矩阵贡献
                ke += B.T @ D @ B * detJ * w
            
            # 将单元刚度矩阵组装到全局矩阵
            for i in range(8):
                for j in range(8):
                    for di in range(2):
                        for dj in range(2):
                            row = 2*elem[i] + di
                            col = 2*elem[j] + dj
                            K[row, col] += ke[2*i+di, 2*j+dj]
    
    # 3. 施加边界条件
    print("Applying boundary conditions...")
    
    if element_order == 1:
        # 线性单元边界条件
        # 找到左端中间节点
        left_nodes = np.where(nodes[:,0] == 0)[0]
        mid_height = H/2
        mid_node = left_nodes[np.argmin(np.abs(nodes[left_nodes,1] - mid_height))]
        
        # 固定中间节点(ux=uy=0)
        fixed_dofs = [2*mid_node, 2*mid_node+1]
        
        # 其他左端节点限制水平位移(ux=0)
        for node in left_nodes:
            if node != mid_node:
                fixed_dofs.append(2*node)
    else:
        # 二次单元边界条件
        # 找到左端所有节点(包括边中点)
        left_nodes = np.where(nodes[:,0] == 0)[0]
        mid_height = H/2
        mid_node = left_nodes[np.argmin(np.abs(nodes[left_nodes,1] - mid_height))]
        
        # 固定中间节点(ux=uy=0)
        fixed_dofs = [2*mid_node, 2*mid_node+1]
        
        # 其他左端节点限制水平位移(ux=0)
        for node in left_nodes:
            if node != mid_node:
                fixed_dofs.append(2*node)
        
        # 对于边中点节点，根据其位置施加约束
        # 顶部和底部的边中点节点只约束ux
        # 中间的边中点节点约束ux和uy
        for node in left_nodes:
            y_pos = nodes[node,1]
            if y_pos == 0 or y_pos == H:  # 底部或顶部边中点
                fixed_dofs.append(2*node)
            elif y_pos == mid_height:  # 中间边中点
                fixed_dofs.extend([2*node, 2*node+1])
    
    # 重新设计边界条件处理
    print(f"\n=== 节点和自由度详细调试 ===")
    print(f"总节点数: {n_nodes}")
    print(f"总自由度: {K.shape[0]} (应为2*{n_nodes}={2*n_nodes})")
    
    # 获取左端所有节点(更可靠的方法)
    left_nodes = []
    for i in range(n_nodes):
        if np.isclose(nodes[i, 0], 0.0):  # x坐标接近0的节点
            left_nodes.append(i)
    
    print(f"左端节点数: {len(left_nodes)}")
    print(f"左端节点索引: {left_nodes}")
    
    # 计算梁高度和中间高度
    min_y = np.min(nodes[:,1])
    max_y = np.max(nodes[:,1])
    H = max_y - min_y
    mid_height = min_y + H/2
    
    # 应用新的边界条件
    fixed_dofs = []
    for node in left_nodes:
        y_pos = nodes[node,1]
        if np.isclose(y_pos, mid_height, atol=1e-6):  # 中间点固定u和v
            fixed_dofs.extend([2*node, 2*node+1])
            print(f"节点{node} (y={y_pos:.4f}) 完全固定")
        else:  # 其他左端点仅限制u
            fixed_dofs.append(2*node)
            print(f"节点{node} (y={y_pos:.4f}) 限制水平位移")
    
    # 去除重复并验证
    fixed_dofs = sorted(set(fixed_dofs))
    print(f"固定自由度索引: {fixed_dofs}")
    
    # 验证自由度范围
    if max(fixed_dofs) >= 2*n_nodes:
        raise ValueError(f"固定自由度索引{max(fixed_dofs)}超出总自由度{2*n_nodes}")
    
    free_dofs = [i for i in range(2*n_nodes) if i not in fixed_dofs]
    
    print(f"有效固定自由度数: {len(fixed_dofs)}")
    print(f"自由自由度数: {len(free_dofs)}")
    
    # 验证自由度编号
    if max(fixed_dofs) >= K.shape[0]:
        raise ValueError(f"固定自由度索引{max(fixed_dofs)}超出总自由度{K.shape[0]}")
    
    # 提取自由度的子矩阵
    K_ff = K[free_dofs, :][:, free_dofs]
    
    # 验证刚度矩阵维度
    if K_ff.shape[0] != len(free_dofs) or K_ff.shape[1] != len(free_dofs):
        raise ValueError(f"刚度矩阵维度不匹配: K_ff形状{K_ff.shape}, 自由自由度数{len(free_dofs)}")
    
    # Debug: 验证刚度矩阵
    print("\n=== DEBUG: Stiffness Matrix Validation ===")
    # 检查对称性
    sym_diff = np.abs(K_ff - K_ff.T).max()
    print(f"Max symmetry difference: {sym_diff:.2e} (should be close to 0)")
    
    # 检查正定性
    try:
        eigvals = np.linalg.eigvalsh(K_ff.toarray())
        min_eigval = np.min(eigvals)
        print(f"Minimum eigenvalue: {min_eigval:.2e} (should be positive)")
        if min_eigval <= 0:
            print("WARNING: Stiffness matrix is not positive definite!")
    except Exception as e:
        print(f"Eigenvalue computation failed: {str(e)}")
    
    # 4. 施加载荷
    print("Applying loads...")
    
    # 创建载荷向量
    if element_order == 1:
        F = np.zeros(2*n_nodes)
    else:
        # 二次单元有额外的边中点节点
        F = np.zeros(2*(n_nodes + len(elements)*4))
    
    # 找到右端节点
    right_nodes = np.where(nodes[:,0] == L)[0]
    n_right_nodes = len(right_nodes)
    
    # 计算每个节点的等效载荷 (总载荷P*H均匀分配到右端节点)
    total_load = P * H  # 总载荷(N)
    
    if element_order == 1:
        # 线性单元载荷施加
        node_load = total_load / n_right_nodes
        for node in right_nodes:
            F[2*node+1] = -node_load  # 负号表示向下
    else:
        # 二次单元载荷施加
        # 均匀分配到所有右端节点（包括边中点节点）
        node_load = total_load / n_right_nodes
        for node in right_nodes:
            F[2*node+1] = -node_load  # 负号表示向下
    
    # 处理载荷向量边界条件
    F_f = F[free_dofs]
    
    # 5. 求解
    print("\n=== DEBUG: Solving System Equations ===")
    print("Converting stiffness matrix to CSR format...")
    K_ff = K_ff.tocsr()  # 转换为CSR格式
    
    print("Solving system...")
    start_time = time.time()
    try:
        u_f = spsolve(K_ff, F_f)
        solve_time = time.time() - start_time
        
        # 检查解的有效性
        if not np.all(np.isfinite(u_f)):
            raise ValueError("求解得到的位移包含非有限值")
            
        # 计算残差
        residual = np.linalg.norm(K_ff.dot(u_f) - F_f)
        
        print(f"System solved in {solve_time:.3f} seconds")
        print(f"Residual norm: {residual:.2e} (should be close to 0)")
        
        # 组装完整位移向量
        u = np.zeros(2*n_nodes)
        u[free_dofs] = u_f
        
    except Exception as e:
        print(f"系统求解失败: {str(e)}")
        print("尝试使用最小二乘解...")
        try:
            u_f = np.linalg.lstsq(K_ff.toarray(), F_f, rcond=None)[0]
            if not np.all(np.isfinite(u_f)):
                raise ValueError("最小二乘解也包含非有限值")
                
            u = np.zeros(2*n_nodes)
            u[free_dofs] = u_f
            print("使用最小二乘解继续计算")
        except Exception as e2:
            print(f"最小二乘解也失败: {str(e2)}")
            print("使用零位移作为默认值")
            u = np.zeros(2*n_nodes)
    
    # 6. 计算应力
    print("Computing stresses...")
    
    # 准备存储应力结果
    nodal_stresses = {}  # 节点应力 (平均来自相邻单元)
    elem_stresses = []   # 单元应力
    
    for elem in elements:
        if element_order == 1:
            node_coords = nodes[elem]
            elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(4)])
            n_points = 4
        else:
            node_coords = nodes[elem[:8]]  # 8节点单元
            elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(8)])
            n_points = 8
        
        # 在积分点计算应力
        gp_stresses = []
        for gp, w in zip(gauss_points, gauss_weights):
            xi, eta = gp
            
            if element_order == 1:
                # 线性单元形函数导数
                dN_dxi = 0.25 * np.array([
                    [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                    [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                ])
            else:
                # 二次单元形函数导数 (修正为2x8矩阵)
                dN_dxi = np.array([
                    [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                     0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                     -xi*(1-eta), 0.5*(1-eta**2), 
                     -xi*(1+eta), -0.5*(1-eta**2)],  # dN/dxi
                     
                    [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                     0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                     -0.5*(1-xi**2), -eta*(1+xi),
                     0.5*(1-xi**2), -eta*(1-xi)]    # dN/deta
                ])
            
            # 计算雅可比矩阵 (显式计算确保维度正确)
            J = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
            
            # 调试信息
            print(f"\n=== 雅可比矩阵详细调试 ===")
            print(f"dN_dxi:\n{dN_dxi}")
            print(f"node_coords:\n{node_coords}")
            print(f"计算得到的雅可比矩阵:\n{J}")
            
            # 检查雅可比矩阵是否可逆
            detJ = np.linalg.det(J)
            if abs(detJ) < 1e-10:
                raise ValueError(f"雅可比矩阵行列式太小({detJ:.2e})，可能导致数值不稳定")
                
            invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            
            B = np.zeros((3, 2*n_points))
            for i in range(n_points):
                B[0, 2*i] = dN_dx[0, i]
                B[1, 2*i+1] = dN_dx[1, i]
                B[2, 2*i] = dN_dx[1, i]
                B[2, 2*i+1] = dN_dx[0, i]
            
            # 计算应力
            try:
                stress = D @ B @ elem_disp
                if not np.all(np.isfinite(stress)):
                    raise ValueError("计算得到的应力包含非有限值")
                gp_stresses.append(stress)
            except Exception as e:
                print(f"应力计算错误: {str(e)}")
                print(f"单元: {elem}")
                print(f"位移: {elem_disp}")
                print(f"B矩阵: {B}")
                raise
        
        # 平均积分点应力作为单元应力，确保结果有效
        if len(gp_stresses) > 0:
            elem_stress = np.mean(gp_stresses, axis=0)
            if np.all(np.isfinite(elem_stress)):
                elem_stresses.append(elem_stress)
            else:
                print(f"警告: 单元{elem}的平均应力包含非有限值")
                elem_stresses.append(np.zeros(3))  # 使用零值替代
        else:
            print(f"警告: 单元{elem}没有有效的应力计算结果")
            elem_stresses.append(np.zeros(3))  # 使用零值替代
        
        # 收集节点应力 (用于后续可视化)
        for i in range(n_points):
            node = elem[i]
            if node not in nodal_stresses:
                nodal_stresses[node] = []
            nodal_stresses[node].append(elem_stress)
    
    # 计算平均节点应力
    avg_nodal_stresses = {}
    for node, stresses in nodal_stresses.items():
        avg_nodal_stresses[node] = np.mean(stresses, axis=0)
    
    # 准备应力数据用于云图
    x_coords = nodes[:,0]
    y_coords = nodes[:,1]
    sxx_values = np.array([avg_nodal_stresses.get(node, [0,0,0])[0] 
                          for node in range(len(nodes))])
    
    # 7. 误差分析和收敛性验证
    print("\n=== Error Analysis ===")
    
    # 计算L2误差和H1半范数误差
    l2_err = 0.0
    h1_err = 0.0
    
    for elem in elements:
        if element_order == 1:
            node_coords = nodes[elem]
            elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(4)])
            n_points = 4
        else:
            node_coords = nodes[elem[:8]]  # 8节点单元
            elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(8)])
            n_points = 8
        
        # 使用高斯积分计算误差
        for gp, w in zip(gauss_points, gauss_weights):
            xi, eta = gp
            
            if element_order == 1:
                # 线性单元形函数
                N = 0.25 * np.array([
                    (1-xi)*(1-eta),
                    (1+xi)*(1-eta),
                    (1+xi)*(1+eta),
                    (1-xi)*(1+eta)
                ])
            else:
                # 二次单元形函数
                N = np.array([
                    0.25*(1-xi)*(1-eta)*(-xi-eta-1),  # N1
                    0.25*(1+xi)*(1-eta)*(xi-eta-1),   # N2
                    0.25*(1+xi)*(1+eta)*(xi+eta-1),   # N3
                    0.25*(1-xi)*(1+eta)*(-xi+eta-1),  # N4
                    0.5*(1-xi**2)*(1-eta),            # N5
                    0.5*(1+xi)*(1-eta**2),            # N6
                    0.5*(1-xi**2)*(1+eta),            # N7
                    0.5*(1-xi)*(1-eta**2)             # N8
                ])
            
            # 计算物理坐标
            x_phys = np.sum(N * node_coords[:,0])
            y_phys = np.sum(N * node_coords[:,1])
            
            # 计算数值解
            u_num = np.zeros(2)
            for i in range(4):
                u_num += N[i] * elem_disp[2*i:2*i+2]
            
            # 计算精确解
            ux_exact, uy_exact = exact_displacement(x_phys, y_phys, L, P, E, H, nu)
            u_exact = np.array([ux_exact, uy_exact])
            
            # 计算精确解的应变
            sxx_exact, _, sxy_exact = exact_stress(x_phys, y_phys, L, P, H)
            strain_exact = np.array([sxx_exact/E, nu*sxx_exact/E, sxy_exact/(E/(2*(1+nu)))])
            
            # 计算数值解的应变
            if element_order == 1:
                # 线性单元形函数导数
                dN_dxi = 0.25 * np.array([
                    [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                    [-(1-xi), -(1+xi), (1+xi), (1-xi)]
                ])
                # 显式计算雅可比矩阵
                J = np.zeros((2, 2))
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
            else:
                # 二次单元形函数导数
                dN_dxi = np.array([
                    [0.25*(1-eta)*(2*xi+eta), 0.25*(1-eta)*(2*xi-eta), 
                     0.25*(1+eta)*(2*xi+eta), 0.25*(1+eta)*(2*xi-eta),
                     -xi*(1-eta), 0.5*(1-eta**2), 
                     -xi*(1+eta), -0.5*(1-eta**2)],  # dN/dxi
                     
                    [0.25*(1-xi)*(xi+2*eta), 0.25*(1+xi)*(xi-2*eta),
                     0.25*(1+xi)*(xi+2*eta), 0.25*(1-xi)*(xi-2*eta),
                     -0.5*(1-xi**2), -eta*(1+xi),
                     0.5*(1-xi**2), -eta*(1-xi)]    # dN/deta
                ])
                # 显式计算雅可比矩阵
                J = np.zeros((2, 2))
                for i in range(2):
                    for j in range(2):
                        J[i, j] = np.sum(dN_dxi[i, :] * node_coords[:, j])
            
            detJ = np.linalg.det(J)
            invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            
            # 应变-位移矩阵B
            if element_order == 1:
                B = np.zeros((3, 8))
                for i in range(4):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
            else:
                B = np.zeros((3, 16))
                for i in range(8):
                    B[0, 2*i] = dN_dx[0, i]
                    B[1, 2*i+1] = dN_dx[1, i]
                    B[2, 2*i] = dN_dx[1, i]
                    B[2, 2*i+1] = dN_dx[0, i]
            
            strain_num = B @ elem_disp
            
            # 累加L2误差
            l2_err += np.sum((u_num - u_exact)**2) * detJ * w
            
            # 累加H1半范数误差
            h1_err += np.sum((strain_num - strain_exact)**2) * detJ * w
    
    # 取平方根得到最终误差
    l2_err = np.sqrt(l2_err)
    h1_err = np.sqrt(h1_err)
    
    print(f"L2 Error: {l2_err:.3e}")
    print(f"H1 Semi-norm Error: {h1_err:.3e}")
    
    # 8. 可视化结果
    print("\nVisualizing results...")
    
    # 放大位移以便观察变形
    scale_factor = 50.0  # 增大变形放大系数
    deformed_nodes = nodes + scale_factor * u.reshape(-1, 2)
    
    # 创建绘图
    fig = plt.figure(figsize=(18, 12))
    
    # 子图1: 变形前后的网格
    ax1 = plt.subplot(2, 2, 1)
    for elem in elements:
        # 绘制原始网格
        x = nodes[elem, 0]
        y = nodes[elem, 1]
        ax1.fill(x, y, edgecolor='gray', fill=False, linestyle=':', linewidth=0.8, alpha=0.7)
        
        # 绘制变形后网格
        xd = deformed_nodes[elem, 0]
        yd = deformed_nodes[elem, 1]
        ax1.fill(xd, yd, edgecolor='red', fill=False, linewidth=1.5)
    
    # 添加位移矢量图
    ax1.quiver(nodes[:,0], nodes[:,1], 
               scale_factor*u[::2], scale_factor*u[1::2],
               color='blue', scale=1.0, scale_units='xy', angles='xy')
    
    ax1.set_title(f'Mesh Deformation (Magnified x{scale_factor})', fontsize=12)
    ax1.set_xlabel('x (m)', fontsize=10)
    ax1.set_ylabel('y (m)', fontsize=10)
    ax1.legend(['Original', 'Deformed', 'Displacement'], loc='upper right')
    ax1.grid(True, linestyle=':', alpha=0.5)
    ax1.axis('equal')
    
    # 子图2: 应力云图
    ax2 = plt.subplot(2, 2, 2)
    # 使用tricontourf创建平滑的应力云图
    levels = 20
    tcf = ax2.tricontourf(nodes[:,0], nodes[:,1], sxx_values, levels=levels, cmap='jet')
    
    # 添加应力云图轮廓线和标注
    cs = ax2.tricontour(nodes[:,0], nodes[:,1], sxx_values, levels=levels, 
                       colors='k', linewidths=0.5, linestyles='-', alpha=0.3)
    ax2.clabel(cs, inline=True, fontsize=8, fmt='%.1f')
    
    # 标记最大/最小应力点
    max_idx = np.argmax(sxx_values)
    min_idx = np.argmin(sxx_values)
    ax2.plot(nodes[max_idx,0], nodes[max_idx,1], 'w*', markersize=10)
    ax2.plot(nodes[min_idx,0], nodes[min_idx,1], 'k*', markersize=10)
    
    # 添加颜色条
    cbar = plt.colorbar(tcf, ax=ax2)
    cbar.set_label(r'Stress $\sigma_{xx}$ (Pa)', fontsize=10)
    
    ax2.set_title(r'Stress $\sigma_{xx}$ Distribution', fontsize=12)
    ax2.set_xlabel('x (m)', fontsize=10)
    ax2.set_ylabel('y (m)', fontsize=10)
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.axis('equal')
    
    # 子图3: 误差分布
    ax3 = plt.subplot(2, 2, 3)
    # 计算每个节点的误差
    error_values = np.zeros(len(nodes))
    for i, (x, y) in enumerate(nodes):
        ux_num = u[2*i]
        uy_num = u[2*i+1]
        ux_exact, uy_exact = exact_displacement(x, y, L, P, E, H, nu)
        error_values[i] = np.sqrt((ux_num-ux_exact)**2 + (uy_num-uy_exact)**2)
    
    # 绘制误差云图
    tcf_err = ax3.tricontourf(nodes[:,0], nodes[:,1], error_values, levels=20, cmap='hot_r')
    cbar_err = plt.colorbar(tcf_err, ax=ax3)
    cbar_err.set_label('Displacement Error (m)', fontsize=10)
    
    ax3.set_title('Displacement Error Distribution', fontsize=12)
    ax3.set_xlabel('x (m)', fontsize=10)
    ax3.set_ylabel('y (m)', fontsize=10)
    ax3.grid(True, linestyle=':', alpha=0.5)
    ax3.axis('equal')
    
    # 子图4: 收敛性分析
    ax4 = plt.subplot(2, 2, 4)
    # 计算理论解和数值解在右端中点的位移
    mid_node = np.argmin(np.abs(nodes[:,0] - L) + np.abs(nodes[:,1] - H/2))
    uy_num = u[2*mid_node+1]
    _, uy_exact = exact_displacement(L, H/2, L, P, E, H, nu)
    
    # 计算相对误差
    rel_error = np.abs(uy_num - uy_exact) / np.abs(uy_exact)
    
    # 显示关键结果
    ax4.text(0.1, 0.8, f"Max Stress: {np.max(sxx_values):.2f} Pa", fontsize=10)
    ax4.text(0.1, 0.7, f"Min Stress: {np.min(sxx_values):.2f} Pa", fontsize=10)
    ax4.text(0.1, 0.6, f"Tip Displacement: {uy_num:.4e} m", fontsize=10)
    ax4.text(0.1, 0.5, f"Exact Tip Displacement: {uy_exact:.4e} m", fontsize=10)
    ax4.text(0.1, 0.4, f"Relative Error: {rel_error*100:.2f}%", fontsize=10)
    
    ax4.set_title('Key Results Summary', fontsize=12)
    ax4.axis('off')
    
    plt.tight_layout(pad=3.0)
    
    # 保存图形
    output_file = 'timoshenko_beam_results.png'
    fig.savefig(output_file, dpi=300)
    print(f"Results visualization saved to {output_file}")
    
    # 关闭图形避免警告
    plt.close(fig)
    
    # 返回结果字典
    return {
        'nx': nx,
        'ny': ny,
        'element_order': element_order,
        'l2_error': l2_err,
        'h1_error': h1_err,
        'num_dofs': len(free_dofs),
        'h': L/nx  # 特征网格尺寸
    }

def main():
    """主函数，运行多组参数对比"""
    # 调试选项
    debug_quad_only = True  # 仅调试二次单元
    
    if debug_quad_only:
        print("\n=== DEBUG MODE: 仅处理二次单元 ===")
        # 仅运行二次单元测试
        order_cases = [
            (20, 4, 2)   # 二次单元
        ]
        # 跳过网格密度测试
        mesh_cases = []
    else:
        # 运行不同网格密度的模拟
        mesh_cases = [
            (20, 4),  # 基础网格
            (40, 4),  # x方向加密
            (40, 8)   # x和y方向都加密
        ]
        
        # 运行不同单元阶次的模拟
        order_cases = [
            (20, 4, 1),  # 线性单元
            (20, 4, 2)   # 二次单元
        ]
    
    # 收集所有结果
    all_results = []
    
    # 运行网格密度对比
    print("\n=== Running Mesh Density Comparison ===")
    for nx, ny in mesh_cases:
        result = run_simulation(nx, ny, element_order=1)
        all_results.append(result)
    
    # 运行单元阶次对比
    print("\n=== Running Element Order Comparison ===")
    for nx, ny, order in order_cases:
        result = run_simulation(nx, ny, element_order=order)
        all_results.append(result)
    
    # 分析并绘制收敛曲线
    analyze_results(all_results)

def analyze_results(results):
    """分析并绘制收敛曲线"""
    print("\n=== Analyzing Results ===")
    
    # 分离网格密度和单元阶次结果
    mesh_results = [r for r in results if r['element_order'] == 1]
    order_results = [r for r in results if r['nx'] == 20 and r['ny'] == 4]
    
    # 绘制网格密度收敛曲线
    plt.figure(figsize=(12, 5))
    
    # L2误差收敛
    plt.subplot(1, 2, 1)
    h = [r['h'] for r in mesh_results]
    l2_err = [r['l2_error'] for r in mesh_results]
    plt.loglog(h, l2_err, 'o-', label='L2 Error')
    
    # 计算收敛速率
    rate = np.polyfit(np.log(h), np.log(l2_err), 1)[0]
    plt.title(f'L2 Error Convergence (rate={rate:.2f})')
    plt.xlabel('Mesh size h')
    plt.ylabel('Error')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    
    # H1误差收敛
    plt.subplot(1, 2, 2)
    h1_err = [r['h1_error'] for r in mesh_results]
    plt.loglog(h, h1_err, 'o-', label='H1 Error')
    
    # 计算收敛速率
    rate = np.polyfit(np.log(h), np.log(h1_err), 1)[0]
    plt.title(f'H1 Error Convergence (rate={rate:.2f})')
    plt.xlabel('Mesh size h')
    plt.ylabel('Error')
    plt.grid(True, which="both", ls="-")
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('convergence_mesh.png')
    plt.close()
    
    # 绘制单元阶次对比
    plt.figure(figsize=(8, 5))
    orders = [r['element_order'] for r in order_results]
    l2_err = [r['l2_error'] for r in order_results]
    plt.semilogy(orders, l2_err, 'o-')
    plt.title('L2 Error vs Element Order')
    plt.xlabel('Element Order')
    plt.ylabel('Error')
    plt.xticks([1, 2], ['Linear', 'Quadratic'])
    plt.grid(True)
    plt.savefig('convergence_order.png')
    plt.close()
    
    # 打印结果表格
    print("\n=== Results Summary ===")
    print("Mesh Density Comparison:")
    print("nx\tny\th\t\tL2 Error\t\tH1 Error")
    for r in mesh_results:
        print(f"{r['nx']}\t{r['ny']}\t{r['h']:.4f}\t{r['l2_error']:.3e}\t{r['h1_error']:.3e}")
    
    print("\nElement Order Comparison:")
    print("Order\tL2 Error\t\tH1 Error")
    for r in order_results:
        print(f"{r['element_order']}\t{r['l2_error']:.3e}\t{r['h1_error']:.3e}")

if __name__ == "__main__":
    main()