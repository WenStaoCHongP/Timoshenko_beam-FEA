import numpy as np

class FEMTools:
    @staticmethod
    def gauss_points(order):
        """返回高斯积分点和权重
        
        参数:
            order: 积分阶数 
                1: 2×2高斯积分(线性单元)
                2: 3×3高斯积分(二次单元)
            
        返回:
            points: 高斯点列表
            weights: 对应权重列表
            
        注意:
            - 2×2积分用于线性单元和剪切项减缩积分
            - 3×3积分用于二次单元弯曲项完整积分
            - 权重总和应为4.0(对应[-1,1]×[-1,1]区域面积)
        """
        if order == 1:
            # 2×2高斯积分 (精度可达3阶多项式)
            a = 1.0 / np.sqrt(3.0)  # ≈0.5773502691896257
            points = [
                [-a, -a], [a, -a], 
                [a, a], [-a, a]
            ]
            weights = [1.0, 1.0, 1.0, 1.0]  # 每个点权重1.0
        else:
            # 3×3高斯积分 (精度可达5阶多项式)
            a = np.sqrt(3.0 / 5.0)  # ≈0.7745966692414834
            w1 = 5.0 / 9.0          # ≈0.5555555555555556
            w2 = 8.0 / 9.0          # ≈0.8888888888888888
            
            # 生成所有组合
            xi_vals = [-a, 0.0, a]
            wx_vals = [w1, w2, w1]
            
            points = []
            weights = []
            for i, xi in enumerate(xi_vals):
                for j, eta in enumerate(xi_vals):
                    points.append([xi, eta])
                    weights.append(wx_vals[i] * wx_vals[j])
                    
            # 验证权重总和≈4.0
            total_weight = sum(weights)
            if not np.isclose(total_weight, 4.0, atol=1e-10):
                raise ValueError(
                    f"高斯积分权重总和{total_weight}≠4.0，积分方案可能有误"
                )
                    
        return points, weights
            
    @staticmethod
    def test_shape_function_derivatives():
        """测试形函数导数是否正确"""
        print("\n=== 形函数导数测试 ===")
        
        # 测试节点处的形函数值 (δ性质)
        print("\n=== 节点值验证 ===")
        nodes = [(-1,-1), (1,-1), (1,1), (-1,1)]  # 角节点
        for i, (xi, eta) in enumerate(nodes):
            N, _ = FEMTools.shape_functions(xi, eta, 2)
            expected = np.zeros(8)
            expected[i] = 1.0  # 在节点i处应为1，其他为0
            if not np.allclose(N, expected, atol=1e-10):
                print(f"节点{i}验证失败: N={N} (应为{expected})")
                return False
        print("所有节点处形函数值验证通过")
        
        # 测试积分点
        test_points = [
            (0.0, 0.0),      # 单元中心
            (1.0, 0.0),      # 边中点
            (0.0, 1.0),      # 边中点
            (-1.0, 0.0),     # 边中点
            (0.0, -1.0),     # 边中点
            (1.0, 1.0),      # 角点
            (-1.0, -1.0),    # 角点
            (1.0, -1.0),     # 角点
            (-1.0, 1.0),     # 角点
            (0.5, 0.5),      # 内部点
            (-0.5, -0.5)     # 内部点
        ]
        
        all_tests_passed = True
        
        for xi, eta in test_points:
            print(f"\n测试点 (xi={xi}, eta={eta}):")
            N, dN_dxi = FEMTools.shape_functions(xi, eta, 2)
            
            # 打印形函数和导数
            print("形函数值 N:", N)
            print("dN/dxi (ξ方向导数):", dN_dxi[0])
            print("dN/deta (η方向导数):", dN_dxi[1])
            
            # 验证导数性质
            sum_dxi = np.sum(dN_dxi[0])
            sum_deta = np.sum(dN_dxi[1])
            print(f"dN/dxi之和: {sum_dxi:.2e} (应为0)")
            print(f"dN/deta之和: {sum_deta:.2e} (应为0)")
            
            # 验证形函数与导数的关系 (有限差分近似)
            eps = 1e-6
            if abs(xi) < 1-eps and abs(eta) < 1-eps:  # 只在内部点验证
                N_xi_eps, _ = FEMTools.shape_functions(xi+eps, eta, 2)
                N_eta_eps, _ = FEMTools.shape_functions(xi, eta+eps, 2)
                fd_dxi = (N_xi_eps - N)/eps
                fd_deta = (N_eta_eps - N)/eps
                
                print("\n有限差分验证:")
                print("解析导数 dN/dxi:", dN_dxi[0])
                print("有限差分 dN/dxi:", fd_dxi)
                print("解析导数 dN/deta:", dN_dxi[1])
                print("有限差分 dN/deta:", fd_deta)
                
                # 检查一致性
                if not np.allclose(dN_dxi[0], fd_dxi, atol=1e-4) or \
                   not np.allclose(dN_dxi[1], fd_deta, atol=1e-4):
                    print("警告: 解析导数与有限差分不匹配!")
                    all_tests_passed = False
            
            if abs(sum_dxi) > 1e-10 or abs(sum_deta) > 1e-10:
                print("警告: 形函数导数和不为零!")
                all_tests_passed = False
        
        if all_tests_passed:
            print("\n所有形函数导数测试通过!")
        else:
            print("\n警告: 部分形函数导数测试失败!")
            
        return all_tests_passed

    @staticmethod
    def shape_functions(xi, eta, order):
        """返回形函数及其导数
        
        参数:
            xi, eta: 自然坐标(ξ,η) ∈ [-1,1]×[-1,1]
            order: 单元阶数 (1=线性, 2=二次)
            
        返回:
            N: 形函数值数组 (n_nodes,)
            dN_dxi: 形函数对自然坐标的导数 (2 x n_nodes)
                - 第一行: ∂N/∂ξ 
                - 第二行: ∂N/∂η
                
        节点顺序:
            - 线性单元(4节点): 
                0: (-1,-1), 1: (1,-1), 2: (1,1), 3: (-1,1)
            - 二次单元(8节点):
                0-3: 角节点(同线性单元顺序)
                4: 底边中点(η=-1), 5: 右边中点(ξ=1)
                6: 顶边中点(η=1), 7: 左边中点(ξ=-1)
                
        异常:
            ValueError: 如果输入无效
        """
        # 验证输入
        if not isinstance(order, int) or order not in [1, 2]:
            raise ValueError(f"无效的单元阶数: {order} (必须是1或2)")
        if abs(xi) > 1.0 + 1e-10 or abs(eta) > 1.0 + 1e-10:
            raise ValueError(f"自然坐标超出范围: xi={xi}, eta={eta} (应在[-1,1]内)")
            
        if order == 1:
            # 线性四边形单元形函数 (4节点)
            N = np.array([
                0.25 * (1 - xi) * (1 - eta),  # 节点0
                0.25 * (1 + xi) * (1 - eta),  # 节点1
                0.25 * (1 + xi) * (1 + eta),  # 节点2
                0.25 * (1 - xi) * (1 + eta)   # 节点3
            ])
            
            # 形函数导数 [∂N/∂ξ; ∂N/∂η]
            dN_dxi = np.array([
                [ -0.25*(1-eta),  0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta)],  # ∂N/∂ξ
                [ -0.25*(1-xi),  -0.25*(1+xi), 0.25*(1+xi),   0.25*(1-xi) ]   # ∂N/∂η
            ])
            
        else:
            # 二次四边形单元形函数 (8节点 Serendipity单元)
            # 角节点形函数 (i=0-3)
            N = np.array([
                0.25 * (1 - xi) * (1 - eta) * (-xi - eta - 1),  # 节点0
                0.25 * (1 + xi) * (1 - eta) * (xi - eta - 1),   # 节点1
                0.25 * (1 + xi) * (1 + eta) * (xi + eta - 1),   # 节点2
                0.25 * (1 - xi) * (1 + eta) * (-xi + eta - 1),  # 节点3
                # 边中点形函数 (i=4-7)
                0.5 * (1 - xi**2) * (1 - eta),    # 节点4 (底边中点)
                0.5 * (1 + xi) * (1 - eta**2),    # 节点5 (右边中点)
                0.5 * (1 - xi**2) * (1 + eta),    # 节点6 (顶边中点)
                0.5 * (1 - xi) * (1 - eta**2)     # 节点7 (左边中点)
            ])
            
            # 形函数导数 [∂N/∂ξ; ∂N/∂η]
            dN_dxi = np.zeros((2, 8))
            
            # ξ方向导数 (∂N/∂ξ)
            # 角节点 (i=0-3)
            xi_i = [-1, 1, 1, -1]  # 各角节点的ξ坐标
            eta_i = [-1, -1, 1, 1] # 各角节点的η坐标
            for i in range(4):
                term1 = xi_i[i] * (1 + eta_i[i]*eta) * (xi_i[i]*xi + eta_i[i]*eta - 1)
                term2 = (1 + xi_i[i]*xi) * (1 + eta_i[i]*eta) * xi_i[i]
                dN_dxi[0,i] = 0.25 * (term1 + term2)
            
            # 边中点节点 (i=4-7)
            dN_dxi[0,4] = -xi * (1 - eta)       # 节点4
            dN_dxi[0,5] = 0.5 * (1 - eta**2)    # 节点5
            dN_dxi[0,6] = -xi * (1 + eta)       # 节点6
            dN_dxi[0,7] = -0.5 * (1 - eta**2)  # 节点7
            
            # η方向导数 (∂N/∂η)
            # 角节点 (i=0-3)
            for i in range(4):
                term1 = eta_i[i] * (1 + xi_i[i]*xi) * (xi_i[i]*xi + eta_i[i]*eta - 1)
                term2 = (1 + xi_i[i]*xi) * (1 + eta_i[i]*eta) * eta_i[i]
                dN_dxi[1,i] = 0.25 * (term1 + term2)
            
            # 边中点节点 (i=4-7)
            dN_dxi[1,4] = -0.5 * (1 - xi**2)    # 节点4
            dN_dxi[1,5] = -eta * (1 + xi)       # 节点5
            dN_dxi[1,6] = 0.5 * (1 - xi**2)    # 节点6
            dN_dxi[1,7] = -eta * (1 - xi)       # 节点7
            
            # 验证形函数导数性质
            sum_dxi = np.sum(dN_dxi[0])
            sum_deta = np.sum(dN_dxi[1])
            if abs(sum_dxi) > 1e-10 or abs(sum_deta) > 1e-10:
                raise ValueError(f"形函数导数和不为零 (dxi={sum_dxi}, deta={sum_deta})")
                
        return N, dN_dxi

    @staticmethod
    def map_to_physical(xi, eta, nodes, order):
        """将自然坐标映射到物理坐标
        
        参数:
            xi, eta: 自然坐标
            nodes: 单元节点坐标数组
            order: 单元阶数
            
        返回:
            (x,y): 物理坐标
            J: 雅可比矩阵
            
        异常:
            ValueError: 如果输入无效或计算失败
        """
        # 验证输入
        if nodes is None or len(nodes) == 0:
            raise ValueError("节点坐标数组为空")
        if nodes.shape[1] != 2:
            raise ValueError("节点坐标维度不正确")
            
        try:
            N, dN_dxi = FEMTools.shape_functions(xi, eta, order)
            
            # 验证形函数
            if len(N) != len(nodes):
                raise ValueError(f"形函数数{len(N)}与节点数{len(nodes)}不匹配")
                
            # 计算物理坐标
            x = np.sum(N * nodes[:,0])
            y = np.sum(N * nodes[:,1])
            
            # 计算雅可比矩阵
            J = np.zeros((2,2))
            for i in range(2):
                for j in range(2):
                    J[i,j] = np.sum(dN_dxi[i] * nodes[:,j])
                    
            # 检查雅可比矩阵
            if np.any(np.isnan(J)) or np.any(np.isinf(J)):
                raise ValueError("雅可比矩阵包含NaN或Inf")
                
            return (x,y), J
            
        except Exception as e:
            raise ValueError(f"物理坐标映射失败: {str(e)}")

    @staticmethod
    def map_to_natural(x, y, nodes, order, tol=1e-6, max_iter=10):
        """将物理坐标映射到自然坐标(迭代法)
        
        参数:
            x, y: 物理坐标
            nodes: 单元节点坐标数组
            order: 单元阶数
            tol: 容差
            max_iter: 最大迭代次数
            
        返回:
            (xi, eta): 自然坐标
            
        异常:
            ValueError: 如果输入无效或迭代不收敛
        """
        # 验证输入
        if nodes is None or len(nodes) == 0:
            raise ValueError("节点坐标数组为空")
        if not isinstance(order, int) or order not in [1, 2]:
            raise ValueError(f"无效的单元阶数: {order} (必须是1或2)")
        if max_iter <= 0:
            raise ValueError(f"最大迭代次数必须为正数: {max_iter}")
            
        # 初始猜测(单元中心)
        xi, eta = 0.0, 0.0
        last_error = float('inf')
        
        for iter in range(max_iter):
            try:
                (x_est, y_est), J = FEMTools.map_to_physical(xi, eta, nodes, order)
                dx = x - x_est
                dy = y - y_est
                current_error = np.sqrt(dx**2 + dy**2)
                
                # 打印迭代信息
                print(f"Iter {iter+1}: xi={xi:.6f}, eta={eta:.6f}, error={current_error:.3e}")
                
                # 检查收敛
                if current_error < tol:
                    # 检查是否在规范单元内
                    if abs(xi) <= 1.0 + tol and abs(eta) <= 1.0 + tol:
                        print(f"收敛于迭代 {iter+1}, 最终误差={current_error:.3e}")
                        return xi, eta
                    else:
                        raise ValueError(f"坐标({xi:.6f},{eta:.6f})不在单元内")
                
                # 检查误差是否减小
                if current_error >= last_error:
                    print(f"警告: 误差未减小 (上次={last_error:.3e}, 当前={current_error:.3e})")
                
                last_error = current_error
                
                # 牛顿迭代
                try:
                    invJ = np.linalg.inv(J)
                    dxi = invJ[0,0]*dx + invJ[0,1]*dy
                    deta = invJ[1,0]*dx + invJ[1,1]*dy
                    
                    # 限制步长
                    max_step = 0.5
                    step = np.sqrt(dxi**2 + deta**2)
                    if step > max_step:
                        scale = max_step / step
                        dxi *= scale
                        deta *= scale
                        print(f"限制步长: {step:.3f} -> {max_step:.3f}")
                        
                    xi += dxi
                    eta += deta
                    
                except np.linalg.LinAlgError as e:
                    raise ValueError(f"雅可比矩阵奇异(行列式={np.linalg.det(J):.3e}): {str(e)}")
                    
            except Exception as e:
                raise ValueError(f"迭代{iter+1}失败: {str(e)}")
            
        raise ValueError(f"在{max_iter}次迭代内未收敛, 最终误差={last_error:.3e}")
        
    @staticmethod
    def calculate_element_area(node_coords):
        """计算单元面积(适用于线性和二次单元)
        
        参数:
            node_coords: 单元节点坐标数组
            
        返回:
            单元面积
            
        异常:
            ValueError: 如果计算失败
        """
        # 使用Shoelace公式计算多边形面积
        x = node_coords[:,0]
        y = node_coords[:,1]
        return 0.5 * np.abs(np.dot(x, np.roll(y,1)) - np.dot(y, np.roll(x,1)))

    @staticmethod 
    def test_isoparametric_mapping():
        """测试等参元映射是否正确"""
        print("\n=== 等参元映射测试 ===")
        
        # 测试线性单元
        print("\n测试线性单元...")
        nodes_linear = np.array([
            [0, 0],  # 节点0: (ξ,η)=(-1,-1)
            [2, 0],  # 节点1: (ξ,η)=(1,-1)
            [2, 1],  # 节点2: (ξ,η)=(1,1)
            [0, 1]   # 节点3: (ξ,η)=(-1,1)
        ])
        
        # 测试二次单元
        print("\n测试二次单元...")
        nodes_quad = np.array([
            [0, 0],    # 节点0: (ξ,η)=(-1,-1)
            [2, 0],    # 节点1: (ξ,η)=(1,-1)
            [2, 1],    # 节点2: (ξ,η)=(1,1)
            [0, 1],    # 节点3: (ξ,η)=(-1,1)
            [1, 0],    # 节点4: (ξ,η)=(0,-1)
            [2, 0.5],  # 节点5: (ξ,η)=(1,0)
            [1, 1],    # 节点6: (ξ,η)=(0,1)
            [0, 0.5]   # 节点7: (ξ,η)=(-1,0)
        ])
        
        test_cases = [
            (1, nodes_linear, [
                (-1, -1, 0, 0),
                (1, -1, 2, 0),
                (1, 1, 2, 1),
                (-1, 1, 0, 1),
                (0, 0, 1, 0.5)
            ]),
            (2, nodes_quad, [
                (-1, -1, 0, 0),
                (1, -1, 2, 0),
                (1, 1, 2, 1),
                (-1, 1, 0, 1),
                (0, -1, 1, 0),
                (1, 0, 2, 0.5),
                (0, 1, 1, 1),
                (-1, 0, 0, 0.5),
                (0, 0, 1.25, 0.625)
            ])
        ]
        
        for order, nodes, cases in test_cases:
            print(f"\n测试{order}阶单元...")
            for xi, eta, x_exp, y_exp in cases:
                # 正向映射测试
                (x_act, y_act), J = FEMTools.map_to_physical(xi, eta, nodes, order)
                print(f"\n自然坐标({xi},{eta}) -> 物理坐标({x_act:.2f},{y_act:.2f})")
                print(f"期望值: ({x_exp:.2f},{y_exp:.2f})")
                print(f"雅可比矩阵:\n{J}")
                
                if not np.isclose(x_act, x_exp, atol=1e-6) or not np.isclose(y_act, y_exp, atol=1e-6):
                    raise ValueError(f"正向映射错误: ({xi},{eta}) -> ({x_act},{y_act}), 应为 ({x_exp},{y_exp})")
                
                # 反向映射测试(仅对物理坐标在单元内的情况)
                if (order == 1 and 0 <= x_exp <= 2 and 0 <= y_exp <= 1) or \
                   (order == 2 and 0 <= x_exp <= 2 and 0 <= y_exp <= 1):
                    xi_back, eta_back = FEMTools.map_to_natural(x_exp, y_exp, nodes, order)
                    print(f"反向映射: ({x_exp:.2f},{y_exp:.2f}) -> ({xi_back:.2f},{eta_back:.2f})")
                    
                    if not np.isclose(xi_back, xi, atol=1e-4) or not np.isclose(eta_back, eta, atol=1e-4):
                        raise ValueError(f"反向映射错误: ({x_exp},{y_exp}) -> ({xi_back},{eta_back}), 应为 ({xi},{eta})")
        
        # 测试错误情况
        print("\n测试错误情况...")
        try:
            FEMTools.map_to_natural(3, 3, nodes_linear, 1)  # 坐标在单元外
            raise AssertionError("未捕获坐标在单元外的错误")
        except ValueError:
            print("成功捕获坐标在单元外的错误")
            
        try:
            bad_nodes = np.array([[0,0],[1,0],[1,1],[0,0]])  # 退化单元
            FEMTools.map_to_natural(0.5, 0.5, bad_nodes, 1)
            raise AssertionError("未捕获退化单元错误")
        except ValueError:
            print("成功捕获退化单元错误")
            
        print("\n所有等参元映射测试通过!")
        return True