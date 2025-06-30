from Model import ElementType
from SetCase import CaseSetup
from SetMesh import MeshGenerator
from Core import FEMCore
from Solver import LinearSolver
from Visualization import ResultVisualizer
import numpy as np
import time
import matplotlib.pyplot as plt
plt.ioff()

def run_analysis(nx=20, ny=4, element_type='quadratic', 
                length=10.0, height=1.0, E=1e9, nu=0.25, P=10.0,
                debug=False):
    """运行完整的有限元分析
    
    参数:
        nx: x方向单元数
        ny: y方向单元数
        element_type: 单元类型('linear'或'quadratic')
        length: 梁长度(m)
        height: 梁高度(m)
        E: 弹性模量(Pa)
        nu: 泊松比
        P: 分布载荷(N/m)
        
    返回:
        results: 包含分析结果的字典
    """
    # 记录开始时间
    start_time = time.time()
    
    try:
        # 1. 设置案例参数
        print("\n=== Setting up case parameters ===")
        if element_type not in ['linear', 'quadratic']:
            raise ValueError(f"无效的单元类型: {element_type} (必须是'linear'或'quadratic')")
            
        element_enum = ElementType.QUADRATIC if element_type == 'quadratic' else ElementType.LINEAR
        case = CaseSetup(length=length, height=height, elastic_modulus=E,
                        poisson_ratio=nu, distributed_load=P, element_type=element_enum,
                        debug_mode=debug)
        case.print_parameters()
        
        # 验证输入参数
        if length <= 0 or height <= 0 or E <= 0 or nu <= 0 or P <= 0:
            raise ValueError("输入参数必须为正数")
        
        # 2. 生成网格 (添加详细错误检查)
        print("\n=== Generating mesh ===")
        mesh_gen = MeshGenerator(case)
        mesh_result = mesh_gen.generate(nx=nx, ny=ny)
        if mesh_result is None or len(mesh_result) != 2:
            raise ValueError("网格生成失败: 返回无效结果")
        nodes, elements = mesh_result
        
        # 验证网格数据
        if nodes is None or elements is None:
            raise ValueError("网格数据为空")
        if len(nodes) == 0 or len(elements) == 0:
            raise ValueError("生成的网格没有节点或单元")
        if nodes.shape[1] != 2:
            raise ValueError("节点坐标维度不正确")
            
        # 3. 组装刚度矩阵 (添加详细错误检查)
        print("\n=== Assembling system ===")
        fem = FEMCore(case, (nodes, elements))
        K = fem.assemble_stiffness_matrix()
        if K is None:
            raise ValueError("刚度矩阵组装失败: 返回None")
        
        # 4. 施加边界条件 (添加详细错误检查)
        print("\n=== Applying boundary conditions ===")
        bc_result = fem.apply_boundary_conditions()
        if bc_result is None or len(bc_result) != 3:
            raise ValueError("边界条件处理失败: 返回无效结果")
        K_ff, F_f, free_dofs = bc_result
        
        # 验证边界条件结果
        if K_ff is None or F_f is None or free_dofs is None:
            raise ValueError("边界条件处理返回空值")
        print(f"F_f (自由自由度载荷向量) = {F_f}")
        print(f"K_ff (自由自由度刚度矩阵) =\n{K_ff.toarray() if hasattr(K_ff, 'toarray') else K_ff}")
        print(f"free_dofs = {free_dofs}")
        if len(free_dofs) == 0:
            print("警告: 没有自由自由度 - 检查边界条件设置")
        
        # 5. 求解线性系统 (添加详细错误检查)
        print("\n=== Solving system ===")
        solver = LinearSolver()
        u_f = solver.solve(K_ff, F_f)
        print(f"u_f (自由自由度位移解) = {u_f}")
        if u_f is None:
            raise ValueError("线性求解器返回空解")
            
        # 组装完整位移向量
        u = np.zeros(K.shape[0])
        u[free_dofs] = u_f
        
        # 6. 计算应力 (添加详细错误检查)
        print("\n=== Post-processing ===")
        stress_result = fem.compute_stresses(u)
        if stress_result is None or len(stress_result) != 2:
            raise ValueError("应力计算失败: 返回无效结果")
        elem_stresses, nodal_stresses = stress_result
        
        # 验证应力结果
        if elem_stresses is None or nodal_stresses is None:
            raise ValueError("应力计算结果为空")
        
        # 处理新的应力数据结构
        nodal_stress_array = np.zeros((len(nodes), 3))
        for node, avg_stress in nodal_stresses.items():
            avg_stress = np.asarray(avg_stress)
            if avg_stress.ndim > 1:
                avg_stress = np.mean(avg_stress, axis=0)
            nodal_stress_array[node] = avg_stress
        
        # 7. 可视化结果 (添加详细错误检查)
        print("\n=== Visualizing results ===")
        try:
            # 准备可视化数据
            print("\n=== 可视化数据准备 ===")
            print(f"节点数: {len(nodes)}, 单元数: {len(elements)}")
            print(f"位移向量长度: {len(u)}, 应力数据节点数: {len(nodal_stresses)}")
            
            # 详细检查应力数据
            print("\n应力数据详细检查:")
            print(f"应力数据类型: {type(nodal_stresses)}")
            if isinstance(nodal_stresses, dict):
                print(f"字典键数量: {len(nodal_stresses)}")
                if nodal_stresses:
                    first_key = next(iter(nodal_stresses))
                    print(f"第一个键值对示例: 节点{first_key} -> {nodal_stresses[first_key]}")
                    print(f"应力分量类型: {type(nodal_stresses[first_key])}")
                    if isinstance(nodal_stresses[first_key], dict):
                        print(f"包含的应力分量: {nodal_stresses[first_key].keys()}")
            else:
                print("警告: 应力数据不是字典类型")
            
            # 验证应力数据
            print("\n验证应力数据...")
            if not nodal_stresses:
                raise ValueError("应力数据字典为空")
                
            print(f"应力数据包含{len(nodal_stresses)}个节点的应力")
            
            # 创建应力数组并验证每个节点
            nodal_stress_array = np.zeros((len(nodes), 3))
            missing_nodes = 0
            
            for node in range(len(nodes)):
                if node not in nodal_stresses:
                    missing_nodes += 1
                    continue
                    
                stresses = nodal_stresses[node]
                if not stresses:
                    raise ValueError(f"节点{node}的应力数据为空")
                    
                # 计算平均应力并验证
                avg_stress = np.mean(stresses, axis=0)
                if avg_stress.shape != (3,):
                    raise ValueError(f"节点{node}的平均应力形状应为(3,), 实际为{avg_stress.shape}")
                    
                nodal_stress_array[node] = avg_stress
                
            if missing_nodes > 0:
                print(f"警告: {missing_nodes}个节点缺少应力数据")
                
            # 验证应力计算结果
            if np.all(nodal_stress_array == 0):
                raise ValueError("所有节点应力计算结果为零，可能计算有误")
            
            # 详细验证输入数据
            print("验证输入数据...")
            if len(nodes) == 0:
                raise ValueError("没有节点数据")
            if len(elements) == 0:
                raise ValueError("没有单元数据")
            if len(u) != 2 * len(nodes):
                raise ValueError(f"位移向量长度{len(u)}应为节点数的2倍(当前节点数:{len(nodes)})")
            if nodal_stress_array.shape != (len(nodes), 3):
                raise ValueError(f"应力数组形状应为({len(nodes)},3), 实际为{nodal_stress_array.shape}")
            
            print("创建可视化器...")
            visualizer = ResultVisualizer(nodes, elements, u, nodal_stress_array)
            
            print("生成结果图...")
            visualizer.plot_all_results(
                scale_factor=50.0,
                element_type=element_type,
                nx=nx,
                ny=ny,
                save_path=None  # 自动命名 results_{element_type}_nx{nx}_ny{ny}.png
            )
            print("可视化完成")
        except Exception as e:
            print("\n=== 可视化错误详情 ===")
            print(f"错误类型: {type(e).__name__}")
            print(f"错误信息: {str(e)}")
            print("\n当前数据状态:")
            print(f"节点数: {len(nodes) if nodes is not None else 'None'}")
            print(f"单元数: {len(elements) if elements is not None else 'None'}")
            print(f"位移向量长度: {len(u) if u is not None else 'None'}")
            if nodal_stresses is not None:
                print(f"应力数据节点数: {len(nodal_stresses)}")
                print(f"应力数组形状: {nodal_stress_array.shape if 'nodal_stress_array' in locals() else '未创建'}")
            raise ValueError("可视化失败") from e
        
        # 计算运行时间
        run_time = time.time() - start_time
        
        # 收集关键结果 (添加验证)
        print("\n收集分析结果...")
        max_stress = np.max(nodal_stress_array[:,0])
        min_stress = np.min(nodal_stress_array[:,0])
        max_disp = np.max(u[1::2])
        min_disp = np.min(u[1::2])
        
        # 验证结果值
        if not np.isfinite(max_stress) or not np.isfinite(min_stress):
            raise ValueError("无效的应力值")
        if not np.isfinite(max_disp) or not np.isfinite(min_disp):
            raise ValueError("无效的位移值")
        
        results = {
            'nodes': nodes,
            'elements': elements,
            'displacements': u,
            'stresses': nodal_stress_array,
            'max_stress': float(max_stress),  # 确保可JSON序列化
            'min_stress': float(min_stress),
            'max_displacement': float(max_disp),
            'min_displacement': float(min_disp),
            'run_time': float(run_time)
        }
        
        # 验证结果字典
        required_keys = ['nodes', 'elements', 'displacements', 'stresses', 
                        'max_stress', 'min_stress', 'max_displacement', 'min_displacement', 'run_time']
        for key in required_keys:
            if key not in results:
                raise ValueError(f"结果字典缺少关键字段: {key}")
        
        print("\n=== Analysis completed successfully ===")
        print(f"Maximum stress: {max_stress:.2f} Pa")
        print(f"Maximum displacement: {max_disp:.2e} m")
        print(f"Total run time: {run_time:.2f} seconds")
        
        
        return results
        
    except Exception as e:
        print(f"\n!!! Analysis failed: {str(e)}")
        raise

def main():
    """主程序入口"""
    print("=== Timoshenko Beam Finite Element Analysis ===")
    
    # 默认运行参数
    params = {
        'nx': 20,
        'ny': 4,
        'element_type': 'linear',  # 可选 'linear' 或 'quadratic'
        'length': 10.0,
        'height': 1.0,
        'E': 1e9,
        'nu': 0.25,
        'P': 10.0,
        'debug': False  # 设置为True以启用调试模式
    }
    
    try:
        # 运行分析
        results = run_analysis(**params)
        
        # 显示结果摘要
        print("\n=== Results Summary ===")
        print(f"Nodes: {len(results['nodes'])}")
        print(f"Elements: {len(results['elements'])}")
        print(f"Max stress (σ_xx): {results['max_stress']:.2f} Pa")
        print(f"Min stress (σ_xx): {results['min_stress']:.2f} Pa")
        print(f"Max displacement: {results['max_displacement']:.2e} m")
        print(f"Run time: {results['run_time']:.2f} seconds")
        
    except Exception as e:
        print("\nProgram terminated due to error.")
        print(f"Error details: {str(e)}")

if __name__ == "__main__":
    main()