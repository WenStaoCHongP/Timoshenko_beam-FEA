"""
Timoshenko悬臂梁有限元分析程序
作业要求：
- 梁长度(l) = 10m
- 梁高度(d) = 1m  
- 弹性模量(E) = 10^9 Pa
- 泊松比(ν) = 0.25
- 均布载荷(P) = 10 N/m (作用在右端)
- 边界条件：左端中间固定，其余位置限制水平位移
"""

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

def generate_mesh(L, H, nx, ny):
    """
    生成矩形区域的四边形网格
    参数:
        L: 长度
        H: 高度
        nx: x方向单元数
        ny: y方向单元数
    返回:
        nodes: 节点坐标数组 (n_nodes x 2)
        elements: 单元连接列表 (n_elements x 4)
    """
    # 生成节点坐标
    x = np.linspace(0, L, nx+1)
    y = np.linspace(0, H, ny+1)
    xx, yy = np.meshgrid(x, y)
    nodes = np.column_stack([xx.ravel(), yy.ravel()])
    
    # 生成单元连接
    elements = []
    for j in range(ny):
        for i in range(nx):
            n0 = i + j*(nx+1)
            n1 = n0 + 1
            n2 = n1 + (nx+1)
            n3 = n0 + (nx+1)
            elements.append([n0, n1, n2, n3])
    
    return nodes, np.array(elements)

# 材料参数
E = 1e9  # 弹性模量(Pa)
nu = 0.25  # 泊松比
P = 10  # 均布载荷(N/m)

# 几何参数
L = 10.0  # 梁长度(m)
H = 1.0   # 梁高度(m)

def main():
    print("Timoshenko Beam Finite Element Analysis - DEBUG MODE")
    
    # Debug: Print input parameters
    print("\n=== DEBUG: Input Parameters ===")
    print(f"Length (L): {L} m")
    print(f"Height (H): {H} m")
    print(f"Elastic modulus (E): {E} Pa")
    print(f"Poisson's ratio (nu): {nu}")
    print(f"Distributed load (P): {P} N/m")
    
    # 1. 网格生成
    print("\nGenerating mesh...")
    nx, ny = 20, 4  # x和y方向的单元数量
    nodes, elements = generate_mesh(L, H, nx, ny)
    
    # 2. 组装刚度矩阵
    print("Assembling stiffness matrix...")
    n_nodes = len(nodes)
    K = lil_matrix((2*n_nodes, 2*n_nodes))  # 每个节点有ux,uy两个自由度
    
    # 材料矩阵D (平面应力)
    D = E/(1-nu**2) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1-nu)/2]
    ])
    
    # 高斯积分点和权重 (4点高斯积分)
    a = 1.0/np.sqrt(3)
    gauss_points = [
        [-a, -a],  # 第一象限点
        [a, -a],   # 第二象限点
        [a, a],    # 第三象限点
        [-a, a]    # 第四象限点
    ]
    gauss_weights = [1.0, 1.0, 1.0, 1.0]  # 权重均为1
    
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
            
            # 计算雅可比矩阵
            J = dN_dxi @ node_coords
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
    
    # 3. 施加边界条件
    print("Applying boundary conditions...")
    
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
    
    # 处理边界条件
    free_dofs = [i for i in range(2*n_nodes) if i not in fixed_dofs]
    K_ff = K[free_dofs, :][:, free_dofs]
    
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
    F = np.zeros(2*n_nodes)
    
    # 找到右端节点
    right_nodes = np.where(nodes[:,0] == L)[0]
    n_right_nodes = len(right_nodes)
    
    # 计算每个节点的等效载荷 (总载荷P*L均匀分配到右端节点)
    total_load = P * H  # 总载荷(N)
    node_load = total_load / n_right_nodes
    
    # 施加垂直向下载荷(负y方向)
    for node in right_nodes:
        F[2*node+1] = -node_load
    
    # 处理载荷向量边界条件
    F_f = F[free_dofs]
    
    # 5. 求解
    print("\n=== DEBUG: Solving System Equations ===")
    print("Converting stiffness matrix to CSR format...")
    K_ff = K_ff.tocsr()  # 转换为CSR格式
    
    print("Solving system...")
    start_time = time.time()
    u_f = spsolve(K_ff, F_f)
    solve_time = time.time() - start_time
    
    # 计算残差
    residual = np.linalg.norm(K_ff.dot(u_f) - F_f)
    
    print(f"System solved in {solve_time:.3f} seconds")
    print(f"Residual norm: {residual:.2e} (should be close to 0)")
    
    # 组装完整位移向量
    u = np.zeros(2*n_nodes)
    u[free_dofs] = u_f
    
    # 6. 计算应力
    print("Computing stresses...")
    stresses = []
    
    for elem in elements:
        node_coords = nodes[elem]
        elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(4)])
        
        # 在单元中心计算应力(xi=0, eta=0)
        xi, eta = 0, 0
        dN_dxi = 0.25 * np.array([
            [-(1-eta), (1-eta), (1+eta), -(1+eta)],
            [-(1-xi), -(1+xi), (1+xi), (1-xi)]
        ])
        J = dN_dxi @ node_coords
        invJ = np.linalg.inv(J)
        dN_dx = invJ @ dN_dxi
        
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2*i] = dN_dx[0, i]
            B[1, 2*i+1] = dN_dx[1, i]
            B[2, 2*i] = dN_dx[1, i]
            B[2, 2*i+1] = dN_dx[0, i]
        
        stress = D @ B @ elem_disp
        stresses.append(stress)
    
    stresses = np.array(stresses)
    
    # 7. 误差分析和收敛性验证
    print("\n=== Error Analysis ===")
    
    # 计算L2误差和H1半范数误差
    l2_err = 0.0
    h1_err = 0.0
    
    for elem in elements:
        node_coords = nodes[elem]
        elem_disp = np.concatenate([u[2*elem[i]:2*elem[i]+2] for i in range(4)])
        
        # 使用高斯积分计算误差
        for gp, w in zip(gauss_points, gauss_weights):
            xi, eta = gp
            # 计算形函数
            N = 0.25 * np.array([
                (1-xi)*(1-eta),
                (1+xi)*(1-eta),
                (1+xi)*(1+eta),
                (1-xi)*(1+eta)
            ])
            
            # 计算物理坐标
            x_phys = N @ node_coords[:,0]
            y_phys = N @ node_coords[:,1]
            
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
            dN_dxi = 0.25 * np.array([
                [-(1-eta), (1-eta), (1+eta), -(1+eta)],
                [-(1-xi), -(1+xi), (1+xi), (1-xi)]
            ])
            J = dN_dxi @ node_coords
            detJ = np.linalg.det(J)
            invJ = np.linalg.inv(J)
            dN_dx = invJ @ dN_dxi
            
            B = np.zeros((3, 8))
            for i in range(4):
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
    scale_factor = 5.0
    deformed_nodes = nodes + scale_factor * u.reshape(-1, 2)
    
    # 创建绘图
    fig = plt.figure(figsize=(15, 5))
    
    # 子图1: 变形前后的网格
    ax1 = plt.subplot(1, 2, 1)
    for elem in elements:
        x = nodes[elem, 0]
        y = nodes[elem, 1]
        ax1.fill(x, y, edgecolor='blue', fill=False, linestyle='--')
        
        xd = deformed_nodes[elem, 0]
        yd = deformed_nodes[elem, 1]
        ax1.fill(xd, yd, edgecolor='red', fill=False)
    
    ax1.set_title('Mesh Deformation')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    ax1.legend(['Original', 'Deformed'])
    ax1.axis('equal')
    
    # 子图2: 应力云图
    ax2 = plt.subplot(1, 2, 2)
    for i, elem in enumerate(elements):
        x = nodes[elem, 0].mean()
        y = nodes[elem, 1].mean()
        sc = ax2.scatter(x, y, c=stresses[i, 0], cmap='jet', s=100)
    
    plt.colorbar(sc, ax=ax2, label='Stress (Pa)')
    ax2.set_title('Stress Distribution (x-direction)')
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.axis('equal')
    
    plt.tight_layout()
    
    # 保存图形
    output_file = 'timoshenko_beam_results.png'
    fig.savefig(output_file)
    print(f"Results visualization saved to {output_file}")
    
    # 关闭图形避免警告
    plt.close(fig)

if __name__ == "__main__":
    main()