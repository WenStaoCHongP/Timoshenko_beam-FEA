import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import matplotlib.tri as mtri

class ResultVisualizer:
    def __init__(self, nodes, elements, displacements=None, stresses=None):
        """初始化可视化器
        
        参数:
            nodes: 节点坐标数组
            elements: 单元连接列表
            displacements: 节点位移数组
            stresses: 节点应力数据数组
        """
        print("\n=== 可视化器初始化 ===")
        print(f"节点数组形状: {nodes.shape if nodes is not None else 'None'}")
        print(f"单元列表长度: {len(elements) if elements is not None else 'None'}")
        print(f"位移数组长度: {len(displacements) if displacements is not None else 'None'}")
        print(f"应力数据类型: {type(stresses)}")
        
        # 详细验证输入数据
        if nodes is None or elements is None:
            raise ValueError("节点或单元数据不能为None")
        if displacements is None:
            raise ValueError("位移数据不能为None")

        self.nodes = nodes
        self.elements = elements
        self.displacements = displacements
        self.stresses = stresses
        
    def plot_deformation(self, scale_factor=30000.0, save_path=None, show_node_numbers=False):
        """绘制变形前后的网格
        
        参数:
            scale_factor: 变形放大系数(默认30000)
            save_path: 图片保存路径(可选)
            show_node_numbers: 是否显示节点编号(默认False)
            
        异常:
            ValueError: 如果输入数据无效或参数不合理
        """
        # 验证输入数据
        if self.displacements is None:
            raise ValueError("没有位移数据可供可视化")
        if self.nodes is None or self.elements is None:
            raise ValueError("缺少节点或单元数据")
        if len(self.displacements) != 2 * len(self.nodes):
            raise ValueError("位移数据与节点数不匹配")
        if scale_factor <= 0:
            raise ValueError("放大系数必须为正数")
            
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 计算变形后节点坐标
        deformed_nodes = self.nodes + scale_factor * self.displacements.reshape(-1, 2)
        
        # 绘制原始网格(带连接顺序)
        for elem in self.elements:
            # 获取节点坐标
            orig_coords = self.nodes[elem]
            def_coords = deformed_nodes[elem]
            
            # 绘制原始单元(虚线)
            orig_poly = plt.Polygon(orig_coords, closed=True, 
                                  edgecolor='gray', fill=False, 
                                  linestyle=':', linewidth=0.8)
            ax.add_patch(orig_poly)
            
            # 绘制变形单元(实线+箭头指示连接顺序)
            def_poly = plt.Polygon(def_coords, closed=True,
                                 edgecolor='red', fill=False,
                                 linewidth=1.5)
            ax.add_patch(def_poly)
            
            # 绘制连接顺序指示线
            for i in range(len(elem)):
                start = def_coords[i]
                end = def_coords[(i+1)%len(elem)]
                dx = end[0] - start[0]
                dy = end[1] - start[1]
                ax.arrow(start[0], start[1], dx*0.9, dy*0.9, 
                        shape='full', lw=0, length_includes_head=True,
                        head_width=0.1, head_length=0.15, fc='blue', ec='blue')
            
            # 标注节点编号
            if show_node_numbers:
                for i, node in enumerate(elem):
                    ax.text(orig_coords[i,0], orig_coords[i,1], f'n{node}', 
                           fontsize=8, color='gray', ha='center', va='center')
                    ax.text(def_coords[i,0], def_coords[i,1], f'n{node}', 
                           fontsize=8, color='red', ha='center', va='center')
        
        # 添加图例和标题
        ax.set_title(f'Enhanced Mesh Deformation (Magnified x{scale_factor})')
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.legend(['Original mesh', 'Deformed mesh', 'Element node order'], 
                 loc='upper right')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.axis('equal')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_stress_contour(self, component='xx', levels=20, save_path=None):
        """绘制应力云图
        
        参数:
            component: 应力分量('xx', 'yy'或'xy')
            levels: 等值线数量
            save_path: 图片保存路径(可选)
        """
        if self.stresses is None:
            raise ValueError("没有应力数据可供可视化")
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 根据分量选择应力数据
        if component == 'xx':
            stress_data = self.stresses[:, 0]
            title = r'Stress $\sigma_{xx}$ Distribution'
        elif component == 'yy':
            stress_data = self.stresses[:, 1]
            title = r'Stress $\sigma_{yy}$ Distribution'
        elif component == 'xy':
            stress_data = self.stresses[:, 2]
            title = r'Stress $\sigma_{xy}$ Distribution'
        else:
            raise ValueError("无效的应力分量，必须是'xx', 'yy'或'xy'")
        
        # 创建三角网格应力云图
        triangulation = plt.tricontourf(
            self.nodes[:, 0], self.nodes[:, 1], stress_data, 
            levels=levels, cmap='jet')
        
        # 添加颜色条
        cbar = fig.colorbar(triangulation, ax=ax)
        cbar.set_label(f'Stress {component} (Pa)')
        
        # 绘制单元边界
        for elem in self.elements:
            x = self.nodes[elem, 0]
            y = self.nodes[elem, 1]
            ax.fill(x, y, edgecolor='k', fill=False, linewidth=0.5, alpha=0.3)
        
        ax.set_title(title)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.axis('equal')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_displacement(self, component='x', save_path=None):
        """绘制位移分布图
        
        参数:
            component: 位移分量('x'或'y')
            save_path: 图片保存路径(可选)
        """
        if self.displacements is None:
            raise ValueError("没有位移数据可供可视化")
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 根据分量选择位移数据
        if component == 'x':
            disp_data = self.displacements[::2]
            title = 'Displacement in X Direction'
        elif component == 'y':
            disp_data = self.displacements[1::2]
            title = 'Displacement in Y Direction'
        else:
            raise ValueError("无效的位移分量，必须是'x'或'y'")
        
        # 创建位移云图
        triangulation = plt.tricontourf(
            self.nodes[:, 0], self.nodes[:, 1], disp_data, 
            levels=20, cmap='coolwarm')
        
        # 添加颜色条
        cbar = fig.colorbar(triangulation, ax=ax)
        cbar.set_label(f'Displacement {component} (m)')
        
        ax.set_title(title)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.grid(True, linestyle=':', alpha=0.5)
        ax.axis('equal')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
    def plot_all_results(self, scale_factor=50.0, save_path=None, element_type='linear', nx=1, ny=1):
        """绘制所有结果综合图
        
        参数:
            scale_factor: 变形放大系数
            save_path: 图片保存路径
            element_type: 单元类型（'linear'或'quadratic'）
            nx: x方向单元数
            ny: y方向单元数
        
        异常:
            ValueError: 如果输入数据无效
        """
        print("\n=== 可视化数据验证 ===")
        print(f"节点数: {len(self.nodes)}")
        print(f"单元数: {len(self.elements)}")
        print(f"位移向量长度: {len(self.displacements)}")
        
        # 详细验证输入数据
        if self.displacements is None:
            raise ValueError("位移数据为空")
        if self.stresses is None:
            raise ValueError("应力数据为空")
        if self.nodes is None:
            raise ValueError("节点坐标数据为空")
        if self.elements is None:
            raise ValueError("单元连接数据为空")
                   
        # 创建图形
        print("\n创建图形...")
        fig = plt.figure(figsize=(18, 12))
        plt.suptitle("Timoshenko Beam FEM Results", fontsize=16)
        
        # 子图1: 变形网格
        ax1 = plt.subplot(2, 2, 1)
        deformed_nodes = self.nodes + scale_factor * self.displacements.reshape(-1, 2)
        for elem in self.elements:
            x = self.nodes[elem, 0]
            y = self.nodes[elem, 1]
            ax1.fill(x, y, edgecolor='gray', fill=False, linestyle=':', linewidth=0.8)
            xd = deformed_nodes[elem, 0]
            yd = deformed_nodes[elem, 1]
            ax1.fill(xd, yd, edgecolor='red', fill=False, linewidth=1.5)
        ax1.set_title(f'Mesh Deformation (Magnified x{scale_factor})')
        ax1.set_xlabel('x (m)')
        ax1.set_ylabel('y (m)')
        ax1.axis('equal')
        
        # 子图2: 应力云图
        print("\n准备绘制应力云图...")
        
        ax2 = plt.subplot(2, 2, 2)
        tcf = ax2.tricontourf(self.nodes[:,0], self.nodes[:,1], 
                             self.stresses[:,0], levels=20, cmap='jet')
        plt.colorbar(tcf, ax=ax2).set_label(r'Stress $\sigma_{xx}$ (Pa)')
        ax2.set_title(r'Stress $\sigma_{xx}$ Distribution')
        ax2.set_xlabel('x (m)')
        ax2.set_ylabel('y (m)')
        ax2.axis('equal')

        # 子图3: 位移云图
        ax3 = plt.subplot(2, 2, 3)
        tcf = ax3.tricontourf(self.nodes[:,0], self.nodes[:,1], 
                             self.displacements[1::2], levels=20, cmap='coolwarm')
        plt.colorbar(tcf, ax=ax3).set_label('Displacement Y (m)')
        ax3.set_title('Displacement in Y Direction')
        ax3.set_xlabel('x (m)')
        ax3.set_ylabel('y (m)')
        ax3.axis('equal')
        
        # 子图4: 结果摘要
        ax4 = plt.subplot(2, 2, 4)
        ax4.text(0.1, 0.8, f"Max Stress: {np.max(self.stresses[:,0]):.2f} Pa", fontsize=10)
        ax4.text(0.1, 0.7, f"Min Stress: {np.min(self.stresses[:,0]):.2f} Pa", fontsize=10)
        ax4.text(0.1, 0.6, f"Max Displacement: {np.max(self.displacements[1::2]):.2e} m", fontsize=10)
        ax4.text(0.1, 0.5, f"Min Displacement: {np.min(self.displacements[1::2]):.2e} m", fontsize=10)
        ax4.axis('off')
        ax4.set_title('Results Summary')
        
        plt.tight_layout()
        if save_path is None:
            save_path = f"results_{element_type}_nx{nx}_ny{ny}.png"
        plt.savefig(save_path, dpi=300)
        plt.close()