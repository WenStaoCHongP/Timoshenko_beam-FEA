from scipy.sparse.linalg import spsolve, cg, minres
import numpy as np
import time
from scipy.sparse import issparse, csr_matrix

class LinearSolver:
    def __init__(self, method='direct', tol=1e-8, max_iter=1000):
        """初始化线性求解器
        
        参数:
            method: 求解方法('direct', 'cg'或'minres')
            tol: 迭代求解器容差
            max_iter: 最大迭代次数
        """
        self.method = method
        self.tol = tol
        self.max_iter = max_iter
        
    def solve(self, K, F):
        """求解线性系统 K*u = F
        
        参数:
            K: 刚度矩阵(稀疏或密集)
            F: 载荷向量
            
        返回:
            u: 位移解向量
            
        异常:
            RuntimeError: 如果所有求解方法都失败
        """
        print("\nSolving linear system...")
        start_time = time.time()
        
        # 验证输入
        if K is None or F is None:
            raise ValueError("刚度矩阵或载荷向量为None")
        if K.shape[0] != K.shape[1]:
            raise ValueError(f"刚度矩阵不是方阵: 形状{K.shape}")
        if K.shape[0] != F.shape[0]:
            raise ValueError(f"维度不匹配: K{K.shape}, F{F.shape}")
            
        # 确保K是稀疏矩阵(如果不是则转换)
        if not issparse(K):
            print("警告: 刚度矩阵不是稀疏格式，将转换为CSR格式")
            K = csr_matrix(K)
            
        # 根据方法选择求解器
        if self.method == 'direct':
            try:
                print("尝试直接求解...")
                u = self._solve_direct(K, F)
                if u is None:
                    raise RuntimeError("直接求解器返回None")
                solve_time = time.time() - start_time
                print(f"直接求解成功, 耗时: {solve_time:.3f}秒")
                return u
            except Exception as e:
                print(f"直接求解失败: {str(e)}")
                if self.method != 'auto':
                    raise RuntimeError(f"直接求解失败: {str(e)}")
                
        if self.method in ['iterative', 'auto']:
            try:
                print(f"尝试迭代求解(方法={self.method})...")
                u = self._solve_iterative(K, F)
                if u is None:
                    raise RuntimeError("迭代求解器返回None")
                solve_time = time.time() - start_time
                print(f"迭代求解成功, 耗时: {solve_time:.3f}秒")
                return u
            except Exception as e:
                print(f"迭代求解失败: {str(e)}")
                if self.method != 'auto':
                    raise RuntimeError(f"迭代求解失败: {str(e)}")
                
        # 尝试最小二乘法作为最后手段
        try:
            print("尝试最小二乘解...")
            u, residuals, rank, s = np.linalg.lstsq(K.toarray(), F, rcond=None)
            if len(u) != F.shape[0]:
                raise RuntimeError("最小二乘解维度不匹配")
            solve_time = time.time() - start_time
            print(f"最小二乘解完成, 耗时: {solve_time:.3f}秒")
            print(f"残差: {residuals[0] if len(residuals) > 0 else 'N/A'}")
            return u
        except Exception as e:
            print(f"最小二乘解失败: {str(e)}")
            raise RuntimeError(f"所有求解方法都失败: {str(e)}")
            
    def _solve_direct(self, K, F):
        """直接求解方法"""
        return spsolve(K, F)
        
    def _solve_iterative(self, K, F):
        """迭代求解方法"""
        if K.shape[0] < 1000:
            print("系统规模较小(<1000)，建议使用直接求解")
            
        # 根据矩阵对称性选择迭代方法
        if np.allclose(K - K.T, 0, atol=1e-8):
            print("使用共轭梯度法(CG)求解对称系统")
            u, info = cg(K, F, tol=self.tol, maxiter=self.max_iter)
        else:
            print("使用最小残差法(MINRES)求解非对称系统")
            u, info = minres(K, F, tol=self.tol, maxiter=self.max_iter)
            
        if info != 0:
            raise RuntimeError(f"迭代求解未收敛 (info={info})")
            
        return u
        
    def check_solution(self, K, u, F):
        """验证解的正确性"""
        residual = K.dot(u) - F
        rel_error = np.linalg.norm(residual) / np.linalg.norm(F)
        print(f"相对残差: {rel_error:.3e}")
        return rel_error < self.tol