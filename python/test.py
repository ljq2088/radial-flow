import sys                                                                                                                                                                                        
sys.path.insert(0, '/home/ljq/code/radial_flow')                                                                                                                                                  
                                                                                                                                                                                                
from python.radial_flow import SchwarzschildParams, WaveParams, RadialFlowSolver                                                                                                                  
                                                                                                                                                                                                
# 设置黑洞参数                                                                                                                                                                                    
bh = SchwarzschildParams(r0=2.0)  # Schwarzschild 半径                                                                                                                                            
                                                                                                                                                                                                
# 设置波参数                                                                                                                                                                                      
# 简化版本 (V = l(l+1)/r², 保留 -1/r 项)                                                                                                                                                          
wave_simple = WaveParams(omega=0.2, l=0, use_full_potential=False)                                                                                                                                
                                                                                                                                                                                                
# 完整版本 (V = l(l+1)/r² + f'/r, 去掉 a 系数的 -1/r 项)                                                                                                                                          
wave_full = WaveParams(omega=0.2, l=0, use_full_potential=True)                                                                                                                                   
                                                                                                                                                                                                
# 创建求解器                                                                                                                                                                                      
solver = RadialFlowSolver(bh, wave_simple)  # 或 wave_full                                                                                                                                        
                                                                                                                                                                                                
# 求解                                                                                                                                                                                            
result = solver.solve(                                                                                                                                                                            
    r_max=200.0,           # 最大积分半径                                                                                                                                                         
    epsilon=1e-6,          # 视界偏移                                                                                                                                                             
    n_points=100000,       # 输出点数                                                                                                                                                             
    convergence_threshold=0.01,  # 收敛阈值                                                                                                                                                       
    n_periods=2,           # 检测周期数                                                                                                                                                           
    max_r_factor=2.0       # 最大扩展倍数                                                                                                                                                         
)                                                                                                                                                                                                 
                                                                                                                                                                                                
# 查看结果                                                                                                                                                                                        
print(f"收敛状态: {result['success']}")                                                                                                                                                           
print(f"反射系数 A/B = {result['reflection_coeff']}")                                                                                                                                             
