import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import interp1d
import time




timelist=np.arange(0, 601) / 12
ZeroRate = np.array(pd.read_excel('ZeroRate.xlsx', sheet_name='Sheet1').iloc[0:62,0:601], dtype=float) * 0.01
#ZeroRate.shape=(63, 241)
    

swaptionrate = np.array(pd.read_excel('swaption2023.xlsx', sheet_name='Sheet1').iloc[1:64,77:85], dtype=float) * 0.01
#插值补全缺失值
x = np.linspace(0, 61, 62)
for i in range(8):
    y = swaptionrate[:, i]
    x_clean = x[~np.isnan(y)]
    y_clean = y[~np.isnan(y)]
    f = interp1d(x_clean, y_clean, kind='linear')
    swaptionrate[:, i] = f(x)
#swaptionrate.shape=(62, 8)

swaptionrateglobal = np.full((62,72),1.0)
col_index = 0
for s in range(72):
    swaptionrateglobal[:, s] = swaptionrate[:, col_index]
    if (s + 1) % 9 == 0:
        col_index += 1
        
swaptionrate=swaptionrateglobal
#swaptionrate.shape=(62, 72)


# Input Market vol
Market_Price=np.array(pd.read_excel('MarketPrice.xlsx', sheet_name='Sheet1').iloc[0:62,0:72], dtype=float)

a_calibrated = np.array(pd.read_excel('a_calibrated.xlsx', sheet_name='Sheet1').iloc[0:62,0], dtype=float)  
#a_calibrated.shape=(62,1)

tenor = np.array([1.0]*9 + [2.0]*9 + [3.0]*9 + [4.0]*9 + [5.0]*9 + [10.0]*9 + [15.0]*9 + [20.0]*9)

maturity = np.array([1,2,3,4,5,6,7,8,10]*8)


def P(rate,t,T):
    return np.exp(-rate * (T-t))
#远期瞬时利率
def f(rate,timelist): 
    return (rate*np.exp(-rate*timelist))   
#f(0,T)=-d(exp(-rT))/dT=r*exp(-rT)
instanforwardrate=f(ZeroRate,timelist)
#instanforwardrate.shape=(62,601)  



def r_star_equation(a,sigma,maturity,tenor,ZeroRate,instanforwardrate,r,swaptionrate):
    #函数r_star_equation中应当输入的参数如下：
    #Scalar：a,sigma,maturity,tenor,swaptionrate
    #Vector:ZeroRate.shape=(1,601),instanforwardrate.shape=(1,601)
    tenor_list=([])
    for num_tenor in range(int(4*tenor)):
        tenor_list=np.append(tenor_list,(num_tenor+1)/4)
        
    B_list=(1/a)*(1-np.exp(-a*tenor_list))
    Vr=((sigma**2)/(2*a))*(1-np.exp(-2*a*maturity))
    
    P_T_Zero=np.exp(-ZeroRate[int(12*maturity)]*maturity)
    ZeroRate_list=([])
    for num_ZeroRate in range(int(4*tenor)):
        ZeroRate_list=np.append(ZeroRate_list,ZeroRate[int(12*maturity+3*(1+num_ZeroRate))])
    P_Tn_list=np.exp(-ZeroRate_list*(maturity+tenor_list))
    
    A_list=np.log(P_Tn_list/P_T_Zero)+B_list*instanforwardrate[int(12*maturity)]-0.5*(B_list**2)*Vr
    X_list=np.exp(A_list-B_list*r)
    
    c_list=np.full(int(tenor*4-1),swaptionrate/4)
    c_list=np.append(c_list,1+swaptionrate/4)
    
    return np.sum(c_list*X_list)
    


    def root_r_star(a,sigma,ZeroRate,maturity,tenor,swaptionrate,instanforwardrate,tolerance=1e-6, max_iterations=50):
    #函数r_star_equation中应当输入的参数如下：
    #Scalar：a,sigma,maturity,tenor,swaptionrate
    #Vector:ZeroRate.shape=(1,601),instanforwardrate.shape=(1,601)
    
    lower_bound = 0.0  # 下界
    upper_bound = 0.5   # 上界
    
    for i in range(max_iterations):
        r_star = (lower_bound + upper_bound) / 2
        result = r_star_equation(a,sigma,maturity,tenor,ZeroRate,instanforwardrate,r_star,swaptionrate)
        
        if abs(result-1)<tolerance:
            
            return r_star
            
        if result > 1:
            lower_bound = r_star
        else:
            upper_bound = r_star  
        
    return r_star



def PSwaption(a,sigma,maturity,tenor,ZeroRate,instanforwardrate,r_star,swaptionrate):
    #函数PSwaption中应当输入的参数如下：
    #Scalar：a,sigma,maturity,tenor,r_star,swaptionrate
    #Vector:ZeroRate.shape=(1,601),instanforwardrate.shape=(1,601)
    
    tenor_list=[]
    for num_tenor in range(int(4*tenor)):
        tenor_list=np.append(tenor_list,(num_tenor+1)/4)
        
    B_list=(1/a)*(1-np.exp(-a*tenor_list))
    Vr=((sigma**2)/(2*a))*(1-np.exp(-2*a*maturity))
    Vp=Vr*(B_list**2)
    Vp_sqrt=Vp**0.5
    
    P_t=np.exp(-ZeroRate[int(12*maturity)]*maturity)
    ZeroRate_list=[]
    for num_ZeroRate in range(int(4*tenor)):
        ZeroRate_list=np.append(ZeroRate_list,ZeroRate[int(12*maturity+3*(1+num_ZeroRate))])
    P_Tn_list=np.exp(-ZeroRate_list*(maturity+tenor_list))
    
    A_list=np.log(P_Tn_list/P_t)+B_list*instanforwardrate[int(12*maturity)]-0.5*(B_list**2)*Vr
    X_list=np.exp(A_list-B_list*r_star)
    
    d_plus_list=np.log(X_list*P_t/P_Tn_list)/Vp_sqrt+0.5*Vp_sqrt
    d_minus_list=np.log(X_list*P_t/P_Tn_list)/Vp_sqrt-0.5*Vp_sqrt
    d_plusnorm_list=norm.cdf(d_plus_list, 0, 1)
    d_minusnorm_list=norm.cdf(d_minus_list, 0, 1)
    
    ZBP_list=X_list*P_t*d_plusnorm_list-P_Tn_list*d_minusnorm_list
    
    c_list=np.full(int(tenor*4-1),swaptionrate/4)
    c_list=np.append(c_list,1+swaptionrate/4)
    
    PriceSwaption=np.sum(c_list*ZBP_list)
    
    return PriceSwaption



def costfunction(a,sigma,maturity,tenor,ZeroRate,instanforwardrate,swaptionrate,Market_Price):
    #函数costfunction中应当输入的参数如下：
    #Scalar：a,sigma,r_star
    #Vector:ZeroRate.shape=(1,601),instanforwardrate.shape=(1,601)
    #Vector:swaptionrate.shape=(1,72),maturity.shape=(1,72),tenor.shape=(1,72),Market_Price.shape=(1,72)
    this_PSwaption_list=np.array([])
    for type_of_swap in range(72):
        this_r_star=root_r_star(a,sigma,ZeroRate,maturity[type_of_swap],tenor[type_of_swap],swaptionrate[type_of_swap],instanforwardrate,tolerance=1e-6, max_iterations=50)
        this_PSwaption=PSwaption(a,sigma,maturity[type_of_swap],tenor[type_of_swap],ZeroRate,instanforwardrate,this_r_star,swaptionrate[type_of_swap])
        this_PSwaption_list=np.append(this_PSwaption_list,this_PSwaption)
        weight=np.reciprocal(np.max(tenor)) 
    cost_list=(this_PSwaption_list*weight-Market_Price)**2
    return np.sum(cost_list)   



def SimulatedAnnealing(a,initial_sigma,maturity,tenor,ZeroRate,instanforwardrate,swaptionrate,Market_Price,mean,std,interation,gama):
    #函数costfunction中应当输入的参数如下：
    #Scalar：a,sigma,r_star
    #Vector:ZeroRate.shape=(1,601),instanforwardrate.shape=(1,601)
    #Vector:swaptionrate.shape=(1,72),maturity.shape=(1,72),tenor.shape=(1,72),Market_Price.shape=(1,72)
    
    sigma_initial=initial_sigma
    sigma_run=sigma_initial
    sigma_min=sigma_initial
    function_run_initial=costfunction(a,sigma_run,maturity,tenor,ZeroRate,instanforwardrate,swaptionrate,Market_Price)
    function_run_min=function_run_initial
    std_initial=std
    std_run=std
    for n_interation in range(interation):
        sigma_run+=np.random.normal(mean, std_run, 1)
        function_run=costfunction(a,sigma_run,maturity,tenor,ZeroRate,instanforwardrate,swaptionrate,Market_Price)
        if(function_run<function_run_min):
            function_run_min=function_run
            sigma_min=sigma_run
            std_run=std_initial*np.exp(-gama*n_interation/(interation-1))
            
    return sigma_min,function_run_min
    
 
