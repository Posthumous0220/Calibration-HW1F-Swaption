import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import interp1d



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
SwaptionRate=swaptionrateglobal

a_calibrated = np.array(pd.read_excel('a_calibrated.xlsx', sheet_name='Sheet1').iloc[0:62,0], dtype=float)  
#a_calibrated.shape=(62,1)

tenor = np.array([1.0]*9 + [2.0]*9 + [3.0]*9 + [4.0]*9 + [5.0]*9 + [10.0]*9 + [15.0]*9 + [20.0]*9)
maturity = np.array([1,2,3,4,5,6,7,8,10]*8)
weight=np.array([1.0]*2 +[5.0]*7+ [10.0]*9 + [15.0]*9 + [20.0]*9 + [20.0]*9 + [25.0]*9 + [22.0]*9 + [20.0]*9)


Market_Price=np.array(pd.read_excel('MarketPrice.xlsx', sheet_name='Sheet1').iloc[0:62,0:72], dtype=float)
listmin_sigma=np.array(pd.read_excel('listmin_sigma_as_initial.xlsx', sheet_name='Sheet1'), dtype=float)  
listmin_cost=np.array(pd.read_excel('listmin_cost_as_initial.xlsx', sheet_name='Sheet1'), dtype=float)  



#每一个类某一天的某类instrument(swaption)

class HW1FSwaption:
    
    def __init__(self,a,sigma,maturity,tenor,ZeroRate,SwaptionRate,TenorLeg):
        self.a=a #初始化mean reversion，scalar
        self.sigma=sigma
        self.maturity=maturity  #scalar
        self.tenor=tenor   #scalar
        self.ZeroRate=ZeroRate  #vector 一维数组  月频
        self.SwaptionRate=SwaptionRate  #scalar
        self.TenorLeg=TenorLeg   #scalar,int,互换频率,例如：每三个月互换一次，则输入为4
        
        #计算远期瞬时利率
        self.forwardrate=self.ZeroRate*np.exp(-self.ZeroRate*(np.arange(0,(self.ZeroRate).shape[0])/12))  
        
        #计算zero bond组合需要的tenor列表
        self.tenor_list=np.array([])   
        for num_tenor in range(int(self.TenorLeg*self.tenor)):
            self.tenor_list=np.append(self.tenor_list,(num_tenor+1)/self.TenorLeg)
            
 
        #计算c_list
        self.c_list=np.full(int(self.tenor*self.TenorLeg-1),self.SwaptionRate/self.TenorLeg)
        self.c_list=np.append(self.c_list,1+self.SwaptionRate/self.TenorLeg)
        
        
        #计算P(0,T0)、P(0,Tn)
        self.P_T_0=np.exp(-self.ZeroRate[int(12*self.maturity)]*self.maturity)
        self.ZeroRate_list=np.array([])
        for num_ZeroRate in range(int(self.TenorLeg*self.tenor)):
            self.ZeroRate_list=np.append(self.ZeroRate_list,self.ZeroRate[int(12*self.maturity+(12/self.TenorLeg)*(1+num_ZeroRate))])
        self.P_Tn_list=np.exp(-self.ZeroRate_list*(self.maturity+self.tenor_list))        
        
 
        
        
    def VarShorRate(self):
        #计算短期利率的方差
        self.Vr=((self.sigma**2)/(2*self.a))*(1-np.exp(-2*self.a*self.maturity))

        
    def AffineB(self):  #计算仿射因子B(t,T)
        self.B_list=(1/self.a)*(1-np.exp(-self.a*self.tenor_list))
        
        
    def AffineA(self):  #计算仿射因子A(t,T)
        self.AffineB()  
        self.VarShorRate()
        self.A_list=np.log(self.P_Tn_list/self.P_T_0)+self.B_list*self.forwardrate[int(12*self.maturity)]-0.5*(self.B_list**2)*self.Vr

    def Root_r_Star(self):   #返回根r*
    
        lower_bound = 0.0  # 下界
        upper_bound = 0.2   # 上界
        
        self.AffineA()
        
        for i in range(30):
            r=(lower_bound + upper_bound) / 2
            result=np.sum(self.c_list*np.exp(self.A_list-self.B_list*r))
            
            if abs(result-1)<1e-6:
                self.r_star=r
                return self.r_star
            if result>1:
                lower_bound=r
            else:
                upper_bound=r
                
        self.r_star=r
        return self.r_star
    
    def PriceSwaption(self):
        self.AffineA()
        self.r_final=self.Root_r_Star()
        
        self.Vp=self.Vr*(self.B_list**2)
        self.Vp_sqrt=self.Vp**0.5
        self.X_list=np.exp(self.A_list-self.B_list*self.r_final)
    
        self.d_plus_list=np.log(self.X_list*self.P_T_0/self.P_Tn_list)/self.Vp_sqrt+0.5*self.Vp_sqrt
        self.d_minus_list=np.log(self.X_list*self.P_T_0/self.P_Tn_list)/self.Vp_sqrt-0.5*self.Vp_sqrt
        self.d_plusnorm_list=norm.cdf(self.d_plus_list, 0, 1)
        self.d_minusnorm_list=norm.cdf(self.d_minus_list, 0, 1)
        
        self.ZBP_list=self.X_list*self.P_T_0*self.d_plusnorm_list-self.P_Tn_list*self.d_minusnorm_list
        self.PSwaption=np.sum(self.c_list*self.ZBP_list)
        


# In[]:  
class SimulatedAnnealing:
    def __init__(self,a,sigma_initial,maturity,tenor,ZeroRate,SwaptionRate,TenorLeg,MarketPrice,weight,mean,std,interation,gama):
        
        self.a=a  #scalar
        self.sigma_initial=sigma_initial   #scalar
        self.maturity=maturity  #vector 一维数组
        self.tenor=tenor   #vector 一维数组
        self.ZeroRate=ZeroRate  #vector 一维数组
        self.SwaptionRate=SwaptionRate  #vector 一维数组
        self.TenorLeg=TenorLeg   #scalar int
        self.MarketPrice=MarketPrice   #vector 一维数组
        self.weight=weight   #vector 一维数组
        self.std=std      #scalar
        self.mean=mean    #scalar
        self.gama=gama    #scalar
        self.interation=interation   #scalar

        
    def SwaptionPriceList(self,this_sigma):
    	#输出一列swaption价格，self.swaption_price_list.shape=tenor.shape，通过循环调用HW1FSwaption类中的PriceSwaption函数
        #这里的this_maturity,this_tenor均为一维数组，需要输入校准当天所有instrument的信息
        self.swaption_price_list=np.array([])
        for num_ins in range(tenor.shape[0]):
            swaption=HW1FSwaption(self.a,this_sigma,self.maturity[num_ins],self.tenor[num_ins],self.ZeroRate,self.SwaptionRate[num_ins],self.TenorLeg)
            swaption.PriceSwaption()
            self.swaption_price_list=np.append(self.swaption_price_list,swaption.PSwaption)
                

    def CostFunction(self,this_sigma):
        self.SwaptionPriceList(this_sigma)          
        self.costvalue=np.sum((self.swaption_price_list-self.weight*self.MarketPrice)**2)
        return self.costvalue
    
    def Optimiser(self):
            sigma_run=self.sigma_initial
            sigma_min=self.sigma_initial
            function_run_initial=self.CostFunction(sigma_run)
            function_run_min=function_run_initial
            std_initial=self.std
            std_run=self.std
            for n_interation in range(self.interation):
                sigma_run+=np.random.normal(self.mean,std_run,1)
                function_run=self.CostFunction(sigma_run)
                if(function_run<function_run_min):
                    function_run_min=function_run
                    sigma_min=sigma_run
                    std_run=std_initial*np.exp(-self.gama*n_interation/(self.interation-1))
                    
            self.sigma_calibrated=sigma_min
            self.cost_value=function_run_min    
