
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import QuantLib as ql
from collections import namedtuple
import math
from pandas import DataFrame
from datetime import datetime


class HW1FSwaption:
    #用于生成HW1F模型下的swaption对象
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
        
        
        








class SimulatedAnnealingSwaption:
    #参数校准
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
        for num_ins in range(self.tenor.shape[0]):
            swaption=HW1FSwaption(self.a,this_sigma,self.maturity[num_ins],self.tenor[num_ins],self.ZeroRate,self.SwaptionRate[num_ins],self.TenorLeg)
            swaption.PriceSwaption()
            self.swaption_price_list=np.append(self.swaption_price_list,swaption.PSwaption)
                

    def CostFunction(self,this_sigma):
        self.SwaptionPriceList(this_sigma)          
        self.costvalue=np.sum((self.weight*self.swaption_price_list-self.MarketPrice)**2)
        return self.costvalue
    
    def Mesh(self):
        sigma_list=np.linspace(0.001,0.05,1000)
        thecost=np.array([])
        for count in range(1000):
            thecost=np.append(thecost,self.CostFunction(sigma_list[count]))
        self.sigma_initial=sigma_list[np.argmin(thecost)]
    
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
            
            
            
            
class HW1FZeroBond:
    #用于生成HW1F模型下的ZeroBond对象
    def __init__(self,a,sigma,ZeroRate,RateInitial,YearLimit):
        self.a=a  #初始化mean reversion，scalar
        self.sigma=sigma
        self.ZeroRate=ZeroRate[:(int(12*YearLimit)+1)]  #vector 一维数组
        self.RateInitial=RateInitial   #scalar
        self.YearLimit=YearLimit  #int scalar,需要校准的年数
        
        #计算time list (1/12,1/6,1/4,1/3,5/12,1/2,7/12,……,)
        self.timelist=np.arange(0,self.YearLimit*12+1)/12
        
        #计算远期瞬时利率
        self.forwardrate=self.ZeroRate#*np.exp(-self.ZeroRate*self.timelist) 
        
        #计算discounting
        self.discouting=np.exp(-self.ZeroRate*self.timelist)
        
        #计算仿射因子B
        self.B_list=(1-np.exp(-self.a*self.timelist))/self.a
        
        
    def SpecialTerm(self):
        term1=self.sigma**2/(4*self.a**3)
        term2=-np.exp(-2*self.a*self.timelist)
        term3=4*np.exp(-self.a*self.timelist)
        term4=2*self.a*self.timelist
        self.specialterm=np.exp(term1*(term2+term3+term4-3))
        
        
    #计算theta scalar
    def Theta(self):
        self.Vst=2*self.sigma**2*(0.5+0.5*np.exp(-2*self.a*self.timelist)-np.exp(-self.a*self.timelist))/self.a**2
        self.Vnd=2*self.sigma**2*np.exp(-self.a*self.timelist)*(1-np.exp(-self.a*self.timelist))/self.a
        self.theta=np.gradient(self.forwardrate,1/12)+self.a*self.forwardrate+0.5*(self.Vnd+self.a*self.Vst)
        self.theta=self.theta/self.a
        
    def ZeroBondPrice(self):
        self.SpecialTerm()
        self.Theta()
        self.zero_bond_price=np.exp((self.theta-self.RateInitial)*self.B_list-self.theta*self.timelist)*self.specialterm





#数据处理，用于输出Matrix(9)*Tenor(8) 矩阵
class Matrixisation:     
    #input请输入一条72元素的一维数组
    def __init__(self,input):
        self.data=input
        
    def Output(self):
        self.matrix=np.full((9,8),1.0)
        for mat in range(9):
            self.matrix[mat,:]=self.data[mat::9]
            
            
            
            
            
            
            
            
# 以下函数用于单独校准均值回归a    
#zero bond
def P(rate,t,T):
    return np.exp(-rate * (T-t))

# To calculate the affine factor B(t,T)
def B_term(a,t,T):
    B=(1/a)*(1-np.exp(-a*(T-t)))
    return B



#SMM Ratio
def ApproximationSMMRatio(a,sigma,maturity,tenor,rate):
    #“rate”处应当输入每日zero_rate列
    #“maturity”、“tenor”处应当输入需要校准的instrument的tenor和maturity矩阵
    ApproSMM=np.full((62,72),1.0)
    for m in range(0,62):
        for n in range(0,72):
            B=B_term(a,maturity[m,n],maturity[m,n]+tenor[m,n])
            Vr=(0.5/a)*sigma**2*(1-np.exp(-2*a*maturity[m,n]))
            Vp=Vr*(B**2)
            P_0=P(rate[m,int(maturity[m,n]*12)],0,maturity[m,n])
            P_n=P(rate[m,int((maturity[m,n]+tenor[m,n])*12)],0,maturity[m,n]+tenor[m,n])
            P_ratio=P_0/(P_0-P_n)
            
            ApproSMM[m,n]=Vp*(P_ratio**2)

    ApproSMMRatio=np.full((62,63),1.0)
    for k in range(0,62):
        for j in range(0,9):
            for i in range(0,7):
                ApproSMMRatio[k,j+9*i]=ApproSMM[k,j+9*(i+1)]/ApproSMM[k,j+9*i]
                
    return ApproSMMRatio



def cost_function_for_a(a,sigma,tenor,maturity,rate,SwaptionRatio):
    ApproRatio=ApproximationSMMRatio(a,sigma,maturity,tenor,rate)
    cost=(SwaptionRatio-ApproRatio**0.5)**2
    cost=np.sum(cost,axis=1, keepdims=True)
    return cost

    
def SimulatedAnnealing_for_a(initial,sigma,tenor,maturity,ZeroRate,SwaptionRatio,mean,std,interation,gama):
    a_calibrated=[]
    cost_list=[]
    
    for day in range(0,62):
        a_initial=initial
        a_run=a_initial
        a_min=a_initial
        fr_initial=cost_function_for_a(a_run,sigma,tenor,maturity,ZeroRate,SwaptionRatio)
        fr_min=fr_initial
        std_initial=std
        std_run=std
        for n_interation in range(interation):
            a_run=a_run+np.random.normal(mean, std_run, 1)
            fr_run=cost_function_for_a(a_run,sigma,tenor,maturity,ZeroRate,SwaptionRatio)
            if(fr_run[day]<fr_min[day]):
                fr_min=fr_run
                a_min=a_run
                std_run=std_initial*np.exp(-gama*n_interation/(interation-1))
        a_calibrated=np.append(a_calibrated,a_min)
        cost_list=np.append(cost_list,fr_min[day])
        
    
    return a_calibrated,cost_list



def create_swaption_helpers(data, index, term_structure, engine):
    swaptions = []
    fixed_leg_tenor = ql.Period(3, ql.Months) 
    fixed_leg_daycounter = ql.Actual360()
    floating_leg_daycounter = ql.Actual360()# 浮动和固定利率的天数计算方式
    for d in data:
        vol_handle = ql.QuoteHandle(ql.SimpleQuote(d.volatility))
        helper = ql.SwaptionHelper(ql.Period(int(d.start), ql.Years),ql.Period(int(d.length), ql.Years),vol_handle,index,fixed_leg_tenor,fixed_leg_daycounter,floating_leg_daycounter,term_structure)
        helper.setPricingEngine(engine)
        swaptions.append(helper)
    return swaptions

def calibration_report(swaptions, data):
    columns = ["Model Price", "Market Price", "Implied Vol", "Market Vol","Rel Error Price", "Rel Error Vols"]
    report_data = []
    cum_err = 0.0
    cum_err2 = 0.0
    for i, s in enumerate(swaptions):
        model_price = s.modelValue()
        market_vol = data[i].volatility
        black_price = s.blackPrice(market_vol)
        rel_error = model_price/black_price - 1.0
        implied_vol = s.impliedVolatility(model_price,1e-5, 50, 0.0, 0.50)
        rel_error2 = implied_vol/market_vol-1.0
        cum_err += rel_error*rel_error
        cum_err2 += rel_error2*rel_error2
        report_data.append((model_price, black_price, implied_vol,market_vol, rel_error, rel_error2))
    print("Cumulative Error Price: %7.5f" % math.sqrt(cum_err))
    print("Cumulative Error Vols : %7.5f" % math.sqrt(cum_err2))
    return DataFrame(report_data, columns=columns,index=['']*len(report_data))




