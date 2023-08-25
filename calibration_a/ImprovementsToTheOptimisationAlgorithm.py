
import numpy as np
import pandas as pd
import time  
from scipy.interpolate import interp1d
from scipy.stats import truncnorm
from scipy.stats import norm
import matplotlib.pyplot as plt
from pandas import DataFrame
import matplotlib.dates as mdates
from datetime import datetime


swaption_vol = np.array(pd.read_excel('swaption2023.xlsx', sheet_name='Sheet1').iloc[1:64,1:73], dtype=float) * 0.0001


#To calculate the Mkt vol ratio
SwaptionRatio=np.full((62,63),1.0)
for k in range(0,62):
    for j in range(0,9):
        for i in range(0,7):
            SwaptionRatio[k,j+9*i]=swaption_vol[k,j+9*(i+1)]/swaption_vol[k,j+9*i]
            
            
#Maturity & Tenor
tenor = np.array([[1.0]*9 + [2.0]*9 + [3.0]*9 + [4.0]*9 + [5.0]*9 + [10.0]*9 + [15.0]*9 + [20.0]*9 for _ in range(62)])
maturity = np.full((62,72),1.0)
maturity = np.tile(np.tile([1,2,3,4,5,6,7,8,10], 8), (62, 1))

sigma = 1.0   #For the convenience of calculation, but not used (eliminated).

ZeroRate = np.array(pd.read_excel('ZeroRate.xlsx', sheet_name='Sheet1').iloc[0:62,0:601], dtype=float) * 0.01
#ZeroRate.shape=(63, 241)




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



def cost_function(a,sigma,tenor,maturity,rate,SwaptionRatio):
    ApproRatio=ApproximationSMMRatio(a,sigma,maturity,tenor,rate)
    cost=(SwaptionRatio-ApproRatio**0.5)**2
    cost=np.sum(cost,axis=1, keepdims=True)
    return cost

    
def SimulatedAnnealing(initial,sigma,tenor,maturity,ZeroRate,SwaptionRatio,mean,std,interation,gama):
    a_calibrated=[]
    cost_list=[]
    
    for day in range(0,62):
        a_initial=initial
        a_run=a_initial
        a_min=a_initial
        fr_initial=cost_function(a_run,sigma,tenor,maturity,ZeroRate,SwaptionRatio)
        fr_min=fr_initial
        std_initial=std
        std_run=std
        for n_interation in range(interation):
            a_run=a_run+np.random.normal(mean, std_run, 1)
            fr_run=cost_function(a_run,sigma,tenor,maturity,ZeroRate,SwaptionRatio)
            if(fr_run[day]<fr_min[day]):
                fr_min=fr_run
                a_min=a_run
                std_run=std_initial*np.exp(-gama*n_interation/(interation-1))
        a_calibrated=np.append(a_calibrated,a_min)
        cost_list=np.append(cost_list,fr_min[day])
        
    
    return a_calibrated,cost_list




initial=0.08
mean=0
std=0.002
interation=100
gama=0.01

a_calibrated,cost_value=SimulatedAnnealing(initial,sigma,tenor,maturity,ZeroRate,SwaptionRatio,mean,std,interation,gama)


#粗糙判断可能达到最优点的位置，以此作为优化算法的初值
a_list=np.linspace(0.06,0.1,400)    
day_list=np.arange(62)
cost_value=np.full((62,400),1.0)

for count in range(400):
        cost_value[:,count]=cost_function(a_list[count],sigma,tenor,maturity,ZeroRate,SwaptionRatio)[:, 0]

listmin_cost=np.array([])
listmin_a=np.array([])
for m in range(62):
    listmin_cost=np.append(listmin_cost,np.min(cost_value[m,:]))
    listmin_a=np.append(listmin_a,a_list[np.argmin(cost_value[m,:])])
    
x=listmin_a  



initial=listmin_a 
mean=0
std=0.0001
interation=50
gama=0.09

a_calibrated_initial=a_calibrated
cost_value_initial=cost_value

a_calibrated_new,cost_value_new=SimulatedAnnealing(initial,sigma,tenor,maturity,ZeroRate,SwaptionRatio,mean,std,interation,gama)
cost_value_run = np.minimum(cost_value_new, cost_value_initial)
a_calibrated_run=np.where(cost_value_new <= cost_value_initial, a_calibrated_new, a_calibrated_initial)

for interation_n in range(interation):
    a_calibrated_new,cost_value_new=SimulatedAnnealing(initial,sigma,tenor,maturity,ZeroRate,SwaptionRatio,mean,std,interation,gama)
    cost_value_run = np.minimum(cost_value_new, cost_value_run)
    a_calibrated_run=np.where(cost_value_run <= cost_value_new, a_calibrated_run, a_calibrated_new)



