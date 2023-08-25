# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 09:24:44 2023

@author: Posthumous0220
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import interp1d
import time
    

timelist=np.arange(0, 601) / 12
ZeroRate = np.array(pd.read_excel('ZeroRate.xlsx', sheet_name='Sheet1').iloc[0:62,0:601], dtype=float) * 0.01
#ZeroRate.shape=(63, 241)
    

#Fixed Swaprate
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
        
SwapRate=swaptionrateglobal
#swaptionrate.shape=(62, 72)


swaption_vol = np.array(pd.read_excel('swaption2023.xlsx', sheet_name='Sheet1').iloc[1:63,1:73], dtype=float) * 0.0001
swaption_vol = np.round(swaption_vol, 8)
#swaption_vol.shape=(62, 72)

tenor = np.array([1.0]*9 + [2.0]*9 + [3.0]*9 + [4.0]*9 + [5.0]*9 + [10.0]*9 + [15.0]*9 + [20.0]*9)
maturity = np.array([1,2,3,4,5,6,8,10]*9)

def P(rate,t,T):
    return np.exp(-rate * (T-t))


def MargrabePrice(ZeroRate,SwapRate,vol,maturity,tenor):
    #Scalar:SwapRate,vol,maturity,tenor
    #Vector:ZeroRate.shape=(1,601)
    float_value=P(ZeroRate[int(12*maturity)],0,maturity)-P(ZeroRate[int(12*(maturity+tenor))],0,maturity+tenor)
    
    tenor_list=([])
    for num_tenor in range(int(4*tenor)):
        tenor_list=np.append(tenor_list,(num_tenor+1)/4)
        
    ZeroRate_list=([])
    for num_ZeroRate in range(int(4*tenor)):
        ZeroRate_list=np.append(ZeroRate_list,ZeroRate[int(12*maturity+3*(1+num_ZeroRate))])
    
    P_Tn_list=np.exp(-ZeroRate_list*(maturity+tenor_list)) 
    
    fixed_value=np.sum(P_Tn_list)*SwapRate/4
    
    d1=(np.log(float_value/fixed_value)+0.5*maturity*vol**2)/(vol*maturity**0.5)
    d2=d1-vol*maturity**0.5
    d1_norm=norm.cdf(d1, 0, 1)
    d2_norm=norm.cdf(d2, 0, 1)
    
    price=float_value*d1_norm-fixed_value*d2_norm
    
    return price
    
