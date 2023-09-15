import QuantLib as ql
from collections import namedtuple
import math
from pandas import DataFrame
import pandas as pd
import numpy as np
from datetime import datetime


#volatility
swaption_vol = np.array(pd.read_excel('swaption2023.xlsx', sheet_name='Sheet1').iloc[4,1:73], dtype=float) * 0.0001
swaption_vol = np.round(swaption_vol, 8)
#swaption_vol.shape=(72,)


#zero rate curve日期
curve_start_date = ql.Date(28, 4, 2023)  
curve_date = [curve_start_date]
current_date = curve_start_date
while len(curve_date) < 601:
    current_date = ql.Date(current_date.dayOfMonth(), current_date.month(), current_date.year())
    current_date += ql.Period(1, ql.Months)
    curve_date.append(current_date)
#len(curve_date)=601


#zero rate
ZeroRate = np.array(pd.read_excel('ZeroRate.xlsx', sheet_name='Sheet1').iloc[0:62,0:601], dtype=float) * 0.01
#ZeroRate.shape=(62, 601)


#Term Structure
today = ql.Date(3, ql.May, 2023)
ql.Settings.instance().evaluationDate = today
term_structure = ql.YieldTermStructureHandle(ql.ForwardCurve(curve_date, ZeroRate[0,:], ql.Actual360()))
index = ql.USDLibor(ql.Period('6M'),term_structure)



CalibrationData = namedtuple("CalibrationData","start, length, volatility")

#Maturity（star）、Tenor（length）序列


Maturity = np.tile(np.array([1,2,3,4,5,6,7,8,10]), 8)

Tenor = np.repeat([1,2,3,4,5,10,15,20], 9)

#这里全局校准9*8=72种instument
data=[]
for i in range(0,72):
    
    calibration_data = CalibrationData(Maturity[i], Tenor[i], swaption_vol[i] )
    data.append(calibration_data)



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


model = ql.HullWhite(term_structure)
engine = ql.JamshidianSwaptionEngine(model)
swaptions = create_swaption_helpers(data, index, term_structure, engine)
optimization_method = ql.LevenbergMarquardt(1.0e-8,1.0e-8,1.0e-8)
end_criteria = ql.EndCriteria(10000, 100, 1e-6, 1e-8, 1e-8)
model.calibrate(swaptions, optimization_method, end_criteria)
a, sigma = model.params()
print("a = %6.5f, sigma = %6.5f" % (a, sigma))


calibration_report(swaptions, data)




