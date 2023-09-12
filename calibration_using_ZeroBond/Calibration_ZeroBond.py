

class HW1FZeroBond:
    def __init__(self,a,sigma,ZeroRate,RateInitial,YearLimit,Step):
        self.a=a  #初始化mean reversion，scalar
        self.sigma=sigma
        self.ZeroRate=ZeroRate[:(int(12*YearLimit)+1)]  #vector 一维数组
        self.RateInitial=RateInitial   #scalar
        self.YearLimit=YearLimit  #int scalar,需要校准的年数
        self.Step=Step  #scalar 优化步长

        
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
  





  class Algorithmboth:
    def __init__(self,ZeroRate,RateInitial,YearLimit,Upper_a,Upper_sigma,Lower_a,Lower_sigma):
        #self.a=a  #初始化mean reversion，scalar
        #self.sigma=sigma
        self.ZeroRate=ZeroRate[:(int(12*YearLimit)+1)]  #vector 一维数组
        self.RateInitial=RateInitial   #scalar
        self.YearLimit=YearLimit  #int scalar,需要校准的年数
        
        #计算time list (1/12,1/6,1/4,1/3,5/12,1/2,7/12,……,)
        self.timelist=np.arange(0,self.YearLimit*12+1)/12

        #参数范围 scalar
        self.Upper_sigma=Upper_sigma
        self.Upper_a=Upper_a
        self.Lower_a=Lower_a
        self.Lower_sigma=Lower_sigma


        
    def CostFunction(self,this_a,this_sigma):
        ZeroBond=HW1FZeroBond(this_a,this_sigma,self.ZeroRate,self.RateInitial,self.YearLimit)
        ZeroBond.ZeroBondPrice()
        self.ModelPrice=ZeroBond.zero_bond_price
        self.costvalue=np.sum((self.ModelPrice-np.exp(-self.ZeroRate*self.timelist))**2)
        return self.costvalue
        
    def Optimiser(self):
        
        a_range = np.linspace(self.Lower_a,self.Upper_a, int(1/self.Step))
        sigma_range = np.linspace(self.Lower_sigma, self.Upper_sigma, int(1/self.Step))
        self.cost_value_matrix=np.full((int(1/self.Step),int(1/self.Step)), 0.1)
        
        #a_calibrated=0
        #sigma_calibrated=0
        
        for m_a in range(int(1/self.Step)):
            for n_sigma in range(int(1/self.Step)):
                cost_value_run=self.CostFunction(a_range[m_a],sigma_range[n_sigma])
                self.cost_value_matrix[m_a,n_sigma]=cost_value_run
                
        self.min_cost_value=np.min(self.cost_value_matrix)
        self.min_index=np.unravel_index(np.argmin(self.cost_value_matrix), self.cost_value_matrix.shape)
        self.a_calibrated=a_range[self.min_index[0]]
        self.sigma_calibrated=sigma_range[self.min_index[1]]