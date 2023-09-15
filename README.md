# Calibration-HW1F-Swaption

本项目复现了论文 ["Calibration Methods of Hull-White Model"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1514192)中实现的Hull-White模型常数参数的校准方法。原论文基于2007年至2009年期间日本市场的交易数据，共选取120种校准工具(Swaption)，分别使用三种校准策略，对常数形式和时间相关的参数进行校准。


## 校准策略

* 策略一： 只对波动率$\sigma$进行校准；
* 策略二： 先对均值回归$a$进行校准，再把校准得到的$a$作为已知量，对波动率$\sigma$进行校准；
* 策略三： 同时校准$a$和$\sigma$。

作者在论文中强调，这三种策略可以同时应用于不同的需求和情景。策略一只对波动率$\sigma$进行校准，具有最大的便捷性和最快的速度，但在适应市场变化方面表现较差。策略三具有最好的拟合质量与稳定性，但校准消耗的时间较长，且不易于实施。策略二可以作为策略三的替代，牺牲了一定的拟合质量，但可以很好地适应市场变化，并缩短运行时间。

本项目选取2023年4月28日至2023年7月24日期间美国市场的Swaption交易数据

## Dynamics
Hull-White模型假设短期利率在风险中性测度下服从:    

$$ dr(t)=\left(\theta(t)-a(t)r(t)\right)dt+\sigma(t)dW^{\mathbb{Q}}(t)$$

	
其中，$\theta(t)$决定了短期利率$r(t)$的长期水平，$a(t)$为短期利率$r(t)$回归到长期水平的速率，简称均值回归，$\sigma(t)$为短期利率的波动率。计算得到短期利率$r(t)$的期望与方差：

$$\mathbb{E}\left[r(t)|\mathcal{F}_{s}\right]=\frac{E(s)}{E(t)}r(s)+\alpha(t)-\frac{E(s)}{E(t)}\alpha(s) $$

$$Var\left[r(t)|\mathcal{F}_{s}\right]=V_{r}(s,t)$$



其中，$E(t)$与$V_{r}(s,t)$具有如下形式：
$$E(t)=e^{\int_{0}^{t}a(u)\,du}  $$
$$V_{r}(s,t)=\frac{1}{E^{2}(t)}\int_{s}^{t}E^{2}(u)\sigma^{2}(u)\,du$$

**零息债券**$P(t,T)$具有仿射结构：
$$P(t,T)=exp\Big(A(t,T)-r(t)B(t,T)\Big)$$
仿射因子$A(t,T)$和$B(t,T)$由以下公式给出：
$$B(t,T)=E(t)\int_{t}^{T}\frac{du}{E(u)}\,du  $$
$$A(t,T)=ln\dfrac{P(0,T)}{P(0,t)} +B(t, T)f(0, t) -\frac{1}{2}B(t,T)^{2}Vr(0, t)  $$

同时，零息债券$P(t,T)$服从对数正态分布：
$$\frac{dP(t,T)}{P(t, T)}= r(t)dt-\sigma(t)B(t, T)dW^{\mathbb{Q}}(t)$$
	
为了得到Swaption价格的闭式解，作者强调零息债券的比率在远期测度下服从：
	
$$d\frac{P(t,T_{F})}{P(t,T_{P})} = \frac{P(t,T_{F})}{P(t,T_{P})}\sigma(t)\big(B(t,T_{P})-B(t,T_{F})\big) dW^{T_{p}}(t)  $$

	
零息债券的比率因此有如下方差：
$$V_{p}(t,T_{F},T_{P})= \int_{t}^{T_{F}}\sigma^{2}(u)\big(B(u,T_{P})-B(u,T_{F})\big)^{2}\,du $$
$$=V_{r}(t,T_{F})B(T_{F},T_{P})^{2}$$
	
## Swaption的解析表达式
### ZBP
首先使用零息债券比率的方差给零息债券看跌期权(zero-bond put option,ZBP)定价:
$$ZBP(T_{F},T_{P},X)=XP(0,T_{F})\mathcal{N}(d_{+})-P(0,T_{P})\mathcal{N}(d_{-})  $$
其中
$$ d_{\pm}=\frac{ln\left(\frac{P(0,T_{F})X}{P(0,T_{P})}\right)}{\sqrt{V_{p}(0,T_{F},T_{P})}} \pm \frac{1}{2}\sqrt{V_{p}(0,T_{F},T_{P})} $$

下面将利用一系列不同期限的零息债券看跌期权给出Swaption的解析表达式。

### Jamshidian分解
考虑一个支付方Swaption，其固定利率为 $K$，到期时间为$T_{0}$，期限为$T_{P}$，以及现金流互换时间$\{T_{i}\}_{i=1,\ldots,n}$，其中$ T_{n}= T_{P}$。可以使用Jamshidian分解将其表示为一系列零息债券看跌期权的加权之和：
$$PSwaption(K,T_{0},T_{P})=\sum_{i=1}^{n}c_{i}ZBP(T_{0},T_{i},X_{i})  $$
$$c_{i}=K\delta(T_{i-1},T_{i})\quad i=1,\ldots,n-1 $$
$$c_{n}=1+K\delta(T_{n-1},T_{n})  $$
$$X_{i}=exp\big(A(T_{0},T_{i})-B(T_{0},T_{i})r^{*}\big)$$

其中，$r^{*}$满足下列等式：
$$\sum_{i=1}^{n}c_{i}exp\big(A(T_{0},T_{i})-B(T_{0},T_{i})r^{*}\big)=1$$


### SMM Approximation
短期利率的波动率有如下近似（SMM估计）：
$$V_{swap}(T_{0},T_{n})=\left[ \frac{P(0,T_{0})}{P(0,T_{0})-P(0,T_{n})}\right]^{2} V_{p}(0,T_{0},T_{n})$$
注意，这个估计是粗糙的，不能代替Jamshidian分解给Swaption定价。但SMM估计有这样的优良性质：具有相同Maturity和不同Tenor之间的SMM估计的比率仅与均值回归$a$有关，这意味着我们可以在波动率$\sigma$未知的情况下校准均值回归$a$。
$$\frac{V_{swap}(M_{i},T_{j})}{V_{swap}(M_{i},T_{k})}=\left[\frac{\Big(P(0,M_{i})-P(0,T_{k})\Big)B(M_{i},T_{j})}{\Big(P(0,M_{i})-P(0,T_{j})\Big)B(M_{i},T_{k})}\right]^{2}$$





### 参考文献
[1] [Calibration Methods of Hull-White Model](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1514192)
[2] [Pricing Swaptions and Coupon Bond Options in Affine Term Structure Models](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=630402)
[3] [The General Hull-White Model and Super Calibration](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1295228)
[4] [F. Jamshidian. An Exact Bond Option Pricing Formula](https://www.jstor.org/stable/2328284)
[5] [Interest rate modelling. Simona Svobod](https://link.springer.com/book/10.1057/9781403946027)
[6] [Leif B. G. Andersen and Vladimir V. Piterbarg: Interest Rate Modeling](https://link.springer.com/article/10.1007/s11408-011-0157-y)
[7] [Using Hull-White Interest Rate Trees](https://www.researchgate.net/publication/228178882_Using_Hull-White_Interest_Rate_Trees)
[8] [Pricing Interest-Rate Derivative Securities](https://www.researchgate.net/publication/5217241_Pricing_Interest-Rate-Derivative_Securities)