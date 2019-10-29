import pandas as pd
df = pd.read_csv('./noload.dat', header=None)
print(df)


from pylab import plt, np
x = np.array(df[3].values)**2
y = df[5].values
# The first two points are not in no-load condition
x = x[2:]
y = y[2:]
print(x,y)
# The first several points in no-load condition are not really in linear region due to leakage inductance
x = x[9:] # 多使用低电压的点，可以对线性部分进行更准确地拟合
y = y[9:] # 多使用低电压的点，可以对线性部分进行更准确地拟合
y += 100 # [W] friction and windage loss
print(x,y)
plt.plot(x, y, 'ko-')


import statsmodels.api as sm
# Ordinary least squares regression
print('Start OLS...')
model_Simple = sm.OLS(y, x).fit()
print(model_Simple.summary())
print('Parameters: ', model_Simple.params)

# Add a constant term like so:
print('\n', '-'*40)
print('Start OLS with a constant term (intercept)...')
model = sm.OLS(y, sm.add_constant(x)).fit()
print(model.summary())
print('Parameters: ', model.params)

x = [0] + x.tolist()
print(x)

plt.plot(x, np.array(x)*model_Simple.params[0], 'bs--')
plt.plot(x, np.array(x)*model.params[1] + model.params[0], 'rs--')

plt.xlabel('Voltage Squared [$\\rm V^2$]')
plt.ylabel('Power [W]')

plt.grid()
plt.show()


