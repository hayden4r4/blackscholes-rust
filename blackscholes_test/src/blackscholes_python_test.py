import blackscholes_python as bs

inputs = bs.Inputs(option_type=bs.OptionType.Put, s=100, k=110, p=10, r=0.05, q=0.05, t=20.0/365.25, sigma=0.2)
price = inputs.calc_price()
delta = inputs.calc_delta()
gamma = inputs.calc_gamma()
vega = inputs.calc_vega()
theta = inputs.calc_theta()
rho = inputs.calc_rho()
iv = inputs.calc_iv(0.0001)

print(f"Price: {price}")
print(f"Delta: {delta}")
print(f"Gamma: {gamma}")
print(f"Vega: {vega}")
print(f"Theta: {theta}")
print(f"Rho: {rho}")
print(f"Implied Volatility: {iv}")

# print("Done: ", price)

# inp = bs.Inputs(bs.OptionType.Call, 100, 100, None, 0.05, 0.02, 50/365.25, 0.2)
# # print(inputscalc_price(inp))
# print(inp)
