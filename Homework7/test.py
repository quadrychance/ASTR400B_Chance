from M33_Analytic import M33AnalyticOrbit
M33 = M33AnalyticOrbit('test')
a = M33.OrbitIntegrator(0,12,.1)

print(a)
