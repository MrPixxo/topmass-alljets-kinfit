import  pyKinFit
inputpt = [[309.6 , 168, 161.8 , 128.4, 82, 76.4]]
inputeta =[[-1.28, -2.09, -1.52, -0.65, -1.55, -1.52]]
inputphi=[[-0.06, -2.83, -2.72, 0.1, 2.93, 0.67]]
inputM=[[24.9, 14, 21.2, 13.2, 12.4, 13.3]]
gen = [[275.2, 194.4, 168.2, 145.1, 91.1 , 36.5],[-1.28, -2.1,-1.51, -0.65, -1.5, -1.47],[-0.06, -2.85, -2.7, 0.08, 2.93, 0.75],[1.5, 4.8, 0.1, 4.8, 0.1, 0.2]]
fitPt, fitEta, fitPhi, fitMass ,indexmask, fitChi2, fitPgof = pyKinFit.setBestCombi(inputpt, inputeta, inputphi, inputM)
print(fitMass[0])
