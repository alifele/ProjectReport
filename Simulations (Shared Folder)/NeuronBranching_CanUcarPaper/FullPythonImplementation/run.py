from pylab import *
from tqdm import tqdm


from main_OLD import BARW

tmax = 200
prob_branch = 0.03
fav = -0.1*(0)## -0.1
fchem = 0.6*0.5

barw = BARW(tmax)
barw.initialize()


for t in tqdm(range(tmax), desc="Simulation Progress"):
    t+=1
    if t<5 and len(barw.node)>1:
        break
    if len(barw.node)!=0:
        barw.tissue1(prob_branch,fav,fchem)
    if len(barw.node)==0:
        break

barw.remove_arcs()
barw.plotter_basic()
barw.plot_graph()




## Useful DEBUG codes
# plt.scatter(self.coordinates[:,0], self.coordinates[:,1])
# plt.show()


# self.coordinates[:52,3] = np.where(self.coordinates[:52,3]==5, 10, self.coordinates[:52,3])