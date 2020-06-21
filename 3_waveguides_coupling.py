#%%
from fdbpm import FdBpm
%matplotlib qt


#%%
# config #1: successful coupling (20s)
fd=FdBpm()
fd.NUM_SAMPLES=101
fd.LENGTH=3E4
fd.dy=1E-6
fd.l_ambda=0.85E-6
fd.L=60E-6
fd.create_space()
plotOn=True
offset=14E-6
# fd.create_source(waist=5E-6,offset=-offset,plotOn=False)
fd.gauss_light(fwhm=8E-6,offset=-offset)
fd.create_guides(width=10E-6,offset=-offset)
fd.create_guides(width=10E-6)
fd.create_guides(width=10E-6,offset=offset,plotOn=True)
# fd.calculate_propagation()
# fd.cmap='prism'

#%%
# config #2: faster coupling (7s)
fd=FdBpm()
fd.NUM_SAMPLES=101
fd.LENGTH=1.2E4
fd.dy=1E-6
fd.l_ambda=0.85E-6
fd.L=60E-6
fd.create_space()
plotOn=True
offset=14E-6
# fd.create_source(waist=5E-6,offset=-offset,plotOn=False)
fd.gauss_light(fwhm=8E-6,offset=-offset)
fd.create_guides(width=10E-6,offset=-offset)
fd.create_guides(width=10E-6)
fd.create_guides(width=10E-6,offset=offset,plotOn=False)
#%%
# config #3: even faster coupling (2s)
fd=FdBpm()
fd.NUM_SAMPLES=101
fd.LENGTH=500
fd.dy=1E-5
fd.l_ambda=1.55E-6
fd.L=60E-6
fd.create_space()
plotOn=True
offset=14E-6
# fd.create_source(waist=5E-6,offset=-offset,plotOn=False)
fd.gauss_light(fwhm=12E-6,offset=-offset+1E-7)
fd.create_guides(width=10E-6,offset=-offset)
fd.create_guides(width=10E-6)
fd.create_guides(width=10E-6,offset=offset,plotOn=False)

#%%
fd.plot_moving_source()

# %%
