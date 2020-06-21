from fdbpm import FdBpm

fd=FdBpm()
fd.NUM_SAMPLES=100
fd.LENGTH=3E4
fd.dy=1E-6
fd.l_ambda=0.85E-6
fd.L=60E-6
fd.create_space()
plotOn=True
offset=14E-6
fd.create_source(waist=5E-6,offset=-offset,plotOn=False)
fd.create_guides(offset=-offset)
fd.create_guides()
fd.create_guides(offset=offset,plotOn=False)
# fd.calculate_propagation()
fd.cmap='prism'
fd.plot_moving_source()