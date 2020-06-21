from fdbpm import FdBpm

fd=FdBpm()

fd.create_space()

fd.create_source(plotOn=False)

fd.cmap='cool_r'
# fd.calculate_propagation()
fd.plot_moving_source()