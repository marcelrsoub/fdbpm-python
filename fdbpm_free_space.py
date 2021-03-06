from fdbpm import FdBpm

fd=FdBpm()
fd.NUM_SAMPLES=101
fd.LENGTH=1E2
fd.dy=1E-4
fd.l_ambda=1.5E-6
fd.L=1000E-6

fd.create_space()

# fd.create_source(plotOn=False)
fd.gauss_light(fwhm=20E-6)
# fd.dn=0.001
# fd.n_env=1
fd.create_guides(width=0)

fd.cmap='gnuplot2'
# fd.calculate_propagation()
fd.plot_moving_source()