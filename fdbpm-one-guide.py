from fdbpm import FdBpm

fd=FdBpm()
fd.NUM_SAMPLES=101
fd.LENGTH=1E2
fd.dy=1E-4
fd.l_ambda=1.5E-6
fd.L=500E-6

fd.create_space()

print(fd.dx)

# fd.create_source(plotOn=False)
fd.gauss_light()
fd.dn=0.058
fd.n_env=1.004
fd.create_guides()

fd.cmap='gnuplot'
# fd.calculate_propagation()
fd.plot_moving_source()