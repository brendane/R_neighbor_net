Neighbor networks are an extension of neighbor-joining trees that have
reticulations -- that is they show relationships that arose from
recombination, not just vertical descent. The software package that
makes neighbor nets, SplitsTree, has limited plotting capabilities, and
a unique nexus output format.

These are some R functions that read in the nexus file created by
SplitsTree and allow you to make plots in R. At present, it uses the
layout from SplitsTree, but allows you to adjust the colors, labels,
etc.

Please keep in mind that this code has barely been tested, so
expectations should be low.
