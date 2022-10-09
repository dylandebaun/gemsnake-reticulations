# gemsnake-reticulations

Classify the quality of the gene tree dataset
Create subclades of tree to identify reticulations (test monophyly of those subdivisions)
  Based on the results of this test, remove taxa or gene trees that are uninformative
Subsample at least one member from each subclade to test for reticulations deeper in the tree

For each subclade and the subsampled backbone:
Run SnaQ for up to 5 reticulations
Run NANUQ, look at strucutre with SplitsTree
Test all reticulation options with MaxQPseudoLikelihood and Goodness of Fit test
Run HyDe on concatonated gene trees (or SNPs) to see if reticulations are supported
