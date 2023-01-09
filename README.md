# gemsnake-reticulations
The following procedure was created in DeBaun et al., 2022 to comprehensively explore the network structure and trends for the evolutionary history of Malagasy gemsnakes.

Classify the quality of the gene tree dataset
Create subclades of tree to identify reticulations (test monophyly of those subdivisions)
  Based on the results of this test, remove taxa or gene trees that are uninformative
Subsample at least one member from each subclade to test for reticulations deeper in the tree

For each subclade and the subsampled backbone:
Run SnaQ for up to 5 reticulations  <br />
Run NANUQ, look at strucutre with SplitsTree  <br />
Test all reticulation options with MaxQPseudoLikelihood and Goodness of Fit test  <br />
Run HyDe on concatonated gene trees (or SNPs) to see if reticulations are supported  <br />




## Classify quality of dataset, identify possible areas of reticulations

We used [PhyParts](https://bitbucket.org/blackrim/phyparts) to carry out this analysis. <br />
GENETREE is the directory to rooted gene trees <br />
SPTREE is the rooted species tree file <br />
out is the outfile <br />

```
java -jar target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -d $GENETREE -m $SPTREE -o $out -v
```
To visualize, I used the python script from mossmatters/PhyPartsPieCharts which uses ete3 (conda installed) to graph the product

```
python phypartspiecharts.py $SPTREE $out $NUM_GENES --show_nodes --svg_name $out.svg
```
Red and green indicates the level of discordance at nodes and could be an indicotor of several processes like ILS or reticulation
Gray indicates the proportion of genes that are uninformative, and could be problematic for a reticulation study. One could remove those uninformative genes when running the analysis, or remove a taxon that seems to be causing the issue. For example, here, we could remove either *Madagascarophis colubrinus B* or *C* to improve the analysis (a caveat being that we cannot infer reticulations among these members)
<img width="468" alt="Screen Shot 2022-10-09 at 9 17 33 PM" src="https://user-images.githubusercontent.com/37371274/194788419-047827dd-346b-4942-b786-21a754138551.png">

## Create subclades and test monophyly
We needed to break the tree up into smaller clades so that the network inference programs could handle them. The largest clade had 19 species, 92 indiviudals. We used [DiscoVista](https://github.com/esayyari/DiscoVista) to analyze if the clades chosen supported monophyly, meaning members in the clades were not jumping into other clades, which could be a sign of reticulation. We wanted each clade we tested to have high support for monophyly.

Note: we used the docker install of DV but it can also be run locally. See their github for more information

First, create a clade definition file so you know which members belong to each clade. This function will contribute the simple tab file (species \t clade_name) into the definition file.
cladeannotation.txt
```
sp1, clade1
sp2, clade1
... ...
```

```
docker run -v ${DV_directory}:/data esayyari/discovista generate_clade-defs.py parameters/cladeannotation_${testname}.txt parameters/clade-defs_${testname}.txt
```
Then we can run discoVista and look at the results.
```
docker run -v /Users/dylandebaun/Desktop/example/1KP:/data esayyari/discovista discoVista.py -c parameters/clade-defs_${testname}.txt -p genes/ -t 75 -m 1 -o genes/gtdiscordresults${testname}/ 
```
High levels of dark and light blue indicate that you are good to go, there is little chance of between clade reticulation. If you see high levels of Red and yellow, you may want to adjust your subclades to comprise more taxa, test out different combinations. 
If the yellow and red proportion goes down signifigantly by adding more members to the clade, you may have a reticulation among those members, and will want to test that entire clade together.
If yellow and red goes down by removing members from a clade, you may have a lot of discordance as the result of a deep divergence, you'll see that Parasteophis and Leioheterodon support monophyly when separated but not when together.
If no grouping metric seems to decrease the proportion of red and yellow, you may have a deep time reticulation present. Here, we couldn't explain the great discordance among Pseudoxyrhopus, Heteroliodon (and Liopholoidophis when added in), this was the result of a deep time reticulation that was found by subsampling.

<img width="500" alt="Screen Shot 2022-10-09 at 9 37 37 PM" src="https://user-images.githubusercontent.com/37371274/194789273-67b90a05-2bf8-41a7-8c2c-ca161c3b2aeb.png">
<img width="296" alt="Screen Shot 2022-10-09 at 9 40 43 PM" src="https://user-images.githubusercontent.com/37371274/194789386-004ffeef-f8aa-4c10-886f-ac1ed7178f91.png">

### You now have the set of clades you will be working with. You should additionally test the backbone of the tree for reticulations. We did this by sampling one or two members from each subclade and running the network analysis on those members.

## Run SnaQ
We will be running SnaQ on a bifurcating tree to identify 1 to (say about) 5 reticulations. 
We need to create all the prep files for a subclade: a file of concordance factors for 

To create for each clade--
x_individuals.txt = list of individuals in clade <br />
x_species.txt = list of species in clade <br />
x_map.csv = mapping of individuals to species (where individuals listed in allele column) <br />
```
allelle, species
... ...
```
**prep_SnaQ_NANUQ.R** 
This script creates the necessary quartet input files for SnaQ and NANUQ

**SnaQ_run.jl**
This script runs SnaQ.

It creates several reticulations, whose results for pseudolikelihood and goodness of fit will be out into a results*.csv file. These topologies can then be compared (see Compare the Networks section).

Note: you can also run bootstraps in SnaQ as an added support metric for your chosen topology

## Run NANUQ

For NANUQ, we will calculate the probabilities that each quartet is a quartet (p_T3) or a star (p_star). you will need to choose your own threshold values for each of thes (alpha for p_T3 and beta for p_star). 

Choice in alpha and beta should be made intentionally, identifying quartet pairs that group together at low values of alpha or high values of beta. This can be difficult to assess when you have many indivduals, so I created a script to concatonate individual information together to look at it at the species level.

**NANUQ.R**
This script runs NANUQ to calculate the p_T3 and p_star. Using that information, create a file called cladename_alpha.cdv; cladename_beta.csv.

The .nex files in NANUQ can be visualized using SplitsTree
'Data > Filter Taxa' can be used to remove problematic individuals (for example those with many uninformative genes as mentioned above) and to simplify the results to a species tree.
<img width="474" alt="Screen Shot 2022-10-09 at 11 20 14 PM" src="https://user-images.githubusercontent.com/37371274/194794970-7cfae374-0021-46cc-b33a-d0fe89543588.png">

Reticulations can be seen in darting behaviors like those mentioned in the ms for NANUQ (Allman et al., 2019). X-0 is the hybrid and the two Xis on either side are the parent lineages. 
<img width="888" alt="Screen Shot 2022-10-09 at 11 22 07 PM" src="https://user-images.githubusercontent.com/37371274/194795076-19fe67de-4ad0-447d-badc-07b0672d8510.png">

Using this information, you can convert the network into Extended Newick format for networks and calculate the pseudolikelihood and goodness of fit for these values.

## Compare the Networks

Metrics to compare these topologies include the [goodness of fit](https://cecileane.github.io/QuartetNetworkGoodnessFit.jl/dev/man/gof/) and pseudolikelihood in [PhyloNetworks](https://crsl4.github.io/PhyloNetworks.jl/latest/). See how these are calculated in **run_snaq.jl** file

Depending on your question you can remove reticulations with certain genomic contributions (i.e. only looking at introgression of >10%) in PhyloNetworks
```
deleteHybridThreshold!(network, 0.1)
```

Calculate the MaxQ Pseudolikihood for a network/tree which will optimizae the branch lengths/genomic contributions and calculate the pseudolikelihood. 
```
topologyMaxQPseudolik!(network, species_level_cf)
```

Calculate the GoF for a network/tree. This will return a p-value, explaining the propotion of quartets that explain the topology as well as an uncorrected z-score and sigma value.
```
quarnetGoFtest!(net1alt, qCF, false; seed=234, nsim=200);
```

## Run HyDe
All of the above methods required gene tree estimation to run. To counteract this, we run [HyDe](https://github.com/pblischak/HyDe) which is a site specific hybrid detection software and explore if our hybrids are signifigant. We ran this on all 109 species at once, but it can be done in broke down groups also.
1. Concatonate all genes together (or better yet: use SNP dataset)
2. create taxon map 
3. Run HyDe
```
run_hyde.py -i concatonated_phylip.txt -m map.txt -o out -n num_indivs -t num_taxa -s num_sites
```
4. Calclulate the proportion of signifigant triplets out of the total triplets
Look at the triplets involving extant parent lineages (those that arised after hybridization event) and extant hybrid lineages. Caclulate the proportion of these triplets that are signifigant.

## Build Reticulation_through_Time Plot

**reticulation_through_time.R** builds the reticulation through time plot. the rate of accumulation was calculated as the slope of the line connecting these points 

<img width="321" alt="image" src="https://user-images.githubusercontent.com/37371274/211403111-2a107ba1-51a5-410e-8664-2172820cf343.png">


## References
PhyParts - https://bitbucket.org/blackrim/phyparts
DV - https://github.com/esayyari/DiscoVista
NANUQ - Allman, E. S., Baños, H., & Rhodes, J. A. (2019). NANUQ: a method for inferring species networks from gene trees under the coalescent model. Algorithms for Molecular Biology, 14(1), 1-25.
SnaQ - https://crsl4.github.io/PhyloNetworks.jl/latest/
GoF - https://cecileane.github.io/QuartetNetworkGoodnessFit.jl/dev/
HyDe - https://github.com/pblischak/HyDe


## Citation
DeBaun et al., 2022. Widespread Reticulate Evolution in an Adaptive Radiation. Evolution. *accepted*

