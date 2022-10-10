# gemsnake-reticulations
The following procedure was created in DeBaun et al., 2022 to comprehensively explore the network structure and trends for the evolutionary history of Malagasy gemsnakes.

Classify the quality of the gene tree dataset
Create subclades of tree to identify reticulations (test monophyly of those subdivisions)
  Based on the results of this test, remove taxa or gene trees that are uninformative
Subsample at least one member from each subclade to test for reticulations deeper in the tree

For each subclade and the subsampled backbone:
Run SnaQ for up to 5 reticulations
Run NANUQ, look at strucutre with SplitsTree
Test all reticulation options with MaxQPseudoLikelihood and Goodness of Fit test
Run HyDe on concatonated gene trees (or SNPs) to see if reticulations are supported




## Classify quality of dataset, identify possible areas of reticulations

We used [PhyParts] (https://bitbucket.org/blackrim/phyparts) to carry out this analysis. <br />
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
We needed to break the tree up into smaller clades so that the network inference programs could handle them. The largest clade had 19 species, 92 indiviudals. We used [DiscoVista] (https://github.com/esayyari/DiscoVista) to analyze if the clades chosen supported monophyly, meaning members in the clades were not jumping into other clades, which could be a sign of reticulation. We wanted each clade we tested to have high support for monophyly.

Note: we used the docker install of DV but it can also be run locally. See their github for more information

First, create a clade definition file so you know which members belong to each clade. This function will contribute the simple tab file (species \t clade_name) into the definition file:

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


##References
PhyParts
DV


## Citation
DeBaun et al., 2022. Widespread Reticulate Evolution in an Adaptive Radiation. Evolution. *in revision*

