############### SNAQ CODE ############### 
############### Written by: Dylan DeBaun ############### 

################################### RUN SNAQ ############################## 
  
using PhyloNetworks, CSV

group = ARGS[1] # outfile prefix, also the name of the directory with all the input files
numretic = ARGS[2] #number of hybridizations to look for

print(group);
print(numretic);

cd(group);

#read in the tree, mapping, and indiv level concordance factor file (made in prep.R)
astraltree=readTopology(string(group,".tre"));
cf = string(group,"_individualsCFs.csv");
mappingFile = string(group,"_map.csv");

#make a species level cf file, write to folder
df_sp = mapAllelesCFtable(mappingFile, cf);
d_sp = readTableCF!(df_sp);
df_sp = writeTableCF(d_sp);
CSV.write(string(group,"CFtable_species.csv"),df_sp)

#run snaq
x=parse(Int64, numretic)
net = snaq!(astraltree,d_sp, hmax=x, filename=string("net",numretic,group), seed=2345);

################################### ANALYZE SNAQ RESULTS ############################## 
########best to run after snaq is run for all desired reticulation numbers ###### 

group = ARGS[1] #outfile prefix, also the name of the directory with all the input files
out = ARGS[2] #outgroup name
cd(group)

#get libraries
using PhyloPlots,QuartetNetworkGoodnessFit, DataFrames, RCall, PhyloNetworks, CSV

#get species level file and remove repeats (1 representative per species)
df_sp = DataFrame(CSV.File(string(group,"CFtable_species.csv")))	
function hasrep(row)
occursin(r"__2$", row[:t1]) || occursin(r"__2$", row[:t2]) ||
occursin(r"__2$", row[:t3]) || occursin(r"__2$", row[:t4])
df_sp_reduced = filter(!hasrep, df_sp)
CSV.write(string(group,"CFtable_species_norep.csv"), df_sp_reduced);
d_sp_reduced = readTableCF(df_sp_reduced);
d_sp = d_sp_reduced;
qCF = DataFrame!(CSV.File(string(group,"CFtable_species_norep.csv")));

#read in tree and calculate maxq pseudolik and goodness of fit
net0 = readTopology(string("net0",group,".out"));
net0alt =topologyMaxQPseudolik!(net0, d_sp);
ressnaq = quarnetGoFtest!(net0, qCF, false; seed=234, nsim=200);
ressnaq1 = quarnetGoFtest!(net0alt, qCF, false; seed=234, nsim=200);
x=writeTopology(net0alt);

#make data frame for the results (network name, num hybrids, snaq found likelihood, max q likelihood, goodness of fit for snaq found tree, goodness of fit for maximized tree)
dffull = DataFrame(network = "0",numhybrids=0,snaq = net0.loglik, maxq = net0alt.loglik, gof = ressnaq[1], gofalt = ressnaq1[1], topology = "x");

#for all possible reticulation searches (1 to 10)
for j in 1:10
	if(isfile(string("net",j,group,".networks"))) #if there is a reticulation search 
		#read in the topologies
		file = string("net",j,group,".networks")
		networks = readMultiTopology(file);
		#get the loglikelihoods
		scoresInString1 = read(`sed -E 's/.+with -loglik ([0-9]+.[0-9]+).+/\1/' $file`, String)
		scores1 = parse.(Float64, split(scoresInString1));

		#plot the trees, rooting at the outgroup
		@rput group
		@rput j
		R"file=paste0(group,j,\"maxq.pdf\")"
		R"pdf(file,width=20)"
		numnets = length(networks)
		@rput numnets
		R"layout=c(1,numnets)" ###
		R"par"(mar=[0,0,3,0])      # for smaller margins

		for i in eachindex(networks)
			#get rid of hybrids <10% contribution, calculate the new maximized pseudologliklihood, and new goodness of fit value

			print(string(j,".",i,"\n"))
			deleteHybridThreshold!(networks[i], 0.1)
    			x= writeTopology(networks[i]);
    			net1alt = try
        			topologyMaxQPseudolik!(networks[i], d_sp);
        			catch
        			0
        			end
    			loglik = try
        			net1alt.loglik
        			catch
        			0
        			end
    			ressnaq = try
        			#quarnetGoFtest!(networks[i], qCF, false; seed=234, nsim=200);
				catch
        			0
        			end
    			resalt = try
        			quarnetGoFtest!(net1alt, qCF, false; seed=234, nsim=200);
        			catch
        			0
				end
    			ressnaq1 = try
        			ressnaq[1]
        			catch
        			0
        			end
    			resalt1 = try
        			resalt[1]
        			catch
        			0
        			end
			#put all of that info in the data frame
    			push!(dffull, (string(j,".",i),networks[i].numHybrids,scores1[i],loglik,ressnaq1,resalt1,x))

			#for plotting
			try
				rootatnode!(networks[i],out);
			catch e
				0
			end
			                        try
                                rootatnode!(net1alt,out);
                        catch e
                                0
                        end
			try
			plot(net1alt, :R, showGamma=true); ##
			num = string(j,".",i," loglik=",scores1[i]," maxqLL=",loglik," p=",ressnaq1," maxqp=",resalt1);
			@rput num
			R"mtext"(num)           # add text annotation: title here
			catch e
				#plot(networks[i], :R, showGamma=true); ##
                        	plot(net0,:R);
				num = string(j,".",i," loglik=",scores1[i]," maxqLL=",loglik," p=",ressnaq1," maxqp=",resalt1);
                        	@rput num
                        	R"mtext"(num)           # add text annotation: title here
			end
		end
		R"dev.off()"
	else
		#break
	end
	CSV.write(string("results",group,".csv"), dffull) #write to file as it updates
end
CSV.write(string("results",group,".csv"), dffull) #write the final to the file

################################### BOOTSTRAP ############################## 
using PhyloNetworks, CSV
  
#read in the name of the group you are working with (group folder has starting tree, list of bootstrap files in it) and read in number of reticulations to look for

group = ARGS[1]
numretic = ARGS[2]
cd(group)
print(group)

#starting tree
net0 = readTopology(string(group,".tre"));

#list of the locations of bootstrap files (1 file/gene with all BSs in it)
file = string(group,"BSlistfiles");
bootTrees = readBootstrapTrees(file);
file = string(group,"bootsnaq",numretic);
#print("done")
x=parse(Int64, numretic)
print(x)

#nrep=number of bootstraps you have,hmax=number of reticulations to search for
bootnet = bootsnaq(net0, bootTrees, hmax=x, nrep=50, runs=5, filename=file, seed=4321)
