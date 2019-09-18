args = commandArgs(trailingOnly=TRUE)

template.path <- args[1]
xml.folder <- args[2]
newick.tree.path <- args[3]
rates.path <- args[4]
n.sim <- as.numeric(args[5])

#template.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_template.xml"
#xml.folder <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/xml/"
#newick.tree.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_trees.txt"
#rates.path <- "/Users/rzha419/Workspace/ConstantDistanceOperator/validation/sample_prior/internal_nodes/test_internalnode_rates.txt"
#n.sim <- "2"

trees = readLines(newick.tree.path)
rates = readLines(rates.path)

template.lines = readLines(template.path)
for(scenario.idx in 1:length(trees)){
  tree = trees[scenario.idx]
  rate = rates[scenario.idx]
   for (sim.idx in 1:n.sim) {
     chain.length = paste0(sim.idx,"0000000")
     for (line in template.lines) {
        line = gsub("\\[TreeHere\\]", tree, line)
        line = gsub("\\[ChainLengthHere\\]", chain.length, line)
        line = gsub("\\[RatesHere\\]", rate, line)
        line = gsub("\\[ScenarioHere\\]",scenario.idx, line)
        line = gsub("\\[nSimHere\\]", sim.idx, line)
        write(line, file=paste0(xml.folder,"internalnode_S",scenario.idx,"_",sim.idx,".xml"), append=TRUE)
      }
    }
}