# I don't use directed diffusion distance since it doesn't
# guarantee that the initial states are the most-uniform

# we might want to use Markov-logic here eventually?
# say: smooth the estimated markov-transition matrix via rank-K
# and then use absorption probabilities? 
chromatin_potential_postprocess <- function(chrom_obj){
  # form the metacell weighted graph where the weights are
  # the exp-correlation-proportions.
  # for any missing edges that do not involve any initial/terminal
  # cells, put exp(-2) [the minimum possible weight]
  
  # under the assumption that cells should be more differented
  # as time proceeds, use the following logic:
  # set the fate-probability vector of all terminal states
  # to be an indicator vector, and the initial state to be
  # 1/K vector, and all other cells NA
  
  # then keep iterating under the two rules:
  # we are trying to find a vector for each metacell such that:
  # 1) each metacell's vector is more-uniform than the
  # weighted sum of its out-arrows
  # 2) each metacell is less-uniform than the weighted sum
  # of its in-arrows
  
  # this proceeds using in the following steps:
  # order all the cells based on the min-distance to any initial/terminal state
  # this is the order we'll deal with cells
  
  # using the posited ordering, initialize via backpropogating, where we 
  # slowly flip the NAs
  # and set each metacell equal to the weighted sum of its outarrows
  # (this is an iterative process)
  
  # then, based on the iteration, repeat the following steps:
  # MODIFYING
  # using the posited ordering, for each compute the weighted
  # average of its in-arrows and average of its out-arrows
  # this gives 2 K-dimension vectors, and the cell has its own K-dimensional
  # vector
  # perform some optimization so that within the values
  # set by the 2-K dimension vectors AND is sandwiched in terms of
  # it's non-uniformity
  
  # TERMINATION?
  # see if our fate-prob vectors satsifies our system-of-equations
  # if not, repeat the above process again
  
}