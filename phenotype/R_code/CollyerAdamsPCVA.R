# The following function estimates vector lengths, contrast, angles, and
# associated p values for differences between traits in multivariate space
# between two pairs of populations. Code adapted from Collyer & Adams (2007)

CollyerAdamsPCVA <- function(y.mat, x.mat, nPermutations = 9999) {
  # first step: finalise the two x matrices (full and reduced)
  x.mat.full<-x.mat
  x.mat.red<-x.mat.full[,-4] # This removes the coding for the interaction effect
  
  # Second step, estimate parameters for full and reduced models
  b.mat.full<-solve(  (t(x.mat.full)%*%x.mat.full)  ) %*% (t(x.mat.full)%*%y.mat)
  b.mat.red<-solve((t(x.mat.red)%*%x.mat.red))%*%(t(x.mat.red)%*%y.mat)

  # Third step, Estimate LS means

  # In this case, but not necessarily always, a = species and b = locality 
  # Vectors will be calculated for 2 species, a1 (P. jordani) and a2 (P. teyahalee)
  # In our example the following dummy variables were used:
  # P. jordani = 1, P. teyahalee = -1, sympatry = 1, allopatry = -1
  # Thus, LS means could be calculated as follows:

  a1b1<-cbind(1,1,1,1)
  a1b2<-cbind(1,1,-1,-1)
  a2b1<-cbind(1,-1,1,-1)
  a2b2<-cbind(1,-1,-1,1)

  x.ls.full<-rbind(a1b1,a1b2,a2b1,a2b2)
  x.ls.red<-x.ls.full[,-4]

  obs.ls.full<-x.ls.full%*%b.mat.full  # Observed ls means (full model)
  obs.ls.red<-x.ls.red%*%b.mat.red     # Observed ls means (reduced model)

  # Fourth Step, vector and statistics calculations

  obs.a1.vect<-obs.ls.full[1,]-obs.ls.full[2,] # These are the phenotypic change vectors
  obs.a2.vect<-obs.ls.full[3,]-obs.ls.full[4,]

  obs.d.a1<-c(sqrt(t(obs.a1.vect)%*%obs.a1.vect)) # These are lengths of vectors
  obs.d.a2<-c(sqrt(t(obs.a2.vect)%*%obs.a2.vect))

  obs.contrast<-abs(obs.d.a1-obs.d.a2)
  obs.angle<-c(acos(t((obs.a1.vect)/obs.d.a1)%*%((obs.a2.vect)/obs.d.a2)))
  obs.angle<-obs.angle*180/pi  # This step is only necessary to convert radians to degrees


  # Fifth Step, set-up permutation procedure

  y.hat<-x.mat.red%*%b.mat.red     # Predicted values from reduced model
  y.res<-y.mat-y.hat               # Resdiuals of reduced mode (these are the permuted units)


  # PERMUTATION PROCEDURE

  # Need to set-up distributions to be generated

	dist.d1<-NULL
	dist.d1<-rbind(dist.d1,obs.d.a1) # Observed value is first random value
	dist.d2<-NULL
	dist.d2<-rbind(dist.d2,obs.d.a2) # Observed value is first random value
	dist.contrast<-NULL
	dist.contrast<-rbind(dist.contrast,obs.contrast) # Observed value is first random value
	dist.angle<-NULL
	dist.angle<-rbind(dist.angle,obs.angle) # Observed value is first random value

	# In addition to saving random values, it is wise to save the outcome of comparisons
	# of observed and random values.
	# This can be done with logical statements (below).
	# Separate distributions are created for these comparisons.
	# The 'p' indicates that these distributions will be used 
	# to calculate empirical probabilities.
	
	pdist.contrast<-NULL
	pdist.contrast<-rbind(pdist.contrast,1) 
	pdist.angle<-NULL
	pdist.angle<-rbind(pdist.angle,1) 


  for(i in 1:nPermutations){
	  y.res.rand <- y.res[sample(1:nrow(y.res)), ]
	
  	# Create random values
  	y.rand<-y.hat+y.res.rand

  	# Estimate parameters
  	b.mat.rand<-solve((t(x.mat.full)%*%x.mat.full))%*%(t(x.mat.full)%*%y.rand)

  	# Calculate LS means
  	rand.ls.full<-x.ls.full%*%b.mat.rand

  	# Repeat fourth step for random data!

  	rand.a1.vect<-rand.ls.full[1,]-rand.ls.full[2,] # These are the phenotypic change vectors
	  rand.a2.vect<-rand.ls.full[3,]-rand.ls.full[4,]

	  rand.d.a1<-c(sqrt(t(rand.a1.vect)%*%rand.a1.vect)) # These are lengths of vectors
	  rand.d.a2<-c(sqrt(t(rand.a2.vect)%*%rand.a2.vect))

	  rand.contrast<-abs(rand.d.a1-rand.d.a2)
	  rand.angle<-c(acos(t((rand.a1.vect)/rand.d.a1)%*%((rand.a2.vect)/rand.d.a2)))
	  rand.angle<-rand.angle*180/pi  # This step is only necessary to convert radians to degrees

	  # Append distributions

	  dist.d1<-rbind(dist.d1,rand.d.a1)
	  dist.d2<-rbind(dist.d2,rand.d.a2)
	  dist.contrast<-rbind(dist.contrast,rand.contrast) 
	  dist.angle<-rbind(dist.angle,rand.angle) 
	
	  aa<-ifelse(rand.contrast>=obs.contrast,1,0)
	  bb<-ifelse(rand.angle>=obs.angle,1,0)
	
	  pdist.contrast<-rbind(pdist.contrast,aa) 
	  pdist.angle<-rbind(pdist.angle,bb) 
	} 

  # Empirical probabilities are calculated as follows

  p.contrast<-sum(pdist.contrast)/(nPermutations+1)
  p.angle<-sum(pdist.angle)/(nPermutations+1)

  return(c(length.v1 = obs.d.a1,
           length.v2 = obs.d.a2,
           contrast = obs.contrast,
           p.contrast = p.contrast,
           angle = obs.angle,
           p.angle = p.angle))
}