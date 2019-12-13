###################################################################################
# PROGRAM: swccrB.R
# PURPOSE: Function that calculates SW covariate-constrained randomization balance 
#		metric B, given 
#		(1) a permutation matrix of possible cluster orders
#		(2) the SW design scheme - I*J matrix of 0/1 cont/trt
#		(3) matrix of Z-scores of the covariates by cluster
#		(4) nC the number of covariates
#		(5) the desired B cutoff for selection
# AUTHOR: Erin Leister Chaussee
#
# MODIFICATIONS: 3DEC2019 - updates to comments, documentation in preparation for sharing
###################################################################################



###################################################################################
# Balancing function
###################################################################################
# INPUT:
# 	swd - SW design - I cluster by J time period matrix with 0/1 entries 
#		representing control/treatment conditions
#	pm - permutations - data frame where each row represents a cluster order
#		data frame will have I (number of clusters) columns
#	zm - matrix of covariate z-scores for each cluster - one column per covariate
#		columns labeled z1-zC where C is the total number of covariates
#	nC - number of covariates (important if 1 vs. greater than 1)
#	w - weights (default value is 1 for everyone - equal weight)
#	Bp - desired balance metric B cutoff in terms of percentile for sampling
#		ex: Bp=0.20 will create a candidate set from the lowest (best) 20th
#			percentile of the B distribution
#		(default is Bp=0.20)
#
###################################################################################
# OUTPUT: is an object of data type "list"
#	pB_df - data frame of length nperms - contains permutation number, B, and the
#		order of clusters labeled "i1" to "iI" that return this value of B
#	rand_perm - randomly selected permutation with value of B, and order of clusters
# 		based on the desired candidate set size of B specified in the input
#
#	The following are also available for documentation or checking purposes
#	swd - returns the swd input by the user 
#	zm - returns the matrix of Z scores input by the user 
#	sm - contains a data frame of length I*nperms 
#		user can see all permutations and how z-scores get re-ordered before
#		applying the Balance metric calculation
#
###################################################################################

swbf<-function(swd,pm,zm,nC,w=1,Bp=0.20){

	# count number of clusters and proportion in each group by cluster
	nclu<-nrow(swd)
	trtbyclv<-rowSums(swd)/ncol(swd)
	contbyclv<-1-trtbyclv
	
	# prepare values for large matrix/data frame of dimension nrows=nperms*nclusters 
	clustnum<-as.vector(t(pm))
	np<-nrow(pm)
	permnum<-rep(1:np,each=nclu)
	clustord<-rep(1:nclu,times=np)
	ptrtbyclv<-rep(trtbyclv,times=np)
	pcontbyclv<-1-ptrtbyclv
	
	# order z-scores by permuted cluster order
	if(nC==1){
		pzs<-zm[t(pm)]
		}
	else{
		pzs<-zm[t(pm),]
		}
	
	# put all z-scores for all cluster permutations in a dataframe, then into a list
	s<-as.data.frame(cbind(permnum,clustord,clustnum,ptrtbyclv,pcontbyclv))
	sm<-as.data.frame(cbind(s,pzs))
	# split into a list - one list item for each permutation number
	lsm<-split(sm,permnum)

	from<-ncol(s)+1
	to<-ncol(s)+nC

	# apply Balance metric calculations for each permutation
	B_list<-lapply(lsm,function(df){
			if(nC==1){
				ztrt<-df$pzs*df$ptrtbyclv
				zcont<-df$pzs*df$pcontbyclv
				ztrts<-sum(ztrt)
				zconts<-sum(zcont)
				}
			else{
				ztrt<-sweep(df[,from:to], MARGIN=1, df[,"ptrtbyclv"], `*`) 
				zcont<-sweep(df[,from:to], MARGIN=1, df[,"pcontbyclv"], `*`) 
				ztrts<-apply(ztrt,MARGIN=2,sum)
				zconts<-apply(zcont,MARGIN=2,sum)
				}	
			zdiff<-ztrts-zconts
			zdiffsq<-zdiff*zdiff
			wzdiffsq<-w*zdiffsq #apply weights
			B<-sum(wzdiffsq)
			return(B)
		})

	B_df<-do.call("rbind",B_list)
	pB_df<-as.data.frame(cbind("permnum"=c(1:np),B_df,pm))
	names(pB_df)[2] <- "B"

	# select the top Bp percentile of B, then select one permutation at random 
	top_Bp<-pB_df[pB_df$B<=quantile(pB_df$B,probs=Bp),]  
	rand_perm<-top_Bp[sample(nrow(top_Bp),1),]

	retlist<-list(swd=swd,zm=zm,sm=sm,pB_df=pB_df,rand_perm=rand_perm)
	return(retlist)
}


###################################################################################
# EXAMPLES OF USING THE FUNCTION
# Need SW design, cluster order permutations, z-score matrix, 
#	and sample weights (weights optional)
###################################################################################

##########################################
# SW example design
##########################################
swd67<-matrix(c(0,1,1,1,1,1,1,
		    0,0,1,1,1,1,1,
		    0,0,0,1,1,1,1,
		    0,0,0,0,1,1,1,
		    0,0,0,0,0,1,1,
		    0,0,0,0,0,0,1),nrow=6,ncol=7,byrow=TRUE)

##########################################
# Matrix of cluster order permutations 
##########################################
library(gtools) # need for 'permutations' function
pt6<-as.data.frame(permutations(n = 6, r = 6, v = 1:6, repeats.allowed=FALSE))
names(pt6)<-c("i1","i2","i3","i4","i5","i6")
head(pt6)

##########################################
# Sample Z-score matrix - 4 covariates
##########################################
z1<-c(0.7293, -1.1024, -0.7971,  1.1363,  0.8311, -0.7971)
z2<-c(-1.1777,  1.2108,  0.8127,  0.6137, -0.8791, -0.5805)
z3<-c(-1.3005, -0.1432,  0.1644, -0.9123,  1.2129,  0.9787)
z4<-c(1.1142, -0.3562, -0.2063, -1.2979, -0.5286,  1.2748)
szmi<-cbind(z1,z2,z3,z4)

##########################################
# Sample weights to use for a 4 covariate example below
##########################################
wts<-c(0.7,0.1,0.1,0.1)



####################################################
# EXAMPLE 1: 1 covariate
####################################################
set.seed(1)
e1<-swbf(swd=swd67,pm=pt6,zm=szmi[,1],nC=1,Bp=0.20)
e1$zm
head(e1$sm,n=12)
head(e1$pB_df)
e1$rand_perm
hist(e1$pB_df$B)


####################################################
# EXAMPLE 2: 4 covariates
####################################################
set.seed(2)
e2<-swbf(swd=swd67,pm=pt6,zm=szmi,nC=4,Bp=0.20)
e2$swd
e2$zm
head(e2$sm,n=12)
head(e2$pB_df)
e2$rand_perm
hist(e2$pB_df$B)


####################################################
# EXAMPLE 3: 4 covariates, weight covariates based
#		on importance in balancing
####################################################
set.seed(3)
e3<-swbf(swd=swd67,pm=pt6,zm=szmi,nC=4,w=wts,Bp=0.20)
e3$swd
e3$zm
head(e3$sm,n=12)
head(e3$pB_df)
e3$rand_perm
hist(e3$pB_df$B)

