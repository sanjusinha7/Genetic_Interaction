#This script predicts genetic interations using GO realtions as features.
#Sanju Sinha

#Loading required Libraries
require(ranger)
require(randomForest)
require(RcppCNPy)

#Functions Required


#Loading Files needed
GO=lapply(readLines('/cbcb/project2-scratch/sanju/project1c-predicting-genetic-interactions-master/data/examples/example-hierarchy-sets.tsv'), function(x) strsplit(x, '\t'))
GI=npyLoad('/cbcb/project2-scratch/sanju/project1c-predicting-genetic-interactions-master/data/examples/example-genetic-interactions.npy')

feature_vector_fora_pair <- function(GA, GB){
	sapply(GO, function(x) sum(!is.na(match(x[[1]], c(paste(GA), paste(GB))))))
}


GeneList=sort(unique(unlist(GO)))
GeneList=GeneList[order(as.numeric(sapply(GeneList, function(x) strsplit(x, '-')[[1]][2] ) ) )]
Pairs=data.frame(GeneA=rep(GeneList, each=100), GeneB=rep(GeneList, 100))
Pairs=Pairs[(Pairs[,1]!=Pairs[,2]),]
Pairs=Pairs[as.numeric(sapply(Pairs[,1], function(x) strsplit(as.character(x), '-')[[1]][2] )) < as.numeric(sapply(Pairs[,2], function(x) strsplit(as.character(x), '-')[[1]][2] )),]

#Feature vector for every pair
Input_RF=data.frame(Score_List=GI[lower.tri(GI, diag=F)], FV=t(apply(Pairs, 1, function(x) feature_vector_fora_pair(x[1], x[2]))))

samp <- sample(nrow(Input_RF), 0.6 * nrow(Input_RF))
train <- Input_RF[samp, ]
test <- Input_RF[-samp, ]
#Number of combination:: Pairs=(99*100/2)= 4950
#RF_Model
RF_Model=randomForest(Score_List~., train)

#Prediction using the above mdoel
pred <- predict(RF_Model, newdata = test)

##Pearson Co-relation
pear.cor= cor.test(test$Score_List, pred, method='pearson')
print('with a statistical significance of', pear.cor$p.value, 'our prediction has pearson co-relation(without any qq-plot check) of', pear.cor$estimate)


#________________Let's use Random forest as Binary classifier for interaction method______________#
## ***Hyperparameter to play around with.****
K=3
##
xtile_Input_RF=data.frame(Score_List=as.factor(xtile(Input_RF$Score_List, K)), Input_RF[,-1])
samp <- sample(nrow(xtile_Input_RF), 0.6 * nrow(xtile_Input_RF))
xtile_train <- xtile_Input_RF[samp, ]
xtile_test <- xtile_Input_RF[-samp, ]
#Number of combination:: Pairs=(99*100/2)= 4950
#RF_Model
xtile_RF_Model=randomForest(Score_List~., xtile_train)

#Prediction using the above mdoel
xtile_pred <- predict(xtile_RF_Model, newdata = xtile_test)

##Truth box
table(xtile_pred, xtile_test$Score_List)

