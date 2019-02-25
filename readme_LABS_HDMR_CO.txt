The function to call is

[score0 score1 fa_vec, miss_vec]=LabsHdmrCo(xtr0,xtr1,xts0,xts1,T1,T2,T3,tvec,GFvec)

xtr0 is train data in class 0
xtr1 is train data in class 1
xts0 is test data with true label being class 0
xts1 is test data with true label being class 1

xtr0, xtr1, xts0, xts1 should contain dosage values, i.e., each element should be 0, 1, or 2

xtr0 is an n0 by F matrix, where n0 is sample size in class 0 and F is the total number of SNPs
xtr1 is an n1 by F matrix, where n1 is sample size in class 1

xts0 is an n0test by F matrix where n0test is test sample size in class 0
xts1 is an n1test by F matrix where n1test is test sample size in class 1

T1, T2, and T3 are threhsolds described in the main paper. In short,
T1 is the threshold for the risk of a single SNP to enter the classification rule
T2 is the threshold for the risk of a snp pair to enter the classification rule
T3 is the threshold for the risk of a snp block to enter the classification rule

tvec is the vector of threhsolds on the total risk function R(X) to obtain the ROC curves

For each element of tvec, tvec(i) yhat=R(X)>tvec(i) is the assigned label.

GFvec is the vector describing the number of SNPs to use in the analysis.

score0 is the score of test points in class 0 (for the last GFvec value).
score1 is the score of test points in class 1 (for the last GFvec value).

fa_vec is a length(GFvec) by length(tvec) matrix where element in row i and column j
is the probability of false alarm (false positive)
when GFvec(i) features and threshold tvec(j) are used

miss_vec is a length(GFvec) by length(tvec) matrix where element in row i and column j
is the probability of miss (false negative)
when GFvec(i) features and threshold tvec(j) are used
