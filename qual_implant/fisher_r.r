
# define the 3x2 contingency table
x = matrix(c(1,7,3,1,4,1),nrow=3,ncol=2)   


# do fisher test
fisher.test(x)