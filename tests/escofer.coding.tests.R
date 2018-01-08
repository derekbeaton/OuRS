rm(list=ls())
data(beer.data)


escofier.transform_v1 <- escofier.coding(beer.data,T,scale="SS1")



therm_v1 <- thermometer.coding(beer.data,norm.to.one = F)

  ## all of these should end up as rowSums = 16
therm_v2 <- thermometer.coding(beer.data,norm.to.one = T)
therm_v3 <- thermometer.coding(beer.data,norm.to.one = T,mins=rep(0,ncol(beer.data)))
therm_v4 <- thermometer.coding(beer.data,norm.to.one = T,maxs=rep(5,ncol(beer.data)))
therm_v5 <- thermometer.coding(beer.data,norm.to.one = T,mins=rep(0,ncol(beer.data)),maxs=rep(5,ncol(beer.data)))


# DATA <- beer.data
#   ## puts everything at 0
#   ### alt: I could put in a vector for the mins.
#   mins <- apply(DATA,2,min)
#   maxs <- apply(DATA,2,max)
#
#    #these.maxs <- rep(5,ncol(beer.data))
#   #these.mins <- rep(0,ncol(beer.data))
#
# from.mins <- sweep(DATA,2,mins,"-")
# from.mins.normed <- sweep(from.mins,2,these.maxs,"/")
#
# from.maxs <- sweep(from.mins,2,these.maxs,"-") * -1
# from.maxs.normed <- sweep(from.maxs,2,these.maxs,"/")
#
# test.dat <- cbind(from.mins,from.maxs)
# test.dat.normed <- cbind(from.mins.normed,from.maxs.normed)
#
#
# rowSums(test.dat.normed)
# rowSums(test.dat)




