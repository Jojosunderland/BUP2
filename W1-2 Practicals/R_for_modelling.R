my.function <- function(x=0, y=0) {
  z <- x*y
  return(z)
  }

my.value <- my.function(x=2)

my.value <- my.function(x=2, y=3)

install.packages("vegan")
library("vegan")

# vectors 
x <- c(1,3,5)
y <- c("a", "b", "c")

x <- 1:5

# matrix and array
# 2D vector is a matrix

X <- matrix(0, nrow=2,ncol=3)

# for more than 2 dimensions use an array

X <- array(0, dim=c(2,3,5))

#data frame

my.data <- data.frame(x=1:5, y=c("a", "b", "c", "d", "e"))
my.list<- list(x=x, X=X, my.data=my.data)


# for loops and if conditions
# a for loop is a way to iterate through a suite of elements

# we can add numbers 1-10 by writing
x <- 0
for(i in 1:10) {
  x <- x+i
}
x

# or we can create words as:
y <- c("a","b","c")
z <- character(3) 
##we create a vector of characters of size 3 but with no values in it
ii <- 0
for(i in y){
  ii <- ii+1
  z[ii] <- paste(i,i)
}

# An if condition is a function that will lead to execute some lines of codes only if a condition
# is met. For example, we can call if() within a for() loop as follows:

x <- 0
for(i in 1:10){
  x <- x + i
  if(x==10){
    print(x)
  }
}

# apply(), lapply() and sapply() can be used instead of for loops
# apply() applies function FUN to input x on dimensions MARGIN, and returns an element of the same type as x
# sum the columns
m1 <- matrix(1:10,nrow=5, ncol=6)
m1
a.m1 <- apply(m1, 2, sum)

#lapply() applies function FUN to input list X, and returns a list
species <- c("DOG","CAT","EAGLE","ADER")
species_lower <-lapply(species, tolower)
str(species_lower)

# sapply() applies function FUN to input list X, and returns a vector.
dt <- cars ## this is a dataset that is included with base R
lmn_cars <- lapply(dt, min)
smn_cars <- sapply(dt, min)
dt
lmn_cars
smn_cars
a.m1