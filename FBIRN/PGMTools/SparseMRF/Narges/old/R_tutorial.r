# You can download R from here :
#   http://cran.cnr.berkeley.edu/
# You might try the Tinn-R editor found here
#   http://sourceforge.net/projects/tinn-r
# The editor will color the text making it easier to read
# and will let

# When you start R you can paste any of the lines
# above and they will do nothing.  This is because
# the character '#' is a comment character and 
# R ignores everything after it

# Now we can assign variables using either the "<-"
# operator or the "=" operator

x <- 5.5
y = 3:10
#length(y) = # elements in the vector y
z <- 1:length(y)

print(y)  # 3 4 5 6 7 8 9 10
print(z)  # 1 2 3 4 5 6 7 8

# --------------------------------------------------------------------
# Selecting elements of a vector
# --------------------------------------------------------------------


# The statement  n:m  creates a vector which
# counts from n to m.  It's useful for selecting
# elements of a vector

print(y[2:4])  # 4 5 6

# We can also select elements using logical expressions
# For example

ind1 <- y >= 4
print(ind1)  # FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

ind2 <- y < 7
print(ind2)  # TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE

print(ind1 && ind2)            # FALSE

print(ind1 & ind2)             # FALSE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE

print(y[ind1 & ind2])          # 4 5 6
print(y[(y >= 4) & (y < 7)])   # 4 5 6

# The other logical command that is useful is the OR, '|'
print( y[(y == 3) | (y >= 6)])   # 3  6  7  8  9 10



# --------------------------------------------------------------------
# Loops and if-else ladders
# --------------------------------------------------------------------

# For more complicated things we may use a 'for' loop

ind1 <- c()     # creates an empty vector
ind2 <- c()

# The vector 'y' consists of length(y) = 8 elements
# we will loop over the integers 1 through 7

for(i in 1:(length(y)-1)){

	# What's inside this 'if' statement only executes
	# if it the condition is true.  

	if( ((y[i]+y[i+1]) %% 3) == 0){

		# This command creates a new vector composed of
		# 'ind1' and 'i'
		ind1 <- c(ind1, i)

		# We are also allowed to 'push back' a new element
		# by assigning to the index one larger than the
		# current length
		ind2[length(ind2)+1] <- i

	}else if( y[i] %% 2 == 0){

		# The condition above is only checked if the 
		# previous one evaluated to FALSE
		print( c(y[i], i) )

	}else{

		# We only reach this point of none of the conditions above
		# evaluated to TRUE
		print(i)
	}
}
print(ind1)    # 3 6
print(ind2)    # 3 6


# We could have also done this in one line with
ind3 <- which( ((y[-length(y)] + y[-1]) %% 3) == 0)

# The command y[-j] returns every element of 'y' except
# the j'th.  The command 'which' takes a logical vector
# and returns the indices that are equal to 'TRUE'. To
# match with the 'for' loop above we have

# For example
(1:10)[-5]    # 1  2  3  4  6  7  8  9 10 


# --------------------------------------------------------------------
# Matrices
# --------------------------------------------------------------------

# We create a matrix of zeroes with 5 rows and 6 columns
m <- matrix(nrow=5, ncol=6, data=0)
print(dim(m))
print(m)

m[1, -1] <- 2:ncol(m)
print(m)

x <- matrix(nrow=ncol(m), ncol=1, data=1:ncol(m))
print(x)

# We perform matrix multiplication as follows
y <- m %*% x

# We can extract elements of a matrix using commands like
z <- m[1:3, 2:4]
print(dim(z))
print(z)

# In R if you extract multiple rows, the object you get
# back is a matrix.  Likewise if you select multiple columns
# If you select just a single row or column, you get back
# a vector instead
z <- m[1, 2:4]
print(dim(z))    # NULL        .... oops  
print(z[1,2])    # Error in z[1, 2] : incorrect number of dimensions
print(z)         #  2    3    4

# If you think this might happen to you, add the 'drop' argument
# at the end.  Sometimes you know you're selecting just one
# row/column and you really do want a vector back
z <- m[1, 2:4, drop=FALSE]
print(dim(z))    # 1 3 
print(z[1,2])    # 3
print(z)         
#         [,1] [,2] [,3]
#   [1,]    2    3    4

# There are other commands

# Extract a square sub matrix
z <- m[1:5, 1:5]
# Set the diagnonal to all 1's
diag(z) <- 1
print(z)

# Compute the inverse
z.inv <- solve(z)
print(z.inv)

z.eig <- eigen(z)
print(z.eig)

z.tall <- rbind(z, z.inv)
print(z.tall)

z.wide <- cbind(z.inv, z)
print(z.wide)




# --------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------

# You can define your own function using the syntax below

my.func <- function(x, y=1, z=1){
	return( x * y / z )
}

# We can call this several ways 
my.func(3)         # 3
my.func(x=3)       # 3
my.func(3, z=2)    # 1.5
my.func(x=3, z=2)  # 1.5

# We can assign functions to variables

g <- my.func

g(x=3,y=2)   # 6

# Finally, the following will not work :
my.func(y=3) # Error in my.func(y = 3) : argument "x" is missing, with no default



# --------------------------------------------------------------------
# Help pages and some examples of plotting
# --------------------------------------------------------------------

# Beyond loops, if-else statements, and functions the next most useful
# thing to know is the help system.  The following commands won't make
# sense until you look at the documentation

x <- runif(30)
y <- rnorm(n=30, sd=0.3) + x

regObj <- lm(y ~ x)
summary(regObj)
#plot the data
plot(x, y, xlab="predictor", ylab="response")
#Add the regression line to the plot
abline(a=regObj$coefficients[1], b=regObj$coefficients[2])

# None of the above should make sense, but if you have a function
# name but you don't know how it works you can get documentation for
# it inside R with the '?' command.  For example
?rnorm

# The help page indicates that this is a function for the normal distribution
# it has documentation for 4 functions : dnorm, rnorm, pnorm, qnorm
# A lot of the random variable functions follow this pattern
# The 'd' prefix means the function will return the density of the distribution
# The 'r' prefix indicates that the function generates random variables
# and the 'p' prefix returns the CDF of the distribution

# Next try the lm command.
?lm

# From the documentation you can see the variables names, but this isn't 
# always useful.  The 'Details' section explains the first argument
# 'formula' and how to use it.  The 'Examples' section of these pages
# is often the best way to understand the functions

# Create a plotting window with 1 row and 3 columns
par(mfrow=c(1,3))

hist(rnorm(n=100, mean=0.5, sd=1), xlab="x", ylab="count", main="a hist" )

# Create a sequence of 100 numbers beginning at -4 and ending at +4
x <- seq(from=-4, to=4, length=100)

plot(x, dnorm(x, mean=0.5, sd=1), xlab="value", ylab="density", main="PDF" )
plot(x, pnorm(x, mean=0.5, sd=1), xlab="val", ylab="P(Z < x)", main="CDF" )

# --------------------------------------------------------------------
# Reading Files
# --------------------------------------------------------------------

# If you are using Linux you can start R from any directory and it
# automatically save/load files from that directory.  If you are using
# Windows, R will look in some obscure, useless place in your C drive
# So set the working directory first

# setwd("F:/classes/202")
dir()  # 'character(0)'

m <- matrix(nrow=500, ncol=4, data=runif(500*4) )

write.table(m, file="m.csv", sep=",", col.names=c("a", "b", "c", "d"),
            row.names=FALSE)
            
dir()    # m.csv

# If you want to save your entire state you can run the function
save.image()    # m.csv

# It looks like nothing was saved.  That's because the 'dir' function
# doesn't print filenames that begin with a period '.'

# Let's look at the documentation
?dir

# There is an argument called all.files
dir(all.files=TRUE)       "."      ".."     ".RData" "m.csv"

# There are other arguments you can try
dir(all.files=TRUE, pattern=".*csv")   # m.csv

# Now we can read the table with the command
m2 <- read.table(file="m.csv", header=TRUE, sep=",")

# When we saved the matrix we gave the write.table function
# some column names.  We can use these to extract columns
z <- m2$a    
  
# We had to know that there was a header row in the 
# table when reading it.  If we had instead written
m3 <- read.table(file="m.csv", header=FALSE, sep=",")
# Then the next command wouldn't work
print(m3$a)     # NULL

x <- 1:4
y <- m %*% x   # This works
y <- m2 %*% x  # This doesn't
# Error in m2 %*% x : requires numeric matrix/vector arguments

# The object m2 is a 'data frame', so we can convert it
# to a numeric matrix using the 'as.matrix' function
m3 <- as.matrix(m2) 
y <- m3 %*% x  # Now it works

# --------------------------------------------------------------------
# Scripts
# --------------------------------------------------------------------

# you can save any sequence of R functions and commands in a script (like
# this file).  You can run it as follows:

# setwd("current_directory")
source("Filename.R");

# for more in-depth programs, this is nice so you don't have to keep copying
# and pasting all the time.

# Open the script editor in R and type the following into a file:

blender = function(x, y, z, w) {
	num1 = x / y - z / w;
	num2 = w / z - y / x;
	num3 = x / w - z / y;
	num4 = y / z - w / x;
	num5 = sqrt(num1 * num2) - (num3 - num4)^2;
	blender = rep(num5, times = num2); 
	return(blender);
}

# now save it as "blender.R"
# at the command prompt type 

# setwd("current_path")
source("blender.R");
my_blender = blender(7, 4, 1, 15);
print(my_blender);

# good...our new function seems to work.  Or does it?

# --------------------------------------------------------------------
# Debugging
# --------------------------------------------------------------------

# try again with different numbers:

my_blender = blender(7, 4, 0, 15);
print(my_blender);

# We get the following error:

# Error in rep(num5, times = num2) : invalid 'times' argument
# In addition: Warning message:
# In blender(7, 4, 0, 15) : NAs introduced by coercion

# First, we look over the code to see what it is doing.  We are 
# creating five numbers num1 - num5 using arithmatic, then we are
# using the rep command, which seems to be causing the error.

# first we look at the help page for rep:

?rep

# we see that rep(x, times = y) repeats the element x y times.  But how can
# we see what value of y was causing the error, since we can't "see" into the
# internal workspace of the function "blender"?

# We can use the following command:

options(error = recover);
my_blender = blender(7, 4, 0, 15);
print(my_blender);

# if you set this option, when an error occurs, you will be directed to 
# a menu that will allow you to enter the work space of any of the
# functions in the trace list that have been called at the time the error 
# occured.  This mode is called browser mode.

# Enter the number of the function that you would like to look inside.  
# When you are inside a function, you can type the name of a variable 
# to view its contents.  You can even execute any R command while in browser
# mode.

# in this case the blender function is the only one on the list.  Type '1'.

num1 # 1.75  
num2   
# we get Inf: this presents a problem since we're trying to repeat by an 
# infinite number of times!
x # 7
y # 4
z # 0: the reason we got Inf is that we were trying to divide by 0.

# To exit browser mode and go back to the menu, type <enter>.  
# Type '0' to exit the menu and return to the command prompt.

# now try

my_blender = blender(7, 4, 100, 15);
print(my_blender);

# we get numeric(0).  What in the world is going on here?  Shouldn't we be getting
# a vector?  Also, since our code didn't actually generate an error, the 
# options(error = recover) command won't help us.  What do we do?

# we can modify our file by putting a browser() command in it an re-save it:

blender = function(x, y, z, w) {
	num1 = x / y - z / w;
	num2 = w / z - y / x;
	num3 = x / w - z / y;
	num4 = y / z - w / x;
	num5 = sqrt(num1 * num2) - (num3 - num4)^2;
	blender = rep(num5, times = num2); 
	browser(); ## ADDED!
	return(blender);
}

# R will direct you to browser mode whenever it reaches a browser()
# statement.  The following commands are useful for navigation:

# 'n': execute one line of code at a time
# 'c' or <enter>: exit browser and continue running (if you are in a loop,
# you will enter browser mode again next time you come to this loop)
# 'Q': exit browser mode and terminate the R script

# now try:  

source("blender.R");
my_blender = blender(7, 4, 100, 15);

# at the Browse[1] prompt, type:

blender # numeric(0), sure enough
num1 #  -4.916667
num5 # -501.6868
num2 #  -0.4214286

# What?  We are trying to repeat -501.6868 by -0.4214286 times, and R is not
# generating an error?  This is an example of quirky behavior by R, and we 
# need to be careful of it.

# what R is most likely doing is rounding -0.4214286 to the nearest integer,
# aka 0, and returning a vector of length 0: numeric(0) is how R prints a
# zero-length integer.

# We may wish to perform additional checks on num2 to make sure it is a
# non-negative integer:


blender = function(x, y, z, w) {
	num1 = x / y - z / w;
	num2 = w / z - y / x;
	num3 = x / w - z / y;
	num4 = y / z - w / x;
	num5 = sqrt(num1 * num2) - (num3 - num4)^2;
	if(num2 > 0 && num2 == floor(num2)) {
		blender = rep(num5, times = num2);
	} else {
		blender = NaN;
		warning("Invalid 'times' argument");
	}
	browser();
	return(blender);
}

source("blender.R");
my_blender = blender(7, 4, 100, 15);
print(my_blender);

# now we get a descriptive "warning" message and my_blender is NaN, which
# is R's code for an invalid numeric quantity.

# you can look up what the floor() function does on your own.

# Since our program now seems to work, we now want to restore R's default 
# error-handling behavior:

options(error = NULL);

# We also can comment out the browser() command so the function runs straight
# through.

