# define packages to be used
 using JuMP
 using Ipopt
 using DataFrames
 using CSV

 my_content ="""Animal type,Colour,Cover
 bunny,white,fur
 dragon,green,scales
 cow,brown,fur
 pigeon,grey,feathers"""

 open("animals.csv", "w") do out_file
    # write will return the number of bytes written to the file
    # that's what that 90 is about
    write(out_file, my_content)
 end
 animals = CSV.read("animals.csv")
 animals = readtable("animals.csv")
 # iris = CSV.read("iris.csv")

 iris=readtable("/Users/Antoine/Documents/McGill Fall 2018/iris.csv") #This works much better as it takes care of spaces

 names(iris)
 names(animals)

 iris_sub = iris[end,2:3]
 animals_sub = animals[end,3]

 iris_vector=convert(DataFrame, iris[end,2:3])
 iris_vector
 iris[:Species]
 iris_vector_magic=iris[:Sepal_Length]

 iris_vector_magic[2]

 bidon=readtable("CSV-Bidon.csv")
 bidon_vector_Watts=bidon[:Amps]



t=[1, 2, 3, 4, 5];
i=[1, 2]

sum(c_g[i] * g[i,t] + c_w * w[t] for i=1:N for t=1:T)
