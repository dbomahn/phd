using Pkg
cd(Pkg.devdir()) #take us to the dev dir

# pkg console
# generate MyPkg
# dev MyPkg
# edit(pathof(MyPkg))


#pkg console
# activate .
import vOptPkg
vOptPkg.greet()
# pkg console
# add JuMP
# add DataFrames
# ...
