from sympy import *
from KratosMultiphysics import *
from sympy_fe_utilities import *



## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("heat_equation_template.cpp")
outstring = templatefile.read()

for dim in dim_vector:
    if(dim == 2):
        nnodes = 3
    elif(dim == 3):
        nnodes = 4

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity) #N[3x1]

    ## Unknown fields definition
    temp = DefineVector('temp',nnodes)                # TEMPERATURE
    tempn = DefineVector('tempn',nnodes)              # Previous step TEMPERATURE
    tempnn = DefineVector('tempnn',nnodes)            # 2 previous step TEMPERATURE

    ## Other simbols definition
    v = DefineMatrix('v',nnodes,dim)            #  [3x2]
    ## k =  DefineMatrix('k',dim,dim) ##
    k = Symbol('k',positive= True)  # thermal conductivity
    cp = Symbol('cp',positive= True)              #  specific heat at constant pressure
    rho = Symbol('rho', positive = True)        # Density
    dt  = Symbol('dt', positive = True)         # Time increment
    tau = Symbol('tau', positive = True)   #  Stabilization parameter
    Q = DefineVector('Q',nnodes)            # Source term

    ## Backward differences coefficients
    bdf0 = Symbol('bdf0')
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')

    ## Test functions definition
    w = DefineVector('w',nnodes)            # TEMPERATURE field test function

    ## Data interpolation to the Gauss points
    w_gauss = w.transpose()*N
    Q_gauss = Q.transpose()*N
    v_gauss = v.transpose()*N   #[2x1]
    temp_gauss = temp.transpose()*N   # T(x_gauss)
    tempder_gauss = (bdf0*temp + bdf1*tempn + bdf2*tempnn).transpose()*N


    ## Gradients computation (fluid dynamics gradient)
    grad_w = DfiDxj(DN,w) #[1x2]
    grad_temp = DfjDxi(DN,temp) #[2x1]
    div_v = div(DN,v)
    print(div_v.shape)
    # print(grad_temp.shape)

    #Terms definition
    w_convective_term = w_gauss*v_gauss.transpose()*grad_temp
    w_diffusive_term = grad_w * k * grad_temp
    w_temporal_term = w_gauss*tempder_gauss
    w_heat_source_term = w_gauss*Q_gauss

    convective_term = v_gauss.transpose()*grad_temp

    ## Compute galerkin functional
    #  rho*cp*w*dT/dt + rho*cp*w*a *dT + dw * k *dT = w*Q
    rv_galerkin =  w_heat_source_term - rho * cp * w_temporal_term - rho * cp * w_convective_term - w_diffusive_term

    ##  Stabilization functional terms
    # HeatEquation residual       R(T) =  rho*cp[dT/dt + v * dT ] - Q  -> 2nd order terms neglected
    Temperature_residual = Q_gauss - rho*cp*(tempder_gauss + convective_term )

    # Compute the SGS stabilization
    # rv_stab =   tau *rho * cp *(v_gauss.transpose()* grad_w.transpose() + div_v*w_gauss.transpose())* Temperature_residual
    rv_stab =   tau *rho * cp *(v_gauss.transpose()* grad_w.transpose() )* Temperature_residual

    ## Add the stabilization terms to the original residual terms
    rv = rv_galerkin  + rv_stab

    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes, 1))
    testfunc = Matrix( zeros(nnodes, 1))

    for i in range(0,nnodes):

        # TEMPERATURE DOFs and test functions
        dofs[i] = temp[i,0]
        testfunc[i] = w[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation).
    rhs = Compute_RHS(rv.copy(), testfunc, False)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)

    lhs = Compute_LHS(rhs, testfunc, dofs, False)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)


    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)


## Write the modified template
out = open("heat_equation.cpp",'w')
out.write(outstring)
out.close()
