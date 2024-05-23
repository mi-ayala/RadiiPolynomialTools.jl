 using FiniteDifferences

function FiniteDiff_forsequences(x, F)
    """
        A function to compute derivatives of sequence functions
        with finite differences.
   
    """
        Space = space(x)
   
        g = y ->  F(Sequence(Space,y) .+ x)
   
        g_r = x -> real(g(x)[:])
        g_i = x -> imag(g(x)[:])
        print("First FD")
        DF_r = jacobian(central_fdm(16, 1), g_r, zeros(dimension(Space)))[1]
        print("Second FD")
        DF_i = jacobian(central_fdm(16, 1), g_i, zeros(dimension(Space)))[1]
   
        return  LinearOperator(Space, Space, DF_r + DF_i*im)
        #return  DF_r + DF_i*im
       
end
