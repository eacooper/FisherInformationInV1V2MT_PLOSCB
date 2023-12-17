function[params,minerr] = fitGeneralizedLaplacian( x, f )
    
    minerr = 1e99;
    for p = 0.01 : 0.05 : 5
        for s = 0.01 : 0.05 : 5
            gl = exp( -(abs(x)/s).^p );
            gl = gl / sum(gl);
            err = mean((gl-f).^2);
            if( err < minerr )
                minerr = err;
                params = [s p];
            end
        end
    end