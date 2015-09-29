f2py -m elem -c ElemIntgl0.f90 mesh.f90 extend_mesh.f90 mod_func.f90 tripole_mod.f90 ./hi_intgl/hi_integral.f90 ./hi_intgl/hi_const.f90 ./hi_intgl/hi_tfunc.f90 ./hi_intgl/hi_funcs.f90 
rm *.so.*
