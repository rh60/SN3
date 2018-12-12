function write_vtkfile(fname::String,dataname::String,
  x::Vector{Float64},y::Vector{Float64},Tri::Matrix{Int},
  U::Vector{Float64})
    nt=size(Tri,1);
    C = Array{MeshCell,1}(undef,nt);
    for i=1:nt
      T=Tri[i,:]
      C[i] = MeshCell(VTKCellTypes.VTK_TRIANGLE, T)
    end
    vtkfile = vtk_grid(fname,x,y,C)
    vtk_point_data(vtkfile,U, dataname)
    vtk_save(vtkfile)
end
