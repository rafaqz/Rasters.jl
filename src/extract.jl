"""
   extract(vsr::VerySimpleRaster, x, y)
   extract(vsr::VerySimpleRaster, points)

Extracts the value of an array at the given points.

TODO: integrate this with GeometryBase.jl or similar
"""
extract(A, args...) = extract(A, (args...))
extract(A, tup) = begin
   x, y = tup[1], tup[2]
   indx, indy = try
      coord_to_cell(vsr, x, y)
   catch
      return missing
   end
   vsr[indx, indy]
end

extract(A, tup::Tuple{T,T}) where T <: AbstractVector = extract.(Ref(vsr), zip(tup...))
extract(A, ar::AbstractVector) = extract.(Ref(vsr), ar)
