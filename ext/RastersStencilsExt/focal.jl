
"""
    focal(f, x::AbstractRaster; stencil=Window{1}())

Broadcast a function of stencils over a raster.

`f` could be, for example `mean`, to implement a `mean` blur.
"""
focal(f, x::AbstractRaster; stencil=Window{1}()) = Stencils.broadcast_stencil(f, x, stencil)
