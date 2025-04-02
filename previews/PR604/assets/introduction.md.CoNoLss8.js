import{_ as e,c as a,o as t,a7 as s}from"./chunks/framework.BC9Fsu8r.js";const m=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"introduction.md","filePath":"introduction.md","lastUpdated":null}'),o={name:"introduction.md"},i=s('<h2 id="Rasters.jl" tabindex="-1">Rasters.jl <a class="header-anchor" href="#Rasters.jl" aria-label="Permalink to &quot;Rasters.jl {#Rasters.jl}&quot;">​</a></h2><p><a href="https://rafaqz.github.io/Rasters.jl/dev" target="_blank" rel="noreferrer">Rasters.jl</a> defines common types and methods for reading, writing and manipulating rasterized spatial data.</p><p>These currently include raster arrays like <code>GeoTIFF</code> and <code>NetCDF</code>, <code>R grd</code> files, multi-layered stacks, and multi-file series of arrays and stacks.</p><div class="info custom-block"><p class="custom-block-title">Data-source abstraction</p><p>Rasters provides a standardised interface that allows many source data types to be used with identical syntax.</p><ul><li>Scripts and packages building on Rasters.jl can treat <code>Raster</code>,</li></ul><p><code>RasterStack</code>, and <code>RasterSeries</code> as black boxes.</p><ul><li><p>The data could hold GeoTiff or NetCDF files, <code>Array</code>s in memory or <code>CuArray</code>s on the GPU - they will all behave in the same way.</p></li><li><p><code>RasterStack</code> can be backed by a Netcdf or HDF5 file, or a <code>NamedTuple</code> of <code>Raster</code> holding <code>.tif</code> files, or all <code>Raster</code> in memory.</p></li><li><p>Users do not have to deal with the specifics of spatial file types.</p></li><li><p><code>Projected</code> lookups with Cylindrical projections can by indexed using other Cylindrical projections</p></li></ul><p>by setting the <code>mappedcrs</code> keyword on construction. You don&#39;t need to know the underlying projection, the conversion is handled automatically. This means lat/lon <code>EPSG(4326)</code> can be used seamlessly if you need that.</p></div><h2 id="Installation" tabindex="-1">Installation <a class="header-anchor" href="#Installation" aria-label="Permalink to &quot;Installation {#Installation}&quot;">​</a></h2><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">] add Rasters</span></span></code></pre></div><h2 id="Packages-extensions" tabindex="-1">Packages extensions <a class="header-anchor" href="#Packages-extensions" aria-label="Permalink to &quot;Packages extensions {#Packages-extensions}&quot;">​</a></h2><div class="tip custom-block"><p class="custom-block-title">Packages extensions and Rasters 0.8 and onwards</p><p>On Julia 1.9 we can put additional packages in extensions, so the code only loads when you load a specific package. Rasters.jl was always intended to work like this, and its finally possible. This reduced package <code>using</code> time from many seconds to well under a second.</p><p>But, it means you have to manually load packages you need for each backend or additional functionality.</p></div><p>For example, to use the GDAL backend, and download files, you now need to do:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Rasters, ArchGDAL, RasterDataSources</span></span></code></pre></div><p>where previously it was just using Rasters.</p><p>Sources and packages needed:</p><ul><li><p>:gdal: <code>using ArchGDAL</code></p></li><li><p>:netcdf: <code>using NCDatasets</code></p></li><li><p>:grd: built-in.</p></li><li><p>:smap: <code>using HDF5</code></p></li><li><p>:grib: not yet finished.</p></li></ul><p>Other functionality in extensions:</p><ul><li><p>Raster data downloads, like Worldclim{Climate}: <code>using RasterDataSources</code></p></li><li><p>Makie plots: <code>using Makie</code></p></li><li><p>Coordinate transformations for gdal rasters: <code>using CoordinateTransformations</code></p></li></ul><h2 id="Bugs-and-errors" tabindex="-1">Bugs and errors <a class="header-anchor" href="#Bugs-and-errors" aria-label="Permalink to &quot;Bugs and errors {#Bugs-and-errors}&quot;">​</a></h2><div class="warning custom-block"><p class="custom-block-title">Bugs, errors and making issues for Rasters.jl</p><p>Raster data is complicated and there are many places for subtle or not-so-subtle bugs to creep in.</p><p>We need bug reports to reduce how often they occur over time. But also, we need issues that are easy to reproduce or it isn&#39;t practically possible to fix them.</p><p>Because there are so many raster file types and variations of them, most of the time we need the <em>exact file</em> that caused your problem to know how to fix it, and be sure that we have actually fixed it when we are done. So fixing a Rasters.jl bug nearly always involves downloading some file and running some code that breaks with it (if you can trigger the bug without a file, thats great! but its not always possible).</p><p>To make an issue we can fix quickly (or at all) there are three key steps:</p><ol><li><p>Include the file in an accessible place on web <em>without autentication</em> or any other work on our part, so we can just get it and find your bug. You can put it on a file hosting platform (e.g. google drive, drop box, whatever you use) and share the url.</p></li><li><p>Add a minimum working example to the issue template that first downloads the file, then runs the function that triggers the bug.</p></li><li><p>Paste the complete stack trace of the error it produces, right to the bottom, into the issue template. Then we can be sure we reproduced the same problem.</p></li></ol><p>Good issues are really appreciated, but they do take just a little extra effort with Rasters.jl because of this need for files.</p></div>',17),n=[i];function l(r,d,c,p,u,h){return t(),a("div",null,n)}const f=e(o,[["render",l]]);export{m as __pageData,f as default};
